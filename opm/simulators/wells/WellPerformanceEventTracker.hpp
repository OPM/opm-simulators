/*
  Copyright 2026 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_WELL_PERFORMANCE_EVENT_TRACKER_HEADER_INCLUDED
#define OPM_WELL_PERFORMANCE_EVENT_TRACKER_HEADER_INCLUDED

#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <opm/output/data/Wells.hpp>

#include <map>
#include <string>
#include <vector>

namespace Opm
{

class Schedule;

/// Dynamic status of the wells and connections owned by the current rank.
///
/// The deck status amended by runtime decisions recorded in WellTestState and
/// WellState.
struct WellStatusSnapshot {
    struct Entry {
        WellStatus status {WellStatus::SHUT};

        /// Completion number of every connection that is currently able to
        /// flow, one entry per connection and sorted numerically for set
        /// operations.  COMPLUMP may map several connections onto one
        /// completion, so this is a multiset: the repeats are what make WPWE1
        /// and WPWE2 count connections rather than completions, while WPWE3
        /// works on the distinct values.  Do not deduplicate.
        std::vector<int> openCompletions {};

        /// Sorted completion numbers currently closed only by the reach of a
        /// '+CON' workover, rather than by a limit of their own.  WPWE2 does
        /// not count them, and a workover that reaches past its offender has
        /// closed the well to the bottom.
        std::vector<int> closedBelowOffender {};

        /// Completion number of the well's topmost connection, whether or not
        /// it is able to flow.  WPWE3 uses this to distinguish bottom-up
        /// closure from a well whose only flowing connection still has closed
        /// connections above it.
        int topCompletion {0};

        bool operator==(const Entry& rhs) const
        {
            return (this->status == rhs.status) && (this->openCompletions == rhs.openCompletions)
                && (this->closedBelowOffender == rhs.closedBelowOffender)
                && (this->topCompletion == rhs.topCompletion);
        }

        template <class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(status);
            serializer(openCompletions);
            serializer(closedBelowOffender);
            serializer(topCompletion);
        }
    };

    std::map<std::string, Entry> wells {};

    bool operator==(const WellStatusSnapshot& rhs) const
    {
        return this->wells == rhs.wells;
    }

    template <class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(wells);
    }
};

/// Accumulates the well and connection status changes that back the WPWE0 to
/// WPWE7 summary vectors.
///
/// The tracker forms the event values by comparing snapshots of dynamic status
/// as the time step progresses.  Connection changes are reported regardless
/// of source: an economic-limit workover, a well test, or a deck keyword such
/// as WELOPEN or COMPDAT.  The exception is connections closed only by the
/// reach of a '+CON' workover below the offending connection; WPWE2 excludes
/// those closures.  Well-status changes are reported only when made by the
/// simulator: a well shut or stopped by the deck is not a WPWE4 or WPWE7
/// event.  Injector/producer conversions come from the schedule and are
/// tracked whether introduced directly by a keyword or through ACTIONX.
///
/// The values cover a single time step -- they are discarded when the next
/// one starts -- so a summary written after a step reports exactly the events
/// of that step, which is what the WPWE indicators denote.  A well test that
/// reopens a well and a limit check that closes it again in the same step
/// are both reported, rather than netted against each other.
///
/// WPWE0, the drilled indicator, is not supported yet and stays zero.  It
/// reports a well taken off the drilling queue defined by QDRILL and WDRILTIM.
/// Because QDRILL is unsupported, no run reaches this tracker with a drilling
/// queue.  Introducing a well through WELSPECS is not a drilling event.
class WellPerformanceEventTracker
{
public:
    /// Start a new time step.
    ///
    /// Discards the counters of the previous step and adopts \p snapshot as
    /// the reference status, reporting the connection changes the schedule
    /// has made since the last step.
    ///
    /// Injector/producer switches (WPWE5, WPWE6) are checked on every step,
    /// including ACTIONX changes within a report step.
    ///
    /// \param[in] schedule Simulation schedule.
    /// \param[in] reportStep Zero-based index of the current report step.
    /// \param[in] snapshot Dynamic status at the start of the time step.
    void beginTimeStep(const Schedule& schedule,
                       int reportStep,
                       WellStatusSnapshot snapshot);

    /// Accept the current snapshot and injector/producer types after a
    /// successful time step.
    /// Call before applying any ACTIONX changes for the next time step.
    /// Failed attempts must not advance this baseline: beginTimeStep() will
    /// then report the same conversion again when the step is retried.
    void commitTimeStep(const Schedule& schedule, int reportStep);

    /// Record the connection changes the deck has just made, and adopt
    /// \p snapshot as the new reference status.
    ///
    /// Well status changes are absorbed rather than reported: WPWE4 and
    /// WPWE7 denote the simulator's own decisions, not the deck's.
    void applyDeckChanges(WellStatusSnapshot snapshot);

    /// Preserve a simulator-driven status change across time-step retries.
    ///
    /// A convergence failure may close a well and commit that decision before
    /// the time step is retried.  The next beginTimeStep() sees the resulting
    /// status in its starting snapshot, where ordinary status differences are
    /// treated as deck changes.  Keep the event separately until the retried
    /// time step is accepted so WPWE4 or WPWE7 is still reported.
    void recordWellStatusChangeForRetry(const std::string& wellName, WellStatus status);

    /// Record the changes between the reference status and \p current, and
    /// adopt \p current as the new reference.
    ///
    /// \param[in] current Dynamic status at the end of a successful time step.
    void accumulate(WellStatusSnapshot current);

    /// Event values accumulated for \p wellName in the current time step.
    ///
    /// Returns a zeroed record for a well that has seen no events.
    const data::WellPerformanceEvents& events(const std::string& wellName) const;

    static WellPerformanceEventTracker serializationTestObject();

    bool operator==(const WellPerformanceEventTracker& rhs) const
    {
        return (this->events_ == rhs.events_) && (this->previous_ == rhs.previous_)
            && (this->lastAccepted_ == rhs.lastAccepted_)
            && (this->hasLastAccepted_ == rhs.hasLastAccepted_)
            && (this->injector_ == rhs.injector_)
            && (this->pendingStatusEvents_ == rhs.pendingStatusEvents_);
    }

    template <class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(events_);
        serializer(previous_);
        serializer(lastAccepted_);
        serializer(hasLastAccepted_);
        serializer(injector_);
        serializer(pendingStatusEvents_);
    }

private:
    /// Event values for the current time step, by well name.
    std::map<std::string, data::WellPerformanceEvents> events_ {};

    /// Most recent status observed during the current time-step attempt.
    WellStatusSnapshot previous_ {};

    /// Status at the end of the last accepted time step.  Failed attempts are
    /// always compared from this baseline again.
    WellStatusSnapshot lastAccepted_ {};

    bool hasLastAccepted_ {false};

    /// Injector/producer flag at the last accepted time step, by well name.
    std::map<std::string, bool> injector_ {};

    /// Simulator-driven status indicators that must survive failed attempts.
    std::map<std::string, data::WellPerformanceEvents> pendingStatusEvents_ {};

    /// Shared implementation of accumulate() and applyDeckChanges().
    void accumulate(WellStatusSnapshot current, bool trackStatus);
};

} // namespace Opm

#endif // OPM_WELL_PERFORMANCE_EVENT_TRACKER_HEADER_INCLUDED
