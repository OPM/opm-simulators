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

#ifndef OPM_WELL_PERFORMANCE_EVENTS_HEADER_INCLUDED
#define OPM_WELL_PERFORMANCE_EVENTS_HEADER_INCLUDED

#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <opm/output/data/Wells.hpp>

#include <map>
#include <string>
#include <vector>

namespace Opm
{

class Schedule;

/// Dynamic open/shut status of the wells owned by the current rank.
///
/// The status a well and its connections have after the simulator has had
/// its say, i.e., the deck status amended by the run time decisions recorded
/// in WellTestState and WellState.
struct WellStatusSnapshot {
    struct Entry {
        WellStatus status {WellStatus::SHUT};

        /// Completion number of every connection that is currently able to
        /// flow, one entry per connection and sorted numerically for set
        /// operations.  Completion numbers need not follow the wellbore
        /// ordering used by +CON, and COMPLUMP may map several connections
        /// onto one completion, so this is a multiset: the repeats are what
        /// make WPWE1 and WPWE2 count connections rather than completions.
        /// Do not deduplicate.
        std::vector<int> openCompletions {};

        /// Sorted completion numbers currently closed by an actual +CON
        /// workover, rather than merely configured with +CON limits.
        std::vector<int> closedByConPlus {};

        /// True when no connection is able to flow and every one of them was
        /// closed by a 'CON' workover -- the manual's first route to WPWE3.
        /// Derived here because only the snapshot sees every connection.
        bool closedToBottomByCon {false};

        bool operator==(const Entry& rhs) const
        {
            return (this->status == rhs.status) && (this->openCompletions == rhs.openCompletions)
                && (this->closedByConPlus == rhs.closedByConPlus)
                && (this->closedToBottomByCon == rhs.closedToBottomByCon);
        }

        template <class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(status);
            serializer(openCompletions);
            serializer(closedByConPlus);
            serializer(closedToBottomByCon);
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
/// Only automatic status and connection changes are reported.  The tracker
/// forms those event counters by comparing a snapshot of the dynamic
/// status taken at the start of a time step -- after the deck has been
/// applied for the step -- against the status at the end of the step.  Deck
/// driven changes (WELOPEN, COMPDAT, WCONPROD, ...) are therefore absorbed
/// by the snapshot rather than reported.  Injector/producer conversions are
/// the exception, and are tracked whether they come from a schedule keyword
/// or from ACTIONX, since they have no other source.  The manual states no
/// such distinction; it has not been checked against a reference run.
///
/// The counters cover a single time step -- they are discarded when the next
/// one starts -- so a summary written after a step reports exactly the events
/// of that step, which is what the WPWE indicators denote.  A well test that
/// re-opens a well and the limit check that closes it again in the same step
/// are both reported, rather than netted against each other.
///
/// WPWE0, the drilled indicator, is always zero: it reports the wells a
/// drilling queue brings on stream, and QDRILL is not supported by Flow.
class WellPerformanceEvents
{
public:
    /// Start a new time step.
    ///
    /// Discards the counters of the previous step and adopts \p snapshot as
    /// the reference status.  Call once the deck has been applied for the
    /// step, so that the deck driven changes are absorbed by the reference
    /// rather than reported as events.
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

    /// Accept the injector/producer types after a successful time step.
    /// Call before applying any ACTIONX changes for the next time step.
    /// Failed attempts must not advance this baseline: beginTimeStep() will
    /// then report the same conversion again when the step is retried.
    void commitTimeStep(const Schedule& schedule, int reportStep);

    /// Adopt \p snapshot as the reference status without recording events.
    ///
    /// Used to step past the deck driven changes of a report step.
    void synchronise(WellStatusSnapshot snapshot);

    /// Record the changes between the reference status and \p current, and
    /// adopt \p current as the new reference.
    ///
    /// \param[in] current Dynamic status at the end of a successful time step.
    void accumulate(WellStatusSnapshot current);

    /// Event counters accumulated for \p wellName in the current time step.
    ///
    /// Returns a zeroed record for a well that has seen no events.
    const data::WellEvents& events(const std::string& wellName) const;

    static WellPerformanceEvents serializationTestObject();

    bool operator==(const WellPerformanceEvents& rhs) const
    {
        return (this->events_ == rhs.events_) && (this->previous_ == rhs.previous_)
            && (this->injector_ == rhs.injector_);
    }

    template <class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(events_);
        serializer(previous_);
        serializer(injector_);
    }

private:
    /// Event counters for the current time step, by well name.
    std::map<std::string, data::WellEvents> events_ {};

    /// Status the wells had at the start of the current time step.
    WellStatusSnapshot previous_ {};

    /// Injector/producer flag at the last accepted time step, by well name.
    std::map<std::string, bool> injector_ {};
};

} // namespace Opm

#endif // OPM_WELL_PERFORMANCE_EVENTS_HEADER_INCLUDED
