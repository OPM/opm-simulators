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

#include <config.h>

#include <opm/simulators/wells/WellPerformanceEventTracker.hpp>

#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <algorithm>
#include <iterator>
#include <utility>

namespace Opm
{

void
WellPerformanceEventTracker::beginTimeStep(const Schedule& schedule,
                                     const int reportStep,
                                     WellStatusSnapshot snapshot)
{
    this->events_ = this->pendingStatusEvents_;

    if (! this->hasLastAccepted_) {
        // First step of a run or restart: there is no earlier dynamic status
        // to compare against, and no recorded well types.  Adopt the wells
        // that exist now, so that a conversion during this step is reported
        // against the type they have here.
        this->previous_ = std::move(snapshot);
        this->commitTimeStep(schedule, reportStep);
        return;
    }

    // Restore the accepted baseline before every attempt.  A failed attempt
    // may have advanced previous_, while the well state itself is rolled back
    // by resetWGState().
    this->previous_ = this->lastAccepted_;

    // A new report step brings its schedule with it, so the snapshot may
    // already hold connections the deck has opened or closed since the last
    // accepted step.  Repeating the comparison from lastAccepted_ preserves
    // those events when the time step is retried.
    this->accumulate(std::move(snapshot), false);

    for (const auto& wellName : schedule.wellNames(reportStep)) {
        const auto pos = this->injector_.find(wellName);
        if (pos == this->injector_.end()) {
            // A well introduced by WELSPECS while the run is under way has no
            // previous type to compare against.  Entering the schedule is not
            // a drilling event: WPWE0 reports the drilling queue, which is not
            // supported yet.
            continue;
        }

        const auto isInjector = schedule.getWell(wellName, reportStep).isInjector();
        if (pos->second != isInjector) {
            auto& events = this->events_[wellName];
            (isInjector ? events.producerToInjector : events.injectorToProducer) = 1;
        }
    }
}

void
WellPerformanceEventTracker::commitTimeStep(const Schedule& schedule, const int reportStep)
{
    this->lastAccepted_ = this->previous_;
    this->hasLastAccepted_ = true;

    for (const auto& wellName : schedule.wellNames(reportStep)) {
        this->injector_.insert_or_assign(wellName, schedule.getWell(wellName, reportStep).isInjector());
    }

    this->pendingStatusEvents_.clear();
}

void
WellPerformanceEventTracker::applyDeckChanges(WellStatusSnapshot snapshot)
{
    this->accumulate(std::move(snapshot), false);
}

void
WellPerformanceEventTracker::recordWellStatusChangeForRetry(const std::string& wellName,
                                                      const WellStatus status)
{
    if (status == WellStatus::STOP) {
        this->pendingStatusEvents_[wellName].stopped = 1;
    }
    else if (status == WellStatus::SHUT) {
        this->pendingStatusEvents_[wellName].shut = 1;
    }
}

void
WellPerformanceEventTracker::accumulate(WellStatusSnapshot current)
{
    this->accumulate(std::move(current), true);
}

void
WellPerformanceEventTracker::accumulate(WellStatusSnapshot current, const bool trackStatus)
{
    for (const auto& [wellName, now] : current.wells) {
        auto prevPos = this->previous_.wells.find(wellName);
        if (prevPos == this->previous_.wells.end()) {
            // No reference status for this well, e.g., because it appeared
            // mid-run through ACTIONX.  Nothing to compare against.
            continue;
        }

        const auto& before = prevPos->second;

        auto& events = this->events_[wellName];

        if (trackStatus) {
            if ((before.status != WellStatus::STOP) && (now.status == WellStatus::STOP)) {
                events.stopped = 1;
            }

            if ((before.status != WellStatus::SHUT) && (now.status == WellStatus::SHUT)) {
                events.shut = 1;
            }
        }

        // Differenced over connections, not completion numbers.  COMPLUMP
        // may put several connections on one completion and may renumber them
        // while the run is under way; neither opens or closes anything.
        auto opened = std::vector<WellStatusSnapshot::ConnectionID> {};
        std::ranges::set_difference(
            now.openConnections, before.openConnections, std::back_inserter(opened));

        auto closed = std::vector<WellStatusSnapshot::ConnectionID> {};
        std::ranges::set_difference(
            before.openConnections, now.openConnections, std::back_inserter(closed));

        if (opened.empty() && closed.empty()) {
            continue;
        }

        // WPWE1 only counts while the well is able to flow.
        if ((now.status != WellStatus::SHUT) && (now.status != WellStatus::STOP)) {
            events.connsOpened += static_cast<int>(opened.size());
        }

        if (closed.empty()) {
            continue;
        }

        // A '+CON' workover reaches past the connection whose limit was
        // violated and closes the rest of the well below it.  Only the
        // violation itself counts towards WPWE2.
        const auto belowOffender = std::ranges::count_if(
            closed, [&now](const WellStatusSnapshot::ConnectionID& connection) {
                return std::ranges::binary_search(now.closedBelowOffender, connection);
            });

        events.connsClosed += static_cast<int>(closed.size()) - static_cast<int>(belowOffender);

        // WPWE3 has two routes to the bottom of the wellbore.  A workover that
        // closes connections below its offending one gets there directly,
        // regardless of what remains open above; a connection already shut,
        // whether by the deck or by an earlier limit, is not closed again and
        // so does not make a reach.  Otherwise, the closures of this time step
        // must leave either nothing able to flow or only the topmost
        // completion.  A well still flowing farther down, with closed
        // completions above it, has not been closed to the bottom.  The cause
        // of the other closures does not matter: if the deck shut them, a
        // workover that takes the final one still reports WPWE3.
        // openCompletions is a sorted multiset; equal endpoints therefore mean
        // that it contains one distinct completion number.
        if ((belowOffender > 0) ||
            now.openCompletions.empty() ||
            ((now.openCompletions.front() == now.openCompletions.back()) &&
             (now.openCompletions.front() == now.topCompletion)))
        {
            events.closedToBottom = 1;
        }
    }

    this->previous_ = std::move(current);
}

WellPerformanceEventTracker
WellPerformanceEventTracker::serializationTestObject()
{
    auto result = WellPerformanceEventTracker{};

    result.events_["W1"] = data::WellPerformanceEvents::serializationTestObject();
    result.previous_.wells["W1"] =
        WellStatusSnapshot::Entry { WellStatus::OPEN, {{0, 10}, {0, 11}, {0, 12}}, {{1, 13}}, {1, 2, 3}, 1 };
    result.previous_.wells["W2"] =
        WellStatusSnapshot::Entry { WellStatus::SHUT, {}, {}, {}, 1 };
    result.lastAccepted_.wells["W1"] =
        WellStatusSnapshot::Entry { WellStatus::OPEN, {{0, 10}, {0, 11}}, {}, {1, 2}, 1 };
    result.hasLastAccepted_ = true;
    result.injector_["W1"] = false;
    result.injector_["W2"] = true;
    result.pendingStatusEvents_["W2"].shut = 1;

    return result;
}

const data::WellPerformanceEvents&
WellPerformanceEventTracker::events(const std::string& wellName) const
{
    static const auto none = data::WellPerformanceEvents {};

    auto pos = this->events_.find(wellName);

    return (pos == this->events_.end()) ? none : pos->second;
}

} // namespace Opm
