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

#include <opm/simulators/wells/WellPerformanceEvents.hpp>

#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <algorithm>
#include <iterator>
#include <utility>

namespace Opm
{

void
WellPerformanceEvents::beginTimeStep(const Schedule& schedule,
                                     const int reportStep,
                                     WellStatusSnapshot snapshot)
{
    this->events_.clear();
    this->previous_ = std::move(snapshot);

    for (const auto& wellName : schedule.wellNames(reportStep)) {
        const auto isInjector = schedule.getWell(wellName, reportStep).isInjector();

        const auto [pos, inserted] = this->injector_.try_emplace(wellName, isInjector);
        if (inserted || (pos->second == isInjector)) {
            continue;
        }

        auto& events = this->events_[wellName];
        (isInjector ? events.producerToInjector : events.injectorToProducer) = 1;
    }
}

void
WellPerformanceEvents::commitTimeStep(const Schedule& schedule, const int reportStep)
{
    for (const auto& wellName : schedule.wellNames(reportStep)) {
        this->injector_.insert_or_assign(wellName, schedule.getWell(wellName, reportStep).isInjector());
    }
}

void
WellPerformanceEvents::synchronise(WellStatusSnapshot snapshot)
{
    this->previous_ = std::move(snapshot);
}

void
WellPerformanceEvents::accumulate(WellStatusSnapshot current)
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

        if ((before.status != WellStatus::STOP) && (now.status == WellStatus::STOP)) {
            events.stopped = 1;
        }

        if ((before.status != WellStatus::SHUT) && (now.status == WellStatus::SHUT)) {
            events.shut = 1;
        }

        auto opened = std::vector<int> {};
        std::set_difference(now.openCompletions.begin(),
                            now.openCompletions.end(),
                            before.openCompletions.begin(),
                            before.openCompletions.end(),
                            std::back_inserter(opened));

        auto closed = std::vector<int> {};
        std::set_difference(before.openCompletions.begin(),
                            before.openCompletions.end(),
                            now.openCompletions.begin(),
                            now.openCompletions.end(),
                            std::back_inserter(closed));

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

        // Several CON and +CON workovers can occur between snapshots.  Only
        // exclude connections actually closed by +CON from the WPWE2 count.
        const auto numConPlus = std::count_if(closed.begin(), closed.end(),
            [&now](const int complnum)
            {
                return std::binary_search(now.closedByConPlus.begin(),
                                          now.closedByConPlus.end(), complnum);
            });

        events.connsClosed += static_cast<int>(closed.size()) - static_cast<int>(numConPlus);

        // WPWE3 has two routes: a '+CON' workover, which closes to the bottom
        // by construction, or a well left with no flowing connection because
        // 'CON' workovers closed all of them.
        if ((numConPlus > 0) || now.closedToBottomByCon) {
            events.closedToBottom = 1;
        }
    }

    this->previous_ = std::move(current);
}

WellPerformanceEvents
WellPerformanceEvents::serializationTestObject()
{
    auto result = WellPerformanceEvents{};

    result.events_["W1"] = data::WellEvents::serializationTestObject();
    result.previous_.wells["W1"] = WellStatusSnapshot::Entry { WellStatus::OPEN, {1, 2, 3}, {4}, false };
    result.previous_.wells["W2"] = WellStatusSnapshot::Entry { WellStatus::SHUT, {}, {}, true };
    result.injector_["W1"] = false;
    result.injector_["W2"] = true;

    return result;
}

const data::WellEvents&
WellPerformanceEvents::events(const std::string& wellName) const
{
    static const auto none = data::WellEvents {};

    auto pos = this->events_.find(wellName);

    return (pos == this->events_.end()) ? none : pos->second;
}

} // namespace Opm
