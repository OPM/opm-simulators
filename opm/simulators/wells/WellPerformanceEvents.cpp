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
#include <opm/input/eclipse/Schedule/Well/Connection.hpp>
#include <opm/input/eclipse/Schedule/Well/ConnectionEconLimits.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#include <opm/input/eclipse/Schedule/Well/WellEconProductionLimits.hpp>

#include <algorithm>
#include <iterator>
#include <utility>

namespace
{

/// Whether the workover procedure that closes \p complnum of \p well also
/// closes every connection below it in the wellbore, i.e., whether it is
/// a '+CON' procedure.
///
/// A connection carries its own CECON limits when the deck defines them.
/// Otherwise the well level WECON procedure applies.
bool
closesToBottom(const Opm::Well& well, const int complnum)
{
    for (const auto& connection : well.getConnections()) {
        if (connection.complnum() != complnum) {
            continue;
        }

        if (connection.hasEconLimits()) {
            return connection.econLimits().workover
                == Opm::ConnectionEconLimits::EconWorkover::CONP;
        }

        break;
    }

    return well.getEconLimits().workover() == Opm::WellEconProductionLimits::EconWorkover::CONP;
}

} // Anonymous namespace

namespace Opm
{

void
WellPerformanceEvents::beginTimeStep(const Schedule& schedule,
                                     const int reportStep,
                                     const bool reportStepStarts,
                                     WellStatusSnapshot snapshot)
{
    this->events_.clear();
    this->previous_ = std::move(snapshot);

    if (!reportStepStarts) {
        return;
    }

    for (const auto& wellName : schedule.wellNames(reportStep)) {
        const auto isInjector = schedule.getWell(wellName, reportStep).isInjector();

        const auto [pos, inserted] = this->injector_.try_emplace(wellName, isInjector);
        if (inserted || (pos->second == isInjector)) {
            continue;
        }

        auto& events = this->events_[wellName];
        (isInjector ? events.producerToInjector : events.injectorToProducer) = 1;

        pos->second = isInjector;
    }
}

void
WellPerformanceEvents::synchronise(WellStatusSnapshot snapshot)
{
    this->previous_ = std::move(snapshot);
}

void
WellPerformanceEvents::accumulate(const Schedule& schedule,
                                  const int reportStep,
                                  WellStatusSnapshot current)
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

        const auto& well = schedule.getWell(wellName, reportStep);

        // A '+CON' workover closes the offending connection and every one
        // below it in a single event, and only the connection that violated
        // its limit carries the procedure.  Recognising the procedure on any
        // of the closed connections therefore attributes the whole event.
        // Such an event is reported through WPWE3 alone.
        const auto toBottom
            = std::any_of(closed.begin(), closed.end(), [&well](const int complnum) {
                  return closesToBottom(well, complnum);
              });

        if (toBottom) {
            events.closedToBottom = 1;
            continue;
        }

        events.connsClosed += static_cast<int>(closed.size());

        // Closing the last connections that were able to flow leaves the well
        // shut off at the bottom just as a '+CON' workover does.
        if (now.openCompletions.empty()) {
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
    result.previous_.wells["W1"] = WellStatusSnapshot::Entry { WellStatus::OPEN, {1, 2, 3} };
    result.previous_.wells["W2"] = WellStatusSnapshot::Entry { WellStatus::SHUT, {} };
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
