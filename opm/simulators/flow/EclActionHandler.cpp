// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/

#include <config.h>
#include <opm/simulators/flow/EclActionHandler.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/utility/TimeService.hpp>
#include <opm/common/TimingMacros.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Schedule/Action/ActionContext.hpp>
#include <opm/input/eclipse/Schedule/Action/Actions.hpp>
#include <opm/input/eclipse/Schedule/Action/ActionX.hpp>
#include <opm/input/eclipse/Schedule/Action/SimulatorUpdate.hpp>
#include <opm/input/eclipse/Schedule/Action/State.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQConfig.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellMatcher.hpp>

#include <opm/simulators/wells/BlackoilWellModelGeneric.hpp>
#include <opm/simulators/utils/ParallelSerialization.hpp>

#include <chrono>
#include <cstddef>
#include <ctime>
#include <string>
#include <unordered_map>
#include <vector>

#include <fmt/chrono.h>
#include <fmt/format.h>

namespace {
    std::string formatActionDate(const Opm::TimeStampUTC& timePoint,
                                 const int                reportStep)
    {
        auto time_point = std::tm{};

        time_point.tm_year = timePoint.year()  - 1900;
        time_point.tm_mon  = timePoint.month() -    1;
        time_point.tm_mday = timePoint.day();

        time_point.tm_hour = timePoint.hour();
        time_point.tm_min  = timePoint.minutes();
        time_point.tm_sec  = timePoint.seconds();

        return fmt::format("{:%d-%b-%Y %H:%M:%S} (report interval {} to {})",
                           time_point, reportStep, reportStep + 1);
    }

    void logActiveAction(const std::string&              actionName,
                         const std::vector<std::string>& matchingWells,
                         const std::string&              timeString)
    {
        const auto wellString = matchingWells.empty()
            ? std::string{}
            : fmt::format(" Well{}: {}",
                          matchingWells.size() != 1 ? "s" : "",
                          fmt::join(matchingWells, ", "));

        const auto message =
            fmt::format("Action {} triggered at {}.{}",
                        actionName, timeString, wellString);

        Opm::OpmLog::info("ACTION_TRIGGERED", message);
    }

    void logInactiveAction(const std::string& actionName,
                           const std::string& timeString)
    {
        const auto message =
            fmt::format("Action {} NOT triggered at {}.",
                        actionName, timeString);

        Opm::OpmLog::debug("NAMED_ACTION_NOT_TRIGGERED", message);
    }

    void logInactiveActions(const int          numInactive,
                            const std::string& timeString)
    {
        const auto message =
            fmt::format("{} action{} NOT triggered at {}.",
                        numInactive,
                        (numInactive != 1) ? "s" : "",
                        timeString);

        Opm::OpmLog::debug("ACTION_NOT_TRIGGERED", message);
    }
} // Anonymous namespace

namespace Opm {

EclActionHandler::EclActionHandler(EclipseState& ecl_state,
                                   Schedule& schedule,
                                   Action::State& actionState,
                                   SummaryState& summaryState,
                                   BlackoilWellModelGeneric& wellModel,
                                   Parallel::Communication comm)
    : ecl_state_(ecl_state)
    , schedule_(schedule)
    , actionState_(actionState)
    , summaryState_(summaryState)
    , wellModel_(wellModel)
    , comm_(comm)
{}

void EclActionHandler::applyActions(const int reportStep,
                                    const double sim_time,
                                    const TransFunc& transUp)
{
    OPM_TIMEBLOCK(applyActions);
    const auto& actions = schedule_[reportStep].actions();
    if (actions.empty()) {
        return;
    }

    const Action::Context context{ summaryState_, schedule_[reportStep].wlist_manager() };

    const auto now = TimeStampUTC{ schedule_.getStartTime() } + std::chrono::duration<double>(sim_time);
    const auto ts  = formatActionDate(now, reportStep);

    bool commit_wellstate = false;
    for (const auto& pyaction : actions.pending_python(actionState_)) {
        auto sim_update = schedule_.runPyAction(reportStep, *pyaction, actionState_,
                                                ecl_state_, summaryState_);
        this->applySimulatorUpdate(reportStep, sim_update, commit_wellstate, transUp);
    }

    auto non_triggered = 0;
    const auto simTime = asTimeT(now);
    for (const auto& action : actions.pending(actionState_, simTime)) {
        const auto actionResult = action->eval(context);
        if (! actionResult) {
            ++non_triggered;
            logInactiveAction(action->name(), ts);
            continue;
        }

        const auto& matching_wells = actionResult.wells();

        logActiveAction(action->name(), matching_wells, ts);

        const auto wellpi = this->fetchWellPI(reportStep, *action, matching_wells);

        const auto sim_update = this->schedule_
            .applyAction(reportStep, *action, matching_wells, wellpi);

        this->applySimulatorUpdate(reportStep, sim_update, commit_wellstate, transUp);
        actionState_.add_run(*action, simTime, std::move(actionResult));
    }

    if (non_triggered > 0) {
        logInactiveActions(non_triggered, ts);
    }

    // The well state has been stored in a previous object when the time
    // step has completed successfully, the action process might have
    // modified the well state, and to be certain that is not overwritten
    // when starting the next timestep we must commit it.
    if (commit_wellstate) {
        this->wellModel_.commitWGState();
    }
}

void EclActionHandler::applySimulatorUpdate(const int report_step,
                                            const SimulatorUpdate& sim_update,
                                            bool& commit_wellstate,
                                            const TransFunc& updateTrans)
{
    OPM_TIMEBLOCK(applySimulatorUpdate);

    this->wellModel_.updateEclWells(report_step, sim_update, this->summaryState_);

    if (!sim_update.affected_wells.empty()) {
        commit_wellstate = true;
    }

    if (sim_update.tran_update) {
        const auto& keywords = schedule_[report_step].geo_keywords();
        ecl_state_.apply_schedule_keywords( keywords );
        eclBroadcast(comm_, ecl_state_.getTransMult() );

        // re-compute transmissibility
        updateTrans(true);
    }
}

std::unordered_map<std::string, double>
EclActionHandler::fetchWellPI(const int reportStep,
                              const Action::ActionX& action,
                              const std::vector<std::string>& matching_wells) const
{

  auto wellpi_wells = action.wellpi_wells(WellMatcher(schedule_[reportStep].well_order(),
                                                      schedule_[reportStep].wlist_manager()),
                                          matching_wells);

  if (wellpi_wells.empty())
      return {};

  const auto num_wells = schedule_[reportStep].well_order().size();
  std::vector<double> wellpi_vector(num_wells);
  for (const auto& wname : wellpi_wells) {
      if (this->wellModel_.hasWell(wname)) {
          const auto& well = schedule_.getWell( wname, reportStep );
          wellpi_vector[well.seqIndex()] = this->wellModel_.wellPI(wname);
      }
  }

  if (comm_.size() > 1) {
      std::vector<double> wellpi_buffer(num_wells * comm_.size());
      comm_.gather( wellpi_vector.data(), wellpi_buffer.data(), num_wells, 0 );
      if (comm_.rank() == 0) {
          for (int rank=1; rank < comm_.size(); rank++) {
              for (std::size_t well_index=0; well_index < num_wells; well_index++) {
                  const auto global_index = rank*num_wells + well_index;
                  const auto value = wellpi_buffer[global_index];
                  if (value != 0)
                      wellpi_vector[well_index] = value;
              }
          }
      }
      comm_.broadcast(wellpi_vector.data(), wellpi_vector.size(), 0);
  }

  std::unordered_map<std::string, double> wellpi;
  for (const auto& wname : wellpi_wells) {
      const auto& well = schedule_.getWell( wname, reportStep );
      wellpi[wname] = wellpi_vector[ well.seqIndex() ];
  }
  return wellpi;
}

void EclActionHandler::evalUDQAssignments(const unsigned episodeIdx,
                                          UDQState& udq_state)
{
    const auto& udq = schedule_[episodeIdx].udq();

    udq.eval_assign(episodeIdx,
                    this->schedule_,
                    this->schedule_.wellMatcher(episodeIdx),
                    this->schedule_.segmentMatcherFactory(episodeIdx),
                    this->summaryState_,
                    udq_state);
}

} // namespace Opm
