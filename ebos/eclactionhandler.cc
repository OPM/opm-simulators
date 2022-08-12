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

#include <ebos/eclactionhandler.hh>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Action/ActionContext.hpp>
#include <opm/input/eclipse/Schedule/Action/ActionX.hpp>
#include <opm/input/eclipse/Schedule/Action/State.hpp>

#include <opm/simulators/utils/ParallelSerialization.hpp>
#include <opm/simulators/wells/BlackoilWellModelGeneric.hpp>

#include <iomanip>
#include <sstream>

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
{
}

void EclActionHandler::applyActions(int reportStep,
                                    double sim_time,
                                    const TransFunc& transUp)
{
    const auto& actions = schedule_[reportStep].actions();
    if (actions.empty())
        return;

    Action::Context context( summaryState_, schedule_[reportStep].wlist_manager() );
    auto now = TimeStampUTC( schedule_.getStartTime() ) + std::chrono::duration<double>(sim_time);
    std::string ts;
    {
        std::ostringstream os;
        os << std::setw(4) <<                      std::to_string(now.year())  << '/'
           << std::setw(2) << std::setfill('0') << std::to_string(now.month()) << '/'
           << std::setw(2) << std::setfill('0') << std::to_string(now.day()) << "  report:" << std::to_string(reportStep);

        ts = os.str();
    }

    bool commit_wellstate = false;
    for (const auto& pyaction : actions.pending_python(actionState_)) {
        auto sim_update = schedule_.runPyAction(reportStep, *pyaction, actionState_,
                                                ecl_state_, summaryState_);
        this->applySimulatorUpdate(reportStep, sim_update, commit_wellstate, transUp);
    }

    auto simTime = asTimeT(now);
    for (const auto& action : actions.pending(actionState_, simTime)) {
        auto actionResult = action->eval(context);
        if (actionResult) {
            std::string wells_string;
            const auto& matching_wells = actionResult.wells();
            if (!matching_wells.empty()) {
                for (std::size_t iw = 0; iw < matching_wells.size() - 1; iw++)
                    wells_string += matching_wells[iw] + ", ";
                wells_string += matching_wells.back();
            }
            std::string msg = "The action: " + action->name() + " evaluated to true at " + ts + " wells: " + wells_string;
            OpmLog::info(msg);

            const auto& wellpi = this->fetchWellPI(reportStep, *action, matching_wells);

            auto sim_update = schedule_.applyAction(reportStep, *action,
                                                    actionResult.wells(), wellpi);
            this->applySimulatorUpdate(reportStep, sim_update,  commit_wellstate, transUp);
            actionState_.add_run(*action, simTime, std::move(actionResult));
        } else {
            std::string msg = "The action: " + action->name() + " evaluated to false at " + ts;
            OpmLog::info(msg);
        }
    }
    /*
      The well state has been stored in a previous object when the time step
      has completed successfully, the action process might have modified the
      well state, and to be certain that is not overwritten when starting
      the next timestep we must commit it.
    */
    if (commit_wellstate)
        this->wellModel_.commitWGState();
}

void EclActionHandler::applySimulatorUpdate(int report_step,
                                            const SimulatorUpdate& sim_update,
                                            bool& commit_wellstate,
                                            const TransFunc& updateTrans)
  {
      this->wellModel_.updateEclWells(report_step, sim_update.affected_wells, summaryState_);
      if (!sim_update.affected_wells.empty())
          commit_wellstate = true;

      if (sim_update.tran_update) {
          const auto& keywords = schedule_[report_step].geo_keywords();
          ecl_state_.apply_schedule_keywords( keywords );
          eclBroadcast(comm_, ecl_state_.getTransMult() );

          // re-compute transmissibility
          updateTrans(true);
      }
  }

std::unordered_map<std::string, double>
EclActionHandler::fetchWellPI(int reportStep,
                              const Action::ActionX& action,
                              const std::vector<std::string>& matching_wells)
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

} // namespace Opm
