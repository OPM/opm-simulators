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

#ifndef OPM_ACTION_HANDLER_HPP
#define OPM_ACTION_HANDLER_HPP

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

namespace Opm {

namespace Action {
class ActionX;
class State;
}

class BlackoilWellModelGeneric;
class EclipseState;
class Schedule;
struct SimulatorUpdate;
class SummaryState;
class UDQState;

//! \brief Class handling Action support in simulator
class ActionHandler
{
public:
    //! \brief Function handle to update transmissiblities.
    using TransFunc = std::function<void(bool)>;

    ActionHandler(EclipseState& ecl_state,
                  Schedule& schedule,
                  Action::State& actionState,
                  SummaryState& summaryState,
                  BlackoilWellModelGeneric& wellModel,
                  Parallel::Communication comm);

    void applyActions(int reportStep,
                      double sim_time,
                      const TransFunc& updateTrans);

    //! \brief Evaluates UDQ assign statements.
    void evalUDQAssignments(const unsigned episodeIdx,
                            UDQState& udq_state);

  private:
    /*
       This function is run after applyAction has been completed in the Schedule
       implementation. The sim_update argument should have members & flags for
       the simulator properties which need to be updated. This functionality is
       probably not complete.
    */
    void applySimulatorUpdate(int report_step,
                              const SimulatorUpdate& sim_update,
                              bool& commit_wellstate,
                              const TransFunc& updateTrans);

    std::unordered_map<std::string, double>
    fetchWellPI(int reportStep,
                const Action::ActionX& action,
                const std::vector<std::string>& matching_wells) const;

    EclipseState& ecl_state_;
    Schedule& schedule_;
    Action::State& actionState_;
    SummaryState& summaryState_;
    BlackoilWellModelGeneric& wellModel_;
    Parallel::Communication comm_;
};

} // namespace Opm

#endif // OPM_ACTION_HANDLER_HPP
