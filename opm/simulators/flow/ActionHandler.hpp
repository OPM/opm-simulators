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

namespace Opm::Action {
    class State;
} // namespace Opm::Action

namespace Opm {
template<class Scalar> class BlackoilWellModelGeneric;
class EclipseState;
class Schedule;
struct SimulatorUpdate;
class SummaryState;
class UDQState;
} // namespace Opm

namespace Opm {

//! \brief Class handling Action support in simulator
template<class Scalar>
class ActionHandler
{
public:
    //! \brief Function handle to update transmissiblities.
    using TransFunc = std::function<void(bool)>;

    /// Constructor.
    ///
    /// \param[in,out] ecl_state Container of static properties such as
    /// permeability and transmissibility.
    ///
    /// \param[in,out] schedule Container of dynamic objects, such as wells.
    ///
    /// \param[in,out] actionState Dynamic state object for all actions.
    ///
    /// \param[in,out] summaryState Dynamic state object for all summary
    /// vectors.
    ///
    /// \param[in,out] wellModel Simulation wells on this rank.
    ///
    /// \param[in] comm MPI communicator object linking all simulation
    /// ranks.
    ActionHandler(EclipseState& ecl_state,
                  Schedule& schedule,
                  Action::State& actionState,
                  SummaryState& summaryState,
                  BlackoilWellModelGeneric<Scalar>& wellModel,
                  Parallel::Communication comm);

    /// Run all pending actions.
    ///
    /// \param[in] reportStep Zero-based report step index.
    ///
    /// \param[in] sim_time Elapsed time since simulation start.
    ///
    /// \param[in] updateTrans Call-back for affecting transmissibility
    /// updates.  Typically invoked if the action triggers a keyword like
    /// MULTZ.
    void applyActions(int reportStep,
                      double sim_time,
                      const TransFunc& updateTrans);

    /// \brief Evaluates UDQ assign statements.
    ///
    /// \param[in] episodeIdx Zero-based report step index.
    ///
    /// \param[in,out] udq_state Dynamic state of all user-defined
    /// quantities.
    void evalUDQAssignments(const unsigned episodeIdx,
                            UDQState& udq_state);

    /// Convey dynamic updates triggered by an action block back to the
    /// running simulator.
    ///
    /// This function is run after applyAction has been completed in the
    /// Schedule implementation.  The sim_update argument should have
    /// members & flags for the simulator properties which need to be
    /// updated.  This functionality is probably not complete.
    ///
    /// \param[in] report_step Zero-based report step index.
    ///
    /// \param[in] sim_update Action's resulting simulator update.
    ///
    /// \param[in] updateTrans Call-back for affecting transmissibility
    /// updates.  Typically invoked if the action triggers a keyword like
    /// MULTZ.
    ///
    /// \param[out] commit_wellstate Whether or not the action affected any
    /// simulation wells which, in turn, may require rebuilding internal
    /// data structures in the simulator and therefore would require
    /// preserving the dynamic well and group states prior to doing so.
    void applySimulatorUpdate(int report_step,
                              const SimulatorUpdate& sim_update,
                              const TransFunc& updateTrans,
                              bool& commit_wellstate);

private:
    /// Static properties such as permeability and transmissibility.
    EclipseState& ecl_state_;

    /// Dynamic objects such as wells.
    Schedule& schedule_;

    /// Dynamic state for all actions--e.g., their run count.
    Action::State& actionState_;

    /// Dynamic state for all user-defined quantities.
    SummaryState& summaryState_;

    /// Simulation wells on this rank.
    BlackoilWellModelGeneric<Scalar>& wellModel_;

    /// MPI communicator object linking all simulation ranks.
    Parallel::Communication comm_;
};

} // namespace Opm

#endif // OPM_ACTION_HANDLER_HPP
