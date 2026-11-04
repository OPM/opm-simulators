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
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace Opm::Action {
    class State;
} // namespace Opm::Action

namespace Opm {
template<typename Scalar, typename IndexTraits> class BlackoilWellModelGeneric;
class EclipseState;
class Schedule;
struct SimulatorUpdate;
class SummaryState;
class UDQState;
} // namespace Opm

namespace Opm {

//! \brief Class handling Action support in simulator
template<typename Scalar, typename IndexTraits>
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
                  BlackoilWellModelGeneric<Scalar, IndexTraits>& wellModel,
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

    //! \brief In a sector run (USEFLUX), trigger the ACTIONX whose
    //!        conditions read quantities the sector does not compute itself
    //!        when the parent run triggered them, instead of evaluating them.
    //!
    //! \details An action is evaluated in this run when every quantity its
    //!   conditions read is one the sector computes: a constant or the date, a
    //!   well, connection or segment quantity of wells inside the sector, a
    //!   group quantity of a group whose wells are all inside it, a block
    //!   quantity of one of its cells, or a UDQ defined from these only.
    //!   Everything else (field and region quantities, wells and cells outside
    //!   the sector, and UDQs of these) is the parent's at the parent's
    //!   sampling, and a condition on it, in particular one the parent met
    //!   between two of its samples or only just, may come out otherwise. Such
    //!   actions are triggered at the end of the time step in which the
    //!   parent triggered them, with the wells it matched. Writes to the PRT
    //!   which actions are replayed and which are not, and why.
    //!
    //! \param[in] triggerNames, triggerTimes, triggerWellCounts, triggerWells
    //!   The parent's triggers, from the FLUX file (FluxFile::ParentActions):
    //!   times in seconds since the start of the run.
    //! \param[in] replay Whether replay is enabled (--flux-replay-actions).
    void setupFluxActionReplay(const std::vector<std::string>& triggerNames,
                               const std::vector<double>& triggerTimes,
                               const std::vector<int>& triggerWellCounts,
                               const std::vector<std::string>& triggerWells,
                               bool replay);

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
    BlackoilWellModelGeneric<Scalar, IndexTraits>& wellModel_;

    /// MPI communicator object linking all simulation ranks.
    Parallel::Communication comm_;

    struct ReplayTrigger
    {
        double time;
        std::vector<std::string> wells;
    };

    /// The parent's triggers of each replayed action, by time.
    std::unordered_map<std::string, std::vector<ReplayTrigger>> replayTriggers_;

    /// Index of the first trigger of each replayed action not replayed yet.
    std::unordered_map<std::string, std::size_t> replayNext_;

    /// The actions that are replayed rather than evaluated.
    std::unordered_set<std::string> replayed_;

    /// Why the sector does not compute a quantity, or nothing when it does.
    std::optional<std::string> parentQuantityReason_(const std::string& quantity,
                                                     const std::vector<std::string>& args,
                                                     const std::unordered_set<std::string>& sectorWells,
                                                     std::unordered_map<std::string, std::optional<std::string>>& udqMemo) const;
};

} // namespace Opm

#endif // OPM_ACTION_HANDLER_HPP
