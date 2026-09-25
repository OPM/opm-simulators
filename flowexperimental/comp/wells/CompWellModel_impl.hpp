/*
  Copyright 2024 SINTEF Digital

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

#ifndef OPM_COMPOSITIONAL_WELL_MODEL_IMPL_HPP
#define OPM_COMPOSITIONAL_WELL_MODEL_IMPL_HPP

// Improve IDE experience
#ifndef OPM_COMPOSITIONAL_WELL_MODEL_HPP
#include <config.h>
#include <flowexperimental/comp/wells/CompWellModel.hpp>
#endif

#include <set>
#include <algorithm>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <opm/models/utils/parametersystem.hpp>

#include <opm/simulators/flow/BlackoilModelParameters.hpp>

namespace Opm {

template <typename TypeTag>
CompWellModel<TypeTag>::CompWellModel(Simulator& simulator, const NewtonIterationContext& /*iter_ctx*/)
    : WellConnectionModule(*this, simulator.gridView().comm())
    , simulator_(simulator)
    , schedule_(simulator.vanguard().schedule())
    , summary_state_(simulator.vanguard().summaryState())
    , ecl_state_(simulator.vanguard().eclState())
    , comm_(simulator.gridView().comm())
    , comp_config_(ecl_state_.compositionalConfig())
    , comp_well_states_(comp_config_)
    , last_valid_comp_well_states_(comp_config_)
    , dwell_fraction_max_(wellNewtonLimit_<Parameters::DwellFractionMax<Scalar>>())
    , dbhp_max_rel_(wellNewtonLimit_<Parameters::DbhpMaxRel<Scalar>>())
{
    local_num_cells_ = simulator.gridView().size(0);
}

template <typename TypeTag>
template <class Param>
typename CompWellModel<TypeTag>::Scalar
CompWellModel<TypeTag>::
wellNewtonLimit_()
{
    return Parameters::IsSet<Param>(/*errorIfNotRegistered=*/false) ? Parameters::Get<Param>()
                                                                    : Param::value;
}

template <typename TypeTag>
void
CompWellModel<TypeTag>::
beginReportStep(unsigned report_step)
{
    // TODO: not considering the parallel running yet
    report_step_start_events_ = schedule_[report_step].wellgroup_events();
    wells_ecl_ = schedule_.getWells(report_step);

    constexpr auto events_mask = ScheduleEvents::WELL_STATUS_CHANGE |
                                 ScheduleEvents::REQUEST_OPEN_WELL |
                                 ScheduleEvents::REQUEST_SHUT_WELL;
    for (const auto& well_ecl : wells_ecl_) {
        if (!well_ecl.hasConnections()) {
            continue;
        }

        if (!report_step_start_events_.hasEvent(well_ecl.name(), events_mask)) {
            continue;
        }

        if (well_ecl.getStatus() == WellStatus::OPEN) {
            well_open_times_.insert_or_assign(well_ecl.name(), simulator_.time());
            well_close_times_.erase(well_ecl.name());
        }
        else if (well_ecl.getStatus() == WellStatus::SHUT) {
            well_close_times_.insert_or_assign(well_ecl.name(), simulator_.time());
            well_open_times_.erase(well_ecl.name());
        }
    }

    initWellConnectionData();
    initWellState();
    // Save the initial accepted state for this report step.
    last_valid_comp_well_states_.copyDynamicStateFrom(comp_well_states_);
}

template <typename TypeTag>
void
CompWellModel<TypeTag>::
beginTimeStep()
{
    createWellContainer();
    initWellContainer();
}

template <typename TypeTag>
void
CompWellModel<TypeTag>::
restoreLastValidState()
{
    comp_well_states_.copyDynamicStateFrom(last_valid_comp_well_states_);
}

template <typename TypeTag>
void
CompWellModel<TypeTag>::
endTimeStep()
{
    // Persist the accepted well state so failed retries restart from the last
    // successful timestep rather than from the beginning of the report step.
    last_valid_comp_well_states_.copyDynamicStateFrom(comp_well_states_);
}

template <typename TypeTag>
void
CompWellModel<TypeTag>::
init()
{
    simulator_.model().addAuxiliaryModule(this);
}

template <typename TypeTag>
void
CompWellModel<TypeTag>::
createWellContainer()
{
    const auto nw = wells_ecl_.size();
    well_container_.clear();
    for (auto w = 0 * nw; w < nw; ++w) {
        const auto& well_name = wells_ecl_[w].name();
        if (!comp_well_states_.has(well_name)
            || comp_well_states_[well_name].status == WellStatus::SHUT
            || well_connection_data_[w].empty()) {
            continue;
        }

        well_container_.emplace_back(std::make_shared<CompWell<TypeTag>>(wells_ecl_[w], w, well_connection_data_[w],
                                                                          dwell_fraction_max_, dbhp_max_rel_));
    }
}

template <typename TypeTag>
void
CompWellModel<TypeTag>::
initWellContainer()
{
    for (auto& well : well_container_) {
        well->init();
    }
}

template <typename TypeTag>
void
CompWellModel<TypeTag>::
initWellConnectionData()
{
    // Rebuild the local perforation data because schedule events can change
    // connections between report steps.
    well_connection_data_.assign(wells_ecl_.size(), {});
    // Serial runs retain state for every schedule well. In parallel, retain
    // state on ranks that own a connection's grid cell, even if the connection
    // is SHUT and may open at a later report step.
    locally_owned_wells_.assign(wells_ecl_.size(), comm_.size() == 1);
    local_well_reference_cells_.assign(wells_ecl_.size(), -1);

    // Use the partitioner's set of wells that are active during the schedule
    // or may be opened by a schedule action.
    std::set<std::string> everActiveWells;
    for (const auto& well : schedule_.getActiveWellsAtEnd()) {
        everActiveWells.insert(well.name());
    }

    int well_index = 0;
    for (const auto& well : wells_ecl_) {
        int connection_index = 0;
        const auto& well_connections = well.getConnections();
        auto& well_connection_data = well_connection_data_[well_index];

        well_connection_data.reserve(well_connections.size());
        for (const auto& connection : well_connections) {
            const int active_index =
                    this->compressedIndexForInterior(connection.global_index());

            const auto connIsOpen =
                    connection.state() == Connection::State::OPEN;

            if (active_index >= 0) {
                locally_owned_wells_[well_index] = true;
                // Prefer the first local OPEN connection; fall back to the
                // first local connection if all are SHUT.
                if (local_well_reference_cells_[well_index] < 0
                    || (connIsOpen && well_connection_data_[well_index].empty())) {
                    local_well_reference_cells_[well_index] = active_index;
                }
            }

            if (connIsOpen && active_index >= 0) {
                auto& pd = well_connection_data_[well_index].emplace_back();

                pd.cell_index = active_index;
                pd.connection_transmissibility_factor = connection.CF();
                pd.satnum_id = connection.satTableId();
                pd.ecl_index = connection_index;
            }
            ++connection_index;
        }

        // Permanently SHUT wells may span ranks because the default
        // partitioner excludes them from its constraints. They have no well
        // equations to distribute.
        const int ranksWithConnections =
            comm_.sum(locally_owned_wells_[well_index] ? 1 : 0);
        if ((ranksWithConnections > 1) && everActiveWells.count(well.name()) > 0) {
            throw std::runtime_error {
                "Distributed compositional wells are not supported: well '" +
                well.name() + "' has connections on multiple MPI ranks"
            };
        }
        ++well_index;
    }

}

template <typename TypeTag>
void
CompWellModel<TypeTag>::
initWellState()
{
    // TODO: the following might need to be adjusted based on understanding
    const auto pressIx = []()
    {
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) ) {
            return FluidSystem::oilPhaseIdx;
        }
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) ) {
            return FluidSystem::gasPhaseIdx;
        }
        assert(false && "the usage of oil and gas phase is not correct");
        return FluidSystem::gasPhaseIdx;
    }();

    auto cell_pressure = std::vector<Scalar>(this->local_num_cells_, Scalar{0.});
    auto cell_mole_fractions = std::vector<std::vector<Scalar>>(this->local_num_cells_,
                                           std::vector<Scalar>(FluidSystem::numComponents, Scalar{0.}));

    auto cell_temperatures = std::vector<Scalar>(this->local_num_cells_, Scalar{0.});

    auto elemCtx = ElementContext { this->simulator_ };
    const auto& gridView = this->simulator_.vanguard().gridView();

    OPM_BEGIN_PARALLEL_TRY_CATCH();
    for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
        elemCtx.updatePrimaryStencil(elem);
        elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

        const auto ix = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
        const auto& fs = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0).fluidState();

        cell_pressure[ix] = fs.pressure(pressIx).value();
        cell_temperatures[ix] = fs.temperature(0).value();
        for (unsigned compIdx = 0; compIdx < FluidSystem::numComponents; ++compIdx) {
            cell_mole_fractions[ix][compIdx] = fs.moleFraction(compIdx).value();
        }
    }

    OPM_END_PARALLEL_TRY_CATCH("ComposotionalWellModel::initializeWellState() failed: ",
                           this->simulator_.vanguard().grid().comm());

    // Use the reservoir temperature at the first local OPEN connection, or
    // the first local connection if all are SHUT. Wells that may become active
    // must have a single owner, so their reference cell is independent of the
    // partition. Retain the first-cell fallback for serial wells without
    // connections.
    auto well_temperatures = std::vector<Scalar>(
        wells_ecl_.size(), cell_temperatures.empty() ? Scalar{0.} : cell_temperatures.front());
    for (std::size_t wellIdx = 0; wellIdx < wells_ecl_.size(); ++wellIdx) {
        const int cellIdx = local_well_reference_cells_[wellIdx];
        if (cellIdx >= 0) {
            well_temperatures[wellIdx] = cell_temperatures[cellIdx];
        }
    }

    // Carry the wellbore inventory (pressure, water fraction and composition)
    // across report steps, so that schedule-derived targets do not replace it.
    // Runs without water still start each report step from the schedule, as
    // carrying the state over changes their regression results.
    const CompWellState<FluidSystem>* prev_well_state = nullptr;
    if constexpr (FluidSystem::waterEnabled) {
        prev_well_state = &this->last_valid_comp_well_states_;
    }
    this->comp_well_states_.init(this->wells_ecl_,
                                 cell_pressure, well_temperatures, cell_mole_fractions, this->well_connection_data_,
                                 this->summary_state_,
                                 this->locally_owned_wells_,
                                 prev_well_state);
}


template <typename TypeTag>
int
CompWellModel<TypeTag>::
compressedIndexForInterior(std::size_t cartesian_cell_idx) const
{
    return simulator_.vanguard().compressedIndexForInterior(cartesian_cell_idx);
}

template <typename TypeTag>
std::vector<int>
CompWellModel<TypeTag>::
getCellsForConnections(const Well& well) const
{
    std::vector<int> wellCells;
    // All possible connections of the well
    const auto& connectionSet = well.getConnections();
    wellCells.reserve(connectionSet.size());

    for (const auto& connection : connectionSet)
    {
        int compressed_idx = compressedIndexForInterior(connection.global_index());
        if (compressed_idx >= 0) { // Ignore connections in inactive/remote cells.
            wellCells.push_back(compressed_idx);
        }
    }

    return wellCells;

}

template <typename TypeTag>
void
CompWellModel<TypeTag>::
beginIteration()
{
    // do we need to do every iteration here?
    const auto& grid = simulator_.vanguard().grid();
    const auto& gridView = grid.leafGridView();
    ElementContext elemCtx(simulator_);
    for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
        elemCtx.updatePrimaryStencil(elem);
        elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
    }

    assemble(simulator_.timeStepSize());
}

template <typename TypeTag>
void
CompWellModel<TypeTag>::
assemble(const double dt)
{
    const int iterationIdx = simulator_.problem().iterationContext().iteration();

    // Well calculations run only on the owner. Propagate failures to every
    // rank before reservoir assembly so all ranks take the same recovery path.
    OPM_BEGIN_PARALLEL_TRY_CATCH();
    if (iterationIdx == 0) {
        this->calculateExplicitQuantities();
    }
    for (auto& well : well_container_) {
        auto& well_state = comp_well_states_[well->name()];
        well->iterateWellEq(simulator_, dt, well_state);
        // currently we use the converged assembly of well equations directly without a new assembling
        // well->assembleWellEq(simulator_, well_state, dt);
    }
    OPM_END_PARALLEL_TRY_CATCH("CompositionalWellModel::assemble() failed: ", comm_);
}

template <typename TypeTag>
void
CompWellModel<TypeTag>::
calculateExplicitQuantities()
{
    for (auto& well : well_container_) {
        const auto& well_state = comp_well_states_[well->name()];
        well->calculateExplicitQuantities(simulator_, well_state);
    }
}

template <typename TypeTag>
void
CompWellModel<TypeTag>::
computeTotalRatesForDof(RateVector& rate,
                        unsigned globalIdx) const {
    for (const auto& well: well_container_) {
        well->addCellRates(rate, globalIdx);
    }
}

template<typename TypeTag>
void
CompWellModel<TypeTag>::
recoverWellSolutionAndUpdateWellState(const BVector& x)
{
    {
        for (const auto& well : well_container_) {
            const auto& cells = well->cells();
            x_local_.resize(cells.size());

            for (size_t i = 0; i < cells.size(); ++i) {
                x_local_[i] = x[cells[i]];
            }
            auto& ws = this->comp_well_states_[well->name()];
            well->recoverWellSolutionAndUpdateWellState(x_local_, ws);
        }
    }
}

template <typename TypeTag>
bool
CompWellModel<TypeTag>::
forceShutWellByName(const std::string& well_name,
                    double simulation_time,
                    bool)
{
    int well_was_shut = 0;

    if (comp_well_states_.has(well_name)) {
        auto& well_state = comp_well_states_[well_name];
        if (well_state.status != WellStatus::SHUT) {
            well_state.status = WellStatus::SHUT;
            well_close_times_.insert_or_assign(well_name, simulation_time);
            well_open_times_.erase(well_name);

            if (last_valid_comp_well_states_.has(well_name)) {
                last_valid_comp_well_states_[well_name].status = WellStatus::SHUT;
            }

            std::erase_if(well_container_,
                          [&well_name](const auto& well)
                          { return well->name() == well_name; });

            well_was_shut = 1;
        }
    }

    well_was_shut = comm_.max(well_was_shut);
    return well_was_shut == 1;
}

template <typename TypeTag>
bool
CompWellModel<TypeTag>::
getWellConvergence() const
{
    int converged = 1;
    for (const auto& well : this->well_container_) {
        converged = converged && well->getConvergence();
    }
    return comm_.min(converged) == 1;
}

template <typename TypeTag>
data::Wells
CompWellModel<TypeTag>::
wellData() const
{
    return this->comp_well_states_.report();
}

} // end of namespace Opm

#endif
