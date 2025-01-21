
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <opm/simulators/utils/phaseUsageFromDeck.hpp>

namespace Opm
{

template <typename TypeTag>
CompositionalWellModel<TypeTag>::
CompositionalWellModel(Simulator& simulator)
    : simulator_(simulator)
    , schedule_(simulator.vanguard().schedule())
    , summary_state_(simulator.vanguard().summaryState())
    , ecl_state_(simulator.vanguard().eclState())
    , comm_(simulator.gridView().comm())
    , phase_usage_(phaseUsageFromDeck(ecl_state_))
    , comp_config_(ecl_state_.compositionalConfig())
    , comp_well_states_(phase_usage_, comp_config_)
{
    local_num_cells_ = simulator.gridView().size(0);
}

template <typename TypeTag>
void
CompositionalWellModel<TypeTag>::
beginReportStep(unsigned report_step)
{
    // TODO: not considering the parallel running yet
    wells_ecl_ = schedule_.getWells(report_step);
    initWellConnectionData();
    initWellState(report_step);
}

template <typename TypeTag>
void CompositionalWellModel<TypeTag>::
beginTimeStep()
{
    createWellContainer();
    initWellContainer();
    for (auto& well : well_container_) {
        well->calculateExplicitQuantities(simulator_, comp_well_states_[well->name()]);
    }
}

template <typename TypeTag>
void CompositionalWellModel<TypeTag>::
createWellContainer()
{
    // const auto& schedule = simulator_.vanguard().schedule();
    const auto nw = wells_ecl_.size(); // not considering the parallel running yet
    well_container_.clear();
    for (auto w = 0 * nw; w < nw; ++w) {
        well_container_.emplace_back(std::make_shared<CompWell<TypeTag>>(wells_ecl_[w], w, well_connection_data_[w]));
    }
}

template <typename TypeTag>
void CompositionalWellModel<TypeTag>::
initWellContainer()
{
    for (auto& well : well_container_) {
        well->init();
    }
}

template <typename TypeTag>
void CompositionalWellModel<TypeTag>::
initWellConnectionData()
{
    // TODO: we need to consider the parallel running
    // we can refer to the BlackoilWellModelGeneric::initializeWellPerfData()
    well_connection_data_.resize(wells_ecl_.size());

    int well_index = 0;
    for (const auto& well : wells_ecl_) {
        int connection_index = 0;
        const auto& well_connections = well.getConnections();
        auto& well_connection_data = well_connection_data_[well_index];

        well_connection_data.reserve(well_connections.size());
        for (const auto& connection : well_connections) {
            const auto active_index =
                    this->compressedIndexForInterior(connection.global_index());

            const auto connIsOpen =
                    connection.state() == Connection::State::OPEN;

            if (connIsOpen && (active_index >= 0)) {
                auto& pd = well_connection_data_[well_index].emplace_back();

                pd.cell_index = active_index;
                pd.connection_transmissibility_factor = connection.CF();
                pd.satnum_id = connection.satTableId();
                pd.ecl_index = connection_index;
            }
            ++connection_index;
        }
        ++well_index;
    }

}

template <typename TypeTag>
void CompositionalWellModel<TypeTag>::
initWellState(const std::size_t report_step)
{
    // TODO: not sure the following is correct
    const auto pressIx = [this]()
    {
        if (phase_usage_.phase_used[FluidSystem::oilPhaseIdx]) {
            return FluidSystem::oilPhaseIdx;
        }
        if (phase_usage_.phase_used[FluidSystem::gasPhaseIdx]) {
            return FluidSystem::gasPhaseIdx;
        }
        assert(false && "the usage of oil and gas phase is not correct");
        return FluidSystem::gasPhaseIdx;
    }();

    auto cell_pressure = std::vector<Scalar>(this->local_num_cells_, Scalar{0.});

    auto elemCtx = ElementContext { this->simulator_ };
    const auto& gridView = this->simulator_.vanguard().gridView();

    OPM_BEGIN_PARALLEL_TRY_CATCH();
    for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
        elemCtx.updatePrimaryStencil(elem);
        elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

        const auto ix = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
        const auto& fs = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0).fluidState();

        cell_pressure[ix] = fs.pressure(pressIx).value();
    }
    OPM_END_PARALLEL_TRY_CATCH("ComposotionalWellModel::initializeWellState() failed: ",
                           this->simulator_.vanguard().grid().comm());

    /* TODO: no prev well state for now */
    this->comp_well_states_.init(this->schedule_, this->wells_ecl_,
                                 cell_pressure, this->well_connection_data_,
                                 this->summary_state_);
}


template <typename TypeTag>
std::size_t CompositionalWellModel<TypeTag>::
compressedIndexForInterior(std::size_t cartesian_cell_idx) const {
    return simulator_.vanguard().compressedIndexForInterior(cartesian_cell_idx);
}

} // end of namespace Opm