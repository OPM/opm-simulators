
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

namespace Opm
{

template <typename TypeTag>
CompositionalWellModel<TypeTag>::
CompositionalWellModel(Simulator& simulator)
    : simulator_(simulator)
    , schedule_(simulator.vanguard().schedule())
    , eclState_(simulator.vanguard().eclState())
    , comm_(simulator.gridView().comm())
      //, schedule_(schedule)
{
}

template <typename TypeTag>
void
CompositionalWellModel<TypeTag>::
beginReportStep(unsigned report_step)
{
    // TODO: not considering the parallel running yet
    wells_ecl_ = schedule_.getWells(report_step);
    initWellConnectionData();
}

template <typename TypeTag>
void CompositionalWellModel<TypeTag>::
beginTimeStep()
{
    createWellContainer();
    initWellContainer();
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
std::size_t CompositionalWellModel<TypeTag>::
compressedIndexForInterior(std::size_t cartesian_cell_idx) const {
    return simulator_.vanguard().compressedIndexForInterior(cartesian_cell_idx);
}

} // end of namespace Opm