#include <opm/grid/utility/cartesianToCompressed.hpp>
namespace Opm
{

template <typename TypeTag>
BlackoilAquiferModel<TypeTag>::BlackoilAquiferModel(Simulator& simulator)
    : simulator_(simulator)
{
    init();
}

template <typename TypeTag>
void
BlackoilAquiferModel<TypeTag>::initialSolutionApplied()
{
    if (aquiferCarterTracyActive()) {
        for (auto& aquifer : aquifers_CarterTracy) {
            aquifer.initialSolutionApplied();
        }
    }
    if (aquiferFetkovichActive()) {
        for (auto& aquifer : aquifers_Fetkovich) {
            aquifer.initialSolutionApplied();
        }
    }
}

template <typename TypeTag>
void
BlackoilAquiferModel<TypeTag>::initFromRestart(const std::vector<data::AquiferData>& aquiferSoln)
{
    if (aquiferCarterTracyActive()) {
        for (auto& aquifer : aquifers_CarterTracy) {
            aquifer.initFromRestart(aquiferSoln);
        }
    }
    if (aquiferFetkovichActive()) {
        for (auto& aquifer : aquifers_Fetkovich) {
            aquifer.initFromRestart(aquiferSoln);
        }
    }
}

template <typename TypeTag>
void
BlackoilAquiferModel<TypeTag>::beginEpisode()
{
}

template <typename TypeTag>
void
BlackoilAquiferModel<TypeTag>::beginIteration()
{
}

template <typename TypeTag>
void
BlackoilAquiferModel<TypeTag>::beginTimeStep()
{
    if (aquiferCarterTracyActive()) {
        for (auto& aquifer : aquifers_CarterTracy) {
            aquifer.beginTimeStep();
        }
    }
    if (aquiferFetkovichActive()) {
        for (auto& aquifer : aquifers_Fetkovich) {
            aquifer.beginTimeStep();
        }
    }
}

template <typename TypeTag>
template <class Context>
void
BlackoilAquiferModel<TypeTag>::addToSource(RateVector& rates,
                                           const Context& context,
                                           unsigned spaceIdx,
                                           unsigned timeIdx) const
{
    if (aquiferCarterTracyActive()) {
        for (auto& aquifer : aquifers_CarterTracy) {
            aquifer.addToSource(rates, context, spaceIdx, timeIdx);
        }
    }
    if (aquiferFetkovichActive()) {
        for (auto& aquifer : aquifers_Fetkovich) {
            aquifer.addToSource(rates, context, spaceIdx, timeIdx);
        }
    }
}

template <typename TypeTag>
void
BlackoilAquiferModel<TypeTag>::endIteration()
{
}

template <typename TypeTag>
void
BlackoilAquiferModel<TypeTag>::endTimeStep()
{
    if (aquiferCarterTracyActive()) {
        for (auto& aquifer : aquifers_CarterTracy) {
            aquifer.endTimeStep();
        }
    }
    if (aquiferFetkovichActive()) {
        for (auto& aquifer : aquifers_Fetkovich) {
            aquifer.endTimeStep();
        }
    }
}
template <typename TypeTag>
void
BlackoilAquiferModel<TypeTag>::endEpisode()
{
}

template <typename TypeTag>
template <class Restarter>
void
BlackoilAquiferModel<TypeTag>::serialize(Restarter& /* res */)
{
    // TODO (?)
    throw std::logic_error("BlackoilAquiferModel::serialize() is not yet implemented");
}

template <typename TypeTag>
template <class Restarter>
void
BlackoilAquiferModel<TypeTag>::deserialize(Restarter& /* res */)
{
    // TODO (?)
    throw std::logic_error("BlackoilAquiferModel::deserialize() is not yet implemented");
}

// Initialize the aquifers in the deck
template <typename TypeTag>
void
BlackoilAquiferModel<TypeTag>::init()
{
    const auto& aquifer = this->simulator_.vanguard().eclState().aquifer();

    if (!aquifer.active()) {
        return;
    }

    const auto& comm = this->simulator_.vanguard().gridView().comm();
    if (comm.size() > 1)
        throw std::runtime_error("Aquifers currently do not work in parallel.");

    // Get all the carter tracy aquifer properties data and put it in aquifers vector
    const auto& ugrid = simulator_.vanguard().grid();
    const int number_of_cells = simulator_.gridView().size(0);

    cartesian_to_compressed_ = cartesianToCompressed(number_of_cells,
                                                     Opm::UgGridHelpers::globalCell(ugrid));

    const auto& connections = aquifer.connections();
    for (const auto& aq : aquifer.ct()) {
        aquifers_CarterTracy.emplace_back(connections[aq.aquiferID],
                                          cartesian_to_compressed_, this->simulator_, aq);
    }

    for (const auto& aq : aquifer.fetp()) {
        aquifers_Fetkovich.emplace_back(connections[aq.aquiferID],
                                        cartesian_to_compressed_, this->simulator_, aq);
    }
}
template <typename TypeTag>
bool
BlackoilAquiferModel<TypeTag>::aquiferActive() const
{
    return (aquiferCarterTracyActive() || aquiferFetkovichActive());
}
template <typename TypeTag>
bool
BlackoilAquiferModel<TypeTag>::aquiferCarterTracyActive() const
{
    return !aquifers_CarterTracy.empty();
}
template <typename TypeTag>
bool
BlackoilAquiferModel<TypeTag>::aquiferFetkovichActive() const
{
    return !aquifers_Fetkovich.empty();
}
} // namespace Opm
