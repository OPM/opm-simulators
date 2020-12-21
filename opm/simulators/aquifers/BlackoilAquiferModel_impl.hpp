/*
  Copyright 2017 TNO - Heat Transfer & Fluid Dynamics, Modelling & Optimization of the Subsurface
  Copyright 2017 Statoil ASA.

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

namespace Opm
{

template <typename TypeTag>
BlackoilAquiferModel<TypeTag>::BlackoilAquiferModel(Simulator& simulator)
    : simulator_(simulator)
{
    // Grid needs to support Facetag
    using Grid = std::remove_const_t<std::remove_reference_t<decltype(simulator.vanguard().grid())>>;
    static_assert(SupportsFaceTag<Grid>::value, "Grid has to support assumptions about face tag.");

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

    // Get all the carter tracy aquifer properties data and put it in aquifers vector
    const auto& connections = aquifer.connections();
    for (const auto& aq : aquifer.ct()) {
        aquifers_CarterTracy.emplace_back(connections[aq.aquiferID],
                                          this->simulator_, aq);
    }

    for (const auto& aq : aquifer.fetp()) {
        aquifers_Fetkovich.emplace_back(connections[aq.aquiferID],
                                        this->simulator_, aq);
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

template<typename TypeTag>
Opm::data::Aquifers BlackoilAquiferModel<TypeTag>::aquiferData() const {
    Opm::data::Aquifers data;
    if (this->aquiferCarterTracyActive()) {
        for (const auto& aqu : aquifers_CarterTracy) {
            Opm::data::AquiferData aqu_data = aqu.aquiferData();
            data[aqu_data.aquiferID] = aqu_data;
        }
    }

    if (this->aquiferFetkovichActive()) {
        for (const auto& aqu : aquifers_Fetkovich) {
            Opm::data::AquiferData aqu_data = aqu.aquiferData();
            data[aqu_data.aquiferID] = aqu_data;
        }
    }
    return data;
}
} // namespace Opm
