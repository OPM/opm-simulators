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
// TODO: remove this include
#include <iostream>

#include <opm/simulators/aquifers/AquiferConstantFlux.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <stdexcept>

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
    for (auto& aquifer : aquifers)
        aquifer->initialSolutionApplied();
}

template <typename TypeTag>
void
BlackoilAquiferModel<TypeTag>::initFromRestart(const data::Aquifers& aquiferSoln)
{
    for (auto& aquifer : this->aquifers)
        aquifer->initFromRestart(aquiferSoln);
}

template <typename TypeTag>
void
BlackoilAquiferModel<TypeTag>::beginEpisode()
{
    // TODO: not totally sure this is the function should be used.
    // basically, we want to update the aquifer related information from SCHEDULE setup in this section
    // it is the beginning of a report step

    const auto& connections = this->simulator_.vanguard().eclState().aquifer().connections();
    const int report_step  = this->simulator_.episodeIndex();
    const auto& aqufluxs = this->simulator_.vanguard().schedule()[report_step].aqufluxs;// .aquflu// simulator.vanguard().schedule()[reportStepIdx].events()
    for (const auto& elem : aqufluxs) {
        const int id = elem.first;
        auto find = std::find_if(begin(this->aquifers), end(this->aquifers), [id](auto& v){ return v->aquiferID() == id; });
        if (find == this->aquifers.end()) {
            // the aquifer id does not exist in aquifers yet
            const auto& aquinfo = elem.second;
            auto aqf = std::make_unique<AquiferConstantFlux<TypeTag>>(aquinfo, connections.getConnections(aquinfo->id), this->simulator_);
            this->aquifers.push_back(std::move(aqf));
        } else {
            const double prev_cumulative_flux = (*find)->cumulativeFlux();
            const auto& aquinfo = elem.second;
            // TODO: it should be improved to be something like a update instead of creating a new one
            auto aqf = std::make_unique<AquiferConstantFlux<TypeTag>>(aquinfo, connections.getConnections(aquinfo->id), this->simulator_, prev_cumulative_flux);
            *find = std::move(aqf);
        }
    }
}

template <typename TypeTag>
void
BlackoilAquiferModel<TypeTag>::beginIteration()
{}

template <typename TypeTag>
void
BlackoilAquiferModel<TypeTag>::beginTimeStep()
{
    for (auto& aquifer : aquifers)
        aquifer->beginTimeStep();
}

template <typename TypeTag>
template <class Context>
void
BlackoilAquiferModel<TypeTag>::addToSource(RateVector& rates,
                                           const Context& context,
                                           unsigned spaceIdx,
                                           unsigned timeIdx) const
{
    for (auto& aquifer : aquifers)
        aquifer->addToSource(rates, context, spaceIdx, timeIdx);
}

template <typename TypeTag>
void
BlackoilAquiferModel<TypeTag>::addToSource(RateVector& rates,
                                           unsigned globalSpaceIdx,
                                           unsigned timeIdx) const
{
    for (auto& aquifer : aquifers)
        aquifer->addToSource(rates, globalSpaceIdx, timeIdx);
}

template <typename TypeTag>
void
BlackoilAquiferModel<TypeTag>::endIteration()
{}

template <typename TypeTag>
void
BlackoilAquiferModel<TypeTag>::endTimeStep()
{
    for (auto& aquifer : aquifers) {
        aquifer->endTimeStep();
        using NumAq = AquiferNumerical<TypeTag>;
        NumAq* num = dynamic_cast<NumAq*>(aquifer.get());
        if (num)
            this->simulator_.vanguard().grid().comm().barrier();
    }
}

template <typename TypeTag>
void
BlackoilAquiferModel<TypeTag>::endEpisode()
{}

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

    const auto& connections = aquifer.connections();
    for (const auto& aq : aquifer.ct()) {
        if (!connections.hasAquiferConnections(aq.aquiferID)) {
            auto msg = fmt::format("No valid connections for Carter-Tracy aquifer {}, aquifer {} will be ignored.",
                                   aq.aquiferID, aq.aquiferID);
            OpmLog::warning(msg);
            continue;
        }
        auto aqf = std::make_unique<AquiferCarterTracy<TypeTag>>(connections.getConnections(aq.aquiferID),
                                                                 this->simulator_, aq);
        aquifers.push_back(std::move(aqf));
    }

    for (const auto& aq : aquifer.fetp()) {
        if (!connections.hasAquiferConnections(aq.aquiferID)) {
            auto msg = fmt::format("No valid connections for Fetkovich aquifer {}, aquifer {} will be ignored.",
                                   aq.aquiferID, aq.aquiferID);
            OpmLog::warning(msg);
            continue;
        }
        auto aqf = std::make_unique<AquiferFetkovich<TypeTag>>(connections.getConnections(aq.aquiferID),
                                                               this->simulator_, aq);
        aquifers.push_back(std::move(aqf));
    }

    if (aquifer.hasNumericalAquifer()) {
        const auto& num_aquifers = aquifer.numericalAquifers().aquifers();
        for ([[maybe_unused]]const auto& [id, aqu] : num_aquifers) {
            auto aqf = std::make_unique<AquiferNumerical<TypeTag>>(aqu, this->simulator_);
            aquifers.push_back(std::move(aqf));
        }
    }

    // first time handle constant flux aquifers, which is stored in the schedule. Other aquifer types might also be refactored later
    // to be able to be updated through SCHEDULE.
}

template<typename TypeTag>
data::Aquifers BlackoilAquiferModel<TypeTag>::aquiferData() const
{
    data::Aquifers data;
    for (const auto& aqu : this->aquifers)
        data.insert_or_assign(aqu->aquiferID(), aqu->aquiferData());

    return data;
}

template<typename TypeTag>
template<class Serializer>
void BlackoilAquiferModel<TypeTag>::
serializeOp(Serializer& serializer)
{
    for (auto& aiPtr : aquifers) {
        auto* ct = dynamic_cast<AquiferCarterTracy<TypeTag>*>(aiPtr.get());
        auto* fetp = dynamic_cast<AquiferFetkovich<TypeTag>*>(aiPtr.get());
        auto* num = dynamic_cast<AquiferNumerical<TypeTag>*>(aiPtr.get());
        if (ct) {
            serializer(*ct);
        } else if (fetp) {
            serializer(*fetp);
        } else if (num) {
            serializer(*num);
        } else {
            OPM_THROW(std::logic_error, "Error serializing BlackoilAquiferModel: unknown aquifer type");
        }
    }
}

} // namespace Opm
