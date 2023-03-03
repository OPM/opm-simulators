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

#include <opm/simulators/aquifers/AquiferConstantFlux.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string_view>

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
    this->computeConnectionAreaFraction();

    for (auto& aquifer : this->aquifers) {
        aquifer->initialSolutionApplied();
    }
}

template <typename TypeTag>
void
BlackoilAquiferModel<TypeTag>::initFromRestart(const data::Aquifers& aquiferSoln)
{
    this->computeConnectionAreaFraction();

    for (auto& aquifer : this->aquifers) {
        aquifer->initFromRestart(aquiferSoln);
    }
}

template <typename TypeTag>
void
BlackoilAquiferModel<TypeTag>::beginEpisode()
{
    // Probably function name beginReportStep() is more appropriate.
    //
    // Basically, we want to update the aquifer related information from
    // SCHEDULE setup in this section it is the beginning of a report step

    this->createDynamicAquifers(this->simulator_.episodeIndex());

    this->computeConnectionAreaFraction();
}

template <typename TypeTag>
void
BlackoilAquiferModel<TypeTag>::beginIteration()
{}

template <typename TypeTag>
void
BlackoilAquiferModel<TypeTag>::beginTimeStep()
{
    for (auto& aquifer : this->aquifers) {
        aquifer->beginTimeStep();
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
    for (auto& aquifer : this->aquifers) {
        aquifer->addToSource(rates, context, spaceIdx, timeIdx);
    }
}

template <typename TypeTag>
void
BlackoilAquiferModel<TypeTag>::addToSource(RateVector& rates,
                                           unsigned globalSpaceIdx,
                                           unsigned timeIdx) const
{
    for (auto& aquifer : this->aquifers) {
        aquifer->addToSource(rates, globalSpaceIdx, timeIdx);
    }
}

template <typename TypeTag>
void
BlackoilAquiferModel<TypeTag>::endIteration()
{}

template <typename TypeTag>
void
BlackoilAquiferModel<TypeTag>::endTimeStep()
{
    using NumAq = AquiferNumerical<TypeTag>;

    for (auto& aquifer : this->aquifers) {
        aquifer->endTimeStep();
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
void BlackoilAquiferModel<TypeTag>::init()
{
    if (this->simulator_.vanguard().eclState().aquifer().active()) {
        this->initializeStaticAquifers();
    }

    if (this->needRestartDynamicAquifers()) {
        this->initializeRestartDynamicAquifers();
    }
}

template<typename TypeTag>
data::Aquifers BlackoilAquiferModel<TypeTag>::aquiferData() const
{
    data::Aquifers data;
    for (const auto& aqu : this->aquifers) {
        data.insert_or_assign(aqu->aquiferID(), aqu->aquiferData());
    }

    return data;
}

template<typename TypeTag>
template<class Serializer>
void BlackoilAquiferModel<TypeTag>::
serializeOp(Serializer& serializer)
{
    for (auto& aiPtr : this->aquifers) {
        auto* ct = dynamic_cast<AquiferCarterTracy<TypeTag>*>(aiPtr.get());
        auto* fetp = dynamic_cast<AquiferFetkovich<TypeTag>*>(aiPtr.get());
        auto* num = dynamic_cast<AquiferNumerical<TypeTag>*>(aiPtr.get());
        auto* flux = dynamic_cast<AquiferConstantFlux<TypeTag>*>(aiPtr.get());
        if (ct) {
            serializer(*ct);
        } else if (fetp) {
            serializer(*fetp);
        } else if (num) {
            serializer(*num);
        } else if (flux) {
            serializer(*flux);
        } else {
            OPM_THROW(std::logic_error, "Error serializing BlackoilAquiferModel: unknown aquifer type");
        }
    }
}

template <typename TypeTag>
void BlackoilAquiferModel<TypeTag>::initializeRestartDynamicAquifers()
{
    const auto rstStep = this->simulator_.vanguard().eclState()
        .getInitConfig().getRestartStep() - 1;

    this->createDynamicAquifers(rstStep);
}

template <typename TypeTag>
void BlackoilAquiferModel<TypeTag>::initializeStaticAquifers()
{
    const auto& aquifer =
        this->simulator_.vanguard().eclState().aquifer();

    for (const auto& aquCT : aquifer.ct()) {
        auto aquCTPtr = this->template createAnalyticAquiferPointer
            <AquiferCarterTracy<TypeTag>>(aquCT, aquCT.aquiferID, "Carter-Tracy");

        if (aquCTPtr != nullptr) {
            this->aquifers.push_back(std::move(aquCTPtr));
        }
    }

    for (const auto& aquFetp : aquifer.fetp()) {
        auto aquFetpPtr = this->template createAnalyticAquiferPointer
            <AquiferFetkovich<TypeTag>>(aquFetp, aquFetp.aquiferID, "Fetkovich");

        if (aquFetpPtr != nullptr) {
            this->aquifers.push_back(std::move(aquFetpPtr));
        }
    }

    for (const auto& [id, aquFlux] : aquifer.aquflux()) {
        // Make sure not dummy constant flux aquifers
        if (! aquFlux.active) { continue; }

        auto aquFluxPtr = this->template createAnalyticAquiferPointer
            <AquiferConstantFlux<TypeTag>>(aquFlux, id, "Constant Flux");

        if (aquFluxPtr != nullptr) {
            this->aquifers.push_back(std::move(aquFluxPtr));
        }
    }

    if (aquifer.hasNumericalAquifer()) {
        for (const auto& aquNum : aquifer.numericalAquifers().aquifers()) {
            auto aquNumPtr = std::make_unique<AquiferNumerical<TypeTag>>
                (aquNum.second, this->simulator_);

            this->aquifers.push_back(std::move(aquNumPtr));
        }
    }
}

template <typename TypeTag>
bool BlackoilAquiferModel<TypeTag>::needRestartDynamicAquifers() const
{
    const auto& initconfig =
        this->simulator_.vanguard().eclState().getInitConfig();

    if (! initconfig.restartRequested()) {
        return false;
    }

    return this->simulator_.vanguard()
        .schedule()[initconfig.getRestartStep() - 1].hasAnalyticalAquifers();
}

template <typename TypeTag>
template <typename AquiferType, typename AquiferData>
std::unique_ptr<AquiferType>
BlackoilAquiferModel<TypeTag>::
createAnalyticAquiferPointer(const AquiferData& aqData,
                             const int          aquiferID,
                             std::string_view   aqType) const
{
    const auto& connections =
        this->simulator_.vanguard().eclState().aquifer().connections();

    if (! connections.hasAquiferConnections(aquiferID)) {
        const auto msg = fmt::format("No valid connections for {} aquifer {}.  "
                                     "Aquifer {} will be ignored.",
                                     aqType, aquiferID, aquiferID);
        OpmLog::warning(msg);

        return {};
    }

    return std::make_unique<AquiferType>
        (connections.getConnections(aquiferID), this->simulator_, aqData);
}

template <typename TypeTag>
void BlackoilAquiferModel<TypeTag>::createDynamicAquifers(const int episode_index)
{
    const auto& sched = this->simulator_.vanguard().schedule()[episode_index];

    for (const auto& [id, aquFlux] : sched.aqufluxs) {
        auto aquPos =
            std::find_if(std::begin(this->aquifers),
                         std::end(this->aquifers),
                [id = id](const auto& aquPtr)
            {
                return aquPtr->aquiferID() == id;
            });

        if (aquPos == std::end(this->aquifers)) {
            // An aquifer with this 'id' does not yet exist in
            // the collection managed by this object.  Create it.
            auto aquFluxPtr = this->template createAnalyticAquiferPointer
                <AquiferConstantFlux<TypeTag>>(aquFlux, id, "Constant Flux");

            if (aquFluxPtr != nullptr) {
                this->aquifers.push_back(std::move(aquFluxPtr));
            }
        }
        else {
            auto aquFluxPtr = dynamic_cast<AquiferConstantFlux<TypeTag>*>(aquPos->get());
            if (aquFluxPtr == nullptr) {
                // If the aquifers can return types easily, we might be able
                // to give a better message with type information.
                const auto msg =
                    fmt::format("Aquifer {} is updated with constant flux "
                                "aquifer keyword AQUFLUX at report step {}, "
                                "while it might be specified to be a "
                                "different type of aquifer before this. "
                                "We do not support the conversion between "
                                "different types of aquifer.\n", id, episode_index);

                OPM_THROW(std::runtime_error, msg);
            }

            aquFluxPtr->updateAquifer(aquFlux);
        }
    }
}

template <typename TypeTag>
void BlackoilAquiferModel<TypeTag>::computeConnectionAreaFraction() const
{
    auto maxAquID =
        std::accumulate(this->aquifers.begin(), this->aquifers.end(), 0,
                        [](const int aquID, const auto& aquifer)
                        { return std::max(aquID, aquifer->aquiferID()); });

    maxAquID = this->simulator_.vanguard().grid().comm().max(maxAquID);

    auto totalConnArea = std::vector<double>(maxAquID, 0.0);
    for (const auto& aquifer : this->aquifers) {
        totalConnArea[aquifer->aquiferID() - 1] += aquifer->totalFaceArea();
    }

    this->simulator_.vanguard().grid().comm().sum(totalConnArea.data(), maxAquID);

    for (auto& aquifer : this->aquifers) {
        aquifer->computeFaceAreaFraction(totalConnArea);
    }
}

} // namespace Opm
