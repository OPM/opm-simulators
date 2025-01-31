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
/*!
 * \file
 * \copydoc Opm::OutputCompositionalModule
 */
#ifndef OPM_OUTPUT_COMPOSITIONAL_MODULE_HPP
#define OPM_OUTPUT_COMPOSITIONAL_MODULE_HPP

#include <dune/grid/common/gridenums.hh>

#include <opm/simulators/utils/moduleVersion.hpp>

#include <opm/common/Exceptions.hpp>
#include <opm/common/TimingMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>

#include <opm/material/common/Valgrind.hpp>
#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/utils/propertysystem.hh>

#include <opm/simulators/flow/FlowBaseVanguard.hpp>
#include <opm/simulators/flow/GenericOutputBlackoilModule.hpp>

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>


namespace Opm
{

// forward declaration
template <class TypeTag>
class EcfvDiscretization;

/*!
 * \ingroup BlackOilSimulator
 *
 * \brief Output module for the results black oil model writing in
 *        ECL binary format.
 */
template <class TypeTag>
class OutputCompositionalModule : public GenericOutputBlackoilModule<GetPropType<TypeTag, Properties::FluidSystem>>
{
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Discretization = GetPropType<TypeTag, Properties::Discretization>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using BaseType = GenericOutputBlackoilModule<FluidSystem>;

    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };

public:
    template <class CollectDataToIORankType>
    OutputCompositionalModule(const Simulator& simulator,
                              const SummaryConfig& smryCfg,
                              const CollectDataToIORankType& collectToIORank)
        : BaseType(simulator.vanguard().eclState(),
                   simulator.vanguard().schedule(),
                   smryCfg,
                   simulator.vanguard().summaryState(),
                   moduleVersionName(),
                   getPropValue<TypeTag, Properties::EnableEnergy>(),
                   getPropValue<TypeTag, Properties::EnableTemperature>(),
                   getPropValue<TypeTag, Properties::EnableMech>(),
                   getPropValue<TypeTag, Properties::EnableSolvent>(),
                   getPropValue<TypeTag, Properties::EnablePolymer>(),
                   getPropValue<TypeTag, Properties::EnableFoam>(),
                   getPropValue<TypeTag, Properties::EnableBrine>(),
                   getPropValue<TypeTag, Properties::EnableSaltPrecipitation>(),
                   getPropValue<TypeTag, Properties::EnableExtbo>(),
                   getPropValue<TypeTag, Properties::EnableMICP>(),
                   true)
        , simulator_(simulator)
    {
        for (auto& region_pair : this->regions_) {
            this->createLocalRegion_(region_pair.second);
        }

        auto isCartIdxOnThisRank = [&collectToIORank](const int idx) {
            return collectToIORank.isCartIdxOnThisRank(idx);
        };

        this->setupBlockData(isCartIdxOnThisRank);

        if (! Parameters::Get<Parameters::OwnerCellsFirst>()) {
            const std::string msg = "The output code does not support --owner-cells-first=false.";
            if (collectToIORank.isIORank()) {
                OpmLog::error(msg);
            }
            OPM_THROW_NOLOG(std::runtime_error, msg);
        }

        if (smryCfg.match("[FB]PP[OGW]") || smryCfg.match("RPP[OGW]*")) {
            auto rset = this->eclState_.fieldProps().fip_regions();
            rset.push_back("PVTNUM");

            // Note: We explicitly use decltype(auto) here because the
            // default scheme (-> auto) will deduce an undesirable type.  We
            // need the "reference to vector" semantics in this instance.
            this->regionAvgDensity_
                .emplace(this->simulator_.gridView().comm(),
                         FluidSystem::numPhases, rset,
                         [fp = std::cref(this->eclState_.fieldProps())]
                         (const std::string& rsetName) -> decltype(auto)
                         { return fp.get().get_int(rsetName); });
        }
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to ECL output files
     */
    void
    allocBuffers(const unsigned bufferSize,
                 const unsigned reportStepNum,
                 const bool     substep,
                 const bool     log,
                 const bool     isRestart)
    {
        if (! std::is_same<Discretization, EcfvDiscretization<TypeTag>>::value) {
            return;
        }

        this->doAllocBuffers(bufferSize,
                              reportStepNum,
                              substep,
                              log,
                              isRestart);
    }

    /*!
     * \brief Modify the internal buffers according to the intensive
     *        quanties relevant for an element
     */
    void processElement(const ElementContext& elemCtx)
    {
        OPM_TIMEBLOCK_LOCAL(processElement);
        if (!std::is_same<Discretization, EcfvDiscretization<TypeTag>>::value)
            return;

        for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
            const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0);
            const auto& fs = intQuants.fluidState();

            const unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
            // const unsigned pvtRegionIdx = 0; // elemCtx.primaryVars(dofIdx, /*timeIdx=*/0).pvtRegionIndex();

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (this->saturation_[phaseIdx].empty())
                    continue;

                this->saturation_[phaseIdx][globalDofIdx] = getValue(fs.saturation(phaseIdx));
                Valgrind::CheckDefined(this->saturation_[phaseIdx][globalDofIdx]);
            }

            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                if (this->moleFractions_[compIdx].empty()) continue;

                this->moleFractions_[compIdx][globalDofIdx] = getValue(fs.moleFraction(compIdx));
            }
            // XMF and YMF
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
                    if (this->phaseMoleFractions_[oilPhaseIdx][compIdx].empty()) continue;
                    this->phaseMoleFractions_[oilPhaseIdx][compIdx][globalDofIdx] = getValue(fs.moleFraction(oilPhaseIdx, compIdx));
                }
                if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
                    if (this->phaseMoleFractions_[gasPhaseIdx][compIdx].empty()) continue;
                    this->phaseMoleFractions_[gasPhaseIdx][compIdx][globalDofIdx] = getValue(fs.moleFraction(gasPhaseIdx, compIdx));
                }
            }

            if (!this->fluidPressure_.empty()) {
                if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
                    // Output oil pressure as default
                    this->fluidPressure_[globalDofIdx] = getValue(fs.pressure(oilPhaseIdx));
                } else if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
                    // Output gas if oil is not present
                    this->fluidPressure_[globalDofIdx] = getValue(fs.pressure(gasPhaseIdx));
                } else {
                    // Output water if neither oil nor gas is present
                    this->fluidPressure_[globalDofIdx] = getValue(fs.pressure(waterPhaseIdx));
                }
                Valgrind::CheckDefined(this->fluidPressure_[globalDofIdx]);
            }

            if (!this->temperature_.empty()) {
                this->temperature_[globalDofIdx] = getValue(fs.temperature(oilPhaseIdx));
                Valgrind::CheckDefined(this->temperature_[globalDofIdx]);
            }
        }
    }

    void processElementFlows(const ElementContext& /* elemCtx */)
    {
        OPM_TIMEBLOCK_LOCAL(processElementBlockData);
        if (!std::is_same_v<Discretization, EcfvDiscretization<TypeTag>>)
            return;
    }

    void processElementBlockData(const ElementContext& /* elemCtx */)
    {
        OPM_TIMEBLOCK_LOCAL(processElementBlockData);
        if (!std::is_same<Discretization, EcfvDiscretization<TypeTag>>::value)
            return;
    }

    /*!
     * \brief Capture connection fluxes, particularly to account for inter-region flows.
     *
     * \tparam ActiveIndex Callable type, typically a lambda, that enables
     *    retrieving the active index, on the local MPI rank, of a
     *    particular cell/element.  Must support a function call operator of
     *    the form
     \code
        int operator()(const Element& elem) const
     \endcode
     *
     * \tparam CartesianIndex Callable type, typically a lambda, that
     *    enables retrieving the globally unique Cartesian index of a
     *    particular cell/element given its active index on the local MPI
     *    rank.  Must support a function call operator of the form
     \code
        int operator()(const int activeIndex) const
     \endcode
     *
     * \param[in] elemCtx Primary lookup structure for per-cell/element
     *    dynamic information.
     *
     * \param[in] activeIndex Mapping from cell/elements to linear indices
     *    on local MPI rank.
     *
     * \param[in] cartesianIndex Mapping from active index on local MPI rank
     *    to globally unique Cartesian cell/element index.
     */
    template <class ActiveIndex, class CartesianIndex>
    void processFluxes(const ElementContext& /* elemCtx */,
                       ActiveIndex&&         /* activeIndex*/,
                       CartesianIndex&&      /* cartesianIndex */)
    {
    }

    /*!
     * \brief Prepare for capturing connection fluxes, particularly to
     *    account for inter-region flows.
     */
    void initializeFluxData()
    {
        // Inter-region flow rates.  Note: ".clear()" prepares to accumulate
        // contributions per bulk connection between FIP regions.
        this->interRegionFlows_.clear();
    }

    /*!
     * \brief Finalize capturing connection fluxes.
     */
    void finalizeFluxData()
    {
        this->interRegionFlows_.compress();
    }

    /*!
     * \brief Get read-only access to collection of inter-region flows.
     */
    const InterRegFlowMap& getInterRegFlows() const
    {
        return this->interRegionFlows_;
    }

    void updateFluidInPlace(const unsigned             /* globalDofIdx */,
                            const IntensiveQuantities& /* intQuants */,
                            const double               /* totVolume */)
    {
        // this->updateFluidInPlace_(globalDofIdx, intQuants, totVolume);
    }

private:
    bool isDefunctParallelWell(std::string wname) const override
    {
        if (simulator_.gridView().comm().size() == 1)
            return false;
        const auto& parallelWells = simulator_.vanguard().parallelWells();
        std::pair<std::string, bool> value {wname, true};
        auto candidate = std::lower_bound(parallelWells.begin(), parallelWells.end(), value);
        return candidate == parallelWells.end() || *candidate != value;
    }

    void createLocalRegion_(std::vector<int>& region)
    {
        std::size_t elemIdx = 0;
        for (const auto& elem : elements(simulator_.gridView())) {
            if (elem.partitionType() != Dune::InteriorEntity) {
                region[elemIdx] = 0;
            }

            ++elemIdx;
        }
    }

    const Simulator& simulator_;
};

} // namespace Opm

#endif // OPM_OUTPUT_COMPOSITIONAL_MODULE_HPP
