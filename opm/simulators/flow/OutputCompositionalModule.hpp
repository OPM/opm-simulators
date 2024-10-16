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
 * \copydoc Opm::OutputBlackOilModule
 */
#ifndef OPM_OUTPUT_COMPOSITIONAL_MODULE_HPP
#define OPM_OUTPUT_COMPOSITIONAL_MODULE_HPP

#include <dune/common/fvector.hh>

#include <opm/simulators/utils/moduleVersion.hpp>

#include <opm/common/Exceptions.hpp>
#include <opm/common/TimingMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/fluidmatrixinteractions/EclEpsScalingPoints.hpp>
// #include <opm/material/fluidstates/BlackOilFluidState.hpp>
// #include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/models/blackoil/blackoilproperties.hh>
#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/utils/propertysystem.hh>

#include <opm/output/data/Cells.hpp>
#include <opm/output/eclipse/EclipseIO.hpp>
#include <opm/output/eclipse/Inplace.hpp>

#include <opm/simulators/flow/FlowBaseVanguard.hpp>
// #include <opm/simulators/flow/GenericOutputCompositionalModule.hpp>
#include <opm/simulators/flow/GenericOutputBlackoilModule.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace Opm::Parameters {

struct ForceDisableFluidInPlaceOutput { static constexpr bool value = false; };
struct ForceDisableResvFluidInPlaceOutput { static constexpr bool value = false; };

} // namespace Opm::Parameters

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
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementIterator = typename GridView::template Codim<0>::Iterator;
    using BaseType = GenericOutputBlackoilModule<FluidSystem>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using Dir = FaceDir::DirEnum;

    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { gasCompIdx = FluidSystem::gasCompIdx };
    enum { oilCompIdx = FluidSystem::oilCompIdx };
    enum { waterCompIdx = FluidSystem::waterCompIdx };
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };

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
                   getPropValue<TypeTag, Properties::EnableMICP>())
        , simulator_(simulator)
    {
        for (auto& region_pair : this->regions_) {
            this->createLocalRegion_(region_pair.second);
        }

        auto isCartIdxOnThisRank = [&collectToIORank](const int idx) {
            return collectToIORank.isCartIdxOnThisRank(idx);
        };

        this->setupBlockData(isCartIdxOnThisRank);

        this->forceDisableFipOutput_ =
            Parameters::Get<Parameters::ForceDisableFluidInPlaceOutput>();

        this->forceDisableFipresvOutput_ =
            Parameters::Get<Parameters::ForceDisableResvFluidInPlaceOutput>();

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
     * \brief Register all run-time parameters for the Vtk output module.
     */
    static void registerParameters()
    {
        Parameters::Register<Parameters::ForceDisableFluidInPlaceOutput>
            ("Do not print fluid-in-place values after each report step "
             "even if requested by the deck.");
        Parameters::Register<Parameters::ForceDisableResvFluidInPlaceOutput>
            ("Do not print reservoir volumes values after each report step "
             "even if requested by the deck.");
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

        const auto& problem = this->simulator_.problem();

        std::cout << " let us do doAllocBuffers " << std::endl;

        this->doAllocBuffers(bufferSize,
                              reportStepNum,
                              substep,
                              log,
                              isRestart);
        //                      problem.vapparsActive(std::max(simulator_.episodeIndex(), 0)),
        //                      problem.materialLawManager()->enablePCHysteresis(),
        //                      problem.materialLawManager()->enableNonWettingHysteresis(),
        //                      problem.materialLawManager()->enableWettingHysteresis(),
        //                      problem.tracerModel().numTracers(),
        //                      problem.tracerModel().enableSolTracers(),
        //                      problem.eclWriter()->getOutputNnc().size());
    }

//    void processElementMech(const ElementContext& elemCtx)
//    {
//        if constexpr (getPropValue<TypeTag, Properties::EnableMech>()) {
//            const auto& problem = elemCtx.simulator().problem();
//            const auto& model = problem.geoMechModel();
//            for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
//                unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
//                if (!this->mechPotentialForce_.empty()) {
//                    // assume all mechanical things should be written
//                    this->mechPotentialForce_[globalDofIdx] = model.mechPotentialForce(globalDofIdx);
//                    this->mechPotentialPressForce_[globalDofIdx] = model.mechPotentialPressForce(globalDofIdx);
//                    this->mechPotentialTempForce_[globalDofIdx] = model.mechPotentialTempForce(globalDofIdx);
//
//                    this->dispX_[globalDofIdx] = model.disp(globalDofIdx, 0);
//                    this->dispY_[globalDofIdx] = model.disp(globalDofIdx, 1);
//                    this->dispZ_[globalDofIdx] = model.disp(globalDofIdx, 2);
//                    this->stressXX_[globalDofIdx] = model.stress(globalDofIdx, 0);
//                    this->stressYY_[globalDofIdx] = model.stress(globalDofIdx, 1);
//                    this->stressZZ_[globalDofIdx] = model.stress(globalDofIdx, 2);
//                    // voight notation
//                    this->stressXY_[globalDofIdx] = model.stress(globalDofIdx, 5);
//                    this->stressXZ_[globalDofIdx] = model.stress(globalDofIdx, 4);
//                    this->stressYZ_[globalDofIdx] = model.stress(globalDofIdx, 3);
//
//                    this->strainXX_[globalDofIdx] = model.strain(globalDofIdx, 0);
//                    this->strainYY_[globalDofIdx] = model.strain(globalDofIdx, 1);
//                    this->strainZZ_[globalDofIdx] = model.strain(globalDofIdx, 2);
//                    // voight notation
//                    this->strainXY_[globalDofIdx] = model.strain(globalDofIdx, 5);
//                    this->strainXZ_[globalDofIdx] = model.strain(globalDofIdx, 4);
//                    this->strainYZ_[globalDofIdx] = model.strain(globalDofIdx, 3);
//
//
//                    this->delstressXX_[globalDofIdx] = model.delstress(globalDofIdx, 0);
//                    this->delstressYY_[globalDofIdx] = model.delstress(globalDofIdx, 1);
//                    this->delstressZZ_[globalDofIdx] = model.delstress(globalDofIdx, 2);
//                    // voight notation
//                    this->delstressXY_[globalDofIdx] = model.delstress(globalDofIdx, 5);
//                    this->delstressXZ_[globalDofIdx] = model.delstress(globalDofIdx, 4);
//                    this->delstressYZ_[globalDofIdx] = model.delstress(globalDofIdx, 3);
//                }
//            }
//        }
//    }

    /*!
     * \brief Modify the internal buffers according to the intensive
     *        quanties relevant for an element
     */
    void processElement(const ElementContext& elemCtx)
    {
        OPM_TIMEBLOCK_LOCAL(processElement);
        if (!std::is_same<Discretization, EcfvDiscretization<TypeTag>>::value)
            return;

        const auto& problem = elemCtx.simulator().problem();
        const auto& modelResid = elemCtx.simulator().model().linearizer().residual();
        for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
            const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0);
            const auto& fs = intQuants.fluidState();

            using FluidState = std::remove_cv_t<std::remove_reference_t<decltype(fs)>>;

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

//            if (this->regionAvgDensity_.has_value()) {
//                // Note: We intentionally exclude effects of rock
//                // compressibility by using referencePorosity() here.
//                const auto porv = 0; // intQuants.referencePorosity()
//                    // * elemCtx.simulator().model().dofTotalVolume(globalDofIdx);
//
//                this->aggregateAverageDensityContributions_(fs, globalDofIdx,
//                                                            static_cast<double>(porv));
//            }

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

    void processElementFlows(const ElementContext& elemCtx)
    {
        OPM_TIMEBLOCK_LOCAL(processElementBlockData);
        if (!std::is_same_v<Discretization, EcfvDiscretization<TypeTag>>)
            return;
    }

    void processElementBlockData(const ElementContext& elemCtx)
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
    void processFluxes(const ElementContext& elemCtx,
                       ActiveIndex&&         activeIndex,
                       CartesianIndex&&      cartesianIndex)
    {
        OPM_TIMEBLOCK_LOCAL(processFluxes);
        const auto identifyCell = [&activeIndex, &cartesianIndex](const Element& elem)
            -> InterRegFlowMap::Cell
        {
            const auto cellIndex = activeIndex(elem);

            return {
                static_cast<int>(cellIndex),
                cartesianIndex(cellIndex),
                elem.partitionType() == Dune::InteriorEntity
            };
        };

        const auto timeIdx = 0u;
        const auto& stencil = elemCtx.stencil(timeIdx);
        const auto numInteriorFaces = elemCtx.numInteriorFaces(timeIdx);

        for (auto scvfIdx = 0 * numInteriorFaces; scvfIdx < numInteriorFaces; ++scvfIdx) {
            const auto& face = stencil.interiorFace(scvfIdx);
            const auto left  = identifyCell(stencil.element(face.interiorIndex()));
            const auto right = identifyCell(stencil.element(face.exteriorIndex()));

            const auto rates = this->
                getComponentSurfaceRates(elemCtx, face.area(), scvfIdx, timeIdx);

            this->interRegionFlows_.addConnection(left, right, rates);
        }
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

    template <class FluidState>
    void assignToFluidState(FluidState& fs, unsigned elemIdx) const
    {
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (this->saturation_[phaseIdx].empty())
                continue;

            fs.setSaturation(phaseIdx, this->saturation_[phaseIdx][elemIdx]);
        }

        if (!this->fluidPressure_.empty()) {
            // this assumes that capillary pressures only depend on the phase saturations
            // and possibly on temperature. (this is always the case for ECL problems.)
            std::array<Scalar, numPhases> pc = {0};
            const MaterialLawParams& matParams = simulator_.problem().materialLawParams(elemIdx);
            MaterialLaw::capillaryPressures(pc, matParams, fs);
            Valgrind::CheckDefined(this->fluidPressure_[elemIdx]);
            Valgrind::CheckDefined(pc);
            assert(FluidSystem::phaseIsActive(oilPhaseIdx));
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;

                fs.setPressure(phaseIdx, this->fluidPressure_[elemIdx] + (pc[phaseIdx] - pc[oilPhaseIdx]));
            }
        }

        if (!this->temperature_.empty())
            fs.setTemperature(this->temperature_[elemIdx]);
        if (!this->rs_.empty())
            fs.setRs(this->rs_[elemIdx]);
        if (!this->rsw_.empty())
            fs.setRsw(this->rsw_[elemIdx]);
        if (!this->rv_.empty())
            fs.setRv(this->rv_[elemIdx]);
        if (!this->rvw_.empty())
            fs.setRvw(this->rvw_[elemIdx]);
    }

    void initHysteresisParams(Simulator& simulator, unsigned elemIdx) const
    {
        if (!this->soMax_.empty())
            simulator.problem().setMaxOilSaturation(elemIdx, this->soMax_[elemIdx]);

        if (simulator.problem().materialLawManager()->enableHysteresis()) {
            auto matLawManager = simulator.problem().materialLawManager();

            if (FluidSystem::phaseIsActive(oilPhaseIdx) 
                && FluidSystem::phaseIsActive(waterPhaseIdx)) {
                    Scalar somax = 2.0;
                    Scalar swmax = -2.0;
                    Scalar swmin = 2.0;

                if (matLawManager->enableNonWettingHysteresis()) {
                    if (!this->soMax_.empty()) {
                        somax = this->soMax_[elemIdx];
                    }
                }
                if (matLawManager->enableWettingHysteresis()) {
                    if (!this->swMax_.empty()) {
                        swmax = this->swMax_[elemIdx];
                    }
                }
                if (matLawManager->enablePCHysteresis()) {
                    if (!this->swmin_.empty()) {
                        swmin = this->swmin_[elemIdx];
                    }
                }
                matLawManager->setOilWaterHysteresisParams(
                        somax, swmax, swmin, elemIdx);
            }
            if (FluidSystem::phaseIsActive(oilPhaseIdx) 
                && FluidSystem::phaseIsActive(gasPhaseIdx)) {
                    Scalar sgmax = 2.0;
                    Scalar shmax = -2.0;
                    Scalar somin = 2.0;

                if (matLawManager->enableNonWettingHysteresis()) {
                    if (!this->sgmax_.empty()) {
                        sgmax = this->sgmax_[elemIdx];
                    }
                }
                if (matLawManager->enableWettingHysteresis()) {
                    if (!this->shmax_.empty()) {
                        shmax = this->shmax_[elemIdx];
                    }
                }
                if (matLawManager->enablePCHysteresis()) {
                    if (!this->somin_.empty()) {
                        somin = this->somin_[elemIdx];
                    }
                }
                matLawManager->setGasOilHysteresisParams(
                        sgmax, shmax, somin, elemIdx);
            }

        }

        if (simulator_.vanguard().eclState().fieldProps().has_double("SWATINIT")) {
            simulator.problem().materialLawManager()
                ->applyRestartSwatInit(elemIdx, this->ppcw_[elemIdx]);
        }
    }

    void updateFluidInPlace(const ElementContext& elemCtx)
    {
        for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
            updateFluidInPlace_(elemCtx, dofIdx);
        }
    }

    void updateFluidInPlace(const unsigned             globalDofIdx,
                            const IntensiveQuantities& intQuants,
                            const double               totVolume)
    {
        this->updateFluidInPlace_(globalDofIdx, intQuants, totVolume);
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

    void updateFluidInPlace_(const ElementContext& elemCtx, const unsigned dofIdx)
    {
        const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0);
        const unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
        const auto totVolume = elemCtx.simulator().model().dofTotalVolume(globalDofIdx);

        this->updateFluidInPlace_(globalDofIdx, intQuants, totVolume);
    }

    void updateFluidInPlace_(const unsigned             globalDofIdx,
                             const IntensiveQuantities& intQuants,
                             const double               totVolume)
    {
        OPM_TIMEBLOCK_LOCAL(updateFluidInPlace);

        this->updateTotalVolumesAndPressures_(globalDofIdx, intQuants, totVolume);

        if (this->computeFip_) {
            this->updatePhaseInplaceVolumes_(globalDofIdx, intQuants, totVolume);
        }
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

    template <typename FluidState>
    void aggregateAverageDensityContributions_(const FluidState&  fs,
                                               const unsigned int globalDofIdx,
                                               const double       porv)
    {
        auto pvCellValue = RegionPhasePoreVolAverage::CellValue{};
        pvCellValue.porv = porv;

        for (auto phaseIdx = 0*FluidSystem::numPhases;
             phaseIdx < FluidSystem::numPhases; ++phaseIdx)
        {
            if (! FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }

            pvCellValue.value = getValue(fs.density(phaseIdx));
            pvCellValue.sat   = getValue(fs.saturation(phaseIdx));

            this->regionAvgDensity_
                ->addCell(globalDofIdx,
                          RegionPhasePoreVolAverage::Phase { phaseIdx },
                          pvCellValue);
        }
    }

    /*!
     * \brief Compute surface level component flow rates across a single
     *   intersection.
     *
     * \param[in] elemCtx Primary lookup structure for per-cell/element
     *    dynamic information.
     *
     * \param[in] scvfIdx Linear index of current interior bulk connection.
     *
     * \param[in] timeIdx Historical time-point at which to evaluate dynamic
     *    quantities (e.g., reciprocal FVF or dissolved gas concentration).
     *    Zero for the current time.
     *
     * \return Surface level component flow rates.
     */
    data::InterRegFlowMap::FlowRates
    getComponentSurfaceRates(const ElementContext& elemCtx,
                             const Scalar          faceArea,
                             const std::size_t     scvfIdx,
                             const std::size_t     timeIdx) const
    {
        using Component = data::InterRegFlowMap::Component;

        auto rates = data::InterRegFlowMap::FlowRates {};

        return rates;
    }

    template <typename FluidState>
    Scalar hydroCarbonFraction(const FluidState& fs) const
    {
        if (this->eclState_.runspec().co2Storage()) {
            // CO2 storage: Hydrocarbon volume is full pore-volume.
            return 1.0;
        }

        // Common case.  Hydrocarbon volume is fraction occupied by actual
        // hydrocarbons.
        auto hydrocarbon = Scalar {0};
        if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
            hydrocarbon += getValue(fs.saturation(oilPhaseIdx));
        }

        if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
            hydrocarbon += getValue(fs.saturation(gasPhaseIdx));
        }

        return hydrocarbon;
    }

    void updateTotalVolumesAndPressures_(const unsigned             globalDofIdx,
                                         const IntensiveQuantities& intQuants,
                                         const double               totVolume)
    {
        const auto& fs = intQuants.fluidState();

        const double pv = totVolume * intQuants.porosity().value();
        const auto hydrocarbon = this->hydroCarbonFraction(fs);

        if (! this->hydrocarbonPoreVolume_.empty()) {
            this->fip_[Inplace::Phase::PoreVolume][globalDofIdx] = 0;
//                 totVolume * intQuants.referencePorosity();

            this->dynamicPoreVolume_[globalDofIdx] = pv;
            this->hydrocarbonPoreVolume_[globalDofIdx] = pv * hydrocarbon;
        }

        if (!this->pressureTimesHydrocarbonVolume_.empty() &&
            !this->pressureTimesPoreVolume_.empty())
        {
            assert(this->hydrocarbonPoreVolume_.size() == this->pressureTimesHydrocarbonVolume_.size());
            assert(this->fip_[Inplace::Phase::PoreVolume].size() == this->pressureTimesPoreVolume_.size());

            if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
                this->pressureTimesPoreVolume_[globalDofIdx] =
                    getValue(fs.pressure(oilPhaseIdx)) * pv;

                this->pressureTimesHydrocarbonVolume_[globalDofIdx] =
                    this->pressureTimesPoreVolume_[globalDofIdx] * hydrocarbon;
            }
            else if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
                this->pressureTimesPoreVolume_[globalDofIdx] =
                    getValue(fs.pressure(gasPhaseIdx)) * pv;

                this->pressureTimesHydrocarbonVolume_[globalDofIdx] =
                    this->pressureTimesPoreVolume_[globalDofIdx] * hydrocarbon;
            }
            else if (FluidSystem::phaseIsActive(waterPhaseIdx)) {
                this->pressureTimesPoreVolume_[globalDofIdx] =
                    getValue(fs.pressure(waterPhaseIdx)) * pv;
            }
        }
    }

    void updatePhaseInplaceVolumes_(const unsigned             globalDofIdx,
                                    const IntensiveQuantities& intQuants,
                                    const double               totVolume)
    {
    }

    template <typename FIPArray>
    void updateInplaceVolumesSurface(const unsigned  globalDofIdx,
                                     const FIPArray& fip)
    {
        if (FluidSystem::phaseIsActive(oilPhaseIdx) &&
            !this->fip_[Inplace::Phase::OIL].empty())
        {
            this->fip_[Inplace::Phase::OIL][globalDofIdx] = fip[oilPhaseIdx];
        }

        if (FluidSystem::phaseIsActive(oilPhaseIdx) &&
            !this->fip_[Inplace::Phase::OilInLiquidPhase].empty())
        {
            this->fip_[Inplace::Phase::OilInLiquidPhase][globalDofIdx] = fip[oilPhaseIdx];
        }

        if (FluidSystem::phaseIsActive(gasPhaseIdx) &&
            !this->fip_[Inplace::Phase::GAS].empty())
        {
            this->fip_[Inplace::Phase::GAS][globalDofIdx] = fip[gasPhaseIdx];
        }

        if (FluidSystem::phaseIsActive(gasPhaseIdx) &&
            !this->fip_[Inplace::Phase::GasInGasPhase].empty())
        {
            this->fip_[Inplace::Phase::GasInGasPhase][globalDofIdx] = fip[gasPhaseIdx];
        }

        if (FluidSystem::phaseIsActive(waterPhaseIdx) &&
            !this->fip_[Inplace::Phase::WATER].empty())
        {
            this->fip_[Inplace::Phase::WATER][globalDofIdx] = fip[waterPhaseIdx];
        }
    }

    template <typename FluidState, typename FIPArray>
    void updateInplaceVolumesReservoir(const unsigned    globalDofIdx,
                                       const FluidState& fs,
                                       const FIPArray&   fipr)
    {
        if (FluidSystem::phaseIsActive(oilPhaseIdx) &&
            !this->fip_[Inplace::Phase::OilResVolume].empty())
        {
            this->fip_[Inplace::Phase::OilResVolume][globalDofIdx] = fipr[oilPhaseIdx];
        }

        if (FluidSystem::phaseIsActive(gasPhaseIdx) &&
            !this->fip_[Inplace::Phase::GasResVolume].empty())
        {
            this->fip_[Inplace::Phase::GasResVolume][globalDofIdx] = fipr[gasPhaseIdx];
        }

        if (FluidSystem::phaseIsActive(waterPhaseIdx) &&
            !this->fip_[Inplace::Phase::WaterResVolume].empty())
        {
            this->fip_[Inplace::Phase::WaterResVolume][globalDofIdx] = fipr[waterPhaseIdx];
        }

        if (FluidSystem::phaseIsActive(waterPhaseIdx) &&
            !this->fip_[Inplace::Phase::SALT].empty())
        {
            this->fip_[Inplace::Phase::SALT][globalDofIdx] =
                fipr[waterPhaseIdx] * fs.saltConcentration().value();
        }
    }

    template <typename FluidState, typename FIPArray>
    void updateOilGasDistribution(const unsigned    globalDofIdx,
                                  const FluidState& fs,
                                  const FIPArray&   fip)
    {
        // Gas dissolved in oil and vaporized oil
        const auto gasInPlaceLiquid = getValue(fs.Rs()) * fip[oilPhaseIdx];
        const auto oilInPlaceGas    = getValue(fs.Rv()) * fip[gasPhaseIdx];

        if (!this->fip_[Inplace::Phase::GasInLiquidPhase].empty()) {
            this->fip_[Inplace::Phase::GasInLiquidPhase][globalDofIdx] = gasInPlaceLiquid;
        }

        if (!this->fip_[Inplace::Phase::OilInGasPhase].empty()) {
            this->fip_[Inplace::Phase::OilInGasPhase][globalDofIdx] = oilInPlaceGas;
        }

        // Add dissolved gas and vaporized oil to total Fip
        if (!this->fip_[Inplace::Phase::OIL].empty()) {
            this->fip_[Inplace::Phase::OIL][globalDofIdx] += oilInPlaceGas;
        }

        if (!this->fip_[Inplace::Phase::GAS].empty()) {
            this->fip_[Inplace::Phase::GAS][globalDofIdx] += gasInPlaceLiquid;
        }
    }

    template <typename FluidState, typename FIPArray>
    void updateGasWaterDistribution(const unsigned    globalDofIdx,
                                    const FluidState& fs,
                                    const FIPArray&   fip)
    {
        // Gas dissolved in water and vaporized water
        const auto gasInPlaceWater = getValue(fs.Rsw()) * fip[waterPhaseIdx];
        const auto waterInPlaceGas = getValue(fs.Rvw()) * fip[gasPhaseIdx];

        if (!this->fip_[Inplace::Phase::WaterInGasPhase].empty()) {
            this->fip_[Inplace::Phase::WaterInGasPhase][globalDofIdx] = waterInPlaceGas;
        }

        if (!this->fip_[Inplace::Phase::WaterInWaterPhase].empty()) {
            this->fip_[Inplace::Phase::WaterInWaterPhase][globalDofIdx] = fip[waterPhaseIdx];
        }

        // For water+gas cases the gas in water is added to the GIPL value
        if (!this->fip_[Inplace::Phase::GasInLiquidPhase].empty() &&
            !FluidSystem::phaseIsActive(oilPhaseIdx))
        {
            this->fip_[Inplace::Phase::GasInLiquidPhase][globalDofIdx] = gasInPlaceWater;
        }

        // Add dissolved gas and vaporized water to total Fip
        if (!this->fip_[Inplace::Phase::WATER].empty()) {
            this->fip_[Inplace::Phase::WATER][globalDofIdx] += waterInPlaceGas;
        }

        if (!this->fip_[Inplace::Phase::GAS].empty()) {
            this->fip_[Inplace::Phase::GAS][globalDofIdx] += gasInPlaceWater;
        }
    }

    template <typename IntensiveQuantities>
    void updateCO2InGas(const unsigned    globalDofIdx,
                        const double      pv,
                        const IntensiveQuantities& intQuants)
    {
        const auto& scaledDrainageInfo = this->simulator_.problem().materialLawManager()
            ->oilWaterScaledEpsInfoDrainage(globalDofIdx);

        const auto& fs = intQuants.fluidState();
        Scalar sgcr = scaledDrainageInfo.Sgcr;
        if (this->simulator_.problem().materialLawManager()->enableHysteresis()) {
            const auto& matParams = simulator_.problem().materialLawParams(globalDofIdx);
            sgcr = MaterialLaw::trappedGasSaturation(matParams, /*maximumTrapping*/false);
        }

        const Scalar sg   = getValue(fs.saturation(gasPhaseIdx));
        const Scalar rhog = getValue(fs.density(gasPhaseIdx));
        const Scalar xgW  = FluidSystem::phaseIsActive(waterPhaseIdx)
            ? FluidSystem::convertRvwToXgW(getValue(fs.Rvw()), fs.pvtRegionIndex())
            : FluidSystem::convertRvToXgO(getValue(fs.Rv()), fs.pvtRegionIndex());

        const Scalar mM = FluidSystem::molarMass(gasCompIdx, fs.pvtRegionIndex());
        const Scalar massGas = (1 - xgW) * pv * rhog;
        if (!this->fip_[Inplace::Phase::CO2Mass].empty()) {
            this->fip_[Inplace::Phase::CO2Mass][globalDofIdx] = massGas * sg;
        }

        if (!this->fip_[Inplace::Phase::CO2MassInGasPhase].empty()) {
            this->fip_[Inplace::Phase::CO2MassInGasPhase][globalDofIdx] = massGas * sg;
        }

        if (!this->fip_[Inplace::Phase::CO2InGasPhaseInMob].empty()) {
            const Scalar imMobileGas = massGas / mM * std::min(sgcr , sg);
            this->fip_[Inplace::Phase::CO2InGasPhaseInMob][globalDofIdx] = imMobileGas;
        }

        if (!this->fip_[Inplace::Phase::CO2InGasPhaseMob].empty()) {
            const Scalar mobileGas = massGas / mM * std::max(Scalar{0.0}, sg - sgcr);
            this->fip_[Inplace::Phase::CO2InGasPhaseMob][globalDofIdx] = mobileGas;
        }

        if (!this->fip_[Inplace::Phase::CO2InGasPhaseInMobKrg].empty()) {
            if (sgcr >= sg) {
                const Scalar imMobileGasKrg = massGas / mM * sg;
                this->fip_[Inplace::Phase::CO2InGasPhaseInMobKrg][globalDofIdx] = imMobileGasKrg;
            } else {
                this->fip_[Inplace::Phase::CO2InGasPhaseInMobKrg][globalDofIdx] = 0;
            }
        }

        if (!this->fip_[Inplace::Phase::CO2InGasPhaseMobKrg].empty()) {
            if (sg > sgcr) {
                const Scalar mobileGasKrg = massGas / mM * sg;
                this->fip_[Inplace::Phase::CO2InGasPhaseMobKrg][globalDofIdx] = mobileGasKrg;
            } else {
                this->fip_[Inplace::Phase::CO2InGasPhaseMobKrg][globalDofIdx] = 0;
            }
        }

        if (!this->fip_[Inplace::Phase::CO2MassInGasPhaseInMob].empty()) {
            const Scalar imMobileMassGas = massGas * std::min(sgcr , sg);
            this->fip_[Inplace::Phase::CO2MassInGasPhaseInMob][globalDofIdx] = imMobileMassGas;
        }

        if (!this->fip_[Inplace::Phase::CO2MassInGasPhaseMob].empty()) {
            const Scalar mobileMassGas = massGas * std::max(Scalar{0.0}, sg - sgcr);
            this->fip_[Inplace::Phase::CO2MassInGasPhaseMob][globalDofIdx] = mobileMassGas;
        }

        if (!this->fip_[Inplace::Phase::CO2MassInGasPhaseInMobKrg].empty()) {
            if (sgcr >= sg) {
                const Scalar imMobileMassGasKrg = massGas * sg;
                this->fip_[Inplace::Phase::CO2MassInGasPhaseInMobKrg][globalDofIdx] = imMobileMassGasKrg;
            } else {
                this->fip_[Inplace::Phase::CO2MassInGasPhaseInMobKrg][globalDofIdx] = 0;
            }
        }

        if (!this->fip_[Inplace::Phase::CO2MassInGasPhaseMobKrg].empty()) {
            if (sg > sgcr) {
                const Scalar mobileMassGasKrg = massGas * sg;
                this->fip_[Inplace::Phase::CO2MassInGasPhaseMobKrg][globalDofIdx] = mobileMassGasKrg;
            } else {
                this->fip_[Inplace::Phase::CO2MassInGasPhaseMobKrg][globalDofIdx] = 0;
            }
        }

        if (!this->fip_[Inplace::Phase::CO2MassInGasPhaseMaximumTrapped].empty() || 
            !this->fip_[Inplace::Phase::CO2MassInGasPhaseMaximumUnTrapped].empty() ) {
            Scalar trappedGasSaturation = scaledDrainageInfo.Sgcr;
            if (this->simulator_.problem().materialLawManager()->enableHysteresis()) {
                const auto& matParams = simulator_.problem().materialLawParams(globalDofIdx);
                // Get the maximum trapped gas saturation 
                trappedGasSaturation = MaterialLaw::trappedGasSaturation(matParams, /*maximumTrapping*/true);
            }
            if (!this->fip_[Inplace::Phase::CO2MassInGasPhaseMaximumTrapped].empty()) {
                const Scalar imMobileMassGas = massGas * std::min(trappedGasSaturation , sg);
                this->fip_[Inplace::Phase::CO2MassInGasPhaseMaximumTrapped][globalDofIdx] = imMobileMassGas;
            }
            if (!this->fip_[Inplace::Phase::CO2MassInGasPhaseMaximumUnTrapped].empty()) {
                const Scalar mobileMassGas = massGas * std::max(Scalar{0.0}, sg - trappedGasSaturation);
                this->fip_[Inplace::Phase::CO2MassInGasPhaseMaximumUnTrapped][globalDofIdx] = mobileMassGas;
            }
        }

        if (!this->fip_[Inplace::Phase::CO2MassInGasPhaseEffectiveTrapped].empty() || 
            !this->fip_[Inplace::Phase::CO2MassInGasPhaseEffectiveUnTrapped].empty()) {
            Scalar trappedGasSaturation = scaledDrainageInfo.Sgcr;
            if (this->simulator_.problem().materialLawManager()->enableHysteresis()) {
                const auto& matParams = simulator_.problem().materialLawParams(globalDofIdx);
                const double krg = getValue(intQuants.relativePermeability(gasPhaseIdx));
                trappedGasSaturation = MaterialLaw::strandedGasSaturation(matParams, sg, krg);
            }
            if (!this->fip_[Inplace::Phase::CO2MassInGasPhaseEffectiveTrapped].empty()) {
                const Scalar imMobileMassGas = massGas * std::min(trappedGasSaturation , sg);
                this->fip_[Inplace::Phase::CO2MassInGasPhaseEffectiveTrapped][globalDofIdx] = imMobileMassGas;
            }
            if (!this->fip_[Inplace::Phase::CO2MassInGasPhaseEffectiveUnTrapped].empty()) {
                const Scalar mobileMassGas = massGas * std::max(Scalar{0.0}, sg - trappedGasSaturation);
                this->fip_[Inplace::Phase::CO2MassInGasPhaseEffectiveUnTrapped][globalDofIdx] = mobileMassGas;
            }
        }
    }

    template <typename FluidState>
    void updateCO2InWater(const unsigned    globalDofIdx,
                          const double      pv,
                          const FluidState& fs)
    {
        const auto co2InWater = FluidSystem::phaseIsActive(oilPhaseIdx)
            ? this->co2InWaterFromOil(fs, pv)
            : this->co2InWaterFromWater(fs, pv);

        const Scalar mM = FluidSystem::molarMass(gasCompIdx, fs.pvtRegionIndex());
        if (!this->fip_[Inplace::Phase::CO2Mass].empty()) {
            this->fip_[Inplace::Phase::CO2Mass][globalDofIdx] += co2InWater  * mM;
        }
        if (!this->fip_[Inplace::Phase::CO2MassInWaterPhase].empty()) {
            this->fip_[Inplace::Phase::CO2MassInWaterPhase][globalDofIdx] = co2InWater  * mM;
        }
        if (!this->fip_[Inplace::Phase::CO2InWaterPhase].empty()) {
            this->fip_[Inplace::Phase::CO2InWaterPhase][globalDofIdx] = co2InWater;
        }
    }

    template <typename FluidState>
    Scalar co2InWaterFromWater(const FluidState& fs, const double pv) const
    {
        const double rhow = getValue(fs.density(waterPhaseIdx));
        const double sw   = getValue(fs.saturation(waterPhaseIdx));
        const double xwG  = FluidSystem::convertRswToXwG(getValue(fs.Rsw()), fs.pvtRegionIndex());

        const Scalar mM = FluidSystem::molarMass(gasCompIdx, fs.pvtRegionIndex());

        return xwG * pv * rhow * sw / mM;
    }

    template <typename FluidState>
    Scalar co2InWaterFromOil(const FluidState& fs, const double pv) const
    {
        const double rhoo = getValue(fs.density(oilPhaseIdx));
        const double so   = getValue(fs.saturation(oilPhaseIdx));
        const double xoG  = FluidSystem::convertRsToXoG(getValue(fs.Rs()), fs.pvtRegionIndex());

        const Scalar mM = FluidSystem::molarMass(gasCompIdx, fs.pvtRegionIndex());

        return xoG * pv * rhoo * so / mM;
    }

    const Simulator& simulator_;
};

} // namespace Opm

#endif // OPM_OUTPUT_BLACK_OIL_MODULE_HPP
