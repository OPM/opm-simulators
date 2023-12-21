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
 * \copydoc Opm::EclOutputBlackOilModule
 */
#ifndef EWOMS_ECL_OUTPUT_BLACK_OIL_MODULE_HH
#define EWOMS_ECL_OUTPUT_BLACK_OIL_MODULE_HH

#include <dune/common/fvector.hh>

#include <ebos/eclgenericoutputblackoilmodule.hh>

#include <opm/common/Exceptions.hpp>
#include <opm/common/TimingMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/fluidmatrixinteractions/EclEpsScalingPoints.hpp>
#include <opm/material/fluidstates/BlackOilFluidState.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/models/blackoil/blackoilproperties.hh>
#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/models/utils/propertysystem.hh>

#include <opm/output/data/Cells.hpp>
#include <opm/output/eclipse/EclipseIO.hpp>
#include <opm/output/eclipse/Inplace.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <initializer_list>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace Opm::Properties
{

// create new type tag for the Ecl-output
namespace TTag
{
    struct EclOutputBlackOil {
    };
} // namespace TTag

template <class TypeTag, class MyTypeTag>
struct ForceDisableFluidInPlaceOutput {
    using type = UndefinedProperty;
};

template <class TypeTag>
struct ForceDisableFluidInPlaceOutput<TypeTag, TTag::EclOutputBlackOil> {
    static constexpr bool value = false;
};


template <class TypeTag, class MyTypeTag>
struct ForceDisableResvFluidInPlaceOutput {
    using type = UndefinedProperty;
};

template <class TypeTag>
struct ForceDisableResvFluidInPlaceOutput<TypeTag, TTag::EclOutputBlackOil> {
    static constexpr bool value = false;
};

} // namespace Opm::Properties

namespace Opm
{

// forward declaration
template <class TypeTag>
class EcfvDiscretization;

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief Output module for the results black oil model writing in
 *        ECL binary format.
 */
template <class TypeTag>
class EclOutputBlackOilModule : public EclGenericOutputBlackoilModule<GetPropType<TypeTag, Properties::FluidSystem>,
                                                                      GetPropType<TypeTag, Properties::Scalar>>
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
    using BaseType = EclGenericOutputBlackoilModule<FluidSystem, Scalar>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using Dir = FaceDir::DirEnum;

    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { numPhases = FluidSystem::numPhases };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { gasCompIdx = FluidSystem::gasCompIdx };
    enum { oilCompIdx = FluidSystem::oilCompIdx };
    enum { waterCompIdx = FluidSystem::waterCompIdx };
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };

public:
    template <class CollectDataToIORankType>
    EclOutputBlackOilModule(const Simulator& simulator,
                            const CollectDataToIORankType& collectToIORank)
        : BaseType(simulator.vanguard().eclState(),
                   simulator.vanguard().schedule(),
                   simulator.vanguard().summaryConfig(),
                   simulator.vanguard().summaryState(),
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
            EWOMS_GET_PARAM(TypeTag, bool, ForceDisableFluidInPlaceOutput);

        this->forceDisableFipresvOutput_ =
            EWOMS_GET_PARAM(TypeTag, bool, ForceDisableResvFluidInPlaceOutput);

        if (! EWOMS_GET_PARAM(TypeTag, bool, OwnerCellsFirst)) {
            const std::string msg = "The output code does not support --owner-cells-first=false.";
            if (collectToIORank.isIORank()) {
                OpmLog::error(msg);
            }
            OPM_THROW_NOLOG(std::runtime_error, msg);
        }
    }

    /*!
     * \brief Register all run-time parameters for the Vtk output module.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(
            TypeTag,
            bool,
            ForceDisableFluidInPlaceOutput,
            "Do not print fluid-in-place values after each report step even if requested by the deck.");
        EWOMS_REGISTER_PARAM(
            TypeTag,
            bool,
            ForceDisableResvFluidInPlaceOutput,
            "Do not print reservoir volumes values after each report step even if requested by the deck.");
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

        this->doAllocBuffers(bufferSize,
                             reportStepNum,
                             substep,
                             log,
                             isRestart,
                             problem.vapparsActive(std::max(simulator_.episodeIndex(), 0)),
                             problem.materialLawManager()->enableHysteresis(),
                             problem.tracerModel().numTracers(),
                             problem.eclWriter()->getOutputNnc().size());
    }

    void processElementMech(const ElementContext& elemCtx)
    {
        if constexpr (getPropValue<TypeTag, Properties::EnableMech>()) {
            const auto& problem = elemCtx.simulator().problem();
            const auto& model = problem.geoMechModel();
            for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
                unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
                if (!this->mechPotentialForce_.empty()) {
                    // assume all mechanical things should be written
                    this->mechPotentialForce_[globalDofIdx] = model.mechPotentialForce(globalDofIdx);
                    this->mechPotentialPressForce_[globalDofIdx] = model.mechPotentialPressForce(globalDofIdx);
                    this->mechPotentialTempForce_[globalDofIdx] = model.mechPotentialTempForce(globalDofIdx);

                    this->dispX_[globalDofIdx] = model.disp(globalDofIdx, 0);
                    this->dispY_[globalDofIdx] = model.disp(globalDofIdx, 1);
                    this->dispZ_[globalDofIdx] = model.disp(globalDofIdx, 2);
                    this->stressXX_[globalDofIdx] = model.stress(globalDofIdx, 0);
                    this->stressYY_[globalDofIdx] = model.stress(globalDofIdx, 1);
                    this->stressZZ_[globalDofIdx] = model.stress(globalDofIdx, 2);
                    // voight notation
                    this->stressXY_[globalDofIdx] = model.stress(globalDofIdx, 5);
                    this->stressXZ_[globalDofIdx] = model.stress(globalDofIdx, 4);
                    this->stressYZ_[globalDofIdx] = model.stress(globalDofIdx, 3);

                    this->strainXX_[globalDofIdx] = model.strain(globalDofIdx, 0);
                    this->strainYY_[globalDofIdx] = model.strain(globalDofIdx, 1);
                    this->strainZZ_[globalDofIdx] = model.strain(globalDofIdx, 2);
                    // voight notation
                    this->strainXY_[globalDofIdx] = model.strain(globalDofIdx, 5);
                    this->strainXZ_[globalDofIdx] = model.strain(globalDofIdx, 4);
                    this->strainYZ_[globalDofIdx] = model.strain(globalDofIdx, 3);


                    this->delstressXX_[globalDofIdx] = model.delstress(globalDofIdx, 0);
                    this->delstressYY_[globalDofIdx] = model.delstress(globalDofIdx, 1);
                    this->delstressZZ_[globalDofIdx] = model.delstress(globalDofIdx, 2);
                    // voight notation
                    this->delstressXY_[globalDofIdx] = model.delstress(globalDofIdx, 5);
                    this->delstressXZ_[globalDofIdx] = model.delstress(globalDofIdx, 4);
                    this->delstressYZ_[globalDofIdx] = model.delstress(globalDofIdx, 3);
                }
            }
        }
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

        const auto& problem = elemCtx.simulator().problem();
        const auto& ebosResid = elemCtx.simulator().model().linearizer().residual();
        for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
            const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0);
            const auto& fs = intQuants.fluidState();

            typedef typename std::remove_const<typename std::remove_reference<decltype(fs)>::type>::type FluidState;
            unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
            unsigned pvtRegionIdx = elemCtx.primaryVars(dofIdx, /*timeIdx=*/0).pvtRegionIndex();

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (this->saturation_[phaseIdx].empty())
                    continue;

                this->saturation_[phaseIdx][globalDofIdx] = getValue(fs.saturation(phaseIdx));
                Valgrind::CheckDefined(this->saturation_[phaseIdx][globalDofIdx]);
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
            if (!this->gasDissolutionFactor_.empty()) {
                Scalar SoMax = elemCtx.problem().maxOilSaturation(globalDofIdx);
                this->gasDissolutionFactor_[globalDofIdx]
                    = FluidSystem::template saturatedDissolutionFactor<FluidState, Scalar>(
                        fs, oilPhaseIdx, pvtRegionIdx, SoMax);
                Valgrind::CheckDefined(this->gasDissolutionFactor_[globalDofIdx]);
            }
            if (!this->oilVaporizationFactor_.empty()) {
                Scalar SoMax = elemCtx.problem().maxOilSaturation(globalDofIdx);
                this->oilVaporizationFactor_[globalDofIdx]
                    = FluidSystem::template saturatedDissolutionFactor<FluidState, Scalar>(
                        fs, gasPhaseIdx, pvtRegionIdx, SoMax);
                Valgrind::CheckDefined(this->oilVaporizationFactor_[globalDofIdx]);
            }
            if (!this->gasFormationVolumeFactor_.empty()) {
                this->gasFormationVolumeFactor_[globalDofIdx] = 1.0
                    / FluidSystem::template inverseFormationVolumeFactor<FluidState, Scalar>(
                                                                    fs, gasPhaseIdx, pvtRegionIdx);
                Valgrind::CheckDefined(this->gasFormationVolumeFactor_[globalDofIdx]);
            }
            if (!this->saturatedOilFormationVolumeFactor_.empty()) {
                this->saturatedOilFormationVolumeFactor_[globalDofIdx] = 1.0
                    / FluidSystem::template saturatedInverseFormationVolumeFactor<FluidState, Scalar>(
                                                                             fs, oilPhaseIdx, pvtRegionIdx);
                Valgrind::CheckDefined(this->saturatedOilFormationVolumeFactor_[globalDofIdx]);
            }
            if (!this->oilSaturationPressure_.empty()) {
                this->oilSaturationPressure_[globalDofIdx]
                    = FluidSystem::template saturationPressure<FluidState, Scalar>(fs, oilPhaseIdx, pvtRegionIdx);
                Valgrind::CheckDefined(this->oilSaturationPressure_[globalDofIdx]);
            }

            if (!this->rs_.empty()) {
                this->rs_[globalDofIdx] = getValue(fs.Rs());
                Valgrind::CheckDefined(this->rs_[globalDofIdx]);
            }
            if (!this->rsw_.empty()) {
                this->rsw_[globalDofIdx] = getValue(fs.Rsw());
                Valgrind::CheckDefined(this->rsw_[globalDofIdx]);
            }

            if (!this->rv_.empty()) {
                this->rv_[globalDofIdx] = getValue(fs.Rv());
                Valgrind::CheckDefined(this->rv_[globalDofIdx]);
            }
            if (!this->pcow_.empty()) {
                this->pcow_[globalDofIdx] = getValue(fs.pressure(oilPhaseIdx)) - getValue(fs.pressure(waterPhaseIdx));
                Valgrind::CheckDefined(this->pcow_[globalDofIdx]);
            }
            if (!this->pcog_.empty()) {
                this->pcog_[globalDofIdx] = getValue(fs.pressure(gasPhaseIdx)) - getValue(fs.pressure(oilPhaseIdx));
                Valgrind::CheckDefined(this->pcog_[globalDofIdx]);
            }

            if (!this->rvw_.empty()) {
                this->rvw_[globalDofIdx] = getValue(fs.Rvw());
                Valgrind::CheckDefined(this->rvw_[globalDofIdx]);
            }

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (this->invB_[phaseIdx].empty())
                    continue;

                this->invB_[phaseIdx][globalDofIdx] = getValue(fs.invB(phaseIdx));
                Valgrind::CheckDefined(this->invB_[phaseIdx][globalDofIdx]);
            }

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (this->density_[phaseIdx].empty())
                    continue;

                this->density_[phaseIdx][globalDofIdx] = getValue(fs.density(phaseIdx));
                Valgrind::CheckDefined(this->density_[phaseIdx][globalDofIdx]);
            }

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (this->viscosity_[phaseIdx].empty())
                    continue;

                if (!this->extboX_.empty() && phaseIdx == oilPhaseIdx)
                    this->viscosity_[phaseIdx][globalDofIdx] = getValue(intQuants.oilViscosity());
                else if (!this->extboX_.empty() && phaseIdx == gasPhaseIdx)
                    this->viscosity_[phaseIdx][globalDofIdx] = getValue(intQuants.gasViscosity());
                else
                    this->viscosity_[phaseIdx][globalDofIdx] = getValue(fs.viscosity(phaseIdx));
                Valgrind::CheckDefined(this->viscosity_[phaseIdx][globalDofIdx]);
            }

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (this->relativePermeability_[phaseIdx].empty())
                    continue;

                this->relativePermeability_[phaseIdx][globalDofIdx]
                    = getValue(intQuants.relativePermeability(phaseIdx));
                Valgrind::CheckDefined(this->relativePermeability_[phaseIdx][globalDofIdx]);
            }

            if (!this->drsdtcon_.empty()) {
                this->drsdtcon_[globalDofIdx] = problem.drsdtcon(globalDofIdx, elemCtx.simulator().episodeIndex());
            }

            if (!this->sSol_.empty()) {
                this->sSol_[globalDofIdx] = intQuants.solventSaturation().value();
            }

            if (!this->rswSol_.empty()) {
                this->rswSol_[globalDofIdx] = intQuants.rsSolw().value();
            }

            if (!this->cPolymer_.empty()) {
                this->cPolymer_[globalDofIdx] = intQuants.polymerConcentration().value();
            }

            if (!this->cFoam_.empty()) {
                this->cFoam_[globalDofIdx] = intQuants.foamConcentration().value();
            }

            if (!this->cSalt_.empty()) {
                this->cSalt_[globalDofIdx] = fs.saltConcentration().value();
            }

            if (!this->pSalt_.empty()) {
                this->pSalt_[globalDofIdx] = intQuants.saltSaturation().value();
            }

            if (!this->permFact_.empty()) {
                this->permFact_[globalDofIdx] = intQuants.permFactor().value();
            }

            if (!this->extboX_.empty()) {
                this->extboX_[globalDofIdx] = intQuants.xVolume().value();
            }

            if (!this->extboY_.empty()) {
                this->extboY_[globalDofIdx] = intQuants.yVolume().value();
            }

            if (!this->extboZ_.empty()) {
                this->extboZ_[globalDofIdx] = intQuants.zFraction().value();
            }

            if (!this->mFracCo2_.empty()) {
                const Scalar stdVolOil = getValue(fs.saturation(oilPhaseIdx)) * getValue(fs.invB(oilPhaseIdx))
                    + getValue(fs.saturation(gasPhaseIdx)) * getValue(fs.invB(gasPhaseIdx)) * getValue(fs.Rv());
                const Scalar stdVolGas = getValue(fs.saturation(gasPhaseIdx)) * getValue(fs.invB(gasPhaseIdx))
                        * (1.0 - intQuants.yVolume().value())
                    + getValue(fs.saturation(oilPhaseIdx)) * getValue(fs.invB(oilPhaseIdx)) * getValue(fs.Rs())
                        * (1.0 - intQuants.xVolume().value());
                const Scalar stdVolCo2 = getValue(fs.saturation(gasPhaseIdx)) * getValue(fs.invB(gasPhaseIdx))
                        * intQuants.yVolume().value()
                    + getValue(fs.saturation(oilPhaseIdx)) * getValue(fs.invB(oilPhaseIdx)) * getValue(fs.Rs())
                        * intQuants.xVolume().value();
                const Scalar rhoO = FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx);
                const Scalar rhoG = FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx);
                const Scalar rhoCO2 = intQuants.zRefDensity();
                const Scalar stdMassTotal = 1.0e-10 + stdVolOil * rhoO + stdVolGas * rhoG + stdVolCo2 * rhoCO2;
                this->mFracOil_[globalDofIdx] = stdVolOil * rhoO / stdMassTotal;
                this->mFracGas_[globalDofIdx] = stdVolGas * rhoG / stdMassTotal;
                this->mFracCo2_[globalDofIdx] = stdVolCo2 * rhoCO2 / stdMassTotal;
            }

            if (!this->cMicrobes_.empty()) {
                this->cMicrobes_[globalDofIdx] = intQuants.microbialConcentration().value();
            }

            if (!this->cOxygen_.empty()) {
                this->cOxygen_[globalDofIdx] = intQuants.oxygenConcentration().value();
            }

            if (!this->cUrea_.empty()) {
                this->cUrea_[globalDofIdx] = 10
                    * intQuants.ureaConcentration()
                          .value(); // Reescaling back the urea concentration (see WellInterface_impl.hpp)
            }

            if (!this->cBiofilm_.empty()) {
                this->cBiofilm_[globalDofIdx] = intQuants.biofilmConcentration().value();
            }

            if (!this->cCalcite_.empty()) {
                this->cCalcite_[globalDofIdx] = intQuants.calciteConcentration().value();
            }

            if (!this->bubblePointPressure_.empty()) {
                try {
                    this->bubblePointPressure_[globalDofIdx]
                        = getValue(FluidSystem::bubblePointPressure(fs, intQuants.pvtRegionIndex()));
                } catch (const NumericalProblem&) {
                    const auto cartesianIdx = elemCtx.simulator().vanguard().cartesianIndex(globalDofIdx);
                    this->failedCellsPb_.push_back(cartesianIdx);
                }
            }

            if (!this->dewPointPressure_.empty()) {
                try {
                    this->dewPointPressure_[globalDofIdx]
                        = getValue(FluidSystem::dewPointPressure(fs, intQuants.pvtRegionIndex()));
                } catch (const NumericalProblem&) {
                    const auto cartesianIdx = elemCtx.simulator().vanguard().cartesianIndex(globalDofIdx);
                    this->failedCellsPd_.push_back(cartesianIdx);
                }
            }

            if (!this->soMax_.empty())
                this->soMax_[globalDofIdx]
                    = std::max(getValue(fs.saturation(oilPhaseIdx)), problem.maxOilSaturation(globalDofIdx));

            if (!this->swMax_.empty())
                this->swMax_[globalDofIdx]
                    = std::max(getValue(fs.saturation(waterPhaseIdx)), problem.maxWaterSaturation(globalDofIdx));

            if (!this->minimumOilPressure_.empty())
                this->minimumOilPressure_[globalDofIdx]
                    = std::min(getValue(fs.pressure(oilPhaseIdx)), problem.minOilPressure(globalDofIdx));

            if (!this->overburdenPressure_.empty())
                this->overburdenPressure_[globalDofIdx] = problem.overburdenPressure(globalDofIdx);

            if (!this->rockCompPorvMultiplier_.empty())
                this->rockCompPorvMultiplier_[globalDofIdx]
                    = problem.template rockCompPoroMultiplier<Scalar>(intQuants, globalDofIdx);

            if (!this->rockCompTransMultiplier_.empty())
                this->rockCompTransMultiplier_[globalDofIdx]
                    = problem.template rockCompTransMultiplier<Scalar>(intQuants, globalDofIdx);

            const auto& matLawManager = problem.materialLawManager();
            if (matLawManager->enableHysteresis()) {
                if (!this->pcSwMdcOw_.empty() && !this->krnSwMdcOw_.empty()) {
                    matLawManager->oilWaterHysteresisParams(
                        this->pcSwMdcOw_[globalDofIdx], this->krnSwMdcOw_[globalDofIdx], globalDofIdx);
                }
                if (!this->pcSwMdcGo_.empty() && !this->krnSwMdcGo_.empty()) {
                    matLawManager->gasOilHysteresisParams(
                        this->pcSwMdcGo_[globalDofIdx], this->krnSwMdcGo_[globalDofIdx], globalDofIdx);
                }
            }

            if (!this->ppcw_.empty()) {
                this->ppcw_[globalDofIdx] = matLawManager->oilWaterScaledEpsInfoDrainage(globalDofIdx).maxPcow;
                // printf("ppcw_[%d] = %lg\n", globalDofIdx, ppcw_[globalDofIdx]);
            }

            // hack to make the intial output of rs and rv Ecl compatible.
            // For cells with swat == 1 Ecl outputs; rs = rsSat and rv=rvSat, in all but the initial step
            // where it outputs rs and rv values calculated by the initialization. To be compatible we overwrite
            // rs and rv with the values computed in the initially.
            // Volume factors, densities and viscosities need to be recalculated with the updated rs and rv values.
            // This can be removed when ebos has 100% controll over output
            if ((elemCtx.simulator().episodeIndex() < 0) &&
                FluidSystem::phaseIsActive(oilPhaseIdx) &&
                FluidSystem::phaseIsActive(gasPhaseIdx))
            {
                const auto& fsInitial = problem.initialFluidState(globalDofIdx);

                // use initial rs and rv values
                if (!this->rv_.empty())
                    this->rv_[globalDofIdx] = fsInitial.Rv();

                if (!this->rs_.empty())
                    this->rs_[globalDofIdx] = fsInitial.Rs();

                if (!this->rsw_.empty())
                    this->rsw_[globalDofIdx] = fsInitial.Rsw();

                if (!this->rvw_.empty())
                    this->rvw_[globalDofIdx] = fsInitial.Rvw();

                // re-compute the volume factors, viscosities and densities if asked for
                if (!this->density_[oilPhaseIdx].empty())
                    this->density_[oilPhaseIdx][globalDofIdx]
                        = FluidSystem::density(fsInitial, oilPhaseIdx, intQuants.pvtRegionIndex());

                if (!this->density_[gasPhaseIdx].empty())
                    this->density_[gasPhaseIdx][globalDofIdx]
                        = FluidSystem::density(fsInitial, gasPhaseIdx, intQuants.pvtRegionIndex());

                if (!this->invB_[oilPhaseIdx].empty())
                    this->invB_[oilPhaseIdx][globalDofIdx]
                        = FluidSystem::inverseFormationVolumeFactor(fsInitial, oilPhaseIdx, intQuants.pvtRegionIndex());

                if (!this->invB_[gasPhaseIdx].empty())
                    this->invB_[gasPhaseIdx][globalDofIdx]
                        = FluidSystem::inverseFormationVolumeFactor(fsInitial, gasPhaseIdx, intQuants.pvtRegionIndex());

                if (!this->viscosity_[oilPhaseIdx].empty())
                    this->viscosity_[oilPhaseIdx][globalDofIdx]
                        = FluidSystem::viscosity(fsInitial, oilPhaseIdx, intQuants.pvtRegionIndex());

                if (!this->viscosity_[gasPhaseIdx].empty())
                    this->viscosity_[gasPhaseIdx][globalDofIdx]
                        = FluidSystem::viscosity(fsInitial, gasPhaseIdx, intQuants.pvtRegionIndex());
            }

            // Adding Well RFT data
            const auto cartesianIdx = elemCtx.simulator().vanguard().cartesianIndex(globalDofIdx);
            if (this->oilConnectionPressures_.count(cartesianIdx) > 0) {
                this->oilConnectionPressures_[cartesianIdx] = getValue(fs.pressure(oilPhaseIdx));
            }
            if (this->waterConnectionSaturations_.count(cartesianIdx) > 0) {
                this->waterConnectionSaturations_[cartesianIdx] = getValue(fs.saturation(waterPhaseIdx));
            }
            if (this->gasConnectionSaturations_.count(cartesianIdx) > 0) {
                this->gasConnectionSaturations_[cartesianIdx] = getValue(fs.saturation(gasPhaseIdx));
            }

            // tracers
            const auto& tracerModel = simulator_.problem().tracerModel();
            if (! this->tracerConcentrations_.empty()) {
                for (int tracerIdx = 0; tracerIdx < tracerModel.numTracers(); ++tracerIdx) {
                    if (this->tracerConcentrations_[tracerIdx].empty()) {
                        continue;
                    }

                    this->tracerConcentrations_[tracerIdx][globalDofIdx] =
                        tracerModel.tracerConcentration(tracerIdx, globalDofIdx);
                }
            }

            // output residual
            for ( int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx )
            {
                if (!this->residual_[phaseIdx].empty()) {
                    const unsigned activeCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
                    this->residual_[phaseIdx][globalDofIdx] = ebosResid[globalDofIdx][activeCompIdx];
                }
            }
        }
    }

    void processElementFlows(const ElementContext& elemCtx)
    {
        OPM_TIMEBLOCK_LOCAL(processElementBlockData);
        if (!std::is_same_v<Discretization, EcfvDiscretization<TypeTag>>)
            return;

        const auto& problem = elemCtx.simulator().problem();
        for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {

            unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
            if (!problem.model().linearizer().getFlowsInfo().empty()) {
                const auto& flowsInf = problem.model().linearizer().getFlowsInfo();
                auto flowsInfos = flowsInf[globalDofIdx];
                for (auto& flowsInfo : flowsInfos) {
                    if (flowsInfo.faceId >= 0) {
                        if (!this->flows_[flowsInfo.faceId][gasCompIdx].empty()) {
                            this->flows_[flowsInfo.faceId][gasCompIdx][globalDofIdx]
                                = flowsInfo.flow[conti0EqIdx + Indices::canonicalToActiveComponentIndex(gasCompIdx)];
                        }
                        if (!this->flows_[flowsInfo.faceId][oilCompIdx].empty()) {
                            this->flows_[flowsInfo.faceId][oilCompIdx][globalDofIdx]
                                = flowsInfo.flow[conti0EqIdx + Indices::canonicalToActiveComponentIndex(oilCompIdx)];
                        }
                        if (!this->flows_[flowsInfo.faceId][waterCompIdx].empty()) {
                            this->flows_[flowsInfo.faceId][waterCompIdx][globalDofIdx]
                                = flowsInfo.flow[conti0EqIdx + Indices::canonicalToActiveComponentIndex(waterCompIdx)];
                        }
                    }
                    if (flowsInfo.faceId == -2) {
                        if (!this->flowsn_[gasCompIdx].second.first.empty()) {
                            this->flowsn_[gasCompIdx].second.first[flowsInfo.nncId] = flowsInfo.nncId;
                            this->flowsn_[gasCompIdx].second.second[flowsInfo.nncId]
                                = flowsInfo.flow[conti0EqIdx + Indices::canonicalToActiveComponentIndex(gasCompIdx)];
                        }
                        if (!this->flowsn_[oilCompIdx].second.first.empty()) {
                            this->flowsn_[oilCompIdx].second.first[flowsInfo.nncId] = flowsInfo.nncId;
                            this->flowsn_[oilCompIdx].second.second[flowsInfo.nncId]
                                = flowsInfo.flow[conti0EqIdx + Indices::canonicalToActiveComponentIndex(oilCompIdx)];
                        }
                        if (!this->flowsn_[waterCompIdx].second.first.empty()) {
                            this->flowsn_[waterCompIdx].second.first[flowsInfo.nncId] = flowsInfo.nncId;
                            this->flowsn_[waterCompIdx].second.second[flowsInfo.nncId]
                                = flowsInfo.flow[conti0EqIdx + Indices::canonicalToActiveComponentIndex(waterCompIdx)];
                        }
                    }
                }
            }

            // flores
            if (!problem.model().linearizer().getFloresInfo().empty()) {
                const auto& floresInf = problem.model().linearizer().getFloresInfo();
                auto floresInfos =floresInf[globalDofIdx];
                for (auto& floresInfo : floresInfos) {
                    if (floresInfo.faceId >= 0) {
                        if (!this->flores_[floresInfo.faceId][gasCompIdx].empty()) {
                            this->flores_[floresInfo.faceId][gasCompIdx][globalDofIdx]
                                = floresInfo.flow[conti0EqIdx + Indices::canonicalToActiveComponentIndex(gasCompIdx)];
                        }
                        if (!this->flores_[floresInfo.faceId][oilCompIdx].empty()) {
                            this->flores_[floresInfo.faceId][oilCompIdx][globalDofIdx]
                                = floresInfo.flow[conti0EqIdx + Indices::canonicalToActiveComponentIndex(oilCompIdx)];
                        }
                        if (!this->flores_[floresInfo.faceId][waterCompIdx].empty()) {
                            this->flores_[floresInfo.faceId][waterCompIdx][globalDofIdx]
                                = floresInfo.flow[conti0EqIdx + Indices::canonicalToActiveComponentIndex(waterCompIdx)];
                        }
                    }
                   
                    if (floresInfo.faceId == -2) {
                        if (!this->floresn_[gasCompIdx].second.first.empty()) {
                            this->floresn_[gasCompIdx].second.first[floresInfo.nncId] = floresInfo.nncId;
                            this->floresn_[gasCompIdx].second.second[floresInfo.nncId]
                                = floresInfo.flow[conti0EqIdx + Indices::canonicalToActiveComponentIndex(gasCompIdx)];
                        }
                        if (!this->floresn_[oilCompIdx].second.first.empty()) {
                            this->floresn_[oilCompIdx].second.first[floresInfo.nncId] = floresInfo.nncId;
                            this->floresn_[oilCompIdx].second.second[floresInfo.nncId]
                                = floresInfo.flow[conti0EqIdx + Indices::canonicalToActiveComponentIndex(oilCompIdx)];
                        }
                        if (!this->floresn_[waterCompIdx].second.first.empty()) {
                            this->floresn_[waterCompIdx].second.first[floresInfo.nncId] = floresInfo.nncId;
                            this->floresn_[waterCompIdx].second.second[floresInfo.nncId]
                                = floresInfo.flow[conti0EqIdx + Indices::canonicalToActiveComponentIndex(waterCompIdx)];
                        }
                    }
                }
            }
        }
    }

    void processElementBlockData(const ElementContext& elemCtx)
    {
        OPM_TIMEBLOCK_LOCAL(processElementBlockData);
        if (!std::is_same<Discretization, EcfvDiscretization<TypeTag>>::value)
            return;

        const auto& problem = elemCtx.simulator().problem();
        for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
            // Adding block data
            const auto globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
            const auto cartesianIdx = elemCtx.simulator().vanguard().cartesianIndex(globalDofIdx);
            const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0);
            const auto& fs = intQuants.fluidState();
            for (auto& val : this->blockData_) {
                const auto& key = val.first;
                assert(key.second > 0);
                unsigned int cartesianIdxBlock = key.second - 1;
                if (cartesianIdx == cartesianIdxBlock) {
                    if ((key.first == "BWSAT") || (key.first == "BSWAT"))
                        val.second = getValue(fs.saturation(waterPhaseIdx));
                    else if ((key.first == "BGSAT") || (key.first == "BSGAS"))
                        val.second = getValue(fs.saturation(gasPhaseIdx));
                    else if ((key.first == "BOSAT") || (key.first == "BSOIL"))
                        val.second = getValue(fs.saturation(oilPhaseIdx));
                    else if (key.first == "BNSAT")
                        val.second = intQuants.solventSaturation().value();
                    else if ((key.first == "BPR") || (key.first == "BPRESSUR")) {
                        if (FluidSystem::phaseIsActive(oilPhaseIdx))
                            val.second = getValue(fs.pressure(oilPhaseIdx));
                        else if (FluidSystem::phaseIsActive(gasPhaseIdx))
                            val.second = getValue(fs.pressure(gasPhaseIdx));
                        else if (FluidSystem::phaseIsActive(waterPhaseIdx))
                            val.second = getValue(fs.pressure(waterPhaseIdx));
                    }
                    else if ((key.first == "BTCNFHEA") || (key.first == "BTEMP")) {
                        if (FluidSystem::phaseIsActive(oilPhaseIdx))
                            val.second = getValue(fs.temperature(oilPhaseIdx));
                        else if (FluidSystem::phaseIsActive(gasPhaseIdx))
                            val.second = getValue(fs.temperature(gasPhaseIdx));
                        else if (FluidSystem::phaseIsActive(waterPhaseIdx))
                            val.second = getValue(fs.temperature(waterPhaseIdx));
                    }
                    else if (key.first == "BWKR" || key.first == "BKRW")
                        val.second = getValue(intQuants.relativePermeability(waterPhaseIdx));
                    else if (key.first == "BGKR" || key.first == "BKRG")
                        val.second = getValue(intQuants.relativePermeability(gasPhaseIdx));
                    else if (key.first == "BOKR" || key.first == "BKRO")
                        val.second = getValue(intQuants.relativePermeability(oilPhaseIdx));
                    else if (key.first == "BKROG") {
                        const auto& materialParams = problem.materialLawParams(elemCtx, dofIdx, /* timeIdx = */ 0);
                        const auto krog
                            = MaterialLaw::template relpermOilInOilGasSystem<Evaluation>(materialParams, fs);
                        val.second = getValue(krog);
                    }
                    else if (key.first == "BKROW") {
                        const auto& materialParams = problem.materialLawParams(elemCtx, dofIdx, /* timeIdx = */ 0);
                        const auto krow
                            = MaterialLaw::template relpermOilInOilWaterSystem<Evaluation>(materialParams, fs);
                        val.second = getValue(krow);
                    }
                    else if (key.first == "BWPC")
                        val.second = getValue(fs.pressure(oilPhaseIdx)) - getValue(fs.pressure(waterPhaseIdx));
                    else if (key.first == "BGPC")
                        val.second = getValue(fs.pressure(gasPhaseIdx)) - getValue(fs.pressure(oilPhaseIdx));
                    else if (key.first == "BWPR")
                        val.second = getValue(fs.pressure(waterPhaseIdx));
                    else if (key.first == "BGPR")
                        val.second = getValue(fs.pressure(gasPhaseIdx));
                    else if (key.first == "BVWAT" || key.first == "BWVIS")
                        val.second = getValue(fs.viscosity(waterPhaseIdx));
                    else if (key.first == "BVGAS" || key.first == "BGVIS")
                        val.second = getValue(fs.viscosity(gasPhaseIdx));
                    else if (key.first == "BVOIL" || key.first == "BOVIS")
                        val.second = getValue(fs.viscosity(oilPhaseIdx));
                    else if ((key.first == "BODEN") || (key.first == "BDENO"))
                        val.second = getValue(fs.density(oilPhaseIdx));
                    else if ((key.first == "BGDEN") || (key.first == "BDENG"))
                        val.second = getValue(fs.density(gasPhaseIdx));
                    else if ((key.first == "BWDEN") || (key.first == "BDENW"))
                        val.second = getValue(fs.density(waterPhaseIdx));
                    else if ((key.first == "BRPV") ||
                             (key.first == "BOPV") ||
                             (key.first == "BWPV") ||
                             (key.first == "BGPV"))
                    {
                        if (key.first == "BRPV") {
                            val.second = 1.0;
                        }
                        else if (key.first == "BOPV") {
                            val.second = getValue(fs.saturation(oilPhaseIdx));
                        }
                        else if (key.first == "BWPV") {
                            val.second = getValue(fs.saturation(waterPhaseIdx));
                        }
                        else {
                            val.second = getValue(fs.saturation(gasPhaseIdx));
                        }

                        // Include active pore-volume.
                        val.second *= elemCtx.simulator().model().dofTotalVolume(globalDofIdx)
                            * getValue(intQuants.porosity());
                    }
                    else if (key.first == "BRS")
                        val.second = getValue(fs.Rs());
                    else if (key.first == "BRV")
                        val.second = getValue(fs.Rv());
                    else if ((key.first == "BOIP") || (key.first == "BOIPL") || (key.first == "BOIPG") ||
                             (key.first == "BGIP") || (key.first == "BGIPL") || (key.first == "BGIPG") ||
                             (key.first == "BWIP"))
                    {
                        if ((key.first == "BOIP") || (key.first == "BOIPL")) {
                            val.second = getValue(fs.invB(oilPhaseIdx)) * getValue(fs.saturation(oilPhaseIdx));

                            if (key.first == "BOIP") {
                                val.second += getValue(fs.Rv()) * getValue(fs.invB(gasPhaseIdx))
                                    * getValue(fs.saturation(gasPhaseIdx));
                            }
                        }
                        else if (key.first == "BOIPG") {
                            val.second = getValue(fs.Rv()) * getValue(fs.invB(gasPhaseIdx))
                                * getValue(fs.saturation(gasPhaseIdx));
                        }
                        else if ((key.first == "BGIP") || (key.first == "BGIPG")) {
                            val.second = getValue(fs.invB(gasPhaseIdx)) * getValue(fs.saturation(gasPhaseIdx));

                            if (key.first == "BGIP") {
                                val.second += getValue(fs.Rs()) * getValue(fs.invB(oilPhaseIdx))
                                    * getValue(fs.saturation(oilPhaseIdx));
                            }
                        }
                        else if (key.first == "BGIPL") {
                            val.second = getValue(fs.Rs()) * getValue(fs.invB(oilPhaseIdx))
                                * getValue(fs.saturation(oilPhaseIdx));
                        }
                        else { // BWIP
                            val.second = getValue(fs.invB(waterPhaseIdx)) * getValue(fs.saturation(waterPhaseIdx));
                        }

                        // Include active pore-volume.
                        val.second *= elemCtx.simulator().model().dofTotalVolume(globalDofIdx)
                            * getValue(intQuants.porosity());
                    }
                    else if (key.first == "BFLOWI")
                        val.second = this->flows_[FaceDir::ToIntersectionIndex(Dir::XPlus)][waterCompIdx][globalDofIdx];
                    else if (key.first == "BFLOWJ")
                        val.second = this->flows_[FaceDir::ToIntersectionIndex(Dir::YPlus)][waterCompIdx][globalDofIdx];
                    else if (key.first == "BFLOWK")
                        val.second = this->flows_[FaceDir::ToIntersectionIndex(Dir::ZPlus)][waterCompIdx][globalDofIdx];
                    else {
                        std::string logstring = "Keyword '";
                        logstring.append(key.first);
                        logstring.append("' is unhandled for output to summary file.");
                        OpmLog::warning("Unhandled output keyword", logstring);
                    }
                }
            }
        }
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
            -> EclInterRegFlowMap::Cell
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
    const EclInterRegFlowMap& getInterRegFlows() const
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

            if (!this->pcSwMdcOw_.empty() && !this->krnSwMdcOw_.empty()) {
                matLawManager->setOilWaterHysteresisParams(
                    this->pcSwMdcOw_[elemIdx], this->krnSwMdcOw_[elemIdx], elemIdx);
            }
            if (!this->pcSwMdcGo_.empty() && !this->krnSwMdcGo_.empty()) {
                matLawManager->setGasOilHysteresisParams(
                    this->pcSwMdcGo_[elemIdx], this->krnSwMdcGo_[elemIdx], elemIdx);
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
            if (elem.partitionType() != Dune::InteriorEntity)
                region[elemIdx] = 0;
            ++elemIdx;
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

        const auto& extQuant = elemCtx.extensiveQuantities(scvfIdx, timeIdx);

        const auto alpha = getValue(extQuant.extrusionFactor()) * faceArea;

        if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
            const auto& up = elemCtx
                .intensiveQuantities(extQuant.upstreamIndex(oilPhaseIdx), timeIdx);

            using FluidState = std::remove_cv_t<std::remove_reference_t<
                decltype(up.fluidState())>>;

            const auto pvtReg = up.pvtRegionIndex();

            const auto bO = getValue(getInvB_<FluidSystem, FluidState, Scalar>
                                     (up.fluidState(), oilPhaseIdx, pvtReg));

            const auto qO = alpha * bO * getValue(extQuant.volumeFlux(oilPhaseIdx));

            rates[Component::Oil] += qO;

            if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
                const auto Rs = getValue(
                    BlackOil::getRs_<FluidSystem, FluidState, Scalar>
                    (up.fluidState(), pvtReg));

                rates[Component::Gas]    += qO * Rs;
                rates[Component::Disgas] += qO * Rs;
            }
        }

        if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
            const auto& up = elemCtx
                .intensiveQuantities(extQuant.upstreamIndex(gasPhaseIdx), timeIdx);

            using FluidState = std::remove_cv_t<std::remove_reference_t<
                decltype(up.fluidState())>>;

            const auto pvtReg = up.pvtRegionIndex();

            const auto bG = getValue(getInvB_<FluidSystem, FluidState, Scalar>
                                     (up.fluidState(), gasPhaseIdx, pvtReg));

            const auto qG = alpha * bG * getValue(extQuant.volumeFlux(gasPhaseIdx));

            rates[Component::Gas] += qG;

            if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
                const auto Rv = getValue(
                    BlackOil::getRv_<FluidSystem, FluidState, Scalar>
                    (up.fluidState(), pvtReg));

                rates[Component::Oil]    += qG * Rv;
                rates[Component::Vapoil] += qG * Rv;
            }
        }

        if (FluidSystem::phaseIsActive(waterPhaseIdx)) {
            const auto& up = elemCtx
                .intensiveQuantities(extQuant.upstreamIndex(waterPhaseIdx), timeIdx);

            using FluidState = std::remove_cv_t<std::remove_reference_t<
                decltype(up.fluidState())>>;

            const auto pvtReg = up.pvtRegionIndex();

            const auto bW = getValue(getInvB_<FluidSystem, FluidState, Scalar>
                                     (up.fluidState(), waterPhaseIdx, pvtReg));

            rates[Component::Water] +=
                alpha * bW * getValue(extQuant.volumeFlux(waterPhaseIdx));
        }

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
            this->fip_[Inplace::Phase::PoreVolume][globalDofIdx] =
                totVolume * intQuants.referencePorosity();

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
        std::array<Scalar, FluidSystem::numPhases> fip {};
        std::array<Scalar, FluidSystem::numPhases> fipr{}; // at reservoir condition

        const auto& fs = intQuants.fluidState();
        const auto  pv = totVolume * intQuants.porosity().value();

        for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }

            const auto b = getValue(fs.invB(phaseIdx));
            const auto s = getValue(fs.saturation(phaseIdx));

            fipr[phaseIdx] = s * pv;
            fip [phaseIdx] = b * fipr[phaseIdx];
        }

        this->updateInplaceVolumesSurface(globalDofIdx, fip);
        this->updateInplaceVolumesReservoir(globalDofIdx, fs, fipr);

        if (FluidSystem::phaseIsActive(oilPhaseIdx) &&
            FluidSystem::phaseIsActive(gasPhaseIdx))
        {
            this->updateOilGasDistribution(globalDofIdx, fs, fip);
        }

        if (FluidSystem::phaseIsActive(waterPhaseIdx) &&
            FluidSystem::phaseIsActive(gasPhaseIdx))
        {
            this->updateGasWaterDistribution(globalDofIdx, fs, fip);
        }

        if (FluidSystem::phaseIsActive(gasPhaseIdx) &&
            (!this->fip_[Inplace::Phase::CO2InGasPhaseInMob].empty() ||
             !this->fip_[Inplace::Phase::CO2InGasPhaseMob].empty()))
        {
            this->updateCO2InGas(globalDofIdx, pv, fs);
        }

        if (!this->fip_[Inplace::Phase::CO2InWaterPhase].empty() &&
            (FluidSystem::phaseIsActive(waterPhaseIdx) ||
             FluidSystem::phaseIsActive(oilPhaseIdx)))
        {
            const auto co2InWater = FluidSystem::phaseIsActive(oilPhaseIdx)
                ? this->co2InWaterFromOil(fs, pv)
                : this->co2InWaterFromWater(fs, pv);

            this->fip_[Inplace::Phase::CO2InWaterPhase][globalDofIdx] = co2InWater;
        }
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

    template <typename FluidState>
    void updateCO2InGas(const unsigned    globalDofIdx,
                        const double      pv,
                        const FluidState& fs)
    {
        const auto& scaledDrainageInfo = this->simulator_.problem().materialLawManager()
            ->oilWaterScaledEpsInfoDrainage(globalDofIdx);

        Scalar sgcr = scaledDrainageInfo.Sgcr;
        if (this->simulator_.problem().materialLawManager()->enableHysteresis()) {
            const auto& matParams = simulator_.problem().materialLawParams(globalDofIdx);
            sgcr = MaterialLaw::trappedGasSaturation(matParams);
        }

        const double sg   = getValue(fs.saturation(gasPhaseIdx));
        const double rhog = getValue(fs.density(gasPhaseIdx));
        const double xgW  = FluidSystem::phaseIsActive(waterPhaseIdx)
            ? FluidSystem::convertRvwToXgW(getValue(fs.Rvw()), fs.pvtRegionIndex())
            : FluidSystem::convertRvToXgO(getValue(fs.Rv()), fs.pvtRegionIndex());

        const Scalar mM = FluidSystem::molarMass(gasCompIdx, fs.pvtRegionIndex());

        if (!this->fip_[Inplace::Phase::CO2InGasPhaseInMob].empty()) {
            const Scalar imMobileGas = (1 - xgW) * pv * rhog / mM * std::min(sgcr , sg);
            this->fip_[Inplace::Phase::CO2InGasPhaseInMob][globalDofIdx] = imMobileGas;
        }

        if (!this->fip_[Inplace::Phase::CO2InGasPhaseMob].empty()) {
            const Scalar mobileGas = (1 - xgW) * pv * rhog / mM * std::max(0.0, sg - sgcr);
            this->fip_[Inplace::Phase::CO2InGasPhaseMob][globalDofIdx] = mobileGas;
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

#endif
