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

#include <array>
#include <numeric>
#include <optional>
#include <stdexcept>

#include <opm/models/blackoil/blackoilproperties.hh>

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>

#include <opm/material/common/Valgrind.hpp>

#include <opm/parser/eclipse/Units/Units.hpp>
#include <opm/parser/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/output/data/Cells.hpp>
#include <opm/output/eclipse/EclipseIO.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/output/eclipse/Inplace.hpp>

#include <ebos/eclgenericoutputblackoilmodule.hh>

#include <dune/common/fvector.hh>

#include <type_traits>

namespace Opm::Properties {

// create new type tag for the Ecl-output
namespace TTag {
struct EclOutputBlackOil {};
}

template<class TypeTag, class MyTypeTag>
struct ForceDisableFluidInPlaceOutput {
    using type = UndefinedProperty;
};

template<class TypeTag>
struct ForceDisableFluidInPlaceOutput<TypeTag, TTag::EclOutputBlackOil> {
    static constexpr bool value = false;
};


template<class TypeTag, class MyTypeTag>
struct ForceDisableResvFluidInPlaceOutput {
    using type = UndefinedProperty;
};

template<class TypeTag>
struct ForceDisableResvFluidInPlaceOutput<TypeTag, TTag::EclOutputBlackOil> {
    static constexpr bool value = false;
};


} // namespace Opm::Properties

namespace Opm {

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
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementIterator = typename GridView::template Codim<0>::Iterator;
    using BaseType = EclGenericOutputBlackoilModule<FluidSystem, Scalar>;

    enum { numPhases = FluidSystem::numPhases };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { gasCompIdx = FluidSystem::gasCompIdx };
    enum { oilCompIdx = FluidSystem::oilCompIdx };
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };

public:
    template<class CollectDataToIORankType>
    EclOutputBlackOilModule(const Simulator& simulator,
                            const std::vector<std::size_t>& wbp_index_list,
                            const CollectDataToIORankType& collectToIORank)
        : BaseType(simulator.vanguard().eclState(),
                   simulator.vanguard().schedule(),
                   simulator.vanguard().summaryConfig(),
                   simulator.vanguard().summaryState(),
                   getPropValue<TypeTag, Properties::EnableEnergy>(),
                   getPropValue<TypeTag, Properties::EnableTemperature>(),
                   getPropValue<TypeTag, Properties::EnableSolvent>(),
                   getPropValue<TypeTag, Properties::EnablePolymer>(),
                   getPropValue<TypeTag, Properties::EnableFoam>(),
                   getPropValue<TypeTag, Properties::EnableBrine>(),
                   getPropValue<TypeTag, Properties::EnableExtbo>())
         , simulator_(simulator)
    {
        const SummaryConfig summaryConfig = simulator_.vanguard().summaryConfig();

        for (auto& region_pair : this->regions_)
            createLocalRegion_(region_pair.second);

        for (const auto& node: summaryConfig) {
            if (node.category() == SummaryConfigNode::Category::Block) {
                if (collectToIORank.isCartIdxOnThisRank(node.number() - 1)) {
                    std::pair<std::string, int> key = std::make_pair(node.keyword(), node.number());
                    this->blockData_[key] = 0.0;
                }
            }
        }

        for (const auto& global_index : wbp_index_list) {
            if (collectToIORank.isCartIdxOnThisRank(global_index - 1))
                this->wbpData_[global_index] = 0.0;
        }

        this->forceDisableFipOutput_ = EWOMS_GET_PARAM(TypeTag, bool, ForceDisableFluidInPlaceOutput);
        this->forceDisableFipresvOutput_ = EWOMS_GET_PARAM(TypeTag, bool, ForceDisableResvFluidInPlaceOutput);
    }

    /*!
     * \brief Register all run-time parameters for the Vtk output module.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, bool, ForceDisableFluidInPlaceOutput,
                             "Do not print fluid-in-place values after each report step even if requested by the deck.");
        EWOMS_REGISTER_PARAM(TypeTag, bool, ForceDisableResvFluidInPlaceOutput,
                             "Do not print reservoir volumes values after each report step even if requested by the deck.");
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to ECL output files
     */
    void allocBuffers(unsigned bufferSize, unsigned reportStepNum, const bool substep, const bool log, const bool isRestart)
    {
        if (!std::is_same<Discretization, EcfvDiscretization<TypeTag> >::value)
            return;

        this->doAllocBuffers(bufferSize,
                             reportStepNum,
                             substep,
                             log,
                             isRestart,
                             simulator_.problem().vapparsActive(std::max(simulator_.episodeIndex(), 0)),
                             simulator_.problem().materialLawManager()->enableHysteresis(),
                             simulator_.problem().tracerModel().numTracers());
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quanties relevant
     *        for an element
     */
    void processElement(const ElementContext& elemCtx)
    {
        if (!std::is_same<Discretization, EcfvDiscretization<TypeTag> >::value)
            return;

        const auto& problem = elemCtx.simulator().problem();
        for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
            const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0);
            const auto& fs = intQuants.fluidState();

            typedef typename std::remove_const<typename std::remove_reference<decltype(fs)>::type>::type FluidState;
            unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
            unsigned pvtRegionIdx = elemCtx.primaryVars(dofIdx, /*timeIdx=*/0).pvtRegionIndex();

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                if (this->saturation_[phaseIdx].empty())
                    continue;

                this->saturation_[phaseIdx][globalDofIdx] = getValue(fs.saturation(phaseIdx));
                Valgrind::CheckDefined(this->saturation_[phaseIdx][globalDofIdx]);
            }

            if (!this->oilPressure_.empty()) {
                if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
                    this->oilPressure_[globalDofIdx] = getValue(fs.pressure(oilPhaseIdx));
                }else{
                    // put pressure in oil pressure for output
                    if (FluidSystem::phaseIsActive(waterPhaseIdx)) {
                        this->oilPressure_[globalDofIdx] = getValue(fs.pressure(waterPhaseIdx));
                    } else {
                        this->oilPressure_[globalDofIdx] = getValue(fs.pressure(gasPhaseIdx));
                    }
                }
                Valgrind::CheckDefined(this->oilPressure_[globalDofIdx]);
            }

            if (!this->temperature_.empty()) {
                this->temperature_[globalDofIdx] = getValue(fs.temperature(oilPhaseIdx));
                Valgrind::CheckDefined(this->temperature_[globalDofIdx]);
            }
            if (!this->gasDissolutionFactor_.empty()) {
                Scalar SoMax = elemCtx.problem().maxOilSaturation(globalDofIdx);
                this->gasDissolutionFactor_[globalDofIdx] =
                    FluidSystem::template saturatedDissolutionFactor<FluidState, Scalar>(fs, oilPhaseIdx, pvtRegionIdx, SoMax);
                Valgrind::CheckDefined(this->gasDissolutionFactor_[globalDofIdx]);

            }
            if (!this->oilVaporizationFactor_.empty()) {
                Scalar SoMax = elemCtx.problem().maxOilSaturation(globalDofIdx);
                this->oilVaporizationFactor_[globalDofIdx] =
                    FluidSystem::template saturatedDissolutionFactor<FluidState, Scalar>(fs, gasPhaseIdx, pvtRegionIdx, SoMax);
                Valgrind::CheckDefined(this->oilVaporizationFactor_[globalDofIdx]);

            }
            if (!this->gasFormationVolumeFactor_.empty()) {
                this->gasFormationVolumeFactor_[globalDofIdx] =
                    1.0/FluidSystem::template inverseFormationVolumeFactor<FluidState, Scalar>(fs, gasPhaseIdx, pvtRegionIdx);
                Valgrind::CheckDefined(this->gasFormationVolumeFactor_[globalDofIdx]);

            }
            if (!this->saturatedOilFormationVolumeFactor_.empty()) {
                this->saturatedOilFormationVolumeFactor_[globalDofIdx] =
                    1.0/FluidSystem::template saturatedInverseFormationVolumeFactor<FluidState, Scalar>(fs, oilPhaseIdx, pvtRegionIdx);
                Valgrind::CheckDefined(this->saturatedOilFormationVolumeFactor_[globalDofIdx]);

            }
            if (!this->oilSaturationPressure_.empty()) {
                this->oilSaturationPressure_[globalDofIdx] =
                    FluidSystem::template saturationPressure<FluidState, Scalar>(fs, oilPhaseIdx, pvtRegionIdx);
                Valgrind::CheckDefined(this->oilSaturationPressure_[globalDofIdx]);

            }

            if (!this->rs_.empty()) {
                this->rs_[globalDofIdx] = getValue(fs.Rs());
                Valgrind::CheckDefined(this->rs_[globalDofIdx]);
            }

            if (!this->rv_.empty()) {
                this->rv_[globalDofIdx] = getValue(fs.Rv());
                Valgrind::CheckDefined(this->rv_[globalDofIdx]);
            }

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                if (this->invB_[phaseIdx].empty())
                    continue;

                this->invB_[phaseIdx][globalDofIdx] = getValue(fs.invB(phaseIdx));
                Valgrind::CheckDefined(this->invB_[phaseIdx][globalDofIdx]);
            }

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                if (this->density_[phaseIdx].empty())
                    continue;

                this->density_[phaseIdx][globalDofIdx] = getValue(fs.density(phaseIdx));
                Valgrind::CheckDefined(this->density_[phaseIdx][globalDofIdx]);
            }

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                if (this->viscosity_[phaseIdx].empty())
                    continue;

                if (!this->extboX_.empty() && phaseIdx==oilPhaseIdx)
                    this->viscosity_[phaseIdx][globalDofIdx] = getValue(intQuants.oilViscosity());
                else if (!this->extboX_.empty() && phaseIdx==gasPhaseIdx)
                    this->viscosity_[phaseIdx][globalDofIdx] = getValue(intQuants.gasViscosity());
                else
                    this->viscosity_[phaseIdx][globalDofIdx] = getValue(fs.viscosity(phaseIdx));
                Valgrind::CheckDefined(this->viscosity_[phaseIdx][globalDofIdx]);
            }

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                if (this->relativePermeability_[phaseIdx].empty())
                    continue;

                this->relativePermeability_[phaseIdx][globalDofIdx] = getValue(intQuants.relativePermeability(phaseIdx));
                Valgrind::CheckDefined(this->relativePermeability_[phaseIdx][globalDofIdx]);
            }

            if (!this->sSol_.empty()) {
                this->sSol_[globalDofIdx] = intQuants.solventSaturation().value();
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
                const Scalar stdVolOil = getValue(fs.saturation(oilPhaseIdx))*getValue(fs.invB(oilPhaseIdx))
                                       + getValue(fs.saturation(gasPhaseIdx))*getValue(fs.invB(gasPhaseIdx))*getValue(fs.Rv());
                const Scalar stdVolGas = getValue(fs.saturation(gasPhaseIdx))*getValue(fs.invB(gasPhaseIdx))*(1.0-intQuants.yVolume().value())
                                       + getValue(fs.saturation(oilPhaseIdx))*getValue(fs.invB(oilPhaseIdx))*getValue(fs.Rs())*(1.0-intQuants.xVolume().value());
                const Scalar stdVolCo2 = getValue(fs.saturation(gasPhaseIdx))*getValue(fs.invB(gasPhaseIdx))*intQuants.yVolume().value()
                                       + getValue(fs.saturation(oilPhaseIdx))*getValue(fs.invB(oilPhaseIdx))*getValue(fs.Rs())*intQuants.xVolume().value();
                const Scalar rhoO= FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx);
                const Scalar rhoG= FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx);
                const Scalar rhoCO2= intQuants.zRefDensity();
                const Scalar stdMassTotal= 1.0e-10 + stdVolOil*rhoO + stdVolGas*rhoG + stdVolCo2*rhoCO2;
                this->mFracOil_[globalDofIdx] = stdVolOil*rhoO/stdMassTotal;
                this->mFracGas_[globalDofIdx] = stdVolGas*rhoG/stdMassTotal;
                this->mFracCo2_[globalDofIdx] = stdVolCo2*rhoCO2/stdMassTotal;
            }

            if (!this->bubblePointPressure_.empty()) {
                try {
                    this->bubblePointPressure_[globalDofIdx] = getValue(FluidSystem::bubblePointPressure(fs, intQuants.pvtRegionIndex()));
                }
                catch (const NumericalIssue&) {
                    const auto cartesianIdx = elemCtx.simulator().vanguard().grid().globalCell()[globalDofIdx];
                    this->failedCellsPb_.push_back(cartesianIdx);
                }
            }
            if (!this->dewPointPressure_.empty()) {
                try {
                    this->dewPointPressure_[globalDofIdx] = getValue(FluidSystem::dewPointPressure(fs, intQuants.pvtRegionIndex()));
                }
                catch (const NumericalIssue&) {
                    const auto cartesianIdx = elemCtx.simulator().vanguard().grid().globalCell()[globalDofIdx];
                    this->failedCellsPd_.push_back(cartesianIdx);
                }
            }

            if (!this->soMax_.empty())
                this->soMax_[globalDofIdx] =
                    std::max(getValue(fs.saturation(oilPhaseIdx)),
                             problem.maxOilSaturation(globalDofIdx));

            if (!this->swMax_.empty())
                this->swMax_[globalDofIdx] =
                    std::max(getValue(fs.saturation(waterPhaseIdx)),
                             problem.maxWaterSaturation(globalDofIdx));

            if (!this->minimumOilPressure_.empty())
                this->minimumOilPressure_[globalDofIdx] =
                    std::min(getValue(fs.pressure(oilPhaseIdx)),
                             problem.minOilPressure(globalDofIdx));

            if (!this->overburdenPressure_.empty())
                this->overburdenPressure_[globalDofIdx] = problem.overburdenPressure(globalDofIdx);

            if (!this->rockCompPorvMultiplier_.empty())
                this->rockCompPorvMultiplier_[globalDofIdx] = problem.template rockCompPoroMultiplier<Scalar>(intQuants, globalDofIdx);

            if (!this->rockCompTransMultiplier_.empty())
                this->rockCompTransMultiplier_[globalDofIdx] = problem.template rockCompTransMultiplier<Scalar>(intQuants, globalDofIdx);

            const auto& matLawManager = problem.materialLawManager();
            if (matLawManager->enableHysteresis()) {
                if (!this->pcSwMdcOw_.empty() && !this->krnSwMdcOw_.empty()) {
                    matLawManager->oilWaterHysteresisParams(
                                this->pcSwMdcOw_[globalDofIdx],
                                this->krnSwMdcOw_[globalDofIdx],
                                globalDofIdx);
                }
                if (!this->pcSwMdcGo_.empty() && !this->krnSwMdcGo_.empty()) {
                    matLawManager->gasOilHysteresisParams(
                                this->pcSwMdcGo_[globalDofIdx],
                                this->krnSwMdcGo_[globalDofIdx],
                                globalDofIdx);
                }
            }


            if (!this->ppcw_.empty()) {
                this->ppcw_[globalDofIdx] = matLawManager->oilWaterScaledEpsInfoDrainage(globalDofIdx).maxPcow;
                //printf("ppcw_[%d] = %lg\n", globalDofIdx, ppcw_[globalDofIdx]);
            }
            // hack to make the intial output of rs and rv Ecl compatible.
            // For cells with swat == 1 Ecl outputs; rs = rsSat and rv=rvSat, in all but the initial step
            // where it outputs rs and rv values calculated by the initialization. To be compatible we overwrite
            // rs and rv with the values computed in the initially.
            // Volume factors, densities and viscosities need to be recalculated with the updated rs and rv values.
            // This can be removed when ebos has 100% controll over output
            if (elemCtx.simulator().episodeIndex() < 0 && FluidSystem::phaseIsActive(oilPhaseIdx) && FluidSystem::phaseIsActive(gasPhaseIdx)) {

                const auto& fsInitial = problem.initialFluidState(globalDofIdx);

                // use initial rs and rv values
                if (!this->rv_.empty())
                    this->rv_[globalDofIdx] = fsInitial.Rv();

                if (!this->rs_.empty())
                    this->rs_[globalDofIdx] = fsInitial.Rs();

                // re-compute the volume factors, viscosities and densities if asked for
                if (!this->density_[oilPhaseIdx].empty())
                    this->density_[oilPhaseIdx][globalDofIdx] = FluidSystem::density(fsInitial,
                                                                                     oilPhaseIdx,
                                                                                     intQuants.pvtRegionIndex());
                if (!this->density_[gasPhaseIdx].empty())
                    this->density_[gasPhaseIdx][globalDofIdx] = FluidSystem::density(fsInitial,
                                                                                     gasPhaseIdx,
                                                                                     intQuants.pvtRegionIndex());

                if (!this->invB_[oilPhaseIdx].empty())
                    this->invB_[oilPhaseIdx][globalDofIdx] = FluidSystem::inverseFormationVolumeFactor(fsInitial,
                                                                                                       oilPhaseIdx,
                                                                                                       intQuants.pvtRegionIndex());
                if (!this->invB_[gasPhaseIdx].empty())
                    this->invB_[gasPhaseIdx][globalDofIdx] = FluidSystem::inverseFormationVolumeFactor(fsInitial,
                                                                                                 gasPhaseIdx,
                                                                                                 intQuants.pvtRegionIndex());
                if (!this->viscosity_[oilPhaseIdx].empty())
                    this->viscosity_[oilPhaseIdx][globalDofIdx] = FluidSystem::viscosity(fsInitial,
                                                                                         oilPhaseIdx,
                                                                                         intQuants.pvtRegionIndex());
                if (!this->viscosity_[gasPhaseIdx].empty())
                    this->viscosity_[gasPhaseIdx][globalDofIdx] = FluidSystem::viscosity(fsInitial,
                                                                                         gasPhaseIdx,
                                                                                         intQuants.pvtRegionIndex());
            }

            // Add fluid in Place values
            updateFluidInPlace_(elemCtx, dofIdx);

            // Adding block data
            const auto cartesianIdx = elemCtx.simulator().vanguard().grid().globalCell()[globalDofIdx];
            for (auto& val: this->blockData_) {
                const auto& key = val.first;
                int cartesianIdxBlock = key.second - 1;
                if (cartesianIdx == cartesianIdxBlock) {
                    if ((key.first == "BWSAT") || (key.first == "BSWAT"))
                        val.second = getValue(fs.saturation(waterPhaseIdx));
                    else if ((key.first == "BGSAT") || (key.first == "BSGAS"))
                        val.second = getValue(fs.saturation(gasPhaseIdx));
                    else if ((key.first == "BOSAT") || (key.first == "BSOIL"))
                        val.second = getValue(fs.saturation(oilPhaseIdx));
                    else if ((key.first == "BPR") || (key.first == "BPRESSUR"))
                        val.second = getValue(fs.pressure(oilPhaseIdx));
                    else if (key.first == "BWKR" || key.first == "BKRW")
                        val.second = getValue(intQuants.relativePermeability(waterPhaseIdx));
                    else if (key.first == "BGKR" || key.first == "BKRG")
                        val.second = getValue(intQuants.relativePermeability(gasPhaseIdx));
                    else if (key.first == "BOKR" || key.first == "BKRO")
                        val.second = getValue(intQuants.relativePermeability(oilPhaseIdx));
                    else if (key.first == "BKROG") {
                        const auto& materialParams = problem.materialLawParams(elemCtx, dofIdx, /* timeIdx = */ 0);
                        const auto krog = MaterialLaw::template relpermOilInOilGasSystem<Evaluation>(materialParams, fs);
                        val.second = getValue(krog);
                    }
                    else if (key.first == "BKROW") {
                        const auto& materialParams = problem.materialLawParams(elemCtx, dofIdx, /* timeIdx = */ 0);
                        const auto krow = MaterialLaw::template relpermOilInOilWaterSystem<Evaluation>(materialParams, fs);
                        val.second = getValue(krow);
                    }
                    else if (key.first == "BWPC")
                        val.second = getValue(fs.pressure(oilPhaseIdx)) - getValue(fs.pressure(waterPhaseIdx));
                    else if (key.first == "BGPC")
                        val.second = getValue(fs.pressure(gasPhaseIdx)) - getValue(fs.pressure(oilPhaseIdx));
                    else if (key.first == "BVWAT" || key.first == "BWVIS")
                        val.second = getValue(fs.viscosity(waterPhaseIdx));
                    else if (key.first == "BVGAS" || key.first == "BGVIS")
                        val.second = getValue(fs.viscosity(gasPhaseIdx));
                    else if (key.first == "BVOIL" || key.first == "BOVIS")
                        val.second = getValue(fs.viscosity(oilPhaseIdx));
                    else {
                        std::string logstring = "Keyword '";
                        logstring.append(key.first);
                        logstring.append("' is unhandled for output to file.");
                        OpmLog::warning("Unhandled output keyword", logstring);
                    }
                }
            }

            // Adding Well RFT data
            if (this->oilConnectionPressures_.count(cartesianIdx) > 0) {
                this->oilConnectionPressures_[cartesianIdx] = getValue(fs.pressure(oilPhaseIdx));
            }
            if (this->waterConnectionSaturations_.count(cartesianIdx) > 0) {
                this->waterConnectionSaturations_[cartesianIdx] = getValue(fs.saturation(waterPhaseIdx));
            }
            if (this->gasConnectionSaturations_.count(cartesianIdx) > 0) {
                this->gasConnectionSaturations_[cartesianIdx] = getValue(fs.saturation(gasPhaseIdx));
            }
            if (this->wbpData_.count(cartesianIdx) > 0)
                this->wbpData_[cartesianIdx] = getValue(fs.pressure(oilPhaseIdx));

            // tracers
            const auto& tracerModel = simulator_.problem().tracerModel();
            if (!this->tracerConcentrations_.empty()) {
                for (int tracerIdx = 0; tracerIdx < tracerModel.numTracers(); tracerIdx++){
                    if (this->tracerConcentrations_[tracerIdx].empty())
                        continue;

                    this->tracerConcentrations_[tracerIdx][globalDofIdx] = tracerModel.tracerConcentration(tracerIdx, globalDofIdx);
                }
            }
        }
    }

    template <class FluidState>
    void assignToFluidState(FluidState& fs, unsigned elemIdx) const
    {
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            if (this->saturation_[phaseIdx].empty())
                continue;

            fs.setSaturation(phaseIdx, this->saturation_[phaseIdx][elemIdx]);
        }

        if (!this->oilPressure_.empty()) {
            // this assumes that capillary pressures only depend on the phase saturations
            // and possibly on temperature. (this is always the case for ECL problems.)
            Dune::FieldVector< Scalar, numPhases > pc(0);
            const MaterialLawParams& matParams = simulator_.problem().materialLawParams(elemIdx);
            MaterialLaw::capillaryPressures(pc, matParams, fs);
            Valgrind::CheckDefined(this->oilPressure_[elemIdx]);
            Valgrind::CheckDefined(pc);
            assert(FluidSystem::phaseIsActive(oilPhaseIdx));
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;

                fs.setPressure(phaseIdx, this->oilPressure_[elemIdx] + (pc[phaseIdx] - pc[oilPhaseIdx]));
            }
        }

        if (!this->temperature_.empty())
            fs.setTemperature(this->temperature_[elemIdx]);
        if (!this->rs_.empty())
           fs.setRs(this->rs_[elemIdx]);
        if (!this->rv_.empty())
           fs.setRv(this->rv_[elemIdx]);
    }

    void initHysteresisParams(Simulator& simulator, unsigned elemIdx) const
    {
        if (!this->soMax_.empty())
            simulator.problem().setMaxOilSaturation(elemIdx, this->soMax_[elemIdx]);

        if (simulator.problem().materialLawManager()->enableHysteresis()) {
            auto matLawManager = simulator.problem().materialLawManager();

            if (!this->pcSwMdcOw_.empty() && !this->krnSwMdcOw_.empty()) {
                matLawManager->setOilWaterHysteresisParams(
                            this->pcSwMdcOw_[elemIdx],
                            this->krnSwMdcOw_[elemIdx],
                            elemIdx);
            }
            if (!this->pcSwMdcGo_.empty() && !this->krnSwMdcGo_.empty()) {
                matLawManager->setGasOilHysteresisParams(
                            this->pcSwMdcGo_[elemIdx],
                            this->krnSwMdcGo_[elemIdx],
                            elemIdx);
            }
        }

        if (simulator_.vanguard().eclState().fieldProps().has_double("SWATINIT")) {
            auto oilWaterScaledEpsInfoDrainage = simulator.problem().materialLawManager()->oilWaterScaledEpsInfoDrainagePointerReferenceHack(elemIdx);
            oilWaterScaledEpsInfoDrainage->maxPcow =  this->ppcw_[elemIdx];
        }

    }

private:
    bool isDefunctParallelWell(std::string wname) const override
    {
        if (simulator_.gridView().comm().size()==1)
            return false;
        const auto& parallelWells = simulator_.vanguard().parallelWells();
        std::pair<std::string,bool> value{wname, true};
        auto candidate = std::lower_bound(parallelWells.begin(), parallelWells.end(),
                                          value);
        return candidate == parallelWells.end() || *candidate != value;
    }

    void updateFluidInPlace_(const ElementContext& elemCtx, unsigned dofIdx)
    {

        const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0);
        const auto& fs = intQuants.fluidState();
        unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);

        // Fluid in Place calculations

        // calculate the pore volume of the current cell. Note that the porosity
        // returned by the intensive quantities is defined as the ratio of pore
        // space to total cell volume and includes all pressure dependent (->
        // rock compressibility) and static modifiers (MULTPV, MULTREGP, NTG,
        // PORV, MINPV and friends). Also note that because of this, the porosity
        // returned by the intensive quantities can be outside of the physical
        // range [0, 1] in pathetic cases.
        const auto totVolume =
            elemCtx.simulator().model().dofTotalVolume(globalDofIdx);
        const double pv = totVolume * intQuants.porosity().value();

        if (!this->pressureTimesHydrocarbonVolume_.empty() && !this->pressureTimesPoreVolume_.empty()) {
            assert(this->hydrocarbonPoreVolume_.size() ==  this->pressureTimesHydrocarbonVolume_.size());
            assert(this->fip_[Inplace::Phase::PoreVolume].size() == this->pressureTimesPoreVolume_.size());

            this->fip_[Inplace::Phase::PoreVolume][globalDofIdx] =
                totVolume * intQuants.referencePorosity();

            this->dynamicPoreVolume_[globalDofIdx] = pv;

            Scalar hydrocarbon = 0.0;
            if (FluidSystem::phaseIsActive(oilPhaseIdx))
                hydrocarbon += getValue(fs.saturation(oilPhaseIdx));
            if (FluidSystem::phaseIsActive(gasPhaseIdx))
                hydrocarbon += getValue(fs.saturation(gasPhaseIdx));

            this->hydrocarbonPoreVolume_[globalDofIdx] = pv * hydrocarbon;

            if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
                this->pressureTimesPoreVolume_[globalDofIdx] = getValue(fs.pressure(oilPhaseIdx)) * pv;
                this->pressureTimesHydrocarbonVolume_[globalDofIdx] = this->pressureTimesPoreVolume_[globalDofIdx] * hydrocarbon;
            }
        }

        if (this->computeFip_) {
            Scalar fip[FluidSystem::numPhases];
            Scalar fipr[FluidSystem::numPhases]; // at reservoir condition
            for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
                fip[phaseIdx] = 0.0;

                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;

                const double b = getValue(fs.invB(phaseIdx));
                const double s = getValue(fs.saturation(phaseIdx));                
                fipr[phaseIdx] = s * pv;
                fip[phaseIdx] = b * fipr[phaseIdx];
            }

            if (FluidSystem::phaseIsActive(oilPhaseIdx) && !this->fip_[Inplace::Phase::OIL].empty())
                this->fip_[Inplace::Phase::OIL][globalDofIdx] = fip[oilPhaseIdx];
            if (FluidSystem::phaseIsActive(gasPhaseIdx) && !this->fip_[Inplace::Phase::GAS].empty())
                this->fip_[Inplace::Phase::GAS][globalDofIdx] = fip[gasPhaseIdx];
            if (FluidSystem::phaseIsActive(waterPhaseIdx) && !this->fip_[Inplace::Phase::WATER].empty())
                this->fip_[Inplace::Phase::WATER][globalDofIdx] = fip[waterPhaseIdx];

            if (FluidSystem::phaseIsActive(oilPhaseIdx) && !this->fip_[Inplace::Phase::OilResVolume].empty())
                this->fip_[Inplace::Phase::OilResVolume][globalDofIdx] = fipr[oilPhaseIdx];
            if (FluidSystem::phaseIsActive(gasPhaseIdx) && !this->fip_[Inplace::Phase::GasResVolume].empty())
                this->fip_[Inplace::Phase::GasResVolume][globalDofIdx] = fipr[gasPhaseIdx];
            if (FluidSystem::phaseIsActive(waterPhaseIdx) && !this->fip_[Inplace::Phase::WaterResVolume].empty())
                this->fip_[Inplace::Phase::WaterResVolume][globalDofIdx] = fipr[waterPhaseIdx];
            
            if (FluidSystem::phaseIsActive(waterPhaseIdx) && !this->fip_[Inplace::Phase::SALT].empty()) 
                this->fip_[Inplace::Phase::SALT][globalDofIdx] = fipr[waterPhaseIdx] * fs.saltConcentration().value();

            // Store the pure oil and gas Fip
            if (FluidSystem::phaseIsActive(oilPhaseIdx) && !this->fip_[Inplace::Phase::OilInLiquidPhase].empty())
                this->fip_[Inplace::Phase::OilInLiquidPhase][globalDofIdx] = fip[oilPhaseIdx];

            if (FluidSystem::phaseIsActive(gasPhaseIdx) && !this->fip_[Inplace::Phase::GasInGasPhase].empty())
                this->fip_[Inplace::Phase::GasInGasPhase][globalDofIdx] = fip[gasPhaseIdx];

            if (FluidSystem::phaseIsActive(oilPhaseIdx) && FluidSystem::phaseIsActive(gasPhaseIdx)) {
                // Gas dissolved in oil and vaporized oil
                Scalar gasInPlaceLiquid = getValue(fs.Rs()) * fip[oilPhaseIdx];
                Scalar oilInPlaceGas = getValue(fs.Rv()) * fip[gasPhaseIdx];
                if (!this->fip_[Inplace::Phase::GasInLiquidPhase].empty())
                    this->fip_[Inplace::Phase::GasInLiquidPhase][globalDofIdx] = gasInPlaceLiquid;
                if (!this->fip_[Inplace::Phase::OilInGasPhase].empty())
                    this->fip_[Inplace::Phase::OilInGasPhase][globalDofIdx] = oilInPlaceGas;

                // Add dissolved gas and vaporized oil to total Fip
                if (!this->fip_[Inplace::Phase::OIL].empty())
                    this->fip_[Inplace::Phase::OIL][globalDofIdx] += oilInPlaceGas;
                if (!this->fip_[Inplace::Phase::GAS].empty())
                    this->fip_[Inplace::Phase::GAS][globalDofIdx] += gasInPlaceLiquid;
            }
        }

    }

    void createLocalRegion_(std::vector<int>& region)
    {
        ElementContext elemCtx(simulator_);
        ElementIterator elemIt = simulator_.gridView().template begin</*codim=*/0>();
        const ElementIterator& elemEndIt = simulator_.gridView().template end</*codim=*/0>();
        size_t elemIdx = 0;
        for (; elemIt != elemEndIt; ++elemIt, ++elemIdx) {
            const Element& elem = *elemIt;
            if (elem.partitionType() != Dune::InteriorEntity)
                region[elemIdx] = 0;
        }
    }

    const Simulator& simulator_;
};

} // namespace Opm

#endif
