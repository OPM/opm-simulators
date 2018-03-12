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
 * \copydoc Ewoms::EclOutputBlackOilModule
 */
#ifndef EWOMS_ECL_OUTPUT_BLACK_OIL_MODULE_HH
#define EWOMS_ECL_OUTPUT_BLACK_OIL_MODULE_HH

#include "eclwriter.hh"

#include <ewoms/models/blackoil/blackoilproperties.hh>

#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>

#include <opm/material/common/Valgrind.hpp>

#if HAVE_ECL_INPUT
#include <opm/parser/eclipse/Units/Units.hpp>
#include <opm/parser/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#endif

#if HAVE_ECL_OUTPUT
#include <opm/output/data/Cells.hpp>
#include <opm/output/eclipse/EclipseIO.hpp>
#endif

#include <opm/common/OpmLog/OpmLog.hpp>

#include <dune/common/fvector.hh>

#include <type_traits>

namespace Ewoms {
namespace Properties {
// create new type tag for the Ecl-output
NEW_TYPE_TAG(EclOutputBlackOil);

} // namespace Properties

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
class EclOutputBlackOilModule
{

    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Discretization) Discretization;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    enum { numPhases = FluidSystem::numPhases };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { gasCompIdx = FluidSystem::gasCompIdx };
    enum { oilCompIdx = FluidSystem::oilCompIdx };

    typedef std::vector<Scalar> ScalarBuffer;

    struct FipDataType {

        enum FipId {
            WaterInPlace = 0, //WIP
            OilInPlace = 1, //OIP
            GasInPlace = 2, //GIP
            OilInPlaceInLiquidPhase = 3, //OIPL
            OilInPlaceInGasPhase = 4, //OIPG
            GasInPlaceInLiquidPhase = 5, //GIPL
            GasInPlaceInGasPhase = 6,  //GIPG
            PoreVolume   = 7, //PV
        };
        static const int numFipValues = PoreVolume + 1 ;
    };

public:
    template<class CollectDataToIORankType>
    EclOutputBlackOilModule(const Simulator& simulator, const CollectDataToIORankType& collectToIORank)
        : simulator_(simulator)
    {
        createLocalFipnum_();

        // Summary output is for all steps
        const Opm::SummaryConfig summaryConfig = simulator_.vanguard().summaryConfig();

        // Initialize block output
        for( const auto& node : summaryConfig ) {
            if (node.type() == ECL_SMSPEC_BLOCK_VAR) {
                if(collectToIORank.isGlobalIdxOnThisRank(node.num() - 1)) {
                    std::pair<std::string, int> key = std::make_pair(node.keyword(), node.num());
                    blockData_[key] = 0.0;
                }
            }
        }
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to ECL output files
     */
    void allocBuffers(unsigned bufferSize, unsigned reportStepNum, const bool substep, const bool log)
    {

        if (!std::is_same<Discretization, Ewoms::EcfvDiscretization<TypeTag> >::value)
            return;

        // Summary output is for all steps
        const Opm::SummaryConfig summaryConfig = simulator_.vanguard().summaryConfig();

        // Only output RESTART_AUXILIARY asked for by the user.
        const Opm::RestartConfig& restartConfig = simulator_.vanguard().eclState().getRestartConfig();
        std::map<std::string, int> rstKeywords = restartConfig.getRestartKeywords(reportStepNum);
        for (auto& keyValue : rstKeywords) {
            keyValue.second = restartConfig.getKeyword(keyValue.first, reportStepNum);
        }

        outputFipRestart_ = false;
        computeFip_ = false;

        // Fluid in place
        for (int i = 0; i<FipDataType::numFipValues; i++) {
            if (!substep || summaryConfig.require3DField(fipEnumToString_(i))) {
                if (rstKeywords["FIP"] > 0) {
                    rstKeywords["FIP"] = 0;
                    outputFipRestart_ = true;
                }
                fip_[i].resize(bufferSize, 0.0);
                computeFip_ = true;
            } else {
                fip_[i].clear();
            }
        }
        if (!substep || summaryConfig.hasKeyword("FPR") || summaryConfig.hasKeyword("FPRP") || summaryConfig.hasKeyword("RPR")) {
            fip_[FipDataType::PoreVolume].resize(bufferSize, 0.0);
            hydrocarbonPoreVolume_.resize(bufferSize, 0.0);
            pressureTimesPoreVolume_.resize(bufferSize, 0.0);
            pressureTimesHydrocarbonVolume_.resize(bufferSize, 0.0);
        } else {
            hydrocarbonPoreVolume_.clear();
            pressureTimesPoreVolume_.clear();
            pressureTimesHydrocarbonVolume_.clear();
        }

        // Well RFT data
        if (!substep) {
            for ( const auto& well : simulator_.vanguard().schedule().getWells( reportStepNum )) {

                // don't bother with wells not on this process
                const auto& defunct_well_names = simulator_.vanguard().defunctWellNames();
                if ( defunct_well_names.find(well->name()) != defunct_well_names.end() ) {
                    continue;
                }

                if( !( well->getRFTActive( reportStepNum )
                       || well->getPLTActive( reportStepNum ) ) )
                    continue;

                for( const auto& completion : well->getCompletions( reportStepNum ) ) {
                    const size_t i = size_t( completion.getI() );
                    const size_t j = size_t( completion.getJ() );
                    const size_t k = size_t( completion.getK() );
                    const size_t index = simulator_.vanguard().eclState().getInputGrid().getGlobalIndex( i, j, k );

                    oilCompletionPressures_.emplace(std::make_pair(index, 0.0));
                    waterCompletionSaturations_.emplace(std::make_pair(index, 0.0));
                    gasCompletionSaturations_.emplace(std::make_pair(index, 0.0));
                }
            }
        }

        // Only provide restart on restart steps
        if (!restartConfig.getWriteRestartFile(reportStepNum) || substep)
            return;

        // always output saturation of active phases
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx))
                continue;

            saturation_[phaseIdx].resize(bufferSize,0.0);
        }
        // and oil pressure
        oilPressure_.resize(bufferSize,0.0);

        if (true) // adding some logic here when the energy module is done
            temperature_.resize(bufferSize,0.0);

        if (FluidSystem::enableDissolvedGas()) {
            rs_.resize(bufferSize,0.0);
        }
        if (FluidSystem::enableVaporizedOil()) {
            rv_.resize(bufferSize,0.0);
        }

        if (GET_PROP_VALUE(TypeTag, EnableSolvent)) {
            sSol_.resize(bufferSize,0.0);
        }
        if (GET_PROP_VALUE(TypeTag, EnablePolymer)) {
            cPolymer_.resize(bufferSize,0.0);
        }

        if (simulator_.problem().vapparsActive())
            soMax_.resize(bufferSize,0.0);

        if (simulator_.problem().materialLawManager()->enableHysteresis()) {
            pcSwMdcOw_.resize(bufferSize,0.0);
            krnSwMdcOw_.resize(bufferSize,0.0);
            pcSwMdcGo_.resize(bufferSize,0.0);
            krnSwMdcGo_.resize(bufferSize,0.0);
        }

        if (FluidSystem::enableDissolvedGas() && rstKeywords["RSSAT"] > 0) {
            rstKeywords["RSSAT"] = 0;
            gasDissolutionFactor_.resize(bufferSize,0.0);
        }
        if (FluidSystem::enableVaporizedOil() && rstKeywords["RVSAT"] > 0) {
            rstKeywords["RVSAT"] = 0;
            oilVaporizationFactor_.resize(bufferSize,0.0);
        }

        if (FluidSystem::phaseIsActive(waterPhaseIdx) && rstKeywords["BW"] > 0)
        {
            rstKeywords["BW"] = 0;
            invB_[waterPhaseIdx].resize(bufferSize,0.0);
        }
        if (FluidSystem::phaseIsActive(oilPhaseIdx) && rstKeywords["BO"] > 0)
        {
            rstKeywords["BO"] = 0;
            invB_[oilPhaseIdx].resize(bufferSize,0.0);
        }
        if (FluidSystem::phaseIsActive(gasPhaseIdx) && rstKeywords["BG"] > 0)
        {
            rstKeywords["BG"] = 0;
            invB_[gasPhaseIdx].resize(bufferSize,0.0);
        }

        if (rstKeywords["DEN"] > 0) {
            rstKeywords["DEN"] = 0;
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;
                density_[phaseIdx].resize(bufferSize,0.0);
            }
        }
        const bool hasVWAT = (rstKeywords["VISC"] > 0) || (rstKeywords["VWAT"] > 0);
        const bool hasVOIL = (rstKeywords["VISC"] > 0) || (rstKeywords["VOIL"] > 0);
        const bool hasVGAS = (rstKeywords["VISC"] > 0) || (rstKeywords["VGAS"] > 0);
        rstKeywords["VISC"] = 0;

        if (FluidSystem::phaseIsActive(waterPhaseIdx) && hasVWAT)
        {
            rstKeywords["VWAT"] = 0;
            viscosity_[waterPhaseIdx].resize(bufferSize,0.0);
        }
        if (FluidSystem::phaseIsActive(oilPhaseIdx) && hasVOIL > 0)
        {
            rstKeywords["VOIL"] = 0;
            viscosity_[oilPhaseIdx].resize(bufferSize,0.0);
        }
        if (FluidSystem::phaseIsActive(gasPhaseIdx) && hasVGAS > 0)
        {
            rstKeywords["VGAS"] = 0;
            viscosity_[gasPhaseIdx].resize(bufferSize,0.0);
        }

        if (FluidSystem::phaseIsActive(waterPhaseIdx) && rstKeywords["KRW"] > 0)
        {
            rstKeywords["KRW"] = 0;
            relativePermeability_[waterPhaseIdx].resize(bufferSize,0.0);
        }
        if (FluidSystem::phaseIsActive(oilPhaseIdx) && rstKeywords["KRO"] > 0)
        {
            rstKeywords["KRO"] = 0;
            relativePermeability_[oilPhaseIdx].resize(bufferSize,0.0);
        }
        if (FluidSystem::phaseIsActive(gasPhaseIdx) && rstKeywords["KRG"] > 0)
        {
            rstKeywords["KRG"] = 0;
            relativePermeability_[gasPhaseIdx].resize(bufferSize,0.0);
        }

        if (rstKeywords["PBPD"] > 0)  {
            rstKeywords["PBPD"] = 0;
            bubblePointPressure_.resize(bufferSize,0.0);
            dewPointPressure_.resize(bufferSize,0.0);
        }

        //Warn for any unhandled keyword
        if (log) {
            for (auto& keyValue : rstKeywords) {
                if (keyValue.second > 0) {
                    std::string logstring = "Keyword '";
                    logstring.append(keyValue.first);
                    logstring.append("' is unhandled for output to file.");
                    Opm::OpmLog::warning("Unhandled output keyword", logstring);
                }
            }
        }

        failedCellsPb_.clear();
        failedCellsPd_.clear();

        // Not supported in flow legacy
        if (false)
            saturatedOilFormationVolumeFactor_.resize(bufferSize,0.0);
        if (false)
            oilSaturationPressure_.resize(bufferSize,0.0);
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quanties relevant
     *        for an element
     */
    void processElement(const ElementContext& elemCtx)
    {

        if (!std::is_same<Discretization, Ewoms::EcfvDiscretization<TypeTag> >::value)
            return;

        for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
            const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0);
            const auto& fs = intQuants.fluidState();

            typedef typename std::remove_const<typename std::remove_reference<decltype(fs)>::type>::type FluidState;
            unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
            unsigned pvtRegionIdx = elemCtx.primaryVars(dofIdx, /*timeIdx=*/0).pvtRegionIndex();

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                if (saturation_[phaseIdx].size() == 0)
                    continue;

                saturation_[phaseIdx][globalDofIdx] = Opm::getValue(fs.saturation(phaseIdx));
                Opm::Valgrind::CheckDefined(saturation_[phaseIdx][globalDofIdx]);
            }

            if (oilPressure_.size() > 0) {
                oilPressure_[globalDofIdx] = Opm::getValue(fs.pressure(oilPhaseIdx));
                Opm::Valgrind::CheckDefined(oilPressure_[globalDofIdx]);
            }

            if (temperature_.size() > 0) {
                temperature_[globalDofIdx] = Opm::getValue(fs.temperature(oilPhaseIdx));
                Opm::Valgrind::CheckDefined(temperature_[globalDofIdx]);
            }
            if (gasDissolutionFactor_.size() > 0) {
                Scalar SoMax = elemCtx.problem().maxOilSaturation(globalDofIdx);
                gasDissolutionFactor_[globalDofIdx] =
                    FluidSystem::template saturatedDissolutionFactor<FluidState, Scalar>(fs, oilPhaseIdx, pvtRegionIdx, SoMax);
                Opm::Valgrind::CheckDefined(gasDissolutionFactor_[globalDofIdx]);

            }
            if (oilVaporizationFactor_.size() > 0) {
                Scalar SoMax = elemCtx.problem().maxOilSaturation(globalDofIdx);
                oilVaporizationFactor_[globalDofIdx] =
                    FluidSystem::template saturatedDissolutionFactor<FluidState, Scalar>(fs, gasPhaseIdx, pvtRegionIdx, SoMax);
                Opm::Valgrind::CheckDefined(oilVaporizationFactor_[globalDofIdx]);

            }
            if (gasFormationVolumeFactor_.size() > 0) {
                gasFormationVolumeFactor_[globalDofIdx] =
                    1.0/FluidSystem::template inverseFormationVolumeFactor<FluidState, Scalar>(fs, gasPhaseIdx, pvtRegionIdx);
                Opm::Valgrind::CheckDefined(gasFormationVolumeFactor_[globalDofIdx]);

            }
            if (saturatedOilFormationVolumeFactor_.size() > 0) {
                saturatedOilFormationVolumeFactor_[globalDofIdx] =
                    1.0/FluidSystem::template saturatedInverseFormationVolumeFactor<FluidState, Scalar>(fs, oilPhaseIdx, pvtRegionIdx);
                Opm::Valgrind::CheckDefined(saturatedOilFormationVolumeFactor_[globalDofIdx]);

            }
            if (oilSaturationPressure_.size() > 0) {
                oilSaturationPressure_[globalDofIdx] =
                    FluidSystem::template saturationPressure<FluidState, Scalar>(fs, oilPhaseIdx, pvtRegionIdx);
                Opm::Valgrind::CheckDefined(oilSaturationPressure_[globalDofIdx]);

            }

            if (rs_.size()) {
                rs_[globalDofIdx] = Opm::getValue(fs.Rs());
                Opm::Valgrind::CheckDefined(rs_[globalDofIdx]);
            }

            if (rv_.size()) {
                rv_[globalDofIdx] = Opm::getValue(fs.Rv());
                Opm::Valgrind::CheckDefined(rv_[globalDofIdx]);
            }

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                if (invB_[phaseIdx].size() == 0)
                    continue;

                invB_[phaseIdx][globalDofIdx] = Opm::getValue(fs.invB(phaseIdx));
                Opm::Valgrind::CheckDefined(invB_[phaseIdx][globalDofIdx]);
            }

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                if (density_[phaseIdx].size() == 0)
                    continue;

                density_[phaseIdx][globalDofIdx] = Opm::getValue(fs.density(phaseIdx));
                Opm::Valgrind::CheckDefined(density_[phaseIdx][globalDofIdx]);
            }

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                if (viscosity_[phaseIdx].size() == 0)
                    continue;

                viscosity_[phaseIdx][globalDofIdx] = Opm::getValue(fs.viscosity(phaseIdx));
                Opm::Valgrind::CheckDefined(viscosity_[phaseIdx][globalDofIdx]);
            }

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                if (relativePermeability_[phaseIdx].size() == 0)
                    continue;

                relativePermeability_[phaseIdx][globalDofIdx] = Opm::getValue(intQuants.relativePermeability(phaseIdx));
                Opm::Valgrind::CheckDefined(relativePermeability_[phaseIdx][globalDofIdx]);
            }

            if (sSol_.size() > 0) {
                sSol_[globalDofIdx] = intQuants.solventSaturation().value();
            }

            if (cPolymer_.size() > 0) {
                cPolymer_[globalDofIdx] = intQuants.polymerConcentration().value();
            }

            if (bubblePointPressure_.size() > 0)
            {
                try {
                    bubblePointPressure_[globalDofIdx] = Opm::getValue(FluidSystem::bubblePointPressure(fs, intQuants.pvtRegionIndex()));
                }
                catch (const Opm::NumericalIssue&) {
                    const auto globalIdx = elemCtx.simulator().vanguard().grid().globalCell()[globalDofIdx];
                    failedCellsPb_.push_back(globalIdx);
                }
            }
            if (dewPointPressure_.size() > 0)
            {
                try {
                    dewPointPressure_[globalDofIdx] = Opm::getValue(FluidSystem::dewPointPressure(fs, intQuants.pvtRegionIndex()));
                }
                catch (const Opm::NumericalIssue&) {
                    const auto globalIdx = elemCtx.simulator().vanguard().grid().globalCell()[globalDofIdx];
                    failedCellsPd_.push_back(globalIdx);
                }
            }

            if (soMax_.size() > 0)
                soMax_[globalDofIdx] = elemCtx.simulator().problem().maxOilSaturation(globalDofIdx);

            const auto& matLawManager = elemCtx.simulator().problem().materialLawManager();
            if (matLawManager->enableHysteresis()) {
                if (pcSwMdcOw_.size() > 0 && krnSwMdcOw_.size() > 0) {
                    matLawManager->oilWaterHysteresisParams(
                                pcSwMdcOw_[globalDofIdx],
                                krnSwMdcOw_[globalDofIdx],
                                globalDofIdx);
                }
                if (pcSwMdcGo_.size() > 0 && krnSwMdcGo_.size() > 0) {
                    matLawManager->gasOilHysteresisParams(
                                pcSwMdcGo_[globalDofIdx],
                                krnSwMdcGo_[globalDofIdx],
                                globalDofIdx);
                }
            }

            // hack to make the intial output of rs and rv Ecl compatible.
            // For cells with swat == 1 Ecl outputs; rs = rsSat and rv=rvSat, in all but the initial step
            // where it outputs rs and rv values calculated by the initialization. To be compatible we overwrite
            // rs and rv with the values computed in the initially.
            // Volume factors, densities and viscosities need to be recalculated with the updated rs and rv values.
            // This can be removed when ebos has 100% controll over output
            if (elemCtx.simulator().episodeIndex() < 0 && FluidSystem::phaseIsActive(oilPhaseIdx) && FluidSystem::phaseIsActive(gasPhaseIdx) ) {

                const auto& fs_initial = elemCtx.simulator().problem().initialFluidState(globalDofIdx);

                // use initial rs and rv values
                if (rv_.size() > 0)
                    rv_[globalDofIdx] = fs_initial.Rv();

                if (rs_.size() > 0)
                    rs_[globalDofIdx] = fs_initial.Rs();

                // re-compute the volume factors, viscosities and densities if asked for
                if (density_[oilPhaseIdx].size() > 0)
                    density_[oilPhaseIdx][globalDofIdx] = FluidSystem::density(fs_initial,
                                                                                    oilPhaseIdx,
                                                                                    intQuants.pvtRegionIndex());
                if (density_[gasPhaseIdx].size() > 0)
                    density_[gasPhaseIdx][globalDofIdx] = FluidSystem::density(fs_initial,
                                                                                    gasPhaseIdx,
                                                                                    intQuants.pvtRegionIndex());

                if (invB_[oilPhaseIdx].size() > 0)
                    invB_[oilPhaseIdx][globalDofIdx] = FluidSystem::inverseFormationVolumeFactor(fs_initial,
                                                                                                      oilPhaseIdx,
                                                                                                      intQuants.pvtRegionIndex());
                if (invB_[gasPhaseIdx].size() > 0)
                    invB_[gasPhaseIdx][globalDofIdx] = FluidSystem::inverseFormationVolumeFactor(fs_initial,
                                                                                                      gasPhaseIdx,
                                                                                                      intQuants.pvtRegionIndex());
                if (viscosity_[oilPhaseIdx].size() > 0)
                    viscosity_[oilPhaseIdx][globalDofIdx] = FluidSystem::viscosity(fs_initial,
                                                                                        oilPhaseIdx,
                                                                                        intQuants.pvtRegionIndex());
                if (viscosity_[gasPhaseIdx].size() > 0)
                    viscosity_[gasPhaseIdx][globalDofIdx] = FluidSystem::viscosity(fs_initial,
                                                                                        gasPhaseIdx,
                                                                                        intQuants.pvtRegionIndex());
            }

            // Add fluid in Place values
            updateFluidInPlace_(elemCtx, dofIdx);

            // Adding block data
            const auto globalIdx = elemCtx.simulator().vanguard().grid().globalCell()[globalDofIdx];
            for( auto& val : blockData_ ) {
                const auto& key = val.first;
                int global_index = key.second - 1;
                if (global_index == globalIdx) {
                    if (key.first == "BWSAT") {
                        val.second = Opm::getValue(fs.saturation(waterPhaseIdx));
                    } else if (key.first == "BGSAT") {
                        val.second = Opm::getValue(fs.saturation(gasPhaseIdx));
                    } else if (key.first == "BPR") {
                        val.second = Opm::getValue(fs.pressure(oilPhaseIdx));
                    } else {
                        std::string logstring = "Keyword '";
                        logstring.append(key.first);
                        logstring.append("' is unhandled for output to file.");
                        Opm::OpmLog::warning("Unhandled output keyword", logstring);
                    }
                }
            }

            // Adding Well RFT data
            if (oilCompletionPressures_.count(globalIdx) > 0) {
                oilCompletionPressures_[globalIdx] = Opm::getValue(fs.pressure(oilPhaseIdx));
            }
            if (waterCompletionSaturations_.count(globalIdx) > 0) {
                waterCompletionSaturations_[globalIdx] = Opm::getValue(fs.saturation(waterPhaseIdx));
            }
            if (gasCompletionSaturations_.count(globalIdx) > 0) {
                gasCompletionSaturations_[globalIdx] = Opm::getValue(fs.saturation(gasPhaseIdx));
            }
        }
    }


    void outputErrorLog()
    {
        const size_t maxNumCellsFaillog = 20;

        int pbSize = failedCellsPb_.size(), pd_size = failedCellsPd_.size();
        std::vector<int> displPb, displPd, recvLenPb, recvLenPd;
        const auto& comm = simulator_.gridView().comm();

        if ( isIORank_() )
        {
            displPb.resize(comm.size()+1, 0);
            displPd.resize(comm.size()+1, 0);
            recvLenPb.resize(comm.size());
            recvLenPd.resize(comm.size());
        }

        comm.gather(&pbSize, recvLenPb.data(), 1, 0);
        comm.gather(&pd_size, recvLenPd.data(), 1, 0);
        std::partial_sum(recvLenPb.begin(), recvLenPb.end(), displPb.begin()+1);
        std::partial_sum(recvLenPd.begin(), recvLenPd.end(), displPd.begin()+1);
        std::vector<int> globalFailedCellsPb, globalFailedCellsPd;

        if ( isIORank_() )
        {
            globalFailedCellsPb.resize(displPb.back());
            globalFailedCellsPd.resize(displPd.back());
        }

        comm.gatherv(failedCellsPb_.data(), static_cast<int>(failedCellsPb_.size()),
                     globalFailedCellsPb.data(), recvLenPb.data(),
                     displPb.data(), 0);
        comm.gatherv(failedCellsPd_.data(), static_cast<int>(failedCellsPd_.size()),
                     globalFailedCellsPd.data(),  recvLenPd.data(),
                     displPd.data(), 0);
        std::sort(globalFailedCellsPb.begin(), globalFailedCellsPb.end());
        std::sort(globalFailedCellsPd.begin(), globalFailedCellsPd.end());

        if (globalFailedCellsPb.size() > 0) {
            std::stringstream errlog;
            errlog << "Finding the bubble point pressure failed for " << globalFailedCellsPb.size() << " cells [";
            errlog << globalFailedCellsPb[0];
            const size_t max_elems = std::min(maxNumCellsFaillog, globalFailedCellsPb.size());
            for (size_t i = 1; i < max_elems; ++i) {
                errlog << ", " << globalFailedCellsPb[i];
            }
            if (globalFailedCellsPb.size() > maxNumCellsFaillog) {
                errlog << ", ...";
            }
            errlog << "]";
            Opm::OpmLog::warning("Bubble point numerical problem", errlog.str());
        }
        if (globalFailedCellsPd.size() > 0) {
            std::stringstream errlog;
            errlog << "Finding the dew point pressure failed for " << globalFailedCellsPd.size() << " cells [";
            errlog << globalFailedCellsPd[0];
            const size_t max_elems = std::min(maxNumCellsFaillog, globalFailedCellsPd.size());
            for (size_t i = 1; i < max_elems; ++i) {
                errlog << ", " << globalFailedCellsPd[i];
            }
            if (globalFailedCellsPd.size() > maxNumCellsFaillog) {
                errlog << ", ...";
            }
            errlog << "]";
            Opm::OpmLog::warning("Dew point numerical problem", errlog.str());
        }
    }

    void addRftDataToWells(Opm::data::Wells& wellDatas, size_t reportStepNum)
    {
        for ( const auto& well : simulator_.vanguard().schedule().getWells( reportStepNum )) {

            // don't bother with wells not on this process
            const auto& defunct_well_names = simulator_.vanguard().defunctWellNames();
            if ( defunct_well_names.find(well->name()) != defunct_well_names.end() ) {
                continue;
            }

            //add data infrastructure for shut wells
            if (!wellDatas.count(well->name())){
                Opm::data::Well wellData;

                if( !( well->getRFTActive( reportStepNum )
                       || well->getPLTActive( reportStepNum ) ) )
                    continue;
                wellData.completions.resize(well->getCompletions( reportStepNum ).size());
                size_t count = 0;
                for( const auto& completion : well->getCompletions( reportStepNum ) ) {

                    const size_t i = size_t( completion.getI() );
                    const size_t j = size_t( completion.getJ() );
                    const size_t k = size_t( completion.getK() );

                    const size_t index = simulator_.vanguard().eclState().getInputGrid().getGlobalIndex( i, j, k );
                    auto& completionData = wellData.completions[ count ];
                    completionData.index = index;
                    count++;
                }
                wellDatas.emplace( std::make_pair(well->name(),wellData) );
            }

            Opm::data::Well& wellData = wellDatas.at(well->name());
            for (auto& completionData : wellData.completions) {
                const auto index = completionData.index;
                if (oilCompletionPressures_.count(index) > 0)
                    completionData.cell_pressure = oilCompletionPressures_.at(index);
                if (waterCompletionSaturations_.count(index) > 0)
                    completionData.cell_saturation_water = waterCompletionSaturations_.at(index);
                if (gasCompletionSaturations_.count(index) > 0)
                    completionData.cell_saturation_gas = gasCompletionSaturations_.at(index);
            }
        }
        oilCompletionPressures_.clear();
        waterCompletionSaturations_.clear();
        gasCompletionSaturations_.clear();
    }

    /*!
     * \brief Move all buffers to data::Solution.
     */
    void assignToSolution(Opm::data::Solution& sol)
    {
        if (!std::is_same<Discretization, Ewoms::EcfvDiscretization<TypeTag> >::value)
            return;

        if ( oilPressure_.size() > 0 ) {
            sol.insert( "PRESSURE", Opm::UnitSystem::measure::pressure, std::move(oilPressure_), Opm::data::TargetType::RESTART_SOLUTION);
        }

        if ( temperature_.size() > 0 ) {
            sol.insert( "TEMP", Opm::UnitSystem::measure::temperature, std::move(temperature_), Opm::data::TargetType::RESTART_SOLUTION);
        }

        if( FluidSystem::phaseIsActive(waterPhaseIdx) && saturation_[waterPhaseIdx].size() > 0 ) {
            sol.insert( "SWAT", Opm::UnitSystem::measure::identity, std::move(saturation_[waterPhaseIdx]), Opm::data::TargetType::RESTART_SOLUTION );
        }
        if( FluidSystem::phaseIsActive(gasPhaseIdx) && saturation_[gasPhaseIdx].size() > 0) {
            sol.insert( "SGAS", Opm::UnitSystem::measure::identity, std::move(saturation_[gasPhaseIdx]), Opm::data::TargetType::RESTART_SOLUTION );
        }

        if ( gasDissolutionFactor_.size() > 0 ) {
            sol.insert( "RSSAT", Opm::UnitSystem::measure::gas_oil_ratio, std::move(gasDissolutionFactor_), Opm::data::TargetType::RESTART_AUXILIARY );

        }
        if ( oilVaporizationFactor_.size() > 0 ) {
            sol.insert( "RVSAT", Opm::UnitSystem::measure::oil_gas_ratio, std::move(oilVaporizationFactor_) , Opm::data::TargetType::RESTART_AUXILIARY );
        }
        if ( rs_.size() > 0 ) {
            sol.insert( "RS", Opm::UnitSystem::measure::gas_oil_ratio, std::move(rs_), Opm::data::TargetType::RESTART_SOLUTION );

        }
        if (rv_.size() > 0 ) {
            sol.insert( "RV", Opm::UnitSystem::measure::oil_gas_ratio, std::move(rv_) , Opm::data::TargetType::RESTART_SOLUTION );
        }
        if (invB_[waterPhaseIdx].size() > 0 ) {
            sol.insert( "1OVERBW", Opm::UnitSystem::measure::water_inverse_formation_volume_factor, std::move(invB_[waterPhaseIdx]), Opm::data::TargetType::RESTART_AUXILIARY );
        }
        if (invB_[oilPhaseIdx].size() > 0 ) {
            sol.insert( "1OVERBO", Opm::UnitSystem::measure::oil_inverse_formation_volume_factor, std::move(invB_[oilPhaseIdx]), Opm::data::TargetType::RESTART_AUXILIARY );
        }
        if (invB_[gasPhaseIdx].size() > 0 ) {
            sol.insert( "1OVERBG", Opm::UnitSystem::measure::gas_inverse_formation_volume_factor, std::move(invB_[gasPhaseIdx]), Opm::data::TargetType::RESTART_AUXILIARY );
        }

        if (density_[waterPhaseIdx].size() > 0 ) {
            sol.insert( "WAT_DEN", Opm::UnitSystem::measure::density, std::move(density_[waterPhaseIdx]), Opm::data::TargetType::RESTART_AUXILIARY );
        }
        if (density_[oilPhaseIdx].size() > 0 ) {
            sol.insert( "OIL_DEN", Opm::UnitSystem::measure::density, std::move(density_[oilPhaseIdx]), Opm::data::TargetType::RESTART_AUXILIARY );
        }
        if (density_[gasPhaseIdx].size() > 0 ) {
            sol.insert( "GAS_DEN", Opm::UnitSystem::measure::density, std::move(density_[gasPhaseIdx]), Opm::data::TargetType::RESTART_AUXILIARY );
        }

        if (viscosity_[waterPhaseIdx].size() > 0 ) {
            sol.insert( "WAT_VISC", Opm::UnitSystem::measure::viscosity, std::move(viscosity_[waterPhaseIdx]), Opm::data::TargetType::RESTART_AUXILIARY);
        }
        if (viscosity_[oilPhaseIdx].size() > 0 ) {
            sol.insert( "OIL_VISC", Opm::UnitSystem::measure::viscosity, std::move(viscosity_[oilPhaseIdx]), Opm::data::TargetType::RESTART_AUXILIARY);
        }
        if (viscosity_[gasPhaseIdx].size() > 0 ) {
            sol.insert( "GAS_VISC", Opm::UnitSystem::measure::viscosity, std::move(viscosity_[gasPhaseIdx]), Opm::data::TargetType::RESTART_AUXILIARY );
        }

        if (relativePermeability_[waterPhaseIdx].size() > 0) {
            sol.insert( "WATKR", Opm::UnitSystem::measure::identity, std::move(relativePermeability_[waterPhaseIdx]), Opm::data::TargetType::RESTART_AUXILIARY);
        }
        if (relativePermeability_[oilPhaseIdx].size() > 0 ) {
            sol.insert( "OILKR", Opm::UnitSystem::measure::identity, std::move(relativePermeability_[oilPhaseIdx]), Opm::data::TargetType::RESTART_AUXILIARY);
        }
        if (relativePermeability_[gasPhaseIdx].size() > 0 ) {
            sol.insert( "GASKR", Opm::UnitSystem::measure::identity, std::move(relativePermeability_[gasPhaseIdx]), Opm::data::TargetType::RESTART_AUXILIARY);
        }

        if (pcSwMdcOw_.size() > 0 )
            sol.insert ("PCSWM_OW", Opm::UnitSystem::measure::identity, std::move(pcSwMdcOw_), Opm::data::TargetType::RESTART_AUXILIARY);

        if (krnSwMdcOw_.size() > 0)
            sol.insert ("KRNSW_OW", Opm::UnitSystem::measure::identity, std::move(krnSwMdcOw_), Opm::data::TargetType::RESTART_AUXILIARY);

        if (pcSwMdcGo_.size() > 0)
            sol.insert ("PCSWM_GO", Opm::UnitSystem::measure::identity, std::move(pcSwMdcGo_), Opm::data::TargetType::RESTART_AUXILIARY);

        if (krnSwMdcGo_.size() > 0)
            sol.insert ("KRNSW_GO", Opm::UnitSystem::measure::identity, std::move(krnSwMdcGo_), Opm::data::TargetType::RESTART_AUXILIARY);

        if (soMax_.size() > 0)
            sol.insert ("SOMAX", Opm::UnitSystem::measure::identity, std::move(soMax_), Opm::data::TargetType::RESTART_SOLUTION);

        if (sSol_.size() > 0)
            sol.insert ("SSOL", Opm::UnitSystem::measure::identity, std::move(sSol_), Opm::data::TargetType::RESTART_SOLUTION);

        if (cPolymer_.size() > 0)
            sol.insert ("POLYMER", Opm::UnitSystem::measure::identity, std::move(cPolymer_), Opm::data::TargetType::RESTART_SOLUTION);

        if (dewPointPressure_.size() > 0)
            sol.insert ("PDEW", Opm::UnitSystem::measure::pressure, std::move(dewPointPressure_), Opm::data::TargetType::RESTART_AUXILIARY);

        if (bubblePointPressure_.size() > 0)
            sol.insert ("PBUB", Opm::UnitSystem::measure::pressure, std::move(bubblePointPressure_), Opm::data::TargetType::RESTART_AUXILIARY);

        // Fluid in place
        for (int i = 0; i<FipDataType::numFipValues; i++) {
            if (outputFipRestart_ && fip_[i].size() > 0) {
                sol.insert(fipEnumToString_(i),
                           Opm::UnitSystem::measure::volume,
                           fip_[i] ,
                           Opm::data::TargetType::SUMMARY);
            }
        }
    }

    // write Fluid In Place to output log
    void outputFipLog(std::map<std::string, double>& miscSummaryData,  std::map<std::string, std::vector<double>>& regionData, const bool substep) {

        const auto& comm = simulator_.gridView().comm();
        size_t ntFip = *std::max_element(fipnum_.begin(), fipnum_.end());
        ntFip = comm.max(ntFip);

        // sum values over each region
        ScalarBuffer regionFipValues[FipDataType::numFipValues];
        for (int i = 0; i<FipDataType::numFipValues; i++) {
            regionFipValues[i] = computeFipForRegions_(fip_[i], fipnum_, ntFip);
            if (isIORank_() && origRegionValues_[i].empty())
                origRegionValues_[i] = regionFipValues[i];
        }

        // sum all region values to compute the field total
        std::vector<int> fieldNum(ntFip, 1);
        ScalarBuffer fieldFipValues(FipDataType::numFipValues,0.0);
        bool comunicateSum = false; // the regionValues are already summed over all ranks.
        for (int i = 0; i<FipDataType::numFipValues; i++) {
            const ScalarBuffer& tmp = computeFipForRegions_(regionFipValues[i], fieldNum, 1, comunicateSum);
            fieldFipValues[i] = tmp[0];
        }

        // compute the hydrocarbon averaged pressure over the regions.
        ScalarBuffer regPressurePv = computeFipForRegions_(pressureTimesPoreVolume_, fipnum_, ntFip);
        ScalarBuffer regPvHydrocarbon = computeFipForRegions_(hydrocarbonPoreVolume_, fipnum_, ntFip);
        ScalarBuffer regPressurePvHydrocarbon = computeFipForRegions_(pressureTimesHydrocarbonVolume_, fipnum_, ntFip);

        ScalarBuffer fieldPressurePv = computeFipForRegions_(regPressurePv, fieldNum, 1, comunicateSum);
        ScalarBuffer fieldPvHydrocarbon = computeFipForRegions_(regPvHydrocarbon, fieldNum, 1, comunicateSum);
        ScalarBuffer fieldPressurePvHydrocarbon = computeFipForRegions_(regPressurePvHydrocarbon, fieldNum, 1, comunicateSum);

        // output on io rank
        // the original Fip values are stored on the first step
        // TODO: Store initial Fip in the init file and restore them
        // and use them here.
        const Opm::SummaryConfig summaryConfig = simulator_.vanguard().summaryConfig();
        if ( isIORank_()) {
            // Field summary output
            for (int i = 0; i<FipDataType::numFipValues; i++) {
                std::string key = "F" + fipEnumToString_(i);
                if (summaryConfig.hasKeyword(key))
                    miscSummaryData[key] = fieldFipValues[i];
            }
            if (summaryConfig.hasKeyword("FOE") && !origTotalValues_.empty())
                miscSummaryData["FOE"] = fieldFipValues[FipDataType::OilInPlace] / origTotalValues_[FipDataType::OilInPlace];

            if (summaryConfig.hasKeyword("FPR"))
                miscSummaryData["FPR"] = pressureAverage_(fieldPressurePvHydrocarbon[0], fieldPvHydrocarbon[0], fieldPressurePv[0], fieldFipValues[FipDataType::PoreVolume], true);

            if (summaryConfig.hasKeyword("FPRP"))
                miscSummaryData["FPR"] = pressureAverage_(fieldPressurePvHydrocarbon[0], fieldPvHydrocarbon[0], fieldPressurePv[0], fieldFipValues[FipDataType::PoreVolume], false);

            // Region summary output
            for (int i = 0; i<FipDataType::numFipValues; i++) {
                std::string key = "R" + fipEnumToString_(i);
                if (summaryConfig.hasKeyword(key))
                    regionData[key] = regionFipValues[i];
            }
            if (summaryConfig.hasKeyword("RPR"))
                regionData["RPR"] = pressureAverage_(regPressurePvHydrocarbon, regPvHydrocarbon, regPressurePv, regionFipValues[FipDataType::PoreVolume], true);

            if (summaryConfig.hasKeyword("RPRP"))
                regionData["RPRP"] = pressureAverage_(regPressurePvHydrocarbon, regPvHydrocarbon, regPressurePv, regionFipValues[FipDataType::PoreVolume], false);

            // Output to log
            if (!substep) {

                fipUnitConvert_(fieldFipValues);
                if (origTotalValues_.empty())
                    origTotalValues_ = fieldFipValues;

                Scalar fieldHydroCarbonPoreVolumeAveragedPressure = pressureAverage_(fieldPressurePvHydrocarbon[0], fieldPvHydrocarbon[0], fieldPressurePv[0], fieldFipValues[FipDataType::PoreVolume], true);
                pressureUnitConvert_(fieldHydroCarbonPoreVolumeAveragedPressure);
                outputRegionFluidInPlace_(origTotalValues_, fieldFipValues, fieldHydroCarbonPoreVolumeAveragedPressure, 0);
                for (size_t reg = 0; reg < ntFip; ++reg ) {
                    ScalarBuffer tmpO(FipDataType::numFipValues,0.0);
                    for (int i = 0; i<FipDataType::numFipValues; i++) {
                        tmpO[i] = origRegionValues_[i][reg];
                    }
                    fipUnitConvert_(tmpO);
                    ScalarBuffer tmp(FipDataType::numFipValues,0.0);
                    for (int i = 0; i<FipDataType::numFipValues; i++) {
                        tmp[i] = regionFipValues[i][reg];
                    }
                    fipUnitConvert_(tmp);
                    Scalar regHydroCarbonPoreVolumeAveragedPressure = pressureAverage_(regPressurePvHydrocarbon[reg], regPvHydrocarbon[reg], regPressurePv[reg], regionFipValues[FipDataType::PoreVolume][reg], true);
                    pressureUnitConvert_(regHydroCarbonPoreVolumeAveragedPressure);
                    outputRegionFluidInPlace_(tmpO, tmp, regHydroCarbonPoreVolumeAveragedPressure, reg + 1);
                }
            }
        }

    }

    void setRestart(const Opm::data::Solution& sol, unsigned elemIdx, unsigned globalDofIndex) 
    {
        Scalar so = 1.0;
        if( saturation_[waterPhaseIdx].size() > 0 && sol.has( "SWAT" ) ) {
            saturation_[waterPhaseIdx][elemIdx] = sol.data("SWAT")[globalDofIndex];
            so -= sol.data("SWAT")[globalDofIndex];
        }

        if( saturation_[gasPhaseIdx].size() > 0 && sol.has( "SGAS" ) ) {
            saturation_[gasPhaseIdx][elemIdx] = sol.data("SGAS")[globalDofIndex];
            so -= sol.data("SGAS")[globalDofIndex];
        }

        assert(saturation_[oilPhaseIdx].size() > 0);
        saturation_[oilPhaseIdx][elemIdx] = so;

        if( oilPressure_.size() > 0 && sol.has( "PRESSURE" ) ) {
            oilPressure_[elemIdx] = sol.data( "PRESSURE" )[globalDofIndex];
        }

        if( temperature_.size() > 0 && sol.has( "TEMP" ) ) {
            temperature_[elemIdx] = sol.data( "TEMP" )[globalDofIndex];
        }

        if( rs_.size() > 0 && sol.has( "RS" ) ) {
            rs_[elemIdx] = sol.data("RS")[globalDofIndex];
        }

        if( rv_.size() > 0 && sol.has( "RV" ) ) {
            rv_[elemIdx] = sol.data("RV")[globalDofIndex];
        }

        if ( sSol_.size() > 0 && sol.has( "SSOL" ) ) {
            sSol_[elemIdx] = sol.data("SSOL")[globalDofIndex];
        }

        if ( cPolymer_.size() > 0 && sol.has("POLYMER" ) ) {
            cPolymer_[elemIdx] = sol.data("POLYMER")[globalDofIndex];
        }

        if ( soMax_.size() > 0 && sol.has("SOMAX" ) ) {
            soMax_[elemIdx] = sol.data("SOMAX")[globalDofIndex];
        }

        if ( pcSwMdcOw_.size() > 0 &&sol.has("PCSWM_OW" ) ) {
            pcSwMdcOw_[elemIdx] = sol.data("PCSWM_OW")[globalDofIndex];
        }

        if ( krnSwMdcOw_.size() > 0 && sol.has("KRNSW_OW" ) ) {
            krnSwMdcOw_[elemIdx] = sol.data("KRNSW_OW")[globalDofIndex];
        }

        if ( pcSwMdcGo_.size() > 0 && sol.has("PCSWM_GO" ) ) {
            pcSwMdcGo_[elemIdx] = sol.data("PCSWM_GO")[globalDofIndex];
        }

        if ( krnSwMdcGo_.size() > 0 && sol.has("KRNSW_GO" ) ) {
            krnSwMdcGo_[elemIdx] = sol.data("KRNSW_GO")[globalDofIndex];
        }
    }


    template <class FluidState>
    void assignToFluidState(FluidState& fs, unsigned elemIdx) const 
    {

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            if (saturation_[phaseIdx].size() == 0)
                continue;

            fs.setSaturation(phaseIdx, saturation_[phaseIdx][elemIdx]);
        }

        if (oilPressure_.size() > 0) {

            // this assumes that capillary pressures only depend on the phase saturations
            // and possibly on temperature. (this is always the case for ECL problems.)
            Dune::FieldVector< Scalar, numPhases > pc( 0 );
            const MaterialLawParams& matParams = simulator_.problem().materialLawParams(elemIdx);
            MaterialLaw::capillaryPressures(pc, matParams, fs);
            Opm::Valgrind::CheckDefined(oilPressure_[elemIdx]);
            Opm::Valgrind::CheckDefined(pc);
            assert(FluidSystem::phaseIsActive(oilPhaseIdx));
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;

                fs.setPressure(phaseIdx, oilPressure_[elemIdx] + (pc[phaseIdx] - pc[oilPhaseIdx]));
            }
        }

        if (temperature_.size() > 0) {
            fs.setTemperature( temperature_[elemIdx]);
        }

        if (rs_.size() > 0) {
           fs.setRs(rs_[elemIdx]);

        }
        if (rv_.size() > 0) {
           fs.setRv(rv_[elemIdx]);
        }
    }

    void initHysteresisParams(Simulator& simulator, unsigned elemIdx) const {

        if (soMax_.size() > 0)
            simulator.problem().setMaxOilSaturation(elemIdx, soMax_[elemIdx]);

        if (simulator.problem().materialLawManager()->enableHysteresis()) {
            auto matLawManager = simulator.problem().materialLawManager();

            if (pcSwMdcOw_.size() > 0 && krnSwMdcOw_.size() > 0) {
                matLawManager->setOilWaterHysteresisParams(
                            pcSwMdcOw_[elemIdx],
                            krnSwMdcOw_[elemIdx],
                            elemIdx);
            }
            if (pcSwMdcGo_.size() > 0 && krnSwMdcGo_.size() > 0) {
                matLawManager->setGasOilHysteresisParams(
                            pcSwMdcGo_[elemIdx],
                            krnSwMdcGo_[elemIdx],
                            elemIdx);
            }
        }

    }

    Scalar getSolventSaturation(unsigned elemIdx) const {
        if(sSol_.size() > 0)
            return sSol_[elemIdx];

        return 0;
    }

    Scalar getPolymerConcentration(unsigned elemIdx) const {
        if(cPolymer_.size() > 0)
            return cPolymer_[elemIdx];

        return 0;
    }

    const std::map<std::pair<std::string, int>, double>& getBlockData() {
        return blockData_;
    }

private:

    bool isIORank_() const
    {
        const auto& comm = simulator_.gridView().comm();
        return comm.rank() == 0;
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
        const double pv =
            elemCtx.simulator().model().dofTotalVolume(globalDofIdx)
            * intQuants.porosity().value();

        if (pressureTimesHydrocarbonVolume_.size() > 0 && pressureTimesPoreVolume_.size() > 0) {
            assert(hydrocarbonPoreVolume_.size() ==  pressureTimesHydrocarbonVolume_.size());
            assert(fip_[FipDataType::PoreVolume].size() == pressureTimesPoreVolume_.size() );

            fip_[FipDataType::PoreVolume][globalDofIdx] = pv;

            Scalar hydrocarbon = 0.0;
            if (FluidSystem::phaseIsActive(oilPhaseIdx))
                hydrocarbon += Opm::getValue(fs.saturation(oilPhaseIdx));
            if (FluidSystem::phaseIsActive(gasPhaseIdx))
                hydrocarbon += Opm::getValue(fs.saturation(gasPhaseIdx));

            hydrocarbonPoreVolume_[globalDofIdx] = pv * hydrocarbon;

            if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
                pressureTimesPoreVolume_[globalDofIdx] = Opm::getValue(fs.pressure(oilPhaseIdx)) * pv;
                pressureTimesHydrocarbonVolume_[globalDofIdx] = pressureTimesPoreVolume_[globalDofIdx] * hydrocarbon;
            }
        }

        if (computeFip_) {
            Scalar fip[FluidSystem::numPhases];
            for (unsigned phase = 0; phase < FluidSystem::numPhases; ++phase) {
                if (!FluidSystem::phaseIsActive(phase)) {
                    continue;
                }

                const double b = Opm::getValue(fs.invB(phase));
                const double s = Opm::getValue(fs.saturation(phase));
                fip[phase] = b * s * pv;
            }

            if (FluidSystem::phaseIsActive(oilPhaseIdx) && fip_[FipDataType::OilInPlace].size() > 0)
                fip_[FipDataType::OilInPlace][globalDofIdx] = fip[oilPhaseIdx];
            if (FluidSystem::phaseIsActive(gasPhaseIdx) && fip_[FipDataType::GasInPlace].size() > 0)
                fip_[FipDataType::GasInPlace][globalDofIdx] = fip[gasPhaseIdx];
            if (FluidSystem::phaseIsActive(waterPhaseIdx) && fip_[FipDataType::WaterInPlace].size() > 0)
                fip_[FipDataType::WaterInPlace][globalDofIdx] = fip[waterPhaseIdx];

            // Store the pure oil and gas Fip
            if (FluidSystem::phaseIsActive(oilPhaseIdx) && fip_[FipDataType::OilInPlaceInLiquidPhase].size() > 0)
                fip_[FipDataType::OilInPlaceInLiquidPhase][globalDofIdx] = fip[oilPhaseIdx];

            if (FluidSystem::phaseIsActive(gasPhaseIdx) && fip_[FipDataType::GasInPlaceInGasPhase].size() > 0)
                fip_[FipDataType::GasInPlaceInGasPhase][globalDofIdx] = fip[gasPhaseIdx];

            if (FluidSystem::phaseIsActive(oilPhaseIdx) && FluidSystem::phaseIsActive(gasPhaseIdx)) {
                // Gas dissolved in oil and vaporized oil
                Scalar gasInPlaceLiquid = Opm::getValue(fs.Rs()) * fip[oilPhaseIdx];
                Scalar oilInPlaceGas = Opm::getValue(fs.Rv()) * fip[gasPhaseIdx];
                if (fip_[FipDataType::GasInPlaceInGasPhase].size() > 0)
                    fip_[FipDataType::GasInPlaceInLiquidPhase][globalDofIdx] = gasInPlaceLiquid;
                if (fip_[FipDataType::OilInPlaceInGasPhase].size() > 0)
                    fip_[FipDataType::OilInPlaceInGasPhase][globalDofIdx] = oilInPlaceGas;

                // Add dissolved gas and vaporized oil to total Fip
                if (fip_[FipDataType::OilInPlace].size() > 0)
                    fip_[FipDataType::OilInPlace][globalDofIdx] += oilInPlaceGas;
                if (fip_[FipDataType::GasInPlace].size() > 0)
                    fip_[FipDataType::GasInPlace][globalDofIdx] += gasInPlaceLiquid;
            }
        }

    }

    void createLocalFipnum_()
    {
        const std::vector<int>& fipnum_global = simulator_.vanguard().eclState().get3DProperties().getIntGridProperty("FIPNUM").getData();
        // Get compressed cell fipnum.
        const auto& gridView = simulator_.vanguard().gridView();
        unsigned numElements = gridView.size(/*codim=*/0);
        fipnum_.resize(numElements, 0.0);
        if (!fipnum_global.empty()) {
            ElementContext elemCtx(simulator_);
            ElementIterator elemIt = gridView.template begin</*codim=*/0>();
            const ElementIterator& elemEndIt = gridView.template end</*codim=*/0>();
            for (; elemIt != elemEndIt; ++elemIt) {
                const Element& elem = *elemIt;
                if (elem.partitionType() != Dune::InteriorEntity)
                    continue; // assign no fipnum regions to ghost elements

                elemCtx.updatePrimaryStencil(elem);
                const unsigned elemIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
                fipnum_[elemIdx] = fipnum_global[simulator_.vanguard().cartesianIndex( elemIdx )];
            }
        }
    }

    // Sum Fip values over regions.
    ScalarBuffer computeFipForRegions_(const ScalarBuffer& fip, std::vector<int>& regionId, size_t maxNumberOfRegions, bool commSum = true)
    {
        ScalarBuffer totals(maxNumberOfRegions, 0.0);

        if (fip.empty())
            return totals;

        assert(regionId.size() == fip.size());
        for (size_t j = 0; j < regionId.size(); ++j) {
            const int regionIdx = regionId[j] - 1;
            // the cell is not attributed to any region. ignore it!
            if (regionIdx < 0)
                continue;

            assert(regionIdx < static_cast<int>(maxNumberOfRegions));
            totals[regionIdx] += fip[j];
        }
        if (commSum) {
            const auto& comm = simulator_.gridView().comm();
            for (size_t i = 0; i < maxNumberOfRegions; ++i)
                totals[i] = comm.sum(totals[i]);
        }

        return totals;
    }

    ScalarBuffer pressureAverage_(const ScalarBuffer& pressurePvHydrocarbon, const ScalarBuffer& pvHydrocarbon, const ScalarBuffer& pressurePv, const ScalarBuffer& pv, bool hydrocarbon) {
        size_t size = pressurePvHydrocarbon.size();
        assert(pvHydrocarbon.size() == size);
        assert(pressurePv.size() == size);
        assert(pv.size() == size);

        ScalarBuffer fraction(size,0.0);
        for (size_t i = 0; i < size; ++i) {
            fraction[i] = pressureAverage_(pressurePvHydrocarbon[i], pvHydrocarbon[i], pressurePv[i], pv[i], hydrocarbon);
        }
        return fraction;
    }

    Scalar pressureAverage_(const Scalar& pressurePvHydrocarbon, const Scalar& pvHydrocarbon, const Scalar& pressurePv, const Scalar& pv, bool hydrocarbon) {
        if (pvHydrocarbon > 1e-10 && hydrocarbon)
            return pressurePvHydrocarbon / pvHydrocarbon;

        return pressurePv / pv;
    }

    void fipUnitConvert_(ScalarBuffer& fip)
    {
        const Opm::UnitSystem& units = simulator_.vanguard().eclState().getUnits();
        if (units.getType() == Opm::UnitSystem::UnitType::UNIT_TYPE_FIELD) {
            fip[FipDataType::WaterInPlace] = Opm::unit::convert::to(fip[FipDataType::WaterInPlace], Opm::unit::stb);
            fip[FipDataType::OilInPlace] = Opm::unit::convert::to(fip[FipDataType::OilInPlace], Opm::unit::stb);
            fip[FipDataType::OilInPlaceInLiquidPhase] = Opm::unit::convert::to(fip[FipDataType::OilInPlaceInLiquidPhase], Opm::unit::stb);
            fip[FipDataType::OilInPlaceInGasPhase] = Opm::unit::convert::to(fip[FipDataType::OilInPlaceInGasPhase], Opm::unit::stb);
            fip[FipDataType::GasInPlace] = Opm::unit::convert::to(fip[FipDataType::GasInPlace], 1000*Opm::unit::cubic(Opm::unit::feet));
            fip[FipDataType::GasInPlaceInLiquidPhase] = Opm::unit::convert::to(fip[FipDataType::GasInPlaceInLiquidPhase], 1000*Opm::unit::cubic(Opm::unit::feet));
            fip[FipDataType::GasInPlaceInGasPhase] = Opm::unit::convert::to(fip[FipDataType::GasInPlaceInGasPhase], 1000*Opm::unit::cubic(Opm::unit::feet));
            fip[FipDataType::PoreVolume] = Opm::unit::convert::to(fip[FipDataType::PoreVolume], Opm::unit::stb);
        }
        else if (units.getType() == Opm::UnitSystem::UnitType::UNIT_TYPE_METRIC) {
            // nothing to do
        }
        else {
            throw std::runtime_error("Unsupported unit type for fluid in place output.");
        }
    }

    void pressureUnitConvert_(Scalar& pav) {
        const Opm::UnitSystem& units = simulator_.vanguard().eclState().getUnits();
        if (units.getType() == Opm::UnitSystem::UnitType::UNIT_TYPE_FIELD) {
            pav = Opm::unit::convert::to(pav, Opm::unit::psia);
        } else if (units.getType() == Opm::UnitSystem::UnitType::UNIT_TYPE_METRIC) {
            pav = Opm::unit::convert::to(pav, Opm::unit::barsa);
        }
        else {
            throw std::runtime_error("Unsupported unit type for fluid in place output.");
        }
    }

    void outputRegionFluidInPlace_(const ScalarBuffer& oip, const ScalarBuffer& cip, const Scalar& pav, const int reg)
    {
        const Opm::UnitSystem& units = simulator_.vanguard().eclState().getUnits();
        std::ostringstream ss;
        if (!reg) {
            ss << "                                                  ===================================================\n"
               << "                                                  :                   Field Totals                  :\n";
        } else {
            ss << "                                                  ===================================================\n"
               << "                                                  :        FIPNUM report region  "
               << std::setw(2) << reg << "                 :\n";
        }
        if (units.getType() == Opm::UnitSystem::UnitType::UNIT_TYPE_METRIC) {
            ss << "                                                  :      PAV  =" << std::setw(14) << pav << " BARSA                 :\n"
               << std::fixed << std::setprecision(0)
               << "                                                  :      PORV =" << std::setw(14) << cip[FipDataType::PoreVolume] << "   RM3                 :\n";
            if (!reg) {
                ss << "                                                  : Pressure is weighted by hydrocarbon pore volume :\n"
                   << "                                                  : Porv volumes are taken at reference conditions  :\n";
            }
            ss << "                         :--------------- Oil    SM3 ---------------:-- Wat    SM3 --:--------------- Gas    SM3 ---------------:\n";
        }
        if (units.getType() == Opm::UnitSystem::UnitType::UNIT_TYPE_FIELD) {
            ss << "                                                  :      PAV  =" << std::setw(14) << pav << "  PSIA                 :\n"
               << std::fixed << std::setprecision(0)
               << "                                                  :      PORV =" << std::setw(14) << cip[FipDataType::PoreVolume] << "   RB                  :\n";
            if (!reg) {
                ss << "                                                  : Pressure is weighted by hydrocarbon pore volume :\n"
                   << "                                                  : Pore volumes are taken at reference conditions  :\n";
            }
            ss << "                         :--------------- Oil    STB ---------------:-- Wat    STB --:--------------- Gas   MSCF ---------------:\n";
        }
        ss << "                         :      Liquid        Vapour        Total   :      Total     :      Free        Dissolved       Total   :" << "\n"
           << ":------------------------:------------------------------------------:----------------:------------------------------------------:" << "\n"
           << ":Currently   in place    :" << std::setw(14) << cip[FipDataType::OilInPlaceInLiquidPhase] << std::setw(14) << cip[FipDataType::OilInPlaceInGasPhase] << std::setw(14) << cip[FipDataType::OilInPlace] << ":"
           << std::setw(13) << cip[FipDataType::WaterInPlace] << "   :" << std::setw(14) << (cip[FipDataType::GasInPlaceInGasPhase]) << std::setw(14) << cip[FipDataType::GasInPlaceInLiquidPhase] << std::setw(14) << cip[FipDataType::GasInPlace] << ":\n"
           << ":------------------------:------------------------------------------:----------------:------------------------------------------:\n"
           << ":Originally  in place    :" << std::setw(14) << oip[FipDataType::OilInPlaceInLiquidPhase] << std::setw(14) << oip[FipDataType::OilInPlaceInGasPhase] << std::setw(14) << oip[FipDataType::OilInPlace] << ":"
           << std::setw(13) << oip[FipDataType::WaterInPlace] << "   :" << std::setw(14) << oip[FipDataType::GasInPlaceInGasPhase] << std::setw(14) << oip[FipDataType::GasInPlaceInLiquidPhase] << std::setw(14) << oip[FipDataType::GasInPlace] << ":\n"
           << ":========================:==========================================:================:==========================================:\n";
        Opm::OpmLog::note(ss.str());
    }

    std::string fipEnumToString_(int i) {
        typedef typename FipDataType::FipId FipId;
        switch( static_cast<FipId>(i) )
        {
        case FipDataType::WaterInPlace: return "WIP";
        case FipDataType::OilInPlace: return "OIP";
        case FipDataType::GasInPlace: return "GIP";
        case FipDataType::OilInPlaceInLiquidPhase: return "OIPL";
        case FipDataType::OilInPlaceInGasPhase: return "OIPG";
        case FipDataType::GasInPlaceInLiquidPhase: return "GIPL";
        case FipDataType::GasInPlaceInGasPhase: return "GIPG";
        case FipDataType::PoreVolume: return "PV";
        }
        return "ERROR";
    }


    const Simulator& simulator_;

    bool outputFipRestart_;
    bool computeFip_;

    ScalarBuffer saturation_[numPhases];
    ScalarBuffer oilPressure_;
    ScalarBuffer temperature_;
    ScalarBuffer gasDissolutionFactor_;
    ScalarBuffer oilVaporizationFactor_;
    ScalarBuffer gasFormationVolumeFactor_;
    ScalarBuffer saturatedOilFormationVolumeFactor_;
    ScalarBuffer oilSaturationPressure_;
    ScalarBuffer rs_;
    ScalarBuffer rv_;
    ScalarBuffer invB_[numPhases];
    ScalarBuffer density_[numPhases];
    ScalarBuffer viscosity_[numPhases];
    ScalarBuffer relativePermeability_[numPhases];
    ScalarBuffer sSol_;
    ScalarBuffer cPolymer_;
    ScalarBuffer soMax_;
    ScalarBuffer pcSwMdcOw_;
    ScalarBuffer krnSwMdcOw_;
    ScalarBuffer pcSwMdcGo_;
    ScalarBuffer krnSwMdcGo_;
    ScalarBuffer bubblePointPressure_;
    ScalarBuffer dewPointPressure_;
    std::vector<int> failedCellsPb_;
    std::vector<int> failedCellsPd_;
    std::vector<int> fipnum_;
    ScalarBuffer fip_[FipDataType::numFipValues];
    ScalarBuffer origTotalValues_;
    ScalarBuffer origRegionValues_[FipDataType::numFipValues];
    ScalarBuffer hydrocarbonPoreVolume_;
    ScalarBuffer pressureTimesPoreVolume_;
    ScalarBuffer pressureTimesHydrocarbonVolume_;
    std::map<std::pair<std::string, int>, double> blockData_;
    std::map<size_t, Scalar> oilCompletionPressures_;
    std::map<size_t, Scalar> waterCompletionSaturations_;
    std::map<size_t, Scalar> gasCompletionSaturations_;
};
} // namespace Ewoms

#endif
