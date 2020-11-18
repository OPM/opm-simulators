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
class EclOutputBlackOilModule
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

    enum { numPhases = FluidSystem::numPhases };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { gasCompIdx = FluidSystem::gasCompIdx };
    enum { oilCompIdx = FluidSystem::oilCompIdx };
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };

    typedef std::vector<Scalar> ScalarBuffer;
    typedef std::vector<std::string> StringBuffer;

    struct FipDataType
    {
        enum FipId
        {
            WaterInPlace = 0, //WIP
            OilInPlace = 1, //OIP
            GasInPlace = 2, //GIP
            OilInPlaceInLiquidPhase = 3, //OIPL
            OilInPlaceInGasPhase = 4, //OIPG
            GasInPlaceInLiquidPhase = 5, //GIPL
            GasInPlaceInGasPhase = 6,  //GIPG
            PoreVolume = 7, //PV
        };
        static const int numFipTypes = PoreVolume + 1 ;


        static std::string EclString(int fip_type) {
            switch(static_cast<FipId>(fip_type)) {
            case FipDataType::WaterInPlace: return "WIP";
            case FipDataType::OilInPlace: return "OIP";
            case FipDataType::GasInPlace: return "GIP";
            case FipDataType::OilInPlaceInLiquidPhase: return "OIPL";
            case FipDataType::OilInPlaceInGasPhase: return "OIPG";
            case FipDataType::GasInPlaceInLiquidPhase: return "GIPL";
            case FipDataType::GasInPlaceInGasPhase: return "GIPG";
            case FipDataType::PoreVolume: return "RPV";
            }
            throw std::logic_error("fip_type:  " + std::to_string(fip_type) + " not recognized");
        }
    };
    struct WellProdDataType
    {
        enum WPId
        {
            WellLocationi = 0, //WLi
            WellLocationj = 1, //WLj
            OilRate = 2, //OR
            WaterRate = 3, //WR
            GasRate = 4, //GR
            FluidResVol = 5, //FRV
            WaterCut = 6, //WC
            GasOilRatio = 7, //GOR
            WatGasRatio = 8, //WGR
            BHP = 9, //BHP
            THP = 10, //THP
            SteadyStatePI = 11, //SteadyStatePI
            WellName = 0, //WName
            CTRLMode = 1, //CTRL
        };
        static const int numWPValues = 12;
        static const int numWPNames = 2;
    };
        struct WellInjDataType
    {
        enum WIId
        {
            WellLocationi = 0, //WLi
            WellLocationj = 1, //WLj
            OilRate = 2, //OR
            WaterRate = 3, //WR
            GasRate = 4, //GR
            FluidResVol = 5, //FRV
            BHP = 6, //BHP
            THP = 7, //THP
            SteadyStateII = 8, //SteadyStateII
            WellName = 0, //WName
            CTRLModeOil = 1, //CTRLo
            CTRLModeWat = 2, //CTRLw
            CTRLModeGas = 3, //CTRLg
        };
        static const int numWIValues = 9;
        static const int numWINames = 4;
    };
        struct WellCumDataType
    {
        enum WCId
        {
            WellLocationi = 0, //WLi
            WellLocationj = 1, //WLj
            OilProd = 2, //OP
            WaterProd = 3, //WP
            GasProd = 4, //GP
            FluidResVolProd = 5, //FRVP
            OilInj = 6, //OI
            WaterInj = 7, //WI
            GasInj = 8, //GI
            FluidResVolInj = 9, //FRVI
            WellName = 0, //WName
            WellType = 1, //WType
            WellCTRL = 2, //WCTRL
        };
        static const int numWCValues = 10;
        static const int numWCNames = 3;
    };

    struct RegionSum {
        std::size_t ntFip;
        std::array<ScalarBuffer, FipDataType::numFipTypes> regFipValues;
        ScalarBuffer regPressurePv;
        ScalarBuffer regPvHydrocarbon;
        ScalarBuffer regPressurePvHydrocarbon;

        std::array<Scalar, FipDataType::numFipTypes> fieldFipValues;
        Scalar fieldPressurePv;
        Scalar fieldPvHydrocarbon;
        Scalar fieldPressurePvHydrocarbon;
    };

public:
    template<class CollectDataToIORankType>
    EclOutputBlackOilModule(const Simulator& simulator, const CollectDataToIORankType& collectToIORank)
        : simulator_(simulator)
    {
        const Opm::SummaryConfig summaryConfig = simulator_.vanguard().summaryConfig();
        const auto& fp = simulator_.vanguard().eclState().fieldProps();

        this->regions_["FIPNUM"] = fp.get_int("FIPNUM");
        for (const auto& region : summaryConfig.fip_regions())
            this->regions_[region] = fp.get_int(region);

        for (auto& region_pair : this->regions_)
            createLocalRegion_(region_pair.second);

        this->RPRNodes_  = summaryConfig.keywords("RPR*");
        this->RPRPNodes_ = summaryConfig.keywords("RPRP*");

        for (int fip_type = 0; fip_type < FipDataType::numFipTypes; fip_type++) {
            std::string key_pattern = "R" + FipDataType::EclString(fip_type) + "*";
            this->regionNodes_[fip_type] = summaryConfig.keywords(key_pattern);
        }

        for (const auto& node: summaryConfig) {
            if (node.category() == SummaryConfigNode::Category::Block) {
                if (collectToIORank.isCartIdxOnThisRank(node.number() - 1)) {
                    std::pair<std::string, int> key = std::make_pair(node.keyword(), node.number());
                    blockData_[key] = 0.0;
                }
            }
        }

        forceDisableFipOutput_ = EWOMS_GET_PARAM(TypeTag, bool, ForceDisableFluidInPlaceOutput);
    }

    /*!
     * \brief Register all run-time parameters for the Vtk output module.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, bool, ForceDisableFluidInPlaceOutput,
                             "Do not print fluid-in-place values after each report step even if requested by the deck.");
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to ECL output files
     */
    void allocBuffers(unsigned bufferSize, unsigned reportStepNum, const bool substep, const bool log, const bool isRestart)
    {
        if (!std::is_same<Discretization, Opm::EcfvDiscretization<TypeTag> >::value)
            return;

        // Summary output is for all steps
        const Opm::SummaryConfig summaryConfig = simulator_.vanguard().summaryConfig();

        // Only output RESTART_AUXILIARY asked for by the user.
        const Opm::RestartConfig& restartConfig = simulator_.vanguard().schedule().restart();
        std::map<std::string, int> rstKeywords = restartConfig.getRestartKeywords(reportStepNum);
        for (auto& [keyword, should_write] : rstKeywords) {
            if (this->isOutputCreationDirective_(keyword)) {
                // 'BASIC', 'FREQ' and similar.  Don't attempt to create
                // cell-based output for these keywords and don't warn about
                // not being able to create such cell-based result vectors.
                should_write = 0;
            }
            else {
                should_write = restartConfig.getKeyword(keyword, reportStepNum);
            }
        }

        outputFipRestart_ = false;
        computeFip_ = false;

        // Fluid in place
        for (int i = 0; i<FipDataType::numFipTypes; i++) {
            if (!substep || summaryConfig.require3DField(FipDataType::EclString(i))) {
                if (rstKeywords["FIP"] > 0) {
                    rstKeywords["FIP"] = 0;
                    outputFipRestart_ = true;
                }
                fip_[i].resize(bufferSize, 0.0);
                computeFip_ = true;
            }
            else
                fip_[i].clear();
        }
        if (!substep || summaryConfig.hasKeyword("FPR") || summaryConfig.hasKeyword("FPRP") || !this->RPRNodes_.empty()) {
            fip_[FipDataType::PoreVolume].resize(bufferSize, 0.0);
            hydrocarbonPoreVolume_.resize(bufferSize, 0.0);
            pressureTimesPoreVolume_.resize(bufferSize, 0.0);
            pressureTimesHydrocarbonVolume_.resize(bufferSize, 0.0);
        }
        else {
            hydrocarbonPoreVolume_.clear();
            pressureTimesPoreVolume_.clear();
            pressureTimesHydrocarbonVolume_.clear();
        }

        // Well RFT data
        if (!substep) {
            const auto& schedule = simulator_.vanguard().schedule();
            const auto& rft_config = schedule.rftConfig();
            for (const auto& well: schedule.getWells(reportStepNum)) {

                // don't bother with wells not on this process
                if (isDefunctParallelWell(well.name())) {
                    continue;
                }

                if (!rft_config.active(reportStepNum))
                    continue;

                for (const auto& connection: well.getConnections()) {
                    const size_t i = size_t(connection.getI());
                    const size_t j = size_t(connection.getJ());
                    const size_t k = size_t(connection.getK());
                    const size_t index = simulator_.vanguard().eclState().gridDims().getGlobalIndex(i, j, k);

                    if (FluidSystem::phaseIsActive(oilPhaseIdx))
                        oilConnectionPressures_.emplace(std::make_pair(index, 0.0));

                    if (FluidSystem::phaseIsActive(waterPhaseIdx))
                        waterConnectionSaturations_.emplace(std::make_pair(index, 0.0));

                    if (FluidSystem::phaseIsActive(gasPhaseIdx))
                        gasConnectionSaturations_.emplace(std::make_pair(index, 0.0));
                }
            }
        }

        // field data should be allocated
        // 1) when we want to restart
        // 2) when it is ask for by the user via restartConfig
        // 3) when it is not a substep
        if (!isRestart && (!restartConfig.getWriteRestartFile(reportStepNum, log) || substep))
            return;

        // always output saturation of active phases
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx))
                continue;

            saturation_[phaseIdx].resize(bufferSize, 0.0);
        }
        // and oil pressure
        oilPressure_.resize(bufferSize, 0.0);
        rstKeywords["PRES"] = 0;
        rstKeywords["PRESSURE"] = 0;

        // allocate memory for temperature
        if (enableEnergy || rstKeywords["TEMP"]) {
            temperature_.resize(bufferSize, 0.0);
            rstKeywords["TEMP"] = 0;
        }

        if (FluidSystem::phaseIsActive(oilPhaseIdx))
            rstKeywords["SOIL"] = 0;
        if (FluidSystem::phaseIsActive(gasPhaseIdx))
            rstKeywords["SGAS"] = 0;
        if (FluidSystem::phaseIsActive(waterPhaseIdx))
            rstKeywords["SWAT"] = 0;

        if (FluidSystem::enableDissolvedGas()) {
            rs_.resize(bufferSize, 0.0);
            rstKeywords["RS"] = 0;
        }
        if (FluidSystem::enableVaporizedOil()) {
            rv_.resize(bufferSize, 0.0);
            rstKeywords["RV"] = 0;
        }

        if (getPropValue<TypeTag, Properties::EnableSolvent>())
            sSol_.resize(bufferSize, 0.0);
        if (getPropValue<TypeTag, Properties::EnablePolymer>())
            cPolymer_.resize(bufferSize, 0.0);
        if (getPropValue<TypeTag, Properties::EnableFoam>())
            cFoam_.resize(bufferSize, 0.0);
        if (getPropValue<TypeTag, Properties::EnableBrine>())
            cSalt_.resize(bufferSize, 0.0);
        if (getPropValue<TypeTag, Properties::EnableExtbo>()) {
            extboX_.resize(bufferSize, 0.0);
            extboY_.resize(bufferSize, 0.0);
            extboZ_.resize(bufferSize, 0.0);
            mFracOil_.resize(bufferSize, 0.0);
            mFracGas_.resize(bufferSize, 0.0);
            mFracCo2_.resize(bufferSize, 0.0);
        }

        if (simulator_.problem().vapparsActive())
            soMax_.resize(bufferSize, 0.0);

        if (simulator_.problem().materialLawManager()->enableHysteresis()) {
            pcSwMdcOw_.resize(bufferSize, 0.0);
            krnSwMdcOw_.resize(bufferSize, 0.0);
            pcSwMdcGo_.resize(bufferSize, 0.0);
            krnSwMdcGo_.resize(bufferSize, 0.0);
        }

        if (simulator_.vanguard().eclState().fieldProps().has_double("SWATINIT")) {
            ppcw_.resize(bufferSize, 0.0);
            rstKeywords["PPCW"] = 0;
        }

        if (FluidSystem::enableDissolvedGas() && rstKeywords["RSSAT"] > 0) {
            rstKeywords["RSSAT"] = 0;
            gasDissolutionFactor_.resize(bufferSize, 0.0);
        }
        if (FluidSystem::enableVaporizedOil() && rstKeywords["RVSAT"] > 0) {
            rstKeywords["RVSAT"] = 0;
            oilVaporizationFactor_.resize(bufferSize, 0.0);
        }

        if (FluidSystem::phaseIsActive(waterPhaseIdx) && rstKeywords["BW"] > 0) {
            rstKeywords["BW"] = 0;
            invB_[waterPhaseIdx].resize(bufferSize, 0.0);
        }
        if (FluidSystem::phaseIsActive(oilPhaseIdx) && rstKeywords["BO"] > 0) {
            rstKeywords["BO"] = 0;
            invB_[oilPhaseIdx].resize(bufferSize, 0.0);
        }
        if (FluidSystem::phaseIsActive(gasPhaseIdx) && rstKeywords["BG"] > 0) {
            rstKeywords["BG"] = 0;
            invB_[gasPhaseIdx].resize(bufferSize, 0.0);
        }

        if (rstKeywords["DEN"] > 0) {
            rstKeywords["DEN"] = 0;
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;
                density_[phaseIdx].resize(bufferSize, 0.0);
            }
        }
        const bool hasVWAT = (rstKeywords["VISC"] > 0) || (rstKeywords["VWAT"] > 0);
        const bool hasVOIL = (rstKeywords["VISC"] > 0) || (rstKeywords["VOIL"] > 0);
        const bool hasVGAS = (rstKeywords["VISC"] > 0) || (rstKeywords["VGAS"] > 0);
        rstKeywords["VISC"] = 0;

        if (FluidSystem::phaseIsActive(waterPhaseIdx) && hasVWAT) {
            rstKeywords["VWAT"] = 0;
            viscosity_[waterPhaseIdx].resize(bufferSize, 0.0);
        }
        if (FluidSystem::phaseIsActive(oilPhaseIdx) && hasVOIL > 0) {
            rstKeywords["VOIL"] = 0;
            viscosity_[oilPhaseIdx].resize(bufferSize, 0.0);
        }
        if (FluidSystem::phaseIsActive(gasPhaseIdx) && hasVGAS > 0) {
            rstKeywords["VGAS"] = 0;
            viscosity_[gasPhaseIdx].resize(bufferSize, 0.0);
        }

        if (FluidSystem::phaseIsActive(waterPhaseIdx) && rstKeywords["KRW"] > 0) {
            rstKeywords["KRW"] = 0;
            relativePermeability_[waterPhaseIdx].resize(bufferSize, 0.0);
        }
        if (FluidSystem::phaseIsActive(oilPhaseIdx) && rstKeywords["KRO"] > 0) {
            rstKeywords["KRO"] = 0;
            relativePermeability_[oilPhaseIdx].resize(bufferSize, 0.0);
        }
        if (FluidSystem::phaseIsActive(gasPhaseIdx) && rstKeywords["KRG"] > 0) {
            rstKeywords["KRG"] = 0;
            relativePermeability_[gasPhaseIdx].resize(bufferSize, 0.0);
        }

        if (rstKeywords["PBPD"] > 0)  {
            rstKeywords["PBPD"] = 0;
            bubblePointPressure_.resize(bufferSize, 0.0);
            dewPointPressure_.resize(bufferSize, 0.0);
        }

        // tracers
        const int numTracers = simulator_.problem().tracerModel().numTracers();
        if (numTracers > 0){
            tracerConcentrations_.resize(numTracers);
            for (int tracerIdx = 0; tracerIdx < numTracers; ++tracerIdx)
            {
                tracerConcentrations_[tracerIdx].resize(bufferSize, 0.0);
            }
        }

        // ROCKC
        if (rstKeywords["ROCKC"] > 0) {
            rstKeywords["ROCKC"] = 0;
            rockCompPorvMultiplier_.resize(bufferSize, 0.0);
            rockCompTransMultiplier_.resize(bufferSize, 0.0);
            swMax_.resize(bufferSize, 0.0);
            minimumOilPressure_.resize(bufferSize, 0.0);
            overburdenPressure_.resize(bufferSize, 0.0);
        }

        //Warn for any unhandled keyword
        if (log) {
            for (auto& keyValue: rstKeywords) {
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
            saturatedOilFormationVolumeFactor_.resize(bufferSize, 0.0);
        if (false)
            oilSaturationPressure_.resize(bufferSize, 0.0);

    }

    /*!
     * \brief Modify the internal buffers according to the intensive quanties relevant
     *        for an element
     */
    void processElement(const ElementContext& elemCtx)
    {
        if (!std::is_same<Discretization, Opm::EcfvDiscretization<TypeTag> >::value)
            return;

        const auto& problem = elemCtx.simulator().problem();
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
                if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
                    oilPressure_[globalDofIdx] = Opm::getValue(fs.pressure(oilPhaseIdx));
                }else{
                    // put pressure in oil pressure for output
                    if (FluidSystem::phaseIsActive(waterPhaseIdx)) {
                        oilPressure_[globalDofIdx] = Opm::getValue(fs.pressure(waterPhaseIdx));
                    } else {
                        oilPressure_[globalDofIdx] = Opm::getValue(fs.pressure(gasPhaseIdx));
                    }
                }
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

                if (extboX_.size() > 0 && phaseIdx==oilPhaseIdx)
                    viscosity_[phaseIdx][globalDofIdx] = Opm::getValue(intQuants.oilViscosity());
                else if (extboX_.size() > 0 && phaseIdx==gasPhaseIdx)
                    viscosity_[phaseIdx][globalDofIdx] = Opm::getValue(intQuants.gasViscosity());
                else
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

            if (cFoam_.size() > 0) {
                cFoam_[globalDofIdx] = intQuants.foamConcentration().value();
            }

            if (cSalt_.size() > 0) {
                cSalt_[globalDofIdx] = fs.saltConcentration().value();
            }

            if (extboX_.size() > 0) {
                extboX_[globalDofIdx] = intQuants.xVolume().value();
            }

            if (extboY_.size() > 0) {
                extboY_[globalDofIdx] = intQuants.yVolume().value();
            }

            if (extboZ_.size() > 0) {
                extboZ_[globalDofIdx] = intQuants.zFraction().value();
            }

            if (mFracCo2_.size() > 0) {
                const Scalar stdVolOil = Opm::getValue(fs.saturation(oilPhaseIdx))*Opm::getValue(fs.invB(oilPhaseIdx))
                                       + Opm::getValue(fs.saturation(gasPhaseIdx))*Opm::getValue(fs.invB(gasPhaseIdx))*Opm::getValue(fs.Rv());
                const Scalar stdVolGas = Opm::getValue(fs.saturation(gasPhaseIdx))*Opm::getValue(fs.invB(gasPhaseIdx))*(1.0-intQuants.yVolume().value())
                                       + Opm::getValue(fs.saturation(oilPhaseIdx))*Opm::getValue(fs.invB(oilPhaseIdx))*Opm::getValue(fs.Rs())*(1.0-intQuants.xVolume().value());
                const Scalar stdVolCo2 = Opm::getValue(fs.saturation(gasPhaseIdx))*Opm::getValue(fs.invB(gasPhaseIdx))*intQuants.yVolume().value()
                                       + Opm::getValue(fs.saturation(oilPhaseIdx))*Opm::getValue(fs.invB(oilPhaseIdx))*Opm::getValue(fs.Rs())*intQuants.xVolume().value();
                const Scalar rhoO= FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx);
                const Scalar rhoG= FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx);
                const Scalar rhoCO2= intQuants.zRefDensity();
                const Scalar stdMassTotal= 1.0e-10 + stdVolOil*rhoO + stdVolGas*rhoG + stdVolCo2*rhoCO2;
                mFracOil_[globalDofIdx] = stdVolOil*rhoO/stdMassTotal;
                mFracGas_[globalDofIdx] = stdVolGas*rhoG/stdMassTotal;
                mFracCo2_[globalDofIdx] = stdVolCo2*rhoCO2/stdMassTotal;
            }

            if (bubblePointPressure_.size() > 0) {
                try {
                    bubblePointPressure_[globalDofIdx] = Opm::getValue(FluidSystem::bubblePointPressure(fs, intQuants.pvtRegionIndex()));
                }
                catch (const Opm::NumericalIssue&) {
                    const auto cartesianIdx = elemCtx.simulator().vanguard().grid().globalCell()[globalDofIdx];
                    failedCellsPb_.push_back(cartesianIdx);
                }
            }
            if (dewPointPressure_.size() > 0) {
                try {
                    dewPointPressure_[globalDofIdx] = Opm::getValue(FluidSystem::dewPointPressure(fs, intQuants.pvtRegionIndex()));
                }
                catch (const Opm::NumericalIssue&) {
                    const auto cartesianIdx = elemCtx.simulator().vanguard().grid().globalCell()[globalDofIdx];
                    failedCellsPd_.push_back(cartesianIdx);
                }
            }

            if (soMax_.size() > 0)
                soMax_[globalDofIdx] =
                    std::max(Opm::getValue(fs.saturation(oilPhaseIdx)),
                             problem.maxOilSaturation(globalDofIdx));

            if (swMax_.size() > 0)
                swMax_[globalDofIdx] =
                    std::max(Opm::getValue(fs.saturation(waterPhaseIdx)),
                             problem.maxWaterSaturation(globalDofIdx));

            if (minimumOilPressure_.size() > 0)
                minimumOilPressure_[globalDofIdx] =
                    std::min(Opm::getValue(fs.pressure(oilPhaseIdx)),
                             problem.minOilPressure(globalDofIdx));

            if (overburdenPressure_.size() > 0)
                overburdenPressure_[globalDofIdx] = problem.overburdenPressure(globalDofIdx);

            if (rockCompPorvMultiplier_.size() > 0)
                rockCompPorvMultiplier_[globalDofIdx] = problem.template rockCompPoroMultiplier<Scalar>(intQuants, globalDofIdx);

            if (rockCompTransMultiplier_.size() > 0)
                rockCompTransMultiplier_[globalDofIdx] = problem.template rockCompTransMultiplier<Scalar>(intQuants, globalDofIdx);

            const auto& matLawManager = problem.materialLawManager();
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


            if (ppcw_.size() > 0) {
                ppcw_[globalDofIdx] = matLawManager->oilWaterScaledEpsInfoDrainage(globalDofIdx).maxPcow;
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
                if (rv_.size() > 0)
                    rv_[globalDofIdx] = fsInitial.Rv();

                if (rs_.size() > 0)
                    rs_[globalDofIdx] = fsInitial.Rs();

                // re-compute the volume factors, viscosities and densities if asked for
                if (density_[oilPhaseIdx].size() > 0)
                    density_[oilPhaseIdx][globalDofIdx] = FluidSystem::density(fsInitial,
                                                                               oilPhaseIdx,
                                                                               intQuants.pvtRegionIndex());
                if (density_[gasPhaseIdx].size() > 0)
                    density_[gasPhaseIdx][globalDofIdx] = FluidSystem::density(fsInitial,
                                                                               gasPhaseIdx,
                                                                               intQuants.pvtRegionIndex());

                if (invB_[oilPhaseIdx].size() > 0)
                    invB_[oilPhaseIdx][globalDofIdx] = FluidSystem::inverseFormationVolumeFactor(fsInitial,
                                                                                                 oilPhaseIdx,
                                                                                                 intQuants.pvtRegionIndex());
                if (invB_[gasPhaseIdx].size() > 0)
                    invB_[gasPhaseIdx][globalDofIdx] = FluidSystem::inverseFormationVolumeFactor(fsInitial,
                                                                                                 gasPhaseIdx,
                                                                                                 intQuants.pvtRegionIndex());
                if (viscosity_[oilPhaseIdx].size() > 0)
                    viscosity_[oilPhaseIdx][globalDofIdx] = FluidSystem::viscosity(fsInitial,
                                                                                   oilPhaseIdx,
                                                                                   intQuants.pvtRegionIndex());
                if (viscosity_[gasPhaseIdx].size() > 0)
                    viscosity_[gasPhaseIdx][globalDofIdx] = FluidSystem::viscosity(fsInitial,
                                                                                   gasPhaseIdx,
                                                                                   intQuants.pvtRegionIndex());
            }

            // Add fluid in Place values
            updateFluidInPlace_(elemCtx, dofIdx);

            // Adding block data
            const auto cartesianIdx = elemCtx.simulator().vanguard().grid().globalCell()[globalDofIdx];
            for (auto& val: blockData_) {
                const auto& key = val.first;
                int cartesianIdxBlock = key.second - 1;
                if (cartesianIdx == cartesianIdxBlock) {
                    if (key.first == "BWSAT")
                        val.second = Opm::getValue(fs.saturation(waterPhaseIdx));
                    else if (key.first == "BGSAT")
                        val.second = Opm::getValue(fs.saturation(gasPhaseIdx));
                    else if (key.first == "BOSAT")
                        val.second = 1. - Opm::getValue(fs.saturation(gasPhaseIdx)) - Opm::getValue(fs.saturation(waterPhaseIdx));
                    else if (key.first == "BPR")
                        val.second = Opm::getValue(fs.pressure(oilPhaseIdx));
                    else if (key.first == "BWKR" || key.first == "BKRW")
                        val.second = Opm::getValue(intQuants.relativePermeability(waterPhaseIdx));
                    else if (key.first == "BGKR" || key.first == "BKRG")
                        val.second = Opm::getValue(intQuants.relativePermeability(gasPhaseIdx));
                    else if (key.first == "BOKR" || key.first == "BKRO")
                        val.second = Opm::getValue(intQuants.relativePermeability(oilPhaseIdx));
                    else if (key.first == "BWPC")
                        val.second = Opm::getValue(fs.pressure(oilPhaseIdx)) - Opm::getValue(fs.pressure(waterPhaseIdx));
                    else if (key.first == "BGPC")
                        val.second = Opm::getValue(fs.pressure(gasPhaseIdx)) - Opm::getValue(fs.pressure(oilPhaseIdx));
                    else if (key.first == "BVWAT" || key.first == "BWVIS")
                        val.second = Opm::getValue(fs.viscosity(waterPhaseIdx));
                    else if (key.first == "BVGAS" || key.first == "BGVIS")
                        val.second = Opm::getValue(fs.viscosity(gasPhaseIdx));
                    else if (key.first == "BVOIL" || key.first == "BOVIS")
                        val.second = Opm::getValue(fs.viscosity(oilPhaseIdx));
                    else {
                        std::string logstring = "Keyword '";
                        logstring.append(key.first);
                        logstring.append("' is unhandled for output to file.");
                        Opm::OpmLog::warning("Unhandled output keyword", logstring);
                    }
                }
            }

            // Adding Well RFT data
            if (oilConnectionPressures_.count(cartesianIdx) > 0) {
                oilConnectionPressures_[cartesianIdx] = Opm::getValue(fs.pressure(oilPhaseIdx));
            }
            if (waterConnectionSaturations_.count(cartesianIdx) > 0) {
                waterConnectionSaturations_[cartesianIdx] = Opm::getValue(fs.saturation(waterPhaseIdx));
            }
            if (gasConnectionSaturations_.count(cartesianIdx) > 0) {
                gasConnectionSaturations_[cartesianIdx] = Opm::getValue(fs.saturation(gasPhaseIdx));
            }

            // tracers
            const auto& tracerModel = simulator_.problem().tracerModel();
            if (tracerConcentrations_.size()>0) {
                for (int tracerIdx = 0; tracerIdx < tracerModel.numTracers(); tracerIdx++){
                    if (tracerConcentrations_[tracerIdx].size() == 0)
                        continue;

                    tracerConcentrations_[tracerIdx][globalDofIdx] = tracerModel.tracerConcentration(tracerIdx, globalDofIdx);
                }
            }
        }
    }


    void outputErrorLog()
    {
        const size_t maxNumCellsFaillog = 20;

        int pbSize = failedCellsPb_.size(), pdSize = failedCellsPd_.size();
        std::vector<int> displPb, displPd, recvLenPb, recvLenPd;
        const auto& comm = simulator_.gridView().comm();

        if (isIORank_()) {
            displPb.resize(comm.size()+1, 0);
            displPd.resize(comm.size()+1, 0);
            recvLenPb.resize(comm.size());
            recvLenPd.resize(comm.size());
        }

        comm.gather(&pbSize, recvLenPb.data(), 1, 0);
        comm.gather(&pdSize, recvLenPd.data(), 1, 0);
        std::partial_sum(recvLenPb.begin(), recvLenPb.end(), displPb.begin()+1);
        std::partial_sum(recvLenPd.begin(), recvLenPd.end(), displPd.begin()+1);
        std::vector<int> globalFailedCellsPb, globalFailedCellsPd;

        if (isIORank_()) {
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
            const size_t maxElems = std::min(maxNumCellsFaillog, globalFailedCellsPb.size());
            for (size_t i = 1; i < maxElems; ++i) {
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
            const size_t maxElems = std::min(maxNumCellsFaillog, globalFailedCellsPd.size());
            for (size_t i = 1; i < maxElems; ++i) {
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
        const auto& schedule = simulator_.vanguard().schedule();
        const auto& rft_config = schedule.rftConfig();
        for (const auto& well: schedule.getWells(reportStepNum)) {

            // don't bother with wells not on this process
            if (isDefunctParallelWell(well.name())) {
                continue;
            }

            //add data infrastructure for shut wells
            if (!wellDatas.count(well.name())) {
                Opm::data::Well wellData;

                if (!rft_config.active(reportStepNum))
                    continue;

                wellData.connections.resize(well.getConnections().size());
                size_t count = 0;
                for (const auto& connection: well.getConnections()) {
                    const size_t i = size_t(connection.getI());
                    const size_t j = size_t(connection.getJ());
                    const size_t k = size_t(connection.getK());

                    const size_t index = simulator_.vanguard().eclState().gridDims().getGlobalIndex(i, j, k);
                    auto& connectionData = wellData.connections[count];
                    connectionData.index = index;
                    count++;
                }
                wellDatas.emplace(std::make_pair(well.name(), wellData));
            }

            Opm::data::Well& wellData = wellDatas.at(well.name());
            for (auto& connectionData: wellData.connections) {
                const auto index = connectionData.index;
                if (oilConnectionPressures_.count(index) > 0)
                    connectionData.cell_pressure = oilConnectionPressures_.at(index);
                if (waterConnectionSaturations_.count(index) > 0)
                    connectionData.cell_saturation_water = waterConnectionSaturations_.at(index);
                if (gasConnectionSaturations_.count(index) > 0)
                    connectionData.cell_saturation_gas = gasConnectionSaturations_.at(index);
            }
        }
        oilConnectionPressures_.clear();
        waterConnectionSaturations_.clear();
        gasConnectionSaturations_.clear();
    }

    /*!
     * \brief Move all buffers to data::Solution.
     */
    void assignToSolution(Opm::data::Solution& sol)
    {
        if (!std::is_same<Discretization, Opm::EcfvDiscretization<TypeTag>>::value)
            return;

        if (oilPressure_.size() > 0) {
            sol.insert("PRESSURE", Opm::UnitSystem::measure::pressure, std::move(oilPressure_), Opm::data::TargetType::RESTART_SOLUTION);
        }

        if (temperature_.size() > 0) {
            sol.insert("TEMP", Opm::UnitSystem::measure::temperature, std::move(temperature_), Opm::data::TargetType::RESTART_SOLUTION);
        }

        if (FluidSystem::phaseIsActive(waterPhaseIdx) && saturation_[waterPhaseIdx].size() > 0) {
            sol.insert("SWAT", Opm::UnitSystem::measure::identity, std::move(saturation_[waterPhaseIdx]), Opm::data::TargetType::RESTART_SOLUTION);
        }
        if (FluidSystem::phaseIsActive(gasPhaseIdx) && saturation_[gasPhaseIdx].size() > 0) {
            sol.insert("SGAS", Opm::UnitSystem::measure::identity, std::move(saturation_[gasPhaseIdx]), Opm::data::TargetType::RESTART_SOLUTION);
        }
        if (ppcw_.size() > 0) {
            sol.insert ("PPCW", Opm::UnitSystem::measure::pressure, std::move(ppcw_), Opm::data::TargetType::RESTART_SOLUTION);
        }

        if (gasDissolutionFactor_.size() > 0) {
            sol.insert("RSSAT", Opm::UnitSystem::measure::gas_oil_ratio, std::move(gasDissolutionFactor_), Opm::data::TargetType::RESTART_AUXILIARY);

        }
        if (oilVaporizationFactor_.size() > 0) {
            sol.insert("RVSAT", Opm::UnitSystem::measure::oil_gas_ratio, std::move(oilVaporizationFactor_), Opm::data::TargetType::RESTART_AUXILIARY);
        }
        if (rs_.size() > 0) {
            sol.insert("RS", Opm::UnitSystem::measure::gas_oil_ratio, std::move(rs_), Opm::data::TargetType::RESTART_SOLUTION);

        }
        if (rv_.size() > 0) {
            sol.insert("RV", Opm::UnitSystem::measure::oil_gas_ratio, std::move(rv_), Opm::data::TargetType::RESTART_SOLUTION);
        }
        if (invB_[waterPhaseIdx].size() > 0) {
            sol.insert("1OVERBW", Opm::UnitSystem::measure::water_inverse_formation_volume_factor, std::move(invB_[waterPhaseIdx]), Opm::data::TargetType::RESTART_AUXILIARY);
        }
        if (invB_[oilPhaseIdx].size() > 0) {
            sol.insert("1OVERBO", Opm::UnitSystem::measure::oil_inverse_formation_volume_factor, std::move(invB_[oilPhaseIdx]), Opm::data::TargetType::RESTART_AUXILIARY);
        }
        if (invB_[gasPhaseIdx].size() > 0) {
            sol.insert("1OVERBG", Opm::UnitSystem::measure::gas_inverse_formation_volume_factor, std::move(invB_[gasPhaseIdx]), Opm::data::TargetType::RESTART_AUXILIARY);
        }

        if (density_[waterPhaseIdx].size() > 0) {
            sol.insert("WAT_DEN", Opm::UnitSystem::measure::density, std::move(density_[waterPhaseIdx]), Opm::data::TargetType::RESTART_AUXILIARY);
        }
        if (density_[oilPhaseIdx].size() > 0) {
            sol.insert("OIL_DEN", Opm::UnitSystem::measure::density, std::move(density_[oilPhaseIdx]), Opm::data::TargetType::RESTART_AUXILIARY);
        }
        if (density_[gasPhaseIdx].size() > 0) {
            sol.insert("GAS_DEN", Opm::UnitSystem::measure::density, std::move(density_[gasPhaseIdx]), Opm::data::TargetType::RESTART_AUXILIARY);
        }

        if (viscosity_[waterPhaseIdx].size() > 0) {
            sol.insert("WAT_VISC", Opm::UnitSystem::measure::viscosity, std::move(viscosity_[waterPhaseIdx]), Opm::data::TargetType::RESTART_AUXILIARY);
        }
        if (viscosity_[oilPhaseIdx].size() > 0) {
            sol.insert("OIL_VISC", Opm::UnitSystem::measure::viscosity, std::move(viscosity_[oilPhaseIdx]), Opm::data::TargetType::RESTART_AUXILIARY);
        }
        if (viscosity_[gasPhaseIdx].size() > 0) {
            sol.insert("GAS_VISC", Opm::UnitSystem::measure::viscosity, std::move(viscosity_[gasPhaseIdx]), Opm::data::TargetType::RESTART_AUXILIARY);
        }

        if (relativePermeability_[waterPhaseIdx].size() > 0) {
            sol.insert("WATKR", Opm::UnitSystem::measure::identity, std::move(relativePermeability_[waterPhaseIdx]), Opm::data::TargetType::RESTART_AUXILIARY);
        }
        if (relativePermeability_[oilPhaseIdx].size() > 0) {
            sol.insert("OILKR", Opm::UnitSystem::measure::identity, std::move(relativePermeability_[oilPhaseIdx]), Opm::data::TargetType::RESTART_AUXILIARY);
        }
        if (relativePermeability_[gasPhaseIdx].size() > 0) {
            sol.insert("GASKR", Opm::UnitSystem::measure::identity, std::move(relativePermeability_[gasPhaseIdx]), Opm::data::TargetType::RESTART_AUXILIARY);
        }

        if (pcSwMdcOw_.size() > 0)
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
            sol.insert ("SSOLVENT", Opm::UnitSystem::measure::identity, std::move(sSol_), Opm::data::TargetType::RESTART_SOLUTION);

        if (extboX_.size() > 0)
            sol.insert ("SS_X", Opm::UnitSystem::measure::identity, std::move(extboX_), Opm::data::TargetType::RESTART_SOLUTION);

        if (extboY_.size() > 0)
            sol.insert ("SS_Y", Opm::UnitSystem::measure::identity, std::move(extboY_), Opm::data::TargetType::RESTART_SOLUTION);

        if (extboZ_.size() > 0)
            sol.insert ("SS_Z", Opm::UnitSystem::measure::identity, std::move(extboZ_), Opm::data::TargetType::RESTART_SOLUTION);

        if (mFracOil_.size() > 0)
            sol.insert ("STD_OIL", Opm::UnitSystem::measure::identity, std::move(mFracOil_), Opm::data::TargetType::RESTART_SOLUTION);

        if (mFracGas_.size() > 0)
            sol.insert ("STD_GAS", Opm::UnitSystem::measure::identity, std::move(mFracGas_), Opm::data::TargetType::RESTART_SOLUTION);

        if (mFracCo2_.size() > 0)
            sol.insert ("STD_CO2", Opm::UnitSystem::measure::identity, std::move(mFracCo2_), Opm::data::TargetType::RESTART_SOLUTION);

        if (cPolymer_.size() > 0)
            sol.insert ("POLYMER", Opm::UnitSystem::measure::identity, std::move(cPolymer_), Opm::data::TargetType::RESTART_SOLUTION);

        if (cFoam_.size() > 0)
            sol.insert ("FOAM", Opm::UnitSystem::measure::identity, std::move(cFoam_), Opm::data::TargetType::RESTART_SOLUTION);

        if (cSalt_.size() > 0)
            sol.insert ("SALT", Opm::UnitSystem::measure::salinity, std::move(cSalt_), Opm::data::TargetType::RESTART_SOLUTION);

        if (dewPointPressure_.size() > 0)
            sol.insert ("PDEW", Opm::UnitSystem::measure::pressure, std::move(dewPointPressure_), Opm::data::TargetType::RESTART_AUXILIARY);

        if (bubblePointPressure_.size() > 0)
            sol.insert ("PBUB", Opm::UnitSystem::measure::pressure, std::move(bubblePointPressure_), Opm::data::TargetType::RESTART_AUXILIARY);


        if (swMax_.size() > 0)
            sol.insert ("SWMAX", Opm::UnitSystem::measure::identity, std::move(swMax_), Opm::data::TargetType::RESTART_SOLUTION);

        if (minimumOilPressure_.size() > 0)
            sol.insert ("PRESROCC", Opm::UnitSystem::measure::pressure, std::move(minimumOilPressure_), Opm::data::TargetType::RESTART_SOLUTION);

        if (overburdenPressure_.size() > 0)
            sol.insert ("PRES_OVB", Opm::UnitSystem::measure::pressure, std::move(overburdenPressure_), Opm::data::TargetType::RESTART_SOLUTION);

        if (rockCompPorvMultiplier_.size() > 0)
            sol.insert ("PORV_RC", Opm::UnitSystem::measure::identity, std::move(rockCompPorvMultiplier_), Opm::data::TargetType::RESTART_SOLUTION);

        if (rockCompTransMultiplier_.size() > 0)
            sol.insert ("TMULT_RC", Opm::UnitSystem::measure::identity, std::move(rockCompTransMultiplier_), Opm::data::TargetType::RESTART_SOLUTION);


        // Fluid in place
        for (int i = 0; i<FipDataType::numFipTypes; i++) {
            if (outputFipRestart_ && fip_[i].size() > 0) {
                sol.insert(FipDataType::EclString(i),
                           Opm::UnitSystem::measure::volume,
                           fip_[i],
                           Opm::data::TargetType::SUMMARY);
            }
        }

        // tracers
        const auto& tracerModel = simulator_.problem().tracerModel();
        if (tracerConcentrations_.size() > 0) {
            for (int tracerIdx = 0; tracerIdx<tracerModel.numTracers(); tracerIdx++){
                std::string tmp = tracerModel.tracerName(tracerIdx) + "F";
                sol.insert(tmp, Opm::UnitSystem::measure::identity, std::move(tracerConcentrations_[tracerIdx]), Opm::data::TargetType::RESTART_SOLUTION);
            }
        }
    }

    int regionMax(const std::vector<int>& region) {
        const auto max_value = region.empty() ? 0 : *std::max_element(region.begin(), region.end());
        return this->simulator_.gridView().comm().max(max_value);
    }


    RegionSum makeRegionSum(const std::vector<int>& region, bool is_fipnum) {
        RegionSum rsum;
        rsum.ntFip = this->regionMax(region);

        // sum values over each region
        rsum.regPressurePv = this->regionSum(this->pressureTimesPoreVolume_, region, rsum.ntFip);
        rsum.regPvHydrocarbon = this->regionSum(this->hydrocarbonPoreVolume_, region, rsum.ntFip);
        rsum.regPressurePvHydrocarbon = this->regionSum(pressureTimesHydrocarbonVolume_, region, rsum.ntFip);
        for (int fip_type = 0; fip_type < FipDataType::numFipTypes; fip_type++)
            rsum.regFipValues[fip_type] = this->regionSum(this->fip_[fip_type], region, rsum.ntFip);



        // sum all region values to compute the field total
        rsum.fieldPressurePv = sum(rsum.regPressurePv);
        rsum.fieldPvHydrocarbon = sum(rsum.regPvHydrocarbon);
        rsum.fieldPressurePvHydrocarbon = sum(rsum.regPressurePvHydrocarbon);
        for (int fip_type = 0; fip_type < FipDataType::numFipTypes; fip_type++) {
            const auto& regionFip = rsum.regFipValues[fip_type];
            rsum.fieldFipValues[fip_type] = sum(regionFip);
        }

        if (is_fipnum) {
            // The first time the outputFipLog function is run we store the inplace values in
            // the initialInplace_ member. This has at least two problems:
            //
            //   o We really want the *initial* value - now we get the value after
            //     the first timestep.
            //
            //   o For restarted runs this is obviously wrong.
            //
            // Finally it is of course not desirable to mutate state in an output
            // routine.
            if (!this->regionInitialInplace_) {
                this->regionInitialInplace_.emplace();
                for (int fip_type = 0; fip_type < FipDataType::numFipTypes; fip_type++)
                    this->regionInitialInplace_.value()[fip_type] = rsum.regFipValues[fip_type];
            }

            if (!this->fieldInitialInplace_)
                this->fieldInitialInplace_ = rsum.fieldFipValues;
        }
        return rsum;
    }


    std::unordered_map<std::string, RegionSum> accumulateRegionSums() {
        std::unordered_map<std::string, RegionSum> rsum_map;
        const Opm::SummaryConfig summaryConfig = simulator_.vanguard().summaryConfig();

        rsum_map.emplace("FIPNUM", makeRegionSum(this->regions_.at("FIPNUM"), true));
        for (const auto& fip_region : summaryConfig.fip_regions()) {
            if (fip_region == "FIPNUM")
                continue;

            const auto& region = this->regions_.at(fip_region);
            rsum_map.emplace(fip_region, makeRegionSum(region, false));
        }
        return rsum_map;
    }

    Scalar sum(const ScalarBuffer& v) {
        return std::accumulate(v.begin(), v.end(), Scalar{0});
    }

    void updateSummaryRegionValues(const std::unordered_map<std::string, RegionSum>& rsum_map,
                                   std::map<std::string, double>& miscSummaryData,
                                   std::map<std::string, std::vector<double>>& regionData) const {

        const Opm::SummaryConfig summaryConfig = simulator_.vanguard().summaryConfig();

        // The field summary vectors should only use the FIPNUM based region sum.
        {
            const auto& rsum = rsum_map.at("FIPNUM");
            for (int fip_type = 0; fip_type < FipDataType::numFipTypes; fip_type++) {
                std::string key = "F" + FipDataType::EclString(fip_type);
                if (summaryConfig.hasKeyword(key))
                    miscSummaryData[key] = rsum.fieldFipValues[fip_type];
            }
            if (summaryConfig.hasKeyword("FOE") && this->fieldInitialInplace_)
                miscSummaryData["FOE"] = rsum.fieldFipValues[FipDataType::OilInPlace]
                    / this->fieldInitialInplace_.value()[FipDataType::OilInPlace];

            if (summaryConfig.hasKeyword("FPR"))
                miscSummaryData["FPR"] = pressureAverage_(rsum.fieldPressurePvHydrocarbon,
                                                          rsum.fieldPvHydrocarbon,
                                                          rsum.fieldPressurePv,
                                                          rsum.fieldFipValues[FipDataType::PoreVolume],
                                                          true);

            if (summaryConfig.hasKeyword("FPRP"))
                miscSummaryData["FPRP"] = pressureAverage_(rsum.fieldPressurePvHydrocarbon,
                                                           rsum.fieldPvHydrocarbon,
                                                           rsum.fieldPressurePv,
                                                           rsum.fieldFipValues[FipDataType::PoreVolume],
                                                           false);

        }

        // The region summary vectors should loop through the FIPxxx regions to
        // support the RPR__xxx summary keywords.
        {
            for (int fip_type = 0; fip_type < FipDataType::numFipTypes; fip_type++) {
                for (const auto& node : this->regionNodes_[fip_type]) {
                    const auto& rsum = rsum_map.at(node.fip_region());
                    regionData[node.keyword()] = rsum.regFipValues[fip_type];
                }
            }

            // The exact same quantity is calculated for RPR and RPRP - is that correct?
            for (const auto& node : this->RPRNodes_) {
                const auto& rsum = rsum_map.at(node.fip_region());
                regionData[node.keyword()] = pressureAverage_(rsum.regPressurePvHydrocarbon,
                                                              rsum.regPvHydrocarbon,
                                                              rsum.regPressurePv,
                                                              rsum.regFipValues[FipDataType::PoreVolume],
                                                              true);
            }

            for (const auto& node : this->RPRPNodes_) {
                const auto& rsum = rsum_map.at(node.fip_region());
                regionData[node.keyword()] = pressureAverage_(rsum.regPressurePvHydrocarbon,
                                                              rsum.regPvHydrocarbon,
                                                              rsum.regPressurePv,
                                                              rsum.regFipValues[FipDataType::PoreVolume],
                                                              false);
            }
        }
    }


    void outputFipLogImpl(const RegionSum& rsum) const {


        {
            Scalar fieldHydroCarbonPoreVolumeAveragedPressure = pressureAverage_(rsum.fieldPressurePvHydrocarbon,
                                                                                 rsum.fieldPvHydrocarbon,
                                                                                 rsum.fieldPressurePv,
                                                                                 rsum.fieldFipValues[FipDataType::PoreVolume],
                                                                                 true);
            std::array<Scalar, FipDataType::numFipTypes> initial_values = *this->fieldInitialInplace_;
            std::array<Scalar, FipDataType::numFipTypes> current_values = rsum.fieldFipValues;
            fipUnitConvert_(initial_values);
            fipUnitConvert_(current_values);

            pressureUnitConvert_(fieldHydroCarbonPoreVolumeAveragedPressure);
            outputRegionFluidInPlace_(initial_values,
                                      current_values,
                                      fieldHydroCarbonPoreVolumeAveragedPressure);
        }

        for (size_t reg = 0; reg < rsum.ntFip; ++reg) {
            std::array<Scalar, FipDataType::numFipTypes> initial_values;
            std::array<Scalar, FipDataType::numFipTypes> current_values;

            for (std::size_t fip_type = 0; fip_type < FipDataType::numFipTypes; fip_type++) {
                initial_values[fip_type] = this->regionInitialInplace_.value()[fip_type][reg];
                current_values[fip_type] = rsum.regFipValues[fip_type][reg];
            }
            fipUnitConvert_(initial_values);
            fipUnitConvert_(current_values);

            Scalar regHydroCarbonPoreVolumeAveragedPressure
                = pressureAverage_(rsum.regPressurePvHydrocarbon[reg],
                                   rsum.regPvHydrocarbon[reg],
                                   rsum.regPressurePv[reg],
                                   rsum.regFipValues[FipDataType::PoreVolume][reg],
                                   true);
            pressureUnitConvert_(regHydroCarbonPoreVolumeAveragedPressure);
            outputRegionFluidInPlace_(initial_values, current_values, regHydroCarbonPoreVolumeAveragedPressure, reg + 1);
        }
    }




    // write Fluid In Place to output log
    void outputFipLog(std::map<std::string, double>& miscSummaryData,  std::map<std::string, std::vector<double>>& regionData, const bool substep)
    {
        auto rsum_map = this->accumulateRegionSums();
        if (!isIORank_())
            return;

        updateSummaryRegionValues(rsum_map,
                                  miscSummaryData,
                                  regionData);

        if (!substep)
            outputFipLogImpl(rsum_map.at("FIPNUM"));
    }


    // write production report to output
    void outputProdLog(size_t reportStepNum, const bool substep, bool forceDisableProdOutput)
    {
                if (!substep) {

                        ScalarBuffer  tmp_values(WellProdDataType::numWPValues, 0.0);
                        StringBuffer  tmp_names(WellProdDataType::numWPNames, "");
                        outputProductionReport_(tmp_values, tmp_names, forceDisableProdOutput);

                        const auto& st = simulator_.vanguard().summaryState();
                        const auto& schedule = simulator_.vanguard().schedule();

                        for (const auto& gname: schedule.groupNames()) {

                                auto gName = static_cast<std::string>(gname);
                                auto get = [&st, &gName](const std::string& vector)
                                {
                                        const auto key = vector + ':' + gName;

                                        return st.has(key) ? st.get(key) : 0.0;
                                };

                                tmp_names[0] = gname;

                                tmp_values[2] = get("GOPR"); //WellProdDataType::OilRate
                                tmp_values[3] = get("GWPR"); //WellProdDataType::WaterRate
                                tmp_values[4] = get("GGPR"); //WellProdDataType::GasRate
                                tmp_values[5] = get("GVPR"); //WellProdDataType::FluidResVol
                                tmp_values[6] = get("GWCT"); //WellProdDataType::WaterCut
                                tmp_values[7] = get("GGOR"); //WellProdDataType::GasOilRatio
                                tmp_values[8] = get("GWPR")/get("GGPR"); //WellProdDataType::WaterGasRatio

                                outputProductionReport_(tmp_values, tmp_names, forceDisableProdOutput);
                        }

                        for (const auto& wname: schedule.wellNames(reportStepNum)) {

                                // don't bother with wells not on this process
                                if (isDefunctParallelWell(wname)) {
                                        continue;
                                }

                                const auto& well = schedule.getWell(wname, reportStepNum);

                                // Ignore injector wells
                                if (well.isInjector()){
                                        continue;
                                }

                                tmp_names[0] = wname;//WellProdDataType::WellName


                                auto wName = static_cast<std::string>(wname);
                                auto get = [&st, &wName](const std::string& vector)
                                {
                                        const auto key = vector + ':' + wName;

                                        return st.has(key) ? st.get(key) : 0.0;
                                };

                                const auto& controls = well.productionControls(st);
                                using CMode = ::Opm::Well::ProducerCMode;

                                auto fctl = [](const auto wmctl) -> std::string
                                {
                                        switch (wmctl) {
                                        case CMode::ORAT: return "ORAT";
                                        case CMode::WRAT: return "WRAT";
                                        case CMode::GRAT: return "GRAT";
                                        case CMode::LRAT: return "LRAT";
                                        case CMode::RESV: return "RESV";
                                        case CMode::THP:  return "THP";
                                        case CMode::BHP:  return "BHP";
                                        case CMode::CRAT: return "CRate";
                                        case CMode::GRUP: return "GRUP";
                                        default:
                                        {
                                                return "none";
                                        }
                                        }
                                };

                                tmp_names[1] = fctl(controls.cmode); //WellProdDataType::CTRLMode

                                tmp_values[0] = well.getHeadI() + 1;//WellProdDataType::WellLocationi
                                tmp_values[1] = well.getHeadJ() + 1;//WellProdDataType::WellLocationj
                                tmp_values[2] = get("WOPR"); //WellProdDataType::OilRate
                                tmp_values[3] = get("WWPR"); //WellProdDataType::WaterRate
                                tmp_values[4] = get("WGPR"); //WellProdDataType::GasRate
                                tmp_values[5] = get("WVPR"); //WellProdDataType::FluidResVol
                                tmp_values[6] = get("WWCT"); //WellProdDataType::WaterCut
                                tmp_values[7] = get("WGOR"); //WellProdDataType::GasOilRatio
                                tmp_values[8] = get("WWPR")/get("WGPR"); //WellProdDataType::WaterGasRatio
                                tmp_values[9] = get("WBHP"); //WellProdDataType::BHP
                                tmp_values[10] = get("WTHP"); //WellProdDataType::THP
                                //tmp_values[11] = 0; //WellProdDataType::SteadyStatePI //

                                outputProductionReport_(tmp_values, tmp_names, forceDisableProdOutput);

                        }
                }
    }


    // write injection report to output
    void outputInjLog(size_t reportStepNum, const bool substep, bool forceDisableInjOutput)
    {
                if (!substep) {
                        ScalarBuffer  tmp_values(WellInjDataType::numWIValues, 0.0);
                        StringBuffer  tmp_names(WellInjDataType::numWINames, "");
                        outputInjectionReport_(tmp_values, tmp_names, forceDisableInjOutput);

                        const auto& st = simulator_.vanguard().summaryState();
                        const auto& schedule = simulator_.vanguard().schedule();
                        for (const auto& gname: schedule.groupNames()) {

                                auto gName = static_cast<std::string>(gname);
                                auto get = [&st, &gName](const std::string& vector)
                                {
                                        const auto key = vector + ':' + gName;

                                        return st.has(key) ? st.get(key) : 0.0;
                                };

                                tmp_names[0] = gname;

                                tmp_values[2] = get("GOIR");//WellInjDataType::OilRate
                                tmp_values[3] = get("GWIR"); //WellInjDataType::WaterRate
                                tmp_values[4] = get("GGIR"); //WellInjDataType::GasRate
                                tmp_values[5] = get("GVIR");//WellInjDataType::FluidResVol

                                outputInjectionReport_(tmp_values, tmp_names, forceDisableInjOutput);
                        }

                        for (const auto& wname: schedule.wellNames(reportStepNum)) {

                                // don't bother with wells not on this process
                                if (isDefunctParallelWell(wname)) {
                                        continue;
                                }

                                const auto& well = schedule.getWell(wname, reportStepNum);

                                // Ignore Producer wells
                                if (well.isProducer()){
                                        continue;
                                }

                                tmp_names[0] = wname;   //WellInjDataType::WellName

                                auto wName = static_cast<std::string>(wname);
                                auto get = [&st, &wName](const std::string& vector)
                                {
                                        const auto key = vector + ':' + wName;

                                        return st.has(key) ? st.get(key) : 0.0;
                                };

                                const auto& controls = well.injectionControls(st);
                                const auto ctlMode = controls.cmode;
                                const auto injType = controls.injector_type;
                                using CMode = ::Opm::Well::InjectorCMode;
                                using WType = ::Opm::InjectorType;

                                auto ftype = [](const auto wtype) -> std::string
                                {
                                        switch (wtype) {
                                        case WType::OIL:   return "Oil";
                                        case WType::WATER: return "Wat";
                                        case WType::GAS:   return "Gas";
                                        case WType::MULTI: return "Multi";
                                        default:
                                        {
                                                return "";
                                        }
                                        }
                                };

                                auto fctl = [](const auto wmctl) -> std::string
                                {
                                        switch (wmctl) {
                                        case CMode::RATE: return "RATE";
                                        case CMode::RESV: return "RESV";
                                        case CMode::THP:  return "THP";
                                        case CMode::BHP:  return "BHP";
                                        case CMode::GRUP: return "GRUP";
                                        default:
                                        {
                                                return "";
                                        }
                                        }
                                };

                                const auto flowtype = ftype(injType);
                                const auto flowctl = fctl(ctlMode);
                                if(flowtype == "Oil") //WellInjDataType::CTRLModeOil
                                {
                                        if (flowctl == "RATE"){ tmp_names[1] = "ORAT"; }
                                        else { tmp_names[1] =  flowctl; }
                                }
                                else if (flowtype == "Wat") //WellInjDataType::CTRLModeWat
                                {
                                        if (flowctl == "RATE"){ tmp_names[3] = "WRAT"; }
                                        else { tmp_names[2] =  flowctl; }
                                }
                                else if (flowtype == "Gas") //WellInjDataType::CTRLModeGas
                                {
                                        if (flowctl == "RATE"){ tmp_names[3] = "GRAT"; }
                                        else { tmp_names[3] =  flowctl; }
                                }

                                tmp_values[0] = well.getHeadI() + 1; //WellInjDataType::wellLocationi
                                tmp_values[1] = well.getHeadJ() + 1; //WellInjDataType::wellLocationj
                                tmp_values[2] = get("WOIR"); //WellInjDataType::OilRate
                                tmp_values[3] = get("WWIR"); //WellInjDataType::WaterRate
                                tmp_values[4] = get("WGIR"); //WellInjDataType::GasRate
                                tmp_values[5] = get("WVIR");//WellInjDataType::FluidResVol
                                tmp_values[6] = get("WBHP"); //WellInjDataType::BHP
                                tmp_values[7] = get("WTHP"); //WellInjDataType::THP
                                //tmp_values[8] = 0; //WellInjDataType::SteadyStateII

                                outputInjectionReport_(tmp_values, tmp_names, forceDisableInjOutput);

                        }
                }
    }

    // write cumulative production and injection reports to output
    void outputCumLog(size_t reportStepNum, const bool substep, bool forceDisableCumOutput)
    {
                if (!substep) {
                        ScalarBuffer  tmp_values(WellCumDataType::numWCValues, 0.0);
                        StringBuffer  tmp_names(WellCumDataType::numWCNames, "");
                        outputCumulativeReport_(tmp_values, tmp_names, forceDisableCumOutput);

                        const auto& st = simulator_.vanguard().summaryState();
                        const auto& schedule = simulator_.vanguard().schedule();

                        for (const auto& gname: schedule.groupNames()) {

                                auto gName = static_cast<std::string>(gname);
                                auto get = [&st, &gName](const std::string& vector)
                                {
                                        const auto key = vector + ':' + gName;

                                        return st.has(key) ? st.get(key) : 0.0;
                                };

                                tmp_names[0] = gname;

                                tmp_values[2] = get("GOPT"); //WellCumDataType::OilProd
                                tmp_values[3] = get("GWPT"); //WellCumDataType::WaterProd
                                tmp_values[4] = get("GGPT"); //WellCumDataType::GasProd
                                tmp_values[5] = get("GVPT");//WellCumDataType::FluidResVolProd
                                tmp_values[6] = get("GOIT"); //WellCumDataType::OilInj
                                tmp_values[7] = get("GWIT"); //WellCumDataType::WaterInj
                                tmp_values[8] = get("GGIT"); //WellCumDataType::GasInj
                                tmp_values[9] = get("GVIT");//WellCumDataType::FluidResVolInj

                                outputCumulativeReport_(tmp_values, tmp_names, forceDisableCumOutput);
                        }

                        for (const auto& wname : schedule.wellNames(reportStepNum))  {

                                // don't bother with wells not on this process
                                if (isDefunctParallelWell(wname)) {
                                        continue;
                                }

                                const auto& well = schedule.getWell(wname, reportStepNum);

                                tmp_names[0] = wname; //WellCumDataType::WellName

                                auto wName = static_cast<std::string>(wname);
                                auto get = [&st, &wName](const std::string& vector)
                                {
                                        const auto key = vector + ':' + wName;

                                        return st.has(key) ? st.get(key) : 0.0;
                                };

                                if (well.isInjector()) {

                                        const auto& controls = well.injectionControls(st);
                                        const auto ctlMode = controls.cmode;
                                        const auto injType = controls.injector_type;
                                        using CMode = ::Opm::Well::InjectorCMode;
                                        using WType = ::Opm::InjectorType;

                                        auto ftype = [](const auto wtype) -> std::string
                                        {
                                            switch (wtype) {
                                            case WType::OIL:   return "Oil";
                                            case WType::WATER: return "Wat";
                                            case WType::GAS:   return "Gas";
                                            case WType::MULTI: return "Multi";
                                            default:
                                            {
                                                    return "";
                                            }
                                            }
                                        };

                                        auto fctl = [](const auto wmctl) -> std::string
                                        {
                                            switch (wmctl) {
                                            case CMode::RATE: return "RATE";
                                            case CMode::RESV: return "RESV";
                                            case CMode::THP:  return "THP";
                                            case CMode::BHP:  return "BHP";
                                            case CMode::GRUP: return "GRUP";
                                            default:
                                            {
                                                    return "";
                                            }
                                            }
                                        };

                                       tmp_names[1] = "INJ"; //WellCumDataType::WellType
                                       const auto flowctl = fctl(ctlMode);
                                       if (flowctl == "RATE") //WellCumDataType::WellCTRL
                                       {
                                                const auto flowtype = ftype(injType);
                                                if(flowtype == "Oil"){ tmp_names[2] = "ORAT"; }
                                                else if(flowtype == "Wat"){ tmp_names[2] = "WRAT"; }
                                                else if(flowtype == "Gas"){ tmp_names[2] = "GRAT"; }
                                        }
                                        else
                                        {
                                                tmp_names[2] = flowctl;
                                        }

                                }
                                else if (well.isProducer()) {

                                        const auto& controls = well.productionControls(st);
                                        using CMode = ::Opm::Well::ProducerCMode;

                                        auto fctl = [](const auto wmctl) -> std::string
                                        {
                                                switch (wmctl) {
                                                        case CMode::ORAT: return "ORAT";
                                                        case CMode::WRAT: return "WRAT";
                                                        case CMode::GRAT: return "GRAT";
                                                        case CMode::LRAT: return "LRAT";
                                                        case CMode::RESV: return "RESV";
                                                        case CMode::THP:  return "THP";
                                                        case CMode::BHP:  return "BHP";
                                                        case CMode::CRAT: return "CRAT";
                                                        case CMode::GRUP: return "GRUP";
                                                        default:
                                                        {
                                                                return "none";
                                                        }
                                                }
                                        };
                                        tmp_names[1] = "PROD"; //WellProdDataType::CTRLMode
                                        tmp_names[2] = fctl(controls.cmode); //WellProdDataType::CTRLMode
                                }

                                        tmp_values[0] = well.getHeadI() + 1; //WellCumDataType::wellLocationi
                                        tmp_values[1] = well.getHeadJ() + 1; //WellCumDataType::wellLocationj
                                        tmp_values[2] = get("WOPT"); //WellCumDataType::OilProd
                                        tmp_values[3] = get("WWPT"); //WellCumDataType::WaterProd
                                        tmp_values[4] = get("WGPT"); //WellCumDataType::GasProd
                                        tmp_values[5] = get("WVPT");//WellCumDataType::FluidResVolProd
                                        tmp_values[6] = get("WOIT"); //WellCumDataType::OilInj
                                        tmp_values[7] = get("WWIT"); //WellCumDataType::WaterInj
                                        tmp_values[8] = get("WGIT"); //WellCumDataType::GasInj
                                        tmp_values[9] = get("WVIT");//WellCumDataType::FluidResVolInj

                                        outputCumulativeReport_(tmp_values, tmp_names, forceDisableCumOutput);

                        }
                }
    }

    void setRestart(const Opm::data::Solution& sol, unsigned elemIdx, unsigned globalDofIndex)
    {
        Scalar so = 1.0;
        if (saturation_[waterPhaseIdx].size() > 0 && sol.has("SWAT")) {
            saturation_[waterPhaseIdx][elemIdx] = sol.data("SWAT")[globalDofIndex];
            so -= sol.data("SWAT")[globalDofIndex];
        }
        if (saturation_[gasPhaseIdx].size() > 0 && sol.has("SGAS")) {
            saturation_[gasPhaseIdx][elemIdx] = sol.data("SGAS")[globalDofIndex];
            so -= sol.data("SGAS")[globalDofIndex];
        }

        if (sSol_.size() > 0) {
            // keep the SSOL option for backward compatibility
            // should be removed after 10.2018 release
            if (sol.has("SSOL"))
                sSol_[elemIdx] = sol.data("SSOL")[globalDofIndex];
            else if (sol.has("SSOLVENT"))
                sSol_[elemIdx] = sol.data("SSOLVENT")[globalDofIndex];

            so -= sSol_[elemIdx];
        }

        assert(saturation_[oilPhaseIdx].size() > 0);
        saturation_[oilPhaseIdx][elemIdx] = so;

        if (oilPressure_.size() > 0 && sol.has("PRESSURE"))
            oilPressure_[elemIdx] = sol.data("PRESSURE")[globalDofIndex];
        if (temperature_.size() > 0 && sol.has("TEMP"))
            temperature_[elemIdx] = sol.data("TEMP")[globalDofIndex];
        if (rs_.size() > 0 && sol.has("RS"))
            rs_[elemIdx] = sol.data("RS")[globalDofIndex];
        if (rv_.size() > 0 && sol.has("RV"))
            rv_[elemIdx] = sol.data("RV")[globalDofIndex];
        if (cPolymer_.size() > 0 && sol.has("POLYMER"))
            cPolymer_[elemIdx] = sol.data("POLYMER")[globalDofIndex];
        if (cFoam_.size() > 0 && sol.has("FOAM"))
            cFoam_[elemIdx] = sol.data("FOAM")[globalDofIndex];
        if (cSalt_.size() > 0 && sol.has("SALT"))
            cSalt_[elemIdx] = sol.data("SALT")[globalDofIndex];
        if (soMax_.size() > 0 && sol.has("SOMAX"))
            soMax_[elemIdx] = sol.data("SOMAX")[globalDofIndex];
        if (pcSwMdcOw_.size() > 0 &&sol.has("PCSWM_OW"))
            pcSwMdcOw_[elemIdx] = sol.data("PCSWM_OW")[globalDofIndex];
        if (krnSwMdcOw_.size() > 0 && sol.has("KRNSW_OW"))
            krnSwMdcOw_[elemIdx] = sol.data("KRNSW_OW")[globalDofIndex];
        if (pcSwMdcGo_.size() > 0 && sol.has("PCSWM_GO"))
            pcSwMdcGo_[elemIdx] = sol.data("PCSWM_GO")[globalDofIndex];
        if (krnSwMdcGo_.size() > 0 && sol.has("KRNSW_GO"))
            krnSwMdcGo_[elemIdx] = sol.data("KRNSW_GO")[globalDofIndex];
        if (ppcw_.size() > 0 && sol.has("PPCW"))
            ppcw_[elemIdx] = sol.data("PPCW")[globalDofIndex];

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
            Dune::FieldVector< Scalar, numPhases > pc(0);
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

        if (temperature_.size() > 0)
            fs.setTemperature(temperature_[elemIdx]);
        if (rs_.size() > 0)
           fs.setRs(rs_[elemIdx]);
        if (rv_.size() > 0)
           fs.setRv(rv_[elemIdx]);
    }

    void initHysteresisParams(Simulator& simulator, unsigned elemIdx) const
    {
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

        if (simulator_.vanguard().eclState().fieldProps().has_double("SWATINIT")) {
            auto oilWaterScaledEpsInfoDrainage = simulator.problem().materialLawManager()->oilWaterScaledEpsInfoDrainagePointerReferenceHack(elemIdx);
            oilWaterScaledEpsInfoDrainage->maxPcow =  ppcw_[elemIdx];
        }

    }

    Scalar getSolventSaturation(unsigned elemIdx) const
    {
        if (sSol_.size() > elemIdx)
            return sSol_[elemIdx];

        return 0;
    }

    Scalar getPolymerConcentration(unsigned elemIdx) const
    {
        if (cPolymer_.size() > elemIdx)
            return cPolymer_[elemIdx];

        return 0;
    }

    Scalar getFoamConcentration(unsigned elemIdx) const
    {
        if (cFoam_.size() > elemIdx)
            return cFoam_[elemIdx];

        return 0;
    }

    Scalar getSaltConcentration(unsigned elemIdx) const
    {
        if (cSalt_.size() > elemIdx)
            return cSalt_[elemIdx];

        return 0;
    }

    const std::map<std::pair<std::string, int>, double>& getBlockData()
    { return blockData_; }

private:

    bool isDefunctParallelWell(std::string wname) const
    {
        if (simulator_.gridView().comm().size()==1)
            return false;
        const auto& parallelWells = simulator_.vanguard().parallelWells();
        std::pair<std::string,bool> value{wname, true};
        auto candidate = std::lower_bound(parallelWells.begin(), parallelWells.end(),
                                          value);
        return candidate == parallelWells.end() || *candidate != value;
    }

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
            assert(fip_[FipDataType::PoreVolume].size() == pressureTimesPoreVolume_.size());

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
            for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
                fip[phaseIdx] = 0.0;

                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;

                const double b = Opm::getValue(fs.invB(phaseIdx));
                const double s = Opm::getValue(fs.saturation(phaseIdx));
                fip[phaseIdx] = b * s * pv;
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
                if (fip_[FipDataType::GasInPlaceInLiquidPhase].size() > 0)
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

    // Sum Fip values over regions.
    ScalarBuffer regionSum(const ScalarBuffer& property, const std::vector<int>& regionId, size_t maxNumberOfRegions) const
    {
        ScalarBuffer totals(maxNumberOfRegions, 0.0);

        if (property.empty())
            return totals;

        assert(regionId.size() == property.size());
        for (size_t j = 0; j < regionId.size(); ++j) {
            const int regionIdx = regionId[j] - 1;
            // the cell is not attributed to any region. ignore it!
            if (regionIdx < 0)
                continue;

            assert(regionIdx < static_cast<int>(maxNumberOfRegions));
            totals[regionIdx] += property[j];
        }

        const auto& comm = simulator_.gridView().comm();
        for (size_t i = 0; i < maxNumberOfRegions; ++i)
            totals[i] = comm.sum(totals[i]);

        return totals;
    }

    ScalarBuffer pressureAverage_(const ScalarBuffer& pressurePvHydrocarbon, const ScalarBuffer& pvHydrocarbon, const ScalarBuffer& pressurePv, const ScalarBuffer& pv, bool hydrocarbon) const
    {
        size_t size = pressurePvHydrocarbon.size();
        assert(pvHydrocarbon.size() == size);
        assert(pressurePv.size() == size);
        assert(pv.size() == size);

        ScalarBuffer fraction(size, 0.0);
        for (size_t i = 0; i < size; ++i) {
            fraction[i] = pressureAverage_(pressurePvHydrocarbon[i], pvHydrocarbon[i], pressurePv[i], pv[i], hydrocarbon);
        }
        return fraction;
    }

    Scalar pressureAverage_(const Scalar& pressurePvHydrocarbon, const Scalar& pvHydrocarbon, const Scalar& pressurePv, const Scalar& pv, bool hydrocarbon) const
    {
        if (pvHydrocarbon > 1e-10 && hydrocarbon)
            return pressurePvHydrocarbon / pvHydrocarbon;

        return pressurePv / pv;
    }

    void fipUnitConvert_(std::array<Scalar, FipDataType::numFipTypes>& fip) const
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
        else if (units.getType() == Opm::UnitSystem::UnitType::UNIT_TYPE_LAB) {
            Scalar scc = Opm::unit::cubic(Opm::prefix::centi * Opm::unit::meter); //standard cubic cm.
            fip[FipDataType::WaterInPlace] = Opm::unit::convert::to(fip[FipDataType::WaterInPlace], scc);
            fip[FipDataType::OilInPlace] = Opm::unit::convert::to(fip[FipDataType::OilInPlace], scc);
            fip[FipDataType::OilInPlaceInLiquidPhase] = Opm::unit::convert::to(fip[FipDataType::OilInPlaceInLiquidPhase], scc);
            fip[FipDataType::OilInPlaceInGasPhase] = Opm::unit::convert::to(fip[FipDataType::OilInPlaceInGasPhase], scc);
            fip[FipDataType::GasInPlace] = Opm::unit::convert::to(fip[FipDataType::GasInPlace], scc);
            fip[FipDataType::GasInPlaceInLiquidPhase] = Opm::unit::convert::to(fip[FipDataType::GasInPlaceInLiquidPhase], scc);
            fip[FipDataType::GasInPlaceInGasPhase] = Opm::unit::convert::to(fip[FipDataType::GasInPlaceInGasPhase], scc);
            fip[FipDataType::PoreVolume] = Opm::unit::convert::to(fip[FipDataType::PoreVolume], scc);
        }
        else if (units.getType() == Opm::UnitSystem::UnitType::UNIT_TYPE_METRIC) {
            // nothing to do
        }
        else {
            throw std::runtime_error("Unsupported unit type for fluid in place output.");
        }
    }

    void pressureUnitConvert_(Scalar& pav) const
    {
        const Opm::UnitSystem& units = simulator_.vanguard().eclState().getUnits();
        if (units.getType() == Opm::UnitSystem::UnitType::UNIT_TYPE_FIELD) {
            pav = Opm::unit::convert::to(pav, Opm::unit::psia);
        }
        else if (units.getType() == Opm::UnitSystem::UnitType::UNIT_TYPE_METRIC) {
            pav = Opm::unit::convert::to(pav, Opm::unit::barsa);
        }
        else if (units.getType() == Opm::UnitSystem::UnitType::UNIT_TYPE_LAB) {
            pav = Opm::unit::convert::to(pav, Opm::unit::atm);

        }
        else {
            throw std::runtime_error("Unsupported unit type for fluid in place output.");
        }
    }

    void outputRegionFluidInPlace_(const std::array<Scalar, FipDataType::numFipTypes>& oip,
                                   const std::array<Scalar, FipDataType::numFipTypes>& cip,
                                   const Scalar& pav, const int reg = 0) const
    {
        if (forceDisableFipOutput_)
            return;

        // don't output FIPNUM report if the region has no porv.
        if (cip[FipDataType::PoreVolume] == 0)
            return;

        const Opm::UnitSystem& units = simulator_.vanguard().eclState().getUnits();
        std::ostringstream ss;
        if (reg == 0) {
            ss << "                                                  ===================================================\n"
               << "                                                  :                   Field Totals                  :\n";
        }
        else {
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

    void outputProductionReport_(const ScalarBuffer& wellProd, const StringBuffer& wellProdNames, const bool forceDisableProdOutput)
    {
                if(forceDisableProdOutput)
                        return;

        const Opm::UnitSystem& units = simulator_.vanguard().eclState().getUnits();
        std::ostringstream ss;
        if (wellProdNames[WellProdDataType::WellName].empty()) {
            ss << "======================================================= PRODUCTION REPORT =======================================================\n"//=================== \n"
               << ":  WELL  :  LOCATION :CTRL:    OIL    :   WATER   :    GAS    :   FLUID   :   WATER   : GAS/OIL  :  WAT/GAS   : BHP OR : THP OR :\n"// STEADY-ST PI       :\n"
               << ":  NAME  :  (I,J,K)  :MODE:    RATE   :   RATE    :    RATE   :  RES.VOL. :    CUT    :  RATIO   :   RATIO    : CON.PR.: BLK.PR.:\n";// OR POTN OF PREF. PH:\n";
            if (units.getType() == Opm::UnitSystem::UnitType::UNIT_TYPE_METRIC) {
                            ss << ":        :           :    :  SCM/DAY  :  SCM/DAY  :  SCM/DAY  :  RCM/DAY  :  SCM/SCM  :  SCM/SCM :  SCM/SCM   :  BARSA :  BARSA :\n";//                    :\n";
                        }
                        if (units.getType() == Opm::UnitSystem::UnitType::UNIT_TYPE_FIELD) {
                            ss << ":        :           :    :  STB/DAY  :  STB/DAY  :  MSCF/DAY :  RB/DAY   :           : MSCF/STB :  STB/MSCF  :  PSIA  :  PSIA  :\n";//                    :\n";
                        }
                    if (units.getType() == Opm::UnitSystem::UnitType::UNIT_TYPE_LAB) {
                                ss << ":        :           :    :  SCC/HR   :  SCC/HR   :  SCC/HR   :    RCC    :  SCC/SCC  :  SCC/SCC :  SCC/SCC   :  ATMA  :  ATMA  :\n";//                    :\n";
                        }
                ss << "=================================================================================================================================\n";//=================== \n";
        }
                else {
            if (wellProd[WellProdDataType::WellLocationi] < 1) {
                ss << std::right << std::fixed << ":" << std::setw (8) << wellProdNames[WellProdDataType::WellName] << ":" << std::setprecision(0) << std::setw(11) << "" << ":" << std::setw(4) << wellProdNames[WellProdDataType::CTRLMode] << ":" << std::setprecision(1) << std::setw(11) << wellProd[WellProdDataType::OilRate] << ":" << std::setw(11) << wellProd[WellProdDataType::WaterRate] << ":" <<  std::setw(11)<< wellProd[WellProdDataType::GasRate] << ":" <<  std::setw(11) << wellProd[WellProdDataType::FluidResVol] << std::setprecision(3) << ":" <<  std::setw(11) << wellProd[WellProdDataType::WaterCut] << std::setprecision(2) << ":" <<  std::setw(10) << wellProd[WellProdDataType::GasOilRatio] << std::setprecision(4) << ":" <<  std::setw(12) << wellProd[WellProdDataType::WatGasRatio] << std::setprecision(1) << ":" <<  std::setw(8) << "" << ":" <<  std::setw(8) << "" << ": \n";//wellProd[WellProdDataType::SteadyStatePI] << std::setw(10) << "\n"
            }
            else {
                ss << std::right << std::fixed << ":" << std::setw (8) << wellProdNames[WellProdDataType::WellName] << ":" << std::setprecision(0) << std::setw(5) << wellProd[WellProdDataType::WellLocationi] << "," << std::setw(5) << wellProd[WellProdDataType::WellLocationj] << ":" << std::setw(4) << wellProdNames[WellProdDataType::CTRLMode] << ":" << std::setprecision(1) << std::setw(11) << wellProd[WellProdDataType::OilRate] << ":" << std::setw(11) << wellProd[WellProdDataType::WaterRate] << ":" <<  std::setw(11)<< wellProd[WellProdDataType::GasRate] << ":" <<  std::setw(11) << wellProd[WellProdDataType::FluidResVol] << std::setprecision(3) << ":" <<  std::setw(11) << wellProd[WellProdDataType::WaterCut] << std::setprecision(2) << ":" <<  std::setw(10) << wellProd[WellProdDataType::GasOilRatio] << std::setprecision(4) << ":" <<  std::setw(12) << wellProd[WellProdDataType::WatGasRatio] << std::setprecision(1) << ":" <<  std::setw(8) << wellProd[WellProdDataType::BHP] << ":" <<  std::setw(8) << wellProd[WellProdDataType::THP] << ": \n";//wellProd[WellProdDataType::SteadyStatePI] << std::setw(10) << "\n"
            }
            ss << ":"<< std::setfill ('-') << std::setw (9) << ":" << std::setfill ('-') << std::setw (12) << ":" << std::setfill ('-') << std::setw (5) << ":" << std::setfill ('-') << std::setw (12) << ":" << std::setfill ('-') << std::setw (12) << ":" << std::setfill ('-') << std::setw (12) << ":" << std::setfill ('-') << std::setw (12) << ":" << std::setfill ('-') << std::setw (12) << ":" << std::setfill ('-') << std::setw (11) << ":" << std::setfill ('-') << std::setw (13) << ":" << std::setfill ('-') << std::setw (9) << ":" << std::setfill ('-') << std::setw (9) << ":" << "\n";
        }
                Opm::OpmLog::note(ss.str());
    }

    void outputInjectionReport_(const ScalarBuffer& wellInj, const StringBuffer& wellInjNames, const bool forceDisableInjOutput)
    {
                if(forceDisableInjOutput)
                        return;

        const Opm::UnitSystem& units = simulator_.vanguard().eclState().getUnits();
        std::ostringstream ss;
        if (wellInjNames[WellInjDataType::WellName].empty()) {
            ss << "=================================================== INJECTION REPORT ========================================\n"//===================== \n"
               << ":  WELL  :  LOCATION : CTRL : CTRL : CTRL :    OIL    :   WATER   :    GAS    :   FLUID   : BHP OR : THP OR :\n"// STEADY-ST II       :\n"
               << ":  NAME  :  (I,J,K)  : MODE : MODE : MODE :    RATE   :   RATE    :    RATE   :  RES.VOL. : CON.PR.: BLK.PR.:\n";// OR POTENTIAL       :\n";
            if (units.getType() == Opm::UnitSystem::UnitType::UNIT_TYPE_METRIC) {
                                ss << ":        :           : OIL  : WAT  : GAS  :  SCM/DAY  :  SCM/DAY  :  SCM/DAY  :  RCM/DAY  :  BARSA :  BARSA :\n";//                    :\n";
                        }
                        if (units.getType() == Opm::UnitSystem::UnitType::UNIT_TYPE_FIELD) {
                                ss << ":        :           : OIL  : WAT  : GAS  :  STB/DAY  :  STB/DAY  :  MSCF/DAY :  RB/DAY   :  PSIA  :  PSIA  :\n";//                    :\n";
                        }
                        if (units.getType() == Opm::UnitSystem::UnitType::UNIT_TYPE_LAB) {
                                ss << ":        :           : OIL  : WAT  : GAS  :   SCC/HR  :  SCC/HR   :  SCC/HR   :  RCC/HR   :  ATMA  :  ATMA  :\n";//                    :\n";
                        }
                ss << "==============================================================================================================\n";//===================== \n";
        }
                else {
            if (wellInj[WellInjDataType::WellLocationi] < 1) {
                ss  << std::right << std::fixed << std::setprecision(0) << ":" << std::setw (8) << wellInjNames[WellInjDataType::WellName] << ":" << std::setw(11) << "" << ":" << std::setw(6) << wellInjNames[WellInjDataType::CTRLModeOil] << ":" << std::setw(6) << wellInjNames[WellInjDataType::CTRLModeWat] << ":" << std::setw(6) << wellInjNames[WellInjDataType::CTRLModeGas] << ":" << std::setprecision(1) << std::setw(11) << wellInj[WellInjDataType::OilRate] << ":" << std::setw(11) << wellInj[WellInjDataType::WaterRate] << ":" << std::setw(11)<< wellInj[WellInjDataType::GasRate] << ":" << std::setw(11) << wellInj[WellInjDataType::FluidResVol] << ":" << std::setw(8)<< "" << ":" << std::setw(8)<< "" << ": \n";//wellInj[WellInjDataType::SteadyStateII] << std::setw(10) << "\n"
            }
            else {
                ss  << std::right << std::fixed << std::setprecision(0) << ":" << std::setw (8) << wellInjNames[WellInjDataType::WellName] << ":" << std::setw(5) << wellInj[WellInjDataType::WellLocationi] << "," << std::setw(5) << wellInj[WellInjDataType::WellLocationj] << ":" << std::setw(6) << wellInjNames[WellInjDataType::CTRLModeOil] << ":" << std::setw(6) << wellInjNames[WellInjDataType::CTRLModeWat] << ":" << std::setw(6) << wellInjNames[WellInjDataType::CTRLModeGas] << ":" << std::setprecision(1) << std::setw(11) << wellInj[WellInjDataType::OilRate] << ":" << std::setw(11) << wellInj[WellInjDataType::WaterRate] << ":" << std::setw(11)<< wellInj[WellInjDataType::GasRate] << ":" << std::setw(11) << wellInj[WellInjDataType::FluidResVol] << ":" << std::setw(8)<< wellInj[WellInjDataType::BHP] << ":" << std::setw(8)<< wellInj[WellInjDataType::THP] << ": \n";//wellInj[WellInjDataType::SteadyStateII] << std::setw(10) << "\n"
            }
            ss << ":--------:-----------:------:------:------:------------:----------:-----------:-----------:--------:--------: \n";//--------------------:\n";
        }
                Opm::OpmLog::note(ss.str());
    }

        void outputCumulativeReport_(const ScalarBuffer& wellCum, const StringBuffer& wellCumNames, const bool forceDisableCumOutput)
    {
                if(forceDisableCumOutput)
                        return;

        const Opm::UnitSystem& units = simulator_.vanguard().eclState().getUnits();
        std::ostringstream ss;
        if (wellCumNames[WellCumDataType::WellName].empty()) {
            ss << "=================================================== CUMULATIVE PRODUCTION/INJECTION REPORT =========================================\n"
               << ":  WELL  :  LOCATION :  WELL  :CTRL:    OIL    :   WATER   :    GAS    :   Prod    :    OIL    :   WATER   :    GAS    :   INJ     :\n"
               << ":  NAME  :  (I,J,K)  :  TYPE  :MODE:    PROD   :   PROD    :    PROD   :  RES.VOL. :    INJ    :   INJ     :    INJ    :  RES.VOL. :\n";
            if (units.getType() == Opm::UnitSystem::UnitType::UNIT_TYPE_METRIC) {
                                ss << ":        :           :        :    :    MSCM   :   MSCM    :    MMSCM  :   MRCM    :    MSCM   :   MSCM    :    MMSCM  :   MRCM    :\n";
                        }
                        if (units.getType() == Opm::UnitSystem::UnitType::UNIT_TYPE_FIELD) {
                                ss << ":        :           :        :    :    MSTB   :   MSTB    :    MMSCF  :   MRB     :    MSTB   :   MSTB    :    MMSCF  :   MRB     :\n";
                        }
                        if (units.getType() == Opm::UnitSystem::UnitType::UNIT_TYPE_LAB) {
                                ss << ":        :           :        :    :     MSCC  :   MSCC    :    MMSCC  :   MRCC    :    MSCC   :   MSCC    :    MMSCC  :   MRCC    :\n";
                        }
                ss << "====================================================================================================================================\n";
        }
                else {
            if (wellCum[WellCumDataType::WellLocationi] < 1) {
                ss  << std::right << std::fixed << std::setprecision(0) << ":" << std::setw (8) << wellCumNames[WellCumDataType::WellName] << ":" << std::setw(11) <<  "" << ":" << std::setw(8) << wellCumNames[WellCumDataType::WellType] << ":" << std::setw(4) << wellCumNames[WellCumDataType::WellCTRL] << ":" << std::setprecision(1) << std::setw(11) << wellCum[WellCumDataType::OilProd]/1000 << ":" << std::setw(11) << wellCum[WellCumDataType::WaterProd]/1000 << ":" << std::setw(11)<< wellCum[WellCumDataType::GasProd]/1000 << ":" << std::setw(11) << wellCum[WellCumDataType::FluidResVolProd]/1000 << ":" << std::setw(11) << wellCum[WellCumDataType::OilInj]/1000 << ":"  << std::setw(11) << wellCum[WellCumDataType::WaterInj]/1000 << ":" << std::setw(11) << wellCum[WellCumDataType::GasInj]/1000 << ":" << std::setw(11) << wellCum[WellCumDataType::FluidResVolInj]/1000 << ": \n";
            }
            else {
                ss  << std::right << std::fixed << std::setprecision(0) << ":" << std::setw (8) << wellCumNames[WellCumDataType::WellName] << ":" << std::setw(5) << wellCum[WellCumDataType::WellLocationi] << "," << std::setw(5) << wellCum[WellCumDataType::WellLocationj] << ":" << std::setw(8) << wellCumNames[WellCumDataType::WellType] << ":" << std::setw(4) << wellCumNames[WellCumDataType::WellCTRL] << ":" << std::setprecision(1) << std::setw(11) << wellCum[WellCumDataType::OilProd]/1000 << ":" << std::setw(11) << wellCum[WellCumDataType::WaterProd]/1000 << ":" << std::setw(11)<< wellCum[WellCumDataType::GasProd]/1000 << ":" << std::setw(11) << wellCum[WellCumDataType::FluidResVolProd]/1000 << ":" << std::setw(11) << wellCum[WellCumDataType::OilInj]/1000 << ":"  << std::setw(11) << wellCum[WellCumDataType::WaterInj]/1000 << ":" << std::setw(11) << wellCum[WellCumDataType::GasInj]/1000 << ":" << std::setw(11) << wellCum[WellCumDataType::FluidResVolInj]/1000 << ": \n";
            }
            ss << ":--------:-----------:--------:----:------------:----------:-----------:-----------:------------:----------:-----------:-----------: \n";
        }
        Opm::OpmLog::note(ss.str());
    }

    std::string WPEnumToString_(int i)
    {
        typedef typename WellProdDataType::WPId WPId;
        switch(static_cast<WPId>(i))
        {
                case WellProdDataType::WellName: return "WName";
                case WellProdDataType::WellLocationi: return "WLi";
                case WellProdDataType::WellLocationj: return "WLj";
                case WellProdDataType::CTRLMode: return "CTRL";
                case WellProdDataType::OilRate: return "OR";
                case WellProdDataType::WaterRate: return "WR";
                case WellProdDataType::GasRate: return "GR";
                case WellProdDataType::FluidResVol: return "FRV";
                case WellProdDataType::WaterCut: return "WC";
                case WellProdDataType::GasOilRatio: return "GOR";
                case WellProdDataType::WatGasRatio: return "WGR";
                case WellProdDataType::BHP: return "BHP";
                case WellProdDataType::THP: return "THP";
                case WellProdDataType::SteadyStatePI: return "SteadyStatePI";
        }
        return "ERROR";
    }
        std::string WIEnumToString_(int i)
    {
        typedef typename WellInjDataType::WIId WIId;
        switch(static_cast<WIId>(i))
        {
                case WellInjDataType::WellName: return "WName";
                case WellInjDataType::WellLocationi: return "WLi";
                case WellInjDataType::WellLocationj: return "WLj";
                case WellInjDataType::CTRLModeOil: return "CTRLo";
                case WellInjDataType::CTRLModeWat: return "CTRLw";
                case WellInjDataType::CTRLModeGas: return "CTRLg";
                case WellInjDataType::OilRate: return "OR";
                case WellInjDataType::WaterRate: return "WR";
                case WellInjDataType::GasRate: return "GR";
                case WellInjDataType::FluidResVol: return "FRV";
                case WellInjDataType::BHP: return "BHP";
                case WellInjDataType::THP: return "THP";
                case WellInjDataType::SteadyStateII: return "SteadyStateII";
        }
        return "ERROR";
    }
        std::string WCEnumToString_(int i)
    {
        typedef typename WellCumDataType::WCId WCId;
        switch(static_cast<WCId>(i))
        {
                case WellCumDataType::WellName: return "WName";
                case WellCumDataType::WellLocationi: return "WLi";
                case WellCumDataType::WellLocationj: return "WLj";
                case WellCumDataType::WellType: return "WType";
                case WellCumDataType::WellCTRL: return "WCTRL";
                case WellCumDataType::OilProd: return "OP";
                case WellCumDataType::WaterProd: return "WP";
                case WellCumDataType::GasProd: return "GP";
                case WellCumDataType::FluidResVolProd: return "FRVP";
                case WellCumDataType::OilInj: return "OI";
                case WellCumDataType::WaterInj: return "WI";
                case WellCumDataType::GasInj: return "GI";
                case WellCumDataType::FluidResVolInj: return "FRVI";
        }
        return "ERROR";
    }

    bool isOutputCreationDirective_(const std::string& keyword) const
    {
        return (keyword == "BASIC") || (keyword == "FREQ")
            || (keyword == "RESTART")                        // From RPTSCHED
            || (keyword == "SAVE")  || (keyword == "SFREQ"); // Not really supported
    }

    const Simulator& simulator_;

    bool outputFipRestart_;
    bool computeFip_;
    bool forceDisableFipOutput_;

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
    ScalarBuffer extboX_;
    ScalarBuffer extboY_;
    ScalarBuffer extboZ_;
    ScalarBuffer mFracOil_;
    ScalarBuffer mFracGas_;
    ScalarBuffer mFracCo2_;
    ScalarBuffer cPolymer_;
    ScalarBuffer cFoam_;
    ScalarBuffer cSalt_;
    ScalarBuffer soMax_;
    ScalarBuffer pcSwMdcOw_;
    ScalarBuffer krnSwMdcOw_;
    ScalarBuffer pcSwMdcGo_;
    ScalarBuffer krnSwMdcGo_;
    ScalarBuffer ppcw_;
    ScalarBuffer bubblePointPressure_;
    ScalarBuffer dewPointPressure_;
    ScalarBuffer rockCompPorvMultiplier_;
    ScalarBuffer rockCompTransMultiplier_;
    ScalarBuffer swMax_;
    ScalarBuffer overburdenPressure_;
    ScalarBuffer minimumOilPressure_;

    std::vector<int> failedCellsPb_;
    std::vector<int> failedCellsPd_;
    std::unordered_map<std::string, std::vector<int>> regions_;
    std::array<ScalarBuffer, FipDataType::numFipTypes> fip_;
    std::optional<std::array<ScalarBuffer, FipDataType::numFipTypes>> regionInitialInplace_;
    std::optional<std::array<Scalar, FipDataType::numFipTypes>> fieldInitialInplace_;

    std::vector<Opm::SummaryConfigNode> RPRNodes_;
    std::vector<Opm::SummaryConfigNode> RPRPNodes_;
    std::array<std::vector<Opm::SummaryConfigNode>, FipDataType::numFipTypes> regionNodes_;



    ScalarBuffer hydrocarbonPoreVolume_;
    ScalarBuffer pressureTimesPoreVolume_;
    ScalarBuffer pressureTimesHydrocarbonVolume_;
    std::map<std::pair<std::string, int>, double> blockData_;
    std::map<size_t, Scalar> oilConnectionPressures_;
    std::map<size_t, Scalar> waterConnectionSaturations_;
    std::map<size_t, Scalar> gasConnectionSaturations_;
    std::vector<ScalarBuffer> tracerConcentrations_;
};
} // namespace Opm

#endif
