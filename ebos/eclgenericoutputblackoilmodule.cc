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

#include <config.h>

#include <ebos/eclgenericoutputblackoilmodule.hh>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidsystems/BlackOilDefaultIndexTraits.hpp>

#include <opm/grid/common/CommunicationUtils.hpp>

#include <opm/output/data/Solution.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/input/eclipse/Schedule/RFTConfig.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#include <opm/input/eclipse/Units/Units.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include <fmt/format.h>

namespace {

std::string EclString(const Opm::Inplace::Phase phase)
{
    switch (phase) {
    case Opm::Inplace::Phase::WATER:
        return "WIP";

    case Opm::Inplace::Phase::OIL:
        return "OIP";

    case Opm::Inplace::Phase::GAS:
        return "GIP";

    case Opm::Inplace::Phase::OilInLiquidPhase:
        return "OIPL";

    case Opm::Inplace::Phase::OilInGasPhase:
        return "OIPG";

    case Opm::Inplace::Phase::GasInLiquidPhase:
        return "GIPL";

    case Opm::Inplace::Phase::GasInGasPhase:
        return "GIPG";

    case Opm::Inplace::Phase::PoreVolume:
        return "RPV";

    case Opm::Inplace::Phase::WaterResVolume:
        return "WIPR";

    case Opm::Inplace::Phase::OilResVolume:
        return "OIPR";

    case Opm::Inplace::Phase::GasResVolume:
        return "GIPR";

    case Opm::Inplace::Phase::SALT:
        return "SIP";

    case Opm::Inplace::Phase::CO2InWaterPhase:
        return "WCD";

    case Opm::Inplace::Phase::CO2InGasPhaseInMob:
        return "GCDI";

    case Opm::Inplace::Phase::CO2InGasPhaseMob:
        return "GCDM";

    case Opm::Inplace::Phase::WaterInGasPhase:
        return "WIPG";

    case Opm::Inplace::Phase::WaterInWaterPhase:
        return "WIPL";

    default:
        throw std::logic_error {
            fmt::format("Phase enum with integer value: "
                        "{} not recognized", static_cast<int>(phase))
        };
    }
}

    std::size_t numCells(const Opm::EclipseState& eclState)
    {
        return eclState.fieldProps().get_int("FIPNUM").size();
    }

    std::vector<Opm::EclInterRegFlowMap::SingleRegion>
    defineInterRegionFlowArrays(const Opm::EclipseState&  eclState,
                                const Opm::SummaryConfig& summaryConfig)
    {
        auto regions = std::vector<Opm::EclInterRegFlowMap::SingleRegion>{};

        const auto& fprops = eclState.fieldProps();
        for (const auto& arrayName : summaryConfig.fip_regions_interreg_flow()) {
            regions.push_back({ arrayName, std::cref(fprops.get_int(arrayName)) });
        }

        return regions;
    }
}

namespace Opm {

template<class FluidSystem, class Scalar>
EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
EclGenericOutputBlackoilModule(const EclipseState& eclState,
                           const Schedule& schedule,
                           const SummaryConfig& summaryConfig,
                           const SummaryState& summaryState,
                           bool enableEnergy,
                           bool enableTemperature,
                           bool enableSolvent,
                           bool enablePolymer,
                           bool enableFoam,
                           bool enableBrine,
                           bool enableSaltPrecipitation,
                           bool enableExtbo,
                           bool enableMICP)
    : eclState_(eclState)
    , schedule_(schedule)
    , summaryConfig_(summaryConfig)
    , summaryState_(summaryState)
    , interRegionFlows_(numCells(eclState), defineInterRegionFlowArrays(eclState, summaryConfig))
    , logOutput_(eclState, schedule, summaryState)
    , enableEnergy_(enableEnergy)
    , enableTemperature_(enableTemperature)
    , enableSolvent_(enableSolvent)
    , enablePolymer_(enablePolymer)
    , enableFoam_(enableFoam)
    , enableBrine_(enableBrine)
    , enableSaltPrecipitation_(enableSaltPrecipitation)
    , enableExtbo_(enableExtbo)
    , enableMICP_(enableMICP)
    , local_data_valid_(false)
{
    const auto& fp = eclState_.fieldProps();

    this->regions_["FIPNUM"] = fp.get_int("FIPNUM");
    for (const auto& region : summaryConfig_.fip_regions())
        this->regions_[region] = fp.get_int(region);

    this->RPRNodes_  = summaryConfig_.keywords("RPR*");
    this->RPRPNodes_ = summaryConfig_.keywords("RPRP*");

    for (const auto& phase : Inplace::phases()) {
        std::string key_pattern = "R" + EclString(phase) + "*";
        this->regionNodes_[phase] = summaryConfig_.keywords(key_pattern);
    }

    // Check if FLORES/FLOWS is set in any RPTRST in the schedule
    anyFlores_ = false;     // Used for the initialization of the sparse table
    anyFlows_ = false;
    enableFlores_ = false;  // Used for the output of i+, j+, k+
    enableFloresn_ = false; // Used for the special case of nnc
    enableFlows_ = false;
    enableFlowsn_ = false;

    for (const auto& block : this->schedule_) { // Uses Schedule::begin() and Schedule::end()
        const auto& rstkw = block.rst_config().keywords;

        if (! anyFlores_) {
            anyFlores_ = rstkw.find("FLORES") != rstkw.end();
        }

        if (! anyFlows_) {
            anyFlows_ = rstkw.find("FLOWS") != rstkw.end();
        }

        if (anyFlores_ && anyFlows_) {
            // Terminate report step loop early if both FLORES and FLOWS
            // have been set at some point as there's no need to search
            // any further in that case.
            break;
        }
    }
}

template<class FluidSystem, class Scalar>
EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
~EclGenericOutputBlackoilModule() = default;

template<class FluidSystem, class Scalar>
void EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
outputCumLog(size_t reportStepNum, const bool substep, bool forceDisableCumOutput)
{
    if (!substep && !forceDisableCumOutput) {
        logOutput_.cumulative(reportStepNum,
                              [this](const std::string& name)
                              { return this->isDefunctParallelWell(name); });
    }
}

template<class FluidSystem,class Scalar>
void EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
outputProdLog(size_t reportStepNum,
              const bool substep,
              bool forceDisableProdOutput)
{
    if (!substep) {
        ScalarBuffer  tmp_values(WellProdDataType::numWPValues, 0.0);
        StringBuffer  tmp_names(WellProdDataType::numWPNames, "");
        outputProductionReport_(tmp_values, tmp_names, forceDisableProdOutput);

        const auto& st = summaryState_;

        for (const auto& gname: schedule_.groupNames()) {

            auto gName = static_cast<std::string>(gname);
            auto get = [&st, &gName](const std::string& vector)
            {
                const auto key = vector + ':' + gName;

                return st.has(key) ? st.get(key) : 0.0;
            };

            tmp_names[0] = gname;

            if (tmp_names[0] == "FIELD"){
                tmp_values[2] = st.get("FOPR", 0.0); //WellProdDataType::OilRate
                tmp_values[3] = st.get("FWPR", 0.0); //WellProdDataType::WaterRate
                tmp_values[4] = st.get("FGPR", 0.0); //WellProdDataType::GasRate
                tmp_values[5] = st.get("FVPR", 0.0); //WellProdDataType::FluidResVol
                tmp_values[6] = st.get("FWCT", 0.0); //WellProdDataType::WaterCut
                tmp_values[7] = st.get("FGOR", 0.0); //WellProdDataType::GasOilRatio
            } else {
                tmp_values[2] = get("GOPR"); //WellProdDataType::OilRate
                tmp_values[3] = get("GWPR"); //WellProdDataType::WaterRate
                tmp_values[4] = get("GGPR"); //WellProdDataType::GasRate
                tmp_values[5] = get("GVPR"); //WellProdDataType::FluidResVol
                tmp_values[6] = get("GWCT"); //WellProdDataType::WaterCut
                tmp_values[7] = get("GGOR"); //WellProdDataType::GasOilRatio
            }

            tmp_values[8] = tmp_values[3]/tmp_values[4]; //WellProdDataType::WaterGasRatio
            if (isnan(tmp_values[8])){
                tmp_values[8] = 0.0;
            }

            outputProductionReport_(tmp_values, tmp_names, forceDisableProdOutput);
        }

        for (const auto& wname: schedule_.wellNames(reportStepNum)) {

            // don't bother with wells not on this process
            if (isDefunctParallelWell(wname)) {
                continue;
            }

            const auto& well = schedule_.getWell(wname, reportStepNum);

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
            using CMode = Well::ProducerCMode;

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
            tmp_values[9] = get("WBHP"); //WellProdDataType::BHP
            tmp_values[10] = get("WTHP"); //WellProdDataType::THP
            //tmp_values[11] = 0; //WellProdDataType::SteadyStatePI //

            tmp_values[8] = tmp_values[3]/tmp_values[4]; //WellProdDataType::WaterGasRatio
            if (isnan(tmp_values[8])){
                tmp_values[8] = 0.0;
            }

            outputProductionReport_(tmp_values, tmp_names, forceDisableProdOutput);
        }
    }
}

template<class FluidSystem,class Scalar>
void EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
outputInjLog(size_t reportStepNum, const bool substep, bool forceDisableInjOutput)
{
    if (!substep) {
        ScalarBuffer  tmp_values(WellInjDataType::numWIValues, 0.0);
        StringBuffer  tmp_names(WellInjDataType::numWINames, "");
        outputInjectionReport_(tmp_values, tmp_names, forceDisableInjOutput);

        const auto& st = summaryState_;
        for (const auto& gname: schedule_.groupNames()) {

            auto gName = static_cast<std::string>(gname);
            auto get = [&st, &gName](const std::string& vector)
            {
                const auto key = vector + ':' + gName;

                return st.has(key) ? st.get(key) : 0.0;
            };

            tmp_names[0] = gname;

            if (tmp_names[0] == "FIELD"){
                tmp_values[2] = st.get("FOIR", 0.0);//WellInjDataType::OilRate
                tmp_values[3] = st.get("FWIR", 0.0); //WellInjDataType::WaterRate
                tmp_values[4] = st.get("FGIR", 0.0); //WellInjDataType::GasRate
                tmp_values[5] = st.get("FVIR", 0.0);//WellInjDataType::FluidResVol
            }
            else {
                tmp_values[2] = get("GOIR");//WellInjDataType::OilRate
                tmp_values[3] = get("GWIR"); //WellInjDataType::WaterRate
                tmp_values[4] = get("GGIR"); //WellInjDataType::GasRate
                tmp_values[5] = get("GVIR");//WellInjDataType::FluidResVol
            }

            outputInjectionReport_(tmp_values, tmp_names, forceDisableInjOutput);
        }

        for (const auto& wname: schedule_.wellNames(reportStepNum)) {

            // don't bother with wells not on this process
            if (isDefunctParallelWell(wname)) {
                continue;
            }

            const auto& well = schedule_.getWell(wname, reportStepNum);

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
            using CMode = Well::InjectorCMode;
            using WType = InjectorType;

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

template<class FluidSystem,class Scalar>
Inplace EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
outputFipLog(std::map<std::string, double>& miscSummaryData,
             std::map<std::string, std::vector<double>>& regionData,
             const bool substep,
             const Parallel::Communication& comm)
{
    auto inplace = this->accumulateRegionSums(comm);
    if (comm.rank() != 0)
        return inplace;

    updateSummaryRegionValues(inplace,
                              miscSummaryData,
                              regionData);

    if (!substep)
        outputFipLogImpl(inplace);

    return inplace;
}

template<class FluidSystem,class Scalar>
Inplace EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
outputFipresvLog(std::map<std::string, double>& miscSummaryData,
             std::map<std::string, std::vector<double>>& regionData,
             const bool substep,
             const Parallel::Communication& comm)
{
    auto inplace = this->accumulateRegionSums(comm);
    if (comm.rank() != 0)
        return inplace;

    updateSummaryRegionValues(inplace,
                              miscSummaryData,
                              regionData);

    if (!substep)
        outputFipresvLogImpl(inplace);

    return inplace;
}

template<class FluidSystem,class Scalar>
void EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
addRftDataToWells(data::Wells& wellDatas, size_t reportStepNum)
{
    const auto& rft_config = schedule_[reportStepNum].rft_config();
    for (const auto& well: schedule_.getWells(reportStepNum)) {

        // don't bother with wells not on this process
        if (isDefunctParallelWell(well.name())) {
            continue;
        }

        //add data infrastructure for shut wells
        if (!wellDatas.count(well.name())) {
            data::Well wellData;

            if (!rft_config.active())
                continue;

            wellData.connections.resize(well.getConnections().size());
            size_t count = 0;
            for (const auto& connection: well.getConnections()) {
                const size_t i = size_t(connection.getI());
                const size_t j = size_t(connection.getJ());
                const size_t k = size_t(connection.getK());

                const size_t index = eclState_.gridDims().getGlobalIndex(i, j, k);
                auto& connectionData = wellData.connections[count];
                connectionData.index = index;
                count++;
            }
            wellDatas.emplace(std::make_pair(well.name(), wellData));
        }

        data::Well& wellData = wellDatas.at(well.name());
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

template<class FluidSystem, class Scalar>
void EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
assignToSolution(data::Solution& sol)
{
    using DataEntry =
        std::tuple<std::string, UnitSystem::measure, const std::vector<Scalar>&>;

    auto doInsert = [&sol](const DataEntry&       entry,
                           const data::TargetType target)
    {
        if (std::get<2>(entry).empty()) {
            return;
        }

        sol.insert(std::get<std::string>(entry),
                   std::get<UnitSystem::measure>(entry),
                   std::move(std::get<2>(entry)),
                   target);
    };

    const auto baseSolutionArrays = std::array {
        DataEntry{"1OVERBG",  UnitSystem::measure::gas_inverse_formation_volume_factor,   invB_[gasPhaseIdx]},
        DataEntry{"1OVERBO",  UnitSystem::measure::oil_inverse_formation_volume_factor,   invB_[oilPhaseIdx]},
        DataEntry{"1OVERBW",  UnitSystem::measure::water_inverse_formation_volume_factor, invB_[waterPhaseIdx]},
        DataEntry{"FLOGASI+", UnitSystem::measure::gas_surface_rate,                      flowsi_[gasCompIdx]},
        DataEntry{"FLOGASJ+", UnitSystem::measure::gas_surface_rate,                      flowsj_[gasCompIdx]},
        DataEntry{"FLOGASK+", UnitSystem::measure::gas_surface_rate,                      flowsk_[gasCompIdx]},
        DataEntry{"FLOOILI+", UnitSystem::measure::liquid_surface_rate,                   flowsi_[oilCompIdx]},
        DataEntry{"FLOOILJ+", UnitSystem::measure::liquid_surface_rate,                   flowsj_[oilCompIdx]},
        DataEntry{"FLOOILK+", UnitSystem::measure::liquid_surface_rate,                   flowsk_[oilCompIdx]},
        DataEntry{"FLOWATI+", UnitSystem::measure::liquid_surface_rate,                   flowsi_[waterCompIdx]},
        DataEntry{"FLOWATJ+", UnitSystem::measure::liquid_surface_rate,                   flowsj_[waterCompIdx]},
        DataEntry{"FLOWATK+", UnitSystem::measure::liquid_surface_rate,                   flowsk_[waterCompIdx]},
        DataEntry{"FLRGASI+", UnitSystem::measure::rate,                                  floresi_[gasCompIdx]},
        DataEntry{"FLRGASJ+", UnitSystem::measure::rate,                                  floresj_[gasCompIdx]},
        DataEntry{"FLRGASK+", UnitSystem::measure::rate,                                  floresk_[gasCompIdx]},
        DataEntry{"FLROILI+", UnitSystem::measure::rate,                                  floresi_[oilCompIdx]},
        DataEntry{"FLROILJ+", UnitSystem::measure::rate,                                  floresj_[oilCompIdx]},
        DataEntry{"FLROILK+", UnitSystem::measure::rate,                                  floresk_[oilCompIdx]},
        DataEntry{"FLRWATI+", UnitSystem::measure::rate,                                  floresi_[waterCompIdx]},
        DataEntry{"FLRWATJ+", UnitSystem::measure::rate,                                  floresj_[waterCompIdx]},
        DataEntry{"FLRWATK+", UnitSystem::measure::rate,                                  floresk_[waterCompIdx]},
        DataEntry{"FOAM",     UnitSystem::measure::identity,                              cFoam_},
        DataEntry{"GASKR",    UnitSystem::measure::identity,                              relativePermeability_[gasPhaseIdx]},
        DataEntry{"GAS_DEN",  UnitSystem::measure::density,                               density_[gasPhaseIdx]},
        DataEntry{"GAS_VISC", UnitSystem::measure::viscosity,                             viscosity_[gasPhaseIdx]},
        DataEntry{"OILKR",    UnitSystem::measure::identity,                              relativePermeability_[oilPhaseIdx]},
        DataEntry{"OIL_DEN",  UnitSystem::measure::density,                               density_[oilPhaseIdx]},
        DataEntry{"OIL_VISC", UnitSystem::measure::viscosity,                             viscosity_[oilPhaseIdx]},
        DataEntry{"PBUB",     UnitSystem::measure::pressure,                              bubblePointPressure_},
        DataEntry{"PCOG",     UnitSystem::measure::pressure,                              pcog_},
        DataEntry{"PCOW",     UnitSystem::measure::pressure,                              pcow_},
        DataEntry{"PDEW",     UnitSystem::measure::pressure,                              dewPointPressure_},
        DataEntry{"POLYMER",  UnitSystem::measure::identity,                              cPolymer_},
        DataEntry{"PPCW",     UnitSystem::measure::pressure,                              ppcw_},
        DataEntry{"PRESROCC", UnitSystem::measure::pressure,                              minimumOilPressure_},
        DataEntry{"PRESSURE", UnitSystem::measure::pressure,                              fluidPressure_},
        DataEntry{"RS",       UnitSystem::measure::gas_oil_ratio,                         rs_},
        DataEntry{"RSSAT",    UnitSystem::measure::gas_oil_ratio,                         gasDissolutionFactor_},
        DataEntry{"RV",       UnitSystem::measure::oil_gas_ratio,                         rv_},
        DataEntry{"RVSAT",    UnitSystem::measure::oil_gas_ratio,                         oilVaporizationFactor_},
        DataEntry{"SALT",     UnitSystem::measure::salinity,                              cSalt_},
        DataEntry{"SOMAX",    UnitSystem::measure::identity,                              soMax_},
        DataEntry{"SSOLVENT", UnitSystem::measure::identity,                              sSol_},
        DataEntry{"SWMAX",    UnitSystem::measure::identity,                              swMax_},
        DataEntry{"WATKR",    UnitSystem::measure::identity,                              relativePermeability_[waterPhaseIdx]},
        DataEntry{"WAT_DEN",  UnitSystem::measure::density,                               density_[waterPhaseIdx]},
        DataEntry{"WAT_VISC", UnitSystem::measure::viscosity,                             viscosity_[waterPhaseIdx]},
    };

    const auto extendedSolutionArrays = std::array {
        DataEntry{"BIOFILM",  UnitSystem::measure::identity,           cBiofilm_},
        DataEntry{"CALCITE",  UnitSystem::measure::identity,           cCalcite_},
        DataEntry{"DRSDTCON", UnitSystem::measure::gas_oil_ratio_rate, drsdtcon_},
        DataEntry{"KRNSW_GO", UnitSystem::measure::identity,           krnSwMdcGo_},
        DataEntry{"KRNSW_OW", UnitSystem::measure::identity,           krnSwMdcOw_},
        DataEntry{"MICROBES", UnitSystem::measure::density,            cMicrobes_},
        DataEntry{"OXYGEN",   UnitSystem::measure::density,            cOxygen_},
        DataEntry{"PCSWM_GO", UnitSystem::measure::identity,           pcSwMdcGo_},
        DataEntry{"PCSWM_OW", UnitSystem::measure::identity,           pcSwMdcOw_},
        DataEntry{"PERMFACT", UnitSystem::measure::identity,           permFact_},
        DataEntry{"PORV_RC",  UnitSystem::measure::identity,           rockCompPorvMultiplier_},
        DataEntry{"PRES_OVB", UnitSystem::measure::pressure,           overburdenPressure_},
        DataEntry{"RSW",      UnitSystem::measure::gas_oil_ratio,      rsw_},
        DataEntry{"RVW",      UnitSystem::measure::oil_gas_ratio,      rvw_},
        DataEntry{"SALTP",    UnitSystem::measure::identity,           pSalt_},
        DataEntry{"SS_X",     UnitSystem::measure::identity,           extboX_},
        DataEntry{"SS_Y",     UnitSystem::measure::identity,           extboY_},
        DataEntry{"SS_Z",     UnitSystem::measure::identity,           extboZ_},
        DataEntry{"STD_CO2",  UnitSystem::measure::identity,           mFracCo2_},
        DataEntry{"STD_GAS",  UnitSystem::measure::identity,           mFracGas_},
        DataEntry{"STD_OIL",  UnitSystem::measure::identity,           mFracOil_},
        DataEntry{"TMULT_RC", UnitSystem::measure::identity,           rockCompTransMultiplier_},
        DataEntry{"UREA",     UnitSystem::measure::density,            cUrea_},
    };

    for (const auto& array : baseSolutionArrays) {
        doInsert(array, data::TargetType::RESTART_SOLUTION);
    }

    for (const auto& array : extendedSolutionArrays) {
        doInsert(array, data::TargetType::RESTART_OPM_EXTENDED);
    }

    if (! this->temperature_.empty())
    {
        sol.insert("TEMP", UnitSystem::measure::temperature,
                   std::move(this->temperature_), data::TargetType::RESTART_SOLUTION);
    }

    if (FluidSystem::phaseIsActive(waterPhaseIdx) &&
        ! this->saturation_[waterPhaseIdx].empty())
    {
        sol.insert("SWAT", UnitSystem::measure::identity,
                   std::move(this->saturation_[waterPhaseIdx]),
                   data::TargetType::RESTART_SOLUTION);
    }

    if (FluidSystem::phaseIsActive(gasPhaseIdx) &&
        ! this->saturation_[gasPhaseIdx].empty())
    {
        sol.insert("SGAS", UnitSystem::measure::identity,
                   std::move(this->saturation_[gasPhaseIdx]),
                   data::TargetType::RESTART_SOLUTION);
    }

    if (eclState_.runspec().co2Storage() && !rsw_.empty()) {
        auto mfrac = std::vector<double>(this->rsw_.size(), 0.0);

        std::transform(this->rsw_.begin(), this->rsw_.end(),
                       this->eclState_.fieldProps().get_int("PVTNUM").begin(),
                       mfrac.begin(),
            [](const auto& rsw, const int pvtReg)
        {
            const auto xwg = FluidSystem::convertRswToXwG(rsw, pvtReg - 1);
            return FluidSystem::convertXwGToxwG(xwg, pvtReg - 1);
        });

        sol.insert("XMFCO2",
                   UnitSystem::measure::identity,
                   std::move(mfrac),
                   data::TargetType::RESTART_OPM_EXTENDED);
    }

    if (eclState_.runspec().co2Storage() && !rvw_.empty()) {
        auto mfrac = std::vector<double>(this->rvw_.size(), 0.0);

        std::transform(this->rvw_.begin(), this->rvw_.end(),
                       this->eclState_.fieldProps().get_int("PVTNUM").begin(),
                       mfrac.begin(),
            [](const auto& rvw, const int pvtReg)
        {
            const auto xgw = FluidSystem::convertRvwToXgW(rvw, pvtReg - 1);
            return FluidSystem::convertXgWToxgW(xgw, pvtReg - 1);
        });

        sol.insert("YMFWAT",
                   UnitSystem::measure::identity,
                   std::move(mfrac),
                   data::TargetType::RESTART_OPM_EXTENDED);
    }

    // Fluid in place
    if (this->outputFipRestart_) {
        for (const auto& phase : Inplace::phases()) {
            if (! this->fip_[phase].empty()) {
                sol.insert(EclString(phase),
                           UnitSystem::measure::volume,
                           this->fip_[phase],
                           data::TargetType::SUMMARY);
            }
        }
    }

    // Tracers
    if (! this->tracerConcentrations_.empty()) {
        const auto& tracers = this->eclState_.tracer();
        for (auto tracerIdx = 0*tracers.size();
             tracerIdx < tracers.size(); ++tracerIdx)
        {
            sol.insert(tracers[tracerIdx].fname(),
                       UnitSystem::measure::identity,
                       std::move(tracerConcentrations_[tracerIdx]),
                       data::TargetType::RESTART_TRACER_SOLUTION);
        }

        // Put tracerConcentrations container into a valid state.  Otherwise
        // we'll move from vectors that have already been moved from if we
        // get here and it's not a restart step.
        this->tracerConcentrations_.clear();
    }
}

template<class FluidSystem,class Scalar>
void EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
setRestart(const data::Solution& sol,
           unsigned elemIdx,
           unsigned globalDofIndex)
{
    Scalar so = 1.0;
    if (!saturation_[waterPhaseIdx].empty() && sol.has("SWAT")) {
        saturation_[waterPhaseIdx][elemIdx] = sol.data("SWAT")[globalDofIndex];
        so -= sol.data("SWAT")[globalDofIndex];
    }
    if (!saturation_[gasPhaseIdx].empty() && sol.has("SGAS")) {
        saturation_[gasPhaseIdx][elemIdx] = sol.data("SGAS")[globalDofIndex];
        so -= sol.data("SGAS")[globalDofIndex];
    }

    if (!sSol_.empty()) {
        // keep the SSOL option for backward compatibility
        // should be removed after 10.2018 release
        if (sol.has("SSOL"))
            sSol_[elemIdx] = sol.data("SSOL")[globalDofIndex];
        else if (sol.has("SSOLVENT"))
            sSol_[elemIdx] = sol.data("SSOLVENT")[globalDofIndex];

        so -= sSol_[elemIdx];
    }

    assert(!saturation_[oilPhaseIdx].empty());
    saturation_[oilPhaseIdx][elemIdx] = so;

    auto assign = [elemIdx, globalDofIndex, &sol](const std::string& name,
                                                  ScalarBuffer& data)

    {
        if (!data.empty() && sol.has(name)) {
            data[elemIdx] = sol.data(name)[globalDofIndex];
        }
    };

    const auto fields = std::array{
        std::pair{"BIOFILM", &cBiofilm_},
        std::pair{"CALCITE",&cCalcite_},
        std::pair{"FOAM", &cFoam_},
        std::pair{"KRNSW_GO", &krnSwMdcGo_},
        std::pair{"KRNSW_OW", &krnSwMdcOw_},
        std::pair{"MICROBES", &cMicrobes_},
        std::pair{"OXYGEN", &cOxygen_},
        std::pair{"PCSWM_GO", &pcSwMdcGo_},
        std::pair{"PCSWM_OW", &pcSwMdcOw_},
        std::pair{"PERMFACT", &permFact_},
        std::pair{"POLYMER", &cPolymer_},
        std::pair{"PPCW", &ppcw_},
        std::pair{"PRESSURE", &fluidPressure_},
        std::pair{"RS", &rs_},
        std::pair{"RSW", &rsw_},
        std::pair{"RV", &rv_},
        std::pair{"RVW", &rvw_},
        std::pair{"SALT", &cSalt_},
        std::pair{"SALTP", &pSalt_},
        std::pair{"SOMAX", &soMax_},
        std::pair{"TEMP", &temperature_},
        std::pair{"UREA", &cUrea_},
    };

    std::for_each(fields.begin(), fields.end(),
                  [&assign](const auto& p)
                  { assign(p.first, *p.second); });
}

template<class FluidSystem,class Scalar>
typename EclGenericOutputBlackoilModule<FluidSystem,Scalar>::ScalarBuffer
EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
regionSum(const ScalarBuffer& property,
          const std::vector<int>& regionId,
          size_t maxNumberOfRegions,
          const Parallel::Communication& comm)
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

        for (size_t i = 0; i < maxNumberOfRegions; ++i)
            totals[i] = comm.sum(totals[i]);

        return totals;
    }

template<class FluidSystem, class Scalar>
void EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
doAllocBuffers(unsigned bufferSize,
               unsigned reportStepNum,
               const bool substep,
               const bool log,
               const bool isRestart,
               const bool vapparsActive,
               const bool enableHysteresis,
               unsigned numTracers,
               unsigned numOutputNnc)
{
    // Output RESTART_OPM_EXTENDED only when explicitly requested by user.
    std::map<std::string, int> rstKeywords = schedule_.rst_keywords(reportStepNum);
    for (auto& [keyword, should_write] : rstKeywords) {
        if (this->isOutputCreationDirective_(keyword)) {
            // 'BASIC', 'FREQ' and similar.  Don't attempt to create
            // cell-based output for these keywords and don't warn about
            // not being able to create such cell-based result vectors.
            should_write = 0;
        }
    }

    if (auto& norst = rstKeywords["NORST"]; norst > 0) {
        // Don't emit diagnostic messages about unsupported 'NORST' key.
        norst = 0;
    }

    this->outputFipRestart_ = false;
    this->computeFip_ = false;

    // Fluid in place
    for (const auto& phase : Inplace::phases()) {
        if (!substep || summaryConfig_.require3DField(EclString(phase))) {
            if (auto& fip = rstKeywords["FIP"]; fip > 0) {
                fip = 0;
                this->outputFipRestart_ = true;
            }

            this->fip_[phase].resize(bufferSize, 0.0);
            this->computeFip_ = true;
        }
        else {
            this->fip_[phase].clear();
        }
    }

    if (!substep ||
        this->summaryConfig_.hasKeyword("FPR") ||
        this->summaryConfig_.hasKeyword("FPRP") ||
        !this->RPRNodes_.empty())
    {
        this->fip_[Inplace::Phase::PoreVolume].resize(bufferSize, 0.0);
        this->dynamicPoreVolume_.resize(bufferSize, 0.0);
        this->hydrocarbonPoreVolume_.resize(bufferSize, 0.0);
        this->pressureTimesPoreVolume_.resize(bufferSize, 0.0);
        this->pressureTimesHydrocarbonVolume_.resize(bufferSize, 0.0);
    }
    else {
        this->dynamicPoreVolume_.clear();
        this->hydrocarbonPoreVolume_.clear();
        this->pressureTimesPoreVolume_.clear();
        this->pressureTimesHydrocarbonVolume_.clear();
    }

    // Well RFT data
    if (!substep) {
        const auto& rft_config = schedule_[reportStepNum].rft_config();
        for (const auto& well: schedule_.getWells(reportStepNum)) {

            // don't bother with wells not on this process
            if (isDefunctParallelWell(well.name())) {
                continue;
            }

            if (!rft_config.active())
                continue;

            for (const auto& connection: well.getConnections()) {
                const size_t i = size_t(connection.getI());
                const size_t j = size_t(connection.getJ());
                const size_t k = size_t(connection.getK());
                const size_t index = eclState_.gridDims().getGlobalIndex(i, j, k);

                if (FluidSystem::phaseIsActive(oilPhaseIdx))
                    oilConnectionPressures_.emplace(std::make_pair(index, 0.0));

                if (FluidSystem::phaseIsActive(waterPhaseIdx))
                    waterConnectionSaturations_.emplace(std::make_pair(index, 0.0));

                if (FluidSystem::phaseIsActive(gasPhaseIdx))
                    gasConnectionSaturations_.emplace(std::make_pair(index, 0.0));
            }
        }
    }

    // Field data should be allocated
    // 1) When we want to restart
    // 2) When it is ask for by the user via restartConfig
    // 3) When it is not a substep
    if (!isRestart && (!schedule_.write_rst_file(reportStepNum) || substep)) {
        return;
    }

    // Always output saturation of active phases
    for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
        if (! FluidSystem::phaseIsActive(phaseIdx)) {
            continue;
        }

        this->saturation_[phaseIdx].resize(bufferSize, 0.0);
    }

    // And oil pressure
    fluidPressure_.resize(bufferSize, 0.0);
    rstKeywords["PRES"] = 0;
    rstKeywords["PRESSURE"] = 0;

    // If TEMP is set in RPTRST we output temperature even if THERMAL
    // is not activated
    if (enableEnergy_ || rstKeywords["TEMP"] > 0) {
        this->temperature_.resize(bufferSize, 0.0);
        rstKeywords["TEMP"] = 0;
    }

    if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
        rstKeywords["SOIL"] = 0;
    }
    if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
        rstKeywords["SGAS"] = 0;
    }
    if (FluidSystem::phaseIsActive(waterPhaseIdx)) {
        rstKeywords["SWAT"] = 0;
    }

    if (FluidSystem::enableDissolvedGas()) {
        rs_.resize(bufferSize, 0.0);
        rstKeywords["RS"] = 0;
    }
    if (FluidSystem::enableDissolvedGasInWater()) {
        rsw_.resize(bufferSize, 0.0);
        rstKeywords["RSW"] = 0;
    }
    if (FluidSystem::enableVaporizedOil()) {
        rv_.resize(bufferSize, 0.0);
        rstKeywords["RV"] = 0;
    }
    if (FluidSystem::enableVaporizedWater()) {
        rvw_.resize(bufferSize, 0.0);
        rstKeywords["RVW"] = 0;
    }

    if (schedule_[reportStepNum].oilvap().drsdtConvective()) {
        drsdtcon_.resize(bufferSize, 0.0);
    }

    if (enableSolvent_) {
        sSol_.resize(bufferSize, 0.0);
    }

    if (enablePolymer_) {
        cPolymer_.resize(bufferSize, 0.0);
    }

    if (enableFoam_) {
        cFoam_.resize(bufferSize, 0.0);
    }

    if (enableBrine_) {
        cSalt_.resize(bufferSize, 0.0);
    }

    if (enableSaltPrecipitation_) {
        pSalt_.resize(bufferSize, 0.0);
        permFact_.resize(bufferSize, 0.0);
    }

    if (enableExtbo_) {
        extboX_.resize(bufferSize, 0.0);
        extboY_.resize(bufferSize, 0.0);
        extboZ_.resize(bufferSize, 0.0);
        mFracOil_.resize(bufferSize, 0.0);
        mFracGas_.resize(bufferSize, 0.0);
        mFracCo2_.resize(bufferSize, 0.0);
    }

    if (enableMICP_) {
        cMicrobes_.resize(bufferSize, 0.0);
        cOxygen_.resize(bufferSize, 0.0);
        cUrea_.resize(bufferSize, 0.0);
        cBiofilm_.resize(bufferSize, 0.0);
        cCalcite_.resize(bufferSize, 0.0);
    }

    if (vapparsActive) {
        soMax_.resize(bufferSize, 0.0);
    }

    if (enableHysteresis) {
        pcSwMdcOw_.resize(bufferSize, 0.0);
        krnSwMdcOw_.resize(bufferSize, 0.0);
        pcSwMdcGo_.resize(bufferSize, 0.0);
        krnSwMdcGo_.resize(bufferSize, 0.0);
    }

    if (eclState_.fieldProps().has_double("SWATINIT")) {
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

    enableFlows_ = false;
    enableFlowsn_ = false;
    if (rstKeywords["FLOWS"] > 0) {
        rstKeywords["FLOWS"] = 0;
        enableFlows_ = true;

        const std::array<int, 3> phaseIdxs = { gasPhaseIdx, oilPhaseIdx, waterPhaseIdx };
        const std::array<int, 3> compIdxs = { gasCompIdx, oilCompIdx, waterCompIdx };
        const auto rstName = std::array { "FLOGASN+", "FLOOILN+", "FLOWATN+" };

        for (unsigned ii = 0; ii < phaseIdxs.size(); ++ii) {
            if (FluidSystem::phaseIsActive(phaseIdxs[ii])) {
                flowsi_[compIdxs[ii]].resize(bufferSize, 0.0);
                flowsj_[compIdxs[ii]].resize(bufferSize, 0.0);
                flowsk_[compIdxs[ii]].resize(bufferSize, 0.0);

                if (numOutputNnc > 0) {
                    enableFlowsn_ = true;

                    flowsn_[compIdxs[ii]].first = rstName[ii];
                    flowsn_[compIdxs[ii]].second.first.resize(numOutputNnc, -1);
                    flowsn_[compIdxs[ii]].second.second.resize(numOutputNnc, 0.0);
                }
            }
        }
    }

    enableFlores_ = false;
    enableFloresn_ = false;
    if (rstKeywords["FLORES"] > 0) {
        rstKeywords["FLORES"] = 0;
        enableFlores_ = true;

        const std::array<int, 3> phaseIdxs = { gasPhaseIdx, oilPhaseIdx, waterPhaseIdx };
        const std::array<int, 3> compIdxs = { gasCompIdx, oilCompIdx, waterCompIdx };
        const auto rstName = std::array{ "FLRGASN+", "FLROILN+", "FLRWATN+" };

        for (unsigned ii = 0; ii < phaseIdxs.size(); ++ii) {
            if (FluidSystem::phaseIsActive(phaseIdxs[ii])) {
                floresi_[compIdxs[ii]].resize(bufferSize, 0.0);
                floresj_[compIdxs[ii]].resize(bufferSize, 0.0);
                floresk_[compIdxs[ii]].resize(bufferSize, 0.0);

                if (numOutputNnc > 0) {
                    enableFloresn_ = true;

                    floresn_[compIdxs[ii]].first = rstName[ii];
                    floresn_[compIdxs[ii]].second.first.resize(numOutputNnc, -1);
                    floresn_[compIdxs[ii]].second.second.resize(numOutputNnc, 0.0);
                }
            }
        }
    }

    if (auto& den = rstKeywords["DEN"]; den > 0) {
        den = 0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }

            this->density_[phaseIdx].resize(bufferSize, 0.0);
        }
    }

    if (auto& deng = rstKeywords["DENG"]; (deng > 0) && FluidSystem::phaseIsActive(gasPhaseIdx)) {
        deng = 0;
        this->density_[gasPhaseIdx].resize(bufferSize, 0.0);
    }

    if (auto& deno = rstKeywords["DENO"]; (deno > 0) && FluidSystem::phaseIsActive(oilPhaseIdx)) {
        deno = 0;
        this->density_[oilPhaseIdx].resize(bufferSize, 0.0);
    }

    if (auto& denw = rstKeywords["DENW"]; (denw > 0) && FluidSystem::phaseIsActive(waterPhaseIdx)) {
        denw = 0;
        this->density_[waterPhaseIdx].resize(bufferSize, 0.0);
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

    if (FluidSystem::phaseIsActive(oilPhaseIdx) && FluidSystem::phaseIsActive(waterPhaseIdx) && rstKeywords["PCOW"] > 0) {
        rstKeywords["PCOW"] = 0;
        pcow_.resize(bufferSize, 0.0);
    }
    if (FluidSystem::phaseIsActive(oilPhaseIdx) && FluidSystem::phaseIsActive(gasPhaseIdx) && rstKeywords["PCOG"] > 0) {
        rstKeywords["PCOG"] = 0;
        pcog_.resize(bufferSize, 0.0);
    }

    if (rstKeywords["PBPD"] > 0)  {
        rstKeywords["PBPD"] = 0;
        bubblePointPressure_.resize(bufferSize, 0.0);
        dewPointPressure_.resize(bufferSize, 0.0);
    }

    // tracers
    if (numTracers > 0) {
        tracerConcentrations_.resize(numTracers);
        for (unsigned tracerIdx = 0; tracerIdx < numTracers; ++tracerIdx)
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
                OpmLog::warning("Unhandled output keyword", logstring);
            }
        }
    }

    failedCellsPb_.clear();
    failedCellsPd_.clear();

    // Not supported in flow legacy
    if (false) {
        saturatedOilFormationVolumeFactor_.resize(bufferSize, 0.0);
    }

    if (false) {
        oilSaturationPressure_.resize(bufferSize, 0.0);
    }
}

template<class FluidSystem, class Scalar>
void EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
fipUnitConvert_(std::unordered_map<Inplace::Phase, Scalar>& fip) const
{
    const UnitSystem& units = eclState_.getUnits();
    using M = UnitSystem::measure;
    const auto unit_map = std::unordered_map<Inplace::Phase, M> {
        {Inplace::Phase::WATER,             M::liquid_surface_volume},
        {Inplace::Phase::OIL,               M::liquid_surface_volume},
        {Inplace::Phase::OilInLiquidPhase,  M::liquid_surface_volume},
        {Inplace::Phase::OilInGasPhase,     M::liquid_surface_volume},
        {Inplace::Phase::GAS,               M::gas_surface_volume},
        {Inplace::Phase::GasInLiquidPhase,  M::gas_surface_volume},
        {Inplace::Phase::GasInGasPhase,     M::gas_surface_volume},
        {Inplace::Phase::PoreVolume,        M::volume},
        {Inplace::Phase::DynamicPoreVolume, M::volume},
        {Inplace::Phase::WaterResVolume,    M::volume},
        {Inplace::Phase::OilResVolume,      M::volume},
        {Inplace::Phase::GasResVolume,      M::volume},
        {Inplace::Phase::SALT,              M::mass},
        {Inplace::Phase::CO2InWaterPhase,   M::moles},
        {Inplace::Phase::CO2InGasPhaseInMob,M::moles},
        {Inplace::Phase::CO2InGasPhaseMob,  M::moles},
        {Inplace::Phase::WaterInWaterPhase, M::liquid_surface_volume},
        {Inplace::Phase::WaterInGasPhase,   M::liquid_surface_volume},
    };

    for (auto& [phase, value] : fip) {
        auto unitPos = unit_map.find(phase);
        if (unitPos != unit_map.end()) {
            value = units.from_si(unitPos->second, value);
        }
    }
}

template<class FluidSystem, class Scalar>
void EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
pressureUnitConvert_(Scalar& pav) const
{
    pav = this->eclState_.getUnits()
        .from_si(UnitSystem::measure::pressure, pav);
}

template<class FluidSystem, class Scalar>
void EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
outputRegionFluidInPlace_(std::unordered_map<Inplace::Phase, Scalar> oip,
                      std::unordered_map<Inplace::Phase, Scalar> cip,
                      const Scalar& pav, const int reg) const
{
    if (forceDisableFipOutput_)
        return;

    // don't output FIPNUM report if the region has no porv.
    if (! (cip[Inplace::Phase::PoreVolume] > Scalar{0}))
        return;

    const UnitSystem& units = eclState_.getUnits();
    std::ostringstream ss;

    ss << '\n';
    if (reg == 0) {
        ss << "Field total";
    }
    else {
        ss << "FIPNUM report region " << reg;
    }

    ss << " pressure dependent pore volume = "
       << std::fixed << std::setprecision(0)
       << cip[Inplace::Phase::DynamicPoreVolume] << ' '
       << units.name(UnitSystem::measure::volume) << "\n\n";

    if (reg == 0) {
        ss << "                                                  ===================================================\n"
           << "                                                  :                   Field Totals                  :\n";
    }
    else {
        ss << "                                                  ===================================================\n"
           << "                                                  :        FIPNUM report region  "
           << std::setw(2) << reg << "                 :\n";
    }
    if (units.getType() == UnitSystem::UnitType::UNIT_TYPE_METRIC) {
        ss << "                                                  :      PAV  =" << std::setw(14) << pav << " BARSA                 :\n"
           << std::fixed << std::setprecision(0)
           << "                                                  :      PORV =" << std::setw(14) << cip[Inplace::Phase::PoreVolume] << "   RM3                 :\n";
        if (!reg) {
            ss << "                                                  : Pressure is weighted by hydrocarbon pore volume :\n"
               << "                                                  : Porv volumes are taken at reference conditions  :\n";
        }
        ss << "                         :--------------- Oil    SM3 ---------------:-- Wat    SM3 --:--------------- Gas    SM3 ---------------:\n";
    }
    if (units.getType() == UnitSystem::UnitType::UNIT_TYPE_FIELD) {
        ss << "                                                  :      PAV  =" << std::setw(14) << pav << "  PSIA                 :\n"
           << std::fixed << std::setprecision(0)
           << "                                                  :      PORV =" << std::setw(14) << cip[Inplace::Phase::PoreVolume] << "   RB                  :\n";
        if (!reg) {
            ss << "                                                  : Pressure is weighted by hydrocarbon pore volume :\n"
               << "                                                  : Pore volumes are taken at reference conditions  :\n";
        }
        ss << "                         :--------------- Oil    STB ---------------:-- Wat    STB --:--------------- Gas   MSCF ---------------:\n";
    }
    ss << "                         :      Liquid        Vapour        Total   :      Total     :      Free        Dissolved       Total   :" << "\n"
       << ":------------------------:------------------------------------------:----------------:------------------------------------------:" << "\n"
       << ":Currently   in place    :" << std::setw(14) << cip[Inplace::Phase::OilInLiquidPhase] << std::setw(14) << cip[Inplace::Phase::OilInGasPhase] << std::setw(14) << cip[Inplace::Phase::OIL] << ":"
       << std::setw(13) << cip[Inplace::Phase::WATER] << "   :" << std::setw(14) << (cip[Inplace::Phase::GasInGasPhase]) << std::setw(14) << cip[Inplace::Phase::GasInLiquidPhase] << std::setw(14) << cip[Inplace::Phase::GAS] << ":\n"
       << ":------------------------:------------------------------------------:----------------:------------------------------------------:\n"
       << ":Originally  in place    :" << std::setw(14) << oip[Inplace::Phase::OilInLiquidPhase] << std::setw(14) << oip[Inplace::Phase::OilInGasPhase] << std::setw(14) << oip[Inplace::Phase::OIL] << ":"
       << std::setw(13) << oip[Inplace::Phase::WATER] << "   :" << std::setw(14) << oip[Inplace::Phase::GasInGasPhase] << std::setw(14) << oip[Inplace::Phase::GasInLiquidPhase] << std::setw(14) << oip[Inplace::Phase::GAS] << ":\n"
       << ":========================:==========================================:================:==========================================:\n";
    OpmLog::note(ss.str());
}

template<class FluidSystem, class Scalar>
void EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
outputResvFluidInPlace_(std::unordered_map<Inplace::Phase, Scalar> cipr, const int reg) const
{
    if (forceDisableFipresvOutput_)
        return;

    // don't output FIPNUM report if the region has no porv.
    if (cipr[Inplace::Phase::PoreVolume] == 0)
        return;
    const UnitSystem& units = eclState_.getUnits();
    std::ostringstream ss;

    if (reg == 0) {
        ss << "                                                     ===================================\n";
        if (units.getType() == UnitSystem::UnitType::UNIT_TYPE_METRIC) {
            ss << "                                                     :  RESERVOIR VOLUMES      M3      :\n";
        }
        if (units.getType() == UnitSystem::UnitType::UNIT_TYPE_FIELD) {
            ss << "                                                     :  RESERVOIR VOLUMES      RB      :\n";
        }
        ss << ":---------:---------------:---------------:---------------:---------------:---------------:\n"
           << ": REGION  :  TOTAL PORE   :  PORE VOLUME  :  PORE VOLUME  : PORE VOLUME   :  PORE VOLUME  :\n"
           << ":         :   VOLUME      :  CONTAINING   :  CONTAINING   : CONTAINING    :  CONTAINING   :\n"
           << ":         :               :     OIL       :    WATER      :    GAS        :  HYDRO-CARBON :\n"
           << ":---------:---------------:---------------:---------------:---------------:---------------\n";
    }
    else {
        ss << std::right << std::fixed << std::setprecision(0) << ":" << std::setw (9) <<  reg << ":" << std::setw(15) << cipr[Inplace::Phase::DynamicPoreVolume] << ":" << std::setw(15) << cipr[Inplace::Phase::OilResVolume] << ":" << std::setw(15) << cipr[Inplace::Phase::WaterResVolume] << ":" << std::setw(15) << cipr[Inplace::Phase::GasResVolume] << ":" << std::setw(15) << cipr[Inplace::Phase::OilResVolume] + cipr[Inplace::Phase::GasResVolume] << ":\n"
        << ":---------:---------------:---------------:---------------:---------------:---------------:\n";
    }
    OpmLog::note(ss.str());
}

template<class FluidSystem, class Scalar>
void EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
outputProductionReport_(const ScalarBuffer& wellProd,
                    const StringBuffer& wellProdNames,
                    const bool forceDisableProdOutput)
{
    if (forceDisableProdOutput)
        return;

    const UnitSystem& units = eclState_.getUnits();
    std::ostringstream ss;
    if (wellProdNames[WellProdDataType::WellName].empty()) {
        ss << "======================================================= PRODUCTION REPORT =======================================================\n"//=================== \n"
           << ":  WELL  :  LOCATION :CTRL:    OIL    :   WATER   :    GAS    :   FLUID   :   WATER   : GAS/OIL  :  WAT/GAS   : BHP OR : THP OR :\n"// STEADY-ST PI       :\n"
           << ":  NAME  :  (I,J,K)  :MODE:    RATE   :   RATE    :    RATE   :  RES.VOL. :    CUT    :  RATIO   :   RATIO    : CON.PR.: BLK.PR.:\n";// OR POTN OF PREF. PH:\n";
        if (units.getType() == UnitSystem::UnitType::UNIT_TYPE_METRIC) {
                        ss << ":        :           :    :  SCM/DAY  :  SCM/DAY  :  SCM/DAY  :  RCM/DAY  :  SCM/SCM  :  SCM/SCM :  SCM/SCM   :  BARSA :  BARSA :\n";//                    :\n";
                    }
                    if (units.getType() == UnitSystem::UnitType::UNIT_TYPE_FIELD) {
                        ss << ":        :           :    :  STB/DAY  :  STB/DAY  :  MSCF/DAY :  RB/DAY   :           : MSCF/STB :  STB/MSCF  :  PSIA  :  PSIA  :\n";//                    :\n";
                    }
                if (units.getType() == UnitSystem::UnitType::UNIT_TYPE_LAB) {
                            ss << ":        :           :    :  SCC/HR   :  SCC/HR   :  SCC/HR   :    RCC    :  SCC/SCC  :  SCC/SCC :  SCC/SCC   :  ATMA  :  ATMA  :\n";//                    :\n";
                    }
            ss << "=================================================================================================================================\n";//=================== \n";
    }
    else {
        if (wellProd[WellProdDataType::WellLocationi] < 1) {
            ss << std::right << std::fixed << ":" << std::setw (8) << wellProdNames[WellProdDataType::WellName] << ":" << std::setprecision(0) << std::setw(11) << "" << ":" << std::setw(4) << wellProdNames[WellProdDataType::CTRLMode] << ":" << std::setprecision(1) << std::setw(11) << wellProd[WellProdDataType::OilRate] << ":" << std::setw(11) << wellProd[WellProdDataType::WaterRate] << ":" <<  std::setw(11)<< wellProd[WellProdDataType::GasRate] << ":" <<  std::setw(11) << wellProd[WellProdDataType::FluidResVol] << std::setprecision(3) << ":" <<  std::setw(11) << wellProd[WellProdDataType::WaterCut] << std::setprecision(2) << ":" <<  std::setw(10) << wellProd[WellProdDataType::GasOilRatio] << std::setprecision(4) << ":" <<  std::setw(12) << wellProd[WellProdDataType::WatGasRatio] << std::setprecision(1) << ":" <<  std::setw(8) << "" << ":" <<  std::setw(8) << "" << ": \n";
        }
        else {
            ss << std::right << std::fixed << ":" << std::setw (8) << wellProdNames[WellProdDataType::WellName] << ":" << std::setprecision(0) << std::setw(5) << wellProd[WellProdDataType::WellLocationi] << "," << std::setw(5) << wellProd[WellProdDataType::WellLocationj] << ":" << std::setw(4) << wellProdNames[WellProdDataType::CTRLMode] << ":" << std::setprecision(1) << std::setw(11) << wellProd[WellProdDataType::OilRate] << ":" << std::setw(11) << wellProd[WellProdDataType::WaterRate] << ":" <<  std::setw(11)<< wellProd[WellProdDataType::GasRate] << ":" <<  std::setw(11) << wellProd[WellProdDataType::FluidResVol] << std::setprecision(3) << ":" <<  std::setw(11) << wellProd[WellProdDataType::WaterCut] << std::setprecision(2) << ":" <<  std::setw(10) << wellProd[WellProdDataType::GasOilRatio] << std::setprecision(4) << ":" <<  std::setw(12) << wellProd[WellProdDataType::WatGasRatio] << std::setprecision(1) << ":" <<  std::setw(8) << wellProd[WellProdDataType::BHP] << ":" <<  std::setw(8) << wellProd[WellProdDataType::THP] << ": \n";
        }
        ss << ":"<< std::setfill ('-') << std::setw (9) << ":" << std::setfill ('-') << std::setw (12) << ":" << std::setfill ('-') << std::setw (5) << ":" << std::setfill ('-') << std::setw (12) << ":" << std::setfill ('-') << std::setw (12) << ":" << std::setfill ('-') << std::setw (12) << ":" << std::setfill ('-') << std::setw (12) << ":" << std::setfill ('-') << std::setw (12) << ":" << std::setfill ('-') << std::setw (11) << ":" << std::setfill ('-') << std::setw (13) << ":" << std::setfill ('-') << std::setw (9) << ":" << std::setfill ('-') << std::setw (9) << ":" << "\n";
    }
    OpmLog::note(ss.str());
}

template<class FluidSystem, class Scalar>
void EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
outputInjectionReport_(const ScalarBuffer& wellInj,
                   const StringBuffer& wellInjNames,
                   const bool forceDisableInjOutput)
{
    if (forceDisableInjOutput)
        return;

    const UnitSystem& units = eclState_.getUnits();
    std::ostringstream ss;
    if (wellInjNames[WellInjDataType::WellName].empty()) {
        ss << "=================================================== INJECTION REPORT ========================================\n"//===================== \n"
           << ":  WELL  :  LOCATION : CTRL : CTRL : CTRL :    OIL    :   WATER   :    GAS    :   FLUID   : BHP OR : THP OR :\n"// STEADY-ST II       :\n"
           << ":  NAME  :  (I,J,K)  : MODE : MODE : MODE :    RATE   :   RATE    :    RATE   :  RES.VOL. : CON.PR.: BLK.PR.:\n";// OR POTENTIAL       :\n";
        if (units.getType() == UnitSystem::UnitType::UNIT_TYPE_METRIC) {
                            ss << ":        :           : OIL  : WAT  : GAS  :  SCM/DAY  :  SCM/DAY  :  SCM/DAY  :  RCM/DAY  :  BARSA :  BARSA :\n";//                    :\n";
                    }
                    if (units.getType() == UnitSystem::UnitType::UNIT_TYPE_FIELD) {
                            ss << ":        :           : OIL  : WAT  : GAS  :  STB/DAY  :  STB/DAY  :  MSCF/DAY :  RB/DAY   :  PSIA  :  PSIA  :\n";//                    :\n";
                    }
                    if (units.getType() == UnitSystem::UnitType::UNIT_TYPE_LAB) {
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
    OpmLog::note(ss.str());
}

template<class FluidSystem,class Scalar>
bool EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
isOutputCreationDirective_(const std::string& keyword)
{
    return (keyword == "BASIC") || (keyword == "FREQ")
        || (keyword == "RESTART")                        // From RPTSCHED
        || (keyword == "SAVE")  || (keyword == "SFREQ"); // Not really supported
}

template<class FluidSystem, class Scalar>
Scalar EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
pressureAverage_(const Scalar& pressurePvHydrocarbon,
                 const Scalar& pvHydrocarbon,
                 const Scalar& pressurePv,
                 const Scalar& pv,
                 const bool    hydrocarbon)
{
    if (hydrocarbon && (pvHydrocarbon > 1e-10))
        return pressurePvHydrocarbon / pvHydrocarbon;

    return pressurePv / pv;
}

template<class FluidSystem,class Scalar>
typename EclGenericOutputBlackoilModule<FluidSystem,Scalar>::ScalarBuffer
EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
pressureAverage_(const ScalarBuffer& pressurePvHydrocarbon,
                 const ScalarBuffer& pvHydrocarbon,
                 const ScalarBuffer& pressurePv,
                 const ScalarBuffer& pv,
                 const bool hydrocarbon)
{
    const std::size_t size = pressurePvHydrocarbon.size();
    assert(pvHydrocarbon.size() == size);
    assert(pressurePv.size() == size);
    assert(pv.size() == size);

    ScalarBuffer fraction(size, 0.0);
    for (std::size_t i = 0; i < size; ++i) {
        fraction[i] = pressureAverage_(pressurePvHydrocarbon[i],
                                       pvHydrocarbon[i],
                                       pressurePv[i],
                                       pv[i],
                                       hydrocarbon);
    }

    return fraction;
}

namespace {
    template <typename IndexVector, typename IJKString>
    void logUniqueFailedCells(const std::string& messageTag,
                              std::string_view   prefix,
                              const std::size_t  maxNumCellsFaillog,
                              IndexVector&&      cells,
                              IJKString&&        ijkString)
    {
        if (cells.empty()) {
            return;
        }

        std::sort(cells.begin(), cells.end());
        auto u = std::unique(cells.begin(), cells.end());

        const auto numFailed = static_cast<std::size_t>
            (std::distance(cells.begin(), u));

        std::ostringstream errlog;
        errlog << prefix << " failed for " << numFailed << " cell"
               << ((numFailed != std::size_t{1}) ? "s" : "")
               << " [" << ijkString(cells[0]);

        const auto maxElems = std::min(maxNumCellsFaillog, numFailed);
        for (auto i = 1 + 0*maxElems; i < maxElems; ++i) {
            errlog << ", " << ijkString(cells[i]);
        }

        if (numFailed > maxNumCellsFaillog) {
            errlog << ", ...";
        }

        errlog << ']';

        OpmLog::warning(messageTag, errlog.str());
    }
} // Namespace anonymous

template<class FluidSystem,class Scalar>
void EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
outputErrorLog(const Parallel::Communication& comm) const
{
    const auto root = 0;
    auto globalFailedCellsPbub = gatherv(this->failedCellsPb_, comm, root);
    auto globalFailedCellsPdew = gatherv(this->failedCellsPd_, comm, root);

    if (std::empty(std::get<0>(globalFailedCellsPbub)) &&
        std::empty(std::get<0>(globalFailedCellsPdew)))
    {
        return;
    }

    auto ijkString = [this](const std::size_t globalIndex)
    {
        const auto ijk = this->eclState_.gridDims().getIJK(globalIndex);

        return fmt::format("({},{},{})", ijk[0] + 1, ijk[1] + 1, ijk[2] + 1);
    };

    const auto maxNumCellsFaillog = static_cast<std::size_t>(20);

    logUniqueFailedCells("Bubble point numerical problem",
                         "Finding the bubble point pressure",
                         maxNumCellsFaillog,
                         std::get<0>(std::move(globalFailedCellsPbub)),
                         ijkString);

    logUniqueFailedCells("Dew point numerical problem",
                         "Finding the dew point pressure",
                         maxNumCellsFaillog,
                         std::get<0>(std::move(globalFailedCellsPdew)),
                         ijkString);
}

template<class FluidSystem,class Scalar>
void EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
outputFipLogImpl(const Inplace& inplace) const
{
    {
        Scalar fieldHydroCarbonPoreVolumeAveragedPressure = pressureAverage_(inplace.get(Inplace::Phase::PressureHydroCarbonPV),
                                                                             inplace.get(Inplace::Phase::HydroCarbonPV),
                                                                             inplace.get(Inplace::Phase::PressurePV),
                                                                             inplace.get(Inplace::Phase::DynamicPoreVolume),
                                                                             true);

        std::unordered_map<Inplace::Phase, Scalar> initial_values;
        std::unordered_map<Inplace::Phase, Scalar> current_values;

        for (const auto& phase : Inplace::phases()) {
            initial_values[phase] = this->initialInplace_->get(phase);
            current_values[phase] = inplace.get(phase);
        }

        current_values[Inplace::Phase::DynamicPoreVolume] =
            inplace.get(Inplace::Phase::DynamicPoreVolume);

        fipUnitConvert_(initial_values);
        fipUnitConvert_(current_values);

        pressureUnitConvert_(fieldHydroCarbonPoreVolumeAveragedPressure);
        outputRegionFluidInPlace_(std::move(initial_values),
                                  std::move(current_values),
                                  fieldHydroCarbonPoreVolumeAveragedPressure);
    }

    for (size_t reg = 1; reg <= inplace.max_region("FIPNUM"); ++reg) {
        std::unordered_map<Inplace::Phase, Scalar> initial_values;
        std::unordered_map<Inplace::Phase, Scalar> current_values;

        for (const auto& phase : Inplace::phases()) {
            initial_values[phase] = this->initialInplace_->get("FIPNUM", phase, reg);
            current_values[phase] = inplace.get("FIPNUM", phase, reg);
        }

        current_values[Inplace::Phase::DynamicPoreVolume] =
            inplace.get("FIPNUM", Inplace::Phase::DynamicPoreVolume, reg);

        fipUnitConvert_(initial_values);
        fipUnitConvert_(current_values);

        Scalar regHydroCarbonPoreVolumeAveragedPressure
                = pressureAverage_(inplace.get("FIPNUM", Inplace::Phase::PressureHydroCarbonPV, reg),
                                   inplace.get("FIPNUM", Inplace::Phase::HydroCarbonPV, reg),
                                   inplace.get("FIPNUM", Inplace::Phase::PressurePV, reg),
                                   inplace.get("FIPNUM", Inplace::Phase::DynamicPoreVolume, reg),
                                   true);
        pressureUnitConvert_(regHydroCarbonPoreVolumeAveragedPressure);
        outputRegionFluidInPlace_(std::move(initial_values),
                                  std::move(current_values),
                                  regHydroCarbonPoreVolumeAveragedPressure, reg);
    }
}

template<class FluidSystem,class Scalar>
void EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
outputFipresvLogImpl(const Inplace& inplace) const
{
    {
        std::unordered_map<Inplace::Phase, Scalar> current_values;

        for (const auto& phase : Inplace::phases()) {
            current_values[phase] = inplace.get(phase);
        }
        fipUnitConvert_(current_values);
        outputResvFluidInPlace_(current_values);
    }

    for (size_t reg = 1; reg <= inplace.max_region("FIPNUM"); ++reg) {
        std::unordered_map<Inplace::Phase, Scalar> current_values;

        for (const auto& phase : Inplace::phases()) {
            current_values[phase] = inplace.get("FIPNUM", phase, reg);
        }
        current_values[Inplace::Phase::DynamicPoreVolume] =
            inplace.get("FIPNUM", Inplace::Phase::DynamicPoreVolume, reg);

        fipUnitConvert_(current_values);
        outputResvFluidInPlace_(current_values, reg);
    }
}

template<class FluidSystem,class Scalar>
int EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
regionMax(const std::vector<int>& region,
          const Parallel::Communication& comm)
{
    const auto max_value = region.empty() ? 0 : *std::max_element(region.begin(), region.end());
    return comm.max(max_value);
}

template<class FluidSystem,class Scalar>
void EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
update(Inplace& inplace,
       const std::string& region_name,
       const Inplace::Phase phase,
       const std::size_t ntFip,
       const ScalarBuffer& values)
{
    double sum = 0.0;
    for (std::size_t region_number = 0; region_number < ntFip; ++region_number) {
        const auto rval = static_cast<double>(values[region_number]);
        inplace.add(region_name, phase, region_number + 1, rval);
        sum += rval;
    }
    inplace.add(phase, sum);
}

template<class FluidSystem,class Scalar>
void EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
makeRegionSum(Inplace& inplace,
              const std::string& region_name,
              const Parallel::Communication& comm) const
{
    const auto& region = this->regions_.at(region_name);
    const std::size_t ntFip = this->regionMax(region, comm);

    auto update_inplace =
        [&inplace, &region, &region_name, &comm, ntFip, this]
        (const Inplace::Phase       phase,
         const std::vector<Scalar>& value)
    {
        update(inplace, region_name, phase, ntFip,
               this->regionSum(value, region, ntFip, comm));
    };

    update_inplace(Inplace::Phase::PressurePV,
                   this->pressureTimesPoreVolume_);

    update_inplace(Inplace::Phase::HydroCarbonPV,
                   this->hydrocarbonPoreVolume_);

    update_inplace(Inplace::Phase::PressureHydroCarbonPV,
                   this->pressureTimesHydrocarbonVolume_);

    update_inplace(Inplace::Phase::DynamicPoreVolume,
                   this->dynamicPoreVolume_);

    for (const auto& phase : Inplace::phases()) {
        auto fipPos = this->fip_.find(phase);
        if (fipPos != this->fip_.end()) {
            update_inplace(phase, fipPos->second);
        }
    }
}

template<class FluidSystem,class Scalar>
Inplace EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
accumulateRegionSums(const Parallel::Communication& comm)
{
    Inplace inplace;

    for (const auto& region : this->regions_) {
        makeRegionSum(inplace, region.first, comm);
    }

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
    if (!this->initialInplace_.has_value())
        this->initialInplace_ = inplace;
    return inplace;
}

template<class FluidSystem,class Scalar>
Scalar EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
sum(const ScalarBuffer& v)
{
    return std::accumulate(v.begin(), v.end(), Scalar{0});
}

template<class FluidSystem,class Scalar>
void EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
updateSummaryRegionValues(const Inplace& inplace,
                          std::map<std::string, double>& miscSummaryData,
                          std::map<std::string, std::vector<double>>& regionData) const
{
    // The field summary vectors should only use the FIPNUM based region sum.
    {
        for (const auto& phase : Inplace::phases()) {
            const std::string key = "F" + EclString(phase);
            if (this->summaryConfig_.hasKeyword(key)) {
                miscSummaryData[key] = inplace.get(phase);
            }
        }

        if (this->summaryConfig_.hasKeyword("FOE") && this->initialInplace_) {
            miscSummaryData["FOE"] = (this->initialInplace_.value().get(Inplace::Phase::OIL) - inplace.get(Inplace::Phase::OIL))
                / this->initialInplace_.value().get(Inplace::Phase::OIL);
        }

        if (this->summaryConfig_.hasKeyword("FPR")) {
            miscSummaryData["FPR"] =
                pressureAverage_(inplace.get(Inplace::Phase::PressureHydroCarbonPV),
                                 inplace.get(Inplace::Phase::HydroCarbonPV),
                                 inplace.get(Inplace::Phase::PressurePV),
                                 inplace.get(Inplace::Phase::DynamicPoreVolume),
                                 true);
        }

        if (this->summaryConfig_.hasKeyword("FPRP")) {
            miscSummaryData["FPRP"] =
                pressureAverage_(inplace.get(Inplace::Phase::PressureHydroCarbonPV),
                                 inplace.get(Inplace::Phase::HydroCarbonPV),
                                 inplace.get(Inplace::Phase::PressurePV),
                                 inplace.get(Inplace::Phase::DynamicPoreVolume),
                                 false);
        }

    }

    // The region summary vectors should loop through the FIPxxx regions to
    // support the RPR__xxx summary keywords.
    {
        auto get_vector = [&inplace]
            (const auto&          node_,
             const Inplace::Phase phase_)
        {
            return inplace.get_vector(node_.fip_region(), phase_);
        };

        for (const auto& phase : Inplace::phases()) {
            for (const auto& node : this->regionNodes_.at(phase))
                regionData[node.keyword()] = get_vector(node, phase);
        }

        for (const auto& node : this->RPRNodes_) {
            regionData[node.keyword()] =
            pressureAverage_(get_vector(node, Inplace::Phase::PressureHydroCarbonPV),
                             get_vector(node, Inplace::Phase::HydroCarbonPV),
                             get_vector(node, Inplace::Phase::PressurePV),
                             get_vector(node, Inplace::Phase::DynamicPoreVolume),
                             true);
        }

        for (const auto& node : this->RPRPNodes_) {
            regionData[node.keyword()] =
            pressureAverage_(get_vector(node, Inplace::Phase::PressureHydroCarbonPV),
                             get_vector(node, Inplace::Phase::HydroCarbonPV),
                             get_vector(node, Inplace::Phase::PressurePV),
                             get_vector(node, Inplace::Phase::DynamicPoreVolume),
                             false);
        }
    }
}

template<class FluidSystem,class Scalar>
void EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
setupBlockData(std::function<bool(int)> isCartIdxOnThisRank)
{
    for (const auto& node : summaryConfig_) {
        if ((node.category() == SummaryConfigNode::Category::Block) &&
            isCartIdxOnThisRank(node.number() - 1))
        {
            this->blockData_.emplace(std::piecewise_construct,
                                     std::forward_as_tuple(node.keyword(),
                                                           node.number()),
                                     std::forward_as_tuple(0.0));
        }
    }
}

template class EclGenericOutputBlackoilModule<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,double>;

} // namespace Opm
