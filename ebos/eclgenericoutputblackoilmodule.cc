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

#include <opm/simulators/utils/PressureAverage.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
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
                           bool enableMech,
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
    , enableMech_(enableMech)
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
    for (const auto& region : fp.fip_regions()) {
        this->regions_[region] = fp.get_int(region);
    }

    this->RPRNodes_  = summaryConfig_.keywords("RPR*");
    this->RPRPNodes_ = summaryConfig_.keywords("RPRP*");

    for (const auto& phase : Inplace::phases()) {
        std::string key_pattern = "R" + EclString(phase) + "*";
        this->regionNodes_[phase] = summaryConfig_.keywords(key_pattern);
    }

    // Check for any BFLOW[I|J|K] summary keys
    blockFlows_ = summaryConfig_.keywords("BFLOW*").size() > 0;

    // Check if FLORES/FLOWS is set in any RPTRST in the schedule
    anyFlores_ = false;     // Used for the initialization of the sparse table
    anyFlows_ = blockFlows_;
    enableFlores_ = false;  // Used for the output of i+, j+, k+
    enableFloresn_ = false; // Used for the special case of nnc
    enableFlows_ = false;
    enableFlowsn_ = false;

    for (const auto& block : this->schedule_) {
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
outputTimeStamp(const std::string& lbl, double elapsed, int rstep, boost::posix_time::ptime currentDate)
{
    logOutput_.timeStamp(lbl, elapsed, rstep, currentDate);
}

template<class FluidSystem, class Scalar>
void EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
outputCumLog(std::size_t reportStepNum)
{
    logOutput_.cumulative(reportStepNum,
                          [this](const std::string& name)
                          { return this->isDefunctParallelWell(name); });
}

template<class FluidSystem,class Scalar>
void EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
outputProdLog(std::size_t reportStepNum)
{
    logOutput_.production(reportStepNum,
                          [this](const std::string& name)
                          { return this->isDefunctParallelWell(name); });
}

template<class FluidSystem,class Scalar>
void EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
outputInjLog(std::size_t reportStepNum)
{
    logOutput_.injection(reportStepNum,
                         [this](const std::string& name)
                         { return this->isDefunctParallelWell(name); });
}


template<class FluidSystem,class Scalar>
Inplace EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
calc_inplace(std::map<std::string, double>& miscSummaryData,
             std::map<std::string, std::vector<double>>& regionData,
             const Parallel::Communication& comm)
{
    auto inplace = this->accumulateRegionSums(comm);
    
    if (comm.rank() != 0)
        return inplace;

    updateSummaryRegionValues(inplace,
                              miscSummaryData,
                              regionData);

    
    return inplace;
}


template<class FluidSystem,class Scalar>
void EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
outputFipAndResvLog(const Inplace& inplace,
                         const std::size_t reportStepNum,
                         double elapsed,
                         boost::posix_time::ptime currentDate,
                         const bool substep,
                         const Parallel::Communication& comm)
{

    if (comm.rank() != 0)
        return;


    // For report step 0 we use the RPTSOL config, else derive from RPTSCHED
    std::unique_ptr<FIPConfig> fipSched;
    if (reportStepNum != 0) {
        const auto& rpt = this->schedule_[reportStepNum].rpt_config.get();
        fipSched = std::make_unique<FIPConfig>(rpt);
    }
    const FIPConfig& fipc = reportStepNum == 0 ? this->eclState_.getEclipseConfig().fip()
                                               : *fipSched;

    if (!substep && !forceDisableFipOutput_ && fipc.output(FIPConfig::OutputField::FIELD)) {

        logOutput_.timeStamp("BALANCE", elapsed, reportStepNum, currentDate);

        logOutput_.fip(inplace, this->initialInplace(), "");  
        
        if (fipc.output(FIPConfig::OutputField::FIPNUM)) { 
            logOutput_.fip(inplace, this->initialInplace(), "FIPNUM");    
            
            if (fipc.output(FIPConfig::OutputField::RESV))
                logOutput_.fipResv(inplace, "FIPNUM"); 
        }
        
        if (fipc.output(FIPConfig::OutputField::FIP)) {
            for (const auto& reg : this->regions_) {
                if (reg.first != "FIPNUM") {
                    std::ostringstream ss;
                    ss << "BAL" << reg.first.substr(3);
                    logOutput_.timeStamp(ss.str(), elapsed, reportStepNum, currentDate);
                    logOutput_.fip(inplace, this->initialInplace(), reg.first);
                    
                    if (fipc.output(FIPConfig::OutputField::RESV))
                        logOutput_.fipResv(inplace, reg.first); 
                }
            }
        }
    }
}


template<class FluidSystem,class Scalar>
void EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
addRftDataToWells(data::Wells& wellDatas, std::size_t reportStepNum)
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
            std::size_t count = 0;
            for (const auto& connection: well.getConnections()) {
                const std::size_t i = std::size_t(connection.getI());
                const std::size_t j = std::size_t(connection.getJ());
                const std::size_t k = std::size_t(connection.getK());

                const std::size_t index = eclState_.gridDims().getGlobalIndex(i, j, k);
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
        DataEntry{"FLRGASI+", UnitSystem::measure::rate,                                  flores_[FaceDir::ToIntersectionIndex(Dir::XPlus)][gasCompIdx]},
        DataEntry{"FLRGASJ+", UnitSystem::measure::rate,                                  flores_[FaceDir::ToIntersectionIndex(Dir::YPlus)][gasCompIdx]},
        DataEntry{"FLRGASK+", UnitSystem::measure::rate,                                  flores_[FaceDir::ToIntersectionIndex(Dir::ZPlus)][gasCompIdx]},
        DataEntry{"FLROILI+", UnitSystem::measure::rate,                                  flores_[FaceDir::ToIntersectionIndex(Dir::XPlus)][oilCompIdx]},
        DataEntry{"FLROILJ+", UnitSystem::measure::rate,                                  flores_[FaceDir::ToIntersectionIndex(Dir::YPlus)][oilCompIdx]},
        DataEntry{"FLROILK+", UnitSystem::measure::rate,                                  flores_[FaceDir::ToIntersectionIndex(Dir::ZPlus)][oilCompIdx]},
        DataEntry{"FLRWATI+", UnitSystem::measure::rate,                                  flores_[FaceDir::ToIntersectionIndex(Dir::XPlus)][waterCompIdx]},
        DataEntry{"FLRWATJ+", UnitSystem::measure::rate,                                  flores_[FaceDir::ToIntersectionIndex(Dir::YPlus)][waterCompIdx]},
        DataEntry{"FLRWATK+", UnitSystem::measure::rate,                                  flores_[FaceDir::ToIntersectionIndex(Dir::ZPlus)][waterCompIdx]},
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
        DataEntry{"RPORV",    UnitSystem::measure::volume,                                rPorV_},
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

    // Separate these as flows*_ may be defined due to BFLOW[I|J|K] even without FLOWS in RPTRST
    const auto flowsSolutionArrays = std::array {
        DataEntry{"FLOGASI+", UnitSystem::measure::gas_surface_rate,                      flows_[FaceDir::ToIntersectionIndex(Dir::XPlus)][gasCompIdx]},
        DataEntry{"FLOGASJ+", UnitSystem::measure::gas_surface_rate,                      flows_[FaceDir::ToIntersectionIndex(Dir::YPlus)][gasCompIdx]},
        DataEntry{"FLOGASK+", UnitSystem::measure::gas_surface_rate,                      flows_[FaceDir::ToIntersectionIndex(Dir::ZPlus)][gasCompIdx]},
        DataEntry{"FLOOILI+", UnitSystem::measure::liquid_surface_rate,                   flows_[FaceDir::ToIntersectionIndex(Dir::XPlus)][oilCompIdx]},
        DataEntry{"FLOOILJ+", UnitSystem::measure::liquid_surface_rate,                   flows_[FaceDir::ToIntersectionIndex(Dir::YPlus)][oilCompIdx]},
        DataEntry{"FLOOILK+", UnitSystem::measure::liquid_surface_rate,                   flows_[FaceDir::ToIntersectionIndex(Dir::ZPlus)][oilCompIdx]},
        DataEntry{"FLOWATI+", UnitSystem::measure::liquid_surface_rate,                   flows_[FaceDir::ToIntersectionIndex(Dir::XPlus)][waterCompIdx]},
        DataEntry{"FLOWATJ+", UnitSystem::measure::liquid_surface_rate,                   flows_[FaceDir::ToIntersectionIndex(Dir::YPlus)][waterCompIdx]},
        DataEntry{"FLOWATK+", UnitSystem::measure::liquid_surface_rate,                   flows_[FaceDir::ToIntersectionIndex(Dir::ZPlus)][waterCompIdx]},
        DataEntry{"FLOGASI-", UnitSystem::measure::gas_surface_rate,                      flows_[FaceDir::ToIntersectionIndex(Dir::XMinus)][gasCompIdx]},
        DataEntry{"FLOGASJ-", UnitSystem::measure::gas_surface_rate,                      flows_[FaceDir::ToIntersectionIndex(Dir::YMinus)][gasCompIdx]},
        DataEntry{"FLOGASK-", UnitSystem::measure::gas_surface_rate,                      flows_[FaceDir::ToIntersectionIndex(Dir::ZMinus)][gasCompIdx]},
        DataEntry{"FLOOILI-", UnitSystem::measure::liquid_surface_rate,                   flows_[FaceDir::ToIntersectionIndex(Dir::XMinus)][oilCompIdx]},
        DataEntry{"FLOOILJ-", UnitSystem::measure::liquid_surface_rate,                   flows_[FaceDir::ToIntersectionIndex(Dir::YMinus)][oilCompIdx]},
        DataEntry{"FLOOILK-", UnitSystem::measure::liquid_surface_rate,                   flows_[FaceDir::ToIntersectionIndex(Dir::ZMinus)][oilCompIdx]},
        DataEntry{"FLOWATI-", UnitSystem::measure::liquid_surface_rate,                   flows_[FaceDir::ToIntersectionIndex(Dir::XMinus)][waterCompIdx]},
        DataEntry{"FLOWATJ-", UnitSystem::measure::liquid_surface_rate,                   flows_[FaceDir::ToIntersectionIndex(Dir::YMinus)][waterCompIdx]},
        DataEntry{"FLOWATK-", UnitSystem::measure::liquid_surface_rate,                   flows_[FaceDir::ToIntersectionIndex(Dir::ZMinus)][waterCompIdx]},
        DataEntry{"FLRGASI-", UnitSystem::measure::rate,                                  flores_[FaceDir::ToIntersectionIndex(Dir::XMinus)][gasCompIdx]},
        DataEntry{"FLRGASJ-", UnitSystem::measure::rate,                                  flores_[FaceDir::ToIntersectionIndex(Dir::YMinus)][gasCompIdx]},
        DataEntry{"FLRGASK-", UnitSystem::measure::rate,                                  flores_[FaceDir::ToIntersectionIndex(Dir::ZMinus)][gasCompIdx]},
        DataEntry{"FLROILI-", UnitSystem::measure::rate,                                  flores_[FaceDir::ToIntersectionIndex(Dir::XMinus)][oilCompIdx]},
        DataEntry{"FLROILJ-", UnitSystem::measure::rate,                                  flores_[FaceDir::ToIntersectionIndex(Dir::YMinus)][oilCompIdx]},
        DataEntry{"FLROILK-", UnitSystem::measure::rate,                                  flores_[FaceDir::ToIntersectionIndex(Dir::ZMinus)][oilCompIdx]},
        DataEntry{"FLRWATI-", UnitSystem::measure::rate,                                  flores_[FaceDir::ToIntersectionIndex(Dir::XMinus)][waterCompIdx]},
        DataEntry{"FLRWATJ-", UnitSystem::measure::rate,                                  flores_[FaceDir::ToIntersectionIndex(Dir::YMinus)][waterCompIdx]},
        DataEntry{"FLRWATK-", UnitSystem::measure::rate,                                  flores_[FaceDir::ToIntersectionIndex(Dir::ZMinus)][waterCompIdx]},
    };

    const auto extendedSolutionArrays = std::array {
        DataEntry{"BIOFILM",  UnitSystem::measure::identity,           cBiofilm_},
        DataEntry{"CALCITE",  UnitSystem::measure::identity,           cCalcite_},
        DataEntry{"DELSTRXX", UnitSystem::measure::pressure,           delstressXX_},
        DataEntry{"DELSTRYY", UnitSystem::measure::pressure,           delstressYY_},
        DataEntry{"DELSTRZZ", UnitSystem::measure::pressure,           delstressZZ_},
        DataEntry{"DELSTRXY", UnitSystem::measure::pressure,           delstressXY_},
        DataEntry{"DELSTRXZ", UnitSystem::measure::pressure,           delstressXZ_},
        DataEntry{"DELSTRYZ", UnitSystem::measure::pressure,           delstressYZ_},
        DataEntry{"DISPX",    UnitSystem::measure::length,             dispX_},
        DataEntry{"DISPY",    UnitSystem::measure::length,             dispY_},
        DataEntry{"DISPZ",    UnitSystem::measure::length,             dispZ_},
        DataEntry{"DRSDTCON", UnitSystem::measure::gas_oil_ratio_rate, drsdtcon_},
        DataEntry{"KRNSW_GO", UnitSystem::measure::identity,           krnSwMdcGo_},
        DataEntry{"KRNSW_OW", UnitSystem::measure::identity,           krnSwMdcOw_},
        DataEntry{"MECHPOTF", UnitSystem::measure::pressure,           mechPotentialForce_},
        DataEntry{"MICROBES", UnitSystem::measure::density,            cMicrobes_},
        DataEntry{"OXYGEN",   UnitSystem::measure::density,            cOxygen_},
        DataEntry{"PCSWM_GO", UnitSystem::measure::identity,           pcSwMdcGo_},
        DataEntry{"PCSWM_OW", UnitSystem::measure::identity,           pcSwMdcOw_},
        DataEntry{"PERMFACT", UnitSystem::measure::identity,           permFact_},
        DataEntry{"PORV_RC",  UnitSystem::measure::identity,           rockCompPorvMultiplier_},
        DataEntry{"PRESPOTF", UnitSystem::measure::pressure,           mechPotentialPressForce_},
        DataEntry{"PRES_OVB", UnitSystem::measure::pressure,           overburdenPressure_},
        DataEntry{"RSW",      UnitSystem::measure::gas_oil_ratio,      rsw_},
        DataEntry{"RSWSOL",   UnitSystem::measure::gas_oil_ratio,      rswSol_},
        DataEntry{"RVW",      UnitSystem::measure::oil_gas_ratio,      rvw_},
        DataEntry{"SALTP",    UnitSystem::measure::identity,           pSalt_},
        DataEntry{"SS_X",     UnitSystem::measure::identity,           extboX_},
        DataEntry{"SS_Y",     UnitSystem::measure::identity,           extboY_},
        DataEntry{"SS_Z",     UnitSystem::measure::identity,           extboZ_},
        DataEntry{"STD_CO2",  UnitSystem::measure::identity,           mFracCo2_},
        DataEntry{"STD_GAS",  UnitSystem::measure::identity,           mFracGas_},
        DataEntry{"STD_OIL",  UnitSystem::measure::identity,           mFracOil_},
        DataEntry{"STRAINXX", UnitSystem::measure::identity,           strainXX_},
        DataEntry{"STRAINYY", UnitSystem::measure::identity,           strainYY_},
        DataEntry{"STRAINZZ", UnitSystem::measure::identity,           strainZZ_},
        DataEntry{"STRAINXY", UnitSystem::measure::identity,           strainXY_},
        DataEntry{"STRAINXZ", UnitSystem::measure::identity,           strainXZ_},
        DataEntry{"STRAINYZ", UnitSystem::measure::identity,           strainYZ_},
        DataEntry{"STRESSXX", UnitSystem::measure::length,             stressXX_},
        DataEntry{"STRESSYY", UnitSystem::measure::length,             stressYY_},
        DataEntry{"STRESSZZ", UnitSystem::measure::length,             stressZZ_},
        DataEntry{"STRESSXY", UnitSystem::measure::length,             stressXY_},
        DataEntry{"STRESSXZ", UnitSystem::measure::length,             stressXZ_},
        DataEntry{"STRESSYZ", UnitSystem::measure::length,             stressYZ_},
        DataEntry{"TEMPPOTF", UnitSystem::measure::pressure,           mechPotentialTempForce_},
        DataEntry{"TMULT_RC", UnitSystem::measure::identity,           rockCompTransMultiplier_},
        DataEntry{"UREA",     UnitSystem::measure::density,            cUrea_},
    };

    for (const auto& array : baseSolutionArrays) {
        doInsert(array, data::TargetType::RESTART_SOLUTION);
    }

    if (this->enableFlows_) {
        for (const auto& array : flowsSolutionArrays) {
            doInsert(array, data::TargetType::RESTART_SOLUTION);
        }
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

    if ((eclState_.runspec().co2Storage() || eclState_.runspec().h2Storage()) && !rsw_.empty()) {
        auto mfrac = std::vector<double>(this->rsw_.size(), 0.0);

        std::transform(this->rsw_.begin(), this->rsw_.end(),
                       this->eclState_.fieldProps().get_int("PVTNUM").begin(),
                       mfrac.begin(),
            [](const auto& rsw, const int pvtReg)
        {
            const auto xwg = FluidSystem::convertRswToXwG(rsw, pvtReg - 1);
            return FluidSystem::convertXwGToxwG(xwg, pvtReg - 1);
        });

        std::string moleFracName = eclState_.runspec().co2Storage() ? "XMFCO2" : "XMFH2";
        sol.insert(moleFracName,
                   UnitSystem::measure::identity,
                   std::move(mfrac),
                   data::TargetType::RESTART_OPM_EXTENDED);
    }

    if ((eclState_.runspec().co2Storage() || eclState_.runspec().h2Storage()) && !rvw_.empty()) {
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

    if (FluidSystem::phaseIsActive(waterPhaseIdx) &&
        ! this->residual_[waterPhaseIdx].empty())
    {
        sol.insert("RES_WAT", UnitSystem::measure::liquid_surface_volume,
                   std::move(this->residual_[waterPhaseIdx]),
                   data::TargetType::RESTART_OPM_EXTENDED);
    }
    if (FluidSystem::phaseIsActive(gasPhaseIdx) &&
        ! this->residual_[gasPhaseIdx].empty())
    {
        sol.insert("RES_GAS", UnitSystem::measure::gas_surface_volume,
                   std::move(this->residual_[gasPhaseIdx]),
                   data::TargetType::RESTART_OPM_EXTENDED);
    }
    if (FluidSystem::phaseIsActive(oilPhaseIdx) &&
        ! this->residual_[oilPhaseIdx].empty())
    {
        sol.insert("RES_OIL", UnitSystem::measure::liquid_surface_volume,
                   std::move(this->residual_[oilPhaseIdx]),
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
        saturation_[waterPhaseIdx][elemIdx] = sol.data<double>("SWAT")[globalDofIndex];
        so -= sol.data<double>("SWAT")[globalDofIndex];
    }
    if (!saturation_[gasPhaseIdx].empty() && sol.has("SGAS")) {
        saturation_[gasPhaseIdx][elemIdx] = sol.data<double>("SGAS")[globalDofIndex];
        so -= sol.data<double>("SGAS")[globalDofIndex];
    }

    if (!sSol_.empty()) {
        // keep the SSOL option for backward compatibility
        // should be removed after 10.2018 release
        if (sol.has("SSOL"))
            sSol_[elemIdx] = sol.data<double>("SSOL")[globalDofIndex];
        else if (sol.has("SSOLVENT"))
            sSol_[elemIdx] = sol.data<double>("SSOLVENT")[globalDofIndex];

        so -= sSol_[elemIdx];
    }

    if (!rswSol_.empty()) {
        if (sol.has("RSWSOL"))
            rswSol_[elemIdx] = sol.data<Scalar>("RSWSOL")[globalDofIndex];

    }

    assert(!saturation_[oilPhaseIdx].empty());
    saturation_[oilPhaseIdx][elemIdx] = so;

    auto assign = [elemIdx, globalDofIndex, &sol](const std::string& name,
                                                  ScalarBuffer& data)

    {
        if (!data.empty() && sol.has(name)) {
            data[elemIdx] = sol.data<double>(name)[globalDofIndex];
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
          std::size_t maxNumberOfRegions,
          const Parallel::Communication& comm)
{
        ScalarBuffer totals(maxNumberOfRegions, 0.0);

        if (property.empty())
            return totals;

        // the regionId contains the ghost cells
        // the property does not contain the ghostcells
        // This code assumes that that the ghostcells are
        // added after the interior cells
        // OwnerCellsFirst = True
        assert(regionId.size() >= property.size());
        for (std::size_t j = 0; j < property.size(); ++j) {
            const int regionIdx = regionId[j] - 1;
            // the cell is not attributed to any region. ignore it!
            if (regionIdx < 0)
                continue;

            assert(regionIdx < static_cast<int>(maxNumberOfRegions));
            totals[regionIdx] += property[j];
        }

        for (std::size_t i = 0; i < maxNumberOfRegions; ++i)
            totals[i] = comm.sum(totals[i]);

        return totals;
    }

template<class FluidSystem, class Scalar>
void EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
doAllocBuffers(const unsigned bufferSize,
               const unsigned reportStepNum,
               const bool     substep,
               const bool     log,
               const bool     isRestart,
               const bool     vapparsActive,
               const bool     enableHysteresis,
               const unsigned numTracers,
               const unsigned numOutputNnc)
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

    const auto needAvgPress = !substep         ||
        !this->RPRNodes_.empty()               ||
        this->summaryConfig_.hasKeyword("FPR") ||
        this->summaryConfig_.hasKeyword("FPRP");

    const auto needPoreVolume = needAvgPress    ||
        this->summaryConfig_.hasKeyword("FHPV") ||
        this->summaryConfig_.match("RHPV*");

    if (needPoreVolume) {
        this->fip_[Inplace::Phase::PoreVolume].resize(bufferSize, 0.0);
        this->dynamicPoreVolume_.resize(bufferSize, 0.0);
        this->hydrocarbonPoreVolume_.resize(bufferSize, 0.0);
    }
    else {
        this->dynamicPoreVolume_.clear();
        this->hydrocarbonPoreVolume_.clear();
    }

    if (needAvgPress) {
        this->pressureTimesPoreVolume_.resize(bufferSize, 0.0);
        this->pressureTimesHydrocarbonVolume_.resize(bufferSize, 0.0);
    }
    else {
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
                const std::size_t i = std::size_t(connection.getI());
                const std::size_t j = std::size_t(connection.getJ());
                const std::size_t k = std::size_t(connection.getK());
                const std::size_t index = eclState_.gridDims().getGlobalIndex(i, j, k);

                if (FluidSystem::phaseIsActive(oilPhaseIdx))
                    oilConnectionPressures_.emplace(std::make_pair(index, 0.0));

                if (FluidSystem::phaseIsActive(waterPhaseIdx))
                    waterConnectionSaturations_.emplace(std::make_pair(index, 0.0));

                if (FluidSystem::phaseIsActive(gasPhaseIdx))
                    gasConnectionSaturations_.emplace(std::make_pair(index, 0.0));
            }
        }
    }

    // Flows may need to be allocated even when there is no restart due to BFLOW* summary keywords
    if (blockFlows_ ) {
        const std::array<int, 3> phaseIdxs = { gasPhaseIdx, oilPhaseIdx, waterPhaseIdx };
        const std::array<int, 3> compIdxs = { gasCompIdx, oilCompIdx, waterCompIdx };

        for (unsigned ii = 0; ii < phaseIdxs.size(); ++ii) {
            if (FluidSystem::phaseIsActive(phaseIdxs[ii])) {
                flows_[FaceDir::ToIntersectionIndex(Dir::XPlus)][compIdxs[ii]].resize(bufferSize, 0.0);
                flows_[FaceDir::ToIntersectionIndex(Dir::YPlus)][compIdxs[ii]].resize(bufferSize, 0.0);
                flows_[FaceDir::ToIntersectionIndex(Dir::ZPlus)][compIdxs[ii]].resize(bufferSize, 0.0);
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

    if (enableMech_ && eclState_.runspec().mech()) {
        this->mechPotentialForce_.resize(bufferSize,0.0);
        rstKeywords["MECHPOTF"] = 0;
        this->mechPotentialTempForce_.resize(bufferSize,0.0);
        rstKeywords["TEMPPOTF"] = 0;
        this->mechPotentialPressForce_.resize(bufferSize,0.0);
        rstKeywords["PRESPOTF"] = 0;

        this->dispX_.resize(bufferSize,0.0);
        rstKeywords["DISPX"] = 0;
        this->dispY_.resize(bufferSize,0.0);
        rstKeywords["DISPY"] = 0;
        this->dispZ_.resize(bufferSize,0.0);
        rstKeywords["DISPZ"] = 0;
        this->stressXX_.resize(bufferSize,0.0);
        rstKeywords["STRESSXX"] = 0;
        this->stressYY_.resize(bufferSize,0.0);
        rstKeywords["STRESSYY"] = 0;
        this->stressZZ_.resize(bufferSize,0.0);
        rstKeywords["STRESSZZ"] = 0;
        this->stressXY_.resize(bufferSize,0.0);
        rstKeywords["STRESSXY"] = 0;
        this->stressXZ_.resize(bufferSize,0.0);
        rstKeywords["STRESSXZ"] = 0;
        this->stressXY_.resize(bufferSize,0.0);
        rstKeywords["STRESSXY"] = 0;
        this->stressYZ_.resize(bufferSize,0.0);
        rstKeywords["STRESSYZ"] = 0;

        this->strainXX_.resize(bufferSize,0.0);
        rstKeywords["STRAINXX"] = 0;
        this->strainYY_.resize(bufferSize,0.0);
        rstKeywords["STRAINYY"] = 0;
        this->strainZZ_.resize(bufferSize,0.0);
        rstKeywords["STRAINZZ"] = 0;
        this->strainXY_.resize(bufferSize,0.0);
        rstKeywords["STRAINXY"] = 0;
        this->strainXZ_.resize(bufferSize,0.0);
        rstKeywords["STRAINXZ"] = 0;
        this->strainXY_.resize(bufferSize,0.0);
        rstKeywords["STRAINXY"] = 0;
        this->strainYZ_.resize(bufferSize,0.0);
        rstKeywords["STRAINYZ"] = 0;

        this->delstressXX_.resize(bufferSize,0.0);
        rstKeywords["DELSTRXX"] = 0;
        this->delstressYY_.resize(bufferSize,0.0);
        rstKeywords["DELSTRYY"] = 0;
        this->delstressZZ_.resize(bufferSize,0.0);
        rstKeywords["DELSTRZZ"] = 0;
        this->delstressXY_.resize(bufferSize,0.0);
        rstKeywords["DELSTRXY"] = 0;
        this->delstressXZ_.resize(bufferSize,0.0);
        rstKeywords["DELSTRXZ"] = 0;
        this->delstressXY_.resize(bufferSize,0.0);
        rstKeywords["DELSTRXY"] = 0;
        this->delstressYZ_.resize(bufferSize,0.0);
        rstKeywords["DELSTRYZ"] = 0;
    }

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
        if (eclState_.getSimulationConfig().hasDISGASW()) {
            rswSol_.resize(bufferSize, 0.0);
        }
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
    if (rstKeywords["RPORV"] > 0) {
        rstKeywords["RPORV"] = 0;
        rPorV_.resize(bufferSize, 0.0);
    }

    enableFlows_ = false;
    enableFlowsn_ = false;
    const bool rstFlows = (rstKeywords["FLOWS"] > 0);
    if (rstFlows) {
        rstKeywords["FLOWS"] = 0;
        enableFlows_ = true;

        const std::array<int, 3> phaseIdxs = { gasPhaseIdx, oilPhaseIdx, waterPhaseIdx };
        const std::array<int, 3> compIdxs = { gasCompIdx, oilCompIdx, waterCompIdx };
        const auto rstName = std::array { "FLOGASN+", "FLOOILN+", "FLOWATN+" };

        for (unsigned ii = 0; ii < phaseIdxs.size(); ++ii) {
            if (FluidSystem::phaseIsActive(phaseIdxs[ii])) {
                if (!blockFlows_) { // Already allocated if summary vectors requested
                    flows_[FaceDir::ToIntersectionIndex(Dir::XPlus)][compIdxs[ii]].resize(bufferSize, 0.0);
                    flows_[FaceDir::ToIntersectionIndex(Dir::YPlus)][compIdxs[ii]].resize(bufferSize, 0.0);
                    flows_[FaceDir::ToIntersectionIndex(Dir::ZPlus)][compIdxs[ii]].resize(bufferSize, 0.0);
                }

                if (rstKeywords["FLOWS-"] > 0) {
                    flows_[FaceDir::ToIntersectionIndex(Dir::XMinus)][compIdxs[ii]].resize(bufferSize, 0.0);
                    flows_[FaceDir::ToIntersectionIndex(Dir::YMinus)][compIdxs[ii]].resize(bufferSize, 0.0);
                    flows_[FaceDir::ToIntersectionIndex(Dir::ZMinus)][compIdxs[ii]].resize(bufferSize, 0.0);  
                }

                if (numOutputNnc > 0) {
                    enableFlowsn_ = true;

                    flowsn_[compIdxs[ii]].first = rstName[ii];
                    flowsn_[compIdxs[ii]].second.first.resize(numOutputNnc, -1);
                    flowsn_[compIdxs[ii]].second.second.resize(numOutputNnc, 0.0);
                }
            }
        }
        if (rstKeywords["FLOWS-"] > 0) {
            rstKeywords["FLOWS-"] = 0;
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
                flores_[FaceDir::ToIntersectionIndex(Dir::XPlus)][compIdxs[ii]].resize(bufferSize, 0.0);
                flores_[FaceDir::ToIntersectionIndex(Dir::YPlus)][compIdxs[ii]].resize(bufferSize, 0.0);
                flores_[FaceDir::ToIntersectionIndex(Dir::ZPlus)][compIdxs[ii]].resize(bufferSize, 0.0);

                if (rstKeywords["FLORES-"] > 0) {
                    flores_[FaceDir::ToIntersectionIndex(Dir::XMinus)][compIdxs[ii]].resize(bufferSize, 0.0);
                    flores_[FaceDir::ToIntersectionIndex(Dir::YMinus)][compIdxs[ii]].resize(bufferSize, 0.0);
                    flores_[FaceDir::ToIntersectionIndex(Dir::ZMinus)][compIdxs[ii]].resize(bufferSize, 0.0);  
                }

                if (numOutputNnc > 0) {
                    enableFloresn_ = true;

                    floresn_[compIdxs[ii]].first = rstName[ii];
                    floresn_[compIdxs[ii]].second.first.resize(numOutputNnc, -1);
                    floresn_[compIdxs[ii]].second.second.resize(numOutputNnc, 0.0);
                }
            }
        }
        if (rstKeywords["FLORES-"] > 0) {
            rstKeywords["FLORES-"] = 0;
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

    if (rstKeywords["RESIDUAL"] > 0) {
        rstKeywords["RESIDUAL"] = 0;
        for (int phaseIdx = 0; phaseIdx <  numPhases; ++phaseIdx)
        {
            if (FluidSystem::phaseIsActive(phaseIdx)) {
                this->residual_[phaseIdx].resize(bufferSize, 0.0);
            }
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
                logstring.append("' is unhandled for output to restart file.");
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

template<class FluidSystem,class Scalar>
bool EclGenericOutputBlackoilModule<FluidSystem,Scalar>::
isOutputCreationDirective_(const std::string& keyword)
{
    return (keyword == "BASIC") || (keyword == "FREQ")
        || (keyword == "RESTART")                        // From RPTSCHED
        || (keyword == "SAVE")  || (keyword == "SFREQ"); // Not really supported
}

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

    logOutput_.error(std::get<0>(globalFailedCellsPbub),
                     std::get<0>(globalFailedCellsPdew));
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
    // the initialInplace_ member. This has a problem:
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

        if (this->summaryConfig_.hasKeyword("FHPV")) {
            miscSummaryData["FHPV"] = inplace.get(Inplace::Phase::HydroCarbonPV);
        }

        if (this->summaryConfig_.hasKeyword("FOE") && this->initialInplace_) {
            miscSummaryData["FOE"] = (this->initialInplace_.value().get(Inplace::Phase::OIL) - inplace.get(Inplace::Phase::OIL))
                / this->initialInplace_.value().get(Inplace::Phase::OIL);
        }

        if (this->summaryConfig_.hasKeyword("FPR")) {
            miscSummaryData["FPR"] =
                detail::pressureAverage(inplace.get(Inplace::Phase::PressureHydroCarbonPV),
                                        inplace.get(Inplace::Phase::HydroCarbonPV),
                                        inplace.get(Inplace::Phase::PressurePV),
                                        inplace.get(Inplace::Phase::DynamicPoreVolume),
                                        true);
        }

        if (this->summaryConfig_.hasKeyword("FPRP")) {
            miscSummaryData["FPRP"] =
                detail::pressureAverage(inplace.get(Inplace::Phase::PressureHydroCarbonPV),
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
                detail::pressureAverage(get_vector(node, Inplace::Phase::PressureHydroCarbonPV),
                                        get_vector(node, Inplace::Phase::HydroCarbonPV),
                                        get_vector(node, Inplace::Phase::PressurePV),
                                        get_vector(node, Inplace::Phase::DynamicPoreVolume),
                                        true);
        }

        for (const auto& node : this->RPRPNodes_) {
            regionData[node.keyword()] =
                detail::pressureAverage(get_vector(node, Inplace::Phase::PressureHydroCarbonPV),
                                        get_vector(node, Inplace::Phase::HydroCarbonPV),
                                        get_vector(node, Inplace::Phase::PressurePV),
                                        get_vector(node, Inplace::Phase::DynamicPoreVolume),
                                        false);
        }

        for (const auto& node : this->summaryConfig_.keywords("RHPV*")) {
            regionData[node.keyword()] =
                get_vector(node, Inplace::Phase::HydroCarbonPV);
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
