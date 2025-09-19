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
#include <opm/simulators/flow/GenericOutputBlackoilModule.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/utility/Visitor.hpp>

#include <opm/grid/common/CommunicationUtils.hpp>

#include <opm/material/fluidmatrixinteractions/EclHysteresisConfig.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>
#include <opm/material/fluidsystems/GenericOilGasWaterFluidSystem.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Runspec.hpp>
#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>

#include <opm/input/eclipse/Schedule/RFTConfig.hpp>
#include <opm/input/eclipse/Schedule/RPTConfig.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/input/eclipse/Schedule/Well/NameOrder.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/models/utils/parametersystem.hpp>

#include <opm/output/data/Solution.hpp>

#include <opm/simulators/utils/PressureAverage.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include <fmt/format.h>

namespace {

    std::size_t numCells(const Opm::EclipseState& eclState)
    {
        return eclState.fieldProps().get_int("FIPNUM").size();
    }

    std::vector<Opm::InterRegFlowMap::SingleRegion>
    defineInterRegionFlowArrays(const Opm::EclipseState&  eclState,
                                const Opm::SummaryConfig& summaryConfig)
    {
        using SRegion = Opm::InterRegFlowMap::SingleRegion;
        auto regions = std::vector<SRegion>{};
        const auto fip_regions = summaryConfig.fip_regions_interreg_flow();
        std::transform(fip_regions.begin(), fip_regions.end(),
                       std::back_inserter(regions),
                       [&fprops = eclState.fieldProps()](const auto& arrayName)
                       { return SRegion{arrayName, std::cref(fprops.get_int(arrayName))}; });

        return regions;
    }

    //! \brief Functor for reordering phase indices for name resolution.
    //! \details This is the identity mapping for BlackoilFluidSystem.
    template<class T>
    struct PhaseReorder
    {
        int operator()(const int phase)
        {
            return phase;
        }
    };

    //! \brief Functor for reordering phase indices for name resolution.
    //! \details Re-order the GenericOilGasWaterFluidSystem phase indices to match
    //!          BlackoilFluidSystem order.
    template<class Scalar, int NumComp, bool enableWater>
    struct PhaseReorder<Opm::GenericOilGasWaterFluidSystem<Scalar,NumComp,enableWater>>
    {
        int operator()(const int phase)
        {
            switch (phase) {
            case 0: return 1;
            case 1: return 2;
            default: return 0;
            }
        }
    };
}

namespace Opm {

template<class FluidSystem>
GenericOutputBlackoilModule<FluidSystem>::
GenericOutputBlackoilModule(const EclipseState& eclState,
                            const Schedule& schedule,
                            const SummaryConfig& summaryConfig,
                            const SummaryState& summaryState,
                            const std::string& moduleVersion,
                            RSTConv::LocalToGlobalCellFunc globalCell,
                            const Parallel::Communication& comm,
                            bool enableEnergy,
                            bool enableTemperature,
                            bool enableMech,
                            bool enableSolvent,
                            bool enablePolymer,
                            bool enableFoam,
                            bool enableBrine,
                            bool enableSaltPrecipitation,
                            bool enableExtbo,
                            bool enableBioeffects)
    : eclState_(eclState)
    , schedule_(schedule)
    , summaryState_(summaryState)
    , summaryConfig_(summaryConfig)
    , interRegionFlows_(numCells(eclState),
                        defineInterRegionFlowArrays(eclState, summaryConfig),
                        declaredMaxRegionID(eclState.runspec()))
    , logOutput_(eclState, schedule, summaryState, moduleVersion)
    , enableEnergy_(enableEnergy)
    , enableTemperature_(enableTemperature)
    , enableMech_(enableMech)
    , enableSolvent_(enableSolvent)
    , enablePolymer_(enablePolymer)
    , enableFoam_(enableFoam)
    , enableBrine_(enableBrine)
    , enableSaltPrecipitation_(enableSaltPrecipitation)
    , enableExtbo_(enableExtbo)
    , enableBioeffects_(enableBioeffects)
    , flowsC_(schedule, summaryConfig)
    , rftC_(eclState_, schedule_,
            [this](const std::string& wname) { return this->isOwnedByCurrentRank(wname); },
            [this](const std::string& wname) { return this->isOnCurrentRank(wname); })
    , rst_conv_(std::move(globalCell), comm)
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
        std::string key_pattern = "R" + Inplace::EclString(phase) + "*";
        this->regionNodes_[phase] = summaryConfig_.keywords(key_pattern);
    }

    forceDisableFipOutput_ =
        Parameters::Get<Parameters::ForceDisableFluidInPlaceOutput>();
    forceDisableFipresvOutput_ =
        Parameters::Get<Parameters::ForceDisableResvFluidInPlaceOutput>();
}

template<class FluidSystem>
GenericOutputBlackoilModule<FluidSystem>::
~GenericOutputBlackoilModule() = default;

template<class FluidSystem>
void GenericOutputBlackoilModule<FluidSystem>::
registerParameters()
{
    Parameters::Register<Parameters::ForceDisableFluidInPlaceOutput>
        ("Do not print fluid-in-place values after each report step "
         "even if requested by the deck.");
    Parameters::Register<Parameters::ForceDisableResvFluidInPlaceOutput>
        ("Do not print reservoir volumes values after each report step "
         "even if requested by the deck.");
}

template<class FluidSystem>
void GenericOutputBlackoilModule<FluidSystem>::
outputTimeStamp(const std::string& lbl,
                const double elapsed,
                const int rstep,
                const boost::posix_time::ptime currentDate)
{
    logOutput_.timeStamp(lbl, elapsed, rstep, currentDate);
}

template<class FluidSystem>
void GenericOutputBlackoilModule<FluidSystem>::
prepareDensityAccumulation()
{
    if (this->regionAvgDensity_.has_value()) {
        this->regionAvgDensity_->prepareAccumulation();
    }
}

template<class FluidSystem>
void GenericOutputBlackoilModule<FluidSystem>::
accumulateDensityParallel()
{
    if (this->regionAvgDensity_.has_value()) {
        this->regionAvgDensity_->accumulateParallel();
    }
}

template<class FluidSystem>
void GenericOutputBlackoilModule<FluidSystem>::
outputCumLog(std::size_t reportStepNum,
             const bool connData)
{
    this->logOutput_.cumulative(reportStepNum, connData);
}

template<class FluidSystem>
void GenericOutputBlackoilModule<FluidSystem>::
outputProdLog(std::size_t reportStepNum,
              const bool connData)
{
    this->logOutput_.production(reportStepNum,
                                connData
                                    ? this->extraBlockData_
                                    : decltype(this->extraBlockData_){});
}

template<class FluidSystem>
void GenericOutputBlackoilModule<FluidSystem>::
outputInjLog(std::size_t reportStepNum,
             const bool connData)
{
    this->logOutput_.injection(reportStepNum,
                               connData
                                   ? this->extraBlockData_
                                   : decltype(this->extraBlockData_){});
}

template<class FluidSystem>
void GenericOutputBlackoilModule<FluidSystem>::
outputMSWLog(std::size_t reportStepNum)
{
    this->logOutput_.msw(reportStepNum);
}

template<class FluidSystem>
void GenericOutputBlackoilModule<FluidSystem>::
calc_initial_inplace(const Parallel::Communication& comm)
{
    if (!this->initialInplace_.has_value()) {
        this->initialInplace_ = this->accumulateRegionSums(comm);
    }
}

template<class FluidSystem>
Inplace GenericOutputBlackoilModule<FluidSystem>::
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

template<class FluidSystem>
void GenericOutputBlackoilModule<FluidSystem>::
outputWellspecReport(const std::vector<std::string>& changedWells,
                     const std::size_t               reportStepNum,
                     const double                    elapsed,
                     boost::posix_time::ptime        currentDate) const
{
    this->logOutput_.timeStamp("WELSPECS", elapsed,
                               static_cast<int>(reportStepNum),
                               currentDate);
    this->logOutput_.wellSpecification(changedWells, reportStepNum);
}

template<class FluidSystem>
void GenericOutputBlackoilModule<FluidSystem>::
assignToSolution(data::Solution& sol)
{
    using DataEntry =
        std::tuple<std::string, UnitSystem::measure, std::vector<Scalar>&>;

    auto doInsert = [&sol](DataEntry&       entry,
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

    // if index not specified, we treat it as valid (>= 0)
    auto addEntry = [](std::vector<DataEntry>& container,
                       const std::string& name,
                       UnitSystem::measure measure,
                       auto& flowArray, int index = 1)
    {
        if (index >= 0) {  // Only add if index is valid
            container.emplace_back(name, measure, flowArray);
        }
    };

    std::vector<DataEntry> baseSolutionVector;
    addEntry(baseSolutionVector, "1OVERBG",  UnitSystem::measure::gas_inverse_formation_volume_factor,   invB_[gasPhaseIdx],                                              gasPhaseIdx);
    addEntry(baseSolutionVector, "1OVERBO",  UnitSystem::measure::oil_inverse_formation_volume_factor,   invB_[oilPhaseIdx],                                              oilPhaseIdx);

    // avoid output with generic fluid system and disabled water phase
    if constexpr (numPhases > 2) {
        addEntry(baseSolutionVector,
                                 "1OVERBW",  UnitSystem::measure::water_inverse_formation_volume_factor, invB_[waterPhaseIdx],                                            waterPhaseIdx);
    }
    addEntry(baseSolutionVector, "FOAM",     UnitSystem::measure::identity,                              cFoam_);
    addEntry(baseSolutionVector, "GASKR",    UnitSystem::measure::identity,                              relativePermeability_[gasPhaseIdx],                              gasPhaseIdx);
    addEntry(baseSolutionVector, "GAS_DEN",  UnitSystem::measure::density,                               density_[gasPhaseIdx],                                           gasPhaseIdx);
    addEntry(baseSolutionVector, "GAS_VISC", UnitSystem::measure::viscosity,                             viscosity_[gasPhaseIdx],                                         gasPhaseIdx);
    addEntry(baseSolutionVector, "OILKR",    UnitSystem::measure::identity,                              relativePermeability_[oilPhaseIdx],                              oilPhaseIdx);
    addEntry(baseSolutionVector, "OIL_DEN",  UnitSystem::measure::density,                               density_[oilPhaseIdx],                                           oilPhaseIdx);
    addEntry(baseSolutionVector, "OIL_VISC", UnitSystem::measure::viscosity,                             viscosity_[oilPhaseIdx],                                         oilPhaseIdx);
    addEntry(baseSolutionVector, "PBUB",     UnitSystem::measure::pressure,                              bubblePointPressure_);
    addEntry(baseSolutionVector, "PCGW",     UnitSystem::measure::pressure,                              pcgw_);
    addEntry(baseSolutionVector, "PCOG",     UnitSystem::measure::pressure,                              pcog_);
    addEntry(baseSolutionVector, "PCOW",     UnitSystem::measure::pressure,                              pcow_);
    addEntry(baseSolutionVector, "PDEW",     UnitSystem::measure::pressure,                              dewPointPressure_);
    addEntry(baseSolutionVector, "POLYMER",  UnitSystem::measure::identity,                              cPolymer_);
    addEntry(baseSolutionVector, "PPCW",     UnitSystem::measure::pressure,                              ppcw_);
    addEntry(baseSolutionVector, "PRESROCC", UnitSystem::measure::pressure,                              minimumOilPressure_);
    addEntry(baseSolutionVector, "PRESSURE", UnitSystem::measure::pressure,                              fluidPressure_);
    addEntry(baseSolutionVector, "RPORV",    UnitSystem::measure::volume,                                rPorV_);
    addEntry(baseSolutionVector, "RS",       UnitSystem::measure::gas_oil_ratio,                         rs_);
    addEntry(baseSolutionVector, "RSSAT",    UnitSystem::measure::gas_oil_ratio,                         gasDissolutionFactor_);
    addEntry(baseSolutionVector, "RV",       UnitSystem::measure::oil_gas_ratio,                         rv_);
    addEntry(baseSolutionVector, "RVSAT",    UnitSystem::measure::oil_gas_ratio,                         oilVaporizationFactor_);
    addEntry(baseSolutionVector, "SALT",     UnitSystem::measure::concentration,                         cSalt_);
    addEntry(baseSolutionVector, "SGMAX",    UnitSystem::measure::identity,                              sgmax_);
    addEntry(baseSolutionVector, "SHMAX",    UnitSystem::measure::identity,                              shmax_);
    addEntry(baseSolutionVector, "SOMAX",    UnitSystem::measure::identity,                              soMax_);
    addEntry(baseSolutionVector, "SOMIN",    UnitSystem::measure::identity,                              somin_);
    addEntry(baseSolutionVector, "SSOLVENT", UnitSystem::measure::identity,                              sSol_);
    addEntry(baseSolutionVector, "SWHY1",    UnitSystem::measure::identity,                              swmin_);
    addEntry(baseSolutionVector, "SWMAX",    UnitSystem::measure::identity,                              swMax_);

    // avoid output with generic fluid system and disabled water phase
    if constexpr (numPhases > 2) {
        addEntry(baseSolutionVector, "WATKR",    UnitSystem::measure::identity,                          relativePermeability_[waterPhaseIdx],                            waterPhaseIdx);
        addEntry(baseSolutionVector, "WAT_DEN",  UnitSystem::measure::density,                           density_[waterPhaseIdx],                                         waterPhaseIdx);
        addEntry(baseSolutionVector, "WAT_VISC", UnitSystem::measure::viscosity,                         viscosity_[waterPhaseIdx],                                       waterPhaseIdx);
    }

    auto extendedSolutionArrays = std::array {
        DataEntry{"DRSDTCON", UnitSystem::measure::gas_oil_ratio_rate, drsdtcon_},
        DataEntry{"PERMFACT", UnitSystem::measure::identity,           permFact_},
        DataEntry{"PORV_RC",  UnitSystem::measure::identity,           rockCompPorvMultiplier_},
        DataEntry{"PRES_OVB", UnitSystem::measure::pressure,           overburdenPressure_},
        DataEntry{"RSW",      UnitSystem::measure::gas_oil_ratio,      rsw_},
        DataEntry{"RSWSAT",   UnitSystem::measure::gas_oil_ratio,      gasDissolutionFactorInWater_},
        DataEntry{"RSWSOL",   UnitSystem::measure::gas_oil_ratio,      rswSol_},
        DataEntry{"RVW",      UnitSystem::measure::oil_gas_ratio,      rvw_},
        DataEntry{"RVWSAT",   UnitSystem::measure::oil_gas_ratio,      waterVaporizationFactor_},
        DataEntry{"SALTP",    UnitSystem::measure::identity,           pSalt_},
        DataEntry{"TMULT_RC", UnitSystem::measure::identity,           rockCompTransMultiplier_},
    };

    for (auto& array : baseSolutionVector) {
        doInsert(array, data::TargetType::RESTART_SOLUTION);
    }

    this->flowsC_.outputRestart(sol);

    if (this->bioeffectsC_.allocated()) {
        // Biofilms for gas-water systems; MICP only for water systems
        this->bioeffectsC_.outputRestart(sol, !FluidSystem::phaseIsActive(gasPhaseIdx));
    }

    for (auto& array : extendedSolutionArrays) {
        doInsert(array, data::TargetType::RESTART_OPM_EXTENDED);
    }

    this->mech_.outputRestart(sol);
    this->extboC_.outputRestart(sol);

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
    this->fipC_.outputRestart(sol);

    // Tracers
    this->tracerC_.outputRestart(sol, eclState_.tracer());
}

template<class FluidSystem>
void GenericOutputBlackoilModule<FluidSystem>::
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
            rswSol_[elemIdx] = sol.data<double>("RSWSOL")[globalDofIndex];

    }
    if (!saturation_[oilPhaseIdx].empty())  {
        saturation_[oilPhaseIdx][elemIdx] = so;
    }

    auto assign = [elemIdx, globalDofIndex, &sol](const std::string& name,
                                                  ScalarBuffer& data)

    {
        if (!data.empty() && sol.has(name)) {
            data[elemIdx] = sol.data<double>(name)[globalDofIndex];
        }
    };

    const auto fields = std::array{
        std::pair{"FOAM",     &cFoam_},
        std::pair{"PERMFACT", &permFact_},
        std::pair{"POLYMER",  &cPolymer_},
        std::pair{"PPCW",     &ppcw_},
        std::pair{"PRESSURE", &fluidPressure_},
        std::pair{"RS",       &rs_},
        std::pair{"RSW",      &rsw_},
        std::pair{"RV",       &rv_},
        std::pair{"RVW",      &rvw_},
        std::pair{"SALT",     &cSalt_},
        std::pair{"SALTP",    &pSalt_},
        std::pair{"SGMAX",    &sgmax_},
        std::pair{"SHMAX",    &shmax_},
        std::pair{"SOMAX",    &soMax_},
        std::pair{"SOMIN",    &somin_},
        std::pair{"SWHY1",    &swmin_},
        std::pair{"SWMAX",    &swMax_},
        std::pair{"TEMP",     &temperature_},
    };

    std::for_each(fields.begin(), fields.end(),
                  [&assign](const auto& p)
                  { assign(p.first, *p.second); });

    if (this->bioeffectsC_.allocated()) {
        // Biofilms for gas-water systems; MICP only for water systems
        this->bioeffectsC_.readRestart(globalDofIndex, elemIdx, sol, !FluidSystem::phaseIsActive(gasPhaseIdx));
    }
}

template<class FluidSystem>
typename GenericOutputBlackoilModule<FluidSystem>::ScalarBuffer
GenericOutputBlackoilModule<FluidSystem>::
regionSum(const ScalarBuffer& property,
          const std::vector<int>& regionId,
          std::size_t maxNumberOfRegions,
          const Parallel::Communication& comm)
{
    ScalarBuffer totals(maxNumberOfRegions, 0.0);

    if (property.empty()) {
        return totals;
    }

    // the regionId contains the ghost cells
    // the property does not contain the ghostcells
    // This code assumes that that the ghostcells are
    // added after the interior cells
    // OwnerCellsFirst = True
    assert(regionId.size() >= property.size());
    for (std::size_t j = 0; j < property.size(); ++j) {
        const int regionIdx = regionId[j] - 1;
        // the cell is not attributed to any region. ignore it!
        if (regionIdx < 0) {
            continue;
        }

        assert(regionIdx < static_cast<int>(maxNumberOfRegions));
        totals[regionIdx] += property[j];
    }

    for (std::size_t i = 0; i < maxNumberOfRegions; ++i) {
        totals[i] = comm.sum(totals[i]);
    }

    return totals;
}

template<class FluidSystem>
void GenericOutputBlackoilModule<FluidSystem>::
doAllocBuffers(const unsigned bufferSize,
               const unsigned reportStepNum,
               const bool     substep,
               const bool     log,
               const bool     isRestart,
               const EclHysteresisConfig* hysteresisConfig,
               const unsigned numOutputNnc,
               std::map<std::string, int> rstKeywords)
{
    if (rstKeywords.empty()) {
        rstKeywords = schedule_.rst_keywords(reportStepNum);
    }

    // Output RESTART_OPM_EXTENDED only when explicitly requested by user.
    for (auto& [keyword, should_write] : rstKeywords) {
        if (this->isOutputCreationDirective_(keyword) || keyword == "WELLS") {
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

    // We always output oil pressure
    rstKeywords["PRES"] = 0;

    // Fluid in place
    this->computeFip_ = this->fipC_.allocate(bufferSize,
                                             summaryConfig_,
                                             !substep,
                                             rstKeywords);

    const auto needAvgPress = !substep         ||
        !this->RPRNodes_.empty()               ||
        this->summaryConfig_.hasKeyword("FPR") ||
        this->summaryConfig_.hasKeyword("FPRP");

    const auto needPoreVolume = needAvgPress    ||
        this->summaryConfig_.hasKeyword("FHPV") ||
        this->summaryConfig_.match("RHPV*");

    if (needPoreVolume) {
        this->fipC_.add(Inplace::Phase::PoreVolume);
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
        this->rftC_.allocate(reportStepNum);
    }

    const bool alloc_fields = isRestart || (schedule_.write_rst_file(reportStepNum) && !substep);
    this->flowsC_.allocate(bufferSize, numOutputNnc, alloc_fields, rstKeywords);

    // Field data should be allocated
    // 1) When we want to restart
    // 2) When it is ask for by the user via restartConfig
    // 3) When it is not a substep
    if (!alloc_fields) {
        return;
    }

    enum class EntryPhaseType {
        // one kw for all fields
        // RESIDUAL
        None,

        // append GWO to kw name
        // KR(G|W|O)
        GWO,

        // Pure kw controls all, additionally append GWO to kw name for each phase
        // DEN, DENG, DENW, DENO
        NGWO,

        // append gas/wat/oil to kw name
        // SGAS, SWAT, SOIL
        GasWatOil,

        // Pure kw controls all, then first letter of kw and apply gas/wat/oil for each phase.
        // VISC, VGAS, VWAT, VOIL
        NGasWatOil,
    };

    struct Entry {
        std::variant<ScalarBuffer*,
                     std::array<ScalarBuffer, numPhases>*,
                     std::vector<ScalarBuffer*>*> data; // Data to resize
        std::string_view kw; // Array keyword
        bool supported; // True if array output is supported
        bool required = true; // True if array keyword is required
        EntryPhaseType phaseType = EntryPhaseType::GasWatOil; // How phase is encoded in keyword
    };

    const auto& oilvap = schedule_[std::max(reportStepNum, 0u)].oilvap();
    const auto& simConfig = eclState_.getSimulationConfig();
    using OilVapP = OilVaporizationProperties::OilVaporization;

    auto pbpd_fields = std::vector{
        &bubblePointPressure_,
        &dewPointPressure_
    };
    auto rockc_fields = std::vector{
        &rockCompPorvMultiplier_,
        &rockCompTransMultiplier_,
        &swMax_,
        &minimumOilPressure_,
        &overburdenPressure_
    };

    const auto entries = std::array{
       Entry{&saturation_,                       "S", true},
       Entry{&fluidPressure_,             "PRESSURE", true},
       // If TEMP is set in RPTRST we output temperature even if THERMAL
       // is not activated
       Entry{&temperature_,                   "TEMP", enableEnergy_ || rstKeywords["TEMP"] > 0},
       Entry{&rs_,                              "RS", FluidSystem::enableDissolvedGas()},
       Entry{&rsw_,                            "RSW", FluidSystem::enableDissolvedGasInWater()},
       Entry{&rv_,                              "RV", FluidSystem::enableVaporizedOil()},
       Entry{&rvw_,                            "RVW", FluidSystem::enableVaporizedWater()},
       Entry{&drsdtcon_,                          "", oilvap.drsdtConvective()},
       Entry{&sSol_,                              "", enableSolvent_},
       Entry{&rswSol_,                            "", enableSolvent_ && simConfig.hasDISGASW()},
       Entry{&cPolymer_,                          "", enablePolymer_},
       Entry{&cFoam_,                             "", enableFoam_},
       Entry{&cSalt_,                             "", enableBrine_},
       Entry{&pSalt_,                             "", enableSaltPrecipitation_},
       Entry{&permFact_,                          "", enableSaltPrecipitation_ || enableBioeffects_},
       Entry{&soMax_,                             "", oilvap.getType() == OilVapP::VAPPARS},
       Entry{&soMax_,                             "", hysteresisConfig &&
                                                      hysteresisConfig->enableNonWettingHysteresis() &&
                                                      FluidSystem::phaseIsActive(oilPhaseIdx) &&
                                                      FluidSystem::phaseIsActive(waterPhaseIdx)},
       Entry{&sgmax_,                             "", hysteresisConfig &&
                                                      hysteresisConfig->enableNonWettingHysteresis() &&
                                                      FluidSystem::phaseIsActive(oilPhaseIdx) &&
                                                      FluidSystem::phaseIsActive(gasPhaseIdx)},
       Entry{&swMax_,                             "", hysteresisConfig &&
                                                      hysteresisConfig->enableWettingHysteresis() &&
                                                      FluidSystem::phaseIsActive(oilPhaseIdx) &&
                                                      FluidSystem::phaseIsActive(waterPhaseIdx)},
       Entry{&shmax_,                             "", hysteresisConfig &&
                                                      hysteresisConfig->enableWettingHysteresis() &&
                                                      FluidSystem::phaseIsActive(oilPhaseIdx) &&
                                                      FluidSystem::phaseIsActive(gasPhaseIdx)},
       Entry{&swmin_,                             "", hysteresisConfig &&
                                                      hysteresisConfig->enablePCHysteresis() &&
                                                      FluidSystem::phaseIsActive(oilPhaseIdx) &&
                                                      FluidSystem::phaseIsActive(waterPhaseIdx)},
       Entry{&somin_,                             "", hysteresisConfig &&
                                                      hysteresisConfig->enablePCHysteresis() &&
                                                      FluidSystem::phaseIsActive(oilPhaseIdx) &&
                                                      FluidSystem::phaseIsActive(gasPhaseIdx)},
       Entry{&ppcw_,                          "PPCW", eclState_.fieldProps().has_double("SWATINIT")},
       Entry{&gasDissolutionFactor_,         "RSSAT", FluidSystem::enableDissolvedGas(), false},
       Entry{&oilVaporizationFactor_,        "RVSAT", FluidSystem::enableVaporizedOil(), false},
       Entry{&gasDissolutionFactorInWater_, "RSWSAT", FluidSystem::enableDissolvedGasInWater(), false},
       Entry{&waterVaporizationFactor_,     "RVWSAT", FluidSystem::enableVaporizedWater(), false},
       Entry{&invB_,                             "B", true, false, EntryPhaseType::GWO},
       Entry{&rPorV_,                        "RPORV", true, false},
       Entry{&density_,                        "DEN", true, false, EntryPhaseType::NGWO},
       Entry{&viscosity_,                     "VISC", true, false, EntryPhaseType::NGasWatOil},
       Entry{&relativePermeability_,            "KR", true, false, EntryPhaseType::GWO},
       Entry{&pcog_,                          "PCOG", FluidSystem::phaseIsActive(oilPhaseIdx) &&
                                                      FluidSystem::phaseIsActive(gasPhaseIdx), false},
       Entry{&pcgw_,                          "PCGW", FluidSystem::phaseIsActive(gasPhaseIdx) &&
                                                      FluidSystem::phaseIsActive(waterPhaseIdx), false},
       Entry{&pcow_,                          "PCOW", FluidSystem::phaseIsActive(oilPhaseIdx) &&
                                                      FluidSystem::phaseIsActive(waterPhaseIdx), false},
       Entry{&pbpd_fields,                    "PBPD", true, false},
       Entry{&residual_,                  "RESIDUAL", true, false, EntryPhaseType::None},
       Entry{&rockc_fields,                  "ROCKC", true, false},
    };

    auto handleScalarEntry =
        [&rstKeywords, bufferSize](ScalarBuffer& data,
                                   const std::string_view& kw,
                                   const bool supported,
                                   const bool required)
    {
        int dummy = 1;
        auto& runtimeRequested = kw.empty() ? dummy : rstKeywords[std::string(kw)];
        if (supported && (required || runtimeRequested > 0)) {
            data.resize(bufferSize, 0.0);
            runtimeRequested = 0;
            return true;
        }
        return false;
    };

    auto getName = [](std::string_view kw, EntryPhaseType type, int phase)
    {
        constexpr auto phaseName = std::array{
            "WAT",
            "OIL",
            "GAS",
        };

        // This method assumes blackoil ordering of phases.
        // Apply a reordering for other fluid systems (compositional).
        PhaseReorder<FluidSystem> reorder;
        phase = reorder(phase);

        switch (type) {
        case EntryPhaseType::None:
            return std::string(kw);

        case EntryPhaseType::NGWO:
        case EntryPhaseType::GWO:
            return std::string(kw) + std::string_view{"WOG"}[phase];

        case EntryPhaseType::NGasWatOil:
            return std::string(1,kw[0]) + phaseName[phase];
        case EntryPhaseType::GasWatOil:
            return std::string(kw) + phaseName[phase];

        default:
            assert(false); // should never be here
            return std::string{};
        }
    };

    using PhaseArray = std::array<ScalarBuffer,numPhases>;
    std::for_each(entries.begin(), entries.end(),
                  [&handleScalarEntry, &getName, &rstKeywords](const auto& entry)
                  {
                      std::visit(VisitorOverloadSet{
                                    // simple scalar entry
                                    [&entry, &handleScalarEntry](ScalarBuffer* v)
                                    {
                                        handleScalarEntry(*v,
                                                          std::string(entry.kw),
                                                          entry.supported,
                                                          entry.required);
                                    },
                                    // multiple outputs controlled by one keyword
                                    [&entry, &handleScalarEntry](std::vector<ScalarBuffer*>* V)
                                    {
                                        auto& v = *V;
                                        if (handleScalarEntry(*v[0],
                                                              std::string(entry.kw),
                                                              entry.supported,
                                                              entry.required))
                                        {
                                            for (std::size_t i = 1; i < v.size(); ++i) {
                                                handleScalarEntry(*v[i], "", true, true);
                                            }
                                        }
                                    },
                                    // one entry per phase
                                    [&entry, &handleScalarEntry, &getName, &rstKeywords](PhaseArray* v)
                                    {
                                       constexpr auto phases = std::array{
                                          waterPhaseIdx,
                                          oilPhaseIdx,
                                          gasPhaseIdx,
                                       };

                                       bool required = entry.required;

                                       // we need to check if the main keyword is defined. If so,
                                       // call the handler with required set to true
                                       if (!required && (entry.phaseType == EntryPhaseType::NGWO ||
                                                         entry.phaseType == EntryPhaseType::NGasWatOil ||
                                                         entry.phaseType == EntryPhaseType::None))
                                       {
                                           auto it = rstKeywords.find(std::string(entry.kw));
                                           if (it != rstKeywords.end() && it->second > 0) {
                                               required = true;
                                               it->second = 0;
                                               if (entry.phaseType != EntryPhaseType::None) {
                                                   // Register phase entries as handled in map
                                                   for (int phase : phases) {
                                                       if (FluidSystem::phaseIsActive(phase)) {
                                                           rstKeywords[getName(entry.kw, entry.phaseType, phase)] = 0;
                                                       }
                                                   }
                                               }
                                           }
                                       }

                                       for (int phase : phases) {
                                           if (FluidSystem::phaseIsActive(phase)) {
                                               handleScalarEntry((*v)[phase],
                                                                 getName(entry.kw, entry.phaseType, phase),
                                                                 entry.supported,
                                                                 required);
                                            }
                                       }
                                    }
                                 }, entry.data);
                  });

    if (enableMech_ && eclState_.runspec().mech()) {
        this->mech_.allocate(bufferSize, rstKeywords);
    }

    if (enableExtbo_) {
        extboC_.allocate(bufferSize);
    }

    if (enableBioeffects_) {
        // Biofilms for gas-water systems; MICP only for water systems  
        this->bioeffectsC_.allocate(bufferSize, !FluidSystem::phaseIsActive(gasPhaseIdx));
    }

    // tracers
    this->tracerC_.allocate(bufferSize, eclState_.tracer());

    //Warn for any unhandled keyword
    if (log) {
        std::for_each(rstKeywords.begin(), rstKeywords.end(),
                      [](const auto& kwpair)
                      {
                          const auto& [kw, value] = kwpair;
                          if (value > 0) {
                              OpmLog::warning(
                                  fmt::format("Keyword '{}' is unhandled for output to restart file.", kw)
                              );
                          }
                      });
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

template<class FluidSystem>
bool GenericOutputBlackoilModule<FluidSystem>::
isOutputCreationDirective_(const std::string& keyword)
{
    return (keyword == "BASIC") || (keyword == "FREQ")
        || (keyword == "RESTART")                        // From RPTSCHED
        || (keyword == "SAVE")  || (keyword == "SFREQ"); // Not really supported
}

template<class FluidSystem>
void GenericOutputBlackoilModule<FluidSystem>::
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

template<class FluidSystem>
int GenericOutputBlackoilModule<FluidSystem>::
regionMax(const std::vector<int>& region,
          const Parallel::Communication& comm)
{
    const auto max_value = region.empty() ? 0 : *std::max_element(region.begin(), region.end());
    return comm.max(max_value);
}

template<class FluidSystem>
void GenericOutputBlackoilModule<FluidSystem>::
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

template<class FluidSystem>
void GenericOutputBlackoilModule<FluidSystem>::
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
        update_inplace(phase, this->fipC_.get(phase));
    }
}

template<class FluidSystem>
Inplace GenericOutputBlackoilModule<FluidSystem>::
accumulateRegionSums(const Parallel::Communication& comm)
{
    Inplace inplace;

    for (const auto& region : this->regions_) {
        makeRegionSum(inplace, region.first, comm);
    }

    return inplace;
}

template<class FluidSystem>
typename GenericOutputBlackoilModule<FluidSystem>::Scalar
GenericOutputBlackoilModule<FluidSystem>::
sum(const ScalarBuffer& v)
{
    return std::accumulate(v.begin(), v.end(), Scalar{0});
}

template<class FluidSystem>
void GenericOutputBlackoilModule<FluidSystem>::
updateSummaryRegionValues(const Inplace& inplace,
                          std::map<std::string, double>& miscSummaryData,
                          std::map<std::string, std::vector<double>>& regionData) const
{
    // The field summary vectors should only use the FIPNUM based region sum.
    {
        for (const auto& phase : Inplace::phases()) {
            const std::string key = "F" + Inplace::EclString(phase);
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

template<class FluidSystem>
void GenericOutputBlackoilModule<FluidSystem>::
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

template<class FluidSystem>
void GenericOutputBlackoilModule<FluidSystem>::
setupExtraBlockData(const std::size_t        reportStepNum,
                    std::function<bool(int)> isCartIdxOnThisRank)
{
    const auto& sched = this->schedule_[reportStepNum - 1];

    for (const auto& wname : sched.well_order()) {
        const auto& well = sched.wells.get(wname);
        for (const auto& connection : well.getConnections()) {
            if (isCartIdxOnThisRank(connection.global_index())) {
                this->extraBlockData_.emplace(std::piecewise_construct,
                                              std::forward_as_tuple("BPR",
                                                                    connection.global_index() + 1),
                                              std::forward_as_tuple(0.0));
            }
        }
    }
}

template<class FluidSystem>
void GenericOutputBlackoilModule<FluidSystem>::
assignGlobalFieldsToSolution(data::Solution& sol)
{
    this->rst_conv_.outputRestart(sol);
}

template<class T> using FS = BlackOilFluidSystem<T, BlackOilDefaultFluidSystemIndices>;

#define INSTANTIATE_TYPE(T) \
    template class GenericOutputBlackoilModule<FS<T>>;

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

#define INSTANTIATE_COMP_THREEPHASE(NUM) \
    template<class T> using FS##NUM = GenericOilGasWaterFluidSystem<T, NUM, true>; \
    template class GenericOutputBlackoilModule<FS##NUM<double>>;

#define INSTANTIATE_COMP_TWOPHASE(NUM) \
    template<class T> using GFS##NUM = GenericOilGasWaterFluidSystem<T, NUM, false>; \
    template class GenericOutputBlackoilModule<GFS##NUM<double>>;

#define INSTANTIATE_COMP(NUM) \
    INSTANTIATE_COMP_THREEPHASE(NUM) \
    INSTANTIATE_COMP_TWOPHASE(NUM)

INSTANTIATE_COMP_THREEPHASE(0)  // \Note: to register the parameter ForceDisableFluidInPlaceOutput
INSTANTIATE_COMP(2)
INSTANTIATE_COMP(3)
INSTANTIATE_COMP(4)
INSTANTIATE_COMP(5)
INSTANTIATE_COMP(6)
INSTANTIATE_COMP(7)

} // namespace Opm
