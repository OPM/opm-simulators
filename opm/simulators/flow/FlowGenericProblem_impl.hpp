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
#ifndef OPM_FLOW_GENERIC_PROBLEM_IMPL_HPP
#define OPM_FLOW_GENERIC_PROBLEM_IMPL_HPP

#ifndef OPM_FLOW_GENERIC_PROBLEM_HPP
#include <config.h>
#include <opm/simulators/flow/FlowGenericProblem.hpp>
#endif

#include <dune/common/parametertree.hh>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Tables/OverburdTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/RockwnodTable.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/models/discretization/common/fvbaseparameters.hh>
#include <opm/models/utils/basicparameters.hh>
#include <opm/models/utils/parametersystem.hpp>

#include <opm/simulators/flow/FlowProblemParameters.hpp>
#include <opm/simulators/flow/SolutionContainers.hpp>

#include <opm/simulators/timestepping/EclTimeSteppingParams.hpp>

#include <boost/date_time.hpp>

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <iostream>
#include <stdexcept>

namespace Opm {

template<class GridView, class FluidSystem>
FlowGenericProblem<GridView,FluidSystem>::
FlowGenericProblem(const EclipseState& eclState,
                   const Schedule& schedule,
                   const GridView& gridView)
    : eclState_(eclState)
    , schedule_(schedule)
    , gridView_(gridView)
    , lookUpData_(gridView)
{
    // we need to update the FluidSystem based on EclipseState before it is passed around
    this->initFluidSystem_();

    enableTuning_ = Parameters::Get<Parameters::EnableTuning>();
    enableDriftCompensation_ = Parameters::Get<Parameters::EnableDriftCompensation>();
    initialTimeStepSize_ = Parameters::Get<Parameters::InitialTimeStepSize<Scalar>>();
    maxTimeStepAfterWellEvent_ = unit::convert::from
        (Parameters::Get<Parameters::TimeStepAfterEventInDays<Scalar>>(), unit::day);

    // The value N for this parameter is defined in the following order of precedence:
    //
    // 1. Command line value (--num-pressure-points-equil=N)
    //
    // 2. EQLDIMS item 2.  Default value from
    //    opm-common/opm/input/eclipse/share/keywords/000_Eclipse100/E/EQLDIMS

    numPressurePointsEquil_ = Parameters::IsSet<Parameters::NumPressurePointsEquil>()
        ? Parameters::Get<Parameters::NumPressurePointsEquil>()
        : eclState.getTableManager().getEqldims().getNumDepthNodesP();

    explicitRockCompaction_ = Parameters::Get<Parameters::ExplicitRockCompaction>();
}

template<class GridView, class FluidSystem>
FlowGenericProblem<GridView,FluidSystem>
FlowGenericProblem<GridView,FluidSystem>::
serializationTestObject(const EclipseState& eclState,
                        const Schedule& schedule,
                        const GridView& gridView)
{
    FlowGenericProblem result(eclState, schedule, gridView);
    result.maxOilSaturation_ = {1.0, 2.0};
    result.maxWaterSaturation_ = {6.0};
    result.minRefPressure_ = {7.0, 8.0, 9.0, 10.0};
    result.overburdenPressure_ = {11.0};
    result.solventSaturation_ = {15.0};
    result.solventRsw_ = {18.0};
    result.polymer_ = PolymerSolutionContainer<Scalar>::serializationTestObject();
    result.micp_ = MICPSolutionContainer<Scalar>::serializationTestObject();

    return result;
}

template<class GridView, class FluidSystem>
std::string
FlowGenericProblem<GridView,FluidSystem>::
helpPreamble(int,
             const char **argv)
{
    std::string desc = FlowGenericProblem::briefDescription();
    if (!desc.empty())
        desc = desc + "\n";

    return
        "Usage: "+std::string(argv[0]) + " [OPTIONS] [ECL_DECK_FILENAME]\n"
        + desc;
}

template<class GridView, class FluidSystem>
std::string
FlowGenericProblem<GridView,FluidSystem>::
briefDescription()
{
    return briefDescription_;
}

template<class GridView, class FluidSystem>
void FlowGenericProblem<GridView,FluidSystem>::
readRockParameters_(const std::vector<Scalar>& cellCenterDepths,
                    std::function<std::array<int,3>(const unsigned)> ijkIndex)
{
    const auto& rock_config = eclState_.getSimulationConfig().rock_config();

    // read the rock compressibility parameters
    {
        const auto& comp = rock_config.comp();
        rockParams_.clear();
        std::transform(comp.begin(), comp.end(), std::back_inserter(rockParams_),
                       [](const auto& c)
                       {
                            return RockParams{static_cast<Scalar>(c.pref),
                                              static_cast<Scalar>(c.compressibility)};
                       });
    }

    // Warn that ROCK and ROCKOPTS item 2 = STORE is used together
    if (rock_config.store()) {
        OpmLog::warning("ROCKOPTS item 2 set to STORE, ROCK item 1 replaced with initial (equilibrated) pressures");
    }

    // read the parameters for water-induced rock compaction
    readRockCompactionParameters_();

    unsigned numElem = gridView_.size(0);
    if (eclState_.fieldProps().has_int(rock_config.rocknum_property())) {
        // Auxiliary function to check rockTableIdx_ values belong to the right range. Otherwise, throws.
        std::function<void(int, int)> valueCheck = [&ijkIndex,&rock_config,this](int fieldPropValue, int coarseElemIdx)
        {
            auto fmtError = [fieldPropValue, coarseElemIdx,&ijkIndex,&rock_config](const char* type, std::size_t size)
            {
                return fmt::format("{} table index {} for elem {} read from {}"
                                   " is out of bounds for number of tables {}",
                                   type,  fieldPropValue,
                                   ijkIndex(coarseElemIdx),
                                   rock_config.rocknum_property(), size);
            };
            if (!rockCompPoroMult_.empty() &&
                fieldPropValue > static_cast<int>(rockCompPoroMult_.size())) {
                throw std::runtime_error(fmtError("Rock compaction",
                                                  rockCompPoroMult_.size()));
            }
            if (!rockCompPoroMultWc_.empty() &&
                fieldPropValue > static_cast<int>(rockCompPoroMultWc_.size())) {
                throw std::runtime_error(fmtError("Rock water compaction",
                                                  rockCompPoroMultWc_.size()));
            }
        };

        rockTableIdx_ = this->lookUpData_.template assignFieldPropsIntOnLeaf<short unsigned int>(eclState_.fieldProps(),
                                                                                                 rock_config.rocknum_property(),
                                                                                                 true /*needsTranslation*/,
                                                                                                 valueCheck);
    }

    // Store overburden pressure pr element
    const auto& overburdTables = eclState_.getTableManager().getOverburdTables();
    if (!overburdTables.empty() && !rock_config.store()) {
        overburdenPressure_.resize(numElem,0.0);
        std::size_t numRocktabTables = rock_config.num_rock_tables();

        if (overburdTables.size() != numRocktabTables)
            throw std::runtime_error(std::to_string(numRocktabTables) +" OVERBURD tables is expected, but " + std::to_string(overburdTables.size()) +" is provided");

        std::vector<Tabulated1DFunction<Scalar>> overburdenTables(numRocktabTables);
        for (std::size_t regionIdx = 0; regionIdx < numRocktabTables; ++regionIdx) {
            const OverburdTable& overburdTable =  overburdTables.template getTable<OverburdTable>(regionIdx);
            overburdenTables[regionIdx].setXYContainers(overburdTable.getDepthColumn(),overburdTable.getOverburdenPressureColumn());
        }

        for (std::size_t elemIdx = 0; elemIdx < numElem; ++ elemIdx) {
            unsigned tableIdx = 0;
            if (!rockTableIdx_.empty()) {
                tableIdx = rockTableIdx_[elemIdx];
            }
            overburdenPressure_[elemIdx] =
                overburdenTables[tableIdx].eval(cellCenterDepths[elemIdx], /*extrapolation=*/true);
        }
    }
    else if (!overburdTables.empty() && rock_config.store()) {
        OpmLog::warning("ROCKOPTS item 2 set to STORE, OVERBURD ignored!");
    } 
}

template<class GridView, class FluidSystem>
void FlowGenericProblem<GridView,FluidSystem>::
readRockCompactionParameters_()
{
    const auto& rock_config = eclState_.getSimulationConfig().rock_config();

    if (!rock_config.active())
        return; // deck does not enable rock compaction

    unsigned numElem = gridView_.size(0);
    switch (rock_config.hysteresis_mode()) {
    case RockConfig::Hysteresis::REVERS:
        break;
    case RockConfig::Hysteresis::IRREVERS:
        // interpolate the porv volume multiplier using the minimum pressure in the cell
        // i.e. don't allow re-inflation.
        minRefPressure_.resize(numElem, 1e99);
        break;
    default:
        throw std::runtime_error("Not support ROCKOMP hysteresis option ");
    }

    std::size_t numRocktabTables = rock_config.num_rock_tables();
    bool waterCompaction = rock_config.water_compaction();

    if (!waterCompaction) {
        const auto& rocktabTables = eclState_.getTableManager().getRocktabTables();
        if (rocktabTables.size() != numRocktabTables)
            throw std::runtime_error("ROCKCOMP is activated." + std::to_string(numRocktabTables)
                                     +" ROCKTAB tables is expected, but " + std::to_string(rocktabTables.size()) +" is provided");

        rockCompPoroMult_.resize(numRocktabTables);
        rockCompTransMult_.resize(numRocktabTables);
        for (std::size_t regionIdx = 0; regionIdx < numRocktabTables; ++regionIdx) {
            const auto& rocktabTable = rocktabTables.template getTable<RocktabTable>(regionIdx);
            const auto& pressureColumn = rocktabTable.getPressureColumn();
            const auto& poroColumn = rocktabTable.getPoreVolumeMultiplierColumn();
            const auto& transColumn = rocktabTable.getTransmissibilityMultiplierColumn();
            rockCompPoroMult_[regionIdx].setXYContainers(pressureColumn, poroColumn);
            rockCompTransMult_[regionIdx].setXYContainers(pressureColumn, transColumn);
        }
    } else {
        const auto& rock2dTables = eclState_.getTableManager().getRock2dTables();
        const auto& rock2dtrTables = eclState_.getTableManager().getRock2dtrTables();
        const auto& rockwnodTables = eclState_.getTableManager().getRockwnodTables();
        maxWaterSaturation_.resize(numElem, 0.0);

        if (rock2dTables.size() != numRocktabTables)
            throw std::runtime_error("Water compation option is selected in ROCKCOMP." + std::to_string(numRocktabTables)
                                     +" ROCK2D tables is expected, but " + std::to_string(rock2dTables.size()) +" is provided");

        if (rockwnodTables.size() != numRocktabTables)
            throw std::runtime_error("Water compation option is selected in ROCKCOMP." + std::to_string(numRocktabTables)
                                     +" ROCKWNOD tables is expected, but " + std::to_string(rockwnodTables.size()) +" is provided");
        //TODO check size match
        rockCompPoroMultWc_.resize(numRocktabTables, TabulatedTwoDFunction(TabulatedTwoDFunction::InterpolationPolicy::Vertical));
        for (std::size_t regionIdx = 0; regionIdx < numRocktabTables; ++regionIdx) {
            const RockwnodTable& rockwnodTable =  rockwnodTables.template getTable<RockwnodTable>(regionIdx);
            const auto& rock2dTable = rock2dTables[regionIdx];

            if (rockwnodTable.getSaturationColumn().size() != rock2dTable.sizeMultValues())
                throw std::runtime_error("Number of entries in ROCKWNOD and ROCK2D needs to match.");

            for (std::size_t xIdx = 0; xIdx < rock2dTable.size(); ++xIdx) {
                rockCompPoroMultWc_[regionIdx].appendXPos(rock2dTable.getPressureValue(xIdx));
                for (std::size_t yIdx = 0; yIdx < rockwnodTable.getSaturationColumn().size(); ++yIdx)
                    rockCompPoroMultWc_[regionIdx].appendSamplePoint(xIdx,
                                                                       rockwnodTable.getSaturationColumn()[yIdx],
                                                                       rock2dTable.getPvmultValue(xIdx, yIdx));
            }
        }

        if (!rock2dtrTables.empty()) {
            rockCompTransMultWc_.resize(numRocktabTables, TabulatedTwoDFunction(TabulatedTwoDFunction::InterpolationPolicy::Vertical));
            for (std::size_t regionIdx = 0; regionIdx < numRocktabTables; ++regionIdx) {
                const RockwnodTable& rockwnodTable =  rockwnodTables.template getTable<RockwnodTable>(regionIdx);
                const auto& rock2dtrTable = rock2dtrTables[regionIdx];

                if (rockwnodTable.getSaturationColumn().size() != rock2dtrTable.sizeMultValues())
                    throw std::runtime_error("Number of entries in ROCKWNOD and ROCK2DTR needs to match.");

                for (std::size_t xIdx = 0; xIdx < rock2dtrTable.size(); ++xIdx) {
                    rockCompTransMultWc_[regionIdx].appendXPos(rock2dtrTable.getPressureValue(xIdx));
                    for (std::size_t yIdx = 0; yIdx < rockwnodTable.getSaturationColumn().size(); ++yIdx)
                        rockCompTransMultWc_[regionIdx].appendSamplePoint(xIdx,
                                                                                 rockwnodTable.getSaturationColumn()[yIdx],
                                                                                 rock2dtrTable.getTransMultValue(xIdx, yIdx));
                }
            }
        }
    }
}

template<class GridView, class FluidSystem>
typename FlowGenericProblem<GridView,FluidSystem>::Scalar
FlowGenericProblem<GridView,FluidSystem>::
rockCompressibility(unsigned globalSpaceIdx) const
{
    if (this->rockParams_.empty())
        return 0.0;

    unsigned tableIdx = 0;
    if (!this->rockTableIdx_.empty()) {
        tableIdx = this->rockTableIdx_[globalSpaceIdx];
    }
    return this->rockParams_[tableIdx].compressibility;
}

template<class GridView, class FluidSystem>
typename FlowGenericProblem<GridView,FluidSystem>::Scalar
FlowGenericProblem<GridView,FluidSystem>::
porosity(unsigned globalSpaceIdx, unsigned timeIdx) const
{
    return this->referencePorosity_[timeIdx][globalSpaceIdx];
}

template<class GridView, class FluidSystem>
typename FlowGenericProblem<GridView,FluidSystem>::Scalar
FlowGenericProblem<GridView,FluidSystem>::
rockFraction(unsigned elementIdx, unsigned timeIdx) const
{
    // the reference porosity is defined as the accumulated pore volume divided by the
    // geometric volume of the element. Note that it can
    // be larger than 1.0 if porevolume multipliers are used
    // to for instance implement larger boundary cells
    auto porosity = this->lookUpData_.fieldPropDouble(eclState_.fieldProps(), "PORO", elementIdx);
    return referencePorosity(elementIdx, timeIdx) / porosity * (1 - porosity);
}

template<class GridView, class FluidSystem>
template<class T>
void FlowGenericProblem<GridView,FluidSystem>::
updateNum(const std::string& name, std::vector<T>& numbers, std::size_t num_regions)
{
    if (!eclState_.fieldProps().has_int(name))
        return;

    std::function<void(T, int)> valueCheck = [num_regions,name](T fieldPropValue, [[maybe_unused]] int fieldPropIdx) {
        if ( fieldPropValue > (int)num_regions) {
            throw std::runtime_error("Values larger than maximum number of regions "
                                     + std::to_string(num_regions) + " provided in " + name);
        }
        if ( fieldPropValue <= 0) {
            throw std::runtime_error("zero or negative values provided for region array: " + name);
        }
    };

    numbers = this->lookUpData_.template assignFieldPropsIntOnLeaf<T>(eclState_.fieldProps(), name,
                                                                      true /*needsTranslation*/, valueCheck);
}

template<class GridView, class FluidSystem>
void FlowGenericProblem<GridView,FluidSystem>::
updatePvtnum_()
{
    const auto num_regions = eclState_.getTableManager().getTabdims().getNumPVTTables();
    updateNum("PVTNUM", pvtnum_, num_regions);
}

template<class GridView, class FluidSystem>
void FlowGenericProblem<GridView,FluidSystem>::
updateSatnum_()
{
    const auto num_regions = eclState_.getTableManager().getTabdims().getNumSatTables();
    updateNum("SATNUM", satnum_, num_regions);
}

template<class GridView, class FluidSystem>
void FlowGenericProblem<GridView,FluidSystem>::
updateMiscnum_()
{
    const auto num_regions = 1; // we only support single region
    updateNum("MISCNUM", miscnum_, num_regions);
}

template<class GridView, class FluidSystem>
void FlowGenericProblem<GridView,FluidSystem>::
updatePlmixnum_()
{
    const auto num_regions = 1; // we only support single region
    updateNum("PLMIXNUM", plmixnum_, num_regions);
}

template<class GridView, class FluidSystem>
bool FlowGenericProblem<GridView,FluidSystem>::
vapparsActive(int episodeIdx) const
{
    const auto& oilVaporizationControl = schedule_[episodeIdx].oilvap();
    return (oilVaporizationControl.getType() == OilVaporizationProperties::OilVaporization::VAPPARS);
}

template<class GridView, class FluidSystem>
bool FlowGenericProblem<GridView,FluidSystem>::
beginEpisode_(bool enableExperiments,
              int episodeIdx)
{
    if (enableExperiments && gridView_.comm().rank() == 0 && episodeIdx >= 0) {
        // print some useful information in experimental mode. (the production
        // simulator does this externally.)
        std::ostringstream ss;
        boost::posix_time::time_facet* facet = new boost::posix_time::time_facet("%d-%b-%Y");
        boost::posix_time::ptime curDateTime =
            boost::posix_time::from_time_t(schedule_.simTime(episodeIdx));
        ss.imbue(std::locale(std::locale::classic(), facet));
        ss << "Report step " << episodeIdx + 1
                  << "/" << schedule_.size() - 1
                  << " at day " << schedule_.seconds(episodeIdx)/(24*3600)
                  << "/" << schedule_.seconds(schedule_.size() - 1)/(24*3600)
                  << ", date = " << curDateTime.date()
                  << "\n ";
        OpmLog::info(ss.str());
    }

    const auto& events = schedule_[episodeIdx].events();

    // react to TUNING changes
    if (episodeIdx > 0 && enableTuning_ && events.hasEvent(ScheduleEvents::TUNING_CHANGE))
    {
        const auto& sched_state = schedule_[episodeIdx];
        const auto& tuning = sched_state.tuning();
        initialTimeStepSize_ = sched_state.max_next_tstep(enableTuning_);
        maxTimeStepAfterWellEvent_ = tuning.TMAXWC;
        return true;
    }

    return false;
}

template<class GridView, class FluidSystem>
void FlowGenericProblem<GridView,FluidSystem>::
beginTimeStep_(bool enableExperiments,
               int episodeIdx,
               int timeStepIndex,
               Scalar startTime,
               Scalar time,
               Scalar timeStepSize,
               Scalar endTime)
{
    if (enableExperiments && gridView_.comm().rank() == 0 && episodeIdx >= 0) {
        std::ostringstream ss;
        boost::posix_time::time_facet* facet = new boost::posix_time::time_facet("%d-%b-%Y");
        boost::posix_time::ptime date = boost::posix_time::from_time_t(startTime) + boost::posix_time::milliseconds(static_cast<long long>(time / prefix::milli));
        ss.imbue(std::locale(std::locale::classic(), facet));
        ss <<"\nTime step " << timeStepIndex << ", stepsize "
               << unit::convert::to(timeStepSize, unit::day) << " days,"
               << " at day " << (double)unit::convert::to(time, unit::day)
               << "/" << (double)unit::convert::to(endTime, unit::day)
               << ", date = " << date;
        OpmLog::info(ss.str());
    }
}

template<class GridView, class FluidSystem>
void FlowGenericProblem<GridView,FluidSystem>::
initFluidSystem_()
{
    FluidSystem::initFromState(eclState_, schedule_);
}

template<class GridView, class FluidSystem>
void FlowGenericProblem<GridView,FluidSystem>::
readBlackoilExtentionsInitialConditions_(std::size_t numDof,
                                         bool enableSolvent,
                                         bool enablePolymer,
                                         bool enablePolymerMolarWeight,
                                         bool enableMICP)
{
    auto getArray = [](const std::vector<double>& input)
    {
        if constexpr (std::is_same_v<Scalar,double>) {
            return input;
        } else {
            return std::vector<Scalar>{input.begin(), input.end()};
        }
    };

    if (enableSolvent) {
        if (eclState_.fieldProps().has_double("SSOL")) {
            solventSaturation_ = getArray(eclState_.fieldProps().get_double("SSOL"));
        } else {
            solventSaturation_.resize(numDof, 0.0);
        }

        solventRsw_.resize(numDof, 0.0);
    }

    if (enablePolymer) {
        if (eclState_.fieldProps().has_double("SPOLY")) {
            polymer_.concentration = getArray(eclState_.fieldProps().get_double("SPOLY"));
        } else {
            polymer_.concentration.resize(numDof, 0.0);
        }
    }

    if (enablePolymerMolarWeight) {
        if (eclState_.fieldProps().has_double("SPOLYMW")) {
            polymer_.moleWeight = getArray(eclState_.fieldProps().get_double("SPOLYMW"));
        } else {
            polymer_.moleWeight.resize(numDof, 0.0);
        }
    }

    if (enableMICP) {
        if (eclState_.fieldProps().has_double("SMICR")) {
            micp_.microbialConcentration = getArray(eclState_.fieldProps().get_double("SMICR"));
        } else {
            micp_.microbialConcentration.resize(numDof, 0.0);
        }
        if (eclState_.fieldProps().has_double("SOXYG")) {
            micp_.oxygenConcentration = getArray(eclState_.fieldProps().get_double("SOXYG"));
        } else {
            micp_.oxygenConcentration.resize(numDof, 0.0);
        }
        if (eclState_.fieldProps().has_double("SUREA")) {
            micp_.ureaConcentration = getArray(eclState_.fieldProps().get_double("SUREA"));
        } else {
            micp_.ureaConcentration.resize(numDof, 0.0);
        }
        if (eclState_.fieldProps().has_double("SBIOF")) {
            micp_.biofilmConcentration = getArray(eclState_.fieldProps().get_double("SBIOF"));
        } else {
            micp_.biofilmConcentration.resize(numDof, 0.0);
        }
        if (eclState_.fieldProps().has_double("SCALC")) {
            micp_.calciteConcentration = getArray(eclState_.fieldProps().get_double("SCALC"));
        } else {
            micp_.calciteConcentration.resize(numDof, 0.0);
        }
    }
}

template<class GridView, class FluidSystem>
typename FlowGenericProblem<GridView,FluidSystem>::Scalar
FlowGenericProblem<GridView,FluidSystem>::
maxWaterSaturation(unsigned globalDofIdx) const
{
    if (maxWaterSaturation_.empty())
        return 0.0;

    return maxWaterSaturation_[globalDofIdx];
}

template<class GridView, class FluidSystem>
typename FlowGenericProblem<GridView,FluidSystem>::Scalar
FlowGenericProblem<GridView,FluidSystem>::
minOilPressure(unsigned globalDofIdx) const
{
    if (minRefPressure_.empty())
        return 0.0;

    return minRefPressure_[globalDofIdx];
}

template<class GridView, class FluidSystem>
typename FlowGenericProblem<GridView,FluidSystem>::Scalar
FlowGenericProblem<GridView,FluidSystem>::
overburdenPressure(unsigned elementIdx) const
{
    if (overburdenPressure_.empty())
        return 0.0;

    return overburdenPressure_[elementIdx];
}

template<class GridView, class FluidSystem>
typename FlowGenericProblem<GridView,FluidSystem>::Scalar
FlowGenericProblem<GridView,FluidSystem>::
solventSaturation(unsigned elemIdx) const
{
    if (solventSaturation_.empty())
        return 0;

    return solventSaturation_[elemIdx];
}

template<class GridView, class FluidSystem>
typename FlowGenericProblem<GridView,FluidSystem>::Scalar
FlowGenericProblem<GridView,FluidSystem>::
solventRsw(unsigned elemIdx) const
{
    if (solventRsw_.empty())
        return 0;

    return solventRsw_[elemIdx];
}



template<class GridView, class FluidSystem>
typename FlowGenericProblem<GridView,FluidSystem>::Scalar
FlowGenericProblem<GridView,FluidSystem>::
polymerConcentration(unsigned elemIdx) const
{
    if (polymer_.concentration.empty()) {
        return 0;
    }

    return polymer_.concentration[elemIdx];
}

template<class GridView, class FluidSystem>
typename FlowGenericProblem<GridView,FluidSystem>::Scalar
FlowGenericProblem<GridView,FluidSystem>::
polymerMolecularWeight(const unsigned elemIdx) const
{
    if (polymer_.moleWeight.empty()) {
        return 0.0;
    }

    return polymer_.moleWeight[elemIdx];
}

template<class GridView, class FluidSystem>
typename FlowGenericProblem<GridView,FluidSystem>::Scalar
FlowGenericProblem<GridView,FluidSystem>::
microbialConcentration(unsigned elemIdx) const
{
    if (micp_.microbialConcentration.empty()) {
        return 0;
    }

    return micp_.microbialConcentration[elemIdx];
}

template<class GridView, class FluidSystem>
typename FlowGenericProblem<GridView,FluidSystem>::Scalar
FlowGenericProblem<GridView,FluidSystem>::
oxygenConcentration(unsigned elemIdx) const
{
    if (micp_.oxygenConcentration.empty()) {
        return 0;
    }

    return micp_.oxygenConcentration[elemIdx];
}

template<class GridView, class FluidSystem>
typename FlowGenericProblem<GridView,FluidSystem>::Scalar
FlowGenericProblem<GridView,FluidSystem>::
ureaConcentration(unsigned elemIdx) const
{
    if (micp_.ureaConcentration.empty()) {
        return 0;
    }

    return micp_.ureaConcentration[elemIdx];
}

template<class GridView, class FluidSystem>
typename FlowGenericProblem<GridView,FluidSystem>::Scalar
FlowGenericProblem<GridView,FluidSystem>::
biofilmConcentration(unsigned elemIdx) const
{
    if (micp_.biofilmConcentration.empty()) {
        return 0;
    }

    return micp_.biofilmConcentration[elemIdx];
}

template<class GridView, class FluidSystem>
typename FlowGenericProblem<GridView,FluidSystem>::Scalar
FlowGenericProblem<GridView,FluidSystem>::
calciteConcentration(unsigned elemIdx) const
{
    if (micp_.calciteConcentration.empty()) {
        return 0;
    }

    return micp_.calciteConcentration[elemIdx];
}

template<class GridView, class FluidSystem>
unsigned FlowGenericProblem<GridView,FluidSystem>::
pvtRegionIndex(unsigned elemIdx) const
{
    if (pvtnum_.empty())
        return 0;

    return pvtnum_[elemIdx];
}

template<class GridView, class FluidSystem>
unsigned FlowGenericProblem<GridView,FluidSystem>::
satnumRegionIndex(unsigned elemIdx) const
{
    if (satnum_.empty())
        return 0;

    return satnum_[elemIdx];
}

template<class GridView, class FluidSystem>
unsigned FlowGenericProblem<GridView,FluidSystem>::
miscnumRegionIndex(unsigned elemIdx) const
{
    if (miscnum_.empty())
        return 0;

    return miscnum_[elemIdx];
}

template<class GridView, class FluidSystem>
unsigned FlowGenericProblem<GridView,FluidSystem>::
plmixnumRegionIndex(unsigned elemIdx) const
{
    if (plmixnum_.empty())
        return 0;

    return plmixnum_[elemIdx];
}

template<class GridView, class FluidSystem>
typename FlowGenericProblem<GridView,FluidSystem>::Scalar
FlowGenericProblem<GridView,FluidSystem>::
maxPolymerAdsorption(unsigned elemIdx) const
{
    if (polymer_.maxAdsorption.empty()) {
        return 0;
    }

    return polymer_.maxAdsorption[elemIdx];
}

template<class GridView, class FluidSystem>
bool FlowGenericProblem<GridView,FluidSystem>::
operator==(const FlowGenericProblem& rhs) const
{
    return this->maxWaterSaturation_ == rhs.maxWaterSaturation_ &&
           this->minRefPressure_ == rhs.minRefPressure_ &&
           this->overburdenPressure_ == rhs.overburdenPressure_ &&
           this->solventSaturation_ == rhs.solventSaturation_ &&
           this->solventRsw_ == rhs.solventRsw_ &&
           this->polymer_ == rhs.polymer_ &&
           this->micp_ == rhs.micp_;
}

} // namespace Opm

#endif // OPM_FLOW_GENERIC_PROBLEM_IMPL_HPP
