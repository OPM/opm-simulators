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
#include <ebos/eclgenericproblem.hh>

#include <ebos/eclalternativeblackoilindices.hh>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/OverburdTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/RockwnodTable.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/polyhedralgrid.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/common/gridpart2gridview.hh>
#include <ebos/femcpgridcompat.hh>
#endif

#include <boost/date_time.hpp>

#include <limits>
#include <stdexcept>
#include <iostream>

namespace Opm {

template<class GridView, class FluidSystem, class Scalar>
EclGenericProblem<GridView,FluidSystem,Scalar>::
EclGenericProblem(const EclipseState& eclState,
                  const Schedule& schedule,
                  const GridView& gridView)
    : eclState_(eclState)
    , schedule_(schedule)
    , gridView_(gridView)
{
}

template<class GridView, class FluidSystem, class Scalar>
std::string
EclGenericProblem<GridView,FluidSystem,Scalar>::
helpPreamble(int,
             const char **argv)
{
    std::string desc = EclGenericProblem::briefDescription();
    if (!desc.empty())
        desc = desc + "\n";

    return
        "Usage: "+std::string(argv[0]) + " [OPTIONS] [ECL_DECK_FILENAME]\n"
        + desc;
}

template<class GridView, class FluidSystem, class Scalar>
std::string
EclGenericProblem<GridView,FluidSystem,Scalar>::
briefDescription()
{
    if (briefDescription_.empty())
        return
            "The Ecl-deck Black-Oil reservoir Simulator (ebos); a hydrocarbon "
            "reservoir simulation program that processes ECL-formatted input "
            "files that is part of the Open Porous Media project "
            "(https://opm-project.org).\n"
            "\n"
            "THE GOAL OF THE `ebos` SIMULATOR IS TO CATER FOR THE NEEDS OF "
            "DEVELOPMENT AND RESEARCH. No guarantees are made for production use!";
    else
        return briefDescription_;
}

template<class GridView, class FluidSystem, class Scalar>
void EclGenericProblem<GridView,FluidSystem,Scalar>::
readRockParameters_(const std::vector<Scalar>& cellCenterDepths)
{
    const auto& rock_config = eclState_.getSimulationConfig().rock_config();

    // read the rock compressibility parameters
    {
        const auto& comp = rock_config.comp();
        rockParams_.clear();
        for (const auto& c : comp)
            rockParams_.push_back( { c.pref, c.compressibility } );
    }

    // read the parameters for water-induced rock compaction
    readRockCompactionParameters_();

    unsigned numElem = gridView_.size(0);
    if (eclState_.fieldProps().has_int(rock_config.rocknum_property())) {
        rockTableIdx_.resize(numElem);
        const auto& num = eclState_.fieldProps().get_int(rock_config.rocknum_property());
        for (size_t elemIdx = 0; elemIdx < numElem; ++ elemIdx) {
            rockTableIdx_[elemIdx] = num[elemIdx] - 1;
        }
    }

    // Store overburden pressure pr element
    const auto& overburdTables = eclState_.getTableManager().getOverburdTables();
    if (!overburdTables.empty()) {
        overburdenPressure_.resize(numElem,0.0);
        size_t numRocktabTables = rock_config.num_rock_tables();

        if (overburdTables.size() != numRocktabTables)
            throw std::runtime_error(std::to_string(numRocktabTables) +" OVERBURD tables is expected, but " + std::to_string(overburdTables.size()) +" is provided");

        std::vector<Tabulated1DFunction<Scalar>> overburdenTables(numRocktabTables);
        for (size_t regionIdx = 0; regionIdx < numRocktabTables; ++regionIdx) {
            const OverburdTable& overburdTable =  overburdTables.template getTable<OverburdTable>(regionIdx);
            overburdenTables[regionIdx].setXYContainers(overburdTable.getDepthColumn(),overburdTable.getOverburdenPressureColumn());
        }

        for (size_t elemIdx = 0; elemIdx < numElem; ++ elemIdx) {
            unsigned tableIdx = 0;
            if (!rockTableIdx_.empty()) {
                tableIdx = rockTableIdx_[elemIdx];
            }
            overburdenPressure_[elemIdx] = overburdenTables[tableIdx].eval(cellCenterDepths[elemIdx], /*extrapolation=*/true);
        }
    }
}

template<class GridView, class FluidSystem, class Scalar>
void EclGenericProblem<GridView,FluidSystem,Scalar>::
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
        minOilPressure_.resize(numElem, 1e99);
        break;
    default:
        throw std::runtime_error("Not support ROCKOMP hysteresis option ");
    }

    size_t numRocktabTables = rock_config.num_rock_tables();
    bool waterCompaction = rock_config.water_compaction();

    if (!waterCompaction) {
        const auto& rocktabTables = eclState_.getTableManager().getRocktabTables();
        if (rocktabTables.size() != numRocktabTables)
            throw std::runtime_error("ROCKCOMP is activated." + std::to_string(numRocktabTables)
                                     +" ROCKTAB tables is expected, but " + std::to_string(rocktabTables.size()) +" is provided");

        rockCompPoroMult_.resize(numRocktabTables);
        rockCompTransMult_.resize(numRocktabTables);
        for (size_t regionIdx = 0; regionIdx < numRocktabTables; ++regionIdx) {
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
        for (size_t regionIdx = 0; regionIdx < numRocktabTables; ++regionIdx) {
            const RockwnodTable& rockwnodTable =  rockwnodTables.template getTable<RockwnodTable>(regionIdx);
            const auto& rock2dTable = rock2dTables[regionIdx];

            if (rockwnodTable.getSaturationColumn().size() != rock2dTable.sizeMultValues())
                throw std::runtime_error("Number of entries in ROCKWNOD and ROCK2D needs to match.");

            for (size_t xIdx = 0; xIdx < rock2dTable.size(); ++xIdx) {
                rockCompPoroMultWc_[regionIdx].appendXPos(rock2dTable.getPressureValue(xIdx));
                for (size_t yIdx = 0; yIdx < rockwnodTable.getSaturationColumn().size(); ++yIdx)
                    rockCompPoroMultWc_[regionIdx].appendSamplePoint(xIdx,
                                                                       rockwnodTable.getSaturationColumn()[yIdx],
                                                                       rock2dTable.getPvmultValue(xIdx, yIdx));
            }
        }

        if (!rock2dtrTables.empty()) {
            rockCompTransMultWc_.resize(numRocktabTables, TabulatedTwoDFunction(TabulatedTwoDFunction::InterpolationPolicy::Vertical));
            for (size_t regionIdx = 0; regionIdx < numRocktabTables; ++regionIdx) {
                const RockwnodTable& rockwnodTable =  rockwnodTables.template getTable<RockwnodTable>(regionIdx);
                const auto& rock2dtrTable = rock2dtrTables[regionIdx];

                if (rockwnodTable.getSaturationColumn().size() != rock2dtrTable.sizeMultValues())
                    throw std::runtime_error("Number of entries in ROCKWNOD and ROCK2DTR needs to match.");

                for (size_t xIdx = 0; xIdx < rock2dtrTable.size(); ++xIdx) {
                    rockCompTransMultWc_[regionIdx].appendXPos(rock2dtrTable.getPressureValue(xIdx));
                    for (size_t yIdx = 0; yIdx < rockwnodTable.getSaturationColumn().size(); ++yIdx)
                        rockCompTransMultWc_[regionIdx].appendSamplePoint(xIdx,
                                                                                 rockwnodTable.getSaturationColumn()[yIdx],
                                                                                 rock2dtrTable.getTransMultValue(xIdx, yIdx));
                }
            }
        }
    }
}

template<class GridView, class FluidSystem, class Scalar>
template<class T>
void EclGenericProblem<GridView,FluidSystem,Scalar>::
updateNum(const std::string& name, std::vector<T>& numbers)
{
    if (!eclState_.fieldProps().has_int(name))
        return;

    const auto& numData = eclState_.fieldProps().get_int(name);

    unsigned numElems = gridView_.size(/*codim=*/0);
    numbers.resize(numElems);
    for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx) {
        numbers[elemIdx] = static_cast<T>(numData[elemIdx]) - 1;
    }
}

template<class GridView, class FluidSystem, class Scalar>
void EclGenericProblem<GridView,FluidSystem,Scalar>::
updatePvtnum_()
{
    updateNum("PVTNUM", pvtnum_);
}

template<class GridView, class FluidSystem, class Scalar>
void EclGenericProblem<GridView,FluidSystem,Scalar>::
updateSatnum_()
{
    updateNum("SATNUM", satnum_);
}

template<class GridView, class FluidSystem, class Scalar>
void EclGenericProblem<GridView,FluidSystem,Scalar>::
updateMiscnum_()
{
    updateNum("MISCNUM", miscnum_);
}

template<class GridView, class FluidSystem, class Scalar>
void EclGenericProblem<GridView,FluidSystem,Scalar>::
updatePlmixnum_()
{
    updateNum("PLMIXNUM", plmixnum_);
}

template<class GridView, class FluidSystem, class Scalar>
bool EclGenericProblem<GridView,FluidSystem,Scalar>::
vapparsActive(int episodeIdx) const
{
    const auto& oilVaporizationControl = schedule_[episodeIdx].oilvap();
    return (oilVaporizationControl.getType() == OilVaporizationProperties::OilVaporization::VAPPARS);
}

template<class GridView, class FluidSystem, class Scalar>
bool EclGenericProblem<GridView,FluidSystem,Scalar>::
drsdtActive_(int episodeIdx) const
{
    const auto& oilVaporizationControl = schedule_[episodeIdx].oilvap();
    const bool bothOilGasActive = FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) &&
                                  FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx);
    return (oilVaporizationControl.drsdtActive() && bothOilGasActive);
}

template<class GridView, class FluidSystem, class Scalar>
bool EclGenericProblem<GridView,FluidSystem,Scalar>::
drvdtActive_(int episodeIdx) const
{
    const auto& oilVaporizationControl = schedule_[episodeIdx].oilvap();
    const bool bothOilGasActive = FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) &&
                                  FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx);
    return (oilVaporizationControl.drvdtActive() && bothOilGasActive);
}

template<class GridView, class FluidSystem, class Scalar>
bool EclGenericProblem<GridView,FluidSystem,Scalar>::
drsdtConvective_(int episodeIdx) const
{
    const auto& oilVaporizationControl = schedule_[episodeIdx].oilvap();
    const bool bothOilGasActive = FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) &&
                                  FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx);
    return (oilVaporizationControl.drsdtConvective() && bothOilGasActive);
}

template<class GridView, class FluidSystem, class Scalar>
bool EclGenericProblem<GridView,FluidSystem,Scalar>::
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

    const auto& debug_config = schedule_[episodeIdx].debug_config();
    OpmLog::updateDebugConfig(debug_config);

    const auto& events = schedule_[episodeIdx].events();
    // react to TUNING changes
    if (episodeIdx > 0 && enableTuning_ && events.hasEvent(ScheduleEvents::TUNING_CHANGE))
    {
        const auto& tuning = schedule_[episodeIdx].tuning();
        initialTimeStepSize_ = tuning.TSINIT;
        maxTimeStepAfterWellEvent_ = tuning.TMAXWC;
        maxTimeStepSize_ = tuning.TSMAXZ;
        restartShrinkFactor_ = 1./tuning.TSFCNV;
        minTimeStepSize_ = tuning.TSMINZ;
        return true;
    }

    return false;
}


template<class GridView, class FluidSystem, class Scalar>
void EclGenericProblem<GridView,FluidSystem,Scalar>::
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

    // update explicit quantities between timesteps.
    const auto& oilVaporizationControl = schedule_[episodeIdx].oilvap();
    if (drsdtActive_(episodeIdx))
        // DRSDT is enabled
        for (size_t pvtRegionIdx = 0; pvtRegionIdx < maxDRs_.size(); ++pvtRegionIdx)
            maxDRs_[pvtRegionIdx] = oilVaporizationControl.getMaxDRSDT(pvtRegionIdx)*timeStepSize;

    if (drvdtActive_(episodeIdx))
        // DRVDT is enabled
        for (size_t pvtRegionIdx = 0; pvtRegionIdx < maxDRv_.size(); ++pvtRegionIdx)
            maxDRv_[pvtRegionIdx] = oilVaporizationControl.getMaxDRVDT(pvtRegionIdx)*timeStepSize;
}

template<class GridView, class FluidSystem, class Scalar>
void EclGenericProblem<GridView,FluidSystem,Scalar>::
checkDeckCompatibility_(const Deck& deck,
                        bool enableApiTracking,
                        bool enableSolvent,
                        bool enablePolymer,
                        bool enableExtbo,
                        bool enableEnergy,
                        int numPhases,
                        bool indicesGasEnabled,
                        bool indicesOilEnabled,
                        bool indicesWaterEnabled) const
{
    if (enableApiTracking)
        throw std::logic_error("API tracking is not yet implemented but requested at compile time.");
    if (!enableApiTracking && deck.hasKeyword("API"))
        throw std::logic_error("The simulator is build with API tracking disabled, but API tracking is requested by the deck.");

    if (enableSolvent && !deck.hasKeyword("SOLVENT"))
        throw std::runtime_error("The simulator requires the solvent option to be enabled, but the deck does not.");
    else if (!enableSolvent && deck.hasKeyword("SOLVENT"))
        throw std::runtime_error("The deck enables the solvent option, but the simulator is compiled without it.");

    if (enablePolymer && !deck.hasKeyword("POLYMER"))
        throw std::runtime_error("The simulator requires the polymer option to be enabled, but the deck does not.");
    else if (!enablePolymer && deck.hasKeyword("POLYMER"))
        throw std::runtime_error("The deck enables the polymer option, but the simulator is compiled without it.");

    if (enableExtbo && !deck.hasKeyword("PVTSOL"))
        throw std::runtime_error("The simulator requires the extendedBO option to be enabled, but the deck does not.");
    else if (!enableExtbo && deck.hasKeyword("PVTSOL"))
        throw std::runtime_error("The deck enables the extendedBO option, but the simulator is compiled without it.");

    if (deck.hasKeyword("TEMP") && deck.hasKeyword("THERMAL"))
        throw std::runtime_error("The deck enables both, the TEMP and the THERMAL options, but they are mutually exclusive.");

    bool deckEnergyEnabled = (deck.hasKeyword("TEMP") || deck.hasKeyword("THERMAL"));
    if (enableEnergy && !deckEnergyEnabled)
        throw std::runtime_error("The simulator requires the TEMP or the THERMAL option to be enabled, but the deck activates neither.");
    else if (!enableEnergy && deckEnergyEnabled)
        throw std::runtime_error("The deck enables the TEMP or the THERMAL option, but the simulator is not compiled to support either.");

    if (deckEnergyEnabled && deck.hasKeyword("TEMP"))
        std::cerr << "WARNING: The deck requests the TEMP option, i.e., treating energy "
                  << "conservation as a post processing step. This is currently unsupported, "
                  << "i.e., energy conservation is always handled fully implicitly." << std::endl;

    int numDeckPhases = FluidSystem::numActivePhases();
    if (numDeckPhases < numPhases)
        std::cerr << "WARNING: The number of active phases specified by the deck ("
                  << numDeckPhases << ") is smaller than the number of compiled-in phases ("
                  << numPhases << "). This usually results in a significant "
                  << "performance degradation compared to using a specialized simulator."  << std::endl;
    else if (numDeckPhases < numPhases)
        throw std::runtime_error("The deck enables "+std::to_string(numDeckPhases)+" phases "
                                 "while this simulator can only handle "+
                                 std::to_string(numPhases)+".");

    // make sure that the correct phases are active
    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && !indicesOilEnabled)
        throw std::runtime_error("The deck enables oil, but this simulator cannot handle it.");
    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && !indicesGasEnabled)
        throw std::runtime_error("The deck enables gas, but this simulator cannot handle it.");
    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) && !indicesWaterEnabled)
        throw std::runtime_error("The deck enables water, but this simulator cannot handle it.");
    // the opposite cases should be fine (albeit a bit slower than what's possible)
}

template<class GridView, class FluidSystem, class Scalar>
void EclGenericProblem<GridView,FluidSystem,Scalar>::
initFluidSystem_()
{
    FluidSystem::initFromState(eclState_, schedule_);
}

template<class GridView, class FluidSystem, class Scalar>
void EclGenericProblem<GridView,FluidSystem,Scalar>::
readBlackoilExtentionsInitialConditions_(size_t numDof,
                                         bool enableSolvent,
                                         bool enablePolymer,
                                         bool enablePolymerMolarWeight)
{
    if (enableSolvent) {
        if (eclState_.fieldProps().has_double("SSOL"))
            solventSaturation_ = eclState_.fieldProps().get_double("SSOL");
        else
            solventSaturation_.resize(numDof, 0.0);
    }

    if (enablePolymer) {
        if (eclState_.fieldProps().has_double("SPOLY"))
            polymerConcentration_ = eclState_.fieldProps().get_double("SPOLY");
        else
            polymerConcentration_.resize(numDof, 0.0);
    }

    if (enablePolymerMolarWeight) {
        if (eclState_.fieldProps().has_double("SPOLYMW"))
            polymerMoleWeight_ = eclState_.fieldProps().get_double("SPOLYMW");
        else
            polymerMoleWeight_.resize(numDof, 0.0);
    }
}


template<class GridView, class FluidSystem, class Scalar>
Scalar EclGenericProblem<GridView,FluidSystem,Scalar>::
maxWaterSaturation(unsigned globalDofIdx) const
{
    if (maxWaterSaturation_.empty())
        return 0.0;

    return maxWaterSaturation_[globalDofIdx];
}

template<class GridView, class FluidSystem, class Scalar>
Scalar EclGenericProblem<GridView,FluidSystem,Scalar>::
minOilPressure(unsigned globalDofIdx) const
{
    if (minOilPressure_.empty())
        return 0.0;

    return minOilPressure_[globalDofIdx];
}

template<class GridView, class FluidSystem, class Scalar>
Scalar EclGenericProblem<GridView,FluidSystem,Scalar>::
overburdenPressure(unsigned elementIdx) const
{
    if (overburdenPressure_.empty())
        return 0.0;

    return overburdenPressure_[elementIdx];
}

template<class GridView, class FluidSystem, class Scalar>
Scalar EclGenericProblem<GridView,FluidSystem,Scalar>::
solventSaturation(unsigned elemIdx) const
{
    if (solventSaturation_.empty())
        return 0;

    return solventSaturation_[elemIdx];
}

template<class GridView, class FluidSystem, class Scalar>
Scalar EclGenericProblem<GridView,FluidSystem,Scalar>::
polymerConcentration(unsigned elemIdx) const
{
    if (polymerConcentration_.empty())
        return 0;

    return polymerConcentration_[elemIdx];
}

template<class GridView, class FluidSystem, class Scalar>
Scalar EclGenericProblem<GridView,FluidSystem,Scalar>::
polymerMolecularWeight(const unsigned elemIdx) const
{
    if (polymerMoleWeight_.empty())
        return 0.0;

    return polymerMoleWeight_[elemIdx];
}

template<class GridView, class FluidSystem, class Scalar>
unsigned EclGenericProblem<GridView,FluidSystem,Scalar>::
pvtRegionIndex(unsigned elemIdx) const
{
    if (pvtnum_.empty())
        return 0;

    return pvtnum_[elemIdx];
}

template<class GridView, class FluidSystem, class Scalar>
unsigned EclGenericProblem<GridView,FluidSystem,Scalar>::
satnumRegionIndex(unsigned elemIdx) const
{
    if (satnum_.empty())
        return 0;

    return satnum_[elemIdx];
}

template<class GridView, class FluidSystem, class Scalar>
unsigned EclGenericProblem<GridView,FluidSystem,Scalar>::
miscnumRegionIndex(unsigned elemIdx) const
{
    if (miscnum_.empty())
        return 0;

    return miscnum_[elemIdx];
}

template<class GridView, class FluidSystem, class Scalar>
unsigned EclGenericProblem<GridView,FluidSystem,Scalar>::
plmixnumRegionIndex(unsigned elemIdx) const
{
    if (plmixnum_.empty())
        return 0;

    return plmixnum_[elemIdx];
}

template<class GridView, class FluidSystem, class Scalar>
Scalar EclGenericProblem<GridView,FluidSystem,Scalar>::
maxPolymerAdsorption(unsigned elemIdx) const
{
    if (maxPolymerAdsorption_.empty())
        return 0;

    return maxPolymerAdsorption_[elemIdx];
}

template<class GridView, class FluidSystem, class Scalar>
void EclGenericProblem<GridView,FluidSystem,Scalar>::
initDRSDT_(size_t numDof,
           int episodeIdx)
{
    // deal with DRSDT
    unsigned ntpvt = eclState_.runspec().tabdims().getNumPVTTables();
    //TODO We may want to only allocate these properties only if active.
    //But since they may be activated at later time we need some more
    //intrastructure to handle it
    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
        maxDRv_.resize(ntpvt, 1e30);
        lastRv_.resize(numDof, 0.0);
        maxDRs_.resize(ntpvt, 1e30);
        dRsDtOnlyFreeGas_.resize(ntpvt, false);
        lastRs_.resize(numDof, 0.0);
        maxDRv_.resize(ntpvt, 1e30);
        lastRv_.resize(numDof, 0.0);
        maxOilSaturation_.resize(numDof, 0.0);
        if (drsdtConvective_(episodeIdx)) {
            convectiveDrs_.resize(numDof, 1.0);
        }
    }
}

#if HAVE_DUNE_FEM
template class EclGenericProblem<Dune::GridView<Dune::Fem::GridPart2GridViewTraits<Dune::Fem::AdaptiveLeafGridPart<Dune::CpGrid, Dune::PartitionIteratorType(4), false>>>,
                                 BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,
                                 double>;
template class EclGenericProblem<Dune::GridView<Dune::Fem::GridPart2GridViewTraits<Dune::Fem::AdaptiveLeafGridPart<Dune::CpGrid, Dune::PartitionIteratorType(4), false>>>,
                                 BlackOilFluidSystem<double,EclAlternativeBlackOilIndexTraits>,
                                 double>;
#else
template class EclGenericProblem<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>,
                                 BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,
                                 double>;
template class EclGenericProblem<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>,
                                 BlackOilFluidSystem<double,EclAlternativeBlackOilIndexTraits>,
                                 double>;
#endif

template class EclGenericProblem<Dune::GridView<Dune::PolyhedralGridViewTraits<3,3,double,Dune::PartitionIteratorType(4)>>,
                                 BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,
                                 double>;

} // namespace Opm
