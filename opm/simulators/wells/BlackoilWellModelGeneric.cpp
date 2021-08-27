/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 - 2017 Statoil ASA.
  Copyright 2017 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2016 - 2018 IRIS AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <config.h>
#include <opm/simulators/wells/BlackoilWellModelGeneric.hpp>

#include <opm/output/data/Groups.hpp>
#include <opm/output/eclipse/RestartValue.hpp>

#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group/GuideRate.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>

#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/GasLiftStage2.hpp>
#include <opm/simulators/wells/VFPProperties.hpp>
#include <opm/simulators/wells/WellGroupHelpers.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <cassert>
#include <stdexcept>

#include <fmt/format.h>

namespace Opm {

BlackoilWellModelGeneric::
BlackoilWellModelGeneric(Schedule& schedule,
                         const SummaryState& summaryState,
                         const EclipseState& eclState,
                         const PhaseUsage& phase_usage,
                         const Comm& comm)
    : schedule_(schedule)
    , summaryState_(summaryState)
    , eclState_(eclState)
    , comm_(comm)
    , phase_usage_(phase_usage)
    , guideRate_(schedule)
    , active_wgstate_(phase_usage)
    , last_valid_wgstate_(phase_usage)
    , nupcol_wgstate_(phase_usage)
{

  const auto numProcs = comm_.size();
  this->not_on_process_ = [this, numProcs](const Well& well) {
      if (numProcs == decltype(numProcs){1})
          return false;

      // Recall: false indicates NOT active!
      const auto value = std::make_pair(well.name(), true);
      auto candidate = std::lower_bound(this->parallel_well_info_.begin(),
                                        this->parallel_well_info_.end(),
                                        value);

      return (candidate == this->parallel_well_info_.end())
              || (*candidate != value);
  };
}

int
BlackoilWellModelGeneric::
numLocalWells() const
{
    return wells_ecl_.size();
}

int
BlackoilWellModelGeneric::
numPhases() const
{
    return phase_usage_.num_phases;
}

bool
BlackoilWellModelGeneric::
hasWell(const std::string& wname)
{
    auto iter = std::find_if(this->wells_ecl_.begin(), this->wells_ecl_.end(),
                             [&wname](const Well& well) { return well.name() == wname; });
    return (iter != this->wells_ecl_.end());
}

bool
BlackoilWellModelGeneric::
wellsActive() const
{
    return wells_active_;
}

bool
BlackoilWellModelGeneric::
localWellsActive() const
{
    return numLocalWells() > 0;
}

bool
BlackoilWellModelGeneric::
anyMSWellOpenLocal() const
{
    for (const auto& well : wells_ecl_) {
        if (well.isMultiSegment()) {
            return true;
        }
    }
    return false;
}

const Well&
BlackoilWellModelGeneric::
getWellEcl(const std::string& well_name) const
{
    // finding the iterator of the well in wells_ecl
    auto well_ecl = std::find_if(wells_ecl_.begin(),
                                 wells_ecl_.end(),
                                 [&well_name](const Well& elem)->bool {
                                     return elem.name() == well_name;
                                 });

    assert(well_ecl != wells_ecl_.end());

    return *well_ecl;
}

void
BlackoilWellModelGeneric::
loadRestartData(const data::Wells& rst_wells,
                const data::GroupAndNetworkValues& grpNwrkValues,
                const PhaseUsage& phases,
                const bool handle_ms_well,
                WellState& well_state)
{
    using GPMode = Group::ProductionCMode;
    using GIMode = Group::InjectionCMode;

    using rt = data::Rates::opt;
    const auto np = phases.num_phases;

    std::vector< rt > phs( np );
    if( phases.phase_used[BlackoilPhases::Aqua] ) {
        phs.at( phases.phase_pos[BlackoilPhases::Aqua] ) = rt::wat;
    }

    if( phases.phase_used[BlackoilPhases::Liquid] ) {
        phs.at( phases.phase_pos[BlackoilPhases::Liquid] ) = rt::oil;
    }

    if( phases.phase_used[BlackoilPhases::Vapour] ) {
        phs.at( phases.phase_pos[BlackoilPhases::Vapour] ) = rt::gas;
    }

    for( std::size_t well_index = 0; well_index < well_state.size(); well_index++) {
        const auto& well_name = well_state.name(well_index);
        const auto& rst_well = rst_wells.at(well_name);
        auto& ws = well_state.well(well_index);
        ws.bhp = rst_well.bhp;
        ws.thp = rst_well.thp;
        ws.temperature = rst_well.temperature;

        if (rst_well.current_control.isProducer) {
            ws.production_cmode = rst_well.current_control.prod;
        }
        else {
            ws.injection_cmode = rst_well.current_control.inj;
        }

        for( size_t i = 0; i < phs.size(); ++i ) {
            assert( rst_well.rates.has( phs[ i ] ) );
            ws.surface_rates[i] = rst_well.rates.get(phs[i]);
        }

        auto& perf_data = ws.perf_data;
        auto& perf_pressure = perf_data.pressure;
        auto& perf_rates = perf_data.rates;
        auto& perf_phase_rates = perf_data.phase_rates;
        const auto& old_perf_data = this->well_perf_data_[well_index];

        for (std::size_t perf_index = 0; perf_index < old_perf_data.size(); perf_index++) {
            const auto& pd = old_perf_data[perf_index];
            const auto& rst_connection = rst_well.connections[pd.ecl_index];
            perf_pressure[perf_index] = rst_connection.pressure;
            perf_rates[perf_index] = rst_connection.reservoir_rate;
            for (int phase_index = 0; phase_index < np; ++phase_index)
                perf_phase_rates[perf_index*np + phase_index] = rst_connection.rates.get(phs[phase_index]);
        }

        if (handle_ms_well && !rst_well.segments.empty()) {
            // we need the well_ecl_ information
            const Well& well_ecl = getWellEcl(well_name);

            const WellSegments& segment_set = well_ecl.getSegments();

            const auto& rst_segments = rst_well.segments;

            // \Note: eventually we need to handle the situations that some segments are shut
            assert(0u + segment_set.size() == rst_segments.size());

            auto& segments = ws.segments;
            auto& segment_pressure = segments.pressure;
            auto& segment_rates  = segments.rates;
            for (const auto& rst_segment : rst_segments) {
                  const int segment_index = segment_set.segmentNumberToIndex(rst_segment.first);

                  // recovering segment rates and pressure from the restart values
                  const auto pres_idx = data::SegmentPressures::Value::Pressure;
                  segment_pressure[segment_index] = rst_segment.second.pressures[pres_idx];

                  const auto& rst_segment_rates = rst_segment.second.rates;
                  for (int p = 0; p < np; ++p) {
                      segment_rates[segment_index * np + p] = rst_segment_rates.get(phs[p]);
                  }
              }
        }
    }

    for (const auto& [group, value] : grpNwrkValues.groupData) {
        const auto cpc = value.currentControl.currentProdConstraint;
        const auto cgi = value.currentControl.currentGasInjectionConstraint;
        const auto cwi = value.currentControl.currentWaterInjectionConstraint;

        if (cpc != GPMode::NONE) {
            this->groupState().production_control(group, cpc);
        }

        if (cgi != GIMode::NONE) {
            this->groupState().injection_control(group, Phase::GAS, cgi);
        }

        if (cwi != GIMode::NONE) {
            this->groupState().injection_control(group, Phase::WATER, cwi);
        }
    }
}

void
BlackoilWellModelGeneric::
initFromRestartFile(const RestartValue& restartValues,
                    const size_t numCells,
                    bool handle_ms_well)
{
    // The restart step value is used to identify wells present at the given
    // time step. Wells that are added at the same time step as RESTART is initiated
    // will not be present in a restart file. Use the previous time step to retrieve
    // wells that have information written to the restart file.
    const int report_step = std::max(eclState_.getInitConfig().getRestartStep() - 1, 0);
    // wells_ecl_ should only contain wells on this processor.
    wells_ecl_ = getLocalWells(report_step);
    local_parallel_well_info_ = createLocalParallelWellInfo(wells_ecl_);

    this->initializeWellProdIndCalculators();
    initializeWellPerfData();

    const int nw = wells_ecl_.size();
    if (nw > 0) {
        handle_ms_well &= anyMSWellOpenLocal();
        this->wellState().resize(wells_ecl_, local_parallel_well_info_, schedule(), handle_ms_well, numCells, well_perf_data_, summaryState_); // Resize for restart step
        loadRestartData(restartValues.wells, restartValues.grp_nwrk, phase_usage_, handle_ms_well, this->wellState());
    }

    this->commitWGState();
    initial_step_ = false;
}

void
BlackoilWellModelGeneric::
setWellsActive(const bool wells_active)
{
    wells_active_ = wells_active;
}

std::vector<Well>
BlackoilWellModelGeneric::
getLocalWells(const int timeStepIdx) const
{
    auto w = schedule().getWells(timeStepIdx);
    w.erase(std::remove_if(w.begin(), w.end(), not_on_process_), w.end());
    return w;
}

std::vector<ParallelWellInfo*>
BlackoilWellModelGeneric::
createLocalParallelWellInfo(const std::vector<Well>& wells)
{
    std::vector<ParallelWellInfo*> local_parallel_well_info;
    local_parallel_well_info.reserve(wells.size());
    for (const auto& well : wells)
    {
        auto wellPair = std::make_pair(well.name(), true);
        auto pwell = std::lower_bound(parallel_well_info_.begin(),
                                      parallel_well_info_.end(),
                                      wellPair);
        assert(pwell != parallel_well_info_.end() &&
               *pwell == wellPair);
        local_parallel_well_info.push_back(&(*pwell));
    }
    return local_parallel_well_info;
}

void
BlackoilWellModelGeneric::
initializeWellProdIndCalculators()
{
    this->prod_index_calc_.clear();
    this->prod_index_calc_.reserve(this->wells_ecl_.size());
    for (const auto& well : this->wells_ecl_) {
        this->prod_index_calc_.emplace_back(well);
    }
}

void
BlackoilWellModelGeneric::
initializeWellPerfData()
{
    well_perf_data_.resize(wells_ecl_.size());
    int well_index = 0;
    for (const auto& well : wells_ecl_) {
        int completion_index = 0;
        // INVALID_ECL_INDEX marks no above perf available
        int completion_index_above = ParallelWellInfo::INVALID_ECL_INDEX;
        well_perf_data_[well_index].clear();
        well_perf_data_[well_index].reserve(well.getConnections().size());
        CheckDistributedWellConnections checker(well, *local_parallel_well_info_[well_index]);
        bool hasFirstPerforation = false;
        bool firstOpenCompletion = true;
        auto& parallelWellInfo = *local_parallel_well_info_[well_index];
        parallelWellInfo.beginReset();

        for (const auto& completion : well.getConnections()) {
            const int active_index =
                cartesian_to_compressed_[completion.global_index()];
            if (completion.state() == Connection::State::OPEN) {
                if (active_index >= 0) {
                    if (firstOpenCompletion)
                    {
                        hasFirstPerforation = true;
                    }
                    checker.connectionFound(completion_index);
                    PerforationData pd;
                    pd.cell_index = active_index;
                    pd.connection_transmissibility_factor = completion.CF();
                    pd.satnum_id = completion.satTableId();
                    pd.ecl_index = completion_index;
                    well_perf_data_[well_index].push_back(pd);
                    parallelWellInfo.pushBackEclIndex(completion_index_above,
                                                      completion_index);
                }
                firstOpenCompletion = false;
                // Next time this index is the one above as each open completion is
                // is stored somehwere.
                completion_index_above = completion_index;
            } else {
                checker.connectionFound(completion_index);
                if (completion.state() != Connection::State::SHUT) {
                    OPM_THROW(std::runtime_error,
                              "Completion state: " << Connection::State2String(completion.state()) << " not handled");
                }
            }
            // Note: we rely on the connections being filtered! I.e. there are only connections
            // to active cells in the global grid.
            ++completion_index;
        }
        parallelWellInfo.endReset();
        checker.checkAllConnectionsFound();
        parallelWellInfo.communicateFirstPerforation(hasFirstPerforation);
        ++well_index;
    }
}

bool
BlackoilWellModelGeneric::
checkGroupConstraints(const Group& group,
                      const int reportStepIdx,
                      DeferredLogger& deferred_logger) const
{
    if (group.isInjectionGroup()) {
        const Phase all[] = {Phase::WATER, Phase::OIL, Phase::GAS};
        for (Phase phase : all) {
            if (!group.hasInjectionControl(phase)) {
                continue;
            }
            const auto& check = checkGroupInjectionConstraints(group, reportStepIdx, phase);
            if (check.first != Group::InjectionCMode::NONE) {
                return true;
            }
        }
    }
    if (group.isProductionGroup()) {
        const auto& check = checkGroupProductionConstraints(group, reportStepIdx, deferred_logger);
        if (check.first != Group::ProductionCMode::NONE)
        {
            return true;
        }
    }

    // call recursively down the group hiearchy
    bool violated = false;
    for (const std::string& groupName : group.groups()) {
        violated = violated || checkGroupConstraints( schedule().getGroup(groupName, reportStepIdx), reportStepIdx, deferred_logger);
    }
    return violated;
}

std::pair<Group::InjectionCMode, double>
BlackoilWellModelGeneric::
checkGroupInjectionConstraints(const Group& group,
                               const int reportStepIdx,
                               const Phase& phase) const
{
    const auto& well_state = this->wellState();

    int phasePos;
    if (phase == Phase::GAS && phase_usage_.phase_used[BlackoilPhases::Vapour] )
        phasePos = phase_usage_.phase_pos[BlackoilPhases::Vapour];
    else if (phase == Phase::OIL && phase_usage_.phase_used[BlackoilPhases::Liquid])
        phasePos = phase_usage_.phase_pos[BlackoilPhases::Liquid];
    else if (phase == Phase::WATER && phase_usage_.phase_used[BlackoilPhases::Aqua] )
        phasePos = phase_usage_.phase_pos[BlackoilPhases::Aqua];
    else
        OPM_THROW(std::runtime_error, "Unknown phase" );

    auto currentControl = this->groupState().injection_control(group.name(), phase);
    if (group.has_control(phase, Group::InjectionCMode::RATE))
    {
        if (currentControl != Group::InjectionCMode::RATE)
        {
            double current_rate = 0.0;
            current_rate += WellGroupHelpers::sumWellSurfaceRates(group, schedule(), well_state, reportStepIdx, phasePos, /*isInjector*/true);

            // sum over all nodes
            current_rate = comm_.sum(current_rate);

            const auto& controls = group.injectionControls(phase, this->summaryState_);
            double target = controls.surface_max_rate;

            if (group.has_gpmaint_control(phase, Group::InjectionCMode::RATE) && this->groupState().has_gpmaint_target(group.name()))
                target = this->groupState().gpmaint_target(group.name());

            if (target < current_rate) {
                double scale = 1.0;
                if (current_rate > 1e-12)
                    scale = target / current_rate;
                return std::make_pair(Group::InjectionCMode::RATE, scale);
            }
        }
    }
    if (group.has_control(phase, Group::InjectionCMode::RESV))
    {
        if (currentControl != Group::InjectionCMode::RESV)
        {
            double current_rate = 0.0;
            current_rate += WellGroupHelpers::sumWellResRates(group, schedule(), well_state, reportStepIdx, phasePos, /*isInjector*/true);
            // sum over all nodes
            current_rate = comm_.sum(current_rate);

            const auto& controls = group.injectionControls(phase, this->summaryState_);
            double target = controls.resv_max_rate;

            if (group.has_gpmaint_control(phase, Group::InjectionCMode::RESV) && this->groupState().has_gpmaint_target(group.name()))
                target = this->groupState().gpmaint_target(group.name());

            if (target < current_rate) {
                double scale = 1.0;
                if (current_rate > 1e-12)
                    scale = target / current_rate;
                return std::make_pair(Group::InjectionCMode::RESV, scale);
            }
        }
    }
    if (group.has_control(phase, Group::InjectionCMode::REIN))
    {
        if (currentControl != Group::InjectionCMode::REIN)
        {
            double production_Rate = 0.0;
            const auto& controls = group.injectionControls(phase, this->summaryState_);
            const Group& groupRein = schedule().getGroup(controls.reinj_group, reportStepIdx);
            production_Rate += WellGroupHelpers::sumWellSurfaceRates(groupRein, schedule(), well_state, reportStepIdx, phasePos, /*isInjector*/false);

            // sum over all nodes
            production_Rate = comm_.sum(production_Rate);

            double current_rate = 0.0;
            current_rate += WellGroupHelpers::sumWellSurfaceRates(group, schedule(), well_state, reportStepIdx, phasePos, /*isInjector*/true);

            // sum over all nodes
            current_rate = comm_.sum(current_rate);

            if (controls.target_reinj_fraction*production_Rate < current_rate) {
                double scale = 1.0;
                if (current_rate > 1e-12)
                    scale = controls.target_reinj_fraction*production_Rate / current_rate;
                return std::make_pair(Group::InjectionCMode::REIN, scale);
            }
        }
    }
    if (group.has_control(phase, Group::InjectionCMode::VREP))
    {
        if (currentControl != Group::InjectionCMode::VREP)
        {
            double voidage_rate = 0.0;
            const auto& controls = group.injectionControls(phase, this->summaryState_);
            const Group& groupVoidage = schedule().getGroup(controls.voidage_group, reportStepIdx);
            voidage_rate += WellGroupHelpers::sumWellResRates(groupVoidage, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Aqua], false);
            voidage_rate += WellGroupHelpers::sumWellResRates(groupVoidage, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Liquid], false);
            voidage_rate += WellGroupHelpers::sumWellResRates(groupVoidage, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Vapour], false);

            // sum over all nodes
            voidage_rate = comm_.sum(voidage_rate);

            double total_rate = 0.0;
            total_rate += WellGroupHelpers::sumWellResRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Aqua], true);
            total_rate += WellGroupHelpers::sumWellResRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Liquid], true);
            total_rate += WellGroupHelpers::sumWellResRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Vapour], true);

            // sum over all nodes
            total_rate = comm_.sum(total_rate);

            if (controls.target_void_fraction*voidage_rate < total_rate) {
                double scale = 1.0;
                if (total_rate > 1e-12)
                    scale = controls.target_void_fraction*voidage_rate / total_rate;
                return std::make_pair(Group::InjectionCMode::VREP, scale);
            }
        }
    }
    return std::make_pair(Group::InjectionCMode::NONE, 1.0);
}

std::pair<Group::ProductionCMode, double>
BlackoilWellModelGeneric::
checkGroupProductionConstraints(const Group& group,
                                const int reportStepIdx,
                                DeferredLogger& deferred_logger) const
{
    const auto& well_state = this->wellState();

    const auto controls = group.productionControls(summaryState_);
    const Group::ProductionCMode& currentControl = this->groupState().production_control(group.name());

    if (group.has_control(Group::ProductionCMode::ORAT))
    {
        if (currentControl != Group::ProductionCMode::ORAT)
        {
            double current_rate = 0.0;
            current_rate += WellGroupHelpers::sumWellSurfaceRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Liquid], false);

            // sum over all nodes
            current_rate = comm_.sum(current_rate);

            if (controls.oil_target < current_rate  ) {
                double scale = 1.0;
                if (current_rate > 1e-12)
                    scale = controls.oil_target / current_rate;
                return std::make_pair(Group::ProductionCMode::ORAT, scale);
            }
        }
    }

    if (group.has_control(Group::ProductionCMode::WRAT))
    {
        if (currentControl != Group::ProductionCMode::WRAT)
        {

            double current_rate = 0.0;
            current_rate += WellGroupHelpers::sumWellSurfaceRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Aqua], false);

            // sum over all nodes
            current_rate = comm_.sum(current_rate);

            if (controls.water_target < current_rate  ) {
                double scale = 1.0;
                if (current_rate > 1e-12)
                    scale = controls.water_target / current_rate;
                return std::make_pair(Group::ProductionCMode::WRAT, scale);
            }
        }
    }
    if (group.has_control(Group::ProductionCMode::GRAT))
    {
        if (currentControl != Group::ProductionCMode::GRAT)
        {
            double current_rate = 0.0;
            current_rate += WellGroupHelpers::sumWellSurfaceRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Vapour], false);

            // sum over all nodes
            current_rate = comm_.sum(current_rate);
            if (controls.gas_target < current_rate  ) {
                double scale = 1.0;
                if (current_rate > 1e-12)
                    scale = controls.gas_target / current_rate;
                return std::make_pair(Group::ProductionCMode::GRAT, scale);
            }
        }
    }
    if (group.has_control(Group::ProductionCMode::LRAT))
    {
        if (currentControl != Group::ProductionCMode::LRAT)
        {
            double current_rate = 0.0;
            current_rate += WellGroupHelpers::sumWellSurfaceRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Liquid], false);
            current_rate += WellGroupHelpers::sumWellSurfaceRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Aqua], false);

            // sum over all nodes
            current_rate = comm_.sum(current_rate);

            if (controls.liquid_target < current_rate  ) {
                double scale = 1.0;
                if (current_rate > 1e-12)
                    scale = controls.liquid_target / current_rate;
                return std::make_pair(Group::ProductionCMode::LRAT, scale);
            }
        }
    }

    if (group.has_control(Group::ProductionCMode::CRAT))
    {
        OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "CRAT control for production groups not implemented" , deferred_logger);
    }
    if (group.has_control(Group::ProductionCMode::RESV))
    {
        if (currentControl != Group::ProductionCMode::RESV)
        {
            double current_rate = 0.0;
            current_rate += WellGroupHelpers::sumWellResRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Aqua], true);
            current_rate += WellGroupHelpers::sumWellResRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Liquid], true);
            current_rate += WellGroupHelpers::sumWellResRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Vapour], true);

            // sum over all nodes
            current_rate = comm_.sum(current_rate);

            double target = controls.resv_target;
            if (group.has_gpmaint_control(Group::ProductionCMode::RESV) && this->groupState().has_gpmaint_target(group.name()))
                target = this->groupState().gpmaint_target(group.name());

            if ( target < current_rate ) {
                double scale = 1.0;
                if (current_rate > 1e-12)
                    scale = target / current_rate;
                return std::make_pair(Group::ProductionCMode::RESV, scale);
            }
        }
    }
    if (group.has_control(Group::ProductionCMode::PRBL))
    {
        OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "PRBL control for production groups not implemented", deferred_logger);
    }
    return std::make_pair(Group::ProductionCMode::NONE, 1.0);
}

void
BlackoilWellModelGeneric::
checkGconsaleLimits(const Group& group,
                    WellState& well_state,
                    const int reportStepIdx,
                    DeferredLogger& deferred_logger)
{
     // call recursively down the group hiearchy
    for (const std::string& groupName : group.groups()) {
        checkGconsaleLimits( schedule().getGroup(groupName, reportStepIdx), well_state, reportStepIdx, deferred_logger);
    }

    // only for groups with gas injection controls
    if (!group.hasInjectionControl(Phase::GAS)) {
        return;
    }

    // check if gconsale is used for this group
    if (!schedule()[reportStepIdx].gconsale().has(group.name()))
        return;

    std::ostringstream ss;

    const auto& gconsale = schedule()[reportStepIdx].gconsale().get(group.name(), summaryState_);
    const Group::ProductionCMode& oldProductionControl = this->groupState().production_control(group.name());

    int gasPos = phase_usage_.phase_pos[BlackoilPhases::Vapour];
    double production_rate = WellGroupHelpers::sumWellSurfaceRates(group, schedule(), well_state, reportStepIdx, gasPos, /*isInjector*/false);
    double injection_rate = WellGroupHelpers::sumWellSurfaceRates(group, schedule(), well_state, reportStepIdx, gasPos, /*isInjector*/true);

    // sum over all nodes
    injection_rate = comm_.sum(injection_rate);
    production_rate = comm_.sum(production_rate);

    double sales_rate = production_rate - injection_rate;
    double production_target = gconsale.sales_target + injection_rate;

    // add import rate and substract consumption rate for group for gas
    if (schedule()[reportStepIdx].gconsump().has(group.name())) {
        const auto& gconsump = schedule()[reportStepIdx].gconsump().get(group.name(), summaryState_);
        if (phase_usage_.phase_used[BlackoilPhases::Vapour]) {
            sales_rate += gconsump.import_rate;
            sales_rate -= gconsump.consumption_rate;
            production_target -= gconsump.import_rate;
            production_target += gconsump.consumption_rate;
        }
    }

    if (sales_rate > gconsale.max_sales_rate) {
        switch(gconsale.max_proc) {
        case GConSale::MaxProcedure::NONE: {
            if (oldProductionControl != Group::ProductionCMode::GRAT && oldProductionControl != Group::ProductionCMode::NONE) {
                ss << "Group sales exceed maximum limit, but the action is NONE for " + group.name() + ". Nothing happens";
            }
            break;
            }
        case GConSale::MaxProcedure::CON: {
            OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "GCONSALE exceed limit CON not implemented", deferred_logger);
            break;
        }
        case GConSale::MaxProcedure::CON_P: {
            OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "GCONSALE exceed limit CON_P not implemented", deferred_logger);
            break;
        }
        case GConSale::MaxProcedure::WELL: {
            OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "GCONSALE exceed limit WELL not implemented", deferred_logger);
            break;
        }
        case GConSale::MaxProcedure::PLUG: {
            OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "GCONSALE exceed limit PLUG not implemented", deferred_logger);
            break;
        }
        case GConSale::MaxProcedure::MAXR: {
            OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "GCONSALE exceed limit MAXR not implemented", deferred_logger);
            break;
        }
        case GConSale::MaxProcedure::END: {
            OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "GCONSALE exceed limit END not implemented", deferred_logger);
            break;
        }
        case GConSale::MaxProcedure::RATE: {
            this->groupState().production_control(group.name(), Group::ProductionCMode::GRAT);
            ss << "Maximum GCONSALE limit violated for " << group.name() << ". The group is switched from ";
            ss << Group::ProductionCMode2String(oldProductionControl) << " to " << Group::ProductionCMode2String(Group::ProductionCMode::GRAT);
            ss << " and limited by the maximum sales rate after consumption and import are considered" ;
            this->groupState().update_grat_sales_target(group.name(), production_target);
            break;
        }
        default:
            throw("Invalid procedure for maximum rate limit selected for group" + group.name());
        }
    }
    if (sales_rate < gconsale.min_sales_rate) {
        const Group::ProductionCMode& currentProductionControl = this->groupState().production_control(group.name());
        if ( currentProductionControl == Group::ProductionCMode::GRAT ) {
            ss << "Group " + group.name() + " has sale rate less then minimum permitted value and is under GRAT control. \n";
            ss << "The GRAT is increased to meet the sales minimum rate. \n";
            this->groupState().update_grat_sales_target(group.name(), production_target);
        //} else if () {//TODO add action for WGASPROD
        //} else if () {//TODO add action for drilling queue
        } else {
            ss << "Group " + group.name() + " has sale rate less then minimum permitted value but cannot increase the group production rate \n";
            ss << "or adjust gas production using WGASPROD or drill new wells to meet the sales target. \n";
            ss << "Note that WGASPROD and drilling queues are not implemented in Flow. No action is taken. \n ";
        }
    }
    if (gconsale.sales_target < 0.0) {
        OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + " has sale rate target less then zero. Not implemented in Flow" , deferred_logger);
    }

    if (!ss.str().empty() && comm_.rank() == 0)
        deferred_logger.info(ss.str());
}

bool
BlackoilWellModelGeneric::
checkGroupHigherConstraints(const Group& group,
                            DeferredLogger& deferred_logger,
                            const int reportStepIdx,
                            std::set<std::string>& switched_groups)
{
    // Set up coefficients for RESV <-> surface rate conversion.
    // Use the pvtRegionIdx from the top cell of the first well.
    // TODO fix this!
    // This is only used for converting RESV rates.
    // What is the proper approach?
    const int fipnum = 0;
    int pvtreg = well_perf_data_.empty() || well_perf_data_[0].empty()
        ? pvt_region_idx_[0]
        : pvt_region_idx_[well_perf_data_[0][0].cell_index];

    bool changed = false;
    if ( comm_.size() > 1)
    {
        // Just like in the sequential case the pvtregion is determined
        // by the first cell of the first well. What is the first well
        // is decided by the order in the Schedule using Well::seqIndex()
        int firstWellIndex = well_perf_data_.empty() ?
            std::numeric_limits<int>::max() : wells_ecl_[0].seqIndex();
        auto regIndexPair = std::make_pair(pvtreg, firstWellIndex);
        std::vector<decltype(regIndexPair)> pairs(comm_.size());
        comm_.allgather(&regIndexPair, 1, pairs.data());
        pvtreg = std::min_element(pairs.begin(), pairs.end(),
                                  [](const auto& p1, const auto& p2){ return p1.second < p2.second;})
            ->first;
    }

    std::vector<double> rates(phase_usage_.num_phases, 0.0);

    const bool skip = switched_groups.count(group.name()) || group.name() == "FIELD";

    if (!skip && group.isInjectionGroup()) {
        // Obtain rates for group.
        std::vector<double> resv_coeff_inj(phase_usage_.num_phases, 0.0);
        calcInjRates(fipnum, pvtreg, resv_coeff_inj);

        for (int phasePos = 0; phasePos < phase_usage_.num_phases; ++phasePos) {
            const double local_current_rate = WellGroupHelpers::sumWellSurfaceRates(group, schedule(), this->wellState(), reportStepIdx, phasePos, /* isInjector */ true);
            // Sum over all processes
            rates[phasePos] = comm_.sum(local_current_rate);
        }
        const Phase all[] = { Phase::WATER, Phase::OIL, Phase::GAS };
        for (Phase phase : all) {
            // Check higher up only if under individual (not FLD) control.
            auto currentControl = this->groupState().injection_control(group.name(), phase);
            if (currentControl != Group::InjectionCMode::FLD && group.injectionGroupControlAvailable(phase)) {
                const Group& parentGroup = schedule().getGroup(group.parent(), reportStepIdx);
                const std::pair<bool, double> changed_this = WellGroupHelpers::checkGroupConstraintsInj(
                    group.name(),
                    group.parent(),
                    parentGroup,
                    this->wellState(),
                    this->groupState(),
                    reportStepIdx,
                    &guideRate_,
                    rates.data(),
                    phase,
                    phase_usage_,
                    group.getGroupEfficiencyFactor(),
                    schedule(),
                    summaryState_,
                    resv_coeff_inj,
                    deferred_logger);
                if (changed_this.first) {
                    switched_groups.insert(group.name());
                    actionOnBrokenConstraints(group, Group::InjectionCMode::FLD, phase, deferred_logger);
                    WellGroupHelpers::updateWellRatesFromGroupTargetScale(changed_this.second, group, schedule(), reportStepIdx, /* isInjector */ true, this->groupState(), this->wellState());
                    changed = true;
                }
            }
        }
    }

    if (!skip && group.isProductionGroup()) {
        // Obtain rates for group.
        for (int phasePos = 0; phasePos < phase_usage_.num_phases; ++phasePos) {
            const double local_current_rate = WellGroupHelpers::sumWellSurfaceRates(group, schedule(), this->wellState(), reportStepIdx, phasePos, /* isInjector */ false);
            // Sum over all processes
            rates[phasePos] = -comm_.sum(local_current_rate);
        }
        std::vector<double> resv_coeff(phase_usage_.num_phases, 0.0);
        calcRates(fipnum, pvtreg, resv_coeff);
        // Check higher up only if under individual (not FLD) control.
        const Group::ProductionCMode& currentControl = this->groupState().production_control(group.name());
            if (currentControl != Group::ProductionCMode::FLD && group.productionGroupControlAvailable()) {
                const Group& parentGroup = schedule().getGroup(group.parent(), reportStepIdx);
                const std::pair<bool, double> changed_this = WellGroupHelpers::checkGroupConstraintsProd(
                    group.name(),
                    group.parent(),
                    parentGroup,
                    this->wellState(),
                    this->groupState(),
                    reportStepIdx,
                    &guideRate_,
                    rates.data(),
                    phase_usage_,
                    group.getGroupEfficiencyFactor(),
                    schedule(),
                    summaryState_,
                    resv_coeff,
                    deferred_logger);
                if (changed_this.first) {
                    switched_groups.insert(group.name());
                    const auto exceed_action = group.productionControls(summaryState_).exceed_action;
                    actionOnBrokenConstraints(group, exceed_action, Group::ProductionCMode::FLD, deferred_logger);
                    WellGroupHelpers::updateWellRatesFromGroupTargetScale(changed_this.second, group, schedule(), reportStepIdx, /* isInjector */ false, this->groupState(), this->wellState());
                    changed = true;
                }
            }
    }

    // call recursively down the group hiearchy
    for (const std::string& groupName : group.groups()) {
        bool changed_this = checkGroupHigherConstraints( schedule().getGroup(groupName, reportStepIdx), deferred_logger, reportStepIdx, switched_groups);
        changed = changed || changed_this;
    }

    return changed;
}

bool
BlackoilWellModelGeneric::
updateGroupIndividualControl(const Group& group,
                             DeferredLogger& deferred_logger,
                             const int reportStepIdx,
                             std::set<std::string>& switched_groups)
{
    bool changed = false;
    const bool skip = switched_groups.count(group.name());
    if (!skip && group.isInjectionGroup())
    {
        const Phase all[] = {Phase::WATER, Phase::OIL, Phase::GAS};
        for (Phase phase : all) {
            if (!group.hasInjectionControl(phase)) {
                continue;
            }
            const auto& changed_this = checkGroupInjectionConstraints(group, reportStepIdx, phase);
            if (changed_this.first != Group::InjectionCMode::NONE)
            {
                switched_groups.insert(group.name());
                actionOnBrokenConstraints(group, changed_this.first, phase, deferred_logger);
                WellGroupHelpers::updateWellRatesFromGroupTargetScale(changed_this.second, group, schedule(), reportStepIdx, /* isInjector */ false, this->groupState(), this->wellState());
                changed = true;
            }
        }
    }
    if (!skip && group.isProductionGroup()) {
        const auto& changed_this = checkGroupProductionConstraints(group, reportStepIdx, deferred_logger);
        const auto controls = group.productionControls(summaryState_);
        if (changed_this.first != Group::ProductionCMode::NONE)
        {
            switched_groups.insert(group.name());
            actionOnBrokenConstraints(group, controls.exceed_action, changed_this.first, deferred_logger);
            WellGroupHelpers::updateWellRatesFromGroupTargetScale(changed_this.second, group, schedule(), reportStepIdx, /* isInjector */ false, this->groupState(), this->wellState());
            changed = true;
        }
    }

    // call recursively down the group hiearchy
    for (const std::string& groupName : group.groups()) {
        bool changed_this = updateGroupIndividualControl( schedule().getGroup(groupName, reportStepIdx), deferred_logger, reportStepIdx, switched_groups);
        changed = changed || changed_this;
    }

    return changed;
}

bool
BlackoilWellModelGeneric::
updateGroupIndividualControls(DeferredLogger& deferred_logger,
                              std::set<std::string>& switched_groups,
                              const int reportStepIdx,
                              const int iterationIdx)
{
    const int nupcol = schedule()[reportStepIdx].nupcol();
    // don't switch group control when iterationIdx > nupcol
    // to avoid oscilations between group controls
    if (iterationIdx > nupcol)
        return false;

    const Group& fieldGroup = schedule().getGroup("FIELD", reportStepIdx);
    return updateGroupIndividualControl(fieldGroup, deferred_logger,
                                        reportStepIdx, switched_groups);
}

bool
BlackoilWellModelGeneric::
updateGroupHigherControls(DeferredLogger& deferred_logger,
                          const int reportStepIdx,
                          std::set<std::string>& switched_groups)
{
    const Group& fieldGroup = schedule().getGroup("FIELD", reportStepIdx);
    return checkGroupHigherConstraints(fieldGroup, deferred_logger, reportStepIdx, switched_groups);
}

void
BlackoilWellModelGeneric::
actionOnBrokenConstraints(const Group& group,
                          const Group::ExceedAction& exceed_action,
                          const Group::ProductionCMode& newControl,
                          DeferredLogger& deferred_logger)
{
    const Group::ProductionCMode oldControl = this->groupState().production_control(group.name());

    std::ostringstream ss;

    switch(exceed_action) {
    case Group::ExceedAction::NONE: {
        if (oldControl != newControl && oldControl != Group::ProductionCMode::NONE) {
            ss << "Group production exceed action is NONE for group " + group.name() + ". Nothing happens.";
        }
        break;
    }
    case Group::ExceedAction::CON: {
        OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "GroupProductionExceedLimit CON not implemented", deferred_logger);
        break;
    }
    case Group::ExceedAction::CON_PLUS: {
        OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "GroupProductionExceedLimit CON_PLUS not implemented", deferred_logger);
        break;
    }
    case Group::ExceedAction::WELL: {
        OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "GroupProductionExceedLimit WELL not implemented", deferred_logger);
        break;
    }
    case Group::ExceedAction::PLUG: {
        OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "GroupProductionExceedLimit PLUG not implemented", deferred_logger);
        break;
    }
    case Group::ExceedAction::RATE: {
        if (oldControl != newControl) {
            this->groupState().production_control(group.name(), newControl);
            ss << "Switching production control mode for group "<< group.name()
               << " from " << Group::ProductionCMode2String(oldControl)
               << " to " << Group::ProductionCMode2String(newControl);
        }
        break;
    }
    default:
        throw("Invalid procedure for maximum rate limit selected for group" + group.name());
    }

    auto cc = Dune::MPIHelper::getCollectiveCommunication();
    if (!ss.str().empty() && cc.rank() == 0)
        deferred_logger.info(ss.str());
}

void
BlackoilWellModelGeneric::
actionOnBrokenConstraints(const Group& group,
                          const Group::InjectionCMode& newControl,
                          const Phase& controlPhase,
                          DeferredLogger& deferred_logger)
{
    auto oldControl = this->groupState().injection_control(group.name(), controlPhase);

    std::ostringstream ss;
    if (oldControl != newControl) {
        const std::string from = Group::InjectionCMode2String(oldControl);
        ss << "Switching injection control mode for group "<< group.name()
           << " from " << Group::InjectionCMode2String(oldControl)
           << " to " << Group::InjectionCMode2String(newControl);
        this->groupState().injection_control(group.name(), controlPhase, newControl);
    }
    auto cc = Dune::MPIHelper::getCollectiveCommunication();
    if (!ss.str().empty() && cc.rank() == 0)
        deferred_logger.info(ss.str());
}

void
BlackoilWellModelGeneric::
updateEclWells(const int timeStepIdx,
               const std::unordered_set<std::string>& wells)
{
    for (const auto& wname : wells) {
        auto well_iter = std::find_if( this->wells_ecl_.begin(), this->wells_ecl_.end(), [wname] (const auto& well) -> bool { return well.name() == wname;});
        if (well_iter != this->wells_ecl_.end()) {
            auto well_index = std::distance( this->wells_ecl_.begin(), well_iter );
            this->wells_ecl_[well_index] = schedule_.getWell(wname, timeStepIdx);

            const auto& well = this->wells_ecl_[well_index];
            auto& pd     = this->well_perf_data_[well_index];
            auto  pdIter = pd.begin();
            for (const auto& conn : well.getConnections()) {
                if (conn.state() != Connection::State::SHUT) {
                    pdIter->connection_transmissibility_factor = conn.CF();
                    ++pdIter;
                }
            }
            this->wellState().updateStatus(well_index, well.getStatus());
            this->wellState().resetConnectionTransFactors(well_index, pd);
            this->prod_index_calc_[well_index].reInit(well);
        }
    }
}

double
BlackoilWellModelGeneric::
wellPI(const int well_index) const
{
    const auto& pu = this->phase_usage_;
    const auto& pi = this->wellState().well(well_index).productivity_index;

    const auto preferred = this->wells_ecl_[well_index].getPreferredPhase();
    switch (preferred) { // Should really have LIQUID = OIL + WATER here too...
    case Phase::WATER:
        return pu.phase_used[BlackoilPhases::PhaseIndex::Aqua]
            ? pi[pu.phase_pos[BlackoilPhases::PhaseIndex::Aqua]]
            : 0.0;

    case Phase::OIL:
        return pu.phase_used[BlackoilPhases::PhaseIndex::Liquid]
            ? pi[pu.phase_pos[BlackoilPhases::PhaseIndex::Liquid]]
            : 0.0;

    case Phase::GAS:
        return pu.phase_used[BlackoilPhases::PhaseIndex::Vapour]
            ? pi[pu.phase_pos[BlackoilPhases::PhaseIndex::Vapour]]
            : 0.0;

    default:
        throw std::invalid_argument {
            "Unsupported preferred phase " +
            std::to_string(static_cast<int>(preferred))
        };
    }
}

double
BlackoilWellModelGeneric::
wellPI(const std::string& well_name) const
{
    auto well_iter = std::find_if(this->wells_ecl_.begin(), this->wells_ecl_.end(),
        [&well_name](const Well& well)
    {
        return well.name() == well_name;
    });

    if (well_iter == this->wells_ecl_.end()) {
        throw std::logic_error { "Could not find well: " + well_name };
    }

    auto well_index = std::distance(this->wells_ecl_.begin(), well_iter);
    return this->wellPI(well_index);
}

bool
BlackoilWellModelGeneric::
wasDynamicallyShutThisTimeStep(const int well_index) const
{
    return this->closed_this_step_.find(this->wells_ecl_[well_index].name()) !=
           this->closed_this_step_.end();
}

void
BlackoilWellModelGeneric::
updateWsolvent(const Group& group,
               const int reportStepIdx,
               const WellState& wellState)
{
    for (const std::string& groupName : group.groups()) {
        const Group& groupTmp = schedule_.getGroup(groupName, reportStepIdx);
        updateWsolvent(groupTmp, reportStepIdx, wellState);
    }

    if (group.isProductionGroup())
        return;

    auto currentGroupControl = this->groupState().injection_control(group.name(), Phase::GAS);
    if( currentGroupControl == Group::InjectionCMode::REIN ) {
        int gasPos = phase_usage_.phase_pos[BlackoilPhases::Vapour];
        const auto& controls = group.injectionControls(Phase::GAS, summaryState_);
        const Group& groupRein = schedule_.getGroup(controls.reinj_group, reportStepIdx);
        double gasProductionRate = WellGroupHelpers::sumWellSurfaceRates(groupRein, schedule_, wellState, reportStepIdx, gasPos, /*isInjector*/false);
        double solventProductionRate = WellGroupHelpers::sumSolventRates(groupRein, schedule_, wellState, reportStepIdx, /*isInjector*/false);

        solventProductionRate = comm_.sum(solventProductionRate);
        gasProductionRate = comm_.sum(gasProductionRate);

        double wsolvent = 0.0;
        if (std::abs(gasProductionRate) > 1e-6)
            wsolvent = solventProductionRate / gasProductionRate;

        setWsolvent(group, reportStepIdx, wsolvent);
    }
}

void
BlackoilWellModelGeneric::
setWsolvent(const Group& group,
            const int reportStepIdx,
            double wsolvent)
{
    for (const std::string& groupName : group.groups()) {
        const Group& groupTmp = schedule_.getGroup(groupName, reportStepIdx);
        setWsolvent(groupTmp, reportStepIdx, wsolvent);
    }

    for (const std::string& wellName : group.wells()) {
        const auto& wellTmp = schedule_.getWell(wellName, reportStepIdx);
        if (wellTmp.getStatus() == Well::Status::SHUT)
            continue;

        getGenWell(wellName)->setWsolvent(wsolvent);
    }
}

data::GuideRateValue
BlackoilWellModelGeneric::
getGuideRateValues(const Well& well) const
{
    auto grval = data::GuideRateValue{};

    const auto& wname = well.name();
    if (!this->wellState().has(wname)) {
        // No flow rates for 'wname' -- might be before well comes
        // online (e.g., for the initial condition before simulation
        // starts).
        return grval;
    }

    if (!this->guideRate_.has(wname)) {
        // No guiderates exist for 'wname'.
        return grval;
    }

    const auto qs = WellGroupHelpers::
        getWellRateVector(this->wellState(), this->phase_usage_, wname);

    this->getGuideRateValues(qs, well.isInjector(), wname, grval);

    return grval;
}

void
BlackoilWellModelGeneric::
getGuideRateValues(const GuideRate::RateVector& qs,
                   const bool                   is_inj,
                   const std::string&           wgname,
                   data::GuideRateValue&        grval) const
{
    auto getGR = [this, &wgname, &qs](const GuideRateModel::Target t)
    {
        return this->guideRate_.getSI(wgname, t, qs);
    };

    // Note: GuideRate does currently (2020-07-20) not support Target::RES.
    grval.set(data::GuideRateValue::Item::Gas,
              getGR(GuideRateModel::Target::GAS));

    grval.set(data::GuideRateValue::Item::Water,
              getGR(GuideRateModel::Target::WAT));

    if (!is_inj) {
        // Producer.  Extract "all" guiderate values.
        grval.set(data::GuideRateValue::Item::Oil,
                  getGR(GuideRateModel::Target::OIL));
    }
}

data::GuideRateValue
BlackoilWellModelGeneric::
getGuideRateValues(const Group& group) const
{
    auto grval = data::GuideRateValue{};
    const auto& gname = group.name();

    if (!this->groupState().has_production_rates(gname)) {
        // No flow rates for production group 'gname' -- might be before group comes
        // online (e.g., for the initial condition before simulation
        // starts).
        return grval;
    }

    if (!this->guideRate_.has(gname)) {
        // No guiderates exist for 'gname'.
        return grval;
    }

    const auto qs = WellGroupHelpers::getProductionGroupRateVector(this->groupState(), this->phase_usage_, gname);

    const auto is_inj = false; // This procedure only applies to G*PGR.
    this->getGuideRateValues(qs, is_inj, gname, grval);

    return grval;
}

data::GuideRateValue
BlackoilWellModelGeneric::
getGuideRateInjectionGroupValues(const Group& group) const
{
    auto grval = data::GuideRateValue{};

    const auto& gname = group.name();
    if (this->guideRate_.has(gname, Phase::GAS)) {
        grval.set(data::GuideRateValue::Item::Gas,
                  this->guideRate_.get(gname, Phase::GAS));
    }
    if (this->guideRate_.has(gname, Phase::WATER)) {
        grval.set(data::GuideRateValue::Item::Water,
                  this->guideRate_.get(gname, Phase::WATER));
    }
    return grval;
}

void
BlackoilWellModelGeneric::
assignWellGuideRates(data::Wells& wsrpt) const
{
    for (const auto& well : this->wells_ecl_) {
        auto xwPos = wsrpt.find(well.name());
        if (xwPos == wsrpt.end()) { // No well results.  Unexpected.
            continue;
        }

        xwPos->second.guide_rates = this->getGuideRateValues(well);
    }
}

void
BlackoilWellModelGeneric::
assignShutConnections(data::Wells& wsrpt,
                      const int reportStepIndex) const
{
    auto wellID = 0;

    for (const auto& well : this->wells_ecl_) {
        auto& xwel = wsrpt[well.name()]; // data::Wells is a std::map<>

        xwel.dynamicStatus = this->schedule()
                             .getWell(well.name(), reportStepIndex).getStatus();

        const auto wellIsOpen = xwel.dynamicStatus == Well::Status::OPEN;
        auto skip = [wellIsOpen](const Connection& conn)
        {
            return wellIsOpen && (conn.state() != Connection::State::SHUT);
        };

        if (this->wellTestState_.hasWellClosed(well.name()) &&
            !this->wasDynamicallyShutThisTimeStep(wellID))
        {
            xwel.dynamicStatus = well.getAutomaticShutIn()
                                 ? Well::Status::SHUT : Well::Status::STOP;
        }

        auto& xcon = xwel.connections;
        for (const auto& conn : well.getConnections()) {
            if (skip(conn)) {
                continue;
            }

            auto& xc = xcon.emplace_back();
            xc.index = conn.global_index();
            xc.pressure = xc.reservoir_rate = 0.0;

            xc.effective_Kh = conn.Kh();
            xc.trans_factor = conn.CF();
        }

        ++wellID;
    }
}

std::unordered_map<std::string, data::GroupGuideRates>
BlackoilWellModelGeneric::
calculateAllGroupGuiderates(const int reportStepIdx) const
{
    auto gr = std::unordered_map<std::string, data::GroupGuideRates>{};
    auto up = std::vector<std::string>{};

    // Start at well level, accumulate contributions towards root of
    // group tree (FIELD group).

    for (const auto& wname : schedule_.wellNames(reportStepIdx)) {
        if (! (this->wellState().has(wname) &&
               this->guideRate_.has(wname)))
        {
            continue;
        }

        const auto& well   = schedule_.getWell(wname, reportStepIdx);
        const auto& parent = well.groupName();

        if (parent == "FIELD") {
            // Well parented directly to "FIELD".  Inadvisable and
            // unexpected, but nothing to do about that here.  Just skip
            // this guide rate contribution.
            continue;
        }

        auto& grval = well.isInjector()
            ? gr[parent].injection
            : gr[parent].production;

        grval += this->getGuideRateValues(well);
        up.push_back(parent);
    }

    // Propagate accumulated guide rates up towards root of group tree.
    // Override accumulation if there is a GUIDERAT specification that
    // applies to a group.
    std::sort(up.begin(), up.end());
    auto start = 0*up.size();
    auto u     = std::unique(up.begin(), up.end());
    auto nu    = std::distance(up.begin(), u);
    while (nu > 0) {
        const auto ntot = up.size();

        for (auto gi = 0*nu; gi < nu; ++gi) {
            const auto& gname = up[start + gi];
            const auto& group = schedule_.getGroup(gname, reportStepIdx);

            if (this->guideRate_.has(gname)) {
                gr[gname].production = this->getGuideRateValues(group);
            }

            if (this->guideRate_.has(gname, Phase::WATER)
                    || this->guideRate_.has(gname, Phase::GAS)) {
                gr[gname].injection = this->getGuideRateInjectionGroupValues(group);
            }

            const auto parent = group.parent();
            if (parent == "FIELD") { continue; }

            gr[parent].injection  += gr[gname].injection;
            gr[parent].production += gr[gname].production;
            up.push_back(parent);
        }

        start = ntot;

        auto begin = up.begin() + ntot;
        std::sort(begin, up.end());
        u  = std::unique(begin, up.end());
        nu = std::distance(begin, u);
    }

    return gr;
}

void
BlackoilWellModelGeneric::
assignGroupControl(const Group& group,
                   data::GroupData& gdata) const
{
    const auto& gname     = group.name();
    const auto  grup_type = group.getGroupType();
    auto&       cgc       = gdata.currentControl;

    cgc.currentProdConstraint = Group::ProductionCMode::NONE;

    cgc.currentGasInjectionConstraint =
    cgc.currentWaterInjectionConstraint = Group::InjectionCMode::NONE;

    if (this->groupState().has_production_control(gname)) {
        cgc.currentProdConstraint = this->groupState().production_control(gname);
    }

    if ((grup_type == ::Opm::Group::GroupType::INJECTION) ||
        (grup_type == ::Opm::Group::GroupType::MIXED))
    {
        if (this->groupState().has_injection_control(gname, Phase::WATER)) {
            cgc.currentWaterInjectionConstraint = this->groupState().injection_control(gname, Phase::WATER);
        }

        if (this->groupState().has_injection_control(gname, Phase::GAS)) {
            cgc.currentGasInjectionConstraint = this->groupState().injection_control(gname, Phase::GAS);
        }
    }
}

void
BlackoilWellModelGeneric::
assignGroupGuideRates(const Group& group,
                      const std::unordered_map<std::string, data::GroupGuideRates>& groupGuideRates,
                      data::GroupData& gdata) const
{
    auto& prod = gdata.guideRates.production;  prod.clear();
    auto& inj  = gdata.guideRates.injection;   inj .clear();

    auto xgrPos = groupGuideRates.find(group.name());
    if ((xgrPos == groupGuideRates.end()) ||
        !this->guideRate_.has(group.name()))
    {
        // No guiderates defined for this group.
        return;
    }

    const auto& xgr = xgrPos->second;
    prod = xgr.production;
    inj  = xgr.injection;
}

void
BlackoilWellModelGeneric::
assignGroupValues(const int                               reportStepIdx,
                  std::map<std::string, data::GroupData>& gvalues) const
{
    const auto groupGuideRates =
        this->calculateAllGroupGuiderates(reportStepIdx);

    for (const auto& gname : schedule_.groupNames(reportStepIdx)) {
        const auto& grup = schedule_.getGroup(gname, reportStepIdx);

        auto& gdata = gvalues[gname];
        this->assignGroupControl(grup, gdata);
        this->assignGroupGuideRates(grup, groupGuideRates, gdata);
    }
}

void
BlackoilWellModelGeneric::
assignNodeValues(std::map<std::string, data::NodeData>& nodevalues) const
{
    nodevalues.clear();
    for (const auto& [node, pressure] : node_pressures_) {
        nodevalues.emplace(node, data::NodeData{pressure});
    }
}

data::GroupAndNetworkValues
BlackoilWellModelGeneric::
groupAndNetworkData(const int reportStepIdx) const
{
    auto grp_nwrk_values = data::GroupAndNetworkValues{};

    this->assignGroupValues(reportStepIdx, grp_nwrk_values.groupData);
    this->assignNodeValues(grp_nwrk_values.nodeData);

    return grp_nwrk_values;
}

void
BlackoilWellModelGeneric::
updateAndCommunicateGroupData(const int reportStepIdx,
                              const int iterationIdx)
{
    const Group& fieldGroup = schedule().getGroup("FIELD", reportStepIdx);
    const int nupcol = schedule()[reportStepIdx].nupcol();

    // This builds some necessary lookup structures, so it must be called
    // before we copy to well_state_nupcol_.
    this->wellState().updateGlobalIsGrup(comm_);

    if (iterationIdx < nupcol) {
        this->updateNupcolWGState();
    }

    auto& well_state = this->wellState();
    const auto& well_state_nupcol = this->nupcolWellState();
    // the group target reduction rates needs to be update since wells may have switched to/from GRUP control
    // The group target reduction does not honor NUPCOL.
    std::vector<double> groupTargetReduction(numPhases(), 0.0);
    WellGroupHelpers::updateGroupTargetReduction(fieldGroup, schedule(), reportStepIdx, /*isInjector*/ false, phase_usage_, guideRate_, well_state, this->groupState(), groupTargetReduction);
    std::vector<double> groupTargetReductionInj(numPhases(), 0.0);
    WellGroupHelpers::updateGroupTargetReduction(fieldGroup, schedule(), reportStepIdx, /*isInjector*/ true, phase_usage_, guideRate_, well_state, this->groupState(), groupTargetReductionInj);

    WellGroupHelpers::updateREINForGroups(fieldGroup, schedule(), reportStepIdx, phase_usage_, summaryState_, well_state_nupcol, this->groupState());
    WellGroupHelpers::updateVREPForGroups(fieldGroup, schedule(), reportStepIdx, well_state_nupcol, this->groupState());

    WellGroupHelpers::updateReservoirRatesInjectionGroups(fieldGroup, schedule(), reportStepIdx, well_state_nupcol, this->groupState());
    WellGroupHelpers::updateSurfaceRatesInjectionGroups(fieldGroup, schedule(), reportStepIdx, well_state_nupcol, this->groupState());

    WellGroupHelpers::updateGroupProductionRates(fieldGroup, schedule(), reportStepIdx, well_state_nupcol, this->groupState());

    // We use the rates from the previous time-step to reduce oscillations
    WellGroupHelpers::updateWellRates(fieldGroup, schedule(), reportStepIdx, this->prevWellState(), well_state);

    // Set ALQ for off-process wells to zero
    for (const auto& wname : schedule().wellNames(reportStepIdx)) {
        const bool is_producer = schedule().getWell(wname, reportStepIdx).isProducer();
        const bool not_on_this_process = !well_state.has(wname);
        if (is_producer && not_on_this_process) {
            well_state.setALQ(wname, 0.0);
        }
    }

    well_state.communicateGroupRates(comm_);
    this->groupState().communicate_rates(comm_);
    // compute wsolvent fraction for REIN wells
    updateWsolvent(fieldGroup, reportStepIdx,  well_state_nupcol);
}

bool
BlackoilWellModelGeneric::
hasTHPConstraints() const
{
    int local_result = false;
    for (const auto& well : well_container_generic_) {
        if (well->wellHasTHPConstraints(summaryState_)) {
            local_result=true;
        }
    }
    return comm_.max(local_result);
}

bool
BlackoilWellModelGeneric::
forceShutWellByNameIfPredictionMode(const std::string& wellname,
                                    const double simulation_time)
{
    // Only add the well to the closed list on the
    // process that owns it.
    int well_was_shut = 0;
    for (const auto& well : well_container_generic_) {
        if (well->name() == wellname && !well->wellIsStopped()) {
            if (well->underPredictionMode()) {
                wellTestState_.closeWell(wellname, WellTestConfig::Reason::PHYSICAL, simulation_time);
                well_was_shut = 1;
            }
            break;
        }
    }

    // Communicate across processes if a well was shut.
    well_was_shut = comm_.max(well_was_shut);

    // Only log a message on the output rank.
    if (terminal_output_ && well_was_shut) {
        const std::string msg = "Well " + wellname
            + " will be shut because it cannot get converged.";
        OpmLog::info(msg);
    }

    return (well_was_shut == 1);
}

void
BlackoilWellModelGeneric::
inferLocalShutWells()
{
    this->local_shut_wells_.clear();

    const auto nw = this->numLocalWells();

    auto used = std::vector<bool>(nw, false);
    for (const auto& wellPtr : this->well_container_generic_) {
        used[wellPtr->indexOfWell()] = true;
    }

    for (auto wellID = 0; wellID < nw; ++wellID) {
        if (! used[wellID]) {
            this->local_shut_wells_.push_back(wellID);
        }
    }
}

void
BlackoilWellModelGeneric::
updateNetworkPressures(const int reportStepIdx)
{
    // Get the network and return if inactive.
    const auto& network = schedule()[reportStepIdx].network();
    if (!network.active()) {
        return;
    }
    node_pressures_ = WellGroupHelpers::computeNetworkPressures(network,
                                                                this->wellState(),
                                                                this->groupState(),
                                                                *(vfp_properties_->getProd()),
                                                                schedule(),
                                                                reportStepIdx);

    // Set the thp limits of wells
    for (auto& well : well_container_generic_) {
        // Producers only, since we so far only support the
        // "extended" network model (properties defined by
        // BRANPROP and NODEPROP) which only applies to producers.
        if (well->isProducer()) {
            const auto it = node_pressures_.find(well->wellEcl().groupName());
            if (it != node_pressures_.end()) {
                // The well belongs to a group with has a network pressure constraint,
                // set the dynamic THP constraint of the well accordingly.
                well->setDynamicThpLimit(it->second);
            }
        }
    }
}

void
BlackoilWellModelGeneric::
calculateEfficiencyFactors(const int reportStepIdx)
{
    if ( !localWellsActive() ) {
        return;
    }

    for (auto& well : well_container_generic_) {
        const Well& wellEcl = well->wellEcl();
        double well_efficiency_factor = wellEcl.getEfficiencyFactor();
        WellGroupHelpers::accumulateGroupEfficiencyFactor(schedule().getGroup(wellEcl.groupName(), reportStepIdx), schedule(), reportStepIdx, well_efficiency_factor);
        well->setWellEfficiencyFactor(well_efficiency_factor);
    }
}

WellInterfaceGeneric*
BlackoilWellModelGeneric::
getGenWell(const std::string& well_name)
{
    // finding the iterator of the well in wells_ecl
    auto well = std::find_if(well_container_generic_.begin(),
                             well_container_generic_.end(),
                                 [&well_name](const WellInterfaceGeneric* elem)->bool {
                                     return elem->name() == well_name;
                                 });

    assert(well != well_container_generic_.end());

    return *well;
}

void
BlackoilWellModelGeneric::
setRepRadiusPerfLength()
{
    for (const auto& well : well_container_generic_) {
        well->setRepRadiusPerfLength(cartesian_to_compressed_);
    }
}

void
BlackoilWellModelGeneric::
gliftDebug(const std::string& msg,
           DeferredLogger& deferred_logger) const
{
    if (this->glift_debug) {
        const std::string message = fmt::format(
            "  GLIFT (DEBUG) : BlackoilWellModel : {}", msg);
        deferred_logger.info(message);
    }
}

void
BlackoilWellModelGeneric::
gliftDebugShowALQ(DeferredLogger& deferred_logger)
{
    for (auto& well : this->well_container_generic_) {
        if (well->isProducer()) {
            auto alq = this->wellState().getALQ(well->name());
            const std::string msg = fmt::format("ALQ_REPORT : {} : {}",
                                                well->name(), alq);
            gliftDebug(msg, deferred_logger);
        }
    }
}

// If a group has any production rate constraints, and/or a limit
// on its total rate of lift gas supply,  allocate lift gas
// preferentially to the wells that gain the most benefit from
// it. Lift gas increments are allocated in turn to the well that
// currently has the largest weighted incremental gradient. The
// procedure takes account of any limits on the group production
// rate or lift gas supply applied to any level of group.
void
BlackoilWellModelGeneric::
gasLiftOptimizationStage2(DeferredLogger& deferred_logger,
                          GLiftProdWells& prod_wells,
                          GLiftOptWells& glift_wells,
                          GLiftWellStateMap& glift_well_state_map,
                          const int episodeIndex)
{
    GasLiftStage2 glift {episodeIndex,
                         comm_,
                         schedule_,
                         summaryState_,
                         deferred_logger,
                         this->wellState(),
                         prod_wells,
                         glift_wells,
                         glift_well_state_map};
    glift.runOptimize();
}

void
BlackoilWellModelGeneric::
updateWellPotentials(const int reportStepIdx,
                     const bool onlyAfterEvent,
                     const SummaryConfig& summaryConfig,
                     DeferredLogger& deferred_logger)
{
    auto well_state_copy = this->wellState();

    const bool write_restart_file = schedule().write_rst_file(reportStepIdx);
    auto exc_type = ExceptionType::NONE;
    std::string exc_msg;
    size_t widx = 0;
    for (const auto& well : well_container_generic_) {
        const bool needed_for_summary =
                ((summaryConfig.hasSummaryKey( "WWPI:" + well->name()) ||
                  summaryConfig.hasSummaryKey( "WOPI:" + well->name()) ||
                  summaryConfig.hasSummaryKey( "WGPI:" + well->name())) && well->isInjector()) ||
                ((summaryConfig.hasKeyword( "GWPI") ||
                  summaryConfig.hasKeyword( "GOPI") ||
                  summaryConfig.hasKeyword( "GGPI")) && well->isInjector()) ||
                ((summaryConfig.hasKeyword( "FWPI") ||
                  summaryConfig.hasKeyword( "FOPI") ||
                  summaryConfig.hasKeyword( "FGPI")) && well->isInjector()) ||
                ((summaryConfig.hasSummaryKey( "WWPP:" + well->name()) ||
                  summaryConfig.hasSummaryKey( "WOPP:" + well->name()) ||
                  summaryConfig.hasSummaryKey( "WGPP:" + well->name())) && well->isProducer()) ||
                ((summaryConfig.hasKeyword( "GWPP") ||
                  summaryConfig.hasKeyword( "GOPP") ||
                  summaryConfig.hasKeyword( "GGPP")) && well->isProducer()) ||
                ((summaryConfig.hasKeyword( "FWPP") ||
                  summaryConfig.hasKeyword( "FOPP") ||
                  summaryConfig.hasKeyword( "FGPP")) && well->isProducer());

        // At the moment, the following events are considered
        // for potentials update
        const uint64_t effective_events_mask = ScheduleEvents::WELL_STATUS_CHANGE
                                             + ScheduleEvents::COMPLETION_CHANGE
                                             + ScheduleEvents::WELL_PRODUCTIVITY_INDEX
                                             + ScheduleEvents::WELL_WELSPECS_UPDATE
                                             + ScheduleEvents::WELLGROUP_EFFICIENCY_UPDATE
                                             + ScheduleEvents::NEW_WELL
                                             + ScheduleEvents::PRODUCTION_UPDATE
                                             + ScheduleEvents::INJECTION_UPDATE;
        const auto& events = schedule()[reportStepIdx].wellgroup_events();
        const bool event = events.hasEvent(well->name(), ScheduleEvents::ACTIONX_WELL_EVENT) ||
                           (report_step_starts_ && events.hasEvent(well->name(), effective_events_mask));
        const bool needPotentialsForGuideRates = well->underPredictionMode() && (!onlyAfterEvent || event);
        const bool needPotentialsForOutput = !onlyAfterEvent && (needed_for_summary || write_restart_file);
        const bool compute_potential = needPotentialsForOutput || needPotentialsForGuideRates;
        if (compute_potential)
        {
            this->computePotentials(widx, well_state_copy, exc_msg, exc_type, deferred_logger);
        }
        ++widx;
    }
    logAndCheckForExceptionsAndThrow(deferred_logger, exc_type,
                                     "computeWellPotentials() failed: " + exc_msg,
                                     terminal_output_);

}

void
BlackoilWellModelGeneric::
runWellPIScaling(const int timeStepIdx,
                 DeferredLogger& local_deferredLogger)
{
    if (this->last_run_wellpi_.has_value() && (*this->last_run_wellpi_ == timeStepIdx)) {
        // We've already run WELPI scaling for this report step.  Most
        // common for the very first report step.  Don't redo WELPI scaling.
        return;
    }

    auto hasWellPIEvent = [this, timeStepIdx](const int well_index) -> bool
    {
        return this->schedule()[timeStepIdx].wellgroup_events()
            .hasEvent(this->wells_ecl_[well_index].name(),
                      ScheduleEvents::Events::WELL_PRODUCTIVITY_INDEX);
    };

    auto updateEclWell = [this, timeStepIdx](const int well_index) -> void
    {
        const auto& schedule = this->schedule();
        const auto& wname = this->wells_ecl_[well_index].name();
        this->wells_ecl_[well_index] = schedule.getWell(wname, timeStepIdx);

        const auto& well = this->wells_ecl_[well_index];
        auto& pd     = this->well_perf_data_[well_index];
        auto  pdIter = pd.begin();
        for (const auto& conn : well.getConnections()) {
            if (conn.state() != Connection::State::SHUT) {
                pdIter->connection_transmissibility_factor = conn.CF();
                ++pdIter;
            }
        }
        this->wellState().resetConnectionTransFactors(well_index, pd);
        this->prod_index_calc_[well_index].reInit(well);
    };


    auto rescaleWellPI =
        [this, timeStepIdx](const int    well_index,
                            const double newWellPI) -> void
    {
        const auto& wname = this->wells_ecl_[well_index].name();

        schedule_.applyWellProdIndexScaling(wname, timeStepIdx, newWellPI);
    };

    // Minimal well setup to compute PI/II values
    {
        auto saved_previous_wgstate = this->prevWGState();
        this->commitWGState();

        this->createWellContainer(timeStepIdx);
        this->inferLocalShutWells();

        this->initWellContainer();

        this->calculateProductivityIndexValues(local_deferredLogger);
        this->calculateProductivityIndexValuesShutWells(timeStepIdx, local_deferredLogger);

        this->commitWGState(std::move(saved_previous_wgstate));
    }

    const auto nw = this->numLocalWells();
    for (auto wellID = 0*nw; wellID < nw; ++wellID) {
        if (hasWellPIEvent(wellID)) {
            rescaleWellPI(wellID, this->wellPI(wellID));
            updateEclWell(wellID);
        }
    }

    this->last_run_wellpi_ = timeStepIdx;
}

}
