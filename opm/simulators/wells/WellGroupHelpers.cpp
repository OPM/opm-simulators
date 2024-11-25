/*
  Copyright 2019 Norce.
  Copyright 2020 Equinor ASA.

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
#include <opm/simulators/wells/WellGroupHelpers.hpp>

#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Group/GConSump.hpp>
#include <opm/input/eclipse/Schedule/Group/GConSale.hpp>
#include <opm/input/eclipse/Schedule/Group/GPMaint.hpp>
#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <opm/input/eclipse/Schedule/Network/ExtNetwork.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <opm/simulators/wells/FractionCalculator.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include <opm/simulators/wells/RegionAverageCalculator.hpp>
#include <opm/simulators/wells/TargetCalculator.hpp>
#include <opm/simulators/wells/VFPProdProperties.hpp>
#include <opm/simulators/wells/WellState.hpp>
#include <opm/simulators/wells/GroupState.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <set>
#include <stack>
#include <stdexcept>

namespace {
    template<class Scalar>
    Opm::GuideRate::RateVector
    getGuideRateVector(const std::vector<Scalar>& rates, const Opm::PhaseUsage& pu)
    {
        using Opm::BlackoilPhases;

        Scalar oilRate = 0.0;
        if (pu.phase_used[BlackoilPhases::Liquid]) {
            oilRate = rates[pu.phase_pos[BlackoilPhases::Liquid]];
        }

        Scalar gasRate = 0.0;
        if (pu.phase_used[BlackoilPhases::Vapour]) {
            gasRate = rates[pu.phase_pos[BlackoilPhases::Vapour]];
        }

        Scalar waterRate = 0.0;
        if (pu.phase_used[BlackoilPhases::Aqua]) {
            waterRate = rates[pu.phase_pos[BlackoilPhases::Aqua]];
        }

        return {oilRate, gasRate, waterRate};
    }

    template<class Scalar>
    Scalar sumWellPhaseRates(bool res_rates,
                             const Opm::Group& group,
                             const Opm::Schedule& schedule,
                             const Opm::WellState<Scalar>& wellState,
                             const int reportStepIdx,
                             const int phasePos,
                             const bool injector)
    {

        Scalar rate = 0.0;
        for (const std::string& groupName : group.groups()) {
            const auto& groupTmp = schedule.getGroup(groupName, reportStepIdx);
            const auto& gefac = groupTmp.getGroupEfficiencyFactor();
            rate += gefac * sumWellPhaseRates(res_rates, groupTmp, schedule, wellState, reportStepIdx, phasePos, injector);
        }

        for (const std::string& wellName : group.wells()) {
            const auto& well_index = wellState.index(wellName);
            if (!well_index.has_value())
                continue;

            if (! wellState.wellIsOwned(well_index.value(), wellName) ) // Only sum once
            {
                continue;
            }

            const auto& wellEcl = schedule.getWell(wellName, reportStepIdx);
            // only count producers or injectors
            if ((wellEcl.isProducer() && injector) || (wellEcl.isInjector() && !injector))
                continue;

            if (wellEcl.getStatus() == Opm::Well::Status::SHUT)
                continue;

            Scalar factor = wellEcl.getEfficiencyFactor();
            const auto& ws = wellState.well(well_index.value());
            if (res_rates) {
                const auto& well_rates = ws.reservoir_rates;
                if (injector)
                    rate += factor * well_rates[phasePos];
                else
                    rate -= factor * well_rates[phasePos];
            } else {
                const auto& well_rates = ws.surface_rates;
                if (injector)
                    rate += factor * well_rates[phasePos];
                else
                    rate -= factor * well_rates[phasePos];
            }
        }
        return rate;
    }
} // namespace Anonymous

namespace Opm {

template<class Scalar>
void WellGroupHelpers<Scalar>::
setCmodeGroup(const Group& group,
              const Schedule& schedule,
              const SummaryState& summaryState,
              const int reportStepIdx,
              GroupState<Scalar>& group_state)
{

    for (const std::string& groupName : group.groups()) {
        setCmodeGroup(schedule.getGroup(groupName, reportStepIdx),
                      schedule, summaryState, reportStepIdx, group_state);
    }

    // use NONE as default control
    const Phase all[] = {Phase::WATER, Phase::OIL, Phase::GAS};
    for (Phase phase : all) {
        if (!group_state.has_injection_control(group.name(), phase)) {
            group_state.injection_control(group.name(), phase, Group::InjectionCMode::NONE);
        }
    }
    if (!group_state.has_production_control(group.name())) {
        group_state.production_control(group.name(), Group::ProductionCMode::NONE);
    }

    const auto& events = schedule[reportStepIdx].wellgroup_events();
    if (group.isInjectionGroup()
        && events.hasEvent(group.name(), ScheduleEvents::GROUP_INJECTION_UPDATE)) {

        for (Phase phase : all) {
            if (!group.hasInjectionControl(phase))
                continue;

            const auto& controls = group.injectionControls(phase, summaryState);
            group_state.injection_control(group.name(), phase, controls.cmode);
        }
    }

    if (group.isProductionGroup()
        && events.hasEvent(group.name(), ScheduleEvents::GROUP_PRODUCTION_UPDATE)) {
        const auto controls = group.productionControls(summaryState);
        group_state.production_control(group.name(), controls.cmode);
    }

    if (group.has_gpmaint_control(Group::ProductionCMode::RESV)) {
        group_state.production_control(group.name(), Group::ProductionCMode::RESV);
    }
    for (Phase phase : all) {
        if (group.has_gpmaint_control(phase, Group::InjectionCMode::RATE)) {
            group_state.injection_control(group.name(), phase, Group::InjectionCMode::RATE);
        } else if (group.has_gpmaint_control(phase, Group::InjectionCMode::RESV)) {
            group_state.injection_control(group.name(), phase, Group::InjectionCMode::RESV);
        }
    }

    if (schedule[reportStepIdx].gconsale().has(group.name())) {
        group_state.injection_control(group.name(), Phase::GAS, Group::InjectionCMode::SALE);
    }
}

template<class Scalar>
void WellGroupHelpers<Scalar>::
accumulateGroupEfficiencyFactor(const Group& group,
                                const Schedule& schedule,
                                const int reportStepIdx,
                                Scalar& factor)
{
    factor *= group.getGroupEfficiencyFactor();
    if (group.parent() != "FIELD" && !group.parent().empty())
        accumulateGroupEfficiencyFactor(
            schedule.getGroup(group.parent(), reportStepIdx), schedule, reportStepIdx, factor);
}

template<class Scalar>
Scalar WellGroupHelpers<Scalar>::
sumWellSurfaceRates(const Group& group,
                    const Schedule& schedule,
                    const WellState<Scalar>& wellState,
                    const int reportStepIdx,
                    const int phasePos,
                    const bool injector)
{
    return sumWellPhaseRates(false, group, schedule, wellState, reportStepIdx, phasePos, injector);
}

template<class Scalar>
Scalar WellGroupHelpers<Scalar>::
sumWellResRates(const Group& group,
                const Schedule& schedule,
                const WellState<Scalar>& wellState,
                const int reportStepIdx,
                const int phasePos,
                const bool injector)
{
    return sumWellPhaseRates(true, group, schedule, wellState, reportStepIdx, phasePos, injector);
}

template<class Scalar>
Scalar WellGroupHelpers<Scalar>::
sumSolventRates(const Group& group,
                const Schedule& schedule,
                const WellState<Scalar>& wellState,
                const int reportStepIdx,
                const bool injector)
{
    Scalar rate = 0.0;
    for (const std::string& groupName : group.groups()) {
        const Group& groupTmp = schedule.getGroup(groupName, reportStepIdx);
        const auto& gefac = groupTmp.getGroupEfficiencyFactor();
        rate += gefac * sumSolventRates(groupTmp, schedule, wellState, reportStepIdx, injector);
    }

    for (const std::string& wellName : group.wells()) {
        const auto& well_index = wellState.index(wellName);
        if (!well_index.has_value())
            continue;

        if (! wellState.wellIsOwned(well_index.value(), wellName) ) // Only sum once
        {
            continue;
        }

        const auto& wellEcl = schedule.getWell(wellName, reportStepIdx);
        // only count producers or injectors
        if ((wellEcl.isProducer() && injector) || (wellEcl.isInjector() && !injector))
            continue;

        if (wellEcl.getStatus() == Well::Status::SHUT)
            continue;

        const auto& ws = wellState.well(well_index.value());
        const Scalar factor = wellEcl.getEfficiencyFactor();
        if (injector)
            rate += factor * ws.sum_solvent_rates();
        else
            rate -= factor * ws.sum_solvent_rates();
    }
    return rate;
}

template<class Scalar>
void WellGroupHelpers<Scalar>::
updateGuideRatesForInjectionGroups(const Group& group,
                                   const Schedule& schedule,
                                   const SummaryState& summaryState,
                                   const PhaseUsage& pu,
                                   const int reportStepIdx,
                                   const WellState<Scalar>& wellState,
                                   const GroupState<Scalar>& group_state,
                                   GuideRate* guideRate,
                                   DeferredLogger& deferred_logger)
{
    for (const std::string& groupName : group.groups()) {
        const Group& groupTmp = schedule.getGroup(groupName, reportStepIdx);
        updateGuideRatesForInjectionGroups(groupTmp, schedule, summaryState, pu, reportStepIdx, wellState, group_state, guideRate, deferred_logger);
    }
    const Phase all[] = {Phase::WATER, Phase::OIL, Phase::GAS};
    for (Phase phase : all) {

        if(!group.hasInjectionControl(phase))
            continue;

        Scalar guideRateValue = 0.0;
        const auto& controls = group.injectionControls(phase, summaryState);
        switch (controls.guide_rate_def){
        case Group::GuideRateInjTarget::RATE:
            break;
        case Group::GuideRateInjTarget::VOID:
        {
            guideRateValue = group_state.injection_vrep_rate(group.name());
            break;
        }
        case Group::GuideRateInjTarget::NETV:
        {
            guideRateValue = group_state.injection_vrep_rate(group.name());
            const std::vector<Scalar>& injRES = group_state.injection_reservoir_rates(group.name());
            if (phase != Phase::OIL && pu.phase_used[BlackoilPhases::Liquid])
                guideRateValue -= injRES[pu.phase_pos[BlackoilPhases::Liquid]];
            if (phase != Phase::GAS && pu.phase_used[BlackoilPhases::Vapour])
                guideRateValue -= injRES[pu.phase_pos[BlackoilPhases::Vapour]];
            if (phase != Phase::WATER && pu.phase_used[BlackoilPhases::Aqua])
                guideRateValue -= injRES[pu.phase_pos[BlackoilPhases::Aqua]];
            break;
        }
        case Group::GuideRateInjTarget::RESV:
            OPM_DEFLOG_THROW(std::runtime_error, "GUIDE PHASE RESV not implemented. Group " + group.name(), deferred_logger);
        case Group::GuideRateInjTarget::POTN:
            break;
        case Group::GuideRateInjTarget::NO_GUIDE_RATE:
            break;
        default:
            OPM_DEFLOG_THROW(std::logic_error,
                             "Invalid GuideRateInjTarget in updateGuideRatesForInjectionGroups",
                             deferred_logger);
        }

        const UnitSystem& unit_system = schedule.getUnits();
        guideRateValue = unit_system.from_si(UnitSystem::measure::rate, guideRateValue);
        guideRate->compute(group.name(), phase, reportStepIdx, guideRateValue);
    }
}

template<class Scalar>
void WellGroupHelpers<Scalar>::
updateGroupTargetReduction(const Group& group,
                           const Schedule& schedule,
                           const int reportStepIdx,
                           const bool isInjector,
                           const PhaseUsage& pu,
                           const GuideRate& guide_rate,
                           const WellState<Scalar>& wellState,
                           GroupState<Scalar>& group_state,
                           std::vector<Scalar>& groupTargetReduction)
{
    const int np = wellState.numPhases();
    for (const std::string& subGroupName : group.groups()) {
        std::vector<Scalar> subGroupTargetReduction(np, 0.0);
        const Group& subGroup = schedule.getGroup(subGroupName, reportStepIdx);
        updateGroupTargetReduction(subGroup,
                                   schedule,
                                   reportStepIdx,
                                   isInjector,
                                   pu,
                                   guide_rate,
                                   wellState,
                                   group_state,
                                   subGroupTargetReduction);

        const Scalar subGroupEfficiency = subGroup.getGroupEfficiencyFactor();

        // accumulate group contribution from sub group
        if (isInjector) {
            const Phase all[] = {Phase::WATER, Phase::OIL, Phase::GAS};
            for (Phase phase : all) {

                int phase_pos = -1;
                switch (phase) {
                case Phase::GAS:
                    if (pu.phase_used[BlackoilPhases::Vapour]) {
                        phase_pos = pu.phase_pos[BlackoilPhases::Vapour];
                    }
                    break;
                case Phase::OIL:
                    if (pu.phase_used[BlackoilPhases::Liquid]) {
                        phase_pos = pu.phase_pos[BlackoilPhases::Liquid];
                    }
                    break;
                case Phase::WATER:
                    if (pu.phase_used[BlackoilPhases::Aqua]) {
                        phase_pos = pu.phase_pos[BlackoilPhases::Aqua];
                    }
                    break;
                default:
                    // just to avoid warning
                    throw std::invalid_argument("unhandled phase enum");
                }

                // the phase is not present
                if (phase_pos == -1)
                    continue;

                const Group::InjectionCMode& currentGroupControl
                        = group_state.injection_control(subGroup.name(), phase);
                const bool individual_control = (currentGroupControl != Group::InjectionCMode::FLD
                        && currentGroupControl != Group::InjectionCMode::NONE);
                const int num_group_controlled_wells
                        = groupControlledWells(schedule, wellState, group_state, reportStepIdx, subGroupName, "", !isInjector, phase);
                if (individual_control || num_group_controlled_wells == 0) {
                    groupTargetReduction[phase_pos]
                        += subGroupEfficiency * sumWellSurfaceRates(subGroup, schedule, wellState, reportStepIdx, phase_pos, isInjector);
                } else {
                    // Accumulate from this subgroup only if no group guide rate is set for it.
                    if (!guide_rate.has(subGroupName, phase)) {
                        groupTargetReduction[phase_pos] += subGroupEfficiency * subGroupTargetReduction[phase_pos];
                    }
                }
            }
        } else {
            const Group::ProductionCMode& currentGroupControl = group_state.production_control(subGroupName);
            const bool individual_control = (currentGroupControl != Group::ProductionCMode::FLD
                                             && currentGroupControl != Group::ProductionCMode::NONE);
            const int num_group_controlled_wells
                = groupControlledWells(schedule, wellState, group_state, reportStepIdx, subGroupName, "", !isInjector, /*injectionPhaseNotUsed*/Phase::OIL);
            if (individual_control || num_group_controlled_wells == 0) {
                for (int phase = 0; phase < np; phase++) {
                    groupTargetReduction[phase]
                        += subGroupEfficiency * sumWellSurfaceRates(subGroup, schedule, wellState, reportStepIdx, phase, isInjector);
                }
            } else {
                // The subgroup may participate in group control.
                if (!guide_rate.has(subGroupName)) {
                    // Accumulate from this subgroup only if no group guide rate is set for it.
                    for (int phase = 0; phase < np; phase++) {
                        groupTargetReduction[phase] += subGroupEfficiency * subGroupTargetReduction[phase];
                    }
                }
            }
        }
    }

    for (const std::string& wellName : group.wells()) {
        const auto& wellTmp = schedule.getWell(wellName, reportStepIdx);

        if (wellTmp.isProducer() && isInjector)
            continue;

        if (wellTmp.isInjector() && !isInjector)
            continue;

        if (wellTmp.getStatus() == Well::Status::SHUT)
            continue;

        const auto& well_index = wellState.index(wellName);
        if (!well_index.has_value())
            continue;

        if (! wellState.wellIsOwned(well_index.value(), wellName) ) // Only sum once
        {
            continue;
        }

        const Scalar efficiency = wellTmp.getEfficiencyFactor();
        // add contributino from wells not under group control
        const auto& ws = wellState.well(well_index.value());
        if (isInjector) {
            if (ws.injection_cmode != Well::InjectorCMode::GRUP)
                for (int phase = 0; phase < np; phase++) {
                    groupTargetReduction[phase] += ws.surface_rates[phase] * efficiency;
                }
        } else {
            if (ws.production_cmode != Well::ProducerCMode::GRUP)
                for (int phase = 0; phase < np; phase++) {
                    groupTargetReduction[phase] -= ws.surface_rates[phase] * efficiency;
                }
        }
    }
    if (isInjector)
        group_state.update_injection_reduction_rates(group.name(), groupTargetReduction);
    else
        group_state.update_production_reduction_rates(group.name(), groupTargetReduction);
}

template<class Scalar>
void WellGroupHelpers<Scalar>::
updateWellRatesFromGroupTargetScale(const Scalar scale,
                                    const Group& group,
                                    const Schedule& schedule,
                                    const int reportStepIdx,
                                    bool isInjector,
                                    const GroupState<Scalar>& group_state,
                                    WellState<Scalar>& wellState)
{
    for (const std::string& groupName : group.groups()) {
        bool individual_control = false;
        if (isInjector) {
            const Phase all[] = {Phase::WATER, Phase::OIL, Phase::GAS};
            for (Phase phase : all) {
                const Group::InjectionCMode& currentGroupControl
                        = group_state.injection_control(groupName, phase);
                individual_control = individual_control || (currentGroupControl != Group::InjectionCMode::FLD
                                                 && currentGroupControl != Group::InjectionCMode::NONE);
            }
        } else {
            const Group::ProductionCMode& currentGroupControl = group_state.production_control(groupName);
            individual_control = (currentGroupControl != Group::ProductionCMode::FLD
                                             && currentGroupControl != Group::ProductionCMode::NONE);
        }
        if (!individual_control) {
            const Group& groupTmp = schedule.getGroup(groupName, reportStepIdx);
            updateWellRatesFromGroupTargetScale(scale, groupTmp, schedule, reportStepIdx, isInjector, group_state, wellState);
        }
    }

    const int np = wellState.numPhases();
    for (const std::string& wellName : group.wells()) {
        const auto& wellTmp = schedule.getWell(wellName, reportStepIdx);

        if (wellTmp.isProducer() && isInjector)
            continue;

        if (wellTmp.isInjector() && !isInjector)
            continue;

        if (wellTmp.getStatus() == Well::Status::SHUT)
            continue;

        const auto& well_index = wellState.index(wellName);
        if (!well_index.has_value())
            continue;

        if (! wellState.wellIsOwned(well_index.value(), wellName) ) // Only sum once
        {
            continue;
        }

        // scale rates
        auto& ws = wellState.well(well_index.value());
        if (isInjector) {
            if (ws.injection_cmode == Well::InjectorCMode::GRUP)
                for (int phase = 0; phase < np; phase++) {
                    ws.surface_rates[phase] *= scale;
                }
        } else {
            if (ws.production_cmode == Well::ProducerCMode::GRUP) {
                for (int phase = 0; phase < np; phase++) {
                    ws.surface_rates[phase] *= scale;
                }
            }
        }
    }


}

template<class Scalar>
void WellGroupHelpers<Scalar>::
updateVREPForGroups(const Group& group,
                    const Schedule& schedule,
                    const int reportStepIdx,
                    const WellState<Scalar>& wellState,
                    GroupState<Scalar>& group_state)
{
    for (const std::string& groupName : group.groups()) {
        const Group& groupTmp = schedule.getGroup(groupName, reportStepIdx);
        updateVREPForGroups(groupTmp, schedule, reportStepIdx, wellState, group_state);
    }
    const int np = wellState.numPhases();
    Scalar resv = 0.0;
    for (int phase = 0; phase < np; ++phase) {
        resv += sumWellPhaseRates(true,
                                  group,
                                  schedule,
                                  wellState,
                                  reportStepIdx,
                                  phase,
                                  /*isInjector*/ false);
    }
    group_state.update_injection_vrep_rate(group.name(), resv);
}

template<class Scalar>
void WellGroupHelpers<Scalar>::
updateReservoirRatesInjectionGroups(const Group& group,
                                    const Schedule& schedule,
                                    const int reportStepIdx,
                                    const WellState<Scalar>& wellState,
                                    GroupState<Scalar>& group_state)
{
    for (const std::string& groupName : group.groups()) {
        const Group& groupTmp = schedule.getGroup(groupName, reportStepIdx);
        updateReservoirRatesInjectionGroups(groupTmp, schedule, reportStepIdx, wellState, group_state);
    }
    const int np = wellState.numPhases();
    std::vector<Scalar> resv(np, 0.0);
    for (int phase = 0; phase < np; ++phase) {
        resv[phase] = sumWellPhaseRates(true,
                                        group,
                                        schedule,
                                        wellState,
                                        reportStepIdx,
                                        phase,
                                        /*isInjector*/ true);
    }
    group_state.update_injection_reservoir_rates(group.name(), resv);
}

template<class Scalar>
void WellGroupHelpers<Scalar>::
updateSurfaceRatesInjectionGroups(const Group& group,
                                  const Schedule& schedule,
                                  const int reportStepIdx,
                                  const WellState<Scalar>& wellState,
                                  GroupState<Scalar>& group_state)
{
    for (const std::string& groupName : group.groups()) {
        const Group& groupTmp = schedule.getGroup(groupName, reportStepIdx);
        updateSurfaceRatesInjectionGroups(groupTmp, schedule, reportStepIdx, wellState, group_state);
    }
    const int np = wellState.numPhases();
    std::vector<Scalar> rates(np, 0.0);
    for (int phase = 0; phase < np; ++phase) {
        rates[phase] = sumWellPhaseRates(false,
                                        group,
                                        schedule,
                                        wellState,
                                        reportStepIdx,
                                        phase,
                                        /*isInjector*/ true);
    }
    group_state.update_injection_surface_rates(group.name(), rates);
}

template<class Scalar>
void WellGroupHelpers<Scalar>::
updateWellRates(const Group& group,
                const Schedule& schedule,
                const int reportStepIdx,
                const WellState<Scalar>& wellStateNupcol,
                WellState<Scalar>& wellState)
{
    for (const std::string& groupName : group.groups()) {
        const Group& groupTmp = schedule.getGroup(groupName, reportStepIdx);
        updateWellRates(groupTmp, schedule, reportStepIdx, wellStateNupcol, wellState);
    }
    const int np = wellState.numPhases();
    for (const std::string& wellName : group.wells()) {
        std::vector<Scalar> rates(np, 0.0);
        const auto& well_index = wellState.index(wellName);
        if (well_index.has_value()) { // the well is found on this node
            const auto& wellTmp = schedule.getWell(wellName, reportStepIdx);
            int sign = 1;
            // production wellRates are negative. The users of currentWellRates uses the convention in
            // opm-common that production and injection rates are positive.
            if (!wellTmp.isInjector())
                sign = -1;
            const auto& ws = wellStateNupcol.well(well_index.value());
            for (int phase = 0; phase < np; ++phase) {
                rates[phase] = sign * ws.surface_rates[phase];
            }
        }
        wellState.setCurrentWellRates(wellName, rates);
    }
}

template<class Scalar>
void WellGroupHelpers<Scalar>::
updateGroupProductionRates(const Group& group,
                           const Schedule& schedule,
                           const int reportStepIdx,
                           const WellState<Scalar>& wellState,
                           GroupState<Scalar>& group_state)
{
    for (const std::string& groupName : group.groups()) {
        const Group& groupTmp = schedule.getGroup(groupName, reportStepIdx);
        updateGroupProductionRates(groupTmp, schedule, reportStepIdx, wellState, group_state);
    }
    const int np = wellState.numPhases();
    std::vector<Scalar> rates(np, 0.0);
    for (int phase = 0; phase < np; ++phase) {
        rates[phase] = sumWellPhaseRates(false, group, schedule, wellState, reportStepIdx, phase, /*isInjector*/ false);
    }
    group_state.update_production_rates(group.name(), rates);
}

template<class Scalar>
void WellGroupHelpers<Scalar>::
updateREINForGroups(const Group& group,
                    const Schedule& schedule,
                    const int reportStepIdx,
                    const PhaseUsage& pu,
                    const SummaryState& st,
                    const WellState<Scalar>& wellState,
                    GroupState<Scalar>& group_state,
                    bool sum_rank)
{
    const int np = wellState.numPhases();
    for (const std::string& groupName : group.groups()) {
        const Group& groupTmp = schedule.getGroup(groupName, reportStepIdx);
        updateREINForGroups(groupTmp, schedule, reportStepIdx, pu, st, wellState, group_state, sum_rank);
    }

    std::vector<Scalar> rein(np, 0.0);
    for (int phase = 0; phase < np; ++phase) {
        rein[phase] = sumWellPhaseRates(false, group, schedule, wellState, reportStepIdx, phase, /*isInjector*/ false);
    }

    // add import rate and subtract consumption rate for group for gas
    if (sum_rank) {
        if (schedule[reportStepIdx].gconsump().has(group.name())) {
            const auto& gconsump = schedule[reportStepIdx].gconsump().get(group.name(), st);
            if (pu.phase_used[BlackoilPhases::Vapour]) {
                rein[pu.phase_pos[BlackoilPhases::Vapour]] += gconsump.import_rate;
                rein[pu.phase_pos[BlackoilPhases::Vapour]] -= gconsump.consumption_rate;
            }
        }
    }

    group_state.update_injection_rein_rates(group.name(), rein);
}

template<class Scalar>
template <class RegionalValues>
void WellGroupHelpers<Scalar>::
updateGpMaintTargetForGroups(const Group& group,
                             const Schedule& schedule,
                             const RegionalValues& regional_values,
                             const int reportStepIdx,
                             const double dt,
                             const WellState<Scalar>& well_state,
                             GroupState<Scalar>& group_state)
{
    for (const std::string& groupName : group.groups()) {
        const Group& groupTmp = schedule.getGroup(groupName, reportStepIdx);
        updateGpMaintTargetForGroups(groupTmp, schedule, regional_values, reportStepIdx, dt, well_state, group_state);
    }
    const auto& gpm = group.gpmaint();
    if (!gpm)
        return;

    const auto& region = gpm->region();
    if (!region)
        return;

    const auto [name, number] = *region;
    const Scalar error = gpm->pressure_target() - regional_values.at(name)->pressure(number);
    Scalar current_rate = 0.0;
    const auto& pu = well_state.phaseUsage();
    bool injection = true;
    Scalar sign = 1.0;
    switch (gpm->flow_target()) {
        case GPMaint::FlowTarget::RESV_PROD:
        {
            current_rate = -group_state.injection_vrep_rate(group.name());
            injection = false;
            sign = -1.0;
            break;
        }
        case GPMaint::FlowTarget::RESV_OINJ:
        {
            if (pu.phase_used[BlackoilPhases::Liquid])
                current_rate = group_state.injection_reservoir_rates(group.name())[pu.phase_pos[BlackoilPhases::Liquid]];

            break;
        }
        case GPMaint::FlowTarget::RESV_WINJ:
        {
            if (pu.phase_used[BlackoilPhases::Aqua])
                current_rate = group_state.injection_reservoir_rates(group.name())[pu.phase_pos[BlackoilPhases::Aqua]];

            break;
        }
        case GPMaint::FlowTarget::RESV_GINJ:
        {
            if (pu.phase_used[BlackoilPhases::Vapour])
                current_rate = group_state.injection_reservoir_rates(group.name())[pu.phase_pos[BlackoilPhases::Vapour]];
            break;
        }
        case GPMaint::FlowTarget::SURF_OINJ:
        {
            if (pu.phase_used[BlackoilPhases::Liquid])
                current_rate = group_state.injection_surface_rates(group.name())[pu.phase_pos[BlackoilPhases::Liquid]];

            break;
        }
        case GPMaint::FlowTarget::SURF_WINJ:
        {
            if (pu.phase_used[BlackoilPhases::Aqua])
                current_rate = group_state.injection_surface_rates(group.name())[pu.phase_pos[BlackoilPhases::Aqua]];

            break;
        }
        case GPMaint::FlowTarget::SURF_GINJ:
        {
            if (pu.phase_used[BlackoilPhases::Vapour])
                current_rate = group_state.injection_surface_rates(group.name())[pu.phase_pos[BlackoilPhases::Vapour]];

            break;
        }
        default:
            throw std::invalid_argument("Invalid Flow target type in GPMAINT");
    }
    auto& gpmaint_state = group_state.gpmaint(group.name());
    // we only activate gpmaint if pressure is lower than the target regional pressure for injectors
    // (i.e. error > 0) and higher for producers.
    bool activate = (injection && error > 0) || (!injection && error < 0);
    Scalar rate = 0.0;
    if (activate) {
        rate = gpm->rate(gpmaint_state, current_rate, error, dt);
    } else {
        gpm->resetState(gpmaint_state);
    }
    group_state.update_gpmaint_target(group.name(), std::max(Scalar{0.0}, sign * rate));
}

template<class Scalar>
std::map<std::string, Scalar>
WellGroupHelpers<Scalar>::
computeNetworkPressures(const Network::ExtNetwork& network,
                        const WellState<Scalar>& well_state,
                        const GroupState<Scalar>& group_state,
                        const VFPProdProperties<Scalar>& vfp_prod_props,
                        const Schedule& schedule,
                        const int report_time_step)
{
    // TODO: Only dealing with production networks for now.

    if (!network.active()) {
        return {};
    }

    std::map<std::string, Scalar> node_pressures;
    auto roots = network.roots();
    for (const auto& root : roots) {
        // Fixed pressure nodes of the network are the roots of trees.
        // Leaf nodes must correspond to groups in the group structure.
        // Let us first find all leaf nodes of the network. We also
        // create a vector of all nodes, ordered so that a child is
        // always after its parent.
        std::stack<std::string> children;
        std::set<std::string> leaf_nodes;
        std::vector<std::string> root_to_child_nodes;
        //children.push(network.root().name());
        children.push(root.get().name());
        while (!children.empty()) {
            const auto node = children.top();
            children.pop();
            root_to_child_nodes.push_back(node);
            auto branches = network.downtree_branches(node);
            if (branches.empty()) {
                leaf_nodes.insert(node);
            }
            for (const auto& branch : branches) {
                children.push(branch.downtree_node());
            }
        }
        assert(children.empty());

        // Starting with the leaf nodes of the network, get the flow rates
        // from the corresponding groups.
        std::map<std::string, std::vector<Scalar>> node_inflows;
        for (const auto& node : leaf_nodes) {
            node_inflows[node] = group_state.production_rates(node);
            // Add the ALQ amounts to the gas rates if requested.
            if (network.node(node).add_gas_lift_gas()) {
                const auto& group = schedule.getGroup(node, report_time_step);
                for (const std::string& wellname : group.wells()) {
                    const Well& well = schedule.getWell(wellname, report_time_step);
                    if (well.isInjector() || !well_state.isOpen(wellname)) continue;
                    // Here we use the efficiency unconditionally, but if WEFAC item 3
                    // for the well is false (it defaults to true) then we should NOT use
                    // the efficiency factor. Fixing this requires not only changing the
                    // code here, but also:
                    //    - Adding a member to the well for this flag, and setting it in Schedule::handleWEFAC().
                    //    - Making the wells' maximum flows (i.e. not time-averaged by using a efficiency factor)
                    //      available and using those (for wells with WEFAC(3) true only) when accumulating group
                    //      rates, but ONLY for network calculations.
                    const Scalar efficiency = well.getEfficiencyFactor();
                    node_inflows[node][BlackoilPhases::Vapour] += well_state.getALQ(wellname) * efficiency;
                }
            }
        }

        // Accumulate in the network, towards the roots. Note that a
        // root (i.e. fixed pressure node) can still be contributing
        // flow towards other nodes in the network, i.e.  a node is
        // the root of a subtree.
        auto child_to_root_nodes = root_to_child_nodes;
        std::reverse(child_to_root_nodes.begin(), child_to_root_nodes.end());
        for (const auto& node : child_to_root_nodes) {
            const auto upbranch = network.uptree_branch(node);
            if (upbranch) {
                // Add downbranch rates to upbranch.
                std::vector<Scalar>& up = node_inflows[(*upbranch).uptree_node()];
                const std::vector<Scalar>& down = node_inflows[node];
                if (up.empty()) {
                    up = down;
                } else {
                    assert (up.size() == down.size());
                    for (std::size_t ii = 0; ii < up.size(); ++ii) {
                        up[ii] += down[ii];
                    }
                }
            }
        }

        // Going the other way (from roots to leafs), calculate the pressure
        // at each node using VFP tables and rates.
        //std::map<std::string, double> node_pressures;
        for (const auto& node : root_to_child_nodes) {
            auto press = network.node(node).terminal_pressure();
            if (press) {
                node_pressures[node] = *press;
            } else {
                const auto upbranch = network.uptree_branch(node);
                assert(upbranch);
                const Scalar up_press = node_pressures[(*upbranch).uptree_node()];
                const auto vfp_table = (*upbranch).vfp_table();
                if (vfp_table) {
                    // The rates are here positive, but the VFP code expects the
                    // convention that production rates are negative, so we must
                    // take a copy and flip signs.
                    auto rates = node_inflows[node];
                    for (auto& r : rates) { r *= -1.0; }
                    assert(rates.size() == 3);
                    // NB! ALQ in extended network is never implicitly the gas lift rate (GRAT), i.e., the
                    //     gas lift rates only enters the network pressure calculations through the rates
                    //     (e.g., in GOR calculations) unless a branch ALQ is set in BRANPROP.
                    //
                    // @TODO: Standard network
                    Scalar alq = (*upbranch).alq_value().value_or(0.0);
                    node_pressures[node] = vfp_prod_props.bhp(*vfp_table,
                                                            rates[BlackoilPhases::Aqua],
                                                            rates[BlackoilPhases::Liquid],
                                                            rates[BlackoilPhases::Vapour],
                                                            up_press,
                                                            alq,
                                                            0.0, //explicit_wfr
                                                            0.0, //explicit_gfr
                                                            false); //use_expvfp we dont support explicit lookup
#define EXTRA_DEBUG_NETWORK 0
#if EXTRA_DEBUG_NETWORK
                    std::ostringstream oss;
                    oss << "parent: " << (*upbranch).uptree_node() << "  child: " << node
                        << "  rates = [ " << rates[0]*86400 << ", " << rates[1]*86400 << ", " << rates[2]*86400 << " ]"
                        << "  p(parent) = " << up_press/1e5 << "  p(child) = " << node_pressures[node]/1e5 << std::endl;
                    OpmLog::debug(oss.str());
#endif
                } else {
                    // Table number specified as 9999 in the deck, no pressure loss.
                    if (network.node(node).as_choke()){
                        // Node pressure is set to the group THP.
                        node_pressures[node] = group_state.well_group_thp(node);
                    } else {
                            node_pressures[node] = up_press;
                    }
                }
            }
        }
    }
    return node_pressures;
}

template<class Scalar>
GuideRate::RateVector
WellGroupHelpers<Scalar>::
getWellRateVector(const WellState<Scalar>& well_state,
                  const PhaseUsage& pu,
                  const std::string& name)
{
    return getGuideRateVector(well_state.currentWellRates(name), pu);
}

template<class Scalar>
GuideRate::RateVector
WellGroupHelpers<Scalar>::
getProductionGroupRateVector(const GroupState<Scalar>& group_state,
                             const PhaseUsage& pu,
                             const std::string& group_name)
{
    return getGuideRateVector(group_state.production_rates(group_name), pu);
}

template<class Scalar>
Scalar WellGroupHelpers<Scalar>::
getGuideRate(const std::string& name,
             const Schedule& schedule,
             const WellState<Scalar>& wellState,
             const GroupState<Scalar>& group_state,
             const int reportStepIdx,
             const GuideRate* guideRate,
             const GuideRateModel::Target target,
             const PhaseUsage& pu)
{
    if (schedule.hasWell(name, reportStepIdx)) {
        if (guideRate->has(name) || guideRate->hasPotentials(name)) {
            return guideRate->get(name, target, getWellRateVector(wellState, pu, name));
        } else {
            return 0.0;
        }
    }

    if (guideRate->has(name)) {
        return guideRate->get(name, target, getProductionGroupRateVector(group_state, pu, name));
    }

    Scalar totalGuideRate = 0.0;
    const Group& group = schedule.getGroup(name, reportStepIdx);

    for (const std::string& groupName : group.groups()) {
        const Group::ProductionCMode& currentGroupControl = group_state.production_control(groupName);
        if (currentGroupControl == Group::ProductionCMode::FLD
            || currentGroupControl == Group::ProductionCMode::NONE) {
            // accumulate from sub wells/groups
            totalGuideRate += getGuideRate(groupName, schedule, wellState, group_state, reportStepIdx, guideRate, target, pu);
        }
    }

    for (const std::string& wellName : group.wells()) {
        const auto& wellTmp = schedule.getWell(wellName, reportStepIdx);

        if (wellTmp.isInjector())
            continue;

        if (wellTmp.getStatus() == Well::Status::SHUT)
            continue;

        // Only count wells under group control or the ru
        if (!wellState.isProductionGrup(wellName))
            continue;

        totalGuideRate += getGuideRate(wellName, schedule, wellState, group_state,
                                       reportStepIdx, guideRate, target, pu);

    }
    return totalGuideRate;
}

template<class Scalar>
Scalar WellGroupHelpers<Scalar>::
getGuideRateInj(const std::string& name,
                const Schedule& schedule,
                const WellState<Scalar>& wellState,
                const GroupState<Scalar>& group_state,
                const int reportStepIdx,
                const GuideRate* guideRate,
                const GuideRateModel::Target target,
                const Phase& injectionPhase,
                const PhaseUsage& pu)
{
    if (schedule.hasWell(name, reportStepIdx)) {
        return getGuideRate(name, schedule, wellState, group_state,
                            reportStepIdx, guideRate, target, pu);
    }

    if (guideRate->has(name, injectionPhase)) {
        return guideRate->get(name, injectionPhase);
    }

    Scalar totalGuideRate = 0.0;
    const Group& group = schedule.getGroup(name, reportStepIdx);

    for (const std::string& groupName : group.groups()) {
        const Group::InjectionCMode& currentGroupControl
            = group_state.injection_control(groupName, injectionPhase);
        if (currentGroupControl == Group::InjectionCMode::FLD
            || currentGroupControl == Group::InjectionCMode::NONE) {
            // accumulate from sub wells/groups
            totalGuideRate += getGuideRateInj(groupName, schedule, wellState, group_state, reportStepIdx, guideRate, target, injectionPhase, pu);
        }
    }

    for (const std::string& wellName : group.wells()) {
        const auto& wellTmp = schedule.getWell(wellName, reportStepIdx);

        if (!wellTmp.isInjector())
            continue;

        if (wellTmp.getStatus() == Well::Status::SHUT)
            continue;

        // Only count wells under group control or the ru
        if (!wellState.isInjectionGrup(wellName))
            continue;

        totalGuideRate += guideRate->get(wellName, target, getWellRateVector(wellState, pu, wellName));
    }
    return totalGuideRate;
}

template<class Scalar>
int WellGroupHelpers<Scalar>::
groupControlledWells(const Schedule& schedule,
                     const WellState<Scalar>& well_state,
                     const GroupState<Scalar>& group_state,
                     const int report_step,
                     const std::string& group_name,
                     const std::string& always_included_child,
                     const bool is_production_group,
                     const Phase injection_phase)
{
    const Group& group = schedule.getGroup(group_name, report_step);
    int num_wells = 0;
    for (const std::string& child_group : group.groups()) {

        bool included = (child_group == always_included_child);
        if (is_production_group) {
            const auto ctrl = group_state.production_control(child_group);
            included = included || (ctrl == Group::ProductionCMode::FLD) || (ctrl == Group::ProductionCMode::NONE);
        } else {
            const auto ctrl = group_state.injection_control(child_group, injection_phase);
            included = included || (ctrl == Group::InjectionCMode::FLD) || (ctrl == Group::InjectionCMode::NONE);
        }

        if (included) {
            num_wells
                += groupControlledWells(schedule, well_state, group_state, report_step, child_group, always_included_child, is_production_group, injection_phase);
        }
    }
    for (const std::string& child_well : group.wells()) {
        bool included = (child_well == always_included_child);
        if (is_production_group) {
            included = included || well_state.isProductionGrup(child_well);
        } else {
            included = included || well_state.isInjectionGrup(child_well);
        }
        if (included) {
            ++num_wells;
        }
    }
    return num_wells;
}

template<class Scalar>
std::vector<std::string>
WellGroupHelpers<Scalar>::
groupChainTopBot(const std::string& bottom,
                 const std::string& top,
                 const Schedule& schedule,
                 const int report_step)
{
    // Get initial parent, 'bottom' can be a well or a group.
    std::string parent;
    if (schedule.hasWell(bottom, report_step)) {
        parent = schedule.getWell(bottom, report_step).groupName();
    } else {
        parent = schedule.getGroup(bottom, report_step).parent();
    }

    // Build the chain from bottom to top.
    std::vector<std::string> chain;
    chain.push_back(bottom);
    chain.push_back(parent);
    while (parent != top) {
        parent = schedule.getGroup(parent, report_step).parent();
        chain.push_back(parent);
    }
    assert(chain.back() == top);

    // Reverse order and return.
    std::reverse(chain.begin(), chain.end());
    return chain;
}

template<class Scalar>
std::pair<bool, Scalar>
WellGroupHelpers<Scalar>::
checkGroupConstraintsProd(const std::string& name,
                          const std::string& parent,
                          const Group& group,
                          const WellState<Scalar>& wellState,
                          const GroupState<Scalar>& group_state,
                          const int reportStepIdx,
                          const GuideRate* guideRate,
                          const Scalar* rates,
                          const PhaseUsage& pu,
                          const Scalar efficiencyFactor,
                          const Schedule& schedule,
                          const SummaryState& summaryState,
                          const std::vector<Scalar>& resv_coeff,
                          DeferredLogger& deferred_logger)
{
    // When called for a well ('name' is a well name), 'parent'
    // will be the name of 'group'. But if we recurse, 'name' and
    // 'parent' will stay fixed while 'group' will be higher up
    // in the group tree.
    // efficiencyfactor is the well efficiency factor for the first group the well is
    // part of. Later it is the accumulated factor including the group efficiency factor
    // of the child of group.

    // NOTE: if name is a group, the rates-argument is the portion of the 'full' group rates that
    // is potentially available for control by an ancestor, i.e. full rates minus reduction rates
    const Group::ProductionCMode& currentGroupControl = group_state.production_control(group.name());

    if (currentGroupControl == Group::ProductionCMode::FLD || currentGroupControl == Group::ProductionCMode::NONE) {
        // Return if we are not available for parent group.
        if (!group.productionGroupControlAvailable()) {
            return std::make_pair(false, 1);
        }
        // Otherwise: check production share of parent's control.
        const auto& parentGroup = schedule.getGroup(group.parent(), reportStepIdx);
        return checkGroupConstraintsProd(name,
                                         parent,
                                         parentGroup,
                                         wellState,
                                         group_state,
                                         reportStepIdx,
                                         guideRate,
                                         rates,
                                         pu,
                                         efficiencyFactor * group.getGroupEfficiencyFactor(),
                                         schedule,
                                         summaryState,
                                         resv_coeff,
                                         deferred_logger);
    }

    // This can be false for FLD-controlled groups, we must therefore
    // check for FLD first (done above).
    if (!group.isProductionGroup()) {
        return std::make_pair(false, Scalar(1.0));
    }

    // If we are here, we are at the topmost group to be visited in the recursion.
    // This is the group containing the control we will check against.

    // gconsale may adjust the grat target.
    // the adjusted rates is send to the targetCalculator
    Scalar gratTargetFromSales = 0.0;
    if (group_state.has_grat_sales_target(group.name()))
        gratTargetFromSales = group_state.grat_sales_target(group.name());

    WGHelpers::TargetCalculator tcalc(currentGroupControl,
                                      pu,
                                      resv_coeff,
                                      gratTargetFromSales,
                                      group.name(),
                                      group_state,
                                      group.has_gpmaint_control(currentGroupControl));

    WGHelpers::FractionCalculator fcalc(schedule,
                                        wellState,
                                        group_state,
                                        reportStepIdx,
                                        guideRate,
                                        tcalc.guideTargetMode(),
                                        pu,
                                        true,
                                        Phase::OIL);

    auto localFraction = [&](const std::string& child) { return fcalc.localFraction(child, name); };

    auto localReduction = [&](const std::string& group_name)
    {
        const std::vector<Scalar>& groupTargetReductions =
            group_state.production_reduction_rates(group_name);
        return tcalc.calcModeRateFromRates(groupTargetReductions);
    };

    auto localCurrentRate = [&](const std::string& group_name)
    {
        const std::vector<Scalar>& groupSurfaceRates =
            group_state.production_rates(group_name);
        return tcalc.calcModeRateFromRates(groupSurfaceRates);
    };

    std::optional<Group::ProductionControls> ctrl;
    if (!group.has_gpmaint_control(currentGroupControl))
        ctrl = group.productionControls(summaryState);

    const Scalar orig_target = tcalc.groupTarget(ctrl, deferred_logger);
    // Assume we have a chain of groups as follows: BOTTOM -> MIDDLE -> TOP.
    // Then ...
    // TODO finish explanation.
    const Scalar current_rate_available
        = -tcalc.calcModeRateFromRates(rates); // Switch sign since 'rates' are negative for producers.
    const auto chain = groupChainTopBot(name, group.name(), schedule, reportStepIdx);
    // Because 'name' is the last of the elements, and not an ancestor, we subtract one below.
    const std::size_t num_ancestors = chain.size() - 1;
    // we need to find out the level where the current well is applied to the local reduction
    std::size_t local_reduction_level = 0;
    for (std::size_t ii = 1; ii < num_ancestors; ++ii) {
        const int num_gr_ctrl = groupControlledWells(schedule,
                                                     wellState,
                                                     group_state,
                                                     reportStepIdx,
                                                     chain[ii],
                                                     "",
                                                     /*is_producer*/ true,
                                                     /*injectionPhaseNotUsed*/ Phase::OIL);
        if (guideRate->has(chain[ii]) && num_gr_ctrl > 0) {
            local_reduction_level = ii;
        }
    }
    // check whether guide rate is violated
    for (std::size_t ii = 1; ii < num_ancestors; ++ii) {
        if (guideRate->has(chain[ii])) {
            const auto& guided_group = chain[ii];
            const Scalar grefficiency
                        = schedule.getGroup(guided_group, reportStepIdx).getGroupEfficiencyFactor();
            const Scalar currentRateFraction = grefficiency * localCurrentRate(guided_group) / (localCurrentRate(chain[ii-1]));
            const Scalar guiderateFraction = localFraction(guided_group);
            // we add a factor here to avoid switching due to numerical instability
            const Scalar factor = 1.01;
            if (currentRateFraction > (guiderateFraction * factor)) {
                return std::make_pair(true, guiderateFraction / currentRateFraction);
            }
        }
    }
    // Compute portion of target corresponding to current_rate_available
    Scalar target = orig_target;
    for (std::size_t ii = 0; ii < num_ancestors; ++ii) {
        if ((ii == 0) || guideRate->has(chain[ii])) {
            // Apply local reductions only at the control level
            // (top) and for levels where we have a specified
            // group guide rate.
            if (local_reduction_level >= ii ) {
                target -= localReduction(chain[ii]);
            }
            // Add my reduction back at the level where it is included in the local reduction
            if (local_reduction_level == ii ) {
                target += current_rate_available * efficiencyFactor;
            }
        }
        target *= localFraction(chain[ii + 1]);
    }
    // Avoid negative target rates comming from too large local reductions.
    const Scalar target_rate_available = std::max(Scalar(1e-12), target / efficiencyFactor);
    Scalar scale = 1.0;
    if (current_rate_available > 1e-12)
        scale = target_rate_available / current_rate_available;

    return std::make_pair(current_rate_available > target_rate_available, scale);
}

template<class Scalar>
std::pair<bool, Scalar>
WellGroupHelpers<Scalar>::
checkGroupConstraintsInj(const std::string& name,
                         const std::string& parent,
                         const Group& group,
                         const WellState<Scalar>& wellState,
                         const GroupState<Scalar>& group_state,
                         const int reportStepIdx,
                         const GuideRate* guideRate,
                         const Scalar* rates,
                         Phase injectionPhase,
                         const PhaseUsage& pu,
                         const Scalar efficiencyFactor,
                         const Schedule& schedule,
                         const SummaryState& summaryState,
                         const std::vector<Scalar>& resv_coeff,
                         DeferredLogger& deferred_logger)
{
    // When called for a well ('name' is a well name), 'parent'
    // will be the name of 'group'. But if we recurse, 'name' and
    // 'parent' will stay fixed while 'group' will be higher up
    // in the group tree.
    // efficiencyfactor is the well efficiency factor for the first group the well is
    // part of. Later it is the accumulated factor including the group efficiency factor
    // of the child of group.

    // NOTE: if name is a group, the rates-argument is the portion of the 'full' group rates that
    // is potentially available for control by an ancestor, i.e. full rates minus reduction rates

    auto currentGroupControl = group_state.injection_control(group.name(), injectionPhase);

    if (currentGroupControl == Group::InjectionCMode::FLD || currentGroupControl == Group::InjectionCMode::NONE) {
        // Return if we are not available for parent group.
        if (!group.injectionGroupControlAvailable(injectionPhase)) {
            return std::make_pair(false, Scalar(1));
        }
        // Otherwise: check production share of parent's control.
        const auto& parentGroup = schedule.getGroup(group.parent(), reportStepIdx);
        return checkGroupConstraintsInj(name,
                                        parent,
                                        parentGroup,
                                        wellState,
                                        group_state,
                                        reportStepIdx,
                                        guideRate,
                                        rates,
                                        injectionPhase,
                                        pu,
                                        efficiencyFactor * group.getGroupEfficiencyFactor(),
                                        schedule,
                                        summaryState,
                                        resv_coeff,
                                        deferred_logger);
    }

    // This can be false for FLD-controlled groups, we must therefore
    // check for FLD first (done above).
    if (!group.isInjectionGroup()) {
        return std::make_pair(false, Scalar(1.0));
    }

    // If we are here, we are at the topmost group to be visited in the recursion.
    // This is the group containing the control we will check against.
    Scalar sales_target = 0;
    if (schedule[reportStepIdx].gconsale().has(group.name())) {
        const auto& gconsale = schedule[reportStepIdx].gconsale().get(group.name(), summaryState);
        sales_target = gconsale.sales_target;
    }
    WGHelpers::InjectionTargetCalculator tcalc(currentGroupControl,
                                               pu,
                                               resv_coeff,
                                               group.name(),
                                               sales_target,
                                               group_state,
                                               injectionPhase,
                                               group.has_gpmaint_control(injectionPhase,
                                                                         currentGroupControl),
                                               deferred_logger);

    WGHelpers::FractionCalculator fcalc(schedule,
                                        wellState,
                                        group_state,
                                        reportStepIdx,
                                        guideRate,
                                        tcalc.guideTargetMode(),
                                        pu,
                                        false,
                                        injectionPhase);

    auto localFraction = [&](const std::string& child) { return fcalc.localFraction(child, name); };

    auto localReduction = [&](const std::string& group_name)
    {
        const std::vector<Scalar>& groupTargetReductions =
            group_state.injection_reduction_rates(group_name);
        return tcalc.calcModeRateFromRates(groupTargetReductions);
    };

    auto localCurrentRate = [&](const std::string& group_name)
    {
        const std::vector<Scalar>& groupSurfaceRates =
            group_state.injection_surface_rates(group_name);
        return tcalc.calcModeRateFromRates(groupSurfaceRates);
    };

    std::optional<Group::InjectionControls> ctrl;
    if (!group.has_gpmaint_control(injectionPhase, currentGroupControl))
        ctrl = group.injectionControls(injectionPhase, summaryState);


    const Scalar orig_target = tcalc.groupTarget(ctrl, deferred_logger);
    // Assume we have a chain of groups as follows: BOTTOM -> MIDDLE -> TOP.
    // Then ...
    // TODO finish explanation.
    const Scalar current_rate_available
        = tcalc.calcModeRateFromRates(rates); // Switch sign since 'rates' are negative for producers.
    const auto chain = groupChainTopBot(name, group.name(), schedule, reportStepIdx);
    // Because 'name' is the last of the elements, and not an ancestor, we subtract one below.
    const std::size_t num_ancestors = chain.size() - 1;
    // we need to find out the level where the current well is applied to the local reduction
    std::size_t local_reduction_level = 0;
    for (std::size_t ii = 1; ii < num_ancestors; ++ii) {
        const int num_gr_ctrl = groupControlledWells(schedule,
                                                     wellState,
                                                             group_state,
                                                             reportStepIdx,
                                                             chain[ii],
                                                             "",
                                                             /*is_producer*/ false,
                                                             injectionPhase);
        if (guideRate->has(chain[ii], injectionPhase) && num_gr_ctrl > 0) {
            local_reduction_level = ii;
        }
    }

    // check whether guide rate is violated
    for (std::size_t ii = 1; ii < num_ancestors; ++ii) {
        if (guideRate->has(chain[ii], injectionPhase)) {
            const auto& guided_group = chain[ii];
            const Scalar grefficiency
                        = schedule.getGroup(guided_group, reportStepIdx).getGroupEfficiencyFactor();
            const Scalar currentRateFraction = grefficiency * localCurrentRate(guided_group) /
                                               localCurrentRate(chain[ii-1]);
            const Scalar guiderateFraction = localFraction(guided_group);
            // we add a factor here to avoid switching due to numerical instability
            const Scalar factor = 1.01;
            if (currentRateFraction > (guiderateFraction * factor)) {
                return std::make_pair(true, guiderateFraction / currentRateFraction);
            }
        }
    }

    // Compute portion of target corresponding to current_rate_available
    Scalar target = orig_target;
    for (std::size_t ii = 0; ii < num_ancestors; ++ii) {
        if ((ii == 0) || guideRate->has(chain[ii], injectionPhase)) {
            // Apply local reductions only at the control level
            // (top) and for levels where we have a specified
            // group guide rate.
            if (local_reduction_level >= ii ) {
                target -= localReduction(chain[ii]);
            }

            // Add my reduction back at the level where it is included in the local reduction
            if (local_reduction_level == ii ) {
                target += current_rate_available * efficiencyFactor;
            }
        }
        target *= localFraction(chain[ii + 1]);
    }
    // Avoid negative target rates comming from too large local reductions.
    const Scalar target_rate_available = std::max(Scalar(1e-12), target / efficiencyFactor);
    Scalar scale = 1.0;
    if (current_rate_available > 1e-12)
        scale = target_rate_available / current_rate_available;

    return std::make_pair(current_rate_available > target_rate_available, scale);
}

template<class Scalar>
std::pair<std::optional<std::string>, Scalar>
WellGroupHelpers<Scalar>::
worstOffendingWell(const Group& group,
                   const Schedule& schedule,
                   const int reportStepIdx,
                   const Group::ProductionCMode& offendedControl,
                   const PhaseUsage& pu,
                   const Parallel::Communication& comm,
                   const WellState<Scalar>& wellState,
                   DeferredLogger& deferred_logger)
{
    std::pair<std::optional<std::string>, Scalar> offending_well {std::nullopt, 0.0};
    for (const std::string& child_group : group.groups()) {
        const auto& this_group = schedule.getGroup(child_group, reportStepIdx);
        const auto & offending_well_this = worstOffendingWell(this_group,
                                                             schedule,
                                                             reportStepIdx,
                                                             offendedControl,
                                                             pu,
                                                             comm,
                                                             wellState,
                                                             deferred_logger);
        if (offending_well_this.second > offending_well.second) {
            offending_well = offending_well_this;
        }
    }

    for (const std::string& child_well : group.wells()) {

        const auto& well_index = wellState.index(child_well);
        Scalar violating_rate = 0.0;
        Scalar prefered_rate = 0.0;
        if (well_index.has_value() && wellState.wellIsOwned(well_index.value(), child_well))
        {
            const auto& ws = wellState.well(child_well);
            switch (offendedControl){
                case Group::ProductionCMode::ORAT:
                    violating_rate = ws.surface_rates[pu.phase_pos[BlackoilPhases::Liquid]];
                    break;
                case Group::ProductionCMode::GRAT:
                    violating_rate = ws.surface_rates[pu.phase_pos[BlackoilPhases::Vapour]];
                    break;
                case Group::ProductionCMode::WRAT:
                    violating_rate = ws.surface_rates[pu.phase_pos[BlackoilPhases::Aqua]];
                    break;
                case Group::ProductionCMode::LRAT:
                    assert(pu.phase_used[BlackoilPhases::Liquid]);
                    assert(pu.phase_used[BlackoilPhases::Aqua]);
                    violating_rate = ws.surface_rates[pu.phase_pos[BlackoilPhases::Liquid]] +
                                     ws.surface_rates[pu.phase_pos[BlackoilPhases::Aqua]];
                    break;
                case Group::ProductionCMode::RESV:
                    for (int p = 0; p < pu.num_phases; ++p) {
                        violating_rate += ws.reservoir_rates[p];
                    }
                    break;
                case Group::ProductionCMode::NONE:
                    break;
                case Group::ProductionCMode::FLD:
                    break;
                case Group::ProductionCMode::PRBL:
                    OPM_DEFLOG_THROW(std::runtime_error,
                                     "Group " + group.name() +
                                     "GroupProductionCMode PRBL not implemented", deferred_logger);
                    break;
                case Group::ProductionCMode::CRAT:
                    OPM_DEFLOG_THROW(std::runtime_error,
                                     "Group " + group.name() +
                                     "GroupProductionCMode CRAT not implemented", deferred_logger);
                    break;
            }
            const auto preferred_phase = schedule.getWell(child_well, reportStepIdx).getPreferredPhase();
             switch (preferred_phase) {
                case Phase::OIL:
                    prefered_rate = ws.surface_rates[pu.phase_pos[BlackoilPhases::Liquid]];
                    break;
                case Phase::GAS:
                    prefered_rate = ws.surface_rates[pu.phase_pos[BlackoilPhases::Vapour]];
                    break;
                case Phase::WATER:
                    prefered_rate = ws.surface_rates[pu.phase_pos[BlackoilPhases::Aqua]];
                    break;
                default:
                    // No others supported.
                    break;
                }
        }
        violating_rate = comm.sum(violating_rate);
        if (violating_rate < 0 ) { // only check producing wells
            prefered_rate = comm.sum(prefered_rate);
            Scalar fraction = prefered_rate < -1e-16 ? violating_rate / prefered_rate : 1.0;
            if ( fraction > offending_well.second) {
                    offending_well = {child_well, fraction};
            }
        }
    }
    return offending_well;
}

template<class Scalar>
template <class AverageRegionalPressureType>
void WellGroupHelpers<Scalar>::
setRegionAveragePressureCalculator(const Group& group,
                                   const Schedule& schedule,
                                   const int reportStepIdx,
                                   const FieldPropsManager& fp,
                                   const PhaseUsage& pu,
                                   std::map<std::string, std::unique_ptr<AverageRegionalPressureType>>& regionalAveragePressureCalculator)
{
    for (const std::string& groupName : group.groups()) {
        setRegionAveragePressureCalculator( schedule.getGroup(groupName, reportStepIdx), schedule,
                                            reportStepIdx, fp, pu, regionalAveragePressureCalculator);
    }
    const auto& gpm = group.gpmaint();
    if (!gpm)
        return;

    const auto& reg = gpm->region();
    if (!reg)
        return;

    if (regionalAveragePressureCalculator.count(reg->first) == 0) {
        const std::string name = (reg->first.rfind("FIP", 0) == 0) ? reg->first : "FIP" + reg->first;
        const auto& fipnum = fp.get_int(name);
        regionalAveragePressureCalculator[reg->first] = std::make_unique<AverageRegionalPressureType>(pu,fipnum);
    }
}

template<class Scalar>
void WellGroupHelpers<Scalar>::
updateGuideRates(const Group& group,
                 const Schedule& schedule,
                 const SummaryState& summary_state,
                 const PhaseUsage& pu,
                 const int report_step,
                 const double sim_time,
                 WellState<Scalar>& well_state,
                 const GroupState<Scalar>& group_state,
                 const Parallel::Communication& comm,
                 GuideRate* guide_rate,
                 std::vector<Scalar>& pot,
                 DeferredLogger& deferred_logger)
{
    guide_rate->updateGuideRateExpiration(sim_time, report_step);
    updateGuideRateForProductionGroups(group, schedule, pu, report_step, sim_time, well_state, group_state, comm, guide_rate, pot);
    updateGuideRatesForInjectionGroups(group, schedule, summary_state, pu, report_step, well_state, group_state, guide_rate, deferred_logger);
    updateGuideRatesForWells(schedule, pu, report_step, sim_time, well_state, comm, guide_rate);
}

template<class Scalar>
void WellGroupHelpers<Scalar>::
updateGuideRateForProductionGroups(const Group& group,
                                   const Schedule& schedule,
                                   const PhaseUsage& pu,
                                   const int reportStepIdx,
                                   const double& simTime,
                                   WellState<Scalar>& wellState,
                                   const GroupState<Scalar>& group_state,
                                   const Parallel::Communication& comm,
                                   GuideRate* guideRate,
                                   std::vector<Scalar>& pot)
{
    const int np = pu.num_phases;
    for (const std::string& groupName : group.groups()) {
        std::vector<Scalar> thisPot(np, 0.0);
        const Group& groupTmp = schedule.getGroup(groupName, reportStepIdx);

        // Note that group effiency factors for groupTmp are applied in updateGuideRateForGroups
        updateGuideRateForProductionGroups(groupTmp, schedule, pu, reportStepIdx, simTime, wellState, group_state, comm, guideRate, thisPot);

        // accumulate group contribution from sub group unconditionally
        const auto currentGroupControl = group_state.production_control(groupName);
        if (currentGroupControl != Group::ProductionCMode::FLD
                && currentGroupControl != Group::ProductionCMode::NONE) {
            continue;
        }

        // Apply group efficiency factor for this goup
        auto gefac = groupTmp.getGroupEfficiencyFactor();

        for (int phase = 0; phase < np; phase++) {
            pot[phase] += gefac*thisPot[phase];
        }

    }
    for (const std::string& wellName : group.wells()) {
        const auto& wellTmp = schedule.getWell(wellName, reportStepIdx);
        const auto wefac = wellTmp.getEfficiencyFactor();

        if (wellTmp.isInjector())
            continue;

        if (wellTmp.getStatus() == Well::Status::SHUT)
            continue;
        const auto& well_index = wellState.index(wellName);
        if (!well_index.has_value()) // the well is not found
            continue;

        if (! wellState.wellIsOwned(well_index.value(), wellName) ) // Only sum once
        {
            continue;
        }

        const auto& ws = wellState.well(well_index.value());
        // add contribution from wells unconditionally
        for (int phase = 0; phase < np; phase++) {
            pot[phase] += wefac * ws.well_potentials[phase];
        }
    }

    std::array<Scalar,3> potentials{};
    auto& [oilPot, gasPot, waterPot] = potentials;
    if (pu.phase_used[BlackoilPhases::Liquid])
        oilPot = pot[pu.phase_pos[BlackoilPhases::Liquid]];

    if (pu.phase_used[BlackoilPhases::Vapour])
        gasPot = pot[pu.phase_pos[BlackoilPhases::Vapour]];

    if (pu.phase_used[BlackoilPhases::Aqua])
        waterPot = pot[pu.phase_pos[BlackoilPhases::Aqua]];

    comm.sum(potentials.data(), potentials.size());
    const UnitSystem& unit_system = schedule.getUnits();
    oilPot = unit_system.from_si(UnitSystem::measure::liquid_surface_rate, oilPot);
    waterPot = unit_system.from_si(UnitSystem::measure::liquid_surface_rate, waterPot);
    gasPot = unit_system.from_si(UnitSystem::measure::gas_surface_rate, gasPot);
    guideRate->compute(group.name(), reportStepIdx, simTime, oilPot, gasPot, waterPot);
}

template<class Scalar>
void WellGroupHelpers<Scalar>::
updateGuideRatesForWells(const Schedule& schedule,
                         const PhaseUsage& pu,
                         const int reportStepIdx,
                         const double& simTime,
                         const WellState<Scalar>& wellState,
                         const Parallel::Communication& comm,
                         GuideRate* guideRate)
{
    for (const auto& well : schedule.getWells(reportStepIdx)) {
        std::array<Scalar,3> potentials{};
        auto& [oilpot, gaspot, waterpot] = potentials;

        const auto& well_index = wellState.index(well.name());
        if (well_index.has_value() && wellState.wellIsOwned(well_index.value(), well.name()))
        {
            // the well is found and owned
            const auto& ws = wellState.well(well_index.value());
            const auto& wpot = ws.well_potentials;
            if (pu.phase_used[BlackoilPhases::Liquid] > 0)
                oilpot = wpot[pu.phase_pos[BlackoilPhases::Liquid]];

            if (pu.phase_used[BlackoilPhases::Vapour] > 0)
                gaspot = wpot[pu.phase_pos[BlackoilPhases::Vapour]];

            if (pu.phase_used[BlackoilPhases::Aqua] > 0)
                waterpot = wpot[pu.phase_pos[BlackoilPhases::Aqua]];
        }
        comm.sum(potentials.data(), potentials.size());
        const UnitSystem& unit_system = schedule.getUnits();
        oilpot = unit_system.from_si(UnitSystem::measure::liquid_surface_rate, oilpot);
        waterpot = unit_system.from_si(UnitSystem::measure::liquid_surface_rate, waterpot);
        gaspot = unit_system.from_si(UnitSystem::measure::gas_surface_rate, gaspot);
        guideRate->compute(well.name(), reportStepIdx, simTime, oilpot, gaspot, waterpot);
    }
}

template<class Scalar>
using AvgP = RegionAverageCalculator::
    AverageRegionalPressure<BlackOilFluidSystem<Scalar>,std::vector<int>>;

template<class Scalar>
using AvgPMap = std::map<std::string, std::unique_ptr<AvgP<Scalar>>>;

#define INSTANTIATE_TYPE(T)                                                   \
    template class WellGroupHelpers<T>;                                       \
    template void WellGroupHelpers<T>::                                       \
        updateGpMaintTargetForGroups<AvgPMap<T>>(const Group&,                \
                                                 const Schedule&,             \
                                                 const AvgPMap<T>&,           \
                                                 int,                         \
                                                 double,                      \
                                                 const WellState<T>&,         \
                                                 GroupState<T>&);             \
    template void WellGroupHelpers<T>::                                       \
        setRegionAveragePressureCalculator<AvgP<T>>(const Group&,             \
                                                    const Schedule&,          \
                                                    const int,                \
                                                    const FieldPropsManager&, \
                                                    const PhaseUsage&,        \
                                                    AvgPMap<T>&);

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} // namespace Opm::WellGroupHelpers
