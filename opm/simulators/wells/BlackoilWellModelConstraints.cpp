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
#include <opm/simulators/wells/BlackoilWellModelConstraints.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <opm/simulators/wells/BlackoilWellModelGeneric.hpp>
#include <opm/simulators/wells/WellGroupHelpers.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>

#include <fmt/format.h>

#include <stdexcept>

namespace Opm {

bool BlackoilWellModelConstraints::
hasTHPConstraints() const
{
    int local_result = false;
    for (const auto& well : wellModel_.genericWells()) {
        if (well->wellHasTHPConstraints(wellModel_.summaryState())) {
            local_result = true;
        }
    }
    return wellModel_.comm().max(local_result);
}

std::pair<Group::InjectionCMode, double>
BlackoilWellModelConstraints::
checkGroupInjectionConstraints(const Group& group,
                               const int reportStepIdx,
                               const Phase& phase) const
{
    const auto& well_state = wellModel_.wellState();
    const auto& pu = wellModel_.phaseUsage();

    int phasePos;
    if (phase == Phase::GAS && pu.phase_used[BlackoilPhases::Vapour] )
        phasePos = pu.phase_pos[BlackoilPhases::Vapour];
    else if (phase == Phase::OIL && pu.phase_used[BlackoilPhases::Liquid])
        phasePos = pu.phase_pos[BlackoilPhases::Liquid];
    else if (phase == Phase::WATER && pu.phase_used[BlackoilPhases::Aqua] )
        phasePos = pu.phase_pos[BlackoilPhases::Aqua];
    else
        OPM_THROW(std::runtime_error, "Unknown phase" );

    auto currentControl = wellModel_.groupState().injection_control(group.name(), phase);
    if (group.has_control(phase, Group::InjectionCMode::RATE))
    {
        if (currentControl != Group::InjectionCMode::RATE)
        {
            double current_rate = 0.0;
            current_rate += WellGroupHelpers::sumWellSurfaceRates(group, wellModel_.schedule(), well_state, reportStepIdx, phasePos, /*isInjector*/true);

            // sum over all nodes
            current_rate = wellModel_.comm().sum(current_rate);

            const auto& controls = group.injectionControls(phase, wellModel_.summaryState());
            double target = controls.surface_max_rate;

            if (group.has_gpmaint_control(phase, Group::InjectionCMode::RATE))
                target = wellModel_.groupState().gpmaint_target(group.name());

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
            current_rate += WellGroupHelpers::sumWellResRates(group, wellModel_.schedule(), well_state, reportStepIdx, phasePos, /*isInjector*/true);
            // sum over all nodes
            current_rate = wellModel_.comm().sum(current_rate);

            const auto& controls = group.injectionControls(phase, wellModel_.summaryState());
            double target = controls.resv_max_rate;

            if (group.has_gpmaint_control(phase, Group::InjectionCMode::RESV))
                target = wellModel_.groupState().gpmaint_target(group.name());

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
            const auto& controls = group.injectionControls(phase, wellModel_.summaryState());
            const Group& groupRein = wellModel_.schedule().getGroup(controls.reinj_group, reportStepIdx);
            production_Rate += WellGroupHelpers::sumWellSurfaceRates(groupRein, wellModel_.schedule(),
                                                                     well_state, reportStepIdx,
                                                                     phasePos, /*isInjector*/false);

            // sum over all nodes
            production_Rate = wellModel_.comm().sum(production_Rate);

            double current_rate = 0.0;
            current_rate += WellGroupHelpers::sumWellSurfaceRates(group, wellModel_.schedule(),
                                                                  well_state, reportStepIdx,
                                                                  phasePos, /*isInjector*/true);

            // sum over all nodes
            current_rate = wellModel_.comm().sum(current_rate);

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
            const auto& controls = group.injectionControls(phase, wellModel_.summaryState());
            const Group& groupVoidage = wellModel_.schedule().getGroup(controls.voidage_group, reportStepIdx);
            voidage_rate += WellGroupHelpers::sumWellResRates(groupVoidage, wellModel_.schedule(),
                                                              well_state, reportStepIdx,
                                                              pu.phase_pos[BlackoilPhases::Aqua], false);
            voidage_rate += WellGroupHelpers::sumWellResRates(groupVoidage, wellModel_.schedule(),
                                                              well_state, reportStepIdx,
                                                              pu.phase_pos[BlackoilPhases::Liquid], false);
            voidage_rate += WellGroupHelpers::sumWellResRates(groupVoidage, wellModel_.schedule(),
                                                              well_state, reportStepIdx,
                                                              pu.phase_pos[BlackoilPhases::Vapour], false);

            // sum over all nodes
            voidage_rate = wellModel_.comm().sum(voidage_rate);

            double total_rate = 0.0;
            total_rate += WellGroupHelpers::sumWellResRates(group, wellModel_.schedule(),
                                                            well_state, reportStepIdx,
                                                            pu.phase_pos[BlackoilPhases::Aqua], true);
            total_rate += WellGroupHelpers::sumWellResRates(group, wellModel_.schedule(),
                                                            well_state, reportStepIdx,
                                                            pu.phase_pos[BlackoilPhases::Liquid], true);
            total_rate += WellGroupHelpers::sumWellResRates(group, wellModel_.schedule(),
                                                            well_state, reportStepIdx,
                                                            pu.phase_pos[BlackoilPhases::Vapour], true);

            // sum over all nodes
            total_rate = wellModel_.comm().sum(total_rate);

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
BlackoilWellModelConstraints::
checkGroupProductionConstraints(const Group& group,
                                const int reportStepIdx,
                                DeferredLogger& deferred_logger) const
{
    const auto& well_state = wellModel_.wellState();
    const auto& pu = wellModel_.phaseUsage();

    const auto controls = group.productionControls(wellModel_.summaryState());
    const Group::ProductionCMode& currentControl = wellModel_.groupState().production_control(group.name());

    if (group.has_control(Group::ProductionCMode::ORAT))
    {
        if (currentControl != Group::ProductionCMode::ORAT)
        {
            double current_rate = 0.0;
            current_rate += WellGroupHelpers::sumWellSurfaceRates(group, wellModel_.schedule(),
                                                                  well_state, reportStepIdx,
                                                                  pu.phase_pos[BlackoilPhases::Liquid], false);

            // sum over all nodes
            current_rate = wellModel_.comm().sum(current_rate);

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
            current_rate += WellGroupHelpers::sumWellSurfaceRates(group, wellModel_.schedule(),
                                                                  well_state, reportStepIdx,
                                                                  pu.phase_pos[BlackoilPhases::Aqua], false);

            // sum over all nodes
            current_rate = wellModel_.comm().sum(current_rate);

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
            current_rate += WellGroupHelpers::sumWellSurfaceRates(group, wellModel_.schedule(),
                                                                  well_state, reportStepIdx,
                                                                  pu.phase_pos[BlackoilPhases::Vapour], false);

            // sum over all nodes
            current_rate = wellModel_.comm().sum(current_rate);
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
            current_rate += WellGroupHelpers::sumWellSurfaceRates(group, wellModel_.schedule(),
                                                                  well_state, reportStepIdx,
                                                                  pu.phase_pos[BlackoilPhases::Liquid], false);
            current_rate += WellGroupHelpers::sumWellSurfaceRates(group, wellModel_.schedule(),
                                                                  well_state, reportStepIdx,
                                                                  pu.phase_pos[BlackoilPhases::Aqua], false);

            // sum over all nodes
            current_rate = wellModel_.comm().sum(current_rate);

            bool skip = false;
            if (controls.liquid_target == controls.oil_target) {
                double current_water_rate = WellGroupHelpers::sumWellSurfaceRates(group, wellModel_.schedule(),
                                                                                  well_state, reportStepIdx,
                                                                                  pu.phase_pos[BlackoilPhases::Aqua], false);
                current_water_rate = wellModel_.comm().sum(current_water_rate);
                if (std::abs(current_water_rate) < 1e-12) {
                    skip = true;
                    deferred_logger.debug("LRAT_ORAT_GROUP", "GROUP " + group.name() + " The LRAT target is equal the ORAT target and the water rate is zero, skip checking LRAT");
                }
            }

            if (!skip && controls.liquid_target < current_rate ) {
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
            current_rate += WellGroupHelpers::sumWellResRates(group, wellModel_.schedule(),
                                                              well_state, reportStepIdx,
                                                              pu.phase_pos[BlackoilPhases::Aqua], false);
            current_rate += WellGroupHelpers::sumWellResRates(group, wellModel_.schedule(),
                                                              well_state, reportStepIdx,
                                                              pu.phase_pos[BlackoilPhases::Liquid], false);
            current_rate += WellGroupHelpers::sumWellResRates(group, wellModel_.schedule(),
                                                              well_state, reportStepIdx,
                                                              pu.phase_pos[BlackoilPhases::Vapour], false);

            // sum over all nodes
            current_rate = wellModel_.comm().sum(current_rate);

            double target = controls.resv_target;
            if (group.has_gpmaint_control(Group::ProductionCMode::RESV))
                target = wellModel_.groupState().gpmaint_target(group.name());

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

bool BlackoilWellModelConstraints::
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
            const auto& check = this->checkGroupInjectionConstraints(group,
                                                                     reportStepIdx,  phase);
            if (check.first != Group::InjectionCMode::NONE) {
                return true;
            }
        }
    }
    if (group.isProductionGroup()) {
        const auto& check = this->checkGroupProductionConstraints(group,
                                                                  reportStepIdx,
                                                                  deferred_logger);
        if (check.first != Group::ProductionCMode::NONE)
        {
            return true;
        }
    }

    // call recursively down the group hierarchy
    bool violated = false;
    for (const std::string& groupName : group.groups()) {
        const auto& grp = wellModel_.schedule().getGroup(groupName, reportStepIdx);
        violated = violated || this->checkGroupConstraints(grp, reportStepIdx, deferred_logger);
    }
    return violated;
}

void BlackoilWellModelConstraints::
actionOnBrokenConstraints(const Group& group,
                          const Group::InjectionCMode& newControl,
                          const Phase& controlPhase,
                          GroupState& group_state,
                          DeferredLogger& deferred_logger) const
{
    auto oldControl = wellModel_.groupState().injection_control(group.name(), controlPhase);

    if (oldControl != newControl) {
        const std::string from = Group::InjectionCMode2String(oldControl);
        group_state.injection_control(group.name(), controlPhase, newControl);
        if (wellModel_.comm().rank() == 0) {
            auto msg = fmt::format("Switching injection control mode for group {} from {} to {}",
                                   group.name(),
                                   Group::InjectionCMode2String(oldControl),
                                   Group::InjectionCMode2String(newControl));
            deferred_logger.debug(msg);
        }
    }
}

void BlackoilWellModelConstraints::
actionOnBrokenConstraints(const Group& group,
                          const Group::GroupLimitAction group_limit_action,
                          const Group::ProductionCMode& newControl,
                          GroupState& group_state,
                          DeferredLogger& deferred_logger) const
{
    const Group::ProductionCMode oldControl = wellModel_.groupState().production_control(group.name());

    std::string ss;
    switch(group_limit_action.allRates) {
    case Group::ExceedAction::NONE: {
        if (oldControl != newControl && oldControl != Group::ProductionCMode::NONE) {
            if ((group_limit_action.water == Group::ExceedAction::RATE && newControl == Group::ProductionCMode::WRAT) ||
                (group_limit_action.gas == Group::ExceedAction::RATE && newControl == Group::ProductionCMode::GRAT) ||
                (group_limit_action.liquid == Group::ExceedAction::RATE && newControl == Group::ProductionCMode::LRAT)) {
                group_state.production_control(group.name(), newControl);
                ss = fmt::format("Switching production control mode for group {} from {} to {}",
                                 group.name(),
                                 Group::ProductionCMode2String(oldControl),
                                 Group::ProductionCMode2String(newControl));
            }
            else {
                ss = fmt::format("Procedure on exceeding {} limit is NONE for group {}. Nothing is done.",
                                 Group::ProductionCMode2String(oldControl),
                                 group.name());
            }
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
            group_state.production_control(group.name(), newControl);
            ss = fmt::format("Switching production control mode for group {} from {} to {}",
                             group.name(),
                             Group::ProductionCMode2String(oldControl),
                             Group::ProductionCMode2String(newControl));
        }
        break;
    }
    default:
        throw("Invalid procedure for maximum rate limit selected for group" + group.name());
    }

    if (!ss.empty() && wellModel_.comm().rank() == 0)
        deferred_logger.debug(ss);
}

bool BlackoilWellModelConstraints::
updateGroupIndividualControl(const Group& group,
                             const int reportStepIdx,
                             std::map<std::pair<std::string,Opm::Phase>,std::string>& switched_inj,
                             std::map<std::string, std::string>& switched_prod,
                             GroupState& group_state,
                             WellState& well_state,
                             DeferredLogger& deferred_logger) const
{
    bool changed = false;
    if (group.isInjectionGroup())
    {
        const Phase all[] = {Phase::WATER, Phase::OIL, Phase::GAS};
        for (Phase phase : all) {
            if (!group.hasInjectionControl(phase)) {
                continue;
            }
            const auto& changed_this = this->checkGroupInjectionConstraints(group,
                                                                            reportStepIdx,
                                                                            phase);
            if (changed_this.first != Group::InjectionCMode::NONE)
            {
                switched_inj.insert_or_assign({group.name(), phase},
                                     Group::InjectionCMode2String(changed_this.first));
                this->actionOnBrokenConstraints(group, changed_this.first, phase,
                                                group_state, deferred_logger);
                WellGroupHelpers::updateWellRatesFromGroupTargetScale(changed_this.second,
                                                                      group,
                                                                      wellModel_.schedule(),
                                                                      reportStepIdx,
                                                                      /* isInjector */ false,
                                                                      wellModel_.groupState(),
                                                                      well_state);
                changed = true;
            }
        }
    }
    if (group.isProductionGroup()) {
        const auto& changed_this = this->checkGroupProductionConstraints(group,
                                                                         reportStepIdx,
                                                                         deferred_logger);
        const auto controls = group.productionControls(wellModel_.summaryState());
        if (changed_this.first != Group::ProductionCMode::NONE)
        {
            switched_prod.insert_or_assign(group.name(),
                                  Group::ProductionCMode2String(changed_this.first));

            this->actionOnBrokenConstraints(group,
                                            controls.group_limit_action,
                                            changed_this.first,
                                            group_state, deferred_logger);
            WellGroupHelpers::updateWellRatesFromGroupTargetScale(changed_this.second,
                                                                  group,
                                                                  wellModel_.schedule(),
                                                                  reportStepIdx,
                                                                  /* isInjector */ false,
                                                                  wellModel_.groupState(),
                                                                  well_state);
            changed = true;
        }
    }

    return changed;
}

}
