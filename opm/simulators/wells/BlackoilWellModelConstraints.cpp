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

#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>

#include <opm/simulators/wells/BlackoilWellModelGeneric.hpp>
#include <opm/simulators/wells/WellGroupHelpers.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>

#include <fmt/format.h>

#include <stdexcept>

namespace Opm {

template<typename Scalar, typename IndexTraits>
std::pair<Group::InjectionCMode, Scalar>
BlackoilWellModelConstraints<Scalar, IndexTraits>::
checkGroupInjectionConstraints(const Group& group,
                               const int reportStepIdx,
                               const Phase& phase) const
{
    const auto& well_state = wellModel_.wellState();
    const auto& pu = wellModel_.phaseUsage();

    int phasePos;
    if (phase == Phase::GAS && pu.phaseIsActive(gasPhaseIdx) )
        phasePos = pu.canonicalToActivePhaseIdx(gasPhaseIdx);
    else if (phase == Phase::OIL && pu.phaseIsActive(oilPhaseIdx) )
        phasePos = pu.canonicalToActivePhaseIdx(oilPhaseIdx);
    else if (phase == Phase::WATER && pu.phaseIsActive(waterPhaseIdx) )
        phasePos = pu.canonicalToActivePhaseIdx(waterPhaseIdx);
    else
        OPM_THROW(std::runtime_error, "Unknown phase" );

    auto currentControl = wellModel_.groupState().injection_control(group.name(), phase);
    if (group.has_control(phase, Group::InjectionCMode::RATE))
    {
        if (currentControl != Group::InjectionCMode::RATE)
        {
            Scalar current_rate = 0.0;
            current_rate += WellGroupHelpersType::sumWellSurfaceRates(group,
                                                                      wellModel_.schedule(),
                                                                      well_state,
                                                                      reportStepIdx,
                                                                      phasePos,
                                                                      /*isInjector*/true,
                                                                      wellModel_.summaryState());

            // sum over all nodes
            current_rate = wellModel_.comm().sum(current_rate);

            const auto& controls = group.injectionControls(phase, wellModel_.summaryState());
            Scalar target = controls.surface_max_rate;

            if (group.has_gpmaint_control(phase, Group::InjectionCMode::RATE))
                target = wellModel_.groupState().gpmaint_target(group.name());

            if (target < current_rate) {
                Scalar scale = 1.0;
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
            Scalar current_rate = 0.0;
            current_rate += WellGroupHelpersType::sumWellResRates(group,
                                                                  wellModel_.schedule(),
                                                                  well_state,
                                                                  reportStepIdx,
                                                                  phasePos,
                                                                  /*isInjector*/true,
                                                                  wellModel_.summaryState());
            // sum over all nodes
            current_rate = wellModel_.comm().sum(current_rate);

            const auto& controls = group.injectionControls(phase, wellModel_.summaryState());
            Scalar target = controls.resv_max_rate;

            if (group.has_gpmaint_control(phase, Group::InjectionCMode::RESV))
                target = wellModel_.groupState().gpmaint_target(group.name());

            if (target < current_rate) {
                Scalar scale = 1.0;
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
            Scalar production_Rate = 0.0;
            const auto& controls = group.injectionControls(phase, wellModel_.summaryState());
            const Group& groupRein = wellModel_.schedule().getGroup(controls.reinj_group, reportStepIdx);
            production_Rate += WellGroupHelpersType::sumWellSurfaceRates(groupRein,
                                                                         wellModel_.schedule(),
                                                                         well_state,
                                                                         reportStepIdx,
                                                                         phasePos,
                                                                         /*isInjector*/false,
                                                                         wellModel_.summaryState());

            // sum over all nodes
            production_Rate = wellModel_.comm().sum(production_Rate);

            Scalar current_rate = 0.0;
            current_rate += WellGroupHelpersType::sumWellSurfaceRates(group,
                                                                      wellModel_.schedule(),
                                                                      well_state,
                                                                      reportStepIdx,
                                                                      phasePos,
                                                                      /*isInjector*/true,
                                                                      wellModel_.summaryState());

            // sum over all nodes
            current_rate = wellModel_.comm().sum(current_rate);

            if (controls.target_reinj_fraction*production_Rate < current_rate) {
                Scalar scale = 1.0;
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
            Scalar voidage_rate = 0.0;
            const auto& controls = group.injectionControls(phase, wellModel_.summaryState());
            const Group& groupVoidage = wellModel_.schedule().getGroup(controls.voidage_group, reportStepIdx);
            voidage_rate += WellGroupHelpersType::sumWellResRates(groupVoidage,
                                                                  wellModel_.schedule(),
                                                                  well_state,
                                                                  reportStepIdx,
                                                                  pu.canonicalToActivePhaseIdx(waterPhaseIdx),
                                                                  false,
                                                                  wellModel_.summaryState());
            voidage_rate += WellGroupHelpersType::sumWellResRates(groupVoidage,
                                                                  wellModel_.schedule(),
                                                                  well_state,
                                                                  reportStepIdx,
                                                                  pu.canonicalToActivePhaseIdx(oilPhaseIdx),
                                                                  false,
                                                                  wellModel_.summaryState());
            voidage_rate += WellGroupHelpersType::sumWellResRates(groupVoidage,
                                                                  wellModel_.schedule(),
                                                                  well_state, reportStepIdx,
                                                                  pu.canonicalToActivePhaseIdx(gasPhaseIdx),
                                                                  false,
                                                                  wellModel_.summaryState());

            // sum over all nodes
            voidage_rate = wellModel_.comm().sum(voidage_rate);

            Scalar total_rate = 0.0;
            total_rate += WellGroupHelpersType::sumWellResRates(group,
                                                                wellModel_.schedule(),
                                                                well_state,
                                                                reportStepIdx,
                                                                pu.canonicalToActivePhaseIdx(waterPhaseIdx),
                                                                true,
                                                                wellModel_.summaryState());
            total_rate += WellGroupHelpersType::sumWellResRates(group,
                                                                wellModel_.schedule(),
                                                                well_state,
                                                                reportStepIdx,
                                                                pu.canonicalToActivePhaseIdx(oilPhaseIdx),
                                                                true,
                                                                wellModel_.summaryState());
            total_rate += WellGroupHelpersType::sumWellResRates(group,
                                                                wellModel_.schedule(),
                                                                well_state,
                                                                reportStepIdx,
                                                                pu.canonicalToActivePhaseIdx(gasPhaseIdx),
                                                                true,
                                                                wellModel_.summaryState());

            // sum over all nodes
            total_rate = wellModel_.comm().sum(total_rate);

            if (controls.target_void_fraction*voidage_rate < total_rate) {
                Scalar scale = 1.0;
                if (total_rate > 1e-12)
                    scale = controls.target_void_fraction*voidage_rate / total_rate;
                return std::make_pair(Group::InjectionCMode::VREP, scale);
            }
        }
    }
    return std::make_pair(Group::InjectionCMode::NONE, 1.0);
}

template<typename Scalar, typename IndexTraits>
std::pair<Group::ProductionCMode, Scalar>
BlackoilWellModelConstraints<Scalar, IndexTraits>::
checkGroupProductionConstraints(const Group& group,
                                const int reportStepIdx,
                                DeferredLogger& deferred_logger) const
{
    const auto& well_state = wellModel_.wellState();

    const auto controls = group.productionControls(wellModel_.summaryState());
    const Group::ProductionCMode& currentControl = wellModel_.groupState().production_control(group.name());
    const auto& pu = wellModel_.phaseUsage();
    if (group.has_control(Group::ProductionCMode::ORAT))
    {
        if (currentControl != Group::ProductionCMode::ORAT)
        {
            Scalar current_rate = 0.0;
            current_rate += WellGroupHelpersType::sumWellSurfaceRates(group,
                                                                      wellModel_.schedule(),
                                                                      well_state,
                                                                      reportStepIdx,
                                                                      pu.canonicalToActivePhaseIdx(oilPhaseIdx),
                                                                      false,
                                                                      wellModel_.summaryState());

            // sum over all nodes
            current_rate = wellModel_.comm().sum(current_rate);

            if (controls.oil_target < current_rate  ) {
                Scalar scale = 1.0;
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
            Scalar current_rate = 0.0;
            current_rate += WellGroupHelpersType::sumWellSurfaceRates(group,
                                                                      wellModel_.schedule(),
                                                                      well_state,
                                                                      reportStepIdx,
                                                                      pu.canonicalToActivePhaseIdx(waterPhaseIdx),
                                                                      false,
                                                                      wellModel_.summaryState());

            // sum over all nodes
            current_rate = wellModel_.comm().sum(current_rate);

            if (controls.water_target < current_rate  ) {
                Scalar scale = 1.0;
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
            Scalar current_rate = 0.0;
            current_rate += WellGroupHelpersType::sumWellSurfaceRates(group,
                                                                      wellModel_.schedule(),
                                                                      well_state,
                                                                      reportStepIdx,
                                                                      pu.canonicalToActivePhaseIdx(gasPhaseIdx),
                                                                      false,
                                                                      wellModel_.summaryState());

            // sum over all nodes
            current_rate = wellModel_.comm().sum(current_rate);
            if (controls.gas_target < current_rate  ) {
                Scalar scale = 1.0;
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
            Scalar current_rate = 0.0;
            current_rate += WellGroupHelpersType::sumWellSurfaceRates(group,
                                                                      wellModel_.schedule(),
                                                                      well_state,
                                                                      reportStepIdx,
                                                                      pu.canonicalToActivePhaseIdx(oilPhaseIdx),
                                                                      false,
                                                                      wellModel_.summaryState());
            current_rate += WellGroupHelpersType::sumWellSurfaceRates(group,
                                                                      wellModel_.schedule(),
                                                                      well_state,
                                                                      reportStepIdx,
                                                                      pu.canonicalToActivePhaseIdx(waterPhaseIdx),
                                                                      false,
                                                                      wellModel_.summaryState());

            // sum over all nodes
            current_rate = wellModel_.comm().sum(current_rate);

            bool skip = false;
            if (controls.liquid_target == controls.oil_target) {
                Scalar current_water_rate =
                    WellGroupHelpersType::sumWellSurfaceRates(group,
                                                              wellModel_.schedule(),
                                                              well_state,
                                                              reportStepIdx,
                                                              pu.canonicalToActivePhaseIdx(waterPhaseIdx),
                                                              false,
                                                              wellModel_.summaryState());
                current_water_rate = wellModel_.comm().sum(current_water_rate);
                if (std::abs(current_water_rate) < 1e-12) {
                    skip = true;
                    deferred_logger.debug("LRAT_ORAT_GROUP", "GROUP " + group.name() + " The LRAT target is equal the ORAT target and the water rate is zero, skip checking LRAT");
                }
            }

            if (!skip && controls.liquid_target < current_rate ) {
                Scalar scale = 1.0;
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
            Scalar current_rate = 0.0;
            current_rate +=
                WellGroupHelpersType::sumWellResRates(group,
                                                      wellModel_.schedule(),
                                                      well_state,
                                                      reportStepIdx,
                                                      pu.canonicalToActivePhaseIdx(waterPhaseIdx),
                                                      false,
                                                      wellModel_.summaryState());
            current_rate +=
                WellGroupHelpersType::sumWellResRates(group,
                                                      wellModel_.schedule(),
                                                      well_state,
                                                      reportStepIdx,
                                                      pu.canonicalToActivePhaseIdx(oilPhaseIdx),
                                                      false,
                                                      wellModel_.summaryState());
            current_rate +=
                WellGroupHelpersType::sumWellResRates(group,
                                                      wellModel_.schedule(),
                                                      well_state,
                                                      reportStepIdx,
                                                      pu.canonicalToActivePhaseIdx(gasPhaseIdx),
                                                      false,
                                                      wellModel_.summaryState());

            // sum over all nodes
            current_rate = wellModel_.comm().sum(current_rate);

            Scalar target = controls.resv_target;
            if (group.has_gpmaint_control(Group::ProductionCMode::RESV))
                target = wellModel_.groupState().gpmaint_target(group.name());

            if ( target < current_rate ) {
                Scalar scale = 1.0;
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
    return std::make_pair(Group::ProductionCMode::NONE, Scalar(1.0));
}

template<typename Scalar, typename IndexTraits>
bool BlackoilWellModelConstraints<Scalar, IndexTraits>::
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

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelConstraints<Scalar, IndexTraits>::
actionOnBrokenConstraints(const Group& group,
                          const Group::InjectionCMode& newControl,
                          const Phase& controlPhase,
                          GroupState<Scalar>& group_state,
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

template<typename Scalar, typename IndexTraits>
bool BlackoilWellModelConstraints<Scalar, IndexTraits>::
actionOnBrokenConstraints(const Group& group,
                          const int reportStepIdx,
                          const Group::GroupLimitAction group_limit_action,
                          const Group::ProductionCMode& newControl,
                          const WellState<Scalar, IndexTraits>& well_state,
                          std::optional<std::string>& worst_offending_well,
                          GroupState<Scalar>& group_state,
                          DeferredLogger& deferred_logger) const
{

    const Group::ProductionCMode oldControl =
        wellModel_.groupState().production_control(group.name());

    // We switch to higher groups independently of the given group limit action in GCONPROD item 7
    if (newControl == Group::ProductionCMode::FLD && oldControl != Group::ProductionCMode::FLD) {
        // If newControl is FLD, the group should be subject to higher order controls
        assert(group.productionGroupControlAvailable());
        group_state.production_control(group.name(), newControl);
        if (wellModel_.comm().rank() == 0) {
            const std::string message = fmt::format("Switching production control mode for group {} from {} to {}",
            group.name(),
            Group::ProductionCMode2String(oldControl),
            Group::ProductionCMode2String(newControl));
            deferred_logger.debug(message);
        }
        return true;
    }

    const auto action = [&group_limit_action](Group::ProductionCMode control) {
        switch (control) {
            case Group::ProductionCMode::ORAT: return group_limit_action.oil;
            case Group::ProductionCMode::WRAT: return group_limit_action.water;
            case Group::ProductionCMode::GRAT: return group_limit_action.gas;
            case Group::ProductionCMode::LRAT: return group_limit_action.liquid;
            default: return Group::ExceedAction::RATE;
        };
    };

    bool changed = false;
    std::string ss;
    switch (action(newControl)) {
    case Group::ExceedAction::NONE: {
        // In GroupKeywordHandlers::handleGCONPROD(), if an action is NONE
        // we do not add the associated control (ORAT, WRAT, GRAT or LRAT) to
        // the set of production controls, so Group::has_control() will return
        // false for such controls.
        // Therefore we should never be here!
        // Nevertheless, to guard somewhat against future changes, we will keep the
        // message here.
        ss = fmt::format("Procedure on exceeding {} limit is NONE for group {}. "
                         "Nothing is done.",
                         Group::ProductionCMode2String(newControl),
                         group.name());
        break;
    }
    case Group::ExceedAction::CON: {
        OPM_DEFLOG_THROW(std::runtime_error,
                         fmt::format("Group {} GroupProductionExceedLimit CON not implemented",
                                     group.name()),
                         deferred_logger);
        break;
    }
    case Group::ExceedAction::CON_PLUS: {
        OPM_DEFLOG_THROW(std::runtime_error,
                         fmt::format("Group {} GroupProductionExceedLimit CON_PLUS not implemented",
                                     group.name()),
                         deferred_logger);
        break;
    }
    case Group::ExceedAction::WELL: {
        std::tie(worst_offending_well, std::ignore) =
            WellGroupHelpersType::worstOffendingWell(group, wellModel_.schedule(), reportStepIdx,
                                                     newControl,
                                                     wellModel_.comm(), well_state, deferred_logger);
        break;
    }
    case Group::ExceedAction::PLUG: {
        OPM_DEFLOG_THROW(std::runtime_error,
                         fmt::format("Group {} GroupProductionExceedLimit PLUG not implemented",
                                     group.name()),
                         deferred_logger);
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
        changed = true;
        break;
    }
    default:
        OPM_THROW(std::runtime_error,
                  "Invalid procedure for maximum rate limit selected for group" + group.name());
    }

    if (!ss.empty() && wellModel_.comm().rank() == 0)
        deferred_logger.debug(ss);

    return changed;
}

template<typename Scalar, typename IndexTraits>
bool BlackoilWellModelConstraints<Scalar, IndexTraits>::
updateGroupIndividualControl(const Group& group,
                             const int reportStepIdx,
                             const int max_number_of_group_switch,
                             const bool update_group_switching_log,
                             std::map<std::string, std::array<std::vector<Group::InjectionCMode>, 3>>& switched_inj,
                             std::map<std::string, std::vector<Group::ProductionCMode>>& switched_prod,
                             std::map<std::string, std::pair<std::string, std::string>>& closed_offending_wells,
                             GroupState<Scalar>& group_state,
                             WellState<Scalar, IndexTraits>& well_state,
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
            bool group_is_oscillating = false;
            const auto currentControl = group_state.injection_control(group.name(), phase);
            if (auto groupPos = switched_inj.find(group.name()); groupPos != switched_inj.end()) {
                auto& ctrls = groupPos->second[static_cast<std::underlying_type_t<Phase>>(phase)];
                const int number_of_switches = std::count(ctrls.begin(), ctrls.end(), currentControl);
                group_is_oscillating = (number_of_switches >= max_number_of_group_switch);
                if (group_is_oscillating) {
                    const bool output_first_time = (number_of_switches == max_number_of_group_switch);
                    if (output_first_time) {
                        if (wellModel_.comm().rank() == 0 ) {
                            std::ostringstream os;
                            os << phase;
                            const std::string msg =
                                fmt::format("Group control for {} injector group {} is oscillating. Group control kept at {}.",
                                            std::move(os).str(),
                                            group.name(),
                                            Group::InjectionCMode2String(currentControl));
                            deferred_logger.info(msg);
                        }
                        ctrls.push_back(currentControl);
                    }
                }
            }

            if (group_is_oscillating) {
                continue;
            }

            const auto& changed_this = this->checkGroupInjectionConstraints(group,
                                                                            reportStepIdx,
                                                                            phase);
            if (changed_this.first != Group::InjectionCMode::NONE)
            {
                auto& group_log = switched_inj[group.name()][static_cast<std::underlying_type_t<Phase>>(phase)];
                if (update_group_switching_log || group_log.empty()) {
                    group_log.push_back(currentControl);
                }
                this->actionOnBrokenConstraints(group, changed_this.first, phase,
                                                group_state, deferred_logger);
                WellGroupHelpersType::updateWellRatesFromGroupTargetScale(changed_this.second,
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

        const Group::ProductionCMode currentControl = group_state.production_control(group.name());
        if (auto groupPos = switched_prod.find(group.name()); groupPos != switched_prod.end()) {
            auto& ctrls = groupPos->second;
            const int number_of_switches = std::count(ctrls.begin(), ctrls.end(), currentControl);
            const bool group_is_oscillating = (number_of_switches >= max_number_of_group_switch);
            if (group_is_oscillating) {
                const bool output_first_time = (number_of_switches == max_number_of_group_switch);
                if (output_first_time) {
                    if (wellModel_.comm().rank() == 0) {
                        const std::string msg =
                        fmt::format("Group control for production group {} is oscillating. Group control kept at {}.",
                                    group.name(),
                                    Group::ProductionCMode2String(currentControl));
                        deferred_logger.info(msg);
                    }
                    ctrls.push_back(currentControl);
                }
                return false;
            }
        }

        const auto& changed_this = this->checkGroupProductionConstraints(group,
                                                                         reportStepIdx,
                                                                         deferred_logger);
        const auto controls = group.productionControls(wellModel_.summaryState());

        if (changed_this.first != Group::ProductionCMode::NONE)
        {
            std::optional<std::string> worst_offending_well = std::nullopt;
            changed = this->actionOnBrokenConstraints(group, reportStepIdx,              
                                            controls.group_limit_action,
                                            changed_this.first, well_state, 
                                            worst_offending_well,
                                            group_state, deferred_logger);

            if(changed) {
                if (update_group_switching_log || switched_prod[group.name()].empty()) {
                    switched_prod[group.name()].push_back(currentControl);
                }

                WellGroupHelpersType::updateWellRatesFromGroupTargetScale(changed_this.second,
                                                                              group,
                                                                              wellModel_.schedule(),
                                                                              reportStepIdx,
                                                                              /* isInjector */ false,
                                                                              wellModel_.groupState(),
                                                                              well_state);
            } else if (worst_offending_well) {
                closed_offending_wells.insert_or_assign(group.name(), 
                            std::make_pair(Group::ProductionCMode2String(changed_this.first), *worst_offending_well));
            }
        }
    }

    return changed;
}

template class BlackoilWellModelConstraints<double, BlackOilDefaultFluidSystemIndices>;

#if FLOW_INSTANTIATE_FLOAT
template class BlackoilWellModelConstraints<float, BlackOilDefaultFluidSystemIndices>;
#endif

}
