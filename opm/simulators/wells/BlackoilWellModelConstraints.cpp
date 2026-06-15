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
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>

#include <fmt/format.h>

#include <stdexcept>

namespace {

template<typename IndexTraits>
int activePhasePos(const Opm::Phase phase,
                   const Opm::PhaseUsageInfo<IndexTraits>& pu)
{
    if (phase == Opm::Phase::GAS && pu.phaseIsActive(IndexTraits::gasPhaseIdx)) {
        return pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx);
    }
    else if (phase == Opm::Phase::OIL && pu.phaseIsActive(IndexTraits::oilPhaseIdx)) {
        return pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
    }
    else if (phase == Opm::Phase::WATER && pu.phaseIsActive(IndexTraits::waterPhaseIdx)) {
        return pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx);
    }

    OPM_THROW(std::runtime_error, "Unknown phase");
}

template<class Scalar, class IndexTraits>
std::optional<std::pair<Opm::Group::InjectionCMode, Scalar>>
checkRateConstraint(const Opm::Group& group,
                    const Opm::Phase& phase,
                    const Opm::Group::InjectionCMode currentControl,
                    const Opm::BlackoilWellModelGeneric<Scalar,IndexTraits>& wellModel)
{
    if (group.has_control(phase, Opm::Group::InjectionCMode::RATE) &&
        currentControl != Opm::Group::InjectionCMode::RATE)
    {
        const int phasePos = activePhasePos(phase, wellModel.phaseUsage());
        Scalar current_rate =
            wellModel.groupStateHelper().sumWellSurfaceRates(group, phasePos, /*isInjector*/true);

        // sum over all nodes
        current_rate = wellModel.comm().sum(current_rate);

        const auto& controls = group.injectionControls(phase, wellModel.summaryState());
        const Scalar target = group.has_gpmaint_control(phase, Opm::Group::InjectionCMode::RATE)
            ? wellModel.groupState().gpmaint_target(group.name())
            : controls.surface_max_rate;

        if (target < current_rate) {
            const Scalar scale = current_rate > 1e-12
                ? target / current_rate
                : 1.0;
            return std::make_pair(Opm::Group::InjectionCMode::RATE, scale);
        }
    }

    return std::nullopt;
}

template<typename Scalar, typename IndexTraits>
std::optional<std::pair<Opm::Group::InjectionCMode, Scalar>>
checkResvConstraint(const Opm::Group& group,
                    const Opm::Phase& phase,
                    const Opm::Group::InjectionCMode currentControl,
                    const Opm::BlackoilWellModelGeneric<Scalar,IndexTraits>& wellModel)
{
    if (group.has_control(phase, Opm::Group::InjectionCMode::RESV) &&
        currentControl != Opm::Group::InjectionCMode::RESV)
    {
        const int phasePos = activePhasePos(phase, wellModel.phaseUsage());
        Scalar current_rate =
            wellModel.groupStateHelper().sumWellResRates(group, phasePos, /*is_injector=*/true);
        // sum over all nodes
        current_rate = wellModel.comm().sum(current_rate);

        const auto& controls = group.injectionControls(phase, wellModel.summaryState());
        const Scalar target = group.has_gpmaint_control(phase, Opm::Group::InjectionCMode::RESV)
            ? wellModel.groupState().gpmaint_target(group.name())
            : controls.resv_max_rate;

        if (target < current_rate) {
            const Scalar scale = current_rate > 1e-12
                ? target / current_rate
                : 1.0;
            return std::make_pair(Opm::Group::InjectionCMode::RESV, scale);
        }
    }

    return std::nullopt;
}

template<typename Scalar, typename IndexTraits>
std::optional<std::pair<Opm::Group::InjectionCMode, Scalar>>
checkReinConstraint(const Opm::Group& group,
                    const Opm::Phase& phase,
                    const Opm::Group::InjectionCMode currentControl,
                    const int reportStepIdx,
                    const Opm::BlackoilWellModelGeneric<Scalar,IndexTraits>& wellModel)
{
    if (group.has_control(phase, Opm::Group::InjectionCMode::REIN) &&
        currentControl != Opm::Group::InjectionCMode::REIN)
    {
        const auto& controls = group.injectionControls(phase, wellModel.summaryState());
        const Opm::Group& groupRein = wellModel.schedule().getGroup(controls.reinj_group, reportStepIdx);
        const int phasePos = activePhasePos(phase, wellModel.phaseUsage());
        const auto& groupStateHelper = wellModel.groupStateHelper();
        Scalar production_Rate =
            groupStateHelper.sumWellSurfaceRates(groupRein, phasePos, /*is_injector=*/false);

        // sum over all nodes
        production_Rate = wellModel.comm().sum(production_Rate);

        Scalar current_rate = 0.0;
        current_rate += groupStateHelper.sumWellSurfaceRates(group, phasePos, /*is_injector=*/true);

        // sum over all nodes
        current_rate = wellModel.comm().sum(current_rate);

        if (controls.target_reinj_fraction * production_Rate < current_rate) {
            const Scalar scale = current_rate > 1e-12
                ? controls.target_reinj_fraction * production_Rate / current_rate
                : 1.0;
            return std::make_pair(Opm::Group::InjectionCMode::REIN, scale);
        }
    }

    return std::nullopt;
}

template<typename Scalar, typename IndexTraits>
std::optional<std::pair<Opm::Group::InjectionCMode, Scalar>>
checkVrepConstraint(const Opm::Group& group,
                    const Opm::Phase& phase,
                    const Opm::Group::InjectionCMode currentControl,
                    const int reportStepIdx,
                    const Opm::BlackoilWellModelGeneric<Scalar,IndexTraits>& wellModel)
{
    if (group.has_control(phase, Opm::Group::InjectionCMode::VREP) &&
        currentControl != Opm::Group::InjectionCMode::VREP)
    {
        const auto& pu = wellModel.phaseUsage();
        const auto& controls = group.injectionControls(phase, wellModel.summaryState());
        const Opm::Group& groupVoidage = wellModel.schedule().getGroup(controls.voidage_group, reportStepIdx);
        const auto& groupStateHelper = wellModel.groupStateHelper();
        Scalar voidage_rate =
            groupStateHelper.sumWellResRates(groupVoidage,
                                             pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx),
                                             /*is_injector=*/false);
        voidage_rate += groupStateHelper.sumWellResRates(groupVoidage,
                                                 pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx),
                                                 /*is_injector=*/false);
        voidage_rate += groupStateHelper.sumWellResRates(groupVoidage,
                                                 pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx),
                                                 /*is_injector=*/false);

        // sum over all nodes
        voidage_rate = wellModel.comm().sum(voidage_rate);

        Scalar total_rate =
            groupStateHelper.sumWellResRates(group,
                                             pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx),
                                             /*is_injector=*/true);
        total_rate += groupStateHelper.sumWellResRates(group,
                                               pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx),
                                               /*is_injector=*/true);
        total_rate += groupStateHelper.sumWellResRates(group,
                                               pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx),
                                               /*is_injector=*/true);

        // sum over all nodes
        total_rate = wellModel.comm().sum(total_rate);

        if (controls.target_void_fraction * voidage_rate < total_rate) {
            const Scalar scale = total_rate > 1e-12
                ? controls.target_void_fraction * voidage_rate / total_rate
                : 1.0;
            return std::make_pair(Opm::Group::InjectionCMode::VREP, scale);
        }
    }

    return std::nullopt;
}

}

namespace Opm {

template<typename Scalar, typename IndexTraits>
std::pair<Group::InjectionCMode, Scalar>
BlackoilWellModelConstraints<Scalar, IndexTraits>::
checkGroupInjectionConstraints(const Group& group,
                               const int reportStepIdx,
                               const Phase& phase) const
{
    const auto currentControl = wellModel_.groupState().injection_control(group.name(), phase);

    if (auto c = checkRateConstraint(group, phase, currentControl, wellModel_); c.has_value()) {
        return *c;
    }
    if (auto c = checkResvConstraint(group, phase, currentControl, wellModel_); c.has_value()) {
        return *c;
    }
    if (auto c = checkReinConstraint(group, phase, currentControl, reportStepIdx, wellModel_); c.has_value()) {
        return *c;
    }
    if (auto c = checkVrepConstraint(group, phase, currentControl, reportStepIdx, wellModel_); c.has_value()) {
        return *c;
    }

    return std::make_pair(Group::InjectionCMode::NONE, 1.0);
}

template<typename Scalar, typename IndexTraits>
std::pair<Group::ProductionCMode, Scalar>
BlackoilWellModelConstraints<Scalar, IndexTraits>::
checkGroupProductionConstraints(const Group& group) const
{
    // NOTE: For reservoir coupling slave groups with no GCONPROD in the slave deck,
    // checkGroupProductionConstraints() below returns NONE because has_control() is
    // false for all modes.  So actionOnBrokenConstraints() will never be called,
    // this is correct since the group only has one mode imposed by the master.
    return groupStateHelper().checkGroupProductionConstraints(group);
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
                          const Group::GroupLimitAction group_limit_action,
                          const Group::ProductionCMode& newControl,
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
            groupStateHelper().worstOffendingWell(group, newControl);
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
updateInjectionGroupControl(const Group& group,
                            const int reportStepIdx,
                            const int max_number_of_group_switch,
                            const bool update_group_switching_log,
                            std::map<std::string, std::array<std::vector<Group::InjectionCMode>, 3>>& switched_inj,
                            GroupState<Scalar>& group_state,
                            WellState<Scalar, IndexTraits>& well_state,
                            DeferredLogger& deferred_logger) const
{
    bool changed = false;
    for (const Phase phase : {Phase::WATER, Phase::OIL, Phase::GAS}) {
        if (!group.hasInjectionControl(phase)) {
            continue;
        }
        bool group_is_oscillating = false;
        const auto currentControl = group_state.injection_control(group.name(), phase);
        if (auto groupPos = switched_inj.find(group.name()); groupPos != switched_inj.end()) {
            auto& ctrls = groupPos->second[static_cast<std::underlying_type_t<Phase>>(phase)];
            const int number_of_switches = std::ranges::count(ctrls, currentControl);
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
        if (changed_this.first != Group::InjectionCMode::NONE) {
            auto& group_log = switched_inj[group.name()][static_cast<std::underlying_type_t<Phase>>(phase)];
            if (update_group_switching_log || group_log.empty()) {
                group_log.push_back(currentControl);
            }
            this->actionOnBrokenConstraints(group, changed_this.first, phase,
                                            group_state, deferred_logger);
            groupStateHelper().updateWellRatesFromGroupTargetScale(changed_this.second,
                                                           group,
                                                           /*is_injector=*/false,
                                                           well_state);
            changed = true;
        }
    }

    return changed;
}


template<typename Scalar, typename IndexTraits>
bool BlackoilWellModelConstraints<Scalar, IndexTraits>::
updateProductionGroupControl(const Group& group,
                             const int max_number_of_group_switch,
                             const bool update_group_switching_log,
                             std::map<std::string, std::vector<Group::ProductionCMode>>& switched_prod,
                             std::map<std::string, std::pair<std::string, std::string>>& closed_offending_wells,
                             GroupState<Scalar>& group_state,
                             WellState<Scalar, IndexTraits>& well_state,
                             DeferredLogger& deferred_logger) const
{
    bool changed = false;
    const Group::ProductionCMode currentControl = group_state.production_control(group.name());
    if (auto groupPos = switched_prod.find(group.name()); groupPos != switched_prod.end()) {
        auto& ctrls = groupPos->second;
        const int number_of_switches = std::ranges::count(ctrls, currentControl);
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

    const auto& changed_this = this->checkGroupProductionConstraints(group);
    const auto controls = group.productionControls(wellModel_.summaryState());

    if (changed_this.first != Group::ProductionCMode::NONE) {
        std::optional<std::string> worst_offending_well = std::nullopt;
        changed = this->actionOnBrokenConstraints(group,
                                                  controls.group_limit_action,
                                                  changed_this.first,
                                                  worst_offending_well,
                                                  group_state, deferred_logger);

        if (changed) {
            if (update_group_switching_log || switched_prod[group.name()].empty()) {
                switched_prod[group.name()].push_back(currentControl);
            }

            groupStateHelper().updateWellRatesFromGroupTargetScale(changed_this.second,
                                                           group,
                                                           /*is_injector=*/false,
                                                           well_state);
        }
        else if (worst_offending_well) {
            closed_offending_wells.insert_or_assign(group.name(),
                        std::make_pair(Group::ProductionCMode2String(changed_this.first), *worst_offending_well));
        }
    }

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
    if (group.isInjectionGroup()) {
        changed = updateInjectionGroupControl(group,
                                              reportStepIdx,
                                              max_number_of_group_switch,
                                              update_group_switching_log,
                                              switched_inj,
                                              group_state,
                                              well_state,
                                              deferred_logger);
    }

    if (group.isProductionGroup()) {
        changed = changed || updateProductionGroupControl(group,
                                                          max_number_of_group_switch,
                                                          update_group_switching_log,
                                                          switched_prod,
                                                          closed_offending_wells,
                                                          group_state,
                                                          well_state,
                                                          deferred_logger);
    }

    return changed;
}

template class BlackoilWellModelConstraints<double, BlackOilDefaultFluidSystemIndices>;

#if FLOW_INSTANTIATE_FLOAT
template class BlackoilWellModelConstraints<float, BlackOilDefaultFluidSystemIndices>;
#endif

}
