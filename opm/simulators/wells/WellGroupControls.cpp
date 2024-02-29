/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2018 IRIS

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
#include <opm/simulators/wells/WellGroupControls.hpp>

#include <opm/input/eclipse/EclipseState/Phase.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Group/GConSale.hpp>
#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <opm/input/eclipse/Schedule/ScheduleTypes.hpp>

#include <opm/material/densead/Evaluation.hpp>

#include <opm/simulators/wells/FractionCalculator.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include <opm/simulators/wells/TargetCalculator.hpp>
#include <opm/simulators/wells/WellGroupHelpers.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <algorithm>
#include <cstddef>
#include <cassert>

namespace Opm
{

template<class EvalWell>
void WellGroupControls::
getGroupInjectionControl(const Group& group,
                         const WellState& well_state,
                         const GroupState& group_state,
                         const Schedule& schedule,
                         const SummaryState& summaryState,
                         const InjectorType& injectorType,
                         const EvalWell& bhp,
                         const EvalWell& injection_rate,
                         const RateConvFunc& rateConverter,
                         double efficiencyFactor,
                         EvalWell& control_eq,
                         DeferredLogger& deferred_logger) const
{
    // Setting some defaults to silence warnings below.
    // Will be overwritten in the switch statement.
    Phase injectionPhase = Phase::WATER;
    switch (injectorType) {
    case InjectorType::WATER:
    {
        injectionPhase = Phase::WATER;
        break;
    }
    case InjectorType::OIL:
    {
        injectionPhase = Phase::OIL;
        break;
    }
    case InjectorType::GAS:
    {
        injectionPhase = Phase::GAS;
        break;
    }
    default:
        // Should not be here.
        assert(false);
    }

    auto currentGroupControl = group_state.injection_control(group.name(), injectionPhase);
    if (currentGroupControl == Group::InjectionCMode::FLD ||
        currentGroupControl == Group::InjectionCMode::NONE) {
        if (!group.injectionGroupControlAvailable(injectionPhase)) {
            // We cannot go any further up the hierarchy. This could
            // be the FIELD group, or any group for which this has
            // been set in GCONINJE or GCONPROD. If we are here
            // anyway, it is likely that the deck set inconsistent
            // requirements, such as GRUP control mode on a well with
            // no appropriate controls defined on any of its
            // containing groups. We will therefore use the wells' bhp
            // limit equation as a fallback.
            const auto& controls = well_.wellEcl().injectionControls(summaryState);
            control_eq = bhp - controls.bhp_limit;
            return;
        } else {
            // Inject share of parents control
            const auto& parent = schedule.getGroup(group.parent(), well_.currentStep());
            efficiencyFactor *= group.getGroupEfficiencyFactor();
            getGroupInjectionControl(parent, well_state,
                                     group_state, schedule,
                                     summaryState, injectorType,
                                     bhp, injection_rate,
                                     rateConverter,
                                     efficiencyFactor,
                                     control_eq,
                                     deferred_logger);
            return;
        }
    }

    const auto& well = well_.wellEcl();
    const auto pu = well_.phaseUsage();

    if (!group.isInjectionGroup()) {
        // use bhp as control eq and let the updateControl code find a valid control
        const auto& controls = well.injectionControls(summaryState);
        control_eq = bhp - controls.bhp_limit;
        return;
    }

    // If we are here, we are at the topmost group to be visited in the recursion.
    // This is the group containing the control we will check against.

    // Make conversion factors for RESV <-> surface rates.
    std::vector<double> resv_coeff(pu.num_phases, 1.0);
    rateConverter(0, well_.pvtRegionIdx(), std::nullopt, resv_coeff); // FIPNUM region 0 here, should use FIPNUM from WELSPECS.

    double sales_target = 0;
    if (schedule[well_.currentStep()].gconsale().has(group.name())) {
        const auto& gconsale = schedule[well_.currentStep()].gconsale().get(group.name(), summaryState);
        sales_target = gconsale.sales_target;
    }
    WellGroupHelpers::InjectionTargetCalculator tcalc(currentGroupControl, pu,
                                                      resv_coeff, group.name(),
                                                      sales_target, group_state,
                                                      injectionPhase,
                                                      group.has_gpmaint_control(injectionPhase, currentGroupControl),
                                                      deferred_logger);
    WellGroupHelpers::FractionCalculator fcalc(schedule, well_state,
                                               group_state, well_.currentStep(),
                                               well_.guideRate(),
                                               tcalc.guideTargetMode(),
                                               pu, false, injectionPhase);

    auto localFraction = [&](const std::string& child) {
        return fcalc.localFraction(child, child);
    };

    auto localReduction = [&](const std::string& group_name) {
        const std::vector<double>& groupTargetReductions = group_state.injection_reduction_rates(group_name);
        return tcalc.calcModeRateFromRates(groupTargetReductions);
    };

    std::optional<Group::InjectionControls> ctrl;
    if (!group.has_gpmaint_control(injectionPhase, currentGroupControl))
        ctrl = group.injectionControls(injectionPhase, summaryState);

    const double orig_target = tcalc.groupTarget(ctrl,
                                                deferred_logger);
    const auto chain = WellGroupHelpers::groupChainTopBot(well_.name(), group.name(),
                                                          schedule, well_.currentStep());
    // Because 'name' is the last of the elements, and not an ancestor, we subtract one below.
    const std::size_t num_ancestors = chain.size() - 1;
    double target = orig_target;
    for (std::size_t ii = 0; ii < num_ancestors; ++ii) {
        if ((ii == 0) || well_.guideRate()->has(chain[ii], injectionPhase)) {
            // Apply local reductions only at the control level
            // (top) and for levels where we have a specified
            // group guide rate.
            target -= localReduction(chain[ii]);
        }
        target *= localFraction(chain[ii+1]);
    }
    // Avoid negative target rates coming from too large local reductions.
    const double target_rate = std::max(0.0, target / efficiencyFactor);
    const auto current_rate = injection_rate;

    control_eq = current_rate - target_rate;
}

std::optional<double>
WellGroupControls::
getGroupInjectionTargetRate(const Group& group,
                            const WellState& well_state,
                            const GroupState& group_state,
                            const Schedule& schedule,
                            const SummaryState& summaryState,
                            const InjectorType& injectorType,
                            const RateConvFunc& rateConverter,
                            double efficiencyFactor,
                            DeferredLogger& deferred_logger) const
{
    // Setting some defaults to silence warnings below.
    // Will be overwritten in the switch statement.
    Phase injectionPhase = Phase::WATER;
    switch (injectorType) {
    case InjectorType::WATER:
    {
        injectionPhase = Phase::WATER;
        break;
    }
    case InjectorType::OIL:
    {
        injectionPhase = Phase::OIL;
        break;
    }
    case InjectorType::GAS:
    {
        injectionPhase = Phase::GAS;
        break;
    }
    default:
        // Should not be here.
        assert(false);
    }
    auto currentGroupControl = group_state.injection_control(group.name(), injectionPhase);
    if (currentGroupControl == Group::InjectionCMode::FLD ||
        currentGroupControl == Group::InjectionCMode::NONE) {
        if (!group.injectionGroupControlAvailable(injectionPhase)) {
            // We cannot go any further up the hierarchy. This could
            // be the FIELD group, or any group for which this has
            // been set in GCONINJE or GCONPROD. If we are here
            // anyway, it is likely that the deck set inconsistent
            // requirements, such as GRUP control mode on a well with
            // no appropriate controls defined on any of its
            // containing groups. We will therefore use the wells' bhp
            // limit equation as a fallback.
            return std::nullopt;
        } else {
            // Inject share of parents control
            const auto& parent = schedule.getGroup( group.parent(), well_.currentStep());
            efficiencyFactor *= group.getGroupEfficiencyFactor();
            return getGroupInjectionTargetRate(parent, well_state, group_state,
                                               schedule, summaryState, injectorType,
                                               rateConverter, efficiencyFactor, deferred_logger);
        }
    }

    const auto pu = well_.phaseUsage();

    if (!group.isInjectionGroup()) {
        return std::nullopt;
    }

    // If we are here, we are at the topmost group to be visited in the recursion.
    // This is the group containing the control we will check against.

    // Make conversion factors for RESV <-> surface rates.
    std::vector<double> resv_coeff(pu.num_phases, 1.0);
    rateConverter(0, well_.pvtRegionIdx(), std::nullopt, resv_coeff); // FIPNUM region 0 here, should use FIPNUM from WELSPECS.

    double sales_target = 0;
    if (schedule[well_.currentStep()].gconsale().has(group.name())) {
        const auto& gconsale = schedule[well_.currentStep()].gconsale().get(group.name(), summaryState);
        sales_target = gconsale.sales_target;
    }
    WellGroupHelpers::InjectionTargetCalculator tcalc(currentGroupControl, pu, resv_coeff,
                                                      group.name(), sales_target, group_state,
                                                      injectionPhase,
                                                      group.has_gpmaint_control(injectionPhase, currentGroupControl),
                                                      deferred_logger);
    WellGroupHelpers::FractionCalculator fcalc(schedule, well_state, group_state,
                                               well_.currentStep(), well_.guideRate(),
                                               tcalc.guideTargetMode(), pu, false, injectionPhase);

    auto localFraction = [&](const std::string& child) {
        return fcalc.localFraction(child, child); //Note child needs to be passed to always include since the global isGrup map is not updated yet.
    };

    auto localReduction = [&](const std::string& group_name) {
        const std::vector<double>& groupTargetReductions = group_state.injection_reduction_rates(group_name);
        return tcalc.calcModeRateFromRates(groupTargetReductions);
    };

    std::optional<Group::InjectionControls> ctrl;
    if (!group.has_gpmaint_control(injectionPhase, currentGroupControl))
        ctrl = group.injectionControls(injectionPhase, summaryState);

    const double orig_target = tcalc.groupTarget(ctrl, deferred_logger);

    const auto chain = WellGroupHelpers::groupChainTopBot(well_.name(), group.name(), schedule, well_.currentStep());
    // Because 'name' is the last of the elements, and not an ancestor, we subtract one below.
    const std::size_t num_ancestors = chain.size() - 1;
    double target = orig_target;
    for (std::size_t ii = 0; ii < num_ancestors; ++ii) {
        if ((ii == 0) || well_.guideRate()->has(chain[ii], injectionPhase)) {
            // Apply local reductions only at the control level
            // (top) and for levels where we have a specified
            // group guide rate.
            target -= localReduction(chain[ii]);
        }
        target *= localFraction(chain[ii+1]);
    }
    return std::max(0.0, target / efficiencyFactor);
}

template<class EvalWell>
void WellGroupControls::getGroupProductionControl(const Group& group,
                                                  const WellState& well_state,
                                                  const GroupState& group_state,
                                                  const Schedule& schedule,
                                                  const SummaryState& summaryState,
                                                  const EvalWell& bhp,
                                                  const std::vector<EvalWell>& rates,
                                                  const RateConvFunc& rateConverter,
                                                  double efficiencyFactor,
                                                  EvalWell& control_eq,
                                                  DeferredLogger& deferred_logger) const
{
    const Group::ProductionCMode& currentGroupControl = group_state.production_control(group.name());
    if (currentGroupControl == Group::ProductionCMode::FLD ||
        currentGroupControl == Group::ProductionCMode::NONE) {
        if (!group.productionGroupControlAvailable()) {
            // We cannot go any further up the hierarchy. This could
            // be the FIELD group, or any group for which this has
            // been set in GCONINJE or GCONPROD. If we are here
            // anyway, it is likely that the deck set inconsistent
            // requirements, such as GRUP control mode on a well with
            // no appropriate controls defined on any of its
            // containing groups. We will therefore use the wells' bhp
            // limit equation as a fallback.
            const auto& controls = well_.wellEcl().productionControls(summaryState);
            control_eq = bhp - controls.bhp_limit;
            return;
        } else {
            // Produce share of parents control
            const auto& parent = schedule.getGroup(group.parent(), well_.currentStep());
            efficiencyFactor *= group.getGroupEfficiencyFactor();
            getGroupProductionControl(parent, well_state, group_state,
                                      schedule, summaryState, bhp,
                                      rates, rateConverter,
                                      efficiencyFactor, control_eq, deferred_logger);
            return;
        }
    }

    const auto& well = well_.wellEcl();
    const auto pu = well_.phaseUsage();

    if (!group.isProductionGroup()) {
        // use bhp as control eq and let the updateControl code find a valid control
        const auto& controls = well.productionControls(summaryState);
        control_eq = bhp - controls.bhp_limit;
        return;
    }

    // If we are here, we are at the topmost group to be visited in the recursion.
    // This is the group containing the control we will check against.

    // Make conversion factors for RESV <-> surface rates.
    std::vector<double> resv_coeff(well_.phaseUsage().num_phases, 1.0);
    rateConverter(0, well_.pvtRegionIdx(), group.name(), resv_coeff); // FIPNUM region 0 here, should use FIPNUM from WELSPECS.

    // gconsale may adjust the grat target.
    // the adjusted rates is send to the targetCalculator
    double gratTargetFromSales = 0.0;
    if (group_state.has_grat_sales_target(group.name()))
        gratTargetFromSales = group_state.grat_sales_target(group.name());

    WellGroupHelpers::TargetCalculator tcalc(currentGroupControl, pu, resv_coeff,
                                             gratTargetFromSales, group.name(),
                                             group_state,
                                             group.has_gpmaint_control(currentGroupControl));
    WellGroupHelpers::FractionCalculator fcalc(schedule, well_state, group_state,
                                               well_.currentStep(),
                                               well_.guideRate(),
                                               tcalc.guideTargetMode(),
                                               pu, true, Phase::OIL);

    auto localFraction = [&](const std::string& child) {
        return fcalc.localFraction(child, child);
    };

    auto localReduction = [&](const std::string& group_name) {
        const std::vector<double>& groupTargetReductions = group_state.production_reduction_rates(group_name);
        return tcalc.calcModeRateFromRates(groupTargetReductions);
    };

    std::optional<Group::ProductionControls> ctrl;
    if (!group.has_gpmaint_control(currentGroupControl))
        ctrl = group.productionControls(summaryState);

    const double orig_target = tcalc.groupTarget(ctrl, deferred_logger);
    const auto chain = WellGroupHelpers::groupChainTopBot(well_.name(), group.name(),
                                                          schedule, well_.currentStep());
    // Because 'name' is the last of the elements, and not an ancestor, we subtract one below.
    const std::size_t num_ancestors = chain.size() - 1;
    double target = orig_target;
    for (std::size_t ii = 0; ii < num_ancestors; ++ii) {
        if ((ii == 0) || well_.guideRate()->has(chain[ii])) {
            // Apply local reductions only at the control level
            // (top) and for levels where we have a specified
            // group guide rate.
            target -= localReduction(chain[ii]);
        }
        target *= localFraction(chain[ii+1]);
    }
    // Avoid negative target rates coming from too large local reductions.
    const double target_rate = std::max(0.0, target / efficiencyFactor);
    const auto current_rate = -tcalc.calcModeRateFromRates(rates); // Switch sign since 'rates' are negative for producers.
    control_eq = current_rate - target_rate;
}

double WellGroupControls::
getGroupProductionTargetRate(const Group& group,
                             const WellState& well_state,
                             const GroupState& group_state,
                             const Schedule& schedule,
                             const SummaryState& summaryState,
                             const RateConvFunc& rateConverter,
                             double efficiencyFactor,
                             DeferredLogger& deferred_logger) const
{
    const Group::ProductionCMode& currentGroupControl = group_state.production_control(group.name());
    if (currentGroupControl == Group::ProductionCMode::FLD ||
        currentGroupControl == Group::ProductionCMode::NONE) {
        if (!group.productionGroupControlAvailable()) {
            return 1.0;
        } else {
            // Produce share of parents control
            const auto& parent = schedule.getGroup(group.parent(), well_.currentStep());
            efficiencyFactor *= group.getGroupEfficiencyFactor();
            return getGroupProductionTargetRate(parent, well_state, group_state,
                                                schedule, summaryState,
                                                rateConverter, efficiencyFactor,
                                                deferred_logger);
        }
    }

    const auto pu = well_.phaseUsage();

    if (!group.isProductionGroup()) {
        return 1.0;
    }

    // If we are here, we are at the topmost group to be visited in the recursion.
    // This is the group containing the control we will check against.

    // Make conversion factors for RESV <-> surface rates.
    std::vector<double> resv_coeff(well_.phaseUsage().num_phases, 1.0);
    rateConverter(0, well_.pvtRegionIdx(), group.name(), resv_coeff); // FIPNUM region 0 here, should use FIPNUM from WELSPECS.

    // gconsale may adjust the grat target.
    // the adjusted rates is send to the targetCalculator
    double gratTargetFromSales = 0.0;
    if (group_state.has_grat_sales_target(group.name()))
        gratTargetFromSales = group_state.grat_sales_target(group.name());

    WellGroupHelpers::TargetCalculator tcalc(currentGroupControl, pu, resv_coeff, gratTargetFromSales, group.name(), group_state, group.has_gpmaint_control(currentGroupControl));
    WellGroupHelpers::FractionCalculator fcalc(schedule, well_state, group_state,
                                               well_.currentStep(),
                                               well_.guideRate(),
                                               tcalc.guideTargetMode(),
                                               pu, true, Phase::OIL);

    auto localFraction = [&](const std::string& child) {
        return fcalc.localFraction(child, child); //Note child needs to be passed to always include since the global isGrup map is not updated yet.
    };

    auto localReduction = [&](const std::string& group_name) {
        const std::vector<double>& groupTargetReductions = group_state.production_reduction_rates(group_name);
        return tcalc.calcModeRateFromRates(groupTargetReductions);
    };

    std::optional<Group::ProductionControls> ctrl;
    if (!group.has_gpmaint_control(currentGroupControl))
        ctrl = group.productionControls(summaryState);

    const double orig_target = tcalc.groupTarget(ctrl, deferred_logger);
    const auto chain = WellGroupHelpers::groupChainTopBot(well_.name(), group.name(),
                                                          schedule, well_.currentStep());
    // Because 'name' is the last of the elements, and not an ancestor, we subtract one below.
    const std::size_t num_ancestors = chain.size() - 1;
    double target = orig_target;
    for (std::size_t ii = 0; ii < num_ancestors; ++ii) {
        if ((ii == 0) || well_.guideRate()->has(chain[ii])) {
            // Apply local reductions only at the control level
            // (top) and for levels where we have a specified
            // group guide rate.
            target -= localReduction(chain[ii]);
        }
        target *= localFraction(chain[ii+1]);
    }
    // Avoid negative target rates coming from too large local reductions.
    const double target_rate = std::max(0.0, target / efficiencyFactor);
    const auto& ws = well_state.well(well_.indexOfWell());
    const auto& rates = ws.surface_rates;
    const auto current_rate = -tcalc.calcModeRateFromRates(rates); // Switch sign since 'rates' are negative for producers.
    double scale = 1.0;
    if (target_rate == 0.0) {
        return 0.0;
    }

    if (current_rate > 1e-14)
        scale = target_rate/current_rate;
    return scale;
}

#define INSTANCE(...) \
template void WellGroupControls:: \
getGroupInjectionControl<__VA_ARGS__>(const Group&, \
                                      const WellState&, \
                                      const GroupState&, \
                                      const Schedule&, \
                                      const SummaryState&, \
                                      const InjectorType&, \
                                      const __VA_ARGS__& bhp, \
                                      const __VA_ARGS__& injection_rate, \
                                      const RateConvFunc& rateConverter, \
                                      double efficiencyFactor, \
                                      __VA_ARGS__& control_eq, \
                                      DeferredLogger& deferred_logger) const; \
template void WellGroupControls:: \
getGroupProductionControl<__VA_ARGS__>(const Group&, \
                                       const WellState&, \
                                       const GroupState&, \
                                       const Schedule&, \
                                       const SummaryState&, \
                                       const __VA_ARGS__& bhp, \
                                       const std::vector<__VA_ARGS__>&, \
                                       const RateConvFunc& rateConverter, \
                                       double efficiencyFactor, \
                                       __VA_ARGS__& control_eq, \
                                       DeferredLogger& deferred_logger) const; \

INSTANCE(DenseAd::Evaluation<double,3,0u>)
INSTANCE(DenseAd::Evaluation<double,4,0u>)
INSTANCE(DenseAd::Evaluation<double,5,0u>)
INSTANCE(DenseAd::Evaluation<double,6,0u>)
INSTANCE(DenseAd::Evaluation<double,7,0u>)
INSTANCE(DenseAd::Evaluation<double,8,0u>)
INSTANCE(DenseAd::Evaluation<double,9,0u>)
INSTANCE(DenseAd::Evaluation<double,10,0u>)
INSTANCE(DenseAd::Evaluation<double,-1,4u>)
INSTANCE(DenseAd::Evaluation<double,-1,5u>)
INSTANCE(DenseAd::Evaluation<double,-1,6u>)
INSTANCE(DenseAd::Evaluation<double,-1,7u>)
INSTANCE(DenseAd::Evaluation<double,-1,8u>)
INSTANCE(DenseAd::Evaluation<double,-1,9u>)
INSTANCE(DenseAd::Evaluation<double,-1,10u>)
INSTANCE(DenseAd::Evaluation<double,-1,11u>)
} // namespace Opm
