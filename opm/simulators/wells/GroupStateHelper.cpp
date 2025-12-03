/*
  Copyright 2025 Equinor ASA

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
#include <opm/simulators/wells/GroupStateHelper.hpp>

#include <opm/common/TimingMacros.hpp>
#include <opm/input/eclipse/Schedule/Group/GConSale.hpp>
#include <opm/input/eclipse/Schedule/Group/GroupSatelliteInjection.hpp>
#include <opm/input/eclipse/Schedule/Network/ExtNetwork.hpp>
#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/FractionCalculator.hpp>
#include <opm/simulators/wells/TargetCalculator.hpp>

#include <array>
#include <cstddef>
#include <stack>
#include <set>

namespace Opm
{

template <typename Scalar, typename IndexTraits>
GroupStateHelper<Scalar, IndexTraits>::GroupStateHelper(WellState<Scalar, IndexTraits>& well_state,
                                                      GroupState<Scalar>& group_state,
                                                      const Schedule& schedule,
                                                      const SummaryState& summary_state,
                                                      const GuideRate& guide_rate,
                                                      const PhaseUsageInfo<IndexTraits>& phase_usage_info)
    : well_state_ {&well_state}
    , group_state_ {&group_state}
    , schedule_ {schedule}
    , summary_state_ {summary_state}
    , guide_rate_ {guide_rate}
    , phase_usage_info_ {phase_usage_info}
{
}


// ---------------------------------------------------------------------
// Public methods
// ---------------------------------------------------------------------

template <typename Scalar, typename IndexTraits>
void
GroupStateHelper<Scalar, IndexTraits>::accumulateGroupEfficiencyFactor(const Group& group,
                                                                      Scalar& factor) const
{
    factor *= group.getGroupEfficiencyFactor();
    if (group.parent() != "FIELD" && !group.parent().empty())
        this->accumulateGroupEfficiencyFactor(this->schedule_.getGroup(group.parent(), this->report_step_),
                                              factor);
}

template <typename Scalar, typename IndexTraits>
std::pair<bool, Scalar>
GroupStateHelper<Scalar, IndexTraits>::checkGroupConstraintsInj(const std::string& name,
                                                               const std::string& parent,
                                                               const Group& group,
                                                               const Scalar* rates,
                                                               const Phase injection_phase,
                                                               const Scalar efficiency_factor,
                                                               const std::vector<Scalar>& resv_coeff,
                                                               const bool check_guide_rate,
                                                               DeferredLogger& deferred_logger) const
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

    auto current_group_control = this->groupState().injection_control(group.name(), injection_phase);

    if (current_group_control == Group::InjectionCMode::FLD
        || current_group_control == Group::InjectionCMode::NONE) {
        // Return if we are not available for parent group.
        if (!group.injectionGroupControlAvailable(injection_phase)) {
            return std::make_pair(false, Scalar(1));
        }
        // Otherwise: check production share of parent's control.
        const auto& parent_group = this->schedule_.getGroup(group.parent(), this->report_step_);
        return this->checkGroupConstraintsInj(name,
                                              parent,
                                              parent_group,
                                              rates,
                                              injection_phase,
                                              efficiency_factor * group.getGroupEfficiencyFactor(),
                                              resv_coeff,
                                              check_guide_rate,
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
    if (this->schedule_[this->report_step_].gconsale().has(group.name())) {
        const auto& gconsale
            = this->schedule_[this->report_step_].gconsale().get(group.name(), this->summary_state_);
        sales_target = gconsale.sales_target;
    }
    GroupStateHelpers::InjectionTargetCalculator<Scalar, IndexTraits> tcalc {
        current_group_control,
        this->phase_usage_info_,
        resv_coeff,
        group.name(),
        sales_target,
        this->groupState(),
        injection_phase,
        group.has_gpmaint_control(injection_phase, current_group_control),
        deferred_logger};

    GroupStateHelpers::FractionCalculator fcalc {this->schedule_,
                                                 *this,
                                                 this->summary_state_,
                                                 this->report_step_,
                                                 &this->guide_rate_,
                                                 tcalc.guideTargetMode(),
                                                 /*is_producer=*/false,
                                                 injection_phase};

    auto local_fraction_lambda = [&](const std::string& child) { return fcalc.localFraction(child, name); };

    auto local_reduction_lambda = [&](const std::string& group_name) {
        const std::vector<Scalar>& group_target_reductions
            = this->groupState().injection_reduction_rates(group_name);
        return tcalc.calcModeRateFromRates(group_target_reductions);
    };

    auto local_current_rate_lambda = [&](const std::string& group_name) {
        const std::vector<Scalar>& group_surface_rates
            = this->groupState().injection_surface_rates(group_name);
        return tcalc.calcModeRateFromRates(group_surface_rates);
    };

    const auto chain = this->groupChainTopBot(name, group.name());
    // Because 'name' is the last of the elements, and not an ancestor, we subtract one below.
    const std::size_t num_ancestors = chain.size() - 1;
    // we need to find out the level where the current well is applied to the local reduction
    std::size_t local_reduction_level = 0;
    for (std::size_t ii = 1; ii < num_ancestors; ++ii) {
        const int num_gr_ctrl = this->groupControlledWells(chain[ii],
                                                           /*always_included_child=*/"",
                                                           /*is_production_group=*/false,
                                                           injection_phase);
        if (this->guide_rate_.has(chain[ii], injection_phase) && num_gr_ctrl > 0) {
            local_reduction_level = ii;
        }
    }

    // check whether guide rate is violated
    if (check_guide_rate) {
        for (std::size_t ii = 1; ii < num_ancestors; ++ii) {
            if (this->guide_rate_.has(chain[ii], injection_phase)) {
                const auto& guided_group = chain[ii];
                const Scalar grefficiency
                    = this->schedule_.getGroup(guided_group, this->report_step_).getGroupEfficiencyFactor();
                const Scalar current_rate_fraction = grefficiency * local_current_rate_lambda(guided_group)
                    / local_current_rate_lambda(chain[ii - 1]);
                const Scalar guiderate_fraction = local_fraction_lambda(guided_group);
                // we add a factor here to avoid switching due to numerical instability
                const Scalar factor = 1.01;
                if (current_rate_fraction > (guiderate_fraction * factor)) {
                    return std::make_pair(true, guiderate_fraction / current_rate_fraction);
                }
            }
        }
    }

    if (this->schedule_.hasWell(name)
        && this->wellState().well(name).group_target) { // for wells we already have computed the target
        Scalar scale = 1.0;
        const auto& group_target = this->wellState().well(name).group_target;
        const Scalar group_target_rate_available = group_target->target_value;
        const Scalar current_well_rate_available
            = tcalc.calcModeRateFromRates(rates); // Switch sign since 'rates' are negative for producers.
        if (current_well_rate_available > 1e-12) {
            scale = group_target_rate_available / current_well_rate_available;
        }
        return std::make_pair(current_well_rate_available > group_target_rate_available, scale);
    }

    std::optional<Group::InjectionControls> ctrl;
    if (!group.has_gpmaint_control(injection_phase, current_group_control))
        ctrl = group.injectionControls(injection_phase, this->summary_state_);


    const Scalar orig_target = tcalc.groupTarget(ctrl, deferred_logger);
    // Assume we have a chain of groups as follows: BOTTOM -> MIDDLE -> TOP.
    // Then ...
    // TODO finish explanation.
    const Scalar current_rate_available
        = tcalc.calcModeRateFromRates(rates); // Switch sign since 'rates' are negative for producers.

    // Compute portion of target corresponding to current_rate_available
    Scalar target = orig_target;
    for (std::size_t ii = 0; ii < num_ancestors; ++ii) {
        if ((ii == 0) || this->guide_rate_.has(chain[ii], injection_phase)) {
            // Apply local reductions only at the control level
            // (top) and for levels where we have a specified
            // group guide rate.
            if (local_reduction_level >= ii) {
                target -= local_reduction_lambda(chain[ii]);
            }

            // Add my reduction back at the level where it is included in the local reduction
            if (local_reduction_level == ii) {
                target += current_rate_available * efficiency_factor;
            }
        }
        target *= local_fraction_lambda(chain[ii + 1]);
    }
    // Avoid negative target rates comming from too large local reductions.
    const Scalar target_rate_available = std::max(Scalar(1e-12), target / efficiency_factor);
    Scalar scale = 1.0;
    if (current_rate_available > 1e-12)
        scale = target_rate_available / current_rate_available;

    return std::make_pair(current_rate_available > target_rate_available, scale);
}

template <typename Scalar, typename IndexTraits>
std::pair<bool, Scalar>
GroupStateHelper<Scalar, IndexTraits>::checkGroupConstraintsProd(const std::string& name,
                                                                const std::string& parent,
                                                                const Group& group,
                                                                const Scalar* rates,
                                                                const Scalar efficiency_factor,
                                                                const std::vector<Scalar>& resv_coeff,
                                                                const bool check_guide_rate,
                                                                DeferredLogger& deferred_logger) const
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
    OPM_TIMEFUNCTION();
    const Group::ProductionCMode& current_group_control = this->groupState().production_control(group.name());

    if (current_group_control == Group::ProductionCMode::FLD
        || current_group_control == Group::ProductionCMode::NONE) {
        // Return if we are not available for parent group.
        if (!group.productionGroupControlAvailable()) {
            return std::make_pair(false, 1);
        }
        // Otherwise: check production share of parent's control.
        const auto& parent_group = this->schedule_.getGroup(group.parent(), this->report_step_);
        return this->checkGroupConstraintsProd(name,
                                               parent,
                                               parent_group,
                                               rates,
                                               efficiency_factor * group.getGroupEfficiencyFactor(),
                                               resv_coeff,
                                               check_guide_rate,
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
    Scalar grat_target_from_sales = 0.0;
    if (this->groupState().has_grat_sales_target(group.name()))
        grat_target_from_sales = this->groupState().grat_sales_target(group.name());

    GroupStateHelpers::TargetCalculator<Scalar, IndexTraits> tcalc {current_group_control,
                                                                    this->phase_usage_info_,
                                                                    resv_coeff,
                                                                    grat_target_from_sales,
                                                                    group.name(),
                                                                    this->groupState(),
                                                                    group.has_gpmaint_control(current_group_control)};

    GroupStateHelpers::FractionCalculator<Scalar, IndexTraits> fcalc {this->schedule_,
                                                                      *this,
                                                                      this->summary_state_,
                                                                      this->report_step_,
                                                                      &this->guide_rate_,
                                                                      tcalc.guideTargetMode(),
                                                                      /*is_producer=*/true,
                                                                      /*injection_phase=*/Phase::OIL};

    auto local_fraction_lambda = [&](const std::string& child) { return fcalc.localFraction(child, name); };

    auto local_reduction_lambda = [&](const std::string& group_name) {
        const std::vector<Scalar>& group_target_reductions
            = this->groupState().production_reduction_rates(group_name);
        return tcalc.calcModeRateFromRates(group_target_reductions);
    };

    auto local_current_rate_lambda = [&](const std::string& group_name) {
        const std::vector<Scalar>& group_surface_rates = this->groupState().production_rates(group_name);
        return tcalc.calcModeRateFromRates(group_surface_rates);
    };

    std::optional<Group::ProductionControls> ctrl;
    if (!group.has_gpmaint_control(current_group_control))
        ctrl = group.productionControls(this->summary_state_);

    const Scalar orig_target = tcalc.groupTarget(ctrl, deferred_logger);
    // Assume we have a chain of groups as follows: BOTTOM -> MIDDLE -> TOP.
    // Then ...
    // TODO finish explanation.
    const Scalar current_rate_available
        = -tcalc.calcModeRateFromRates(rates); // Switch sign since 'rates' are negative for producers.
    const auto chain = this->groupChainTopBot(name, group.name());
    // Because 'name' is the last of the elements, and not an ancestor, we subtract one below.
    const std::size_t num_ancestors = chain.size() - 1;

    // check whether guide rate is violated
    if (check_guide_rate) {
        for (std::size_t ii = 1; ii < num_ancestors; ++ii) {
            if (this->guide_rate_.has(chain[ii])) {
                const auto& guided_group = chain[ii];
                const Scalar grefficiency
                    = this->schedule_.getGroup(guided_group, this->report_step_).getGroupEfficiencyFactor();
                const Scalar current_rate_fraction = grefficiency * local_current_rate_lambda(guided_group)
                    / (local_current_rate_lambda(chain[ii - 1]));
                const Scalar guiderate_fraction = local_fraction_lambda(guided_group);
                // we add a factor here to avoid switching due to numerical instability
                const Scalar factor = 1.01;
                if (current_rate_fraction > (guiderate_fraction * factor)) {
                    return std::make_pair(true, guiderate_fraction / current_rate_fraction);
                }
            }
        }
    }

    if (this->schedule_.hasWell(name) && this->wellState().well(name).group_target) {
        // for wells we already have computed the target

        // Switch sign since 'rates' are negative for producers.
        const Scalar current_well_rate_available = -tcalc.calcModeRateFromRates(rates);
        const auto& group_target = this->wellState().well(name).group_target;
        const Scalar group_target_rate_available = group_target->target_value;
        Scalar scale = 1.0;
        if (current_well_rate_available > 1e-12) {
            scale = group_target_rate_available / current_well_rate_available;
        }

        return std::make_pair(current_well_rate_available > group_target_rate_available, scale);
    }

    // we need to find out the level where the current well is applied to the local reduction
    std::size_t local_reduction_level = 0;
    for (std::size_t ii = 1; ii < num_ancestors; ++ii) {
        const int num_gr_ctrl = this->groupControlledWells(chain[ii],
                                                           /*always_included_child=*/"",
                                                           /*is_production_group=*/true,
                                                           /*injection_phase=*/Phase::OIL/*not used*/);
        if (this->guide_rate_.has(chain[ii]) && num_gr_ctrl > 0) {
            local_reduction_level = ii;
        }
    }

    // Compute portion of target corresponding to current_rate_available
    Scalar target = orig_target;
    for (std::size_t ii = 0; ii < num_ancestors; ++ii) {
        if ((ii == 0) || this->guide_rate_.has(chain[ii])) {
            // Apply local reductions only at the control level
            // (top) and for levels where we have a specified
            // group guide rate.
            if (local_reduction_level >= ii) {
                target -= local_reduction_lambda(chain[ii]);
            }
            // Add my reduction back at the level where it is included in the local reduction
            if (local_reduction_level == ii) {
                target += current_rate_available * efficiency_factor;
            }
        }
        target *= local_fraction_lambda(chain[ii + 1]);
    }
    // Avoid negative target rates comming from too large local reductions.
    const Scalar target_rate_available = std::max(Scalar(1e-12), target / efficiency_factor);

    Scalar scale = 1.0;
    if (current_rate_available > 1e-12)
        scale = target_rate_available / current_rate_available;

    return std::make_pair(current_rate_available > target_rate_available, scale);
}

template <typename Scalar, typename IndexTraits>
Scalar
GroupStateHelper<Scalar, IndexTraits>::getGuideRate(const std::string& name,
                                                   const GuideRateModel::Target target) const
{
    if (this->schedule_.hasWell(name, this->report_step_)) {
        if (this->guide_rate_.has(name) || this->guide_rate_.hasPotentials(name)) {
            return this->guide_rate_.get(name, target, this->getWellRateVector(name));
        } else {
            return 0.0;
        }
    }

    if (this->guide_rate_.has(name)) {
        return this->guide_rate_.get(name, target, this->getProductionGroupRateVector(name));
    }

    Scalar total_guide_rate = 0.0;
    const Group& group = this->schedule_.getGroup(name, this->report_step_);

    for (const std::string& group_name : group.groups()) {
        const Group::ProductionCMode& current_group_control
            = this->groupState().production_control(group_name);
        if (current_group_control == Group::ProductionCMode::FLD
            || current_group_control == Group::ProductionCMode::NONE) {
            // accumulate from sub wells/groups
            total_guide_rate += this->getGuideRate(group_name, target);
        }
    }

    for (const std::string& well_name : group.wells()) {
        const auto& well_tmp = this->schedule_.getWell(well_name, this->report_step_);

        if (well_tmp.isInjector())
            continue;

        const auto well_index = this->wellState().index(well_name);
        if (!well_index.has_value())
            continue;

        const auto& ws = this->wellState().well(well_index.value());
        if (ws.status == Well::Status::SHUT)
            continue;

        if (!this->wellState().isProductionGrup(well_name))
            continue;

        // Only count wells under group control or the ru
        if (!this->wellState().isProductionGrup(well_name))
            continue;

        total_guide_rate += this->getGuideRate(well_name, target);
    }
    return total_guide_rate;
}

template <typename Scalar, typename IndexTraits>
GuideRate::RateVector
GroupStateHelper<Scalar, IndexTraits>::getProductionGroupRateVector(const std::string& group_name) const
{
    return this->getGuideRateVector_(this->groupState().production_rates(group_name));
}

template <typename Scalar, typename IndexTraits>
std::optional<typename SingleWellState<Scalar, IndexTraits>::GroupTarget>
GroupStateHelper<Scalar, IndexTraits>::getWellGroupTargetInjector(const std::string& name,
                                                                  const std::string& parent,
                                                                  const Group& group,
                                                                  const Scalar* rates,
                                                                  Phase injection_phase,
                                                                  const Scalar efficiency_factor,
                                                                  const std::vector<Scalar>& resv_coeff,
                                                                  DeferredLogger& deferred_logger) const
{
    // This function computes a wells group target.
    // 'parent' will be the name of 'group'. But if we recurse, 'name' and
    // 'parent' will stay fixed while 'group' will be higher up
    // in the group tree.
    // Efficiency factor is the well efficiency factor for the first group the well is
    // part of. Later it is the accumulated factor including the group efficiency factor
    // of the child of group.
    auto current_group_control = this->groupState().injection_control(group.name(), injection_phase);
    if (current_group_control == Group::InjectionCMode::FLD
        || current_group_control == Group::InjectionCMode::NONE) {
        // Return if we are not available for parent group.
        if (!group.injectionGroupControlAvailable(injection_phase)) {
            return std::nullopt;
        }
        // Otherwise: check production share of parent's control.
        const auto& parent_group = this->schedule_.getGroup(group.parent(), this->report_step_);
        return this->getWellGroupTargetInjector(name,
                                                parent,
                                                parent_group,
                                                rates,
                                                injection_phase,
                                                efficiency_factor * group.getGroupEfficiencyFactor(),
                                                resv_coeff,
                                                deferred_logger);
    }

    // This can be false for FLD-controlled groups, we must therefore
    // check for FLD first (done above).
    if (!group.isInjectionGroup()) {
        return std::nullopt;
    }

    // If we are here, we are at the topmost group to be visited in the recursion.
    // This is the group containing the control we will check against.
    Scalar sales_target = 0;
    if (this->schedule_[this->report_step_].gconsale().has(group.name())) {
        const auto& gconsale
            = this->schedule_[this->report_step_].gconsale().get(group.name(), this->summary_state_);
        sales_target = gconsale.sales_target;
    }
    GroupStateHelpers::InjectionTargetCalculator<Scalar, IndexTraits> tcalc {
        current_group_control,
        this->phase_usage_info_,
        resv_coeff,
        group.name(),
        sales_target,
        this->groupState(),
        injection_phase,
        group.has_gpmaint_control(injection_phase, current_group_control),
        deferred_logger};

    GroupStateHelpers::FractionCalculator<Scalar, IndexTraits> fcalc {this->schedule_,
                                                                      *this,
                                                                      this->summary_state_,
                                                                      this->report_step_,
                                                                      &this->guide_rate_,
                                                                      tcalc.guideTargetMode(),
                                                                      /*is_producer=*/false,
                                                                      injection_phase};

    auto local_fraction_lambda = [&](const std::string& child, const std::string& always_incluced_name) {
        return fcalc.localFraction(child, always_incluced_name);
    };

    auto local_reduction_lambda = [&](const std::string& group_name) {
        const std::vector<Scalar>& group_target_reductions
            = this->groupState().injection_reduction_rates(group_name);
        return tcalc.calcModeRateFromRates(group_target_reductions);
    };

    std::optional<Group::InjectionControls> ctrl;
    if (!group.has_gpmaint_control(injection_phase, current_group_control))
        ctrl = group.injectionControls(injection_phase, this->summary_state_);


    const Scalar orig_target = tcalc.groupTarget(ctrl, deferred_logger);
    // Assume we have a chain of groups as follows: BOTTOM -> MIDDLE -> TOP.
    // Then ...
    // TODO finish explanation.
    const Scalar current_rate_available
        = tcalc.calcModeRateFromRates(rates); // Switch sign since 'rates' are negative for producers.
    const auto chain = this->groupChainTopBot(name, group.name());
    // Because 'name' is the last of the elements, and not an ancestor, we subtract one below.
    const std::size_t num_ancestors = chain.size() - 1;
    // we need to find out the level where the local reduction is applied
    std::size_t local_reduction_level = 0;
    for (std::size_t ii = 1; ii < num_ancestors; ++ii) {
        const int num_gr_ctrl = this->groupControlledWells(chain[ii],
                                                           /*always_included_child=*/"",
                                                           /*is_production_group=*/false,
                                                           injection_phase);
        if (this->guide_rate_.has(chain[ii], injection_phase) && num_gr_ctrl > 0) {
            local_reduction_level = ii;
        }
    }

    // Compute portion of target corresponding to current_rate_available
    Scalar target = orig_target;
    for (std::size_t ii = 0; ii < num_ancestors; ++ii) {
        if ((ii == 0) || this->guide_rate_.has(chain[ii], injection_phase)) {
            // Apply local reductions only at the control level
            // (top) and for levels where we have a specified
            // group guide rate.
            if (local_reduction_level >= ii) {
                target -= local_reduction_lambda(chain[ii]);
            }

            // If we are under individual control we need to add the wells rate back at the level where it is
            // included in the local reduction
            if (local_reduction_level == ii && !this->wellState().isInjectionGrup(name)) {
                target += current_rate_available * efficiency_factor;
            }
        }
        if (!this->wellState().isInjectionGrup(name)) {
            target *= local_fraction_lambda(chain[ii + 1], name);
        } else {
            target *= local_fraction_lambda(chain[ii + 1], chain[ii + 1]);
        }
    }
    // Avoid negative target rates coming from too large local reductions.
    const auto target_value = std::max(Scalar(0.0), target / efficiency_factor);
    return GroupTarget::injectionGroupTarget(group.name(), current_group_control, target_value);
}

template <typename Scalar, typename IndexTraits>
std::optional<typename SingleWellState<Scalar, IndexTraits>::GroupTarget>
GroupStateHelper<Scalar, IndexTraits>::getWellGroupTargetProducer(const std::string& name,
                                                                  const std::string& parent,
                                                                  const Group& group,
                                                                  const Scalar* rates,
                                                                  const Scalar efficiency_factor,
                                                                  const std::vector<Scalar>& resv_coeff,
                                                                  DeferredLogger& deferred_logger) const
{
    // This function computes a wells group target.
    // 'parent' will be the name of 'group'. But if we recurse, 'name' and
    // 'parent' will stay fixed while 'group' will be higher up
    // in the group tree.
    // Eficiencyfactor is the well efficiency factor for the first group the well is
    // part of. Later it is the accumulated factor including the group efficiency factor
    // of the child of group.
    OPM_TIMEFUNCTION();
    const Group::ProductionCMode& current_group_control = this->groupState().production_control(group.name());

    if (current_group_control == Group::ProductionCMode::FLD
        || current_group_control == Group::ProductionCMode::NONE) {
        // Return if we are not available for parent group.
        if (!group.productionGroupControlAvailable()) {
            return std::nullopt;
        }
        // Otherwise: check production share of parent's control.
        const auto& parent_group = this->schedule_.getGroup(group.parent(), this->report_step_);
        return this->getWellGroupTargetProducer(name,
                                                parent,
                                                parent_group,
                                                rates,
                                                efficiency_factor * group.getGroupEfficiencyFactor(),
                                                resv_coeff,
                                                deferred_logger);
    }

    // This can be false for FLD-controlled groups, we must therefore
    // check for FLD first (done above).
    if (!group.isProductionGroup()) {
        return std::nullopt;
    }

    // If we are here, we are at the topmost group to be visited in the recursion.
    // This is the group containing the control we will check against.

    // gconsale may adjust the grat target.
    // the adjusted rates is sent to the targetCalculator
    Scalar grat_target_from_sales = 0.0;
    if (this->groupState().has_grat_sales_target(group.name()))
        grat_target_from_sales = this->groupState().grat_sales_target(group.name());

    GroupStateHelpers::TargetCalculator<Scalar, IndexTraits> tcalc {current_group_control,
                                                                    this->phase_usage_info_,
                                                                    resv_coeff,
                                                                    grat_target_from_sales,
                                                                    group.name(),
                                                                    this->groupState(),
                                                                    group.has_gpmaint_control(current_group_control)};

    GroupStateHelpers::FractionCalculator<Scalar, IndexTraits> fcalc {this->schedule_,
                                                                      *this,
                                                                      this->summary_state_,
                                                                      this->report_step_,
                                                                      &this->guide_rate_,
                                                                      tcalc.guideTargetMode(),
                                                                      true,
                                                                      Phase::OIL};
    auto local_fraction_lambda = [&](const std::string& child, const std::string& always_incluced_name) {
        return fcalc.localFraction(child, always_incluced_name);
    };

    auto local_reduction_lambda = [&](const std::string& group_name) {
        const std::vector<Scalar>& group_target_reductions
            = this->groupState().production_reduction_rates(group_name);
        return tcalc.calcModeRateFromRates(group_target_reductions);
    };

    std::optional<Group::ProductionControls> ctrl;
    if (!group.has_gpmaint_control(current_group_control))
        ctrl = group.productionControls(this->summary_state_);

    const Scalar orig_target = tcalc.groupTarget(ctrl, deferred_logger);
    // Assume we have a chain of groups as follows: BOTTOM -> MIDDLE -> TOP.
    // Then ...
    // TODO finish explanation.
    const Scalar current_rate_available
        = -tcalc.calcModeRateFromRates(rates); // Switch sign since 'rates' are negative for producers.
    const auto chain = this->groupChainTopBot(name, group.name());
    // Because 'name' is the last of the elements, and not an ancestor, we subtract one below.
    const std::size_t num_ancestors = chain.size() - 1;
    // we need to find out the level where the local reduction is applied
    std::size_t local_reduction_level = 0;
    for (std::size_t ii = 1; ii < num_ancestors; ++ii) {
        const int num_gr_ctrl = this->groupControlledWells(chain[ii],
                                                           /*always_included_child=*/"",
                                                           /*is_production_group=*/true,
                                                           /*injection_phase=*/Phase::OIL);
        if (this->guide_rate_.has(chain[ii]) && num_gr_ctrl > 0) {
            local_reduction_level = ii;
        }
    }
    // Compute portion of target corresponding to current_rate_available
    Scalar target = orig_target;
    for (std::size_t ii = 0; ii < num_ancestors; ++ii) {
        if ((ii == 0) || this->guide_rate_.has(chain[ii])) {
            // Apply local reductions only at the control level
            // (top) and for levels where we have a specified
            // group guide rate.
            if (local_reduction_level >= ii) {
                target -= local_reduction_lambda(chain[ii]);
            }
            // If we are under individual control we need to add the wells rate back at the level where it is
            // included in the local reduction
            if (local_reduction_level == ii && !this->wellState().isProductionGrup(name)) {
                target += current_rate_available * efficiency_factor;
            }
        }
        if (this->wellState().isProductionGrup(name)) {
            target *= local_fraction_lambda(chain[ii + 1], chain[ii + 1]);
        } else {
            target *= local_fraction_lambda(chain[ii + 1], name);
        }
    }
    // Avoid negative target rates coming from too large local reductions.
    const auto target_value = std::max(Scalar(0.0), target / efficiency_factor);
    return GroupTarget::productionGroupTarget(group.name(), current_group_control, target_value);
}

template <typename Scalar, typename IndexTraits>
GuideRate::RateVector
GroupStateHelper<Scalar, IndexTraits>::getWellRateVector(const std::string& name) const
{
    return this->getGuideRateVector_(this->wellState().currentWellRates(name));
}

template <typename Scalar, typename IndexTraits>
std::vector<std::string>
GroupStateHelper<Scalar, IndexTraits>::groupChainTopBot(const std::string& bottom,
                                                       const std::string& top) const
{
    // Get initial parent, 'bottom' can be a well or a group.
    std::string parent;
    if (this->schedule_.hasWell(bottom, this->report_step_)) {
        parent = this->schedule_.getWell(bottom, this->report_step_).groupName();
    } else {
        parent = this->schedule_.getGroup(bottom, this->report_step_).parent();
    }

    // Build the chain from bottom to top.
    std::vector<std::string> chain;
    chain.push_back(bottom);
    chain.push_back(parent);
    while (parent != top) {
        parent = this->schedule_.getGroup(parent, this->report_step_).parent();
        chain.push_back(parent);
    }
    assert(chain.back() == top);

    // Reverse order and return.
    std::reverse(chain.begin(), chain.end());
    return chain;
}

template <typename Scalar, typename IndexTraits>
int
GroupStateHelper<Scalar, IndexTraits>::groupControlledWells(const std::string& group_name,
                                                           const std::string& always_included_child,
                                                           const bool is_production_group,
                                                           const Phase injection_phase) const
{
    auto num_wells = is_production_group
        ? this->groupState().number_of_wells_under_group_control(group_name)
        : this->groupState().number_of_wells_under_inj_group_control(group_name, injection_phase);
    if (this->schedule_.hasWell(always_included_child, this->report_step_)) {
        const bool isInGroup = this->isInGroupChainTopBot_(always_included_child, group_name);
        const bool already_included = is_production_group
            ? this->wellState().isProductionGrup(always_included_child)
            : this->wellState().isInjectionGrup(always_included_child);
        if (!already_included && isInGroup) {
            num_wells++;
        }
    }
    return num_wells;
}

template <typename Scalar, typename IndexTraits>
void
GroupStateHelper<Scalar, IndexTraits>::setCmodeGroup(const Group& group)
{

    for (const std::string& group_name : group.groups()) {
        this->setCmodeGroup(this->schedule_.getGroup(group_name, this->report_step_));
    }

    // use NONE as default control
    const Phase all[] = {Phase::WATER, Phase::OIL, Phase::GAS};
    for (Phase phase : all) {
        if (!this->groupState().has_injection_control(group.name(), phase)) {
            this->groupState().injection_control(group.name(), phase, Group::InjectionCMode::NONE);
        }
    }
    if (!this->groupState().has_production_control(group.name())) {
        this->groupState().production_control(group.name(), Group::ProductionCMode::NONE);
    }

    const auto& events = this->schedule_[this->report_step_].wellgroup_events();
    if (group.isInjectionGroup() && events.hasEvent(group.name(), ScheduleEvents::GROUP_INJECTION_UPDATE)) {

        for (Phase phase : all) {
            if (!group.hasInjectionControl(phase))
                continue;

            const auto& controls = group.injectionControls(phase, this->summary_state_);
            this->groupState().injection_control(group.name(), phase, controls.cmode);
        }
    }

    if (group.isProductionGroup() && events.hasEvent(group.name(), ScheduleEvents::GROUP_PRODUCTION_UPDATE)) {
        const auto controls = group.productionControls(this->summary_state_);
        this->groupState().production_control(group.name(), controls.cmode);
    }

    if (group.has_gpmaint_control(Group::ProductionCMode::RESV)) {
        this->groupState().production_control(group.name(), Group::ProductionCMode::RESV);
    }
    for (Phase phase : all) {
        if (group.has_gpmaint_control(phase, Group::InjectionCMode::RATE)) {
            this->groupState().injection_control(group.name(), phase, Group::InjectionCMode::RATE);
        } else if (group.has_gpmaint_control(phase, Group::InjectionCMode::RESV)) {
            this->groupState().injection_control(group.name(), phase, Group::InjectionCMode::RESV);
        }
    }

    if (this->schedule_[this->report_step_].gconsale().has(group.name())) {
        this->groupState().injection_control(group.name(), Phase::GAS, Group::InjectionCMode::SALE);
    }
}

// NOTE: setRegionAveragePressureCalculator() is a template member function and must be
//       defined in the header file (see GroupStateHelper.hpp).
//
// See: updateGpMaintTargetForGroups() for detailed rationale.

template <typename Scalar, typename IndexTraits>
Scalar
GroupStateHelper<Scalar, IndexTraits>::sumSolventRates(const Group& group, const bool is_injector) const
{
    Scalar rate = 0.0;
    for (const std::string& group_name : group.groups()) {
        const Group& group_tmp = this->schedule_.getGroup(group_name, this->report_step_);
        const auto& gefac = group_tmp.getGroupEfficiencyFactor();
        rate += gefac * this->sumSolventRates(group_tmp, is_injector);
    }

    for (const std::string& well_name : group.wells()) {
        const auto well_index = this->wellState().index(well_name);
        if (!well_index.has_value())
            continue;

        if (!this->wellState().wellIsOwned(well_index.value(), well_name)) // Only sum once
        {
            continue;
        }

        const auto& well_ecl = this->schedule_.getWell(well_name, this->report_step_);
        // only count producers or injectors
        if ((well_ecl.isProducer() && is_injector) || (well_ecl.isInjector() && !is_injector))
            continue;

        const auto& ws = this->wellState().well(well_index.value());

        if (ws.status == Well::Status::SHUT)
            continue;

        const Scalar factor
            = well_ecl.getEfficiencyFactor() * this->wellState()[well_ecl.name()].efficiency_scaling_factor;
        if (is_injector)
            rate += factor * ws.sum_solvent_rates();
        else
            rate -= factor * ws.sum_solvent_rates();
    }
    return rate;
}

template <typename Scalar, typename IndexTraits>
Scalar
GroupStateHelper<Scalar, IndexTraits>::sumWellPhaseRates(bool res_rates,
                                                        const Opm::Group& group,
                                                        const int phase_pos,
                                                        const bool injector,
                                                        const bool network) const
{
    // Only obtain satellite rates once (on rank 0)
    Scalar rate = 0.0;
    if (this->wellState().isRank0() && (group.hasSatelliteProduction() || group.hasSatelliteInjection())) {
        if (injector) {
            rate = this->satelliteInjectionRate_(
                this->schedule_[this->report_step_], group, phase_pos, res_rates);
        } else {
            const auto rate_comp = this->selectRateComponent_(phase_pos);
            if (rate_comp.has_value()) {
                rate = this->satelliteProductionRate_(
                    this->schedule_[this->report_step_], group, *rate_comp, res_rates);
            }
        }
        // Satellite groups have no sub groups/wells so we're done
        return rate;
    }

    for (const std::string& group_name : group.groups()) {
        const auto& group_tmp = this->schedule_.getGroup(group_name, this->report_step_);
        const auto& gefac = group_tmp.getGroupEfficiencyFactor(network);
        rate += gefac * this->sumWellPhaseRates(res_rates, group_tmp, phase_pos, injector, network);
    }

    for (const std::string& well_name : group.wells()) {
        const auto well_index = this->wellState().index(well_name);
        if (!well_index.has_value())
            continue;

        if (!this->wellState().wellIsOwned(well_index.value(), well_name)) // Only sum once
        {
            continue;
        }

        const auto& well_ecl = this->schedule_.getWell(well_name, this->report_step_);
        // only count producers or injectors
        if ((well_ecl.isProducer() && injector) || (well_ecl.isInjector() && !injector))
            continue;

        const auto& ws = this->wellState().well(well_index.value());
        if (ws.status == Opm::Well::Status::SHUT)
            continue;

        const Scalar factor = well_ecl.getEfficiencyFactor(network)
            * this->wellState().well(well_index.value()).efficiency_scaling_factor;
        if (res_rates) {
            const auto& well_rates = ws.reservoir_rates;
            if (injector)
                rate += factor * well_rates[phase_pos];
            else
                rate -= factor * well_rates[phase_pos];
        } else {
            const auto& well_rates = ws.surface_rates;
            if (injector)
                rate += factor * well_rates[phase_pos];
            else
                rate -= factor * well_rates[phase_pos];
        }
    }
    return rate;
}

template <typename Scalar, typename IndexTraits>
Scalar
GroupStateHelper<Scalar, IndexTraits>::sumWellResRates(const Group& group,
                                                      const int phasePos,
                                                      const bool injector) const
{
    return this->sumWellPhaseRates(/*res_rates=*/true, group, phasePos, injector);
}

template <typename Scalar, typename IndexTraits>
Scalar
GroupStateHelper<Scalar, IndexTraits>::sumWellSurfaceRates(const Group& group,
                                                          const int phase_pos,
                                                          const bool injector) const
{
    return this->sumWellPhaseRates(/*res_rates=*/false, group, phase_pos, injector);
}

// NOTE: updateGpMaintTargetForGroups() is a template member function and must be
//       defined in the header file (see GroupStateHelper.hpp).
//
// Rationale:
//   - The RegionalValues template parameter depends on BlackoilWellModel's FluidSystem,
//     which varies by derived class and is not available in this compilation unit.
//   - Template functions must be visible at their point of instantiation, so the
//     implementation cannot be in this .cpp file.
//   - We cannot use explicit template instantiation because the specific RegionalValues
//     type depends on types from derived classes (e.g., AverageRegionalPressure<FluidSystem>).
//
// Historical Note:
//   - The original implementation in WellGroupHelpers.cpp (static utility class) could keep
//     the implementation in the .cpp file because it hard-coded BlackOilFluidSystem<Scalar>
//     and used explicit template instantiation for that specific type.
//   - As a member-based class, GroupStateHelper must support whatever FluidSystem the TypeTag
//     provides, so we cannot hard-code the type and must use header-based implementation.
//
// If moved here, the linker will report:
//   undefined reference to
//   `GroupStateHelper<...>::updateGpMaintTargetForGroups<std::map<...>>(const Group&, ...)`

template <typename Scalar, typename IndexTraits>
int
GroupStateHelper<Scalar, IndexTraits>::updateGroupControlledWells(const bool is_production_group,
                                                                 const Phase injection_phase,
                                                                 DeferredLogger& deferred_logger)
{
    OPM_TIMEFUNCTION();
    const auto& group_name = "FIELD";
    return this->updateGroupControlledWellsRecursive_(
        group_name, is_production_group, injection_phase, deferred_logger);
}

template <typename Scalar, typename IndexTraits>
void
GroupStateHelper<Scalar, IndexTraits>::updateGroupProductionRates(const Group& group)
{
    OPM_TIMEFUNCTION();
    for (const std::string& group_name : group.groups()) {
        const Group& group_tmp = this->schedule_.getGroup(group_name, this->report_step_);
        this->updateGroupProductionRates(group_tmp);
    }
    const int np = this->wellState().numPhases();
    std::vector<Scalar> rates(np, 0.0);
    for (int phase = 0; phase < np; ++phase) {
        rates[phase] = this->sumWellPhaseRates(
            /*res_rates=*/false, group, phase, /*injector=*/false, /*network=*/false);
    }
    this->groupState().update_production_rates(group.name(), rates);
}

template <typename Scalar, typename IndexTraits>
void
GroupStateHelper<Scalar, IndexTraits>::updateGroupTargetReduction(const Group& group,
                                                                 const bool is_injector)
{
    OPM_TIMEFUNCTION();
    std::vector<Scalar> group_target_reduction(this->wellState().numPhases(), 0.0);
    this->updateGroupTargetReductionRecursive_(group, is_injector, group_target_reduction);
}

template <typename Scalar, typename IndexTraits>
void
GroupStateHelper<Scalar, IndexTraits>::updateNetworkLeafNodeProductionRates()
{
    const auto& network = this->schedule_[this->report_step_].network();
    if (network.active()) {
        const int np = this->wellState().numPhases();
        for (const auto& group_name : network.leaf_nodes()) {
            std::vector<Scalar> network_rates(np, 0.0);
            if (this->schedule_[this->report_step_].groups.has(
                    group_name)) { // Allow empty leaf nodes that are not groups
                const auto& group = this->schedule_[this->report_step_].groups.get(group_name);
                if (group.numWells() > 0) {
                    for (int phase = 0; phase < np; ++phase) {
                        network_rates[phase] = this->sumWellPhaseRates(
                            /*res_rates=*/false, group, phase, /*injector=*/false, /*network=*/true);
                    }
                }
            }
            this->groupState().update_network_leaf_node_production_rates(group_name, network_rates);
        }
    }
}

template <typename Scalar, typename IndexTraits>
void
GroupStateHelper<Scalar, IndexTraits>::updateREINForGroups(const Group& group,
                                                          bool sum_rank)
{
    OPM_TIMEFUNCTION();
    const int np = this->wellState().numPhases();
    for (const std::string& groupName : group.groups()) {
        const Group& group_tmp = this->schedule_.getGroup(groupName, this->report_step_);
        this->updateREINForGroups(group_tmp, sum_rank);
    }

    std::vector<Scalar> rein(np, 0.0);
    for (int phase = 0; phase < np; ++phase) {
        rein[phase] = this->sumWellPhaseRates(
            /*res_rates=*/false, group, phase, /*injector=*/false, /*network=*/false);
    }

    // add import rate and subtract consumption rate for group for gas
    if (sum_rank) {
        const auto& pu = this->wellState().phaseUsageInfo();
        if (pu.phaseIsActive(IndexTraits::gasPhaseIdx)) {
            const auto& [consumption_rate, import_rate] = this->groupState().gconsump_rates(group.name());
            const auto ig = pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx);
            rein[ig] += import_rate;
            rein[ig] -= consumption_rate;
        }
    }

    this->groupState().update_injection_rein_rates(group.name(), rein);
}

template <typename Scalar, typename IndexTraits>
void
GroupStateHelper<Scalar, IndexTraits>::updateReservoirRatesInjectionGroups(const Group& group)
{
    OPM_TIMEFUNCTION();
    for (const std::string& group_name : group.groups()) {
        const Group& group_tmp = this->schedule_.getGroup(group_name, this->report_step_);
        this->updateReservoirRatesInjectionGroups(group_tmp);
    }
    const int np = this->wellState().numPhases();
    std::vector<Scalar> resv(np, 0.0);
    for (int phase = 0; phase < np; ++phase) {
        resv[phase] = this->sumWellPhaseRates(
            /*res_rates=*/true, group, phase, /*injector=*/true, /*network=*/false);
    }
    this->groupState().update_injection_reservoir_rates(group.name(), resv);
}

template <typename Scalar, typename IndexTraits>
void
GroupStateHelper<Scalar, IndexTraits>::updateState(WellState<Scalar, IndexTraits>& well_state,
                                                  GroupState<Scalar>& group_state)
{
    this->well_state_ = &well_state;
    this->group_state_ = &group_state;
}

template <typename Scalar, typename IndexTraits>
void
GroupStateHelper<Scalar, IndexTraits>::updateSurfaceRatesInjectionGroups(const Group& group)
{
    OPM_TIMEFUNCTION();
    for (const std::string& group_name : group.groups()) {
        const Group& group_tmp = this->schedule_.getGroup(group_name, this->report_step_);
        this->updateSurfaceRatesInjectionGroups(group_tmp);
    }
    const int np = this->wellState().numPhases();
    std::vector<Scalar> rates(np, 0.0);
    for (int phase = 0; phase < np; ++phase) {
        rates[phase] = this->sumWellPhaseRates(
            /*res_rates=*/false,
            group,
            phase,
            /*injector=*/true,
            /*network=*/false);
    }
    this->groupState().update_injection_surface_rates(group.name(), rates);
}

template <typename Scalar, typename IndexTraits>
void
GroupStateHelper<Scalar, IndexTraits>::updateVREPForGroups(const Group& group)
{
    OPM_TIMEFUNCTION();
    for (const std::string& group_name : group.groups()) {
        const Group& group_tmp = this->schedule_.getGroup(group_name, this->report_step_);
        this->updateVREPForGroups(group_tmp);
    }
    const int np = this->wellState().numPhases();
    Scalar resv = 0.0;
    for (int phase = 0; phase < np; ++phase) {
        resv += this->sumWellPhaseRates(
            /*res_rates=*/true, group, phase, /*injector=*/false, /*network=*/false);
    }
    this->groupState().update_injection_vrep_rate(group.name(), resv);
}

template <typename Scalar, typename IndexTraits>
void
GroupStateHelper<Scalar, IndexTraits>::updateWellRates(const Group& group,
                                                      const WellState<Scalar, IndexTraits>& well_state_nupcol,
                                                      WellState<Scalar, IndexTraits>& well_state) const
{
    OPM_TIMEFUNCTION();
    for (const std::string& group_name : group.groups()) {
        const Group& group_tmp = this->schedule_.getGroup(group_name, this->report_step_);
        this->updateWellRates(group_tmp, well_state_nupcol, well_state);
    }
    const int np = well_state.numPhases();
    for (const std::string& well_name : group.wells()) {
        std::vector<Scalar> rates(np, 0.0);
        const auto well_index = well_state.index(well_name);
        if (well_index.has_value()) { // the well is found on this node
            const auto& well_tmp = this->schedule_.getWell(well_name, this->report_step_);
            int sign = 1;
            // production wellRates are negative. The users of currentWellRates uses the convention in
            // opm-common that production and injection rates are positive.
            if (!well_tmp.isInjector())
                sign = -1;
            const auto& ws = well_state_nupcol.well(well_index.value());
            for (int phase = 0; phase < np; ++phase) {
                rates[phase] = sign * ws.surface_rates[phase];
            }
        }
        well_state.setCurrentWellRates(well_name, rates);
    }
}

template <typename Scalar, typename IndexTraits>
void
GroupStateHelper<Scalar, IndexTraits>::updateWellRatesFromGroupTargetScale(
    const Scalar scale,
    const Group& group,
    bool is_injector,
    WellState<Scalar, IndexTraits>& well_state) const
{
    OPM_TIMEFUNCTION();
    for (const std::string& group_name : group.groups()) {
        bool individual_control = false;
        if (is_injector) {
            const Phase all[] = {Phase::WATER, Phase::OIL, Phase::GAS};
            for (Phase phase : all) {
                const Group::InjectionCMode& current_group_control
                    = this->groupState().injection_control(group_name, phase);
                individual_control = individual_control
                    || (current_group_control != Group::InjectionCMode::FLD
                        && current_group_control != Group::InjectionCMode::NONE);
            }
        } else {
            const Group::ProductionCMode& current_group_control
                = this->groupState().production_control(group_name);
            individual_control = (current_group_control != Group::ProductionCMode::FLD
                                  && current_group_control != Group::ProductionCMode::NONE);
        }
        if (!individual_control) {
            const Group& group_tmp = this->schedule_.getGroup(group_name, this->report_step_);
            this->updateWellRatesFromGroupTargetScale(scale, group_tmp, is_injector, well_state);
        }
    }

    const int np = well_state.numPhases();
    for (const std::string& well_name : group.wells()) {
        const auto& well_tmp = this->schedule_.getWell(well_name, this->report_step_);

        if (well_tmp.isProducer() && is_injector)
            continue;

        if (well_tmp.isInjector() && !is_injector)
            continue;

        const auto well_index = well_state.index(well_name);
        if (!well_index.has_value())
            continue;

        auto& ws = well_state.well(well_index.value());
        if (ws.status == Well::Status::SHUT)
            continue;

        // scale rates
        if (is_injector) {
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

template <typename Scalar, typename IndexTraits>
std::pair<std::optional<std::string>, Scalar>
GroupStateHelper<Scalar, IndexTraits>::worstOffendingWell(const Group& group,
                                                         const Group::ProductionCMode& offended_control,
                                                         const Parallel::Communication& comm,
                                                         DeferredLogger& deferred_logger) const
{
    std::pair<std::optional<std::string>, Scalar> offending_well {std::nullopt, 0.0};
    for (const std::string& child_group : group.groups()) {
        const auto& this_group = this->schedule_.getGroup(child_group, this->report_step_);
        const auto& offending_well_this = this->worstOffendingWell(
            this_group, offended_control, comm, deferred_logger
        );
        if (offending_well_this.second > offending_well.second) {
            offending_well = offending_well_this;
        }
    }

    for (const std::string& child_well : group.wells()) {

        const auto& well_index = this->wellState().index(child_well);
        Scalar violating_rate = 0.0;
        Scalar prefered_rate = 0.0;
        const auto& pu = this->phase_usage_info_;
        if (well_index.has_value() && this->wellState().wellIsOwned(well_index.value(), child_well)) {
            const auto& ws = this->wellState().well(child_well);
            switch (offended_control) {
            case Group::ProductionCMode::ORAT: {
                const int oil_pos = pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
                violating_rate = ws.surface_rates[oil_pos];
                break;
            }
            case Group::ProductionCMode::GRAT: {
                const int gas_pos = pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx);
                violating_rate = ws.surface_rates[gas_pos];
                break;
            }
            case Group::ProductionCMode::WRAT: {
                const int water_pos = pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx);
                violating_rate = ws.surface_rates[water_pos];
                break;
            }
            case Group::ProductionCMode::LRAT: {
                assert(pu.phaseIsActive(IndexTraits::oilPhaseIdx));
                assert(pu.phaseIsActive(IndexTraits::waterPhaseIdx));
                const int oil_pos = pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
                const int water_pos = pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx);
                violating_rate = ws.surface_rates[oil_pos] + ws.surface_rates[water_pos];
                break;
            }
            case Group::ProductionCMode::RESV:
                for (int p = 0; p < this->wellState().numPhases(); ++p) {
                    violating_rate += ws.reservoir_rates[p];
                }
                break;
            case Group::ProductionCMode::NONE:
                break;
            case Group::ProductionCMode::FLD:
                break;
            case Group::ProductionCMode::PRBL:
                OPM_DEFLOG_THROW(std::runtime_error,
                                 "Group " + group.name() + " GroupProductionCMode PRBL not implemented",
                                 deferred_logger);
                break;
            case Group::ProductionCMode::CRAT:
                OPM_DEFLOG_THROW(std::runtime_error,
                                 "Group " + group.name() + " GroupProductionCMode CRAT not implemented",
                                 deferred_logger);
                break;
            }
            const auto preferred_phase
                = this->schedule_.getWell(child_well, this->report_step_).getPreferredPhase();
            switch (preferred_phase) {
            case Phase::OIL: {
                const int oil_pos
                    = this->phase_usage_info_.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
                prefered_rate = ws.surface_rates[oil_pos];
                break;
            }
            case Phase::GAS: {
                const int gas_pos
                    = this->phase_usage_info_.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx);
                prefered_rate = ws.surface_rates[gas_pos];
                break;
            }
            case Phase::WATER: {
                const int water_pos
                    = this->phase_usage_info_.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx);
                prefered_rate = ws.surface_rates[water_pos];
                break;
            }
            default:
                // No others supported.
                break;
            }
        }
        violating_rate = comm.sum(violating_rate);
        if (violating_rate < 0) { // only check producing wells
            prefered_rate = comm.sum(prefered_rate);
            Scalar fraction = prefered_rate < -1e-16 ? violating_rate / prefered_rate : 1.0;
            if (fraction > offending_well.second) {
                offending_well = {child_well, fraction};
            }
        }
    }
    return offending_well;
}

// ---------------------------------------------------------------------
// Private methods
// ---------------------------------------------------------------------

template <typename Scalar, typename IndexTraits>
std::string
GroupStateHelper<Scalar, IndexTraits>::controlGroup_(const Group& group) const
{
    const Group::ProductionCMode& currentGroupControl = this->groupState().production_control(group.name());

    if (currentGroupControl == Group::ProductionCMode::FLD
        || currentGroupControl == Group::ProductionCMode::NONE) {
        const auto& parent_name = group.control_group();
        if (parent_name) {
            const auto& parent = this->schedule_.getGroup(parent_name.value(), this->report_step_);
            return this->controlGroup_(parent);
        }
    }

    return group.name();
}

template <class Scalar, typename IndexTraits>
Opm::GuideRate::RateVector
GroupStateHelper<Scalar, IndexTraits>::getGuideRateVector_(const std::vector<Scalar>& rates) const
{
    const auto& pu = this->phase_usage_info_;
    Scalar oilRate = 0.0;
    if (pu.phaseIsActive(IndexTraits::oilPhaseIdx)) {
        oilRate = rates[pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx)];
    }

    Scalar gasRate = 0.0;
    if (pu.phaseIsActive(IndexTraits::gasPhaseIdx)) {
        gasRate = rates[pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx)];
    }

    Scalar waterRate = 0.0;
    if (pu.phaseIsActive(IndexTraits::waterPhaseIdx)) {
        waterRate = rates[pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx)];
    }

    return {oilRate, gasRate, waterRate};
}

template <typename Scalar, typename IndexTraits>
bool
GroupStateHelper<Scalar, IndexTraits>::isInGroupChainTopBot_(const std::string& bottom,
                                                            const std::string& top) const
{
    // Get initial parent, 'bottom' can be a well or a group.
    std::string parent;
    if (this->schedule_.hasWell(bottom, this->report_step_)) {
        parent = this->schedule_.getWell(bottom, this->report_step_).groupName();
    } else {
        parent = this->schedule_.getGroup(bottom, this->report_step_).parent();
    }

    while (parent != top) {
        parent = this->schedule_.getGroup(parent, this->report_step_).parent();
        if (parent == top) {
            return true;
        } else if (parent == "FIELD") {
            return false;
        }
    }
    return true;
}

template <typename Scalar, typename IndexTraits>
int
GroupStateHelper<Scalar, IndexTraits>::phaseToActivePhaseIdx_(const Phase phase) const
{
    const auto& pu = this->phase_usage_info_;
    int phase_pos = -1;
    switch (phase) {
    case Phase::GAS:
        if (pu.phaseIsActive(IndexTraits::gasPhaseIdx)) {
            phase_pos = pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx);
        }
        break;
    case Phase::OIL:
        if (pu.phaseIsActive(IndexTraits::oilPhaseIdx)) {
            phase_pos = pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
        }
        break;
    case Phase::WATER:
        if (pu.phaseIsActive(IndexTraits::waterPhaseIdx)) {
            phase_pos = pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx);
        }
        break;
    default:
        // just to avoid warning
        throw std::invalid_argument("unhandled phase enum");
    }
    return phase_pos;
}

template <typename Scalar, typename IndexTraits>
Scalar
GroupStateHelper<Scalar, IndexTraits>::satelliteInjectionRate_(const ScheduleState& sched,
                                                              const Group& group,
                                                              const int phase_pos,
                                                              bool res_rates) const
{
    Scalar rate = 0.0;
    const auto& pu = this->phase_usage_info_;
    if (group.hasSatelliteInjection()) {
        std::optional<Phase> ph;
        for (const auto& [bo_phase, phase] : std::array {std::pair(IndexTraits::waterPhaseIdx, Phase::WATER),
                                                         std::pair(IndexTraits::oilPhaseIdx, Phase::OIL),
                                                         std::pair(IndexTraits::gasPhaseIdx, Phase::GAS)}) {
            if (pu.phaseIsActive(bo_phase) && (pu.canonicalToActivePhaseIdx(bo_phase) == phase_pos)) {
                ph = phase;
            }
        }
        if (ph.has_value()) {
            const auto& satellite_inj = sched.satelliteInjection(group.name());
            const auto& rate_ix = satellite_inj.rateIndex(ph.value());
            if (rate_ix.has_value()) {
                const auto& satrates = satellite_inj[*rate_ix];
                if (!res_rates) { // surface rates
                    if (const auto& qs = satrates.surface(); qs.has_value()) {
                        rate = *qs;
                    }
                    // We don't support reservoir rates for satellite injection groups
                }
            }
        }
    }
    return rate;
}

template <typename Scalar, typename IndexTraits>
Scalar
GroupStateHelper<Scalar, IndexTraits>::satelliteProductionRate_(
    const ScheduleState& sched,
    const Group& group,
    const GSatProd::GSatProdGroupProp::Rate rate_comp,
    bool res_rates) const
{
    Scalar rate = 0.0;
    if (group.hasSatelliteProduction()) {
        const auto& gsat_prod = sched.gsatprod();
        if (!res_rates) {
            rate = gsat_prod.get(group.name(), this->summary_state_).rate[rate_comp];
        }
        // We don't support reservoir rates for satellite production groups
    }
    return rate;
}

template <typename Scalar, typename IndexTraits>
std::optional<GSatProd::GSatProdGroupProp::Rate>
GroupStateHelper<Scalar, IndexTraits>::selectRateComponent_(const int phase_pos) const
{
    // TODO: this function can be wrong, phase_pos is not used anymore, this function requries checking and
    // refactoring.
    const auto& pu = this->phase_usage_info_;
    using Rate = GSatProd::GSatProdGroupProp::Rate;

    for (const auto& [phase, rate_comp] : std::array {std::pair {IndexTraits::waterPhaseIdx, Rate::Water},
                                                      std::pair {IndexTraits::oilPhaseIdx, Rate::Oil},
                                                      std::pair {IndexTraits::gasPhaseIdx, Rate::Gas}}) {
        if (pu.phaseIsActive(phase)) {
            const auto active_phase_pos = pu.canonicalToActivePhaseIdx(phase);
            if (active_phase_pos == phase_pos) {
                return rate_comp;
            }
        }
    }

    return std::nullopt;
}

template <typename Scalar, typename IndexTraits>
int
GroupStateHelper<Scalar, IndexTraits>::updateGroupControlledWellsRecursive_(
    const std::string& group_name,
    const bool is_production_group,
    const Phase injection_phase,
    DeferredLogger& deferred_logger)
{
    const Group& group = this->schedule_.getGroup(group_name, this->report_step_);
    int num_wells = 0;
    for (const std::string& child_group : group.groups()) {

        bool included = false;
        if (is_production_group) {
            const auto ctrl = this->groupState().production_control(child_group);
            included = (ctrl == Group::ProductionCMode::FLD || ctrl == Group::ProductionCMode::NONE);
        } else {
            const auto ctrl = this->groupState().injection_control(child_group, injection_phase);
            included = (ctrl == Group::InjectionCMode::FLD || ctrl == Group::InjectionCMode::NONE);
        }

        if (included) {
            num_wells += this->updateGroupControlledWellsRecursive_(
                child_group, is_production_group, injection_phase, deferred_logger);
        } else {
            this->updateGroupControlledWellsRecursive_(
                child_group, is_production_group, injection_phase, deferred_logger);
        }
    }
    for (const std::string& child_well : group.wells()) {
        bool included = false;
        const Well& well = this->schedule_.getWell(child_well, this->report_step_);
        if (is_production_group && well.isProducer()) {
            included = (this->wellState().isProductionGrup(child_well) || group.as_choke());
        } else if (!is_production_group && !well.isProducer()) {
            const auto& well_controls = well.injectionControls(this->summary_state_);
            auto injectorType = well_controls.injector_type;
            if ((injection_phase == Phase::WATER && injectorType == InjectorType::WATER)
                || (injection_phase == Phase::OIL && injectorType == InjectorType::OIL)
                || (injection_phase == Phase::GAS && injectorType == InjectorType::GAS)) {
                included = this->wellState().isInjectionGrup(child_well);
            }
        }
        const auto ctrl1 = this->groupState().production_control(group.name());
        if (group.as_choke()
            && ((ctrl1 == Group::ProductionCMode::FLD) || (ctrl1 == Group::ProductionCMode::NONE))) {
            // The auto choke group has not own group control but inherits control from an ancestor group.
            // Number of wells should be calculated as zero when wells of auto choke group do not deliver
            // target. This behaviour is then similar to no-autochoke group with wells not on GRUP control.
            // The rates of these wells are summed up. The parent group target is reduced with this rate.
            // This reduced target becomes the target of the other child group of this parent.
            const auto num_phases = this->wellState().numPhases();
            std::vector<Scalar> rates(num_phases, 0.0);
            for (int phase_pos = 0; phase_pos < num_phases; ++phase_pos) {
                rates[phase_pos] = this->sumWellSurfaceRates(group, phase_pos, /*injector=*/false);
            }

            // Get the ancestor of the auto choke group that has group control (cmode != FLD, NONE)
            const auto& control_group_name = this->controlGroup_(group);
            const auto& control_group = this->schedule_.getGroup(control_group_name, this->report_step_);
            const auto& ctrl = control_group.productionControls(this->summary_state_);
            const auto& control_group_cmode = ctrl.cmode;

            const auto& group_guide_rate = group.productionControls(this->summary_state_).guide_rate;

            if (group_guide_rate > 0) {
                // Guide rate is not default for the auto choke group
                Scalar grat_target_from_sales = 0.0;
                if (this->groupState().has_grat_sales_target(control_group_name))
                    grat_target_from_sales = this->groupState().grat_sales_target(control_group_name);

                std::vector<Scalar> resv_coeff(num_phases, 1.0);
                GroupStateHelpers::TargetCalculator<Scalar, IndexTraits> tcalc {
                    control_group_cmode,
                    this->phase_usage_info_,
                    resv_coeff,
                    grat_target_from_sales,
                    group.name(),
                    this->groupState(),
                    group.has_gpmaint_control(control_group_cmode)};
                const auto& control_group_target = tcalc.groupTarget(ctrl, deferred_logger);

                // Calculates the guide rate of the parent group with control.
                // It is allowed that the guide rate of this group is defaulted. The guide rate will be
                // derived from the children groups
                const auto& control_group_guide_rate
                    = this->getGuideRate(control_group_name, tcalc.guideTargetMode());

                if (control_group_guide_rate > 0) {
                    // Target rate for the auto choke group
                    const Scalar target_rate
                        = control_group_target * group_guide_rate / control_group_guide_rate;
                    const Scalar current_rate = tcalc.calcModeRateFromRates(rates);

                    if (current_rate < target_rate)
                        included = false;
                }
            }
        }

        if (included) {
            ++num_wells;
        }
    }
    if (is_production_group) {
        this->groupState().update_number_of_wells_under_group_control(group_name, num_wells);
    } else {
        this->groupState().update_number_of_wells_under_inj_group_control(group_name, injection_phase, num_wells);
    }

    return num_wells;
}

template <typename Scalar, typename IndexTraits>
void
GroupStateHelper<Scalar, IndexTraits>::updateGroupTargetReductionRecursive_(
    const Group& group,
    const bool is_injector,
    std::vector<Scalar>& group_target_reduction)
{
    const int np = this->wellState().numPhases();
    for (const std::string& sub_group_name : group.groups()) {
        std::vector<Scalar> sub_group_target_reduction(np, 0.0);
        const Group& sub_group = this->schedule_.getGroup(sub_group_name, this->report_step_);
        this->updateGroupTargetReductionRecursive_(
            sub_group, is_injector, sub_group_target_reduction);

        const Scalar sub_group_efficiency = sub_group.getGroupEfficiencyFactor();

        // accumulate group contribution from sub group
        if (is_injector) {
            const Phase all[] = {Phase::WATER, Phase::OIL, Phase::GAS};
            for (Phase phase : all) {
                const auto phase_pos = this->phaseToActivePhaseIdx_(phase);
                // the phase is not present
                if (phase_pos == -1)
                    continue;
                const Group::InjectionCMode& current_group_control
                    = this->groupState().injection_control(sub_group.name(), phase);
                const bool individual_control = (current_group_control != Group::InjectionCMode::FLD
                                                 && current_group_control != Group::InjectionCMode::NONE);
                const int num_group_controlled_wells
                    = this->groupControlledWells(sub_group.name(),
                                                 /*always_included_child=*/"",
                                                 !is_injector,
                                                 phase);
                if (individual_control || num_group_controlled_wells == 0) {
                    group_target_reduction[phase_pos] += sub_group_efficiency
                        * this->sumWellSurfaceRates(sub_group, phase_pos, is_injector);
                } else {
                    // Accumulate from this subgroup only if no group guide rate is set for it.
                    if (!this->guide_rate_.has(sub_group.name(), phase)) {
                        group_target_reduction[phase_pos]
                            += sub_group_efficiency * sub_group_target_reduction[phase_pos];
                    }
                }
            }
        } else {
            const Group::ProductionCMode& current_group_control
                = this->groupState().production_control(sub_group.name());
            const bool individual_control = (current_group_control != Group::ProductionCMode::FLD
                                             && current_group_control != Group::ProductionCMode::NONE);
            const int num_group_controlled_wells
                = this->groupControlledWells(sub_group.name(),
                                             /*always_included_child=*/"",
                                             !is_injector,
                                             /*injection_phase=*/Phase::OIL);
            if (individual_control || num_group_controlled_wells == 0) {
                for (int phase = 0; phase < np; phase++) {
                    group_target_reduction[phase]
                        += sub_group_efficiency * this->sumWellSurfaceRates(sub_group, phase, is_injector);
                }
            } else {
                // The subgroup may participate in group control.
                if (!this->guide_rate_.has(sub_group.name())) {
                    // Accumulate from this subgroup only if no group guide rate is set for it.
                    for (int phase = 0; phase < np; phase++) {
                        group_target_reduction[phase]
                            += sub_group_efficiency * sub_group_target_reduction[phase];
                    }
                }
            }
        }
    }

    for (const std::string& well_name : group.wells()) {
        const auto& well_tmp = this->schedule_.getWell(well_name, this->report_step_);

        if (well_tmp.isProducer() && is_injector)
            continue;

        if (well_tmp.isInjector() && !is_injector)
            continue;

        const auto well_index = this->wellState().index(well_name);
        if (!well_index.has_value())
            continue;

        if (!this->wellState().wellIsOwned(well_index.value(), well_name)) // Only sum once
        {
            continue;
        }

        const auto& ws = this->wellState().well(well_index.value());
        if (ws.status == Well::Status::SHUT)
            continue;

        const Scalar efficiency = well_tmp.getEfficiencyFactor()
            * this->wellState().well(well_index.value()).efficiency_scaling_factor;

        // add contribution from wells not under group control
        if (is_injector) {
            if (ws.injection_cmode != Well::InjectorCMode::GRUP)
                for (int phase = 0; phase < np; phase++) {
                    group_target_reduction[phase] += ws.surface_rates[phase] * efficiency;
                }
        } else {
            if ((ws.production_cmode != Well::ProducerCMode::GRUP)) {
                if (!group.as_choke()) {
                    for (int phase = 0; phase < np; phase++) {
                        group_target_reduction[phase] -= ws.surface_rates[phase] * efficiency;
                    }
                }
            }
        }
    }
    if (is_injector) {
        this->groupState().update_injection_reduction_rates(group.name(), group_target_reduction);
    } else {
        this->groupState().update_production_reduction_rates(group.name(), group_target_reduction);
    }
}


template class GroupStateHelper<double, BlackOilDefaultFluidSystemIndices>;
#ifdef FLOW_INSTANTIATE_FLOAT
template class GroupStateHelper<float, BlackOilDefaultFluidSystemIndices>;
#endif
} // namespace Opm
