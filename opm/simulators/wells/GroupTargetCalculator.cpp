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
#include <opm/input/eclipse/Schedule/Group/GuideRateModel.hpp>
#include <opm/simulators/wells/GroupTargetCalculator.hpp>

#include <tuple>  // for std::tie
#include <fmt/format.h>

namespace Opm {

// ----------------------------------------------------
// Constructor for the GuideRateTargetCalculator class
// ----------------------------------------------------
template<class Scalar>
GroupTargetCalculator<Scalar>::
GroupTargetCalculator(
    const BlackoilWellModelGeneric<Scalar>& well_model,
    const WellState<Scalar>& well_state,
    const GroupState<Scalar>& group_state,
    const Schedule& schedule,
    const SummaryState& summary_state,
    const PhaseUsage& phase_usage,
    const GuideRate& guide_rate,
    const int report_step_idx,
    DeferredLogger& deferred_logger
) :
    well_model_{well_model},
    well_state_{well_state},
    group_state_{group_state},
    schedule_{schedule},
    summary_state_{summary_state},
    phase_usage_{phase_usage},
    guide_rate_{guide_rate},
    report_step_idx_{report_step_idx},
    deferred_logger_{deferred_logger},
    resv_coeffs_inj_(phase_usage.num_phases, 0.0)
{
    std::tie(this->fipnum_, this->pvtreg_) = this->well_model_.getGroupFipnumAndPvtreg();
    // Get the reservoir coefficients for injection
    // NOTE: The production reservoir coefficients depend on the group production rates, so we
    //  cannot precompute them here.
    this->well_model_.calcInjResvCoeff(this->fipnum_, this->pvtreg_, this->resv_coeffs_inj_);
}

template<class Scalar>
std::optional<typename GroupTargetCalculator<Scalar>::InjectionTargetInfo>
GroupTargetCalculator<Scalar>::
groupInjectionTarget(const Group& group, Phase injection_phase) {
    GeneralCalculator general_calculator{*this, group, injection_phase};
    auto target_info = general_calculator.calculateGroupTarget();
    if (!target_info) {
        // No target was calculated, return an empty optional.
        return std::nullopt;
    }
    return std::make_optional<InjectionTargetInfo>(
        InjectionTargetInfo{target_info->target, std::get<Group::InjectionCMode>(target_info->cmode)}
    );
}

template<class Scalar>
std::optional<typename GroupTargetCalculator<Scalar>::ProductionTargetInfo>
GroupTargetCalculator<Scalar>::
groupProductionTarget(const Group& group) {
    GeneralCalculator general_calculator{*this, group};
    auto target_info = general_calculator.calculateGroupTarget();
    if (!target_info) {
        // No target was calculated, return an empty optional.
        return std::nullopt;
    }
    return std::make_optional<ProductionTargetInfo>(
        ProductionTargetInfo{
            target_info->target,
            std::get<Group::ProductionCMode>(target_info->cmode)
        }
    );
}

// ----------------------------------------------------
// Constructor for inner class GeneralCalculator
// ----------------------------------------------------

template<class Scalar>
GroupTargetCalculator<Scalar>::
GeneralCalculator::
GeneralCalculator(
    GroupTargetCalculator<Scalar>& parent_calculator,
    const Group& original_group,
    std::optional<Phase> injection_phase  // Only used for injectors
) :
    parent_calculator_{parent_calculator},
    original_group_{original_group},
    injection_phase_{injection_phase},
    resv_coeffs_prod_(this->phaseUsage().num_phases, 0.0)
{
    this->wellModel().calcResvCoeff(
        this->fipnum(),
        this->pvtreg(),
        this->groupState().production_rates(original_group.name()),
        this->resv_coeffs_prod_
    );
}

// -------------------------------------------------------
// Public methods for the GeneralCalculator class
// -------------------------------------------------------

template<class Scalar>
std::optional<typename GroupTargetCalculator<Scalar>::TargetInfo>
GroupTargetCalculator<Scalar>::
GeneralCalculator::
calculateGroupTarget()
{
    const auto& group = this->original_group_;  // The bottom group we want to calculate the target for.
    const auto efficiency_factor = group.getGroupEfficiencyFactor();
    return this->calculateGroupTargetRecursive_(group, efficiency_factor);
}

// -------------------------------------------------------
// Private methods for the GeneralCalculator class
// -------------------------------------------------------


// Calculates the target for a group. If the group is a production group, it is
//   assumed that Group::ExceedAction::RATE (GCONPROD item 7 is "RATE").
// The group control mode can be:
// (a) FLD (then item 8 of GCONPROD/GCONINJE is ignored), but a guide rate must be defined
//   in item 9 and 10 of GCONPROD/GCONINJE.
// (b) NONE (then item 8 of GCONPROD/GCONINJE must be YES and a guide rate must be defined.
// - ORAT, WRAT, GRAT, LRAT, CRAT, RESV,... and item 8 of GCONPROD/GCONINJE is "NO", then
//   the target (item 3, 4, 5 or 6) of GCONPROD is used, or for an injector the target defined
//   in GCONINJE is used.
// (c) ORAT, WRAT, GRAT, LRAT, CRAT, RESV,... and item 8 of GCONPROD/GCONINJE is "YES", then
//   if a higher level group target is found, and a guide rate is defined, then this target is used
//   instead of the target defined in GCONPROD/GCONINJE.
//  
// - If the group has control different from FLD or NONE:
//   - If the group still is available for control (item 8 of GCONPROD/GCONINJE is "YES" or):
//     - Use procedure A
//   - else:
//     -

// with control FLD or control NONE. The group is assumed to
// have a guide rate defined. Also, a chain of parent groups with similar controls must exist until
// a group with a control mode different from FLD or NONE (e.g. ORAT) is reached. The corresponding
// target of this group is then distributed down the chain of parent groups by subtracting local
// reductions and multiplying by local fractions. The target is then returned.
template<class Scalar>
std::optional<typename GroupTargetCalculator<Scalar>::TargetInfo>
GroupTargetCalculator<Scalar>::
GeneralCalculator::
calculateGroupTargetRecursive_(const Group& group, const Scalar efficiency_factor)
{
    if (this->hasFldOrNoneControl_(group)) {
        if (!this->parentGroupControlAvailable_(group)) {
            // If we are here, the group has control mode NONE, and either
            //   - item 8 of GCONPROD/GCONINJE is "NO" or,
            //   - group.is_field() is true.
            // NOTE: it cannot have control mode FLD, since FLD will ovverride item 8 of
            //   GCONPROD/GCONINJE, see GroupKeywordHandlers.cpp
            return std::nullopt;
        } else {
            auto new_efficiency_factor = efficiency_factor * group.getGroupEfficiencyFactor();
            return this->calculateGroupTargetRecursive_(this->parentGroup(group), new_efficiency_factor);
        }
    }

    if ((this->isInjector() && !group.isInjectionGroup()) ||
            (this->isProducer() && !group.isProductionGroup()))
    {
        // Controlling group that is not of the type we are interested in.
        // TODO: When does this happen?
        return std::nullopt;
    }

    // If we are here, we are at the topmost group to be visited in the recursion.
    // This is the group containing the control we will check against.

    TopToBottomCalculator top_bottom_calc{*this, group, efficiency_factor};
    return top_bottom_calc.calculateGroupTarget();
}

template<class Scalar>
bool
GroupTargetCalculator<Scalar>::
GeneralCalculator::
hasFldOrNoneControl_(const Group& group) const
{
    // Check if the group has control FLD or NONE.
    const auto& name = group.name();
    if (this->isInjector()) {
        return this->groupState().has_field_or_none_control(name, this->injection_phase_.value());
    }
    else {
        return this->groupState().has_field_or_none_control(name);
    }
}

template<class Scalar>
bool
GroupTargetCalculator<Scalar>::
GeneralCalculator::
parentGroupControlAvailable_(const Group& group) const
{
    if (this->isInjector()) {
        return group.injectionGroupControlAvailable(this->injection_phase_.value());
    } else {
        return group.productionGroupControlAvailable();
    }
}

// ----------------------------------------------------
// Constructor for inner class TopToBottomCalculator
// ----------------------------------------------------

template<class Scalar>
GroupTargetCalculator<Scalar>::
TopToBottomCalculator::
TopToBottomCalculator(
    GeneralCalculator& parent_calculator,
    const Group& group,
    Scalar efficiency_factor
) :
    parent_calculator_{parent_calculator},
    group_{group},
    efficiency_factor_{efficiency_factor}
{
    if (this->isInjector()) {
        this->initForInjector_();
    }
    else {
        this->initForProducer_();
    }
}

// -------------------------------------------------------
// Public methods for the TopToBottomCalculator class
// -------------------------------------------------------

template<class Scalar>
std::optional<typename GroupTargetCalculator<Scalar>::TargetInfo>
GroupTargetCalculator<Scalar>::
TopToBottomCalculator::
calculateGroupTarget()
{
    auto orig_target = this->getTopLevelTarget_();
    const auto chain = this->getGroupChainTopBot_();
    // Because the bottom group (the original group) is the last of the elements,
    //   and not an ancestor, we subtract one:
    const std::size_t num_ancestors = chain.size() - 1;
    Scalar target = orig_target;
    for (std::size_t i = 0; i < num_ancestors; ++i) {
        if ((i == 0) || this->guideRate().has(chain[i])) {
            // Apply local reductions only at the control level
            // (top) and for levels where we have a specified
            // group guide rate.
            target -= this->localReduction_(chain[i]);
        }
        target *= this->localFraction_(chain[i + 1]);
    }
    // Avoid negative target rates coming from too large local reductions.
    target = std::max(Scalar(0.0), target / this->efficiency_factor_);
    return TargetInfo{target, this->toplevel_control_mode_};
}

// -------------------------------------------------------
// Private methods for the TopToBottomCalculator class
// -------------------------------------------------------

template<class Scalar>
void
GroupTargetCalculator<Scalar>::
TopToBottomCalculator::
initForInjector_()
{
    const auto& toplevel_group_control_mode = this->groupState().injection_control(
        this->group_.name(), this->injectionPhase()
    );
    // Save the toplevel control mode for later use.
    this->toplevel_control_mode_ = toplevel_group_control_mode;
    Scalar sales_target = this->getGratSalesInjectionTarget_();
    bool has_gpmaint_control = this->group_.has_gpmaint_control(
        this->injectionPhase(), toplevel_group_control_mode
    );
    this->target_calculator_.template emplace<InjectionTargetCalculator>(
        toplevel_group_control_mode,
        this->phaseUsage(),
        this->resvCoeffsInj(),
        this->group_.name(),
        sales_target,
        this->groupState(),
        this->injectionPhase(),
        has_gpmaint_control,
        this->deferredLogger()
    );
    auto guide_target_mode =
        std::get<InjectionTargetCalculator>(this->target_calculator_).guideTargetMode();
    this->fraction_calculator_.emplace(
        this->schedule(),
        this->wellState(),
        this->groupState(),
        this->summaryState(),
        this->reportStepIdx(),
        &this->guideRate(),
        GuideRateModel::Target{guide_target_mode},
        this->phaseUsage(),
        /*is_injector=*/true,
        this->injectionPhase()
    );
}

template<class Scalar>
void
GroupTargetCalculator<Scalar>::
TopToBottomCalculator::
initForProducer_()
{
    const auto& toplevel_group_control_mode =
        this->groupState().production_control(this->group_.name());
    // Save the toplevel control mode for later use.
    this->toplevel_control_mode_ = toplevel_group_control_mode;
    Scalar grat_sales_target = this->getGratSalesProductionTarget_();
    bool has_gpmaint_control = this->group_.has_gpmaint_control(toplevel_group_control_mode);
    this->target_calculator_.template emplace<TargetCalculator>(
        toplevel_group_control_mode,
        this->phaseUsage(),
        this->resvCoeffsProd(),
        grat_sales_target,
        this->group_.name(),
        this->groupState(),
        has_gpmaint_control
    );
    auto guide_target_mode =
        std::get<TargetCalculator>(this->target_calculator_).guideTargetMode();
    auto dummy_phase = Phase::OIL; // Dummy phase, not used for producers.
    this->fraction_calculator_.emplace(
        this->schedule(),
        this->wellState(),
        this->groupState(),
        this->summaryState(),
        this->reportStepIdx(),
        &this->guideRate(),
        GuideRateModel::Target{guide_target_mode},
        this->phaseUsage(),
        /*is_injector=*/true,
        dummy_phase
    );
}


template<class Scalar>
Scalar
GroupTargetCalculator<Scalar>::
TopToBottomCalculator::
getGratSalesInjectionTarget_() const
{
    const auto& gconsale = this->gconsale();
    const auto& group_name = this->group_.name();
    if (gconsale.has(group_name)) {
        return gconsale.get(group_name, this->summaryState()).sales_target;
    }
    return 0.0;

}

template<class Scalar>
Scalar
GroupTargetCalculator<Scalar>::
TopToBottomCalculator::
getGratSalesProductionTarget_() const
{
    if (this->groupState().has_grat_sales_target(this->group_.name()))
        return this->groupState().grat_sales_target(this->group_.name());
    return 0.0;

}

template<class Scalar>
std::vector<std::string>
GroupTargetCalculator<Scalar>::
TopToBottomCalculator::
getGroupChainTopBot_() const
{
    return WellGroupHelpers<Scalar>::groupChainTopBot(
        this->bottomGroup().name(), this->group_.name(), this->schedule(), this->reportStepIdx()
    );
}

template<class Scalar>
Scalar
GroupTargetCalculator<Scalar>::
TopToBottomCalculator::
getTopLevelTarget_()
{
    if (this->isInjector()) {
        auto control_mode = this->groupState().injection_control(
            this->group_.name(), this->injectionPhase()
        );
        std::optional<Group::InjectionControls> ctrl;
        if (!this->group_.has_gpmaint_control(this->injectionPhase(), control_mode))
            ctrl = this->group_.injectionControls(this->injectionPhase(), this->summaryState());
        return std::get<InjectionTargetCalculator>(this->target_calculator_).groupTarget(
            ctrl, this->deferredLogger()
        );
    }
    else {
        auto control_mode = this->groupState().production_control(this->group_.name());
        std::optional<Group::ProductionControls> ctrl;
        if (!this->group_.has_gpmaint_control(control_mode))
            ctrl = this->group_.productionControls(this->summaryState());
        return std::get<TargetCalculator>(this->target_calculator_).groupTarget(
            ctrl, this->deferredLogger()
        );
    }
}

template<class Scalar>
Scalar
GroupTargetCalculator<Scalar>::
TopToBottomCalculator::
localFraction_(const std::string& group_name)
{
    const auto& always_included_child = this->bottomGroup().name();
    return this->fraction_calculator_->localFraction(group_name, always_included_child);
}

template<class Scalar>
Scalar
GroupTargetCalculator<Scalar>::
TopToBottomCalculator::
localReduction_(const std::string& group_name)
{
    if (this->isInjector()) {
        const std::vector<Scalar>& group_target_reductions =
            this->groupState().injection_reduction_rates(group_name);
        return std::get<InjectionTargetCalculator>(this->target_calculator_).calcModeRateFromRates(
            group_target_reductions
        );
    }
    else {
        const std::vector<Scalar>& group_target_reductions =
            this->groupState().production_reduction_rates(group_name);
        return std::get<TargetCalculator>(this->target_calculator_).calcModeRateFromRates(
            group_target_reductions
        );

    }
}


template class GroupTargetCalculator<double>;

#if FLOW_INSTANTIATE_FLOAT
template class GroupTargetCalculator<float>;
#endif

} // namespace Opm
