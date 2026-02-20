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
#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>

#include <tuple>  // for std::tie
#include <fmt/format.h>

namespace Opm {

// ----------------------------------------------------
// Constructor for the GroupTargetCalculator class
// ----------------------------------------------------
template<class Scalar, class IndexTraits>
GroupTargetCalculator<Scalar, IndexTraits>::
GroupTargetCalculator(
    const BlackoilWellModelGeneric<Scalar, IndexTraits>& well_model,
    const GroupStateHelperType& group_state_helper
) :
    well_model_{well_model},
    group_state_helper_{group_state_helper},
    well_state_{group_state_helper.wellState()},
    group_state_{group_state_helper.groupState()},
    schedule_{group_state_helper.schedule()},
    summary_state_{group_state_helper.summaryState()},
    phase_usage_{group_state_helper.phaseUsage()},
    guide_rate_{group_state_helper.guideRate()},
    report_step_idx_{group_state_helper.reportStepIdx()},
    resv_coeffs_inj_(group_state_helper.phaseUsage().numPhases, 0.0)
{
    std::tie(this->fipnum_, this->pvtreg_) = this->well_model_.getGroupFipnumAndPvtreg();
    // Get the reservoir coefficients for injection
    // NOTE: The production reservoir coefficients depend on the group production rates, so we
    //  cannot precompute them here.
    this->well_model_.calcInjResvCoeff(this->fipnum_, this->pvtreg_, this->resv_coeffs_inj_);
}

template<class Scalar, class IndexTraits>
std::optional<typename GroupTargetCalculator<Scalar, IndexTraits>::InjectionTargetInfo>
GroupTargetCalculator<Scalar, IndexTraits>::
groupInjectionTarget(const Group& group, ReservoirCoupling::Phase injection_phase)
{
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

template<class Scalar, class IndexTraits>
std::optional<typename GroupTargetCalculator<Scalar, IndexTraits>::ProductionTargetInfo>
GroupTargetCalculator<Scalar, IndexTraits>::
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

template<class Scalar, class IndexTraits>
GroupTargetCalculator<Scalar, IndexTraits>::
GeneralCalculator::
GeneralCalculator(
    GroupTargetCalculator<Scalar, IndexTraits>& parent_calculator,
    const Group& original_group,
    std::optional<ReservoirCoupling::Phase> injection_phase  // Only used for injectors
) :
    parent_calculator_{parent_calculator},
    original_group_{original_group},
    injection_phase_{injection_phase},
    resv_coeffs_prod_(this->phaseUsage().numPhases, 0.0)
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

template<class Scalar, class IndexTraits>
std::optional<typename GroupTargetCalculator<Scalar, IndexTraits>::TargetInfo>
GroupTargetCalculator<Scalar, IndexTraits>::
GeneralCalculator::
calculateGroupTarget()
{
    // TODO: For now this is adapted to the reservoir coupling implementation where "group" below will
    //   correspond to a master group. In the future we might want to port the logic to include
    //   e.g. checkGroupConstraintsProd() and checkGroupConstraintsInj() in GroupStateHelper.cpp.
    const auto& group = this->original_group_;  // The bottom group we want to calculate the target for.
    // For injection, check if the group has injection control defined for this phase.
    // If not, there is no target to calculate.
    if (this->targetType() == TargetType::Injection) {
        if (!group.hasInjectionControl(this->injectionPhase_())) {
            return std::nullopt;
        }
    }
    if (group.is_field() || !this->parentGroupControlAvailable_(group)) {
        return this->getGroupTargetNoGuideRate(group);
    }
    assert(this->parentGroupControlAvailable_(group));
    if (!this->hasGuideRate(group)) {
        if (this->hasFldOrNoneControl(group)) {
            // Parent control is available, but no guide rate is defined. This is illegal for a master group
            // under FLD or NONE control.
            OPM_DEFLOG_THROW(
                std::runtime_error,
                fmt::format("No guide rate defined for master group {}", group.name()),
                this->deferredLogger());
        }
        // A master group with:
        //   - individual (not FLD or NONE) control,
        //   - parent control available,
        //   - but no guide rate
        // this could be considered an error, but we can also fall back to use its own target.
        return this->getGroupTargetNoGuideRate(group);
    }
    const auto efficiency_factor = group.getGroupEfficiencyFactor();
    auto target = this->calculateGroupTargetRecursive_(this->parentGroup(group), efficiency_factor);
    if (target) {
        // TODO: We could probably switch the group to FLD control mode now. However, this call is coming
        //    from beginTimeStep() in BlackoilWellModel_impl.hpp, but the regular higher constraints
        //    checks (for all groups, not just the master group) will first be done during the assemble() step,
        //    at which point the group should be switched to FLD control mode if necessary.
        return target;
    }
    // If a higher level target was not found or not violated, use the group's own target.
    return this->getGroupTargetNoGuideRate(group);
}

template<class Scalar, class IndexTraits>
typename GroupTargetCalculator<Scalar, IndexTraits>::GeneralCalculator::TargetCalculatorType
GroupTargetCalculator<Scalar, IndexTraits>::
GeneralCalculator::
getInjectionTargetCalculator(const Group& /*group*/)
{
    return InjectionTargetCalculator{
        this->groupStateHelper(),
        this->injectionPhase_()
    };
}

template<class Scalar, class IndexTraits>
typename GroupTargetCalculator<Scalar, IndexTraits>::GeneralCalculator::TargetCalculatorType
GroupTargetCalculator<Scalar, IndexTraits>::
GeneralCalculator::
getProductionTargetCalculator(const Group& group) const
{
    return TargetCalculator{
        this->groupStateHelper(),
        this->resvCoeffsProd(),
        group
    };
}

template<class Scalar, class IndexTraits>
typename GroupTargetCalculator<Scalar, IndexTraits>::GeneralCalculator::TargetCalculatorType
GroupTargetCalculator<Scalar, IndexTraits>::
GeneralCalculator::
getTargetCalculator(const Group& group)
{
    if (this->targetType() == TargetType::Injection) {
        return this->getInjectionTargetCalculator(group);
    }
    else {
        return this->getProductionTargetCalculator(group);
    }
}

template<class Scalar, class IndexTraits>
typename GroupTargetCalculator<Scalar, IndexTraits>::TargetInfo
GroupTargetCalculator<Scalar, IndexTraits>::
GeneralCalculator::
getGroupTargetNoGuideRate(const Group& group)
{
    // Efficiency factor is not applied here, since it only applies to child rates, not to the group's
    // own target.
    // NOTE: If a master group is under individual control and has a RESV rate target, it should be
    //   transferred into a RESV target for the slave group. But it is not obvious how to transform it
    //   since the master and slave groups in general have different PVT properties.
    // TODO: For now we assume the RESV rate target is in slave reservoir units, and just send it as is
    //   to the slave group. Alternatively, we could also throw an error here for that case
    if (this->targetType() == TargetType::Injection) {
        const auto& control_mode = this->groupState().injection_control(group.name(), this->injectionPhase_());
        return TargetInfo{
            this->groupStateHelper().getInjectionGroupTarget(group, this->injectionPhase_(), this->resvCoeffsInj()),
            control_mode
        };
    }
    else {
        const auto& control_mode = this->groupState().production_control(group.name());
        return TargetInfo{
            this->groupStateHelper().getProductionGroupTarget(group),
            control_mode
        };
    }
}

template<class Scalar, class IndexTraits>
bool
GroupTargetCalculator<Scalar, IndexTraits>::
GeneralCalculator::
hasFldOrNoneControl(const Group& group)
{
    // Check if the group (constrained to production or injection) has control FLD or NONE.
    // For example, a pure injection group will have production control NONE.
    const auto& name = group.name();
    if (this->targetType() == TargetType::Injection) {
        return this->groupState().has_field_or_none_control(name, this->injectionPhase_());
    }
    else {
        return this->groupState().has_field_or_none_control(name);
    }
}


// -------------------------------------------------------
// Private methods for the GeneralCalculator class
// -------------------------------------------------------

// See comments above RescoupTargetCalculator::calculateMasterGroupTargetsAndSendToSlaves() about how
// the target is calculated.
template<class Scalar, class IndexTraits>
std::optional<typename GroupTargetCalculator<Scalar, IndexTraits>::TargetInfo>
GroupTargetCalculator<Scalar, IndexTraits>::
GeneralCalculator::
calculateGroupTargetRecursive_(const Group& group, const Scalar efficiency_factor)
{
    if (this->hasFldOrNoneControl(group)) {
        // NOTE: A pure injection group is assumed to have production control NONE and
        //   a pure production group is assumed to have injection control NONE, see
        //   GroupStateHelper.cpp::setCmodeGroup() for details.
        if (!this->parentGroupControlAvailable_(group)) {
            // If we are here, the group has control mode NONE, and either
            //   - item 8 of GCONPROD/GCONINJE is "NO" or,
            //   - group.is_field() is true.
            // NOTE: it cannot have control mode FLD, since FLD will ovverride item 8 of
            //   GCONPROD/GCONINJE, see GroupKeywordHandlers.cpp
            // TODO: We could throw an error here. Or try to fall back to use the group's own target.
            // Fall back to use the group's own target.
            return std::nullopt;
        } else {
            assert(!group.is_field());  // The field group should not have parent control.
            auto new_efficiency_factor = efficiency_factor * group.getGroupEfficiencyFactor();
            return this->calculateGroupTargetRecursive_(this->parentGroup(group), new_efficiency_factor);
        }
    }

    if ((this->targetType() == TargetType::Injection && !group.isInjectionGroup()) ||
            (this->targetType() == TargetType::Production && !group.isProductionGroup()))
    {
        // Controlling group that is not of the type we are interested in.
        // This should never happen, since any group that is not an injector should have injection
        // control NONE. And any group that is not a producer should have production control NONE.
        // See GroupStateHelper.cpp::setCmodeGroup() for details.
        // .. and therefore should be caught by the hasFldOrNoneControl() check above.
        OPM_DEFLOG_THROW(
            std::runtime_error,
            fmt::format("Controlling group that is not of the type we are interested in: {}", group.name()),
            this->deferredLogger()
        );
    }

    // If we are here, we are at the topmost group to be visited in the recursion.
    // This is the group containing the control we will check against.

    const auto& top_group = group;
    const auto& bottom_group = this->original_group_;
    TopToBottomCalculator top_bottom_calc{*this, top_group, bottom_group, efficiency_factor};
    return top_bottom_calc.calculateGroupTarget();
}

template<class Scalar, class IndexTraits>
bool
GroupTargetCalculator<Scalar, IndexTraits>::
GeneralCalculator::
parentGroupControlAvailable_(const Group& group)
{
    if (this->targetType() == TargetType::Injection) {
        return group.injectionGroupControlAvailable(this->injectionPhase_());
    } else {
        return group.productionGroupControlAvailable();
    }
}

template<class Scalar, class IndexTraits>
Phase
GroupTargetCalculator<Scalar, IndexTraits>::
GeneralCalculator::
injectionPhase_() const
{
    if (this->injection_phase_.has_value()) {
        switch (this->injection_phase_.value()) {
            case ReservoirCoupling::Phase::Water:
                return Phase::WATER;
            case ReservoirCoupling::Phase::Oil:
                return Phase::OIL;
            case ReservoirCoupling::Phase::Gas:
                return Phase::GAS;
            default:
                OPM_DEFLOG_THROW(std::runtime_error, "Invalid injection phase", this->deferredLogger());
                return Phase::WATER; // Unreachable, but satisfies compiler
        }
    }
    else {
        return Phase::OIL;   // Dummy phase, not used for producers.
    }
}

// ----------------------------------------------------
// Constructor for inner class TopToBottomCalculator
// ----------------------------------------------------

template<class Scalar, class IndexTraits>
GroupTargetCalculator<Scalar, IndexTraits>::
TopToBottomCalculator::
TopToBottomCalculator(
    GeneralCalculator& parent_calculator,
    const Group& top_group,
    const Group& bottom_group,
    Scalar chain_efficiency_factor
) :
    parent_calculator_{parent_calculator},
    top_group_{top_group},
    bottom_group_{bottom_group},
    chain_efficiency_factor_{chain_efficiency_factor}
{
    if (this->targetType() == TargetType::Injection) {
        this->initForInjector_();
    }
    else {
        this->initForProducer_();
    }
}

// -------------------------------------------------------
// Public methods for the TopToBottomCalculator class
// -------------------------------------------------------

template<class Scalar, class IndexTraits>
std::optional<typename GroupTargetCalculator<Scalar, IndexTraits>::TargetInfo>
GroupTargetCalculator<Scalar, IndexTraits>::
TopToBottomCalculator::
calculateGroupTarget()
{
    const auto chain = this->getGroupChainTopBot_();
    // NOTE: For reservoir coupling: we are called from beginTimeStep() in BlackoilWellModel_impl.hpp,
    //   so we do not check for guide rate violations here, compare with
    //  GroupStateHelper.cpp::checkGroupConstraintsProd() and
    //  GroupStateHelper.cpp::checkGroupConstraintsInj() since the guide rate will be checked during the
    //   assemble() step.
    const auto local_reduction_level = this->getLocalReductionLevel_(chain);
    Scalar top_level_target = this->getTopLevelTarget_();
    Scalar bottom_group_current_rate_available = this->getBottomGroupCurrentRateAvailable_();
    // Because the bottom group (the original group) is the last of the elements,
    //   and not an ancestor, we subtract one:
    const std::size_t num_ancestors = chain.size() - 1;
    Scalar target = top_level_target;
    for (std::size_t i = 0; i < num_ancestors; ++i) {
        if ((i == 0) || this->hasGuideRate_(chain[i])) {
            // Apply local reductions only at the control level
            // (top) and for levels where we have a specified
            // group guide rate.
            if (i <= local_reduction_level) {
                target -= this->localReduction_(chain[i]);
            }
            // Add my reduction back at the level where it is included in the local reduction
            if (this->bottomGroupHasIndividualControl_() && i == local_reduction_level) {
                const Scalar addback_efficiency
                    = this->computeAddbackEfficiency_(chain, local_reduction_level);
                target += bottom_group_current_rate_available * addback_efficiency;
            }
        }
        target *= this->localFraction_(chain[i + 1]);
    }
    // Avoid negative target rates coming from too large local reductions.
    target = std::max(Scalar(0.0), target);
    // Divide by the accumulated efficiency factor along the chain from top to bottom, excluding the top group.
    // This is the full target rate that should be assigned to the bottom group, to be seen as the unscaled
    // target rate for the top group.
    const Scalar full_target = target / this->chain_efficiency_factor_;
    if (bottom_group_current_rate_available > TARGET_RATE_TOLERANCE) {  // > 1e-12
        const auto toplevel_control_mode = this->getToplevelControlMode_();
        if ((bottom_group_current_rate_available > full_target)) {
            // The bottom group is producing too much according to the top level group target,
            // so we need to scale down the target rate.
            if (this->isProducerAndRESVControl_(this->top_group_)) {
                // For RESV mode with reservoir coupling master groups, we scale the slave's
                // actual reservoir rate to get a target in slave's reservoir units.
                const Scalar slave_resv_rate = this->getSlaveGroupReservoirRate_(this->bottom_group_);
                // Scale the slave's reservoir rate by the same factor we would scale surface rates
                const Scalar scale = full_target / bottom_group_current_rate_available;
                return TargetInfo{scale * slave_resv_rate, toplevel_control_mode};
            }
            else {
                return TargetInfo{full_target, toplevel_control_mode};
            }
        }
        else if (this->hasFLDControl_(this->bottom_group_)) {
            // The bottom group is under FLD control, so we return the top level target as is.
            return TargetInfo{full_target, toplevel_control_mode};
        }
    }
    // Return the bottom group's target rate as is.
    return this->getGroupTargetNoGuideRate_(this->bottom_group_);
}

// -------------------------------------------------------
// Private methods for the TopToBottomCalculator class
// -------------------------------------------------------

template <class Scalar, class IndexTraits>
bool
GroupTargetCalculator<Scalar, IndexTraits>::
TopToBottomCalculator::
bottomGroupHasIndividualControl_()
{
    return !this->hasFldOrNoneControl_(this->bottom_group_);
}

template <class Scalar, class IndexTraits>
Scalar
GroupTargetCalculator<Scalar, IndexTraits>::
TopToBottomCalculator::
computeAddbackEfficiency_(
    const std::vector<std::string>& chain,
    const std::size_t local_reduction_level) const
{
    // Compute partial efficiency factor from local_reduction_level down to the entity.
    // Chain is ordered [control_group, ..., local_reduction_level, ..., entity].
    // We multiply efficiency factors from index (local_reduction_level + 1) to the entity.
    // Note: The entity at chain[num_ancestors] may be a well or a group.
    // Both wells (WEFAC) and groups (GEFAC) have efficiency factors that must be included,
    // since the local_reduction rates include these efficiency factors.
    const std::size_t num_ancestors = chain.size() - 1;
    Scalar efficiency = 1.0;
    for (std::size_t jj = local_reduction_level + 1; jj <= num_ancestors; ++jj) {
        const std::string& name = chain[jj];
        if (this->schedule().hasGroup(name, this->reportStepIdx())) {
            const auto& grp = this->schedule().getGroup(name, this->reportStepIdx());
            efficiency *= grp.getGroupEfficiencyFactor();
        } else if (this->schedule().hasWell(name, this->reportStepIdx())) {
            const auto& well = this->schedule().getWell(name, this->reportStepIdx());
            efficiency *= well.getEfficiencyFactor()
                * this->wellState()[name].efficiency_scaling_factor;
        }
    }
    return efficiency;
}

template<class Scalar, class IndexTraits>
Scalar
GroupTargetCalculator<Scalar, IndexTraits>::
TopToBottomCalculator::
getBottomGroupCurrentRateAvailable_() const
{
    const auto& group = this->bottom_group_;
    if (this->targetType() == TargetType::Injection) {
        return std::get<InjectionTargetCalculator>(this->target_calculator_).calcModeRateFromRates(
            this->groupStateHelper().getGroupRatesAvailableForHigherLevelControl(group, /*is_injector=*/true)
        );
    }
    else {
        const auto& tcalc = std::get<TargetCalculator>(this->target_calculator_);
        return tcalc.calcModeRateFromRates(
            this->groupStateHelper().getGroupRatesAvailableForHigherLevelControl(group, /*is_injector=*/false)
        );
    }
}

template<class Scalar, class IndexTraits>
Scalar
GroupTargetCalculator<Scalar, IndexTraits>::
TopToBottomCalculator::
getSlaveGroupReservoirRate_(const Group& group)
{
    if (!this->groupStateHelper().isReservoirCouplingMasterGroup(group)) {
        OPM_DEFLOG_THROW(std::runtime_error, "Group is not a reservoir coupling master group", this->deferredLogger());
    }
    // Sum reservoir rates for all phases from the slave
    Scalar total_resv_rate = 0.0;
    const auto& master = this->groupStateHelper().reservoirCouplingMaster();
    for (auto phase : {ReservoirCoupling::Phase::Oil, ReservoirCoupling::Phase::Gas, ReservoirCoupling::Phase::Water}) {
        total_resv_rate += master.getMasterGroupRate(group.name(), phase, ReservoirCoupling::RateKind::ProductionReservoir);
    }
    return total_resv_rate;
}

template<class Scalar, class IndexTraits>
std::vector<std::string>
GroupTargetCalculator<Scalar, IndexTraits>::
TopToBottomCalculator::
getGroupChainTopBot_() const
{
    return this->groupStateHelper().groupChainTopBot(this->bottom_group_.name(), this->top_group_.name());
}

template<class Scalar, class IndexTraits>
std::size_t
GroupTargetCalculator<Scalar, IndexTraits>::
TopToBottomCalculator::
getLocalReductionLevel_(const std::vector<std::string>& chain)
{
    const std::size_t num_ancestors = chain.size() - 1;
    std::size_t local_reduction_level = 0;
    for (std::size_t ii = 1; ii < num_ancestors; ++ii) {
        const int num_gr_ctrl = this->groupStateHelper().groupControlledWells(
            chain[ii],
            /*always_included_child=*/"",
            /*is_producer=*/this->targetType() == TargetType::Production,
            /*injection_phase=*/this->injectionPhase_());
        if (this->hasGuideRate_(chain[ii]) && num_gr_ctrl > 0) {
            local_reduction_level = ii;
        }
    }
    return local_reduction_level;
}

template<class Scalar, class IndexTraits>
typename GroupTargetCalculator<Scalar, IndexTraits>::ControlMode
GroupTargetCalculator<Scalar, IndexTraits>::
TopToBottomCalculator::
getToplevelControlMode_() const
{
    if (this->targetType() == TargetType::Injection) {
        return this->groupState().injection_control(this->top_group_.name(), this->injectionPhase_());
    }
    else {
        return this->groupState().production_control(this->top_group_.name());
    }
}

template<class Scalar, class IndexTraits>
Scalar
GroupTargetCalculator<Scalar, IndexTraits>::
TopToBottomCalculator::
getTopLevelTarget_()
{
    if (this->targetType() == TargetType::Injection) {
        return this->groupStateHelper().getInjectionGroupTarget(
            this->top_group_, this->injectionPhase_(), this->resvCoeffsInj());
    }
    else {
        return this->groupStateHelper().getProductionGroupTarget(this->top_group_);
    }
}

template<class Scalar, class IndexTraits>
bool
GroupTargetCalculator<Scalar, IndexTraits>::
TopToBottomCalculator::
hasFLDControl_(const Group& group) const
{
    if (this->targetType() == TargetType::Injection) {
        return this->groupState().injection_control(group.name(), this->injectionPhase_()) == Group::InjectionCMode::FLD;
    }
    else {
        return this->groupState().production_control(group.name()) == Group::ProductionCMode::FLD;
    }
}

template<class Scalar, class IndexTraits>
void
GroupTargetCalculator<Scalar, IndexTraits>::
TopToBottomCalculator::
initForInjector_()
{
    this->target_calculator_.template emplace<InjectionTargetCalculator>(
        this->groupStateHelper(),
        this->injectionPhase_()
    );
    const auto guide_target_mode = this->groupStateHelper().getInjectionGuideTargetMode(this->injectionPhase_());
    this->fraction_calculator_.emplace(
        this->schedule(),
        this->groupStateHelper(),
        this->summaryState(),
        this->reportStepIdx(),
        &this->guideRate(),
        GuideRateModel::Target{guide_target_mode},
        /*is_producer=*/false,
        this->injectionPhase_()
    );
}

template<class Scalar, class IndexTraits>
void
GroupTargetCalculator<Scalar, IndexTraits>::
TopToBottomCalculator::
initForProducer_()
{
    this->target_calculator_.template emplace<TargetCalculator>(
        this->groupStateHelper(),
        this->resvCoeffsProd(),
        this->top_group_
    );
    const auto guide_target_mode = this->groupStateHelper().getProductionGuideTargetMode(this->top_group_);
    const auto dummy_phase = Phase::OIL; // Dummy phase, not used for producers.
    this->fraction_calculator_.emplace(
        this->schedule(),
        this->groupStateHelper(),
        this->summaryState(),
        this->reportStepIdx(),
        &this->guideRate(),
        GuideRateModel::Target{guide_target_mode},
        /*is_producer=*/true,
        dummy_phase
    );
}

template<class Scalar, class IndexTraits>
bool
GroupTargetCalculator<Scalar, IndexTraits>::
TopToBottomCalculator::
isProducerAndRESVControl_(const Group& group) const
{
    if (this->targetType() == TargetType::Production) {
        return this->groupState().production_control(group.name()) == Group::ProductionCMode::RESV;
    }
    else {
        return false;
    }
}

template<class Scalar, class IndexTraits>
Scalar
GroupTargetCalculator<Scalar, IndexTraits>::
TopToBottomCalculator::
localFraction_(const std::string& group_name)
{
    const auto& always_included_child = this->bottom_group_.name();
    return this->fraction_calculator_->localFraction(group_name, always_included_child);
}

template<class Scalar, class IndexTraits>
Scalar
GroupTargetCalculator<Scalar, IndexTraits>::
TopToBottomCalculator::
localReduction_(const std::string& group_name)
{
    if (this->targetType() == TargetType::Injection) {
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


template class GroupTargetCalculator<double, BlackOilDefaultFluidSystemIndices>;

#if FLOW_INSTANTIATE_FLOAT
template class GroupTargetCalculator<float, BlackOilDefaultFluidSystemIndices>;
#endif

} // namespace Opm
