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
#include <opm/input/eclipse/Schedule/ResCoup/ReservoirCouplingInfo.hpp>
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
                                                      const PhaseUsageInfo<IndexTraits>& phase_usage_info,
                                                      const Parallel::Communication& comm,
                                                      bool terminal_output)
    : well_state_ {&well_state}
    , group_state_ {&group_state}
    , schedule_ {schedule}
    , summary_state_ {summary_state}
    , guide_rate_ {guide_rate}
    , phase_usage_info_ {phase_usage_info}
    , comm_ {comm}
    , terminal_output_ {terminal_output}
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
                                                               const bool check_guide_rate) const
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
                                              check_guide_rate);
    }

    // This can be false for FLD-controlled groups, we must therefore
    // check for FLD first (done above).
    if (!group.isInjectionGroup()) {
        return std::make_pair(false, Scalar(1.0));
    }

    // If we are here, we are at the topmost group to be visited in the recursion.
    // This is the group containing the control we will check against.
    GroupStateHelpers::InjectionTargetCalculator<Scalar, IndexTraits> tcalc {*this,
                                                                             injection_phase};

    GroupStateHelpers::FractionCalculator fcalc {this->schedule_,
                                                 *this,
                                                 this->summary_state_,
                                                 this->report_step_,
                                                 &this->guide_rate_,
                                                 this->getInjectionGuideTargetMode(injection_phase),
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

    const Scalar orig_target = this->getInjectionGroupTarget(group, injection_phase, resv_coeff);
    const Scalar current_rate_available = tcalc.calcModeRateFromRates(rates);
    const std::size_t local_reduction_level = this->getLocalReductionLevel_(
        chain, /*is_production_group=*/false, injection_phase);

    const Scalar target = this->applyReductionsAndFractions_(chain,
                                                             orig_target,
                                                             current_rate_available,
                                                             local_reduction_level,
                                                             /*is_production_group=*/false,
                                                             injection_phase,
                                                             local_reduction_lambda,
                                                             local_fraction_lambda,
                                                             /*do_addback=*/true);

    // Avoid negative target rates coming from too large local reductions.
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
                                                                const bool check_guide_rate) const
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
                                               check_guide_rate);
    }

    // This can be false for FLD-controlled groups, we must therefore
    // check for FLD first (done above).
    if (!group.isProductionGroup()) {
        return std::make_pair(false, Scalar(1.0));
    }

    // If we are here, we are at the topmost group to be visited in the recursion.
    // This is the group containing the control we will check against.

    GroupStateHelpers::TargetCalculator<Scalar, IndexTraits> tcalc {*this,
                                                                    resv_coeff,
                                                                    group};

    GroupStateHelpers::FractionCalculator<Scalar, IndexTraits> fcalc {this->schedule_,
                                                                      *this,
                                                                      this->summary_state_,
                                                                      this->report_step_,
                                                                      &this->guide_rate_,
                                                                      this->getProductionGuideTargetMode(group),
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

    const Scalar orig_target = this->getProductionGroupTarget(group);
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

    const std::size_t local_reduction_level = this->getLocalReductionLevel_(
        chain, /*is_production_group=*/true, /*injection_phase=*/Phase::OIL);

    const Scalar target = this->applyReductionsAndFractions_(chain,
                                                             orig_target,
                                                             current_rate_available,
                                                             local_reduction_level,
                                                             /*is_production_group=*/true,
                                                             /*injection_phase=*/Phase::OIL,
                                                             local_reduction_lambda,
                                                             local_fraction_lambda,
                                                             /*do_addback=*/true);

    // Avoid negative target rates coming from too large local reductions.
    const Scalar target_rate_available = std::max(Scalar(1e-12), target / efficiency_factor);

    Scalar scale = 1.0;
    if (current_rate_available > 1e-12)
        scale = target_rate_available / current_rate_available;

    return std::make_pair(current_rate_available > target_rate_available, scale);
}

template<typename Scalar, typename IndexTraits>
std::pair<Group::ProductionCMode, Scalar>
GroupStateHelper<Scalar, IndexTraits>::
checkGroupProductionConstraints(const Group& group) const
{
    const auto controls = group.productionControls(this->summary_state_);
    const auto currentControl = this->groupState().production_control(group.name());

    for (const auto cmode : {
        Group::ProductionCMode::ORAT,
        Group::ProductionCMode::WRAT,
        Group::ProductionCMode::GRAT,
        Group::ProductionCMode::LRAT,
        Group::ProductionCMode::RESV})
    {
        if (!group.has_control(cmode) || currentControl == cmode) {
            continue;
        }

        Scalar current_rate = this->sumProductionRate_(group, cmode);
        Scalar target = this->getProductionConstraintTarget_(group, cmode, controls);

        // LRAT skip heuristic: if liquid and oil targets are equal
        // and water rate is ~0, skip the LRAT check.
        if (cmode == Group::ProductionCMode::LRAT
            && target == controls.oil_target)
        {
            const auto& pu = this->phaseUsage();
            Scalar water_rate = this->sumWellSurfaceRates(group,
                pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx),
                false);
            water_rate = this->comm().sum(water_rate);
            if (std::abs(water_rate) < 1e-12) {
                this->deferredLogger().debug(
                    "LRAT_ORAT_GROUP",
                    "GROUP " + group.name()
                    + " The LRAT target is equal the ORAT target"
                      " and the water rate is zero, skip checking LRAT");
                continue;
            }
        }

        auto result = this->checkProductionRateConstraint_(
            group, cmode, currentControl, target, current_rate);
        if (result.first != Group::ProductionCMode::NONE) {
            return result;
        }
    }

    if (group.has_control(Group::ProductionCMode::CRAT)) {
        OPM_DEFLOG_THROW(std::runtime_error,
            "Group " + group.name()
            + "CRAT control for production groups not implemented",
            this->deferredLogger());
    }
    if (group.has_control(Group::ProductionCMode::PRBL)) {
        OPM_DEFLOG_THROW(std::runtime_error,
            "Group " + group.name()
            + "PRBL control for production groups not implemented",
            this->deferredLogger());
    }

    return {Group::ProductionCMode::NONE, Scalar(1.0)};
}

template <typename Scalar, typename IndexTraits>
Scalar
GroupStateHelper<Scalar, IndexTraits>::
getInjectionGroupTarget(
    const Group& group,
    const Phase& injection_phase,
    const std::vector<Scalar>& resv_coeff) const
{
    const auto cmode = this->groupState().injection_control(group.name(), injection_phase);

    // Check for master target override (reservoir coupling slave)
#ifdef RESERVOIR_COUPLING_ENABLED
    if (this->isReservoirCouplingSlave()
        && this->reservoirCouplingSlave().hasMasterInjectionTarget(group.name(), injection_phase))
    {
        auto [master_target, master_cmode] =
            this->reservoirCouplingSlave().masterInjectionTarget(group.name(), injection_phase);
        auto filter = this->getInjectionFilterFlag_(group.name(), injection_phase);
        using FilterFlag = ReservoirCoupling::GrupSlav::FilterFlag;
        if (filter == FilterFlag::MAST) {
            return master_target;
        }
        if (filter == FilterFlag::BOTH) {
            Scalar slave_target = this->getInjectionGroupTargetForMode_(
                group, injection_phase, resv_coeff, cmode);
            return std::min(master_target, slave_target);
        }
        // FilterFlag::SLAV: fall through to original (non-RC) logic below
    }
#endif

    return this->getInjectionGroupTargetForMode_(group, injection_phase, resv_coeff, cmode);
}

template<typename Scalar, typename IndexTraits>
Scalar
GroupStateHelper<Scalar, IndexTraits>::
getInjectionGroupTargetForMode_(
    const Group& group,
    const Phase& injection_phase,
    const std::vector<Scalar>& resv_coeff,
    Group::InjectionCMode cmode) const
{
    const auto& pu = this->phaseUsage();
    const int pos = this->phaseToActivePhaseIdx(injection_phase);
    Group::InjectionControls ctrl = group.injectionControls(injection_phase, this->summary_state_);
    bool use_gpmaint = group.has_gpmaint_control(injection_phase, cmode)
                                && this->groupState().has_gpmaint_target(group.name());
    switch (cmode) {
    case Group::InjectionCMode::RATE:
        if (use_gpmaint) {
            return this->groupState().gpmaint_target(group.name());
        }
        return ctrl.surface_max_rate;
    case Group::InjectionCMode::RESV:
        if (use_gpmaint)
            return this->groupState().gpmaint_target(group.name()) / resv_coeff[pos];

        return ctrl.resv_max_rate / resv_coeff[pos];
    case Group::InjectionCMode::REIN: {
        Scalar production_rate = this->groupState().injection_rein_rates(ctrl.reinj_group)[pos];
        return ctrl.target_reinj_fraction * production_rate;
    }
    case Group::InjectionCMode::VREP: {
        // We use the injection_reservoir_rates directly instead of the reduction rates here to account for the
        // possibility that the group in question has both a VREP control and another injection control for a different phase.
        const std::vector<Scalar>& group_injection_reservoir_rates =
                                this->groupState().injection_reservoir_rates(group.name());
        Scalar voidage_rate = this->groupState().injection_vrep_rate(ctrl.voidage_group) * ctrl.target_void_fraction;
        if (ctrl.phase != Phase::WATER && pu.phaseIsActive(IndexTraits::waterPhaseIdx)) {
            const int water_pos = pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx);
            voidage_rate -= group_injection_reservoir_rates[water_pos];
        }
        if (ctrl.phase != Phase::OIL && pu.phaseIsActive(IndexTraits::oilPhaseIdx)) {
            const int oil_pos = pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
            voidage_rate -= group_injection_reservoir_rates[oil_pos];
        }
        if (ctrl.phase != Phase::GAS && pu.phaseIsActive(IndexTraits::gasPhaseIdx)) {
            const int gas_pos = pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx);
            voidage_rate -= group_injection_reservoir_rates[gas_pos];
        }
        return voidage_rate / resv_coeff[pos];
    }
    case Group::InjectionCMode::SALE: {
        assert(pos == pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx) );
        Scalar sales_target = 0;
        if (this->schedule_[this->report_step_].gconsale().has(group.name())) {
            const auto& gconsale
                = this->schedule_[this->report_step_].gconsale().get(group.name(), this->summary_state_);
            sales_target = gconsale.sales_target;
        }
        // Gas injection rate = Total gas production rate + gas import rate - gas consumption rate - sales rate;
        // Gas import and consumption is already included in the REIN rates
        Scalar inj_rate = this->groupState().injection_rein_rates(group.name())[pos];
        inj_rate -= sales_target;
        return inj_rate;
    }
    default:
        OPM_DEFLOG_THROW(std::logic_error,
                         "Invalid Group::InjectionCMode in getInjectionGroupTargetForMode_",
                         this->deferredLogger());
        return 0.0;
    }
}

template <typename Scalar, typename IndexTraits>
GuideRate::RateVector
GroupStateHelper<Scalar, IndexTraits>::getProductionGroupRateVector(const std::string& group_name) const
{
    return this->getGuideRateVector_(this->groupState().production_rates(group_name));
}

template<typename Scalar, typename IndexTraits>
Scalar
GroupStateHelper<Scalar, IndexTraits>::
getProductionGroupTarget(const Group& group) const
{
    const auto cmode = this->groupState().production_control(group.name());

    // Check for master target override (reservoir coupling slave)
#ifdef RESERVOIR_COUPLING_ENABLED
    if (this->isReservoirCouplingSlave()
        && this->reservoirCouplingSlave().hasMasterProductionTarget(group.name()))
    {
        auto [master_target, master_cmode] =
            this->reservoirCouplingSlave().masterProductionTarget(group.name());
        auto filter = this->getProductionFilterFlag_(group.name(), master_cmode);
        using FilterFlag = ReservoirCoupling::GrupSlav::FilterFlag;
        if (filter == FilterFlag::MAST) {
            return master_target;
        }
        if (filter == FilterFlag::BOTH) {
            Scalar slave_target = this->getProductionGroupTargetForMode_(group, cmode);
            return std::min(master_target, slave_target);
        }
        // FilterFlag::SLAV: fall through to original (non-RC) logic below
    }
#endif

    return this->getProductionGroupTargetForMode_(group, cmode);
}

template<typename Scalar, typename IndexTraits>
Scalar
GroupStateHelper<Scalar, IndexTraits>::
getProductionGroupTargetForMode(const Group& group,
                                 Group::ProductionCMode cmode) const
{
    return this->getProductionGroupTargetForMode_(group, cmode);
}

template<typename Scalar, typename IndexTraits>
GuideRateModel::Target
GroupStateHelper<Scalar, IndexTraits>::
getProductionGuideTargetMode(const Group& group) const
{
    const auto cmode = this->groupState().production_control(group.name());
    switch (cmode) {
    case Group::ProductionCMode::ORAT:
        return GuideRateModel::Target::OIL;
    case Group::ProductionCMode::WRAT:
        return GuideRateModel::Target::WAT;
    case Group::ProductionCMode::GRAT:
        return GuideRateModel::Target::GAS;
    case Group::ProductionCMode::LRAT:
        return GuideRateModel::Target::LIQ;
    case Group::ProductionCMode::RESV:
        return GuideRateModel::Target::RES;
    default:
        return GuideRateModel::Target::NONE;
    }
}

template<typename Scalar, typename IndexTraits>
GuideRateModel::Target
GroupStateHelper<Scalar, IndexTraits>::
getInjectionGuideTargetMode(Phase injection_phase) const
{
    switch (injection_phase) {
    case Phase::WATER:
        return GuideRateModel::Target::WAT;
    case Phase::OIL:
        return GuideRateModel::Target::OIL;
    case Phase::GAS:
        return GuideRateModel::Target::GAS;
    default:
        OPM_DEFLOG_THROW(std::logic_error,
                         "Invalid injection phase in getInjectionGuideTargetMode",
                         this->deferredLogger());
        return GuideRateModel::Target::NONE;
    }
}

template<typename Scalar, typename IndexTraits>
std::pair<Scalar, Group::ProductionCMode>
GroupStateHelper<Scalar, IndexTraits>::
getAutoChokeGroupProductionTargetRate(const Group& bottom_group,
                                      const Group& group,
                                      const std::vector<Scalar>& resv_coeff,
                                      Scalar efficiencyFactor) const
{
    const Group::ProductionCMode& currentGroupControl
        = this->groupState().production_control(group.name());
    if (currentGroupControl == Group::ProductionCMode::FLD ||
        currentGroupControl == Group::ProductionCMode::NONE) {
        if (!group.productionGroupControlAvailable()) {
            return std::make_pair(1.0, currentGroupControl);
        } else {
            // Produce share of parents control
            const auto& parent = this->schedule_.getGroup(group.parent(), this->report_step_);
            efficiencyFactor *= group.getGroupEfficiencyFactor();
            return this->getAutoChokeGroupProductionTargetRate(bottom_group, parent,
                                                               resv_coeff, efficiencyFactor);
        }
    }

    if (!group.isProductionGroup()) {
        return std::make_pair(1.0, currentGroupControl);
    }

    // If we are here, we are at the topmost group to be visited in the recursion.
    // This is the group containing the control we will check against.

    GroupStateHelpers::TargetCalculator<Scalar, IndexTraits> tcalc{*this,
                                                                    resv_coeff,
                                                                    group};

    GroupStateHelpers::FractionCalculator<Scalar, IndexTraits> fcalc(this->schedule_,
                                                                      *this,
                                                                      this->summary_state_,
                                                                      this->report_step_,
                                                                      &this->guide_rate_,
                                                                      this->getProductionGuideTargetMode(group),
                                                                      true,
                                                                      Phase::OIL);

    auto localFraction = [&](const std::string& child) {
        return fcalc.localFraction(child, child); //Note child needs to be passed to always include since the global isGrup map is not updated yet.
    };

    auto localReduction = [&](const std::string& group_name) {
        const std::vector<Scalar>& groupTargetReductions
            = this->groupState().production_reduction_rates(group_name);
        return tcalc.calcModeRateFromRates(groupTargetReductions);
    };

    const Scalar orig_target = this->getProductionGroupTarget(group);
    const auto chain = this->groupChainTopBot(bottom_group.name(), group.name());
    const std::size_t local_reduction_level = this->getLocalReductionLevel_(
        chain, /*is_production_group=*/true, Phase::OIL);

    // Delegate to applyReductionsAndFractions_ with no addback.  Addback is not
    // needed because autochoke wells are excluded from reductions via the
    // !group.as_choke() check in updateGroupTargetReductionRecursive_().
    const Scalar target = this->applyReductionsAndFractions_(chain,
                                                              orig_target,
                                                              /*current_rate_available=*/Scalar(0.0),
                                                              local_reduction_level,
                                                              /*is_production_group=*/true,
                                                              /*injection_phase=*/Phase::OIL,
                                                              localReduction,
                                                              localFraction,
                                                              /*do_addback=*/false);
    // Avoid negative target rates coming from too large local reductions.
    const Scalar target_rate = std::max(Scalar(0.0), target / efficiencyFactor);

    return std::make_pair(target_rate, currentGroupControl);
}

template <typename Scalar, typename IndexTraits>
std::vector<Scalar>
GroupStateHelper<Scalar, IndexTraits>::
getGroupRatesAvailableForHigherLevelControl(const Group& group, const bool is_injector) const
{
    std::vector<Scalar> rates(this->numPhases(), 0.0);
    if (is_injector) {
        const std::vector<Scalar> reduction_rates = this->groupState().injection_reduction_rates(group.name());
        for (int phasePos = 0; phasePos < this->numPhases(); ++phasePos) {
            const Scalar local_current_rate = this->sumWellPhaseRates(
                /*res_rates=*/false, group, phasePos, /*is_injector=*/true);
            rates[phasePos] = this->comm_.sum(local_current_rate) - reduction_rates[phasePos];
        }
    }
    else {
        const std::vector<Scalar> reduction_rates = this->groupState().production_reduction_rates(group.name());
        for (int phasePos = 0; phasePos < this->numPhases(); ++phasePos) {
            const Scalar local_current_rate = this->sumWellPhaseRates(
                /*res_rates=*/false, group, phasePos, /*is_injector=*/false);
            rates[phasePos] = -this->comm_.sum(local_current_rate) - reduction_rates[phasePos];
        }
    }
    return rates;
}

template <typename Scalar, typename IndexTraits>
std::optional<typename SingleWellState<Scalar, IndexTraits>::GroupTarget>
GroupStateHelper<Scalar, IndexTraits>::getWellGroupTargetInjector(const std::string& name,
                                                                 const std::string& parent,
                                                                 const Group& group,
                                                                 const Scalar* rates,
                                                                 Phase injection_phase,
                                                                 const Scalar efficiency_factor,
                                                                 const std::vector<Scalar>& resv_coeff) const
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
                                                resv_coeff);
    }

    // This can be false for FLD-controlled groups, we must therefore
    // check for FLD first (done above).
    if (!group.isInjectionGroup()) {
        return std::nullopt;
    }

    // If we are here, we are at the topmost group to be visited in the recursion.
    // This is the group containing the control we will check against.
    GroupStateHelpers::InjectionTargetCalculator<Scalar, IndexTraits> tcalc{*this,
                                                                             injection_phase};

    GroupStateHelpers::FractionCalculator<Scalar, IndexTraits> fcalc {this->schedule_,
                                                                      *this,
                                                                      this->summary_state_,
                                                                      this->report_step_,
                                                                      &this->guide_rate_,
                                                                      this->getInjectionGuideTargetMode(injection_phase),
                                                                      /*is_producer=*/false,
                                                                      injection_phase};

    auto local_reduction_lambda = [&](const std::string& group_name) {
        const std::vector<Scalar>& group_target_reductions
            = this->groupState().injection_reduction_rates(group_name);
        return tcalc.calcModeRateFromRates(group_target_reductions);
    };

    // For wells under individual control (not GRUP), include the well name as always_included_child.
    // For wells under GRUP control, use the child group itself as always_included_child.
    const bool is_grup_control = this->wellState().isInjectionGrup(name);
    auto local_fraction_lambda = [&](const std::string& child) {
        const std::string& always_included = is_grup_control ? child : name;
        return fcalc.localFraction(child, always_included);
    };

    const Scalar orig_target = this->getInjectionGroupTarget(group, injection_phase, resv_coeff);
    const Scalar current_rate_available = tcalc.calcModeRateFromRates(rates);
    const auto chain = this->groupChainTopBot(name, group.name());

    const std::size_t local_reduction_level = this->getLocalReductionLevel_(
        chain, /*is_production_group=*/false, injection_phase);

    // Add-back is only performed when the well is under individual control (not GRUP)
    const bool do_addback = !is_grup_control;

    const Scalar target = this->applyReductionsAndFractions_(chain,
                                                             orig_target,
                                                             current_rate_available,
                                                             local_reduction_level,
                                                             /*is_production_group=*/false,
                                                             injection_phase,
                                                             local_reduction_lambda,
                                                             local_fraction_lambda,
                                                             do_addback);

    // Avoid negative target rates coming from too large local reductions.
    return GroupTarget{group.name(), std::max(Scalar(0.0), target / efficiency_factor)};
}

template <typename Scalar, typename IndexTraits>
std::optional<typename SingleWellState<Scalar, IndexTraits>::GroupTarget>
GroupStateHelper<Scalar, IndexTraits>::getWellGroupTargetProducer(const std::string& name,
                                                                 const std::string& parent,
                                                                 const Group& group,
                                                                 const Scalar* rates,
                                                                 const Scalar efficiency_factor,
                                                                 const std::vector<Scalar>& resv_coeff) const
{
    // This function computes a wells group target.
    // 'parent' will be the name of 'group'. But if we recurse, 'name' and
    // 'parent' will stay fixed while 'group' will be higher up
    // in the group tree.
    // Efficiency factor is the well efficiency factor for the first group the well is
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
                                                resv_coeff);
    }

    // This can be false for FLD-controlled groups, we must therefore
    // check for FLD first (done above).
    if (!group.isProductionGroup()) {
        return std::nullopt;
    }

    // If we are here, we are at the topmost group to be visited in the recursion.
    // This is the group containing the control we will check against.

    GroupStateHelpers::TargetCalculator<Scalar, IndexTraits> tcalc{*this,
                                                                    resv_coeff,
                                                                    group};

    GroupStateHelpers::FractionCalculator<Scalar, IndexTraits> fcalc {this->schedule_,
                                                                      *this,
                                                                      this->summary_state_,
                                                                      this->report_step_,
                                                                      &this->guide_rate_,
                                                                      this->getProductionGuideTargetMode(group),
                                                                      /*is_producer=*/true,
                                                                      /*injection_phase=*/Phase::OIL};

    auto local_reduction_lambda = [&](const std::string& group_name) {
        const std::vector<Scalar>& group_target_reductions
            = this->groupState().production_reduction_rates(group_name);
        return tcalc.calcModeRateFromRates(group_target_reductions);
    };

    // For wells under individual control (not GRUP), include the well name as always_included_child.
    // For wells under GRUP control, use the child group itself as always_included_child.
    const bool is_grup_control = this->wellState().isProductionGrup(name);
    auto local_fraction_lambda = [&](const std::string& child) {
        const std::string& always_included = is_grup_control ? child : name;
        return fcalc.localFraction(child, always_included);
    };

    const Scalar orig_target = this->getProductionGroupTarget(group);
    // Switch sign since 'rates' are negative for producers.
    const Scalar current_rate_available = -tcalc.calcModeRateFromRates(rates);
    const auto chain = this->groupChainTopBot(name, group.name());

    const std::size_t local_reduction_level = this->getLocalReductionLevel_(
        chain, /*is_production_group=*/true, /*injection_phase=*/Phase::OIL);

    // Add-back is only performed when the well is under individual control (not GRUP)
    const bool do_addback = !is_grup_control;

    const Scalar target = this->applyReductionsAndFractions_(chain,
                                                             orig_target,
                                                             current_rate_available,
                                                             local_reduction_level,
                                                             /*is_production_group=*/true,
                                                             /*injection_phase=*/Phase::OIL,
                                                             local_reduction_lambda,
                                                             local_fraction_lambda,
                                                             do_addback);

    // Avoid negative target rates coming from too large local reductions.
    return GroupTarget{group.name(), std::max(Scalar(0.0), target / efficiency_factor)};
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
int
GroupStateHelper<Scalar, IndexTraits>::phaseToActivePhaseIdx(const Phase phase) const
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
                                                        const bool is_injector,
                                                        const bool network) const
{
    if (this->isRank0()) {
        // Only obtain satellite rates once (on rank 0)
        if (this->isSatelliteGroup_(group)) {
            return this->getSatelliteRate_(group, phase_pos, res_rates, is_injector);
        }
#ifdef RESERVOIR_COUPLING_ENABLED
        if (this->isReservoirCouplingMasterGroup(group)) {
            using RateKind = ReservoirCoupling::RateKind;
            RateKind kind;
            if (is_injector) {
                kind = res_rates ? RateKind::InjectionReservoir : RateKind::InjectionSurface;
            } else {
                if (res_rates)
                    kind = RateKind::ProductionReservoir;
                else if (network)
                    kind = RateKind::ProductionNetworkSurface;
                else
                    kind = RateKind::ProductionSurface;
            }
            return this->getReservoirCouplingMasterGroupRate_(group, phase_pos, kind);
        }
#endif
    }
    Scalar rate = 0.0;
    for (const std::string& group_name : group.groups()) {
        const auto& group_tmp = this->schedule_.getGroup(group_name, this->report_step_);
        const auto& gefac = group_tmp.getGroupEfficiencyFactor(network);
        rate += gefac * this->sumWellPhaseRates(res_rates, group_tmp, phase_pos, is_injector, network);
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
        if (ws.status == Opm::Well::Status::SHUT)
            continue;

        const Scalar factor = well_ecl.getEfficiencyFactor(network)
            * this->wellState().well(well_index.value()).efficiency_scaling_factor;
        if (res_rates) {
            const auto& well_rates = ws.reservoir_rates;
            if (is_injector)
                rate += factor * well_rates[phase_pos];
            else
                rate -= factor * well_rates[phase_pos];
        } else {
            const auto& well_rates = ws.surface_rates;
            if (is_injector)
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
//   - The original implementation in GroupStateHelpers.cpp (static utility class) could keep
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
                                                                  const Phase injection_phase)
{
    OPM_TIMEFUNCTION();
    const auto& group_name = "FIELD";
    return this->updateGroupControlledWellsRecursive_(group_name, is_production_group, injection_phase);
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
    const int np = this->numPhases();
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
    std::vector<Scalar> group_target_reduction(this->numPhases(), 0.0);
    this->updateGroupTargetReductionRecursive_(group, is_injector, group_target_reduction);
}

template <typename Scalar, typename IndexTraits>
void
GroupStateHelper<Scalar, IndexTraits>::updateNetworkLeafNodeProductionRates()
{
    const auto& network = this->schedule_[this->report_step_].network();
    if (network.active()) {
        const int np = this->numPhases();
        for (const auto& group_name : network.leaf_nodes()) {
            std::vector<Scalar> network_rates(np, 0.0);
            if (this->schedule_[this->report_step_].groups.has(
                    group_name)) { // Allow empty leaf nodes that are not groups
                const auto& group = this->schedule_[this->report_step_].groups.get(group_name);
                for (int phase = 0; phase < np; ++phase) {
                    network_rates[phase] = this->sumWellPhaseRates(
                        /*res_rates=*/false, group, phase, /*injector=*/false, /*network=*/true);
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
    const int np = this->numPhases();
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
        const auto& pu = this->phaseUsage();
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
    const int np = this->numPhases();
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
    const int np = this->numPhases();
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
    const int np = this->numPhases();
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
    const int np = this->numPhases();
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

    const int np = this->numPhases();
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
                                                          const Group::ProductionCMode& offended_control) const
{
    std::pair<std::optional<std::string>, Scalar> offending_well {std::nullopt, 0.0};
    for (const std::string& child_group : group.groups()) {
        const auto& this_group = this->schedule_.getGroup(child_group, this->report_step_);
        const auto& offending_well_this = this->worstOffendingWell(this_group, offended_control);
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
                for (int p = 0; p < this->numPhases(); ++p) {
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
                                 this->deferredLogger());
                break;
            case Group::ProductionCMode::CRAT:
                OPM_DEFLOG_THROW(std::runtime_error,
                                 "Group " + group.name() + " GroupProductionCMode CRAT not implemented",
                                 this->deferredLogger());
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
        violating_rate = this->comm_.sum(violating_rate);
        if (violating_rate < 0) { // only check producing wells
            prefered_rate = this->comm_.sum(prefered_rate);
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

#ifdef RESERVOIR_COUPLING_ENABLED
template <typename Scalar, typename IndexTraits>
ReservoirCoupling::Phase
GroupStateHelper<Scalar, IndexTraits>::
activePhaseIdxToRescoupPhase_(int phase_pos) const
{
    const auto& pu = this->phase_usage_info_;
    // Map active phase index back to canonical phase, then to ReservoirCoupling::Phase
    if (pu.phaseIsActive(IndexTraits::oilPhaseIdx) &&
        phase_pos == pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx)) {
        return ReservoirCoupling::Phase::Oil;
    }
    if (pu.phaseIsActive(IndexTraits::gasPhaseIdx) &&
        phase_pos == pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx)) {
        return ReservoirCoupling::Phase::Gas;
    }
    if (pu.phaseIsActive(IndexTraits::waterPhaseIdx) &&
        phase_pos == pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx)) {
        return ReservoirCoupling::Phase::Water;
    }
    // TODO: We would like to use OPM_DEFLOG_THROW() here, but it is requires the deferred logger to be passed
    //   as an argument to all the sumWellPhaseRates() calls, which is a major refactoring effort.
    //   Alternatively, we could store a reference to the deferred logger in the GroupStateHelper class,
    //   which would also require a major refactoring effort.
    throw std::logic_error("Invalid phase_pos in activePhaseIdxToRescoupPhase");
    return ReservoirCoupling::Phase::Oil; // just to avoid warning
}
#endif

template <typename Scalar, typename IndexTraits>
template <typename ReductionLambda, typename FractionLambda>
Scalar
GroupStateHelper<Scalar, IndexTraits>::applyReductionsAndFractions_(const std::vector<std::string>& chain,
                                                                    const Scalar orig_target,
                                                                    const Scalar current_rate_available,
                                                                    const std::size_t local_reduction_level,
                                                                    const bool is_production_group,
                                                                    const Phase injection_phase,
                                                                    ReductionLambda&& local_reduction_lambda,
                                                                    FractionLambda&& local_fraction_lambda,
                                                                    const bool do_addback) const
{
    // Compute portion of target corresponding to current_rate_available.
    // The chain is ordered [control_group (top), ..., bottom_group].
    // num_ancestors excludes the bottom group itself. For the bottom group, the always_included_child argument
    //   (encapsulated into the local_fraction_lambda) to FractionCalculator::localFraction() is used to calculate
    // the correct guide rate sum for a group under individual control.
    const std::size_t num_ancestors = chain.size() - 1;
    Scalar target = orig_target;
    for (std::size_t ii = 0; ii < num_ancestors; ++ii) {
        // Check if this level has a guide rate (or is the top level)
        const bool has_guide_rate = is_production_group
            ? this->guide_rate_.has(chain[ii])
            : this->guide_rate_.has(chain[ii], injection_phase);
        if ((ii == 0) || has_guide_rate) {
            // Apply local reductions only at the control level (top)
            // and for levels where we have a specified group guide rate.
            if (local_reduction_level >= ii) {
                target -= local_reduction_lambda(chain[ii]);
            }
            // Add the bottom group's rate back at the level where it is included
            // in the local reduction (if add-back is requested).
            if (do_addback && local_reduction_level == ii) {
                const Scalar addback_efficiency
                    = this->computeAddbackEfficiency_(chain, local_reduction_level);
                target += current_rate_available * addback_efficiency;
            }
        }
        // Apply guide rate fraction to get portion for next level
        target *= local_fraction_lambda(chain[ii + 1]);
    }
    return target;
}

template <typename Scalar, typename IndexTraits>
Scalar
GroupStateHelper<Scalar, IndexTraits>::computeAddbackEfficiency_(const std::vector<std::string>& chain,
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
        if (this->schedule_.hasGroup(name, this->report_step_)) {
            const auto& grp = this->schedule_.getGroup(name, this->report_step_);
            efficiency *= grp.getGroupEfficiencyFactor();
        } else if (this->schedule_.hasWell(name, this->report_step_)) {
            const auto& well = this->schedule_.getWell(name, this->report_step_);
            efficiency *= well.getEfficiencyFactor()
                * this->wellState()[name].efficiency_scaling_factor;
        }
    }
    return efficiency;
}

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
std::size_t
GroupStateHelper<Scalar, IndexTraits>::getLocalReductionLevel_(const std::vector<std::string>& chain,
                                                               const bool is_production_group,
                                                               const Phase injection_phase) const
{
    // The local reduction level is the deepest level in the chain (starting from level 1)
    // where a group has both a guide rate and group-controlled wells (GCW > 0).
    // Level 0 is the top controlling group, so we start from level 1.
    const std::size_t num_ancestors = chain.size() - 1;
    std::size_t local_reduction_level = 0;
    for (std::size_t ii = 1; ii < num_ancestors; ++ii) {
        const int num_gr_ctrl = this->groupControlledWells(chain[ii],
                                                           /*always_included_child=*/"",
                                                           is_production_group,
                                                           injection_phase);
        const bool has_guide_rate = is_production_group
            ? this->guide_rate_.has(chain[ii])
            : this->guide_rate_.has(chain[ii], injection_phase);
        if (has_guide_rate && num_gr_ctrl > 0) {
            local_reduction_level = ii;
        }
    }
    return local_reduction_level;
}

template<typename Scalar, typename IndexTraits>
Scalar
GroupStateHelper<Scalar, IndexTraits>::
getProductionConstraintTarget_(const Group& group,
                                Group::ProductionCMode cmode,
                                const Group::ProductionControls& controls) const
{
    Scalar target = 0.0;
    switch (cmode) {
    case Group::ProductionCMode::ORAT: target = controls.oil_target; break;
    case Group::ProductionCMode::WRAT: target = controls.water_target; break;
    case Group::ProductionCMode::GRAT: target = controls.gas_target; break;
    case Group::ProductionCMode::LRAT: target = controls.liquid_target; break;
    case Group::ProductionCMode::RESV: {
        if (group.has_gpmaint_control(Group::ProductionCMode::RESV))
            target = this->groupState().gpmaint_target(group.name());
        else
            target = controls.resv_target;
        break;
    }
    default: return 0.0;
    }
#ifdef RESERVOIR_COUPLING_ENABLED
    if (this->isReservoirCouplingSlave()
        && this->isReservoirCouplingSlaveGroup(group))
    {
        target = this->getEffectiveProductionLimit(group.name(), cmode, target);
    }
#endif
    return target;
}

template<typename Scalar, typename IndexTraits>
Scalar
GroupStateHelper<Scalar, IndexTraits>::
getProductionGroupTargetForMode_(const Group& group,
                                  Group::ProductionCMode cmode) const
{
    Group::ProductionControls ctrl = group.productionControls(this->summary_state_);
    switch (cmode) {
    case Group::ProductionCMode::ORAT:
        return ctrl.oil_target;
    case Group::ProductionCMode::WRAT:
        return ctrl.water_target;
    case Group::ProductionCMode::GRAT:
    {
        Scalar grat_target_from_sales = 0.0;
        if (this->groupState().has_grat_sales_target(group.name())) {
            grat_target_from_sales = this->groupState().grat_sales_target(group.name());
        }
        // gas target may have been adjusted by GCONSALE
        if (grat_target_from_sales > 0)
            return grat_target_from_sales;

        return ctrl.gas_target;
    }
    case Group::ProductionCMode::LRAT:
        return ctrl.liquid_target;
    case Group::ProductionCMode::RESV:
    {
        bool use_gpmaint = group.has_gpmaint_control(cmode);
        if (use_gpmaint && this->groupState().has_gpmaint_target(group.name()))
            return this->groupState().gpmaint_target(group.name());

        return ctrl.resv_target;
    }
    default:
        OPM_DEFLOG_THROW(std::logic_error,
                         "Invalid Group::ProductionCMode in getProductionGroupTargetForMode_",
                         this->deferredLogger());
        return 0.0;
    }
}

template<typename Scalar, typename IndexTraits>
std::pair<Group::ProductionCMode, Scalar>
GroupStateHelper<Scalar, IndexTraits>::
checkProductionRateConstraint_(const Group& group,
                                Group::ProductionCMode cmode,
                                Group::ProductionCMode currentControl,
                                Scalar target,
                                Scalar current_rate) const
{
    if (!group.has_control(cmode)) {
        return {Group::ProductionCMode::NONE, Scalar(1.0)};
    }
    if (currentControl == cmode) {
        return {Group::ProductionCMode::NONE, Scalar(1.0)};
    }
    if (target < current_rate) {
        Scalar scale = 1.0;
        if (current_rate > 1e-12)
            scale = target / current_rate;
        return {cmode, scale};
    }
    return {Group::ProductionCMode::NONE, Scalar(1.0)};
}

#ifdef RESERVOIR_COUPLING_ENABLED
template <typename Scalar, typename IndexTraits>
Scalar
GroupStateHelper<Scalar, IndexTraits>::
getEffectiveProductionLimit(
    const std::string& gname,
    Group::ProductionCMode rate_type,
    Scalar slave_local_target) const
{
    auto& slave = this->reservoirCouplingSlave();
    if (!slave.hasMasterProductionLimits(gname)) {
        return slave_local_target;
    }
    const auto& limits = slave.masterProductionLimits(gname);
    Scalar master_limit = -1;
    switch (rate_type) {
    case Group::ProductionCMode::ORAT: master_limit = limits.oil_limit; break;
    case Group::ProductionCMode::WRAT: master_limit = limits.water_limit; break;
    case Group::ProductionCMode::GRAT: master_limit = limits.gas_limit; break;
    case Group::ProductionCMode::LRAT: master_limit = limits.liquid_limit; break;
    case Group::ProductionCMode::RESV: master_limit = limits.resv_limit; break;
    default: return slave_local_target;
    }
    if (master_limit < 0) {
        return slave_local_target;  // no master limit for this rate type
    }
    auto filter = this->getProductionFilterFlag_(gname, rate_type);
    using FilterFlag = ReservoirCoupling::GrupSlav::FilterFlag;
    if (filter == FilterFlag::MAST) {
        return master_limit;
    }
    if (filter == FilterFlag::BOTH) {
        return std::min(master_limit, slave_local_target);
    }
    // FilterFlag::SLAV
    return slave_local_target;
}

template <typename Scalar, typename IndexTraits>
ReservoirCoupling::GrupSlav::FilterFlag
GroupStateHelper<Scalar, IndexTraits>::
getInjectionFilterFlag_(const std::string& group_name,
                        Phase injection_phase) const
{
    const auto& rescoup_info = this->schedule_[this->report_step_].rescoup();
    if (!rescoup_info.hasGrupSlav(group_name)) {
        return ReservoirCoupling::GrupSlav::FilterFlag::MAST;
    }
    const auto& grup_slav = rescoup_info.grupSlav(group_name);
    switch (injection_phase) {
    case Phase::OIL:
        return grup_slav.oilInjFlag();
    case Phase::WATER:
        return grup_slav.waterInjFlag();
    case Phase::GAS:
        return grup_slav.gasInjFlag();
    default:
        return ReservoirCoupling::GrupSlav::FilterFlag::MAST;
    }
}

template <typename Scalar, typename IndexTraits>
ReservoirCoupling::GrupSlav::FilterFlag
GroupStateHelper<Scalar, IndexTraits>::
getProductionFilterFlag_(const std::string& group_name,
                         Group::ProductionCMode cmode) const
{
    const auto& rescoup_info = this->schedule_[this->report_step_].rescoup();
    if (!rescoup_info.hasGrupSlav(group_name)) {
        return ReservoirCoupling::GrupSlav::FilterFlag::MAST;
    }
    const auto& grup_slav = rescoup_info.grupSlav(group_name);
    switch (cmode) {
    case Group::ProductionCMode::ORAT:
        return grup_slav.oilProdFlag();
    case Group::ProductionCMode::WRAT:
    case Group::ProductionCMode::LRAT:
        return grup_slav.liquidProdFlag();
    case Group::ProductionCMode::GRAT:
        return grup_slav.gasProdFlag();
    case Group::ProductionCMode::RESV:
        return grup_slav.fluidVolumeProdFlag();
    default:
        return ReservoirCoupling::GrupSlav::FilterFlag::MAST;
    }
}

template <typename Scalar, typename IndexTraits>
Scalar
GroupStateHelper<Scalar, IndexTraits>::
getReservoirCouplingMasterGroupRate_(const Group& group,
                                     const int phase_pos,
                                     ReservoirCoupling::RateKind kind) const
{
    if (this->isReservoirCouplingMaster()) {
        ReservoirCoupling::Phase rescoup_phase = this->activePhaseIdxToRescoupPhase_(phase_pos);
        return this->reservoirCouplingMaster().getMasterGroupRate(group.name(), rescoup_phase, kind);
    }
    else {
        return 0.0;
    }
}
#endif // RESERVOIR_COUPLING_ENABLED

template <typename Scalar, typename IndexTraits>
Scalar
GroupStateHelper<Scalar, IndexTraits>::getSatelliteRate_(const Group& group,
                                                               const int phase_pos,
                                                               const bool res_rates,
                                                               const bool is_injector) const
{
    // Only obtain satellite rates once (on rank 0)
    assert(this->isRank0() && this->isSatelliteGroup_(group));
    Scalar rate = 0.0;
    if (is_injector) {
        return this->satelliteInjectionRate_(this->schedule_[this->report_step_], group, phase_pos, res_rates);
    } else {
        const auto rate_comp = this->selectRateComponent_(phase_pos);
        if (rate_comp.has_value()) {
            return this->satelliteProductionRate_(this->schedule_[this->report_step_], group, *rate_comp, res_rates);
        }
    }
    return rate;
}

template <typename Scalar, typename IndexTraits>
bool
GroupStateHelper<Scalar, IndexTraits>::isAutoChokeGroupUnderperforming_(const Group& group) const
{
    // Only applies to production auto choke groups
    if (!group.as_choke()) {
        return false;
    }

    // Auto choke groups inherit control from an ancestor (cmode == FLD or NONE)
    const auto ctrl = this->groupState().production_control(group.name());
    if (ctrl != Group::ProductionCMode::FLD && ctrl != Group::ProductionCMode::NONE) {
        return false;
    }

    // Check if guide rate is set for this group
    const auto& group_guide_rate = group.productionControls(this->summary_state_).guide_rate;
    if (group_guide_rate <= 0) {
        return false;
    }

    // Find the ancestor group that has actual control
    const auto& control_group_name = this->controlGroup_(group);
    const auto& control_group = this->schedule_.getGroup(control_group_name, this->report_step_);

    // Calculate target rate for this auto choke group using a guide rate ratio
    // formula. This is an approximation that ignores intermediate reductions.
    // A more accurate calculation would use getAutoChokeGroupProductionTargetRate(),
    // but that requires reduction rates not yet available during
    // updateGroupControlledWells() (circular dependency). The efficiency factor is however accounted for
    // below meaning that the current target rate will be correct for the special case where the control
    // group has no reductions and the auto choke group is a direct child of the control group.
    // TODO: The issue with missing handling of reductions should be fixed by refactoring the auto choke
    // target calculation below to use the more accurate
    // getAutoChokeGroupProductionTargetRate() function.
    const auto num_phases = this->numPhases();
    std::vector<Scalar> resv_coeff(num_phases, 1.0);
    GroupStateHelpers::TargetCalculator<Scalar, IndexTraits> tcalc{*this,
                                                                   resv_coeff,
                                                                   control_group};
    const Scalar control_group_target = this->getProductionGroupTarget(control_group);

    // Sum guide rates of the control group's FLD children. The control group's
    // own guide rate (from GCONPROD) is NOT used here — it governs the control
    // group's share of its parent's target, not how its own target is distributed.
    Scalar control_group_guide_rate = 0.0;
    for (const std::string& child : control_group.groups()) {
        const auto child_ctrl = this->groupState().production_control(child);
        if (child_ctrl == Group::ProductionCMode::FLD
            || child_ctrl == Group::ProductionCMode::NONE) {
            if (this->guide_rate_.has(child)) {
                control_group_guide_rate += this->guide_rate_.get(
                    child, this->getProductionGuideTargetMode(control_group),
                    this->getProductionGroupRateVector(child));
            }
        }
    }

    if (control_group_guide_rate <= 0) {
        return false;
    }

    // Accumulate group efficiency factors from the autochoke group up to (but not
    // including) the control group. This matches what
    // getAutoChokeGroupProductionTargetRate() does during its recursion.
    Scalar accumulated_efficiency = 1.0;
    {
        std::string current = group.name();
        while (current != control_group_name) {
            const auto& g = this->schedule_.getGroup(current, this->report_step_);
            accumulated_efficiency *= g.getGroupEfficiencyFactor();
            current = g.parent();
        }
    }

    const Scalar target_rate = control_group_target * group_guide_rate
                             / control_group_guide_rate / accumulated_efficiency;

    // Calculate current rate for this group
    std::vector<Scalar> rates(num_phases, 0.0);
    for (int phase_pos = 0; phase_pos < num_phases; ++phase_pos) {
        rates[phase_pos] = this->sumWellSurfaceRates(group, phase_pos, /*injector=*/false);
    }
    this->comm_.sum(rates.data(), rates.size());
    const Scalar current_rate = tcalc.calcModeRateFromRates(rates);

    // If underperforming, wells should be excluded from GCW count
    return current_rate < target_rate;
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
bool
GroupStateHelper<Scalar, IndexTraits>::isSatelliteGroup_(const Group& group) const
{
    return group.hasSatelliteProduction() || group.hasSatelliteInjection();
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

template<typename Scalar, typename IndexTraits>
Scalar
GroupStateHelper<Scalar, IndexTraits>::
sumProductionRate_(const Group& group,
                    Group::ProductionCMode cmode) const
{
    const auto& pu = this->phaseUsage();
    Scalar rate = 0.0;
    switch (cmode) {
    case Group::ProductionCMode::ORAT:
        rate = this->sumWellSurfaceRates(group,
            pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx), false);
        break;
    case Group::ProductionCMode::WRAT:
        rate = this->sumWellSurfaceRates(group,
            pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx), false);
        break;
    case Group::ProductionCMode::GRAT:
        rate = this->sumWellSurfaceRates(group,
            pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx), false);
        break;
    case Group::ProductionCMode::LRAT:
        rate = this->sumWellSurfaceRates(group,
            pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx), false);
        rate += this->sumWellSurfaceRates(group,
            pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx), false);
        break;
    case Group::ProductionCMode::RESV:
        rate = this->sumWellResRates(group,
            pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx), false);
        rate += this->sumWellResRates(group,
            pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx), false);
        rate += this->sumWellResRates(group,
            pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx), false);
        break;
    default:
        break;
    }
    return this->comm().sum(rate);
}

template <typename Scalar, typename IndexTraits>
int
GroupStateHelper<Scalar, IndexTraits>::updateGroupControlledWellsRecursive_(
    const std::string& group_name,
    const bool is_production_group,
    const Phase injection_phase)
{
    // NOTE: The number of group controlled wells (GCW) is only relevant for groups that are NOT under
    //   individual control. However, it is collected for all groups here.
    const Group& group = this->schedule_.getGroup(group_name, this->report_step_);
    int num_wells = 0;
    // NOTE: If "group" is a reservoir coupling master group, it should have no child groups or wells.
    //   As a convention, we will assign GCW = 1 (GCW=Group Controlled Wells). This is consistent with
    //   the usage of GCW elsewhere in the code: GCW is only checked for GCW>0, the number of wells if
    //   greater than zero is irrelevant.
    //   1) for target reductions, GCW is used to check if a group not under individual control and with a guide rate
    //      set should participate in target reduction, and
    //   2) for guide rate control, GCW is used to check if a group not under individual control should (or
    //      under individual control if "always_included_child" in localFraction():FractionCalculator.cpp is true)
    //      should participate in guide rate control.
    if (this->isReservoirCouplingMasterGroup(group)) {
        num_wells = 1;
    }
    else {
        // NOTE: A group with sub groups cannot also have direct wells (one level below) under its control.
        // So a group (either under individual control or not) having all its direct child groups (one level below)
        // under individual control will have GCW = 0 (GCW=Group Controlled Wells).
        // I.e. sub groups that are under individual control do not contribute to the group's GCW. (Another way
        // to phrase it: It will only have GCW=0 if all child groups not under individual control have GCW=0.)
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
                    child_group, is_production_group, injection_phase);
            } else {
                this->updateGroupControlledWellsRecursive_(
                    child_group, is_production_group, injection_phase);
            }
        }
        // For production auto choke groups: check if we should exclude all wells from GCW count.
        // If the group is underperforming its target, wells are not counted as group-controlled,
        // effectively excluding this group from guide rate distribution at the parent level.
        const bool exclude_for_auto_choke = is_production_group
            && this->isAutoChokeGroupUnderperforming_(group);

        // Below loop is only entered for well groups (i.e. groups with only wells as direct children).
        // For such a group (assuming it is not an auto choke group), if all its wells are under individual control,
        //   then the group will have GCW = 0 (GCW=Group Controlled Wells). This is true whether the group itself is
        //   under individual control or not.
        for (const std::string& child_well : group.wells()) {
            bool included = false;
            const Well& well = this->schedule_.getWell(child_well, this->report_step_);
            if (is_production_group && well.isProducer()) {
                included = (this->wellState().isProductionGrup(child_well) || group.as_choke());
                // Auto choke groups: exclude all wells if group is underperforming its target
                if (exclude_for_auto_choke) {
                    included = false;
                }
            } else if (!is_production_group && !well.isProducer()) {
                const auto& well_controls = well.injectionControls(this->summary_state_);
                auto injectorType = well_controls.injector_type;
                if ((injection_phase == Phase::WATER && injectorType == InjectorType::WATER)
                    || (injection_phase == Phase::OIL && injectorType == InjectorType::OIL)
                    || (injection_phase == Phase::GAS && injectorType == InjectorType::GAS)) {
                    included = this->wellState().isInjectionGrup(child_well);
                }
            }
            if (included) {
                ++num_wells;
            }
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
    const int np = this->numPhases();
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
                const auto phase_pos = this->phaseToActivePhaseIdx(phase);
                // the phase is not present
                if (phase_pos == -1)
                    continue;
                const Group::InjectionCMode& current_group_control
                    = this->groupState().injection_control(sub_group.name(), phase);
                const bool individual_control = (current_group_control != Group::InjectionCMode::FLD
                                                 && current_group_control != Group::InjectionCMode::NONE);
            // NOTE: A reservoir coupling master group: will have GCW (group controlled wells) set to 1 by convention.
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
                    // NOTE: For reservoir coupling master groups that are not under individual control
                    //   it is required that they have a guide rate set, see GroupTargetCalculator.cpp.
                    if (!this->guide_rate_.has(sub_group.name(), phase)) {
                        group_target_reduction[phase_pos]
                            += sub_group_efficiency * sub_group_target_reduction[phase_pos];
                    }
                    // else: if [(NOT indivdual control) AND (guide rate is set) AND (GCW > 0)]
                    //   it means that the current target will be distributed according to guide rate fraction,
                    //   and then the current (higher level) group target should not be reduced at this level.
                }
            }
        } else {
            const Group::ProductionCMode& current_group_control
                = this->groupState().production_control(sub_group.name());
            const bool individual_control = (current_group_control != Group::ProductionCMode::FLD
                                             && current_group_control != Group::ProductionCMode::NONE);
            // NOTE: A reservoir coupling master group: will have GCW (group controlled wells) set to 1 by convention.
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
                // NOTE: Reservoir coupling master groups that are not under individual control
                //   must have a guide rate set, see GroupTargetCalculator.cpp.
                if (!this->guide_rate_.has(sub_group.name())) {
                    // Accumulate from this subgroup only if no group guide rate is set for it.
                    for (int phase = 0; phase < np; phase++) {
                        group_target_reduction[phase]
                            += sub_group_efficiency * sub_group_target_reduction[phase];
                    }
                }
                // else: if [(NOT indivdual control) AND (guide rate is set) AND (GCW > 0)]
                //   it means that the current target will be distributed according to guide rate fraction,
                //   and then the current (higher level) group target should not be reduced at this level.
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

#ifdef RESERVOIR_COUPLING_ENABLED
template <typename Scalar, typename IndexTraits>
void
GroupStateHelper<Scalar, IndexTraits>::
updateSlaveGroupCmodesFromMaster()
{
    using FilterFlag = ReservoirCoupling::GrupSlav::FilterFlag;
    auto& slave = this->reservoirCouplingSlave();
    for (std::size_t i = 0; i < slave.numSlaveGroups(); ++i) {
        const auto& gname = slave.slaveGroupIdxToGroupName(i);
        // Production cmode
        if (slave.hasMasterProductionTarget(gname)) {
            auto [master_target, master_cmode] = slave.masterProductionTarget(gname);
            auto filter = this->getProductionFilterFlag_(gname, master_cmode);
            if (filter != FilterFlag::SLAV) {
                this->groupState().production_control(gname, master_cmode);
            }
        }
        // Injection cmode — check each phase
        for (const auto phase : {Phase::WATER, Phase::OIL, Phase::GAS}) {
            if (slave.hasMasterInjectionTarget(gname, phase)) {
                auto [master_target, master_cmode] = slave.masterInjectionTarget(gname, phase);
                auto filter = this->getInjectionFilterFlag_(gname, phase);
                if (filter != FilterFlag::SLAV) {
                    this->groupState().injection_control(gname, phase, master_cmode);
                }
            }
        }
    }
}
#endif // RESERVOIR_COUPLING_ENABLED

template class GroupStateHelper<double, BlackOilDefaultFluidSystemIndices>;
#ifdef FLOW_INSTANTIATE_FLOAT
template class GroupStateHelper<float, BlackOilDefaultFluidSystemIndices>;
#endif
} // namespace Opm
