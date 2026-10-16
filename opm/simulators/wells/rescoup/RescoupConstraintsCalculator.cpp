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
#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>
#include <opm/simulators/wells/rescoup/RescoupConstraintsCalculator.hpp>

#include <array>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <fmt/format.h>

namespace Opm {

// -------------------------------------------------------
// Constructor for the RescoupConstraintsCalculator class
// -------------------------------------------------------
template <class Scalar, class IndexTraits>
RescoupConstraintsCalculator<Scalar, IndexTraits>::
RescoupConstraintsCalculator(
    GuideRateHandler<Scalar, IndexTraits>& guide_rate_handler,
    GroupStateHelper<Scalar, IndexTraits>& group_state_helper
)
    : guide_rate_handler_{guide_rate_handler}
    , group_state_helper_{group_state_helper}
    , well_state_{group_state_helper.wellState()}
    , group_state_{group_state_helper.groupState()}
    , report_step_idx_{group_state_helper.reportStepIdx()}
    , schedule_{group_state_helper.schedule()}
    , summary_state_{group_state_helper.summaryState()}
    , deferred_logger_{guide_rate_handler.deferredLogger()}
    , reservoir_coupling_master_{group_state_helper.reservoirCouplingMaster()}
    , well_model_{guide_rate_handler.wellModel()}
    , phase_usage_{group_state_helper.phaseUsage()}
{
}

// Calculates the constraints (target and per-rate-type limits) for each master group.
//
// Runs at the start of each sync step from the slave rates reported before the
// slaves solve their wells, and again from the rates reported after the
// slaves' initial well solve and at each network iteration; see
// BlackoilWellModelRescoup::refreshAndSendGroupConstraints().
//
// A master group defines both an active control-mode target and limits for other rate types
// (ORAT, WRAT, GRAT, LRAT, RESV). The active target alone is not sufficient for reservoir-coupling slaves:
// a slave group must know every effective limit so that it can enforce the most restrictive
// constraint for each rate type independently (the active cmode on the master side may differ
// from what is binding on the slave side).
//
// Overview of the phases:
// -----------------------
// The constraint calculation is split into three phases.  All phases run on
// every rank of the master communicator (collective MPI ops are required by
// the underlying GroupStateHelper / GroupConstraintCalculator); the
// slave-facing MPI sends in Phase 3 are rank-0-only inside the send helpers.
//
//  A note on "effective GCW" for master groups:
//    GCW (group controlled wells) gates whether a group participates in its
//    parent's guide-rate distribution (GCW>0) or instead has its own rate
//    subtracted as a parent target reduction (GCW=0).  For ordinary groups GCW
//    is derived from the control mode.  A master group that is available for higher-level
//    control takes part in its parent's distribution whatever its cmode (it can be on its own
//    ORAT limit, from the deck or by group control), so its cmode cannot be used to decide
//    participation. The master therefore maintains a per-master-group "effective GCW"
//    (ReservoirCouplingMaster::effectiveGCW), set in this routine and read by
//    GroupStateHelper::updateGroupControlledWellsRecursive_.  It separates two
//    concerns:
//      - static eligibility: GCONPROD item 8 RESPOND_TO_PARENT = YES
//        (Group::productionGroupControlAvailable) — a deck property; from
//      - dynamic weight this sync step: 1 if the group participates, 0 if
//        none of its slave producers are under group control (Phase 2) or if
//        its slave is inactive (pre-phase).
//    An eligible group with effective GCW = 0 is excluded from the guide-rate
//    sum and contributes its rate as a parent target reduction instead — which
//    is precisely how an excluded group's shortfall is redistributed to its
//    siblings.  The control mode of such a group is not touched by this routine:
//    it is the Schedule's, or the one group control has switched it to, and the
//    calculator and the master's own group control read it.
//
//  Pre-phase: reset the effective GCW (participating master groups default to
//    1), set the effective GCW of inactive slaves' master groups to 0 (and put
//    those not available for higher-level control on individual control) so
//    they don't consume any share of the parent's target, then recompute GCW +
//    FIELD target reductions with the resulting state.
//
//  Phase 1 (per-slave per-master-group): compute initial guide-rate-
//    distributed targets and per-rate-type limits via
//    GroupConstraintCalculator.  This walks up the group hierarchy from each
//    master group, applying local reductions and guide-rate fractions, and
//    returns the more restrictive of the higher-level distributed target and
//    the master group's own GCONPROD constraint.  For limits, the same
//    traversal is reused once per non-active rate type (ORAT/WRAT/GRAT/LRAT/
//    RESV) with an explicit cmode parameter; see the Phase 1 details below.
//
//  Phase 2 (exclude-and-redistribute): the master does for a master group what
//    Flow's group control does for a well or group that cannot meet its share:
//    it runs at its own limit, drops out of the guide-rate distribution, and its
//    current rate becomes a reduction on the parent's target.  A master group
//    none of whose slave producers are under group control (e.g. all on their
//    own THP or BHP limit) is handled in three steps:
//      1. Set its effective GCW to 0, so that it drops out of the parent's
//         guide-rate sum.  (Its control mode is left as it is: it is the
//         effective GCW, not the cmode, that drives the exclusion; see the
//         effective-GCW note above.)
//      2. Recompute the GCW and the target reductions, so that its current rate,
//         the one the slave last reported, is subtracted from the parent's target
//         as a reduction.  The remaining groups absorb the shortfall through the
//         standard localReduction / FractionCalculator machinery.
//      3. Recompute the targets of all groups.  The excluded group's rate is in
//         the parent's reduction and is added back for its own target, and its
//         own guide rate is counted in its fraction (see
//         FractionCalculator::groupControlledWells()), so it is sent the target
//         it would get on returning to group control.  That
//         leaves its wells at their own limits; if they can meet that target,
//         they return to group control and the group takes part in the
//         distribution again from the next recalculation (see
//         BlackoilWellModelRescoup::refreshAndSendGroupConstraints()).
//    Finally the master completes its own time step assuming slave rates remain
//    constant; the effective GCW is reset at the start of the next calculation.
//    See excludeAndRedistributeProductionTargets_().
//
//  Phase 2b (injection exclude-and-redistribute): the same three steps per
//    injection phase.  A master group none of whose slave injectors are under
//    group control for a phase (e.g. all on their own BHP limit) has its
//    effective injection GCW for the phase set to 0.  Master groups available
//    for higher-level injection control (GCONINJE item 8 = YES) take part in
//    their parent's injection guide-rate distribution whatever their injection
//    control mode, exactly as GCONPROD item 8 does for production.  See
//    excludeAndRedistributeInjectionTargets_().
//
//  Phase 3: send the resulting (target, cmode, per-rate-type limits)
//    triples to each activated slave via MPI.  Only rank 0 actually sends;
//    other master ranks reach the send helpers but the underlying
//    sendNum/Injection/Production functions in
//    ReservoirCouplingMasterReportStep are gated by `if (comm.rank() == 0)`.
//
//  Consistency across master ranks: the target reductions are per-rank
//    partial sums (a master group's rate is added on rank 0 only).  Each time
//    this calculator recomputes them (updateGCWAndTargetReductions_(),
//    updateInjectionGCWAndTargetReductions_()) it sums them across the master's
//    ranks at once, so every rank computes the same targets and exclusions,
//    and so makes the same collective calls.  All other group rates were
//    already summed by the group-data update that precedes this calculation.
//
// Details on the target calculation (Phase 1):
// --------------------------------------------
// - If the group is a production group:
//  * it is assumed that the action on exceeding the limit (GCONPROD item 7) is "RATE", other
//    actions are not implemented.
//  * the control mode can be:
//   (a) FLD : then item 8 of GCONPROD (available for higher control) is ignored (if it is "NO"),
//       and a guide rate must be defined in item 9 and 10 of GCONPROD.
//   (b) NONE : then item 8 of GCONPROD (available for higher control) must be YES and a guide rate
//       must be defined (or else it will not be possible to distribute a higher level target to the master
//       group since it has no knowledge about the slave group's guide rates).
//   (c) ORAT, WRAT, GRAT, LRAT, CRAT, RESV,... :
//      - if item 8 of GCONPROD is "NO", then the target of the group itself (GCONPROD item 3, 4, ...)
//        is used,
//      - if item 8 of GCONPROD is "YES", then 1) if a higher level group target is found, and a guide
//        rate is defined for the master group, it will be used. 2) If a higher level target is not found,
//        or guiderates are not defined for the master group, then the target of the master group itself
//        is used.
//   - NOTE: If the guiderate definition in the master group (GCONPROD item 10) is different
//           from the control mode of the higher level group target (item 2) in (a), (b), or (c), above
//           then the guide rate is transformed into a guide rate for the phase of the higher level
//           using the production rates of the slave groups as communicated from the slave process at
//           the beginning of the time step.
//  * NOTE: If the group is available for higher level control (item 8 is "YES") and a guide rate
//      is required as noted above, then either:
//      (a)  item 9 of GCONPROD must be set to a positive guide rate value, or
//      (b)  item 9 must be defaulted and item 10 is set to FORM.
//    - In case (b), the formula defined in GUIDERAT will use the communicated slave group
//      potentials to calculate the master group guide rates
//
// - If the group is an injection group:
//  * the control mode for a given phase (OIL, WATER, GAS) can be:
//   (a) FLD : then item 8 (available for higher control) of GCONINJE is ignored (if it is "NO"),
//       and a guide rate must be defined in item 9 and 10 of GCONINJE. Also, a higher level group
//       target for the same phase must be defined.
//   (b) NONE : then item 8 of GCONINJE must be YES and a guide rate must be defined as for (a) above.
//   (c) RATE, RESV, REIN, VREP: then:
//     - if item 8 of GCONINJE is "NO", then the target of the master group itself (GCONINJE item 4
//       or item 5) is used,
//     - if item 8 of GCONINJE is "YES", then 1) if a higher level group target for the same phase is found,
//       and a guide rate is defined for the master group, then the higher level group target will
//       be used. 2) If a higher level target is not found, or a guiderate is not defined
//       for the master group, then the target of the master group itself is used.
//   * NOTE: The RESV, REIN, and VREP targets for a master group depend on slave group reservoir
//       injection rates, surface production rates, or voidage production rate as communicated at
//       the beginning of the time step. See more details in RescoupSendSlaveGroupData.cpp.
//
// Details on the per-rate-type limit calculation (Phase 1):
// ---------------------------------------------------------
// For each non-active rate type (ORAT, WRAT, GRAT, LRAT, RESV) that is not the active
// cmode, the same hierarchy-traversal logic described above is reused with an "explicit
// cmode" parameter.  The difference is in the recursion stopping criterion:
//  - Active target: stops at an ancestor whose production_control() != FLD/NONE.
//  - Per-rate-type limit: stops at an ancestor that has_control(rate_type), i.e. that
//    defines a GCONPROD limit for that specific rate type.
// The guide-rate fraction calculation (FractionCalculator) is identical in both cases,
// since fractions represent proportional capacity allocation independent of rate type.
// If no ancestor defines a limit for a given rate type, the limit is set to -1 (undefined).
// See GroupConstraintCalculator::groupProductionConstraints() for the implementation.
template <class Scalar, class IndexTraits>
void
RescoupConstraintsCalculator<Scalar, IndexTraits>::
calculateMasterGroupConstraintsAndSendToSlaves()
{
    const auto constraints = this->calculateMasterGroupConstraints();

    // Phase 3: send to slaves.  The send functions are rank-0-only internally.
    const auto& rescoup_master = this->reservoir_coupling_master_;
    const auto num_slaves = rescoup_master.numSlaves();
    for (std::size_t slave_idx = 0; slave_idx < num_slaves; ++slave_idx) {
        if (rescoup_master.slaveIsCoupled(slave_idx)) {
            this->sendSlaveGroupConstraintsToSlave_(
                rescoup_master, slave_idx,
                constraints.injection_targets[slave_idx],
                constraints.production_constraints[slave_idx]
            );
        }
    }
}

// The pre-phase and Phases 1, 2 and 2b of
// calculateMasterGroupConstraintsAndSendToSlaves(), see the overview above.
// Separate from the send so that unit tests, which have no slaves to send to,
// can check the constraints (see tests/rescoup/test_rescoup_constraints.cpp).
template <class Scalar, class IndexTraits>
typename RescoupConstraintsCalculator<Scalar, IndexTraits>::MasterGroupConstraints
RescoupConstraintsCalculator<Scalar, IndexTraits>::
calculateMasterGroupConstraints()
{
    // NOTE: Since this object can only be constructed for a master process,
    //   we can be sure that if we are here, we are running as master.
    //
    // The body below must run on all ranks to ensure correct behavior:
    //   - updateGroupControlledWells() uses a collective comm_.sum()
    //   - GroupConstraintCalculator and updateGroupTargetReduction depend on
    //     GroupStateHelper which performs collective operations on the
    //     master's MPI communicator.
    auto& rescoup_master = this->reservoir_coupling_master_;
    GroupConstraintCalculator calculator{
        this->well_model_,
        this->group_state_helper_
    };
    // Reset effective-GCW entries from the previous sync step to GCW=1 before
    // excludeInactiveSlaveMasterGroupsFromDistribution_() (and later the
    // Phase 2 exclusion) repopulate the 0-entries for this step.
    rescoup_master.resetEffectiveGCW();
    this->excludeInactiveSlaveMasterGroupsFromDistribution_();
    // Recompute GCW and reduction rates after the control changes above.
    // The earlier updateAndCommunicateGroupData() in beginTimeStep() may
    // have computed reductions with different controls. Inactive slaves'
    // master groups are now on individual control with zero rates → excluded
    // from guide rate fractions via GCW=0.
    this->updateGCWAndTargetReductions_();
    this->updateInjectionGCWAndTargetReductions_();

    // Phase 1: compute initial targets for all slaves
    const auto num_slaves = rescoup_master.numSlaves();
    std::vector<std::vector<InjectionGroupTarget>> all_injection_targets(num_slaves);
    std::vector<std::vector<ProductionGroupConstraints>> all_production_constraints(num_slaves);
    for (std::size_t slave_idx = 0; slave_idx < num_slaves; ++slave_idx) {
        if (rescoup_master.slaveIsCoupled(slave_idx)) {
            auto [inj, prod] = this->calculateSlaveGroupConstraints_(slave_idx, calculator);
            all_injection_targets[slave_idx] = std::move(inj);
            all_production_constraints[slave_idx] = std::move(prod);
        }
    }

    // Phase 2: exclude master groups without producers under group control and
    // redistribute the shortfall to sibling master groups via the standard
    // localReduction / FractionCalculator machinery.
    this->excludeAndRedistributeProductionTargets_(calculator, all_production_constraints);

    // Phase 2b: exclude master groups without injectors under group control
    // and redistribute the surplus to sibling master groups, per phase.
    this->excludeAndRedistributeInjectionTargets_(calculator, all_injection_targets);

    return {std::move(all_injection_targets), std::move(all_production_constraints)};
}

// ----------------------------------------------------------------------
// Private methods alphabetically for class RescoupConstraintsCalculator
// ----------------------------------------------------------------------

template <class Scalar, class IndexTraits>
std::tuple<
  std::vector<typename RescoupConstraintsCalculator<Scalar, IndexTraits>::InjectionGroupTarget>,
  std::vector<typename RescoupConstraintsCalculator<Scalar, IndexTraits>::ProductionGroupConstraints>
>
RescoupConstraintsCalculator<Scalar, IndexTraits>::
calculateSlaveGroupConstraints_(std::size_t slave_idx, GroupConstraintCalculator<Scalar, IndexTraits>& calculator) const
{
    std::vector<InjectionGroupTarget> injection_targets =
        this->calculateSlaveGroupInjectionTargets_(slave_idx, calculator);
    std::vector<ProductionGroupConstraints> production_constraints;
    auto& rescoup_master = this->reservoir_coupling_master_;
    const auto& master_groups = rescoup_master.getMasterGroupNamesForSlave(slave_idx);
    for (std::size_t group_idx = 0; group_idx < master_groups.size(); ++group_idx) {
        const auto& group_name = master_groups[group_idx];
        const Group& group = this->schedule_.getGroup(group_name, this->report_step_idx_);
        if (group.isProductionGroup()) {
            auto constraints = calculator.groupProductionConstraints(group);
            if (constraints.has_value()) {
                production_constraints.push_back(
                    makeProductionGroupConstraints_(group_idx, *constraints));
            }
        }
    }
    return {injection_targets, production_constraints};
}

template <class Scalar, class IndexTraits>
std::vector<typename RescoupConstraintsCalculator<Scalar, IndexTraits>::InjectionGroupTarget>
RescoupConstraintsCalculator<Scalar, IndexTraits>::
calculateSlaveGroupInjectionTargets_(std::size_t slave_idx, GroupConstraintCalculator<Scalar, IndexTraits>& calculator) const
{
    std::vector<InjectionGroupTarget> injection_targets;
    auto& rescoup_master = this->reservoir_coupling_master_;
    static const std::array<ReservoirCoupling::Phase, 3> phases = {
        ReservoirCoupling::Phase::Water, ReservoirCoupling::Phase::Oil, ReservoirCoupling::Phase::Gas
    };
    const auto& master_groups = rescoup_master.getMasterGroupNamesForSlave(slave_idx);
    for (std::size_t group_idx = 0; group_idx < master_groups.size(); ++group_idx) {
        const auto& group_name = master_groups[group_idx];
        const Group& group = this->schedule_.getGroup(group_name, this->report_step_idx_);
        if (!group.isInjectionGroup()) {
            continue;
        }
        for (ReservoirCoupling::Phase phase : phases) {
            auto target_info = calculator.groupInjectionTarget(group, phase);
            if (target_info.has_value()) {
                // Always send injection targets as RATE. The numeric value is
                // already a surface rate for all modes (RATE, REIN, RESV, VREP),
                // and the slave cannot evaluate derived modes (REIN, VREP, RESV)
                // because it lacks the master's schedule data (reinj_group,
                // voidage_group, GCONSUMP, resv_coeff, etc.).
                injection_targets.push_back(
                    InjectionGroupTarget{
                        group_idx, target_info->constraint,
                        Group::InjectionCMode::RATE, phase
                    }
                );
            }
        }
    }
    return injection_targets;
}

template <class Scalar, class IndexTraits>
void
RescoupConstraintsCalculator<Scalar, IndexTraits>::
excludeAndRedistributeInjectionTargets_(
    GroupConstraintCalculator<Scalar, IndexTraits>& calculator,
    std::vector<std::vector<InjectionGroupTarget>>& all_injection_targets
)
{
    // Drop each master group none of whose slave injectors follow the group's
    // target for a phase from the guide-rate distribution for that phase, and
    // redistribute the surplus to the sibling groups using the standard
    // localReduction / FractionCalculator machinery.
    //
    // This is the master-level counterpart of what happens to a well under
    // group control that cannot meet its share: it switches to its own limit
    // (e.g. BHP), drops out of the guide-rate distribution, and its rate
    // becomes a reduction on the parent's target instead.  On the slave this
    // shows as a group with no injectors under group control, which the slave
    // reports to the master.  Without the exclusion such a master group keeps a
    // guide-rate share it does not inject, the shortfall is never taken up by
    // anyone, and the parent's target is under-delivered.
    //
    // The decision is not based on the slave's injection potentials: for an
    // injector with cross-flow between its connections the potential can be
    // far from what the well injects (see the potential calculation in
    // StandardWell::computeWellPotentialsImplicit()).
    //
    // The steps, and the target the excluded group is sent, are described in
    // the overview above, Phase 2.
    auto& rescoup_master = this->reservoir_coupling_master_;
    const auto num_slaves = all_injection_targets.size();

    // Step 1: exclude the groups without injectors under group control.
    bool newly_excluded = false;
    for (std::size_t slave_idx = 0; slave_idx < num_slaves; ++slave_idx) {
        const auto& master_groups = rescoup_master.getMasterGroupNamesForSlave(slave_idx);
        for (const auto& it : all_injection_targets[slave_idx]) {
            const auto& group_name = master_groups[it.group_name_idx];
            const Phase phase = ReservoirCoupling::convertToOpmPhase(it.phase);
            // A group not available for higher-level control (GCONINJE item 8
            // = NO) takes no part in the distribution in the first place.
            if (!this->group_state_helper_.isMasterGroupEligibleForInjectionGuideRate(group_name, phase)
                || rescoup_master.getSlaveGroupNumGroupControlledInjectors(group_name, it.phase) > 0) {
                continue;
            }
            if (rescoup_master.effectiveInjectionGCW(group_name, phase) != 0) {
                this->deferred_logger_.debug(fmt::format(
                    "RC injection redistribution: {} phase {} has no injectors under group "
                    "control, excluding from guide-rate distribution",
                    group_name, static_cast<int>(it.phase)));
                newly_excluded = true;
            }
            // Drop the group from its siblings' injection guide-rate sum for
            // this phase, so its rate becomes a parent target reduction on the
            // recompute below.
            rescoup_master.setEffectiveInjectionGCW(group_name, phase, 0);
        }
    }

    // If no group was excluded for the first time, the targets already
    // reflect every exclusion in force.
    if (!newly_excluded) {
        return;
    }

    // NOTE: updateGroupTargetReduction() below reads an excluded master group's
    // rate from the slave-reported injection surface rates, so the siblings
    // absorb the surplus over what the group injects now.

    // Step 2: recompute the injection GCW and reductions with the excluded
    // groups excluded.
    this->updateInjectionGCWAndTargetReductions_();

    // Step 3: recompute the targets of all groups.
    for (std::size_t slave_idx = 0; slave_idx < num_slaves; ++slave_idx) {
        const auto& master_groups = rescoup_master.getMasterGroupNamesForSlave(slave_idx);
        for (auto& it : all_injection_targets[slave_idx]) {
            const auto& group_name = master_groups[it.group_name_idx];
            const Scalar old_target = it.target;
            const Group& group = this->schedule_.getGroup(group_name, this->report_step_idx_);
            const auto target_info = calculator.groupInjectionTarget(group, it.phase);
            if (target_info.has_value()) {
                it.target = target_info->constraint;
            }
            this->deferred_logger_.debug(fmt::format(
                "RC injection redistribution: {} phase {} old_target={:.4f} new_target={:.4f}",
                group_name, static_cast<int>(it.phase), old_target, it.target));
        }
    }
}

template <class Scalar, class IndexTraits>
void
RescoupConstraintsCalculator<Scalar, IndexTraits>::
excludeAndRedistributeProductionTargets_(
    GroupConstraintCalculator<Scalar, IndexTraits>& calculator,
    std::vector<std::vector<ProductionGroupConstraints>>& all_production_constraints
)
{
    // Drop each master group none of whose slave producers follow the group's
    // target from the guide-rate distribution, and redistribute the shortfall
    // to the sibling groups using the standard localReduction /
    // FractionCalculator machinery.
    //
    // On the slave a group none of whose producers follow its target shows as
    // a group with no producers under group control, which the slave reports to
    // the master.  The steps, and the target the excluded group is sent, are
    // described in the overview above, Phase 2.
    //
    // The excluded group's rate is the one the slave last reported.  The
    // constraints are recalculated after the slaves' initial well solve of the
    // sync step and at each network iteration (see
    // BlackoilWellModelRescoup::refreshAndSendGroupConstraints()), but not
    // over the Newton iterations, so if an excluded group's rate changes later
    // in the step, e.g. while gas lift adds lift gas, the parent's target is
    // met only to within that change until the next sync step.
    //
    // Implementation: set the effective GCW of the excluded groups to 0 so
    // that updateGroupTargetReduction includes their rates as reduction for the
    // parent group, then recompute the targets of all groups.

    auto& rescoup_master = this->reservoir_coupling_master_;
    const auto num_slaves = all_production_constraints.size();

    // Step 1: exclude the groups without producers under group control
    std::set<std::string> excluded_groups;
    for (std::size_t slave_idx = 0; slave_idx < num_slaves; ++slave_idx) {
        const auto& master_groups = rescoup_master.getMasterGroupNamesForSlave(slave_idx);
        for (const auto& pc : all_production_constraints[slave_idx]) {
            const auto& group_name = master_groups[pc.group_name_idx];
            // A group not available for higher-level control (GCONPROD item 8
            // = NO) takes no part in the distribution in the first place.
            if (!this->group_state_helper_.isMasterGroupEligibleForGuideRate(group_name)
                || rescoup_master.getSlaveGroupNumGroupControlledProducers(group_name) > 0) {
                continue;
            }
            this->deferred_logger_.debug(fmt::format(
                "RC redistribution: {} has no producers under group control, "
                "excluding from guide-rate distribution", group_name));
            excluded_groups.insert(group_name);
            // Drop the group from guide-rate distribution to its siblings: set
            // effective GCW=0 so it is excluded from
            // FractionCalculator::guideRateSum, and its rate is included in the
            // parent's target reduction, on the recompute below.  This must
            // happen before updateGCWAndTargetReductions_().  The control mode
            // is left alone: for an eligible master group the effective GCW,
            // not the control mode, decides both (see
            // GroupStateHelper::updateGroupTargetReductionRecursive_()), and
            // individual control would make the recompute fall back to the
            // group's own GCONPROD target.
            rescoup_master.setEffectiveGCW(group_name, 0);
        }
    }

    // If no groups are excluded, no further action is needed.
    if (excluded_groups.empty()) return;

    // Step 2: recompute GCW and reduction rates with the excluded groups now
    // dropped from the distribution.
    this->updateGCWAndTargetReductions_();

    // Step 3: recompute the targets of all groups using the updated reductions
    for (std::size_t slave_idx = 0; slave_idx < num_slaves; ++slave_idx) {
        const auto& master_groups = rescoup_master.getMasterGroupNamesForSlave(slave_idx);
        for (auto& pc : all_production_constraints[slave_idx]) {
            const auto& group_name = master_groups[pc.group_name_idx];
            const Scalar old_target = pc.target;
            const Group& group = this->schedule_.getGroup(group_name, this->report_step_idx_);
            // The limits are shares of the parent's limits too, and change with
            // the reductions, so replace them together with the target.
            auto constraints = calculator.groupProductionConstraints(group);
            if (constraints.has_value()) {
                pc = makeProductionGroupConstraints_(pc.group_name_idx, *constraints);
            }
            this->deferred_logger_.debug(fmt::format(
                "RC redistribution: {} old_target={:.4f} new_target={:.4f}",
                group_name, old_target, pc.target));
        }
    }

    // Step 3b: share the unused part of each controlling target among the
    // excluded groups below it, instead of offering all of it to each of them.
    this->shareUnusedTargetAmongExcludedGroups_(excluded_groups, all_production_constraints);

    // Step 4: the master completes its own time step, assuming that the
    // production and injection rates of the slave groups remain constant over
    // the time step.  A master group not available for higher-level control
    // (GCONPROD item 8 = NO) is put on the control mode of its allocated target,
    // which is what excludes it from the parent's distribution.  The control
    // mode of the other master groups is left as it is: their share of the
    // parent's target is decided by their effective GCW, and the mode is what
    // the calculator and the master's own group control read, so it must stay
    // the one the Schedule or group control gives, not a marker of this step.
    for (std::size_t slave_idx = 0; slave_idx < num_slaves; ++slave_idx) {
        const auto& master_groups = rescoup_master.getMasterGroupNamesForSlave(slave_idx);
        for (auto& pc : all_production_constraints[slave_idx]) {
            const auto& group_name = master_groups[pc.group_name_idx];
            if (this->group_state_helper_.isMasterGroupEligibleForGuideRate(group_name)) {
                continue;
            }
            this->group_state_helper_.groupState().production_control(
                group_name, pc.cmode);
        }
    }
    // Recompute GCW and reduction with the final state: the groups excluded
    // above have GCW=0, and their slave-reported rates are in the reductions.
    this->updateGCWAndTargetReductions_();
}

// Switch the master groups associated with currently-uncoupled slaves to
// individual control so they are excluded from guide-rate distribution.
// An uncoupled slave contributes zero rate and zero potential, so
// its master groups should not consume any share of the parent's target.
//
// A slave is uncoupled either because it has not activated yet, or because it
// has reached the end of its own schedule while the master run continues.
//
// TODO: Implement GECON item 8, which lets a master deck ask for the master run
//   to stop when one of its slaves finishes, instead of continuing without it.
template <class Scalar, class IndexTraits>
void
RescoupConstraintsCalculator<Scalar, IndexTraits>::
excludeInactiveSlaveMasterGroupsFromDistribution_()
{
    auto& rescoup_master = this->reservoir_coupling_master_;
    const auto num_slaves = rescoup_master.numSlaves();
    // The choice of ORAT here is arbitrary, any non-FLD/NONE control mode works
    // because GroupStateHelper::updateGroupControlledWells() assigns GCW=0 for
    // any master group whose control is not FLD/NONE. ORAT is just a conventional
    // rate mode used here as a marker.  It is only needed for a master group not
    // available for higher-level control (GCONPROD item 8 = NO); the others are
    // excluded through their effective GCW below, and their control mode is left
    // as it is, see Step 4 of excludeAndRedistributeProductionTargets_().
    const Group::ProductionCMode individual_cmode = Group::ProductionCMode::ORAT;
    for (std::size_t slave_idx = 0; slave_idx < num_slaves; ++slave_idx) {
        if (!rescoup_master.slaveIsCoupled(slave_idx)) {
            const auto& master_groups = rescoup_master.getMasterGroupNamesForSlave(slave_idx);
            for (const auto& group_name : master_groups) {
                if (!this->group_state_helper_.isMasterGroupEligibleForGuideRate(group_name)) {
                    this->group_state_helper_.groupState().production_control(
                        group_name, individual_cmode);
                }
                // For a participating master group (GCONPROD item 8 = YES) the GCW
                // comes from effectiveGCW(), not from the control mode.  Force it
                // to 0 here so an inactive slave's groups are excluded from
                // guide-rate distribution regardless of item 8.
                rescoup_master.setEffectiveGCW(group_name, 0);
            }
        }
    }
    // Likewise for the injection guide-rate distribution of every phase.
    this->excludeInactiveSlaveMasterGroupsFromInjectionDistribution_();
}

// Exclude the master groups of currently-uncoupled slaves from the injection
// guide-rate distribution of every phase, by setting their effective injection
// GCW to 0.  The injection part of
// excludeInactiveSlaveMasterGroupsFromDistribution_().
template <class Scalar, class IndexTraits>
void
RescoupConstraintsCalculator<Scalar, IndexTraits>::
excludeInactiveSlaveMasterGroupsFromInjectionDistribution_()
{
    auto& rescoup_master = this->reservoir_coupling_master_;
    const auto num_slaves = rescoup_master.numSlaves();
    for (std::size_t slave_idx = 0; slave_idx < num_slaves; ++slave_idx) {
        if (!rescoup_master.slaveIsCoupled(slave_idx)) {
            for (const auto& group_name : rescoup_master.getMasterGroupNamesForSlave(slave_idx)) {
                for (const Phase phase : {Phase::WATER, Phase::OIL, Phase::GAS}) {
                    rescoup_master.setEffectiveInjectionGCW(group_name, phase, 0);
                }
            }
        }
    }
}

// The constraints sent to a slave for one master group, from the target and
// limits the calculator found for it.
template <class Scalar, class IndexTraits>
typename RescoupConstraintsCalculator<Scalar, IndexTraits>::ProductionGroupConstraints
RescoupConstraintsCalculator<Scalar, IndexTraits>::
makeProductionGroupConstraints_(
    std::size_t group_idx,
    const typename GroupConstraintCalculator<Scalar, IndexTraits>::ProductionConstraintResult& constraints)
{
    return ProductionGroupConstraints{
        group_idx,
        constraints.active_target,
        constraints.active_cmode,
        constraints.oil_limit,
        constraints.water_limit,
        constraints.gas_limit,
        constraints.liquid_limit,
        constraints.resv_limit
    };
}

// Step 3b of excludeAndRedistributeProductionTargets_().
//
// Step 3 gives an excluded group the target it would get on returning to group
// control on its own: the controlling ancestor's target less the current rates
// of all the other groups below it.  That is not a consistent allocation once
// two or more groups are excluded: each of them is offered the whole unused
// part of the ancestor's target.  As long as their wells sit at their own (THP
// or BHP) limits this does no harm, but anything that lets them produce more --
// gas lift, or just their rates rising over a long sync step -- then works
// towards a sum of targets above the ancestor's target.  Without reservoir
// coupling the ancestor's own constraint check would catch this within the
// time step; here the master only sees it at the next sync step.
//
// So, for each controlling ancestor and rate type, the unused part of its
// target is worked out -- its target less its current rate, less the room the
// participating master groups below it still have up to their targets -- and
// shared among the excluded groups below it.  Each excluded group gets its
// current rate plus its share, and never more than the target from Step 3.  The
// shares are in proportion to how far each group's Step 3 target exceeds its
// rate.  A negative unused part (the ancestor produces above its target) is
// shared in proportion to the groups' rates.  When participating groups are
// present they take up the unused part, so the excluded groups are held at
// their current rates; only when none is present is the room shared among the
// excluded groups.
//
// The price is that an excluded group's target is often its current rate, so a
// small rise in its rate switches its wells to group control, and from the next
// recalculation it takes part in the distribution again.  If its wells cannot
// produce the share it then gets, the parent's target is under-delivered until
// they are back on their own limits and the group is excluded again.
//
// The same is done for the per-rate-type limits (ORAT, WRAT, GRAT, LRAT), each
// with the ancestor that defines that limit.  RESV targets and limits are left
// as Step 3 computed them, as are the master's own wells: their room up to
// their targets is not deducted from the unused part.
template <class Scalar, class IndexTraits>
void
RescoupConstraintsCalculator<Scalar, IndexTraits>::
shareUnusedTargetAmongExcludedGroups_(
    const std::set<std::string>& excluded_groups,
    std::vector<std::vector<ProductionGroupConstraints>>& all_production_constraints)
{
    using CMode = Group::ProductionCMode;
    auto& rescoup_master = this->reservoir_coupling_master_;
    const auto& group_state = this->group_state_helper_.groupState();

    auto rate_for_mode = [this](const std::string& name, const CMode mode) -> Scalar {
        const auto rates = this->group_state_helper_.getProductionGroupRateVector(name);
        switch (mode) {
        case CMode::ORAT: return rates.oil_rat;
        case CMode::WRAT: return rates.wat_rat;
        case CMode::GRAT: return rates.gas_rat;
        case CMode::LRAT: return rates.oil_rat + rates.wat_rat;
        default: return Scalar{0};
        }
    };
    // The value of a constraint for a rate type: the active target, or one of
    // the per-rate-type limits (-1 when there is none).
    auto value_for_mode = [](ProductionGroupConstraints& pc, const CMode mode) -> Scalar* {
        if (pc.cmode == mode) {
            return &pc.target;
        }
        switch (mode) {
        case CMode::ORAT: return &pc.oil_limit;
        case CMode::WRAT: return &pc.water_limit;
        case CMode::GRAT: return &pc.gas_limit;
        case CMode::LRAT: return &pc.liquid_limit;
        default: return nullptr;
        }
    };
    // The ancestor whose target or limit for the rate type the group's
    // constraint is a share of, as GroupConstraintCalculator finds it, together
    // with the product of the efficiency factors from the group up to it.
    auto controlling_ancestor = [this, &group_state](const std::string& name, const CMode mode,
                                                     const bool is_target)
        -> std::optional<std::pair<std::string, Scalar>>
    {
        const Group* group = &this->schedule_.getGroup(name, this->report_step_idx_);
        Scalar efficiency = 1.0;
        while (!group->is_field() && group->productionGroupControlAvailable()) {
            efficiency *= group->getGroupEfficiencyFactor();
            group = &this->schedule_.getGroup(group->parent(), this->report_step_idx_);
            const bool controls = is_target
                ? !group_state.has_field_or_none_control(group->name())
                : group->has_control(mode);
            if (controls) {
                return std::make_pair(group->name(), efficiency);
            }
        }
        return std::nullopt;
    };
    auto is_below = [this](std::string name, const std::string& ancestor) {
        while (name != "FIELD") {
            name = this->schedule_.getGroup(name, this->report_step_idx_).parent();
            if (name == ancestor) {
                return true;
            }
        }
        return false;
    };

    struct Member {
        Scalar* value;       // the constraint to adjust
        Scalar step3_value;  // its value from Step 3
        Scalar rate;         // the group's current rate for the rate type
        Scalar efficiency;   // efficiency factor from the group up to the ancestor
    };
    for (const CMode mode : {CMode::ORAT, CMode::WRAT, CMode::GRAT, CMode::LRAT}) {
        // Excluded groups by controlling ancestor
        std::map<std::string, std::vector<Member>> excluded_by_ancestor;
        for (std::size_t slave_idx = 0; slave_idx < all_production_constraints.size(); ++slave_idx) {
            const auto& master_groups = rescoup_master.getMasterGroupNamesForSlave(slave_idx);
            for (auto& pc : all_production_constraints[slave_idx]) {
                const auto& name = master_groups[pc.group_name_idx];
                if (excluded_groups.count(name) == 0) {
                    continue;
                }
                Scalar* value = value_for_mode(pc, mode);
                if (value == nullptr || *value < 0) {
                    continue;
                }
                const auto ancestor = controlling_ancestor(name, mode, /*is_target=*/pc.cmode == mode);
                // A group with an efficiency factor of 0 takes nothing of what the ancestor
                // has left, and the shares below are divided by the factor.
                if (!ancestor.has_value() || !(ancestor->second > 0)) {
                    continue;
                }
                excluded_by_ancestor[ancestor->first].push_back(
                    Member{value, *value, rate_for_mode(name, mode), ancestor->second});
            }
        }
        for (auto& [ancestor, members] : excluded_by_ancestor) {
            const Group& top = this->schedule_.getGroup(ancestor, this->report_step_idx_);
            // How far the ancestor is below its target; negative when it is above it.
            const Scalar slack = this->group_state_helper_.getProductionGroupTargetForMode(top, mode)
                - rate_for_mode(ancestor, mode);
            // The room the participating master groups below the ancestor still
            // have up to their own targets or limits is not unused.
            Scalar room = 0.0;
            for (std::size_t slave_idx = 0; slave_idx < all_production_constraints.size(); ++slave_idx) {
                const auto& master_groups = rescoup_master.getMasterGroupNamesForSlave(slave_idx);
                for (auto& pc : all_production_constraints[slave_idx]) {
                    const auto& name = master_groups[pc.group_name_idx];
                    if (excluded_groups.count(name) > 0
                        || !this->group_state_helper_.isMasterGroupEligibleForGuideRate(name)
                        || !is_below(name, ancestor)) {
                        continue;
                    }
                    const Scalar* value = value_for_mode(pc, mode);
                    if (value == nullptr || *value < 0) {
                        continue;
                    }
                    const auto chain = controlling_ancestor(name, mode, /*is_target=*/pc.cmode == mode);
                    const Scalar efficiency = chain.has_value() ? chain->second : Scalar{1};
                    room += efficiency * std::max(Scalar{0}, *value - rate_for_mode(name, mode));
                }
            }
            // An ancestor below its target hands what the participating groups do not
            // take to the excluded groups, and never takes anything from them: they are
            // held at their current rates when the participating groups take it all.  Only
            // an ancestor above its target takes from the excluded groups.
            const Scalar unused = (slack >= 0) ? std::max(Scalar{0}, slack - room) : slack;
            Scalar weight_sum = 0.0;
            for (const auto& m : members) {
                weight_sum += (unused >= 0)
                    ? m.efficiency * std::max(Scalar{0}, m.step3_value - m.rate)
                    : m.efficiency * m.rate;
            }
            for (auto& m : members) {
                Scalar share = 0.0;
                if (weight_sum > 0) {
                    const Scalar weight = (unused >= 0)
                        ? m.efficiency * std::max(Scalar{0}, m.step3_value - m.rate)
                        : m.efficiency * m.rate;
                    share = unused * weight / weight_sum / m.efficiency;
                }
                *m.value = std::clamp(m.rate + share, Scalar{0}, m.step3_value);
            }
            this->deferred_logger_.debug(fmt::format(
                "RC redistribution: {} unused {} {:.4f} shared among {} excluded groups",
                ancestor, Group::ProductionCMode2String(mode), unused, members.size()));
        }
    }
}

template <class Scalar, class IndexTraits>
void
RescoupConstraintsCalculator<Scalar, IndexTraits>::
sendSlaveGroupConstraintsToSlave_(
    const ReservoirCouplingMaster<Scalar>& rescoup_master,
    std::size_t slave_idx,
    const std::vector<InjectionGroupTarget>& injection_targets,
    const std::vector<ProductionGroupConstraints>& production_constraints
) const
{
    auto num_injection_targets = injection_targets.size();
    auto num_production_constraints = production_constraints.size();
    // First, send the number of constraints such that the slave can know if it can expect none
    // or more constraints.
    rescoup_master.sendNumGroupConstraintsToSlave(slave_idx, num_injection_targets, num_production_constraints);
    if (num_injection_targets > 0) {
        rescoup_master.sendInjectionTargetsToSlave(slave_idx, injection_targets);
    }
    if (num_production_constraints > 0) {
        rescoup_master.sendProductionConstraintsToSlave(slave_idx, production_constraints);
    }
}

// Recompute GCW (Group Controlled Wells) and the FIELD-level target reduction
// after master-group control modes have been changed.  GCW must be updated
// before the reduction since updateGroupTargetReduction() depends on it.
template <class Scalar, class IndexTraits>
void
RescoupConstraintsCalculator<Scalar, IndexTraits>::
updateGCWAndTargetReductions_()
{
    this->group_state_helper_.updateGroupControlledWells(
        /*is_production_group=*/true, /*dummy_injection_phase=*/Phase::OIL);
    const Group& fieldGroup = this->schedule_.getGroup("FIELD", this->report_step_idx_);
    this->group_state_helper_.updateGroupTargetReduction(fieldGroup, /*is_injector=*/false);
    // The reduction rates are per-rank partial sums (a master group's rate is
    // added on rank 0 only); sum them so every rank computes the same targets
    // and exclusions, see updateInjectionGCWAndTargetReductions_().
    this->group_state_helper_.groupState().communicate_reduction_rates(
        this->group_state_helper_.comm(), /*is_injector=*/false);
}

// Recompute the injection GCW for every phase and the FIELD-level injection
// target reduction after effective injection GCW entries have changed.
template <class Scalar, class IndexTraits>
void
RescoupConstraintsCalculator<Scalar, IndexTraits>::
updateInjectionGCWAndTargetReductions_()
{
    for (const Phase phase : {Phase::WATER, Phase::OIL, Phase::GAS}) {
        this->group_state_helper_.updateGroupControlledWells(
            /*is_production_group=*/false, phase);
    }
    const Group& fieldGroup = this->schedule_.getGroup("FIELD", this->report_step_idx_);
    this->group_state_helper_.updateGroupTargetReduction(fieldGroup, /*is_injector=*/true);
    // updateGroupTargetReduction() leaves per-rank partial sums: a master
    // group's slave-reported rate is added on rank 0 only.  Sum them across the
    // master's ranks before any target is computed from them.  Otherwise the
    // targets, and with them the exclusions, can differ between ranks, and
    // the ranks then make different numbers of the collective calls inside the
    // target computation (GroupStateHelper::getGroupRatesAvailableForHigherLevelControl()).
    this->group_state_helper_.groupState().communicate_reduction_rates(
        this->group_state_helper_.comm(), /*is_injector=*/true);
}

template class RescoupConstraintsCalculator<double, BlackOilDefaultFluidSystemIndices>;

#if FLOW_INSTANTIATE_FLOAT
template class RescoupConstraintsCalculator<float, BlackOilDefaultFluidSystemIndices>;
#endif

}// namespace Opm
