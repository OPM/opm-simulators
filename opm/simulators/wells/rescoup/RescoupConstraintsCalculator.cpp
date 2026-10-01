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
//    is derived from the control mode.  Master groups, however, carry the individual control mode
//    purely as a slave-communication signal (see the Phase 2 finalize), so their cmode cannot
//    be used to decide participation. The master therefore maintains a per-master-group "effective GCW"
//    (ReservoirCouplingMaster::effectiveGCW), set in this routine and read by
//    GroupStateHelper::updateGroupControlledWellsRecursive_.  It separates two
//    concerns:
//      - static eligibility: GCONPROD item 8 RESPOND_TO_PARENT = YES
//        (Group::productionGroupControlAvailable) — a deck property; from
//      - dynamic weight this sync step: 1 if the group participates and is
//        uncapped, 0 if capped at its slave potential (Phase 2) or if its
//        slave is inactive (pre-phase).
//    An eligible group with effective GCW = 0 is excluded from the guide-rate
//    sum and contributes its rate as a parent target reduction instead — which
//    is precisely how a capped group's shortfall is redistributed to its
//    siblings.  This decoupling is needed because the cmode (the natural
//    dynamic control channel) is already taken by the slave-communication
//    marker for master groups.
//
//  Pre-phase: restore master-group control modes from the Schedule, reset the
//    effective GCW (participating master groups default to 1), switch inactive
//    slaves' master groups to individual control and set their effective GCW
//    to 0 so they don't consume any share of the parent's target, then
//    recompute GCW + FIELD target reductions with the resulting state.
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
//  Phase 2 (cap-and-redistribute): if any master group's Phase-1 target
//    exceeds its slave's reported production potential, cap the target at
//    that potential and set the group's effective GCW to 0.  With effective
//    GCW = 0 the group drops out of the parent's guide-rate sum, and on the
//    recompute below its (capped) rate is subtracted as a FIELD target
//    reduction instead; the remaining uncapped groups then absorb the
//    shortfall through the standard localReduction / FractionCalculator
//    machinery.  (The capped group is also switched to individual control,
//    but for master groups the cmode is only a slave-communication marker —
//    it is the effective GCW, not the cmode, that drives the guide-rate
//    exclusion; see the effective-GCW note above.)  Finally all master groups
//    are switched to individual control so the master completes its own time
//    step assuming slave rates remain constant; the controls and effective
//    GCW are reset at the start of the next sync step.  See
//    capAndRedistributeProductionTargets_().
//
//  Phase 2b (injection cap-and-redistribute): the same idea per injection
//    phase.  A master group whose injection target for a phase exceeds the
//    injection potential its slave reported for that phase is capped at the
//    potential and its effective injection GCW for the phase is set to 0, so
//    its rate is subtracted from the parent's target instead of it holding a
//    guide-rate share, and the remaining groups absorb the surplus.  Master
//    groups available for higher-level injection control (GCONINJE item 8 =
//    YES) take part in their parent's injection guide-rate distribution
//    whatever their injection control mode, exactly as GCONPROD item 8 does
//    for production.  See capAndRedistributeInjectionTargets_().
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
//    ranks at once, so every rank computes the same targets and cap decisions,
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
    // Phase 2 cap) repopulate the 0-entries for this step.
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

    // Phase 2: cap production targets at slave potentials and redistribute
    // shortfall to sibling master groups via the standard localReduction /
    // FractionCalculator machinery.
    this->capAndRedistributeProductionTargets_(calculator, all_production_constraints);

    // Phase 2b: cap injection targets at slave injection potentials and
    // redistribute the surplus to sibling master groups, per phase.
    this->capAndRedistributeInjectionTargets_(calculator, all_injection_targets);

    // Phase 3: send to slaves.  The send functions are rank-0-only internally.
    for (std::size_t slave_idx = 0; slave_idx < num_slaves; ++slave_idx) {
        if (rescoup_master.slaveIsCoupled(slave_idx)) {
            this->sendSlaveGroupConstraintsToSlave_(
                rescoup_master, slave_idx,
                all_injection_targets[slave_idx],
                all_production_constraints[slave_idx]
            );
        }
    }
}

// Recompute the injection targets against the slave rates the master holds now
// and ship them to the slaves, replacing the targets sent earlier in this sync
// step.  See the declaration in the header for why this second send exists.
template <class Scalar, class IndexTraits>
void
RescoupConstraintsCalculator<Scalar, IndexTraits>::
recalculateInjectionTargetsAndSendToSlaves()
{
    // As in calculateMasterGroupConstraintsAndSendToSlaves(), the body must run
    // on every rank of the master communicator: GroupConstraintCalculator relies
    // on GroupStateHelper, which performs collective operations.  The MPI sends
    // themselves are rank-0-only inside the send helpers.
    auto& rescoup_master = this->reservoir_coupling_master_;
    GroupConstraintCalculator calculator{
        this->well_model_,
        this->group_state_helper_
    };
    // Decide the injection caps afresh from the latest slave data, exactly as
    // the first calculation of the sync step does.  A cap kept from earlier in
    // the sync step would leave its group with effective injection GCW 0, and
    // such a group drops out of its own guide-rate share: its target would come
    // out as 0, which the cap test (target > potential) cannot detect.  Only the
    // injection caps are reset; the production caps decided earlier in the sync
    // step stay in force.
    rescoup_master.resetEffectiveInjectionGCW();
    this->excludeInactiveSlaveMasterGroupsFromInjectionDistribution_();
    // The slave injection rates have also changed since the targets were first
    // computed this sync step, so refresh the injection target reductions.
    this->updateInjectionGCWAndTargetReductions_();
    const auto num_slaves = rescoup_master.numSlaves();
    std::vector<std::vector<InjectionGroupTarget>> all_injection_targets(num_slaves);
    for (std::size_t slave_idx = 0; slave_idx < num_slaves; ++slave_idx) {
        if (rescoup_master.slaveIsCoupled(slave_idx)) {
            all_injection_targets[slave_idx] =
                this->calculateSlaveGroupInjectionTargets_(slave_idx, calculator);
        }
    }
    this->capAndRedistributeInjectionTargets_(calculator, all_injection_targets);
    for (std::size_t slave_idx = 0; slave_idx < num_slaves; ++slave_idx) {
        if (!rescoup_master.slaveIsCoupled(slave_idx)) {
            continue;
        }
        // An empty production-constraint list tells the slave that no production
        // constraints follow; the ones it received earlier this sync step stay in
        // force.
        this->sendSlaveGroupConstraintsToSlave_(
            rescoup_master, slave_idx, all_injection_targets[slave_idx],
            /*production_constraints=*/{}
        );
    }
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
                    ProductionGroupConstraints{
                        group_idx,
                        constraints->active_target,
                        constraints->active_cmode,
                        constraints->oil_limit,
                        constraints->water_limit,
                        constraints->gas_limit,
                        constraints->liquid_limit,
                        constraints->resv_limit
                    }
                );
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
capAndRedistributeInjectionTargets_(
    GroupConstraintCalculator<Scalar, IndexTraits>& calculator,
    std::vector<std::vector<InjectionGroupTarget>>& all_injection_targets
)
{
    // Cap each master group's injection target for a phase at its slave's
    // reported injection potential for that phase, then redistribute the
    // surplus to the sibling groups using the standard localReduction /
    // FractionCalculator machinery.
    //
    // Without the cap a master group keeps a guide-rate share it cannot inject:
    // the slave's wells run into their own limits, the shortfall is never taken
    // up by anyone, and the parent's target is under-delivered.  This is the
    // master-level counterpart of what happens to a well under group control
    // that cannot meet its share: it drops out of the guide-rate distribution
    // and its rate becomes a reduction on the parent's target instead.
    //
    // The steps are repeated until a pass caps nothing new: capping one group
    // raises its siblings' shares, which can push another sibling over its own
    // potential.  Every further pass caps at least one more (group, phase)
    // pair, so the loop ends after at most as many passes as there are targets.
    auto& rescoup_master = this->reservoir_coupling_master_;
    const auto num_slaves = all_injection_targets.size();
    std::set<std::pair<std::string, ReservoirCoupling::Phase>> capped;
    while (true) {
        // Step 1: identify targets that exceed the potential and cap them.
        bool newly_capped = false;
        for (std::size_t slave_idx = 0; slave_idx < num_slaves; ++slave_idx) {
            const auto& master_groups = rescoup_master.getMasterGroupNamesForSlave(slave_idx);
            for (auto& it : all_injection_targets[slave_idx]) {
                const auto& group_name = master_groups[it.group_name_idx];
                if (capped.count({group_name, it.phase}) > 0) {
                    continue;  // already capped in an earlier pass
                }
                const auto& potentials = rescoup_master.getSlaveGroupInjectionPotentials(group_name);
                const Scalar pot = this->potentialForInjectionPhase_(potentials, it.phase);
                if (it.target > pot) {
                    this->deferred_logger_.debug(fmt::format(
                        "RC injection redistribution: {} phase {} target={:.4f} exceeds "
                        "potential={:.4f}, capping",
                        group_name, static_cast<int>(it.phase), it.target, pot));
                    it.target = pot;
                    capped.emplace(group_name, it.phase);
                    newly_capped = true;
                    // Drop the capped group from its siblings' injection guide-rate
                    // sum for this phase, so its rate becomes a parent target
                    // reduction on the recompute below.
                    rescoup_master.setEffectiveInjectionGCW(
                        group_name, ReservoirCoupling::convertToOpmPhase(it.phase), 0);
                }
            }
        }

        // If no group was capped in this pass, every target is within its
        // potential.
        if (!newly_capped) {
            return;
        }

        // NOTE: updateGroupTargetReduction() below reads a capped master group's
        // rate from the slave-reported injection surface rates.  As for production
        // (see capAndRedistributeProductionTargets_()), a group whose target exceeds
        // its potential has its injectors at their own limits, so the reported rate
        // is close to the potential the target was capped at.

        // Step 2: recompute the injection GCW and reductions with the capped groups
        // excluded.
        this->updateInjectionGCWAndTargetReductions_();

        // Step 3: recompute the targets of the uncapped groups.
        for (std::size_t slave_idx = 0; slave_idx < num_slaves; ++slave_idx) {
            const auto& master_groups = rescoup_master.getMasterGroupNamesForSlave(slave_idx);
            for (auto& it : all_injection_targets[slave_idx]) {
                const auto& group_name = master_groups[it.group_name_idx];
                if (capped.count({group_name, it.phase}) > 0) {
                    continue;  // keep the capped target
                }
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
}

template <class Scalar, class IndexTraits>
void
RescoupConstraintsCalculator<Scalar, IndexTraits>::
capAndRedistributeProductionTargets_(
    GroupConstraintCalculator<Scalar, IndexTraits>& calculator,
    std::vector<std::vector<ProductionGroupConstraints>>& all_production_constraints
)
{
    // Cap each master group's production target at its slave's reported
    // potential, then redistribute the shortfall to sibling groups using
    // the standard localReduction / FractionCalculator machinery.
    //
    // Implementation: switch capped groups to individual control so that
    // updateGroupTargetReduction includes their rates as reduction for the
    // parent group.  Then recompute targets for uncapped groups — they get
    // a larger share because the parent's available target increased.

    auto& rescoup_master = this->reservoir_coupling_master_;
    const auto num_slaves = all_production_constraints.size();

    // Step 1: identify groups where target > potential and cap them
    std::set<std::string> capped_groups;
    for (std::size_t slave_idx = 0; slave_idx < num_slaves; ++slave_idx) {
        const auto& master_groups = rescoup_master.getMasterGroupNamesForSlave(slave_idx);
        for (auto& pc : all_production_constraints[slave_idx]) {
            const auto& group_name = master_groups[pc.group_name_idx];
            const auto& potentials = rescoup_master.getSlaveGroupPotentials(group_name);
            const auto pot = this->potentialForProductionCmode_(
                potentials, pc.cmode);
            if (pot >= 0 && pc.target > pot) {
                pc.target = pot;
                capped_groups.insert(group_name);
                // Switch to individual control so updateGroupTargetReduction()
                // includes this group's rate as reduction for the parent.
                this->group_state_helper_.groupState().production_control(
                    group_name, pc.cmode);
                // Drop the capped group from guide-rate distribution to its
                // siblings: set effective GCW=0 so it is excluded from
                // FractionCalculator::guideRateSum on the recompute below.  This
                // must happen before updateGCWAndTargetReductions_().
                rescoup_master.setEffectiveGCW(group_name, 0);
            }
        }
    }

    // If no groups are capped, no further action is needed.
    if (capped_groups.empty()) return;

    // Note: updateGroupTargetReduction() below calls sumWellSurfaceRates(), which
    // for master groups reads from the MPI-received surface_rates in
    // slave_group_production_data_ rather than from groupState().production_rates().
    // No override is needed here: the cap condition (target > potential) implies
    // the slave cannot meet the target, so its wells saturate at BHP and the
    // group produces at its full potential.  Both surface_rates and potentials
    // are reported by the slave from the same reservoir-state snapshot, so
    // surface_rates ≈ potential for any capped group, i.e. the reduction computed
    // below already reflects the capped rate.

    // Step 2: recompute GCW and reduction rates with capped groups now
    // individually controlled.
    this->updateGCWAndTargetReductions_();

    // Step 3: recompute targets for uncapped groups using the updated reductions
    for (std::size_t slave_idx = 0; slave_idx < num_slaves; ++slave_idx) {
        const auto& master_groups = rescoup_master.getMasterGroupNamesForSlave(slave_idx);
        for (auto& pc : all_production_constraints[slave_idx]) {
            const auto& group_name = master_groups[pc.group_name_idx];
            if (capped_groups.count(group_name) > 0) {
                this->deferred_logger_.debug(fmt::format(
                    "RC redistribution: {} capped at potential, target={:.4f}",
                    group_name, pc.target));
                continue;  // keep the capped target
            }
            const Scalar old_target = pc.target;
            const Group& group = this->schedule_.getGroup(group_name, this->report_step_idx_);
            auto constraints = calculator.groupProductionConstraints(group);
            if (constraints.has_value()) {
                pc.target = constraints->active_target;
            }
            this->deferred_logger_.debug(fmt::format(
                "RC redistribution: {} old_target={:.4f} new_target={:.4f}",
                group_name, old_target, pc.target));
        }
    }

    // Step 4: switch ALL master groups to individual control with their
    // allocated targets. I.e. the master completes its own time step, assuming
    // that the production and injection rates of the slave groups remain constant
    // over the time step.
    for (std::size_t slave_idx = 0; slave_idx < num_slaves; ++slave_idx) {
        const auto& master_groups = rescoup_master.getMasterGroupNamesForSlave(slave_idx);
        for (auto& pc : all_production_constraints[slave_idx]) {
            const auto& group_name = master_groups[pc.group_name_idx];
            if (capped_groups.count(group_name) == 0) {
                // Uncapped groups: switch to individual control too
                this->group_state_helper_.groupState().production_control(
                    group_name, pc.cmode);
            }
        }
    }
    // Recompute GCW and reduction: all master groups are now individually
    // controlled, so GCW=0 for all of them and FIELD's reduction sums every
    // master group's slave-reported rate.
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
    // rate mode used here as a marker.
    const Group::ProductionCMode individual_cmode = Group::ProductionCMode::ORAT;
    for (std::size_t slave_idx = 0; slave_idx < num_slaves; ++slave_idx) {
        if (!rescoup_master.slaveIsCoupled(slave_idx)) {
            const auto& master_groups = rescoup_master.getMasterGroupNamesForSlave(slave_idx);
            for (const auto& group_name : master_groups) {
                this->group_state_helper_.groupState().production_control(
                    group_name, individual_cmode);
                // For a participating master group (GCONPROD item 8 = YES) the
                // ORAT marker above no longer forces GCW=0 — its GCW now comes
                // from effectiveGCW().  Force it to 0 here so an inactive slave's
                // groups are excluded from guide-rate distribution regardless of
                // item 8.
                rescoup_master.setEffectiveGCW(group_name, 0);
            }
        }
    }
    // Likewise for the injection guide-rate distribution of every phase.
    this->excludeInactiveSlaveMasterGroupsFromInjectionDistribution_();
}

// Exclude the master groups of currently-uncoupled slaves from the injection
// guide-rate distribution of every phase, by setting their effective injection
// GCW to 0.  Split out from excludeInactiveSlaveMasterGroupsFromDistribution_()
// so the injection caps can be recomputed within a sync step without touching
// the production state.
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

template <class Scalar, class IndexTraits>
Scalar
RescoupConstraintsCalculator<Scalar, IndexTraits>::
potentialForInjectionPhase_(
    const ReservoirCoupling::Potentials<Scalar>& potentials,
    ReservoirCoupling::Phase phase) const
{
    // The potentials are stored in display units by GuideRateHandler, so
    // convert to SI for comparison with the (SI) target from
    // GroupConstraintCalculator.
    const auto& units = this->schedule_.getUnits();
    const auto measure = (phase == ReservoirCoupling::Phase::Gas)
        ? UnitSystem::measure::gas_surface_rate
        : UnitSystem::measure::liquid_surface_rate;
    return units.to_si(measure, potentials[phase]);
}

template <class Scalar, class IndexTraits>
Scalar
RescoupConstraintsCalculator<Scalar, IndexTraits>::
potentialForProductionCmode_(
    const ReservoirCoupling::Potentials<Scalar>& potentials,
    Group::ProductionCMode cmode) const
{
    // The potentials are stored in display units (e.g., SM3/DAY for METRIC)
    // by GuideRateHandler, so we convert to SI (m3/s) for comparison with
    // the target from GroupConstraintCalculator which is in SI.
    using RcPhase = ReservoirCoupling::Phase;
    const auto& units = this->schedule_.getUnits();
    switch (cmode) {
    case Group::ProductionCMode::ORAT:
        return units.to_si(UnitSystem::measure::liquid_surface_rate,
                           potentials[RcPhase::Oil]);
    case Group::ProductionCMode::GRAT:
        return units.to_si(UnitSystem::measure::gas_surface_rate,
                           potentials[RcPhase::Gas]);
    case Group::ProductionCMode::WRAT:
        return units.to_si(UnitSystem::measure::liquid_surface_rate,
                           potentials[RcPhase::Water]);
    case Group::ProductionCMode::LRAT:
        return units.to_si(UnitSystem::measure::liquid_surface_rate,
                           potentials[RcPhase::Oil] + potentials[RcPhase::Water]);
    default:
        return -1;  // No capping for RESV, FLD, NONE, etc.
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
    // and cap decisions, see updateInjectionGCWAndTargetReductions_().
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
    // targets, and with them the cap decisions, can differ between ranks, and
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
