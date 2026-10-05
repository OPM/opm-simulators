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

#ifndef OPM_RESCOUP_CONSTRAINTS_CALCULATOR_HPP
#define OPM_RESCOUP_CONSTRAINTS_CALCULATOR_HPP
#include <opm/input/eclipse/Schedule/Group/GuideRate.hpp>
#include <opm/material/fluidsystems/PhaseUsageInfo.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCoupling.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingMaster.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/BlackoilWellModelGeneric.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include <opm/simulators/wells/GroupStateHelper.hpp>
#include <opm/simulators/wells/GroupConstraintCalculator.hpp>
#include <opm/simulators/wells/GuideRateHandler.hpp>
#include <opm/simulators/wells/WellState.hpp>

namespace Opm {

/// @brief Computes per-master-group production targets and per-rate-type
///   limits for reservoir coupling, and sends them to the slaves.
///
/// @details Constructed once per sync step on the master process, drives
///   the full target-distribution flow when its single public entry point
///   `calculateMasterGroupConstraintsAndSendToSlaves()` is called.  The
///   flow is structured as a pre-phase (control restore + inactive-slave
///   handling + GCW/reduction recompute) followed by three phases:
///   target computation (Phase 1), cap-and-redistribute (Phase 2), and
///   slave send (Phase 3).  See the function-level comment on
///   `calculateMasterGroupConstraintsAndSendToSlaves()` in the
///   corresponding .cpp file for the per-phase narrative.
///
/// @tparam Scalar Floating-point type used for rates and targets.
/// @tparam IndexTraits Phase-index traits.
template<class Scalar, class IndexTraits>
class RescoupConstraintsCalculator {
public:
    using InjectionGroupTarget = ReservoirCoupling::InjectionGroupTarget<Scalar>;
    using ProductionGroupConstraints = ReservoirCoupling::ProductionGroupConstraints<Scalar>;

    /// @brief Construct a calculator bound to the master's per-sync-step
    ///   state.
    /// @param guide_rate_handler Source of guide rates and the deferred
    ///   logger used during the constraint calculation.
    /// @param group_state_helper Provides access to the group state, the
    ///   schedule, the well state, and the reservoir-coupling master
    ///   facade.  All intermediate state changes (control modes, GCW,
    ///   target reductions) are written through this helper.
    RescoupConstraintsCalculator(
        GuideRateHandler<Scalar, IndexTraits>& guide_rate_handler,
        GroupStateHelper<Scalar, IndexTraits>& group_state_helper
    );

    /// @brief Run the full master-side target distribution and dispatch
    ///   the resulting constraints to each activated slave.
    /// @details Runs on every rank of the master communicator; the
    ///   slave-facing MPI sends are rank-0-only inside the underlying
    ///   helpers.  The driver writes group-state side effects (cmodes,
    ///   GCW, reductions) and finishes with a
    ///   `GroupState::communicate_rates(comm)` call so the
    ///   post-redistribution state is consistent across master ranks.
    ///   See the function-level comment on the implementation in the
    ///   corresponding .cpp file for the per-phase walkthrough and the
    ///   collective-call invariants.
    void calculateMasterGroupConstraintsAndSendToSlaves();

private:
    /// @brief Phase 1: compute initial guide-rate-distributed targets and
    ///   per-rate-type limits for one slave's master groups.
    /// @details Invoked once per activated slave by
    ///   `calculateMasterGroupConstraintsAndSendToSlaves()`.  Walks the
    ///   master groups attached to this slave and asks
    ///   `GroupConstraintCalculator` for the active target plus the four
    ///   non-active rate-type limits, returning them in the
    ///   `(injection_targets, production_constraints)` pair.  Injection
    ///   targets are always emitted in surface-rate units (cmode `RATE`)
    ///   because the slave cannot evaluate derived modes (REIN, VREP,
    ///   RESV) on its own.
    /// @param slave_idx Zero-based index of the activated slave.
    /// @param calculator Group-constraint calculator bound to the current
    ///   group/well state.  Stateful caches inside the calculator are
    ///   reused across slaves.
    /// @return `(injection_targets, production_constraints)` for this
    ///   slave's master groups, ready for Phase 2 / Phase 3.
    std::tuple<std::vector<InjectionGroupTarget>, std::vector<ProductionGroupConstraints>>
        calculateSlaveGroupConstraints_(std::size_t slave_idx, GroupConstraintCalculator<Scalar, IndexTraits>& calculator) const;

    /// @brief Compute the per-phase injection targets for one slave's
    ///   master groups.
    /// @details The injection half of `calculateSlaveGroupConstraints_()`.
    /// @param slave_idx Zero-based index of the activated slave.
    /// @param calculator Group-constraint calculator bound to the current
    ///   group/well state.
    /// @return One entry per (master group, phase) that has a target.
    std::vector<InjectionGroupTarget>
        calculateSlaveGroupInjectionTargets_(std::size_t slave_idx, GroupConstraintCalculator<Scalar, IndexTraits>& calculator) const;

    /// @brief Phase 2: exclude master groups none of whose slave producers
    ///   are under group control from the guide-rate distribution, and
    ///   redistribute the shortfall to sibling groups.
    /// @details Switches excluded groups to individual control so their
    ///   current rate becomes a reduction on the parent group, recomputes
    ///   the reductions, then re-evaluates all groups via
    ///   `GroupConstraintCalculator`; an excluded group gets the target it
    ///   would get on returning to group control.  Finally switches every
    ///   master group to individual control with its final allocated
    ///   target so the master completes its own time step assuming slave
    ///   rates remain constant.  Injection targets are handled separately,
    ///   see capAndRedistributeInjectionTargets_().
    /// @param calculator Group-constraint calculator (same instance as
    ///   used in Phase 1).
    /// @param all_production_constraints In/out: production constraints
    ///   from Phase 1, modified in place.
    void capAndRedistributeProductionTargets_(
        GroupConstraintCalculator<Scalar, IndexTraits>& calculator,
        std::vector<std::vector<ProductionGroupConstraints>>& all_production_constraints);

    /// @brief Phase 2b: exclude master groups none of whose slave injectors
    ///   are under group control for a phase from the guide-rate
    ///   distribution, and redistribute the surplus to sibling master groups.
    /// @details The injection counterpart of
    ///   capAndRedistributeProductionTargets_(), but decided from the number
    ///   of the slave group's injectors under group control (reported by the
    ///   slave) rather than from potentials, which for injectors with
    ///   cross-flow can be far from what the wells inject.  An excluded group
    ///   has its effective injection GCW for the phase set to 0, which drops
    ///   it from the parent's injection guide-rate sum and makes its rate a
    ///   parent target reduction instead.  The injection GCW and the
    ///   injection target reductions are then recomputed and all targets are
    ///   re-evaluated, so the remaining groups absorb the surplus.  The
    ///   excluded group gets the target it would get on returning to group
    ///   control, as a well that cannot meet its share does.
    /// @param calculator Group-constraint calculator (same instance as
    ///   used for the initial targets).
    /// @param all_injection_targets In/out: per-slave injection targets,
    ///   modified in place.
    void capAndRedistributeInjectionTargets_(
        GroupConstraintCalculator<Scalar, IndexTraits>& calculator,
        std::vector<std::vector<InjectionGroupTarget>>& all_injection_targets);

    /// @brief Pre-phase: switch master groups that route to currently-
    ///   inactive slaves to individual control so they are excluded from
    ///   guide-rate distribution.
    /// @details An inactive slave contributes zero rate and zero
    ///   potential and must not consume any share of the parent's
    ///   target.  Currently handles the not-yet-activated case; the
    ///   slave-finished-before-master case deferred to a follow-up PR.
    void excludeInactiveSlaveMasterGroupsFromDistribution_();

    /// @brief Exclude the master groups of currently-inactive slaves from
    ///   the injection guide-rate distribution of every phase.
    /// @details Sets their effective injection GCW to 0.  The injection
    ///   part of excludeInactiveSlaveMasterGroupsFromDistribution_().
    void excludeInactiveSlaveMasterGroupsFromInjectionDistribution_();

    /// @brief The constraints sent to a slave for one master group: the
    ///   active target and mode and the per-rate-type limits the calculator
    ///   found for it.
    static ProductionGroupConstraints makeProductionGroupConstraints_(
        std::size_t group_idx,
        const typename GroupConstraintCalculator<Scalar, IndexTraits>::ProductionConstraintResult& constraints);

    /// @brief Pre-phase: restore the production control mode the Schedule
    ///   gives each master group.
    /// @details The Phase 2 finalize of the previous calculation leaves every
    ///   master group on individual control for the master's own time step.
    ///   Left in place, a master group on FLD or NONE in the deck would be
    ///   treated as being on individual control, and fall back to its own,
    ///   undefined, GCONPROD target where its share of the parent's target is
    ///   not binding.
    void restoreMasterGroupProductionControls_();

    /// @brief Phase 3: send the computed constraints for one slave to
    ///   that slave over MPI.
    /// @details The underlying MPI sends in
    ///   `ReservoirCouplingMasterReportStep` are rank-0-only, so it is
    ///   safe to invoke this from every master rank inside the Phase-3
    ///   loop.
    /// @param rescoup_master Reservoir-coupling master facade owning the
    ///   slave communicators.
    /// @param slave_idx Zero-based index of the activated slave.
    /// @param injection_targets Per-phase injection rate targets to send.
    /// @param production_constraints Per-master-group production
    ///   constraint triples (target, cmode, per-rate-type limits) to
    ///   send.
    void sendSlaveGroupConstraintsToSlave_(
        const ReservoirCouplingMaster<Scalar>& rescoup_master,
        std::size_t slave_idx,
        const std::vector<InjectionGroupTarget>& injection_targets,
        const std::vector<ProductionGroupConstraints>& production_constraints
    ) const;

    /// @brief Recompute the Group-Controlled-Wells count and the
    ///   FIELD-level production target reduction.
    /// @details Called from the pre-phase and from Phase 2 (twice: after
    ///   the cap loop and after the finalize-to-individual loop).  The
    ///   ordering matters: `updateGroupControlledWells()` must run before
    ///   `updateGroupTargetReduction()` because the reduction depends on
    ///   GCW.  The reduction rates are then summed across the master's
    ///   ranks (they are per-rank partial sums), so that every rank
    ///   computes the same targets from them.  MPI-collective on the master
    ///   communicator.
    void updateGCWAndTargetReductions_();

    /// @brief Recompute the injection Group-Controlled-Wells count for
    ///   every phase and the FIELD-level injection target reduction.
    /// @details The injection counterpart of
    ///   updateGCWAndTargetReductions_(), needed after effective injection
    ///   GCW entries change.  The injection reduction rates are then summed
    ///   across the master's ranks, as for production.  MPI-collective on the
    ///   master communicator.
    void updateInjectionGCWAndTargetReductions_();

    GuideRateHandler<Scalar, IndexTraits>& guide_rate_handler_;
    GroupStateHelper<Scalar, IndexTraits>& group_state_helper_;
    const WellState<Scalar, IndexTraits>& well_state_;
    const GroupState<Scalar>& group_state_;
    const int report_step_idx_;
    const Schedule& schedule_;
    const SummaryState& summary_state_;
    DeferredLogger& deferred_logger_;
    ReservoirCouplingMaster<Scalar>& reservoir_coupling_master_;
    BlackoilWellModelGeneric<Scalar, IndexTraits>& well_model_;
    const PhaseUsageInfo<IndexTraits>& phase_usage_;
};

}  // namespace Opm
#endif // OPM_RESCOUP_CONSTRAINTS_CALCULATOR_HPP
