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

    /// @brief Phase 2: cap each capacity-limited master group's
    ///   production target at the slave's reported potential and
    ///   redistribute the surplus to sibling groups.
    /// @details Switches capped groups to individual control so their
    ///   capped rate becomes a reduction on the parent group, recomputes
    ///   the FIELD-level reduction, then re-evaluates the uncapped
    ///   groups via `GroupConstraintCalculator`.  Finally switches every
    ///   master group to individual control with its final allocated
    ///   target so the master completes its own time step assuming slave
    ///   rates remain constant.  Injection targets are not capped here;
    ///   deferred to future PR.
    /// @param calculator Group-constraint calculator (same instance as
    ///   used in Phase 1).
    /// @param all_production_constraints In/out: production constraints
    ///   from Phase 1, modified in place.
    void capAndRedistributeProductionTargets_(
        GroupConstraintCalculator<Scalar, IndexTraits>& calculator,
        std::vector<std::vector<ProductionGroupConstraints>>& all_production_constraints);

    /// @brief Pre-phase: switch master groups that route to currently-
    ///   inactive slaves to individual control so they are excluded from
    ///   guide-rate distribution.
    /// @details An inactive slave contributes zero rate and zero
    ///   potential and must not consume any share of the parent's
    ///   target.  Currently handles the not-yet-activated case; the
    ///   slave-finished-before-master case deferred to a follow-up PR.
    void excludeInactiveSlaveMasterGroupsFromDistribution_();

    /// @brief Project a slave-reported potentials triple onto a single
    ///   value comparable to a target in the given control mode.
    /// @param potentials Slave-reported potentials (oil, gas, water).
    /// @param cmode Active production control mode the cap test is
    ///   expressed in (`ORAT`, `WRAT`, `GRAT`, `LRAT`, ...).
    /// @return Potential value in SI units suitable for direct
    ///   comparison against the target, or a negative sentinel when the
    ///   cmode does not map to a single phase / linear combination
    ///   (e.g. `RESV`, `FLD`, `NONE`).  In that case Phase 2 skips the
    ///   cap test for the group.
    Scalar potentialForProductionCmode_(
        const ReservoirCoupling::Potentials<Scalar>& potentials,
        Group::ProductionCMode cmode) const;

    /// @brief Pre-phase: restore each master group's `production_control`
    ///   cmode to its Schedule-defined GCONPROD value.
    /// @details The previous sync step's Phase 2 finalize switched these
    ///   to individual control; this restore lets the upcoming
    ///   distribution see the user-specified cmode.  Reads from the
    ///   Schedule directly (rather than caching the pre-sync state) so
    ///   new GCONPROD records at report-step boundaries are picked up
    ///   automatically.
    void restoreMasterGroupControlsFromSchedule_();

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
    ///   GCW.  Both delegated calls are MPI-collective on the master
    ///   communicator.
    void updateGCWAndTargetReductions_();

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
