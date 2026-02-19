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
#include <opm/simulators/wells/rescoup/RescoupTargetCalculator.hpp>

#include <array>
#include <string>
#include <vector>
#include <tuple>

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
// Details on the target calculation:
// ----------------------------------
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
// Details on the per-rate-type limit calculation:
// -----------------------------------------------
// For each non-active rate type (ORAT, WRAT, GRAT, LRAT, RESV) that is not the active
// cmode, the same hierarchy-traversal logic described above is reused with an "explicit
// cmode" parameter.  The difference is in the recursion stopping criterion:
//  - Active target: stops at an ancestor whose production_control() != FLD/NONE.
//  - Per-rate-type limit: stops at an ancestor that has_control(rate_type), i.e. that
//    defines a GCONPROD limit for that specific rate type.
// The guide-rate fraction calculation (FractionCalculator) is identical in both cases,
// since fractions represent proportional capacity allocation independent of rate type.
// If no ancestor defines a limit for a given rate type, the limit is set to -1 (undefined).
// See GroupTargetCalculator::groupProductionConstraints() for the implementation.
template <class Scalar, class IndexTraits>
void
RescoupConstraintsCalculator<Scalar, IndexTraits>::
calculateMasterGroupConstraintsAndSendToSlaves()
{
    // NOTE: Since this object can only be constructed for a master process, we can be
    //   sure that if we are here, we are running as master
    auto& rescoup_master = this->reservoir_coupling_master_;
    const auto& comm = rescoup_master.getComm();
    if (comm.rank() == 0) {
        GroupTargetCalculator calculator{
            this->well_model_,
            this->group_state_helper_
        };
        auto num_slaves = rescoup_master.numSlaves();
        for (std::size_t slave_idx = 0; slave_idx < num_slaves; ++slave_idx) {
            if (rescoup_master.slaveIsActivated(slave_idx)) {
                auto [injection_targets, production_constraints] =
                    this->calculateSlaveGroupConstraints_(slave_idx, calculator);
                this->sendSlaveGroupConstraintsToSlave_(
                    rescoup_master, slave_idx, injection_targets, production_constraints
                );
            }
        }
    }
}

// NOTE on reuse and future refactor:
// This method relies on GroupTargetCalculator to compute group-level targets
// for reservoir coupling. Similar target logic exists in GroupStateHelper
// (checkGroupContraintsProd/getWellGroupTargetProducer and
// checkGroupContraintsInj/getWellGroupTargetInjector). The plan is to make
// GroupTargetCalculator the general implementation and refactor the existing
// helpers to reuse it, removing duplication. We will address that in a follow-up
// PR to keep this PR scoped to reservoir coupling behavior.
template <class Scalar, class IndexTraits>
std::tuple<
  std::vector<typename RescoupConstraintsCalculator<Scalar, IndexTraits>::InjectionGroupTarget>,
  std::vector<typename RescoupConstraintsCalculator<Scalar, IndexTraits>::ProductionGroupConstraints>
>
RescoupConstraintsCalculator<Scalar, IndexTraits>::
calculateSlaveGroupConstraints_(std::size_t slave_idx, GroupTargetCalculator<Scalar, IndexTraits>& calculator) const
{
    std::vector<InjectionGroupTarget> injection_targets;
    std::vector<ProductionGroupConstraints> production_constraints;
    auto& rescoup_master = this->reservoir_coupling_master_;
    static const std::array<ReservoirCoupling::Phase, 3> phases = {
        ReservoirCoupling::Phase::Water, ReservoirCoupling::Phase::Oil, ReservoirCoupling::Phase::Gas
    };
    const auto& master_groups = rescoup_master.getMasterGroupNamesForSlave(slave_idx);
    for (std::size_t group_idx = 0; group_idx < master_groups.size(); ++group_idx) {
        const auto& group_name = master_groups[group_idx];
        const Group& group = this->schedule_.getGroup(group_name, this->report_step_idx_);
        if (group.isInjectionGroup()) {
            for (ReservoirCoupling::Phase phase : phases) {
                auto target_info = calculator.groupInjectionTarget(group, phase);
                if (target_info.has_value()) {
                    injection_targets.push_back(
                        InjectionGroupTarget{
                            group_idx, target_info->target, target_info->cmode, phase
                        }
                    );
                }
            }
        }
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

template class RescoupConstraintsCalculator<double, BlackOilDefaultFluidSystemIndices>;

#if FLOW_INSTANTIATE_FLOAT
template class RescoupConstraintsCalculator<float, BlackOilDefaultFluidSystemIndices>;
#endif

}// namespace Opm
