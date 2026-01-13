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
// Constructor for the RescoupTargetCalculator class
// -------------------------------------------------------
template <class Scalar, class IndexTraits>
RescoupTargetCalculator<Scalar, IndexTraits>::
RescoupTargetCalculator(
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

// Calculates the target of each master group.
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
template <class Scalar, class IndexTraits>
void
RescoupTargetCalculator<Scalar, IndexTraits>::
calculateMasterGroupTargetsAndSendToSlaves()
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
                auto [injection_targets, production_targets] =
                    this->calculateSlaveGroupTargets_(slave_idx, calculator);
                this->sendSlaveGroupTargetsToSlave_(
                    rescoup_master, slave_idx, injection_targets, production_targets
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
  std::vector<typename RescoupTargetCalculator<Scalar, IndexTraits>::InjectionGroupTarget>,
  std::vector<typename RescoupTargetCalculator<Scalar, IndexTraits>::ProductionGroupTarget>
>
RescoupTargetCalculator<Scalar, IndexTraits>::
calculateSlaveGroupTargets_(std::size_t slave_idx, GroupTargetCalculator<Scalar, IndexTraits>& calculator) const
{
    std::vector<InjectionGroupTarget> injection_targets;
    std::vector<ProductionGroupTarget> production_targets;
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
            auto target_info = calculator.groupProductionTarget(group);
            if (target_info.has_value()) {
                production_targets.push_back(
                    ProductionGroupTarget{
                        group_idx, target_info->target, target_info->cmode
                    }
                );
            }
        }
    }
    return {injection_targets, production_targets};
}

template <class Scalar, class IndexTraits>
void
RescoupTargetCalculator<Scalar, IndexTraits>::
sendSlaveGroupTargetsToSlave_(
    const ReservoirCouplingMaster<Scalar>& rescoup_master,
    std::size_t slave_idx,
    const std::vector<InjectionGroupTarget>& injection_targets,
    const std::vector<ProductionGroupTarget>& production_targets
) const
{
    auto num_injection_targets = injection_targets.size();
    auto num_production_targets = production_targets.size();
    // First, send the number of targets such that the slave can know if it can expect none
    // or more targets.
    rescoup_master.sendNumGroupTargetsToSlave(slave_idx, num_injection_targets, num_production_targets);
    if (num_injection_targets > 0) {
        rescoup_master.sendInjectionTargetsToSlave(slave_idx, injection_targets);
    }
    if (num_production_targets > 0) {
        rescoup_master.sendProductionTargetsToSlave(slave_idx, production_targets);
    }
}

template class RescoupTargetCalculator<double, BlackOilDefaultFluidSystemIndices>;

#if FLOW_INSTANTIATE_FLOAT
template class RescoupTargetCalculator<float, BlackOilDefaultFluidSystemIndices>;
#endif

}// namespace Opm
