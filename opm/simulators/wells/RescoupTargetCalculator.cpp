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
#include <opm/simulators/wells/GroupTargetCalculator.hpp>
#include <opm/simulators/wells/RescoupTargetCalculator.hpp>

#include <array>
#include <string>
#include <vector>
#include <tuple>

#include <fmt/format.h>

namespace Opm {

// -------------------------------------------------------
// Constructor for the RescoupTargetCalculator class
// -------------------------------------------------------
template <class Scalar>
RescoupTargetCalculator<Scalar>::
RescoupTargetCalculator(
    GuideRateHandler<Scalar>& guide_rate_handler,
    const WellState<Scalar>& well_state,
    const GroupState<Scalar>& group_state,
    const int report_step_idx
)
    : guide_rate_handler_{guide_rate_handler}
    , well_state_{well_state}
    , group_state_{group_state}
    , report_step_idx_{report_step_idx}
    , schedule_{guide_rate_handler.schedule()}
    , summary_state_{guide_rate_handler.summaryState()}
    , deferred_logger_{guide_rate_handler.deferredLogger()}
    , reservoir_coupling_master_{guide_rate_handler.reservoirCouplingMaster()}
    , well_model_{guide_rate_handler.wellModel()}
    , phase_usage_{guide_rate_handler.phaseUsage()}
{
}

// The master group targets are calculated based on their guide rates and a
//  fixed target of a parent group that is distributed down to the master groups.
// This means that the master groups must have their guide rates defined, i.e. either
//      (a)  item 9 of GCONPROD is set to a positive guide rate value, or
//      (b)  item 9 is defaulted and item 10 is set to FORM.
//   In case (b), the the formula defined in GUIDERAT will use the communicated slave group
//   potentials to calculate the master group guide rates
template <class Scalar>
void
RescoupTargetCalculator<Scalar>::
calculateMasterGroupTargetsAndSendToSlaves()
{
    // NOTE: Since this object can only be constructed for a master process, we can be
    //   sure that if we are here, we are running as master
    auto& rescoup_master = this->reservoir_coupling_master_;
    const auto& comm = rescoup_master.getComm();
    if (comm.rank() == 0) {
        GroupTargetCalculator calculator{
            this->well_model_,
            this->well_state_,
            this->group_state_,
            this->schedule_,
            this->summary_state_,
            this->guide_rate_handler_.phaseUsage(),
            this->guide_rate_handler_.guideRate(),
            this->report_step_idx_,
            this->deferred_logger_
        };
        auto num_slaves = rescoup_master.numSlaves();
        static const std::array<Phase, 3> phases = { Phase::WATER, Phase::OIL, Phase::GAS };
        for (std::size_t slave_idx = 0; slave_idx < num_slaves; ++slave_idx) {
            std::vector<InjectionGroupTarget> injection_targets;
            std::vector<ProductionGroupTarget> production_targets;
            const auto& master_groups = rescoup_master.getMasterGroupNamesForSlave(slave_idx);
            for (std::size_t group_idx = 0; group_idx < master_groups.size(); ++group_idx) {
                const auto& group_name = master_groups[group_idx];
                const Group& group = this->schedule_[this->report_step_idx_].groups.get(group_name);
                if (group.isInjectionGroup()) {
                    for (Phase phase : phases) {
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
            this->sendSlaveGroupTargetsToSlave_(
                rescoup_master, slave_idx, injection_targets, production_targets
            );
        }
    }
}

template <class Scalar>
void
RescoupTargetCalculator<Scalar>::
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

template class RescoupTargetCalculator<double>;

#if FLOW_INSTANTIATE_FLOAT
template class RescoupTargetCalculator<float>;
#endif

}// namespace Opm
