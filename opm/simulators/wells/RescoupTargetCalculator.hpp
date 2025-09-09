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

#ifndef OPM_RESCOUP_TARGET_CALCULATOR_HPP
#define OPM_RESCOUP_TARGET_CALCULATOR_HPP
#include <opm/simulators/flow/ReservoirCoupling.hpp>
#include <opm/simulators/flow/ReservoirCouplingMaster.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRate.hpp>
#include <opm/simulators/utils/BlackoilPhases.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/BlackoilWellModelGeneric.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include <opm/simulators/wells/GroupTargetCalculator.hpp>
#include <opm/simulators/wells/GuideRateHandler.hpp>
#include <opm/simulators/wells/WellState.hpp>

namespace Opm {

template<class Scalar>
class RescoupTargetCalculator {
public:
    using InjectionGroupTarget = ReservoirCoupling::InjectionGroupTarget<Scalar>;
    using ProductionGroupTarget = ReservoirCoupling::ProductionGroupTarget<Scalar>;
    RescoupTargetCalculator(
        GuideRateHandler<Scalar>& guide_rate_handler,
        const WellState<Scalar>& well_state,
        const GroupState<Scalar>& group_state,
        const int report_step_idx
    );

    void calculateMasterGroupTargetsAndSendToSlaves();
private:
    std::tuple<std::vector<InjectionGroupTarget>, std::vector<ProductionGroupTarget>>
        calculateSlaveGroupTargets_(std::size_t slave_idx, GroupTargetCalculator<Scalar>& calculator) const;
    void sendSlaveGroupTargetsToSlave_(
        const ReservoirCouplingMaster<Scalar>& rescoup_master,
        std::size_t slave_idx,
        const std::vector<InjectionGroupTarget>& injection_targets,
        const std::vector<ProductionGroupTarget>& production_targets
    ) const;

    GuideRateHandler<Scalar>& guide_rate_handler_;
    const WellState<Scalar>& well_state_;
    const GroupState<Scalar>& group_state_;
    const int report_step_idx_;
    const Schedule& schedule_;
    const SummaryState& summary_state_;
    DeferredLogger& deferred_logger_;
    ReservoirCouplingMaster<Scalar>& reservoir_coupling_master_;
    BlackoilWellModelGeneric<Scalar>& well_model_;
    const PhaseUsage& phase_usage_;
};

}  // namespace Opm
#endif // OPM_RESCOUP_TARGET_CALCULATOR_HPP
