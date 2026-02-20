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
#include <opm/simulators/wells/GroupTargetCalculator.hpp>
#include <opm/simulators/wells/GuideRateHandler.hpp>
#include <opm/simulators/wells/WellState.hpp>

namespace Opm {

template<class Scalar, class IndexTraits>
class RescoupConstraintsCalculator {
public:
    using InjectionGroupTarget = ReservoirCoupling::InjectionGroupTarget<Scalar>;
    using ProductionGroupConstraints = ReservoirCoupling::ProductionGroupConstraints<Scalar>;
    RescoupConstraintsCalculator(
        GuideRateHandler<Scalar, IndexTraits>& guide_rate_handler,
        GroupStateHelper<Scalar, IndexTraits>& group_state_helper
    );

    void calculateMasterGroupConstraintsAndSendToSlaves();
private:
    std::tuple<std::vector<InjectionGroupTarget>, std::vector<ProductionGroupConstraints>>
        calculateSlaveGroupConstraints_(std::size_t slave_idx, GroupTargetCalculator<Scalar, IndexTraits>& calculator) const;
    void sendSlaveGroupConstraintsToSlave_(
        const ReservoirCouplingMaster<Scalar>& rescoup_master,
        std::size_t slave_idx,
        const std::vector<InjectionGroupTarget>& injection_targets,
        const std::vector<ProductionGroupConstraints>& production_constraints
    ) const;

    GuideRateHandler<Scalar, IndexTraits>& guide_rate_handler_;
    const GroupStateHelper<Scalar, IndexTraits>& group_state_helper_;
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
