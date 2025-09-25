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

#ifndef OPM_RESCOUP_RECEIVE_SLAVE_GROUP_DATA_HPP
#define OPM_RESCOUP_RECEIVE_SLAVE_GROUP_DATA_HPP

#include <opm/simulators/flow/rescoup/ReservoirCoupling.hpp>
#include <opm/simulators/wells/GuideRateHandler.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include <opm/simulators/wells/WellGroupHelper.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>

namespace Opm {

template<class Scalar, class IndexTraits>
class RescoupReceiveSlaveGroupData {
public:
    using SlaveGroupProductionData  = ReservoirCoupling::SlaveGroupProductionData<Scalar>;
    using SlaveGroupInjectionData = ReservoirCoupling::SlaveGroupInjectionData<Scalar>;
    using Potentials = ReservoirCoupling::Potentials<Scalar>;
    using ProductionRates = ReservoirCoupling::ProductionRates<Scalar>;
    using WellGroupHelperType = WellGroupHelper<Scalar, IndexTraits>;

    RescoupReceiveSlaveGroupData(WellGroupHelperType& wg_helper);

    void receiveSlaveGroupData();
private:

    WellGroupHelperType& wg_helper_;
    ReservoirCouplingMaster<Scalar>& reservoir_coupling_master_;
    const Schedule& schedule_;
    const GroupState<Scalar>& group_state_;
    const PhaseUsageInfo<IndexTraits>& phase_usage_;
    const int report_step_idx_;
};

} // namespace Opm
#endif // OPM_RESCOUP_RECEIVE_SLAVE_GROUP_DATA_HPP
