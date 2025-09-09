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
#include <opm/simulators/wells/RescoupReceiveGroupTargets.hpp>

#include <fmt/format.h>

namespace Opm {

template <class Scalar, class IndexTraits>
RescoupReceiveGroupTargets<Scalar, IndexTraits>::
RescoupReceiveGroupTargets(
    GuideRateHandler<Scalar, IndexTraits>& guide_rate_handler,
    const WellState<Scalar, IndexTraits>& well_state,
    const GroupState<Scalar>& group_state,
    const int report_step_idx
)
    : guide_rate_handler_{guide_rate_handler}
    , well_state_{well_state}
    , group_state_{group_state}
    , report_step_idx_{report_step_idx}
    , reservoir_coupling_slave_{guide_rate_handler.reservoirCouplingSlave()}
{
}

template <class Scalar, class IndexTraits>
void
RescoupReceiveGroupTargets<Scalar, IndexTraits>::
receiveGroupTargetsFromMaster()
{
    // NOTE: Since this object can only be constructed for a slave process, we can be
    //   sure that if we are here, we are running as a slave
    auto& rescoup_slave = this->reservoir_coupling_slave_;
    const auto& comm = rescoup_slave.getComm();
    if (comm.rank() == 0) {
        auto [num_inj_targets, num_prod_targets] = rescoup_slave.receiveNumGroupTargetsFromMaster();
        if (num_inj_targets > 0) {
            rescoup_slave.receiveInjectionGroupTargetsFromMaster(num_inj_targets);
        }
        if (num_prod_targets > 0) {
            rescoup_slave.receiveProductionGroupTargetsFromMaster(num_prod_targets);
        }
    }
}

template class RescoupReceiveGroupTargets<double, BlackOilDefaultFluidSystemIndices>;

#if FLOW_INSTANTIATE_FLOAT
template class RescoupReceiveGroupTargets<float, BlackOilDefaultFluidSystemIndices>;
#endif

} // namespace Opm
