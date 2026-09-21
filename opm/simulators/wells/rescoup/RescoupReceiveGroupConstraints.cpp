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
#include <opm/simulators/wells/rescoup/RescoupReceiveGroupConstraints.hpp>

#include <fmt/format.h>

namespace Opm {

template <class Scalar, class IndexTraits>
RescoupReceiveGroupConstraints<Scalar, IndexTraits>::
RescoupReceiveGroupConstraints(
    GuideRateHandler<Scalar, IndexTraits>& guide_rate_handler,
    GroupStateHelper<Scalar, IndexTraits>& group_state_helper
)
    : guide_rate_handler_{guide_rate_handler}
    , group_state_helper_{group_state_helper}
    , reservoir_coupling_slave_{guide_rate_handler.reservoirCouplingSlave()}
{
}

template <class Scalar, class IndexTraits>
void
RescoupReceiveGroupConstraints<Scalar, IndexTraits>::
receiveGroupConstraintsFromMaster(const ReservoirCoupling::GroupConstraintsSend send)
{
    // NOTE: All ranks must call these functions because they contain broadcasts.
    //   The MPI_Recv parts inside the functions have their own rank 0 checks.
    auto& rescoup_slave = this->reservoir_coupling_slave_;
    auto [num_inj_targets, num_prod_constraints] = rescoup_slave.receiveNumGroupConstraintsFromMaster();
    if (num_inj_targets > 0) {
        rescoup_slave.receiveInjectionGroupTargetsFromMaster(num_inj_targets);
    }
    if (num_prod_constraints > 0) {
        rescoup_slave.receiveProductionGroupConstraintsFromMaster(num_prod_constraints);
    }
    else if (send == ReservoirCoupling::GroupConstraintsSend::Handshake) {
        // The handshake's production list is complete, so none means none:
        // the master has withdrawn every production constraint it imposed on
        // this slave, and the ones held from an earlier step must go.  A
        // refresh sends an empty production list to mean "unchanged", so
        // there the held ones stay.
        rescoup_slave.clearMasterProductionConstraints();
    }
}

template class RescoupReceiveGroupConstraints<double, BlackOilDefaultFluidSystemIndices>;

#if FLOW_INSTANTIATE_FLOAT
template class RescoupReceiveGroupConstraints<float, BlackOilDefaultFluidSystemIndices>;
#endif

} // namespace Opm
