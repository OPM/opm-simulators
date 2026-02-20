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
#include <opm/simulators/wells/rescoup/RescoupReceiveSlaveGroupData.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingMaster.hpp>

#include <array>
#include <string>
#include <vector>
#include <tuple>

#include <fmt/format.h>

namespace Opm {

// -------------------------------------------------------
// Constructor for the RescoupReceiveSlaveGroupData class
// -------------------------------------------------------
template <class Scalar, class IndexTraits>
RescoupReceiveSlaveGroupData<Scalar, IndexTraits>::
RescoupReceiveSlaveGroupData(
    GroupStateHelperType& groupStateHelper
)
    : groupStateHelper_{groupStateHelper}
    , reservoir_coupling_master_{groupStateHelper.reservoirCouplingMaster()}
    , schedule_{groupStateHelper.schedule()}
    , group_state_{groupStateHelper.groupState()}
    , phase_usage_{groupStateHelper.phaseUsage()}
    , report_step_idx_{groupStateHelper.reportStepIdx()}
{
}

template <class Scalar, class IndexTraits>
void
RescoupReceiveSlaveGroupData<Scalar, IndexTraits>::
receiveSlaveGroupData()
{
    auto& rescoup_master = this->reservoir_coupling_master_;
    rescoup_master.receiveProductionDataFromSlaves();
    rescoup_master.receiveInjectionDataFromSlaves();
}


template class RescoupReceiveSlaveGroupData<double, BlackOilDefaultFluidSystemIndices>;

#if FLOW_INSTANTIATE_FLOAT
template class RescoupReceiveSlaveGroupData<float, BlackOilDefaultFluidSystemIndices>;
#endif

} // namespace Opm
