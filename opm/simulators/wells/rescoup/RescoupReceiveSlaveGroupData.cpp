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
// Constructor for the RescoupTargetCalculator class
// -------------------------------------------------------
template <class Scalar, class IndexTraits>
RescoupReceiveSlaveGroupData<Scalar, IndexTraits>::
RescoupReceiveSlaveGroupData(
    WellGroupHelperType& wg_helper
)
    : wg_helper_{wg_helper}
    , reservoir_coupling_master_{wg_helper.reservoirCouplingMaster()}
    , schedule_{wg_helper.schedule()}
    , group_state_{wg_helper.groupState()}
    , phase_usage_{wg_helper.phaseUsage()}
    , report_step_idx_{wg_helper.reportStepIdx()}
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
