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
#include <opm/simulators/flow/rescoup/ReservoirCoupling.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingMpiTraits.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingSlave.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingSlaveReportStep.hpp>

#include <opm/input/eclipse/Schedule/ResCoup/ReservoirCouplingInfo.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/MasterGroup.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/Slaves.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <dune/common/parallel/mpitraits.hh>

#include <vector>
#include <fmt/format.h>

namespace Opm {

template <class Scalar>
ReservoirCouplingSlaveReportStep<Scalar>::
ReservoirCouplingSlaveReportStep(
    ReservoirCouplingSlave<Scalar> &slave
) :
    slave_{slave}
{
}

// ------------------
// Public methods
// ------------------

template <class Scalar>
void
ReservoirCouplingSlaveReportStep<Scalar>::
sendInjectionDataToMaster(
    const std::vector<ReservoirCoupling::SlaveGroupInjectionData<Scalar>> &injection_data
) const
{
    sendDataToMaster_(injection_data, MessageTag::SlaveInjectionData, "injection data");
}

template <class Scalar>
void
ReservoirCouplingSlaveReportStep<Scalar>::
sendProductionDataToMaster(
    const std::vector<ReservoirCoupling::SlaveGroupProductionData<Scalar>> &production_data
) const
{
    sendDataToMaster_(production_data, MessageTag::SlaveProductionData, "production data");
}

// ------------------
// Private methods
// ------------------

template <class Scalar>
template <class DataType>
void
ReservoirCouplingSlaveReportStep<Scalar>::
sendDataToMaster_(
    const std::vector<DataType>& data,
    MessageTag tag,
    const std::string& data_type_name
) const
{
    // NOTE: The master can determine from the ordering of the data in the vector
    //   which slave group for a given slave name the given data belong to,
    //   so we do not need to send the slave group names also.
    if (this->comm().rank() == 0) {
        auto num_groups = data.size();
        auto MPI_DATA_TYPE = Dune::MPITraits<DataType>::getType();
        MPI_Send(
            data.data(),
            /*count=*/num_groups,
            /*datatype=*/MPI_DATA_TYPE,
            /*dest_rank=*/0,
            /*tag=*/static_cast<int>(tag),
            this->getSlaveMasterComm()
        );
        this->logger().info(
            fmt::format("Sent {} for {} groups to master process from rank 0, slave name: {}",
                data_type_name, num_groups, this->slaveName())
        );
    }
}

// Explicit instantiations
template class ReservoirCouplingSlaveReportStep<double>;
#if FLOW_INSTANTIATE_FLOAT
template class ReservoirCouplingSlaveReportStep<float>;
#endif

} // namespace Opm
