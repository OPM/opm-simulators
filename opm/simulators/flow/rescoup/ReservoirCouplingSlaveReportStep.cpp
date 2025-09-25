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
    // NOTE: The master can determine from the ordering of the injection data in the vector
    //   which slave group for a given slave name the given injection data belong to,
    //   so we do not need to send the slave group names also.
    if (this->comm().rank() == 0) {
        auto num_groups = injection_data.size();
        auto MPI_INJECTION_DATA_TYPE = Dune::MPITraits<SlaveGroupInjectionData>::getType();
        MPI_Send(
            injection_data.data(),
            /*count=*/num_groups,
            /*datatype=*/MPI_INJECTION_DATA_TYPE,
            /*dest_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::SlaveInjectionData),
            this->getSlaveMasterComm()
        );
        this->logger().info(
            fmt::format("Sent injection data for {} groups to master process from rank 0, slave name: {}",
                num_groups, this->slaveName())
        );
    }
}

template <class Scalar>
void
ReservoirCouplingSlaveReportStep<Scalar>::
sendProductionDataToMaster(
    const std::vector<ReservoirCoupling::SlaveGroupProductionData<Scalar>> &production_data
) const
{
    // NOTE: The master can determine from the ordering of the production data in the vector
    //   which slave group for a given slave name the given production data belong to,
    //   so we do not need to send the slave group names also.
    if (this->comm().rank() == 0) {
        auto num_groups = production_data.size();
        auto MPI_PRODUCTION_DATA_TYPE = Dune::MPITraits<SlaveGroupProductionData>::getType();
        MPI_Send(
            production_data.data(),
            /*count=*/num_groups,
            /*datatype=*/MPI_PRODUCTION_DATA_TYPE,
            /*dest_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::SlaveProductionData),
            this->getSlaveMasterComm()
        );
        this->logger().info(
            fmt::format("Sent production data for {} groups to master process from rank 0, slave name: {}",
                num_groups, this->slaveName())
        );
    }
}

// Explicit instantiations
template class ReservoirCouplingSlaveReportStep<double>;
#if FLOW_INSTANTIATE_FLOAT
template class ReservoirCouplingSlaveReportStep<float>;
#endif

} // namespace Opm
