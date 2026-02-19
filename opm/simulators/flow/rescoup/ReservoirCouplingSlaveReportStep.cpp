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

template <class Scalar>
std::pair<std::size_t, std::size_t>
ReservoirCouplingSlaveReportStep<Scalar>::
receiveNumGroupConstraintsFromMaster() const {
    std::vector<std::size_t> num_group_targets(2);
    if (this->comm().rank() == 0) {
        auto MPI_SIZE_T_TYPE = Dune::MPITraits<std::size_t>::getType();
        // NOTE: See comment about error handling at the top of this file.
        MPI_Recv(
            num_group_targets.data(),
            /*count=*/2,
            /*datatype=*/MPI_SIZE_T_TYPE,
            /*source_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::NumSlaveGroupConstraints),
            this->getSlaveMasterComm(),
            MPI_STATUS_IGNORE
        );
        this->logger().info("Received number of slave group constraints from master process rank 0");
    }
    this->comm().broadcast(num_group_targets.data(), /*count=*/2, /*emitter_rank=*/0);
    auto num_injection_targets = num_group_targets[0];
    auto num_production_constraints = num_group_targets[1];
    this->logger().info(fmt::format("Received number of injection targets: {} and "
                             "production constraints: {} from master process",
                             num_injection_targets, num_production_constraints));
    return std::make_pair(num_injection_targets, num_production_constraints);
}

template <class Scalar>
void
ReservoirCouplingSlaveReportStep<Scalar>::
receiveInjectionGroupTargetsFromMaster(std::size_t num_targets) const
{
    std::vector<InjectionGroupTarget> injection_targets(num_targets);
    if (this->comm().rank() == 0) {
        auto MPI_INJECTION_GROUP_TARGET_TYPE = Dune::MPITraits<InjectionGroupTarget>::getType();
        // NOTE: See comment about error handling at the top of this file.
        MPI_Recv(
            injection_targets.data(),
            /*count=*/num_targets,
            /*datatype=*/MPI_INJECTION_GROUP_TARGET_TYPE,
            /*source_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::InjectionGroupTargets),
            this->getSlaveMasterComm(),
            MPI_STATUS_IGNORE
        );
        this->logger().info(fmt::format(
            "Received injection {} group targets from master process rank 0", num_targets
        ));
    }
    this->comm().broadcast(
        injection_targets.data(), num_targets, /*emitter_rank=*/0
    );
}

template <class Scalar>
void
ReservoirCouplingSlaveReportStep<Scalar>::
receiveProductionGroupConstraintsFromMaster(std::size_t num_targets) const
{
    std::vector<ProductionGroupConstraints> production_constraints(num_targets);
    if (this->comm().rank() == 0) {
        auto MPI_PRODUCTION_GROUP_CONSTRAINTS_TYPE = Dune::MPITraits<ProductionGroupConstraints>::getType();
        // NOTE: See comment about error handling at the top of this file.
        MPI_Recv(
            production_constraints.data(),
            /*count=*/num_targets,
            /*datatype=*/MPI_PRODUCTION_GROUP_CONSTRAINTS_TYPE,
            /*source_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::ProductionGroupConstraints),
            this->getSlaveMasterComm(),
            MPI_STATUS_IGNORE
        );
        this->logger().info(fmt::format(
            "Received production {} group constraints from master process rank 0", num_targets
        ));
    }
    this->comm().broadcast(
        production_constraints.data(), num_targets, /*emitter_rank=*/0
    );
}

// Explicit instantiations
template class ReservoirCouplingSlaveReportStep<double>;
#if FLOW_INSTANTIATE_FLOAT
template class ReservoirCouplingSlaveReportStep<float>;
#endif

} // namespace Opm
