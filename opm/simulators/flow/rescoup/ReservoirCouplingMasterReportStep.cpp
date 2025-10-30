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
#include <opm/simulators/flow/rescoup/ReservoirCouplingMaster.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingMasterReportStep.hpp>

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
ReservoirCouplingMasterReportStep<Scalar>::
ReservoirCouplingMasterReportStep(
    ReservoirCouplingMaster<Scalar> &master
) :
    master_{master}
{
}

// ------------------
// Public methods
// ------------------

template <class Scalar>
const ReservoirCoupling::Potentials<Scalar>&
ReservoirCouplingMasterReportStep<Scalar>::
getSlaveGroupPotentials(const std::string &master_group_name) const
{
    auto it = this->getMasterGroupToSlaveNameMap().find(master_group_name);
    if (it != this->getMasterGroupToSlaveNameMap().end()) {
        auto& slave_name = it->second;
        auto group_idx = this->getMasterGroupCanonicalIdx(slave_name, master_group_name);
        return this->slave_group_production_data_.at(slave_name)[group_idx].potentials;
    }
    else {
        OPM_THROW(
            std::runtime_error,
            fmt::format(
                "Master group name {} not found in master-to-slave-group-name mapping",
                master_group_name
            )
        );
    }
}

template <class Scalar>
void
ReservoirCouplingMasterReportStep<Scalar>::
receiveInjectionDataFromSlaves()
{
    auto num_slaves = this->numSlaves();
    this->logger().info("Receiving injection data from slave processes");
    for (unsigned int i = 0; i < num_slaves; i++) {
        auto num_slave_groups = this->numSlaveGroups(i);
        assert( num_slave_groups > 0 );
        std::vector<SlaveGroupInjectionData> injection_data(num_slave_groups);
        if (this->comm().rank() == 0) {
            if (this->slaveIsActivated(i)) {
                // NOTE: See comment about error handling at the top of this file.
                auto MPI_INJECTION_DATA_TYPE = Dune::MPITraits<SlaveGroupInjectionData>::getType();
                MPI_Recv(
                    injection_data.data(),
                    /*count=*/num_slave_groups,
                    /*datatype=*/MPI_INJECTION_DATA_TYPE,
                    /*source_rank=*/0,
                    /*tag=*/static_cast<int>(MessageTag::SlaveInjectionData),
                    this->getSlaveComm(i),
                    MPI_STATUS_IGNORE
                );
                this->logger().info(
                    fmt::format(
                        "Received injection data for {} groups from slave process with name: {}. "
                        "Number of slave groups: {}", num_slave_groups, this->slaveName(i), num_slave_groups
                    )
                );
            }
            else {
                this->logger().info(fmt::format(
                    "Slave {} has not activated yet, skipping receiving injection data from slave",
                        this->slaveName(i)
                    )
                );
                injection_data.assign(num_slave_groups, SlaveGroupInjectionData{}); // Set to zero injection data
            }
        }
        // NOTE: The dune broadcast() below will do something like:
        //    MPI_Bcast(inout,len,MPITraits<SlaveGroupInjectionData>::getType(),root,communicator)
        //  so it should use the custom SlaveGroupInjectionData MPI type that we defined in
        //  ReservoirCouplingMpiTraits.hpp
        this->comm().broadcast(injection_data.data(), /*count=*/num_slave_groups, /*emitter_rank=*/0);
        this->slave_group_injection_data_[this->slaveName(i)] = injection_data;
    }
}

template <class Scalar>
void
ReservoirCouplingMasterReportStep<Scalar>::
receiveProductionDataFromSlaves()
{
    auto num_slaves = this->numSlaves();
    this->logger().info("Receiving production data from slave processes");
    for (unsigned int i = 0; i < num_slaves; i++) {
        auto num_slave_groups = this->numSlaveGroups(i);
        assert( num_slave_groups > 0 );
        std::vector<SlaveGroupProductionData> production_data(num_slave_groups);
        if (this->comm().rank() == 0) {
            if (this->slaveIsActivated(i)) {
                // NOTE: See comment about error handling at the top of this file.
                auto MPI_PRODUCTION_DATA_TYPE = Dune::MPITraits<SlaveGroupProductionData>::getType();
                MPI_Recv(
                    production_data.data(),
                    /*count=*/num_slave_groups,
                    /*datatype=*/MPI_PRODUCTION_DATA_TYPE,
                    /*source_rank=*/0,
                    /*tag=*/static_cast<int>(MessageTag::SlaveProductionData),
                    this->getSlaveComm(i),
                    MPI_STATUS_IGNORE
                );
                this->logger().info(
                    fmt::format(
                        "Received production data for {} groups from slave process with name: {}. "
                        "Number of slave groups: {}", num_slave_groups, this->slaveName(i), num_slave_groups
                    )
                );
            }
            else {
                this->logger().info(fmt::format(
                    "Slave {} has not activated yet, skipping receiving production data from slave",
                        this->slaveName(i)
                    )
                );
                production_data.assign(num_slave_groups, SlaveGroupProductionData{}); // Set to zero production data
            }
        }
        // NOTE: The dune broadcast() below will do something like:
        //    MPI_Bcast(inout,len,MPITraits<SlaveGroupProductionData>::getType(),root,communicator)
        //  so it should use the custom SlaveGroupProductionData MPI type that we defined in
        //  ReservoirCouplingMpiTraits.hpp
        this->comm().broadcast(production_data.data(), /*count=*/num_slave_groups, /*emitter_rank=*/0);
        this->slave_group_production_data_[this->slaveName(i)] = production_data;
    }
}

template <class Scalar>
void
ReservoirCouplingMasterReportStep<Scalar>::
setReportStepIdx(int report_step_idx)
{
    this->report_step_idx_ = report_step_idx;
}


// ------------------
// Private methods
// ------------------


// Explicit instantiations
template class ReservoirCouplingMasterReportStep<double>;
#if FLOW_INSTANTIATE_FLOAT
template class ReservoirCouplingMasterReportStep<float>;
#endif

} // namespace Opm
