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
#include <opm/simulators/flow/rescoup/ReservoirCouplingErrorMacros.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingMpiTraits.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingMaster.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingMasterReportStep.hpp>

#include <opm/input/eclipse/Schedule/ResCoup/ReservoirCouplingInfo.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/MasterGroup.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/Slaves.hpp>
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
std::size_t
ReservoirCouplingMasterReportStep<Scalar>::
getMasterGroupCanonicalIdx(const std::string &slave_name, const std::string &master_group_name) const
{
    return this->master_.getMasterGroupCanonicalIdx(slave_name, master_group_name);
}

template <class Scalar>
Scalar
ReservoirCouplingMasterReportStep<Scalar>::
getMasterGroupInjectionSurfaceRate(const std::string &group_name, ReservoirCoupling::Phase phase) const
{
    return this->getMasterGroupRate_(group_name, phase, /*reservoir_rates=*/false, /*is_injection=*/true);
}

template <class Scalar>
Scalar
ReservoirCouplingMasterReportStep<Scalar>::
getMasterGroupInjectionReservoirRate(const std::string &group_name, ReservoirCoupling::Phase phase) const
{
    return this->getMasterGroupRate_(group_name, phase, /*reservoir_rates=*/true, /*is_injection=*/true);
}

template <class Scalar>
Scalar
ReservoirCouplingMasterReportStep<Scalar>::
getMasterGroupProductionSurfaceRate(const std::string &group_name, ReservoirCoupling::Phase phase) const
{
    return this->getMasterGroupRate_(group_name, phase, /*reservoir_rates=*/false, /*is_injection=*/false);
}

template <class Scalar>
Scalar
ReservoirCouplingMasterReportStep<Scalar>::
getMasterGroupNetworkProductionSurfaceRate(
    const std::string &group_name, ReservoirCoupling::Phase phase
) const
{
    return this->getMasterGroupRate_(
        group_name, phase, /*reservoir_rates=*/false, /*is_injection=*/false, /*network=*/true
    );
}

template <class Scalar>
Scalar
ReservoirCouplingMasterReportStep<Scalar>::
getMasterGroupProductionReservoirRate(const std::string &group_name, ReservoirCoupling::Phase phase) const
{
    return this->getMasterGroupRate_(group_name, phase, /*reservoir_rates=*/true, /*is_injection=*/false);
}

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
        RCOUP_LOG_THROW(
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
    this->logger().debug("Receiving injection data from slave processes");
    for (unsigned int i = 0; i < num_slaves; i++) {
        auto num_slave_groups = this->numSlaveGroups(i);
        if (num_slave_groups == 0) {
            // History mode: no slave groups defined, skip data exchange
            this->logger().debug(fmt::format(
                "Slave {} has no slave groups (history mode), skipping injection data exchange",
                this->slaveName(i)
            ));
            continue;
        }
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
                this->logger().debug(fmt::format(
                    "Received injection data for {} groups from {}",
                    num_slave_groups, this->slaveName(i)
                ));
            }
            else {
                this->logger().debug(fmt::format(
                    "Slave {} has not activated yet, skipping injection data",
                    this->slaveName(i)
                ));
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
    this->logger().debug("Receiving production data from slave processes");
    for (unsigned int i = 0; i < num_slaves; i++) {
        auto num_slave_groups = this->numSlaveGroups(i);
        if (num_slave_groups == 0) {
            // History mode: no slave groups defined, skip data exchange
            this->logger().debug(fmt::format(
                "Slave {} has no slave groups (history mode), skipping production data exchange",
                this->slaveName(i)
            ));
            continue;
        }
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
                this->logger().debug(fmt::format(
                    "Received production data for {} groups from {}",
                    num_slave_groups, this->slaveName(i)
                ));
            }
            else {
                this->logger().debug(fmt::format(
                    "Slave {} has not activated yet, skipping production data",
                    this->slaveName(i)
                ));
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
sendInjectionTargetsToSlave(std::size_t slave_idx,
                            const std::vector<InjectionGroupTarget>& injection_targets) const
{
    // Only rank 0 sends data to slaves. Other ranks in the master's MPI communicator
    // do not participate in master-slave communication (no else branch needed).
    if (this->comm().rank() == 0) {
        auto num_injection_targets = injection_targets.size();
        auto MPI_INJECTION_TARGETS_TYPE = Dune::MPITraits<InjectionGroupTarget>::getType();
        // NOTE: See comment about error handling at the top of this file.
        MPI_Send(
            injection_targets.data(),
            /*count=*/num_injection_targets,
            /*datatype=*/MPI_INJECTION_TARGETS_TYPE,
            /*dest_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::InjectionGroupTargets),
            this->getSlaveComm(slave_idx)
        );
        this->logger().debug(fmt::format(
            "Sent {} injection targets to {}",
            num_injection_targets, this->slaveName(slave_idx)
        ));
    }
}

// The slave process will use this information to determine how many targets to expect
template <class Scalar>
void
ReservoirCouplingMasterReportStep<Scalar>::
sendNumGroupTargetsToSlave(std::size_t slave_idx,
                           std::size_t num_injection_targets,
                           std::size_t num_production_targets) const
{
    // Only rank 0 sends data to slaves. Other ranks in the master's MPI communicator
    // do not participate in master-slave communication (no else branch needed).
    if (this->comm().rank() == 0) {
        std::vector<std::size_t> num_targets(2);
        num_targets[0] = num_injection_targets;
        num_targets[1] = num_production_targets;
        auto MPI_SIZE_T_TYPE = Dune::MPITraits<std::size_t>::getType();
        // NOTE: See comment about error handling at the top of this file.
        MPI_Send(
            num_targets.data(),
            /*count=*/2,
            /*datatype=*/MPI_SIZE_T_TYPE,
            /*dest_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::NumSlaveGroupTargets),
            this->getSlaveComm(slave_idx)
        );
        this->logger().debug(fmt::format(
            "Sent target counts (inj={}, prod={}) to {}",
            num_injection_targets, num_production_targets, this->slaveName(slave_idx)
        ));
    }
}

template <class Scalar>
void
ReservoirCouplingMasterReportStep<Scalar>::
sendProductionTargetsToSlave(std::size_t slave_idx,
                             const std::vector<ProductionGroupTarget>& production_targets) const
{
    // Only rank 0 sends data to slaves. Other ranks in the master's MPI communicator
    // do not participate in master-slave communication (no else branch needed).
    if (this->comm().rank() == 0) {
        auto num_production_targets = production_targets.size();
        auto MPI_PRODUCTION_TARGETS_TYPE = Dune::MPITraits<ProductionGroupTarget>::getType();
        // NOTE: See comment about error handling at the top of this file.
        MPI_Send(
            production_targets.data(),
            /*count=*/num_production_targets,
            /*datatype=*/MPI_PRODUCTION_TARGETS_TYPE,
            /*dest_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::ProductionGroupTargets),
            this->getSlaveComm(slave_idx)
        );
        this->logger().debug(fmt::format(
            "Sent {} production targets to {}",
            num_production_targets, this->slaveName(slave_idx)
        ));
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

template <class Scalar>
Scalar
ReservoirCouplingMasterReportStep<Scalar>::
getMasterGroupRate_(const std::string &group_name, ReservoirCoupling::Phase phase,
                    bool reservoir_rates, bool is_injection, bool network) const
{
    auto it = this->getMasterGroupToSlaveNameMap().find(group_name);
    if (it != this->getMasterGroupToSlaveNameMap().end()) {
        auto& slave_name = it->second;
        auto group_idx = this->getMasterGroupCanonicalIdx(slave_name, group_name);
        if (is_injection) {
            const auto& rates = reservoir_rates
                ? this->slave_group_injection_data_.at(slave_name)[group_idx].reservoir_rates
                : this->slave_group_injection_data_.at(slave_name)[group_idx].surface_rates;
            return rates[phase];
        }
        else {
            const auto& prod_data = this->slave_group_production_data_.at(slave_name)[group_idx];
            if (reservoir_rates) {
                return prod_data.reservoir_rates[phase];
            }
            else if (network) {
                return prod_data.network_surface_rates[phase];
            }
            else {
                return prod_data.surface_rates[phase];
            }
        }
    }
    else {
        RCOUP_LOG_THROW(
            std::runtime_error,
            fmt::format(
                "Master group name {} not found in master-to-slave-group-name mapping",
                group_name
            )
        );
    }
}

// Explicit instantiations
template class ReservoirCouplingMasterReportStep<double>;
#if FLOW_INSTANTIATE_FLOAT
template class ReservoirCouplingMasterReportStep<float>;
#endif

} // namespace Opm
