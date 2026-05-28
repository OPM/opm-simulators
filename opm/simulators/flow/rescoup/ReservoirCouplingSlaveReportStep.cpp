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
bool
ReservoirCouplingSlaveReportStep<Scalar>::
hasMasterGroupNodePressure(const std::string& gname) const
{
    return this->master_group_node_pressures_.count(gname) > 0;
}

template <class Scalar>
bool
ReservoirCouplingSlaveReportStep<Scalar>::
hasMasterInjectionTarget(const std::string& gname, const Phase phase) const
{
    return this->master_injection_targets_.count({phase, gname}) > 0;
}

template <class Scalar>
bool
ReservoirCouplingSlaveReportStep<Scalar>::
hasMasterProductionLimits(const std::string& gname) const
{
    return this->master_production_limits_.count(gname) > 0;
}

template <class Scalar>
bool
ReservoirCouplingSlaveReportStep<Scalar>::
hasMasterProductionTarget(const std::string& gname) const
{
    return this->master_production_targets_.count(gname) > 0;
}

template <class Scalar>
Scalar
ReservoirCouplingSlaveReportStep<Scalar>::
masterGroupNodePressure(const std::string& gname) const
{
    return this->master_group_node_pressures_.at(gname);
}

template <class Scalar>
std::pair<Scalar, Group::InjectionCMode>
ReservoirCouplingSlaveReportStep<Scalar>::
masterInjectionTarget(const std::string& gname, const Phase phase) const
{
    return this->master_injection_targets_.at({phase, gname});
}

template <class Scalar>
const typename ReservoirCouplingSlaveReportStep<Scalar>::MasterProductionLimits&
ReservoirCouplingSlaveReportStep<Scalar>::
masterProductionLimits(const std::string& gname) const
{
    return this->master_production_limits_.at(gname);
}

template <class Scalar>
std::pair<Scalar, Group::ProductionCMode>
ReservoirCouplingSlaveReportStep<Scalar>::
masterProductionTarget(const std::string& gname) const
{
    return this->master_production_targets_.at(gname);
}

// Receive the master's single-bool flag telling the slave whether the master
// will iterate the cross-rescoup network exchange this sync timestep.
// Mirrored into `last_received_master_group_node_pressures_is_final_`:
// is_final = !active, so when the master is active the slave's outer-loop
// gate sees "not final" and participates in the per-iteration node-pressure
// receive in updateWellControlsAndNetworkIteration.
template <class Scalar>
void
ReservoirCouplingSlaveReportStep<Scalar>::
receiveCoupledNetworkActiveStatusFromMaster()
{
    std::size_t payload = 0u;
    if (this->comm().rank() == 0) {
        auto MPI_SIZE_T_TYPE = Dune::MPITraits<std::size_t>::getType();
        // NOTE: See comment about error handling at the top of this file.
        MPI_Recv(
            &payload,
            /*count=*/1,
            /*datatype=*/MPI_SIZE_T_TYPE,
            /*source_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::CoupledNetworkActiveStatus),
            this->getSlaveMasterComm(),
            MPI_STATUS_IGNORE
        );
    }
    this->comm().broadcast(&payload, /*count=*/1, /*emitter_rank=*/0);
    // The received flag is now per-slave: "is this slave connected to the
    // master's cross-rescoup network this sync step?".
    const bool active = payload != 0u;
    this->connected_to_master_coupled_network_ = active;
    this->last_received_master_group_node_pressures_is_final_ = !active;
    // No coupled-network iteration this sync step: the per-iteration pressure
    // receive (which clears the map) is not entered, so drop any pressures
    // left over from a previous sync step to avoid re-applying stale THP limits.
    if (!active) {
        this->master_group_node_pressures_.clear();
    }
    this->logger().debug(fmt::format(
        "Received coupled-network active status (active={}) from master process",
        active));
}

template <class Scalar>
void
ReservoirCouplingSlaveReportStep<Scalar>::
receiveInjectionGroupTargetsFromMaster(std::size_t num_targets)
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
        this->logger().debug(fmt::format(
            "Received injection {} group targets from master process rank 0", num_targets
        ));
    }
    this->comm().broadcast(
        injection_targets.data(), num_targets, /*emitter_rank=*/0
    );
    // Clear old targets before storing new ones from master
    this->master_injection_targets_.clear();
    for (const auto& target : injection_targets) {
        const auto& group_name = this->slave_.slaveGroupIdxToGroupName(target.group_name_idx);
        // Convert ReservoirCoupling::Phase to Opm::Phase
        const auto opm_phase = ReservoirCoupling::convertToOpmPhase(target.phase);
        this->setMasterInjectionTarget(
            group_name, opm_phase, target.target, target.cmode
        );
        this->logger().debug(fmt::format(
            "Stored master injection target for group '{}': target={}, phase={}",
            group_name, target.target, static_cast<int>(target.phase)
        ));
    }
}

template <class Scalar>
void
ReservoirCouplingSlaveReportStep<Scalar>::
receiveMasterGroupNodePressuresFromMaster(std::size_t num_pressures)
{
    std::vector<MasterGroupNodePressure> pressures(num_pressures);
    if (this->comm().rank() == 0) {
        auto MPI_MASTER_GROUP_NODE_PRESSURE_TYPE =
            Dune::MPITraits<MasterGroupNodePressure>::getType();
        // NOTE: See comment about error handling at the top of this file.
        MPI_Recv(
            pressures.data(),
            /*count=*/num_pressures,
            /*datatype=*/MPI_MASTER_GROUP_NODE_PRESSURE_TYPE,
            /*source_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::MasterGroupNodePressures),
            this->getSlaveMasterComm(),
            MPI_STATUS_IGNORE
        );
        this->logger().debug(fmt::format(
            "Received {} master group node pressures from master process rank 0",
            num_pressures
        ));
    }
    this->comm().broadcast(
        pressures.data(), num_pressures, /*emitter_rank=*/0
    );
    // Clear old pressures before storing new ones from master
    this->master_group_node_pressures_.clear();
    for (const auto& entry : pressures) {
        const auto& group_name = this->slave_.slaveGroupIdxToGroupName(entry.group_name_idx);
        this->setMasterGroupNodePressure(group_name, entry.pressure);
        this->logger().debug(fmt::format(
            "Stored master group node pressure for group '{}': pressure={}",
            group_name, entry.pressure
        ));
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
        this->logger().debug("Received number of slave group constraints from master process rank 0");
    }
    this->comm().broadcast(num_group_targets.data(), /*count=*/2, /*emitter_rank=*/0);
    auto num_injection_targets = num_group_targets[0];
    auto num_production_constraints = num_group_targets[1];
    this->logger().debug(fmt::format("Received number of injection targets: {} and "
                             "production constraints: {} from master process",
                             num_injection_targets, num_production_constraints));
    return std::make_pair(num_injection_targets, num_production_constraints);
}

template <class Scalar>
std::pair<std::size_t, bool>
ReservoirCouplingSlaveReportStep<Scalar>::
receiveNumMasterGroupNodePressuresFromMaster()
{
    std::vector<std::size_t> header(2);
    if (this->comm().rank() == 0) {
        auto MPI_SIZE_T_TYPE = Dune::MPITraits<std::size_t>::getType();
        // NOTE: See comment about error handling at the top of this file.
        MPI_Recv(
            header.data(),
            /*count=*/2,
            /*datatype=*/MPI_SIZE_T_TYPE,
            /*source_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::NumMasterGroupNodePressures),
            this->getSlaveMasterComm(),
            MPI_STATUS_IGNORE
        );
        this->logger().debug("Received master-group-node-pressures header from master process rank 0");
    }
    this->comm().broadcast(header.data(), /*count=*/2, /*emitter_rank=*/0);
    const auto num_pressures = header[0];
    const bool is_final = header[1] != 0u;
    this->last_received_master_group_node_pressures_is_final_ = is_final;
    // When the master sends count == 0 the payload receive (which clears and
    // repopulates the map) is skipped by the caller. Clear here so stale
    // pressures from a previous iteration are not re-applied as THP limits.
    if (num_pressures == 0) {
        this->master_group_node_pressures_.clear();
    }
    this->logger().debug(fmt::format(
        "Received master-group-node-pressures header (count={}, is_final={}) from master process",
        num_pressures, is_final
    ));
    return {num_pressures, is_final};
}

template <class Scalar>
void
ReservoirCouplingSlaveReportStep<Scalar>::
receiveProductionGroupConstraintsFromMaster(std::size_t num_targets)
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
        this->logger().debug(fmt::format(
            "Received production {} group constraints from master process rank 0", num_targets
        ));
    }
    this->comm().broadcast(
        production_constraints.data(), num_targets, /*emitter_rank=*/0
    );
    // Clear old targets and limits before storing new ones from master
    this->master_production_targets_.clear();
    this->master_production_limits_.clear();
    for (const auto& target : production_constraints) {
        const auto& group_name = this->slave_.slaveGroupIdxToGroupName(target.group_name_idx);
        this->setMasterProductionTarget(
            group_name, target.target, target.cmode
        );
        MasterProductionLimits limits;
        limits.oil_limit = target.oil_limit;
        limits.water_limit = target.water_limit;
        limits.gas_limit = target.gas_limit;
        limits.liquid_limit = target.liquid_limit;
        limits.resv_limit = target.resv_limit;
        this->setMasterProductionLimits(group_name, limits);
        this->logger().debug(fmt::format(
            "Stored master production target for group '{}': target={}, cmode={}",
            group_name, target.target, static_cast<int>(target.cmode)
        ));
    }
}

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

template <class Scalar>
void
ReservoirCouplingSlaveReportStep<Scalar>::
setMasterGroupNodePressure(const std::string& gname, const Scalar pressure)
{
    this->master_group_node_pressures_[gname] = pressure;
}

template <class Scalar>
void
ReservoirCouplingSlaveReportStep<Scalar>::
setMasterInjectionTarget(
    const std::string& gname, const Phase phase, const Scalar target, const Group::InjectionCMode cmode
)
{
    this->master_injection_targets_[{phase, gname}] = {target, cmode};
}

template <class Scalar>
void
ReservoirCouplingSlaveReportStep<Scalar>::
setMasterProductionLimits(const std::string& gname, const MasterProductionLimits& limits)
{
    this->master_production_limits_[gname] = limits;
}

template <class Scalar>
void
ReservoirCouplingSlaveReportStep<Scalar>::
setMasterProductionTarget(const std::string& gname, const Scalar target, const Group::ProductionCMode cmode)
{
    this->master_production_targets_[gname] = {target, cmode};
}

template <class Scalar>
void
ReservoirCouplingSlaveReportStep<Scalar>::
markSlaveGroupsInSchedule(Schedule& schedule, const int report_step_idx)
{
    for (std::size_t i = 0; i < this->slave_.numSlaveGroups(); ++i) {
        const auto& gname = this->slave_.slaveGroupIdxToGroupName(i);
        if (this->hasMasterProductionTarget(gname)) {
            schedule.markSlaveProductionGroup(report_step_idx, gname);
        }
        for (const auto phase : {Phase::WATER, Phase::OIL, Phase::GAS}) {
            if (this->hasMasterInjectionTarget(gname, phase)) {
                schedule.markSlaveInjectionGroup(report_step_idx, gname);
                break;  // only need to mark once per group
            }
        }
    }
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
        this->logger().debug(
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
