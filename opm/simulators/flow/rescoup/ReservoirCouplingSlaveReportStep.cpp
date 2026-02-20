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
hasMasterInjectionTarget(const std::string& gname, Phase phase) const
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
std::pair<Scalar, Group::InjectionCMode>
ReservoirCouplingSlaveReportStep<Scalar>::
masterInjectionTarget(const std::string& gname, Phase phase) const
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
setMasterInjectionTarget(const std::string& gname, Phase phase, Scalar target, Group::InjectionCMode cmode)
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
setMasterProductionTarget(const std::string& gname, Scalar target, Group::ProductionCMode cmode)
{
    this->master_production_targets_[gname] = {target, cmode};
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
