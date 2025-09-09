/*
  Copyright 2024 Equinor ASA

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

#include <opm/simulators/flow/ReservoirCouplingMaster.hpp>
#include <opm/simulators/flow/ReservoirCouplingMpiTraits.hpp>
#include <opm/simulators/flow/ReservoirCouplingSpawnSlaves.hpp>

#include <opm/input/eclipse/Schedule/ResCoup/ReservoirCouplingInfo.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/MasterGroup.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/Slaves.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <opm/input/eclipse/Schedule/Schedule.hpp>


#include <filesystem>
#include <vector>

#include <fmt/format.h>

namespace Opm {
// NOTE: All slave-master communicators have set a custom error handler, which eventually
//   will call MPI_Abort() so there is no need to check the return value of any MPI_Recv()
//   or MPI_Send() calls below.

template <class Scalar>
ReservoirCouplingMaster<Scalar>::
ReservoirCouplingMaster(
    const Parallel::Communication &comm,
    const Schedule &schedule,
    int argc, char **argv
) :
    comm_{comm},
    schedule_{schedule},
    argc_{argc},
    argv_{argv}
{
    this->activation_date_ = this->getMasterActivationDate_();
}

// ------------------
// Public methods
// ------------------

template <class Scalar>
bool
ReservoirCouplingMaster<Scalar>::
isMasterGroup(const std::string &group_name) const
{
    return this->master_group_slave_names_.find(group_name) !=
           this->master_group_slave_names_.end();
}

template <class Scalar>
const std::vector<std::string>&
ReservoirCouplingMaster<Scalar>::
getMasterGroupNamesForSlave(std::size_t slave_idx) const
{
    if (slave_idx >= this->slave_idx_to_master_groups_.size()) {
        OPM_THROW(
            std::runtime_error,
            fmt::format(
                "Slave index {} out of bounds. Valid range is [0, {})",
                slave_idx, this->slave_idx_to_master_groups_.size()
            )
        );
    }

    // Direct O(1) vector access for performance-critical usage
    return this->slave_idx_to_master_groups_[slave_idx];
}

template <class Scalar>
std::size_t
ReservoirCouplingMaster<Scalar>::
getMasterGroupPotIdx(
    const std::string &slave_name, const std::string &master_group_name) const
{
    return this->master_group_name_order_.at(slave_name).at(master_group_name);
}

template <class Scalar>
const ReservoirCoupling::Potentials<Scalar>&
ReservoirCouplingMaster<Scalar>::
getSlaveGroupPotentials(const std::string &master_group_name)
{
    auto it = this->master_group_slave_names_.find(master_group_name);
    if (it != this->master_group_slave_names_.end()) {
        auto& slave_name = it->second;
        auto pot_idx = this->getMasterGroupPotIdx(slave_name, master_group_name);
        return this->slave_group_potentials_[slave_name][pot_idx];
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
int
ReservoirCouplingMaster<Scalar>::
getSlaveIdx(const std::string &slave_name) const
{
    auto it = std::find(this->slave_names_.begin(), this->slave_names_.end(), slave_name);
    if (it != this->slave_names_.end()) {
        return std::distance(this->slave_names_.begin(), it);
    }
    OPM_THROW(std::runtime_error, "Slave name not found: " + slave_name);
}

template <class Scalar>
void
ReservoirCouplingMaster<Scalar>::
maybeActivate(int report_step) {
    if (!this->activated()) {
        double start_date = this->schedule_.getStartTime();
        auto current_time = start_date + this->schedule_.seconds(report_step);
        if (Seconds::compare_gt_or_eq(current_time, this->activation_date_)) {
            this->activated_ = true;
        }
    }
}

template <class Scalar>
void
ReservoirCouplingMaster<Scalar>::
maybeSpawnSlaveProcesses(int report_step)
{
    if (this->numSlavesStarted() > 0) {  // We have already spawned the slave processes
        return;
    }
    const auto& rescoup = this->schedule_[report_step].rescoup();
    auto slave_count = rescoup.slaveCount();
    auto master_group_count = rescoup.masterGroupCount();
    if (slave_count > 0 && master_group_count > 0) {
        ReservoirCouplingSpawnSlaves<Scalar> spawn_slaves{*this, rescoup};
        spawn_slaves.spawn();
    }
}

template <class Scalar>
double
ReservoirCouplingMaster<Scalar>::
maybeChopSubStep(double suggested_timestep_original, double elapsed_time) const
{
    // Check if the suggested timestep needs to be adjusted based on the slave processes'
    // next report step, or if the slave process has not started yet: the start of a slave process.
    // NOTE: getStartTime() returns a std::time_t value, which is typically a long integer. It should
    //     be possible to represent reasonable epoch values within a double. See comment for
    //     getMasterActivationDate_() for more information.
    double start_date = this->schedule_.getStartTime();
    double step_start_date{start_date + elapsed_time};
    double step_end_date{step_start_date + suggested_timestep_original};
    double suggested_timestep{suggested_timestep_original};
    auto num_slaves = this->numSlavesStarted();
    // Determine the minimum step_end_date and the corresponding suggested_timestep such that no
    // slave process will report or start during the timestep [step_start_date, step_end_date]
    // where suggested_timestep = step_end_date - step_start_date
    for (std::size_t i = 0; i < num_slaves; i++) {
        double slave_start_date = this->slave_start_dates_[i];
        double slave_activation_date = this->slave_activation_dates_[i];
        double slave_next_report_date{this->slave_next_report_time_offsets_[i] + slave_start_date};
        if (Seconds::compare_gt_or_eq(slave_activation_date, step_end_date)) {
            // The slave process has not activated yet, and will not activate during this timestep
            continue;
        }
        double slave_elapsed_time;
        if (Seconds::compare_lt_or_eq(slave_activation_date,step_start_date)) {
            // The slave process has already activated, and will continue during this timestep
            if (Seconds::compare_gt(slave_next_report_date, step_end_date)) {
                // The slave process will not report during this timestep
                continue;
            }
            // The slave process will report during this timestep
            slave_elapsed_time = slave_next_report_date - step_start_date;
        }
        else {
            // The slave process will activate during the timestep, but not at the beginning
            // Ensure that we synchronize the slave activation report date with the
            //   start of a master timestep by chopping the current master timestep.
            slave_elapsed_time = slave_activation_date - step_start_date;
        }
        suggested_timestep = slave_elapsed_time;
        step_end_date = step_start_date + suggested_timestep;
    }
    return suggested_timestep;
}

template <class Scalar>
void
ReservoirCouplingMaster<Scalar>::
maybeReceiveActivationHandshakeFromSlaves(double current_time)
{
    // Initialize on first call
    if (this->slave_activation_status_.empty()) {
        this->slave_activation_status_.resize(this->numSlavesStarted(), false);
    }

    if (this->comm_.rank() == 0) {
        auto current_date = this->schedule_.getStartTime() + current_time;
        for (unsigned int i = 0; i < this->numSlavesStarted(); i++) {
            // Skip if already activated
            if (this->slaveIsActivated(i)) {
                continue;
            }
            // Check if slave should activate during this timestep
            double slave_activation_date = this->slave_activation_dates_[i];
            // NOTE: The master will adjust its time stepping to ensure that its step will always
            //    conincide with the slave report steps (and hence the slave activation date)
            assert(Seconds::compare_gt_or_eq(slave_activation_date, current_date));
            if (Seconds::compare_eq(slave_activation_date, current_date)) {
                // Use non-blocking probe first to check for handshake
                int flag;
                MPI_Status status;
                MPI_Iprobe(0, static_cast<int>(MessageTag::SlaveActivationHandshake),
                          this->master_slave_comm_[i], &flag, &status);

                if (!flag) {
                    // Handshake not received yet, for debugging reasons it can be useful to know which slave
                    //   is waiting for the handshake
                    OpmLog::info(fmt::format("Waiting for activation handshake from slave {}",
                        this->slave_names_[i]));
                }
                std::uint8_t handshake;
                auto MPI_UINT8_T_TYPE = Dune::MPITraits<std::uint8_t>::getType();
                MPI_Recv(&handshake, 1, MPI_UINT8_T_TYPE, 0,
                        static_cast<int>(MessageTag::SlaveActivationHandshake),
                        this->master_slave_comm_[i], MPI_STATUS_IGNORE);
                this->slave_activation_status_[i] = handshake;
                OpmLog::info(fmt::format("Received activation handshake from slave {}",
                                        this->slave_names_[i]));

            }
        }
    }
    // Broadcast status to all ranks
    comm_.broadcast(
        this->slave_activation_status_.data(), /*count=*/this->numSlavesStarted(), /*emitter_rank=*/0
    );
}

template <class Scalar>
std::size_t
ReservoirCouplingMaster<Scalar>::
numSlaveGroups(unsigned int index)
{
    return this->master_group_name_order_[this->slave_names_[index]].size();
}

template <class Scalar>
std::size_t
ReservoirCouplingMaster<Scalar>::
numSlavesStarted() const
{
    return this->slave_names_.size();
}

template <class Scalar>
void
ReservoirCouplingMaster<Scalar>::
sendNextTimeStepToSlaves(double timestep)
{
    if (this->comm_.rank() == 0) {
        for (unsigned int i = 0; i < this->master_slave_comm_.size(); i++) {
            if (!this->slaveIsActivated(i)) {
                OpmLog::info(fmt::format(
                    "Slave {} has not activated yet, skipping sending next time step to slave",
                    this->slave_names_[i]
                ));
                continue;
            }
            // NOTE: See comment about error handling at the top of this file.
            MPI_Send(
                &timestep,
                /*count=*/1,
                /*datatype=*/MPI_DOUBLE,
                /*dest_rank=*/0,
                /*tag=*/static_cast<int>(MessageTag::SlaveNextTimeStep),
                this->getSlaveComm(i)
            );
            OpmLog::info(fmt::format(
                "Sent next time step {} from master process rank 0 to slave process "
                "rank 0 with name: {}", timestep, this->slave_names_[i])
            );
        }
   }
}

template <class Scalar>
void
ReservoirCouplingMaster<Scalar>::
rebuildSlaveIdxToMasterGroupsVector()
{
    // Rebuild the index-based vector from the name-based map
    // This ensures consistent ordering with slave_names_ vector
    this->slave_idx_to_master_groups_.clear();
    this->slave_idx_to_master_groups_.resize(this->slave_names_.size());

    for (std::size_t slave_idx = 0; slave_idx < this->slave_names_.size(); ++slave_idx) {
        const auto& slave_name = this->slave_names_[slave_idx];
        auto it = this->slave_name_to_master_groups_map_.find(slave_name);
        if (it != this->slave_name_to_master_groups_map_.end()) {
            this->slave_idx_to_master_groups_[slave_idx] = it->second;
        }
        // If slave_name not found in map, leave the vector empty for that index
    }
}

template <class Scalar>
void
ReservoirCouplingMaster<Scalar>::
receiveNextReportDateFromSlaves()
{
    auto num_slaves = this->numSlavesStarted();
    if (this->comm_.rank() == 0) {
        OpmLog::info("Receiving next report dates from slave processes");
        for (unsigned int i = 0; i < num_slaves; i++) {
            if (!this->slaveIsActivated(i)) {
                // Set to zero to indicate that the slave has not activated yet
                this->slave_next_report_time_offsets_[i] = 0.0;
                OpmLog::info(fmt::format("Slave {} has not activated yet, setting next report date to 0.0",
                                        this->slave_names_[i]));
                continue;
            }
            double slave_next_report_time_offset; // Elapsed time from the beginning of the simulation
            // NOTE: See comment about error handling at the top of this file.
            MPI_Recv(
                &slave_next_report_time_offset,
                /*count=*/1,
                /*datatype=*/MPI_DOUBLE,
                /*source_rank=*/0,
                /*tag=*/static_cast<int>(MessageTag::SlaveNextReportDate),
                this->getSlaveComm(i),
                MPI_STATUS_IGNORE
            );
            this->slave_next_report_time_offsets_[i] = slave_next_report_time_offset;
            OpmLog::info(
                fmt::format(
                    "Received simulation slave next report date from slave process with name: {}. "
                    "Next report date: {}", this->slave_names_[i], slave_next_report_time_offset
                )
            );
        }
    }
    this->comm_.broadcast(
        this->slave_next_report_time_offsets_.data(), /*count=*/num_slaves, /*emitter_rank=*/0
    );
    OpmLog::info("Broadcasted slave next report dates to all ranks");
}

template <class Scalar>
void
ReservoirCouplingMaster<Scalar>::
receivePotentialsFromSlaves()
{
    auto num_slaves = this->numSlavesStarted();
    this->logger_.info("Receiving potentials from slave processes");
    for (unsigned int i = 0; i < num_slaves; i++) {
        auto num_slave_groups = this->numSlaveGroups(i);
        assert( num_slave_groups > 0 );
        std::vector<Potentials> potentials(num_slave_groups);
        if (this->comm_.rank() == 0) {
            if (this->slaveIsActivated(i)) {
                // NOTE: See comment about error handling at the top of this file.
                auto MPI_POTENTIALS_TYPE = Dune::MPITraits<Potentials>::getType();
                MPI_Recv(
                    potentials.data(),
                    /*count=*/num_slave_groups,
                    /*datatype=*/MPI_POTENTIALS_TYPE,
                    /*source_rank=*/0,
                    /*tag=*/static_cast<int>(MessageTag::Potentials),
                    this->getSlaveComm(i),
                    MPI_STATUS_IGNORE
                );
                this->logger_.info(
                    fmt::format(
                        "Received potentials from slave process with name: {}. "
                        "Number of slave groups: {}", this->slave_names_[i], num_slave_groups
                    )
                );
            }
            else {
                this->logger_.info(fmt::format(
                    "Slave {} has not activated yet, skipping receiving potentials from slave",
                        this->slave_names_[i]
                    )
                );
                potentials.assign(num_slave_groups, Potentials{}); // Set to zero potentials
            }
        }
        // NOTE: The dune broadcast() below will do something like:
        //    MPI_Bcast(inout,len,MPITraits<Potentials>::getType(),root,communicator)
        //  so it should use the custom Potentials MPI type that we defined in
        //  ReservoirCouplingMpiTraits.hpp
        this->comm_.broadcast(potentials.data(), /*count=*/num_slave_groups, /*emitter_rank=*/0);
        this->slave_group_potentials_[this->slave_names_[i]] = potentials;
    }
}

template <class Scalar>
void
ReservoirCouplingMaster<Scalar>::
sendInjectionTargetsToSlave(std::size_t slave_idx,
                            const std::vector<InjectionGroupTarget>& injection_targets) const
{
    // Only rank 0 sends data to slaves. Other ranks in the master's MPI communicator
    // do not participate in master-slave communication (no else branch needed).
    if (this->comm_.rank() == 0) {
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
        this->logger_.info(fmt::format(
            "Sent {} injection targets to slave process with name: {}",
            num_injection_targets, this->getSlaveName(slave_idx)
        ));
    }
}

// The slave process will use this information to determine how many targets to expect
template <class Scalar>
void
ReservoirCouplingMaster<Scalar>::
sendNumGroupTargetsToSlave(std::size_t slave_idx,
                           std::size_t num_injection_targets,
                           std::size_t num_production_targets) const
{
    // Only rank 0 sends data to slaves. Other ranks in the master's MPI communicator
    // do not participate in master-slave communication (no else branch needed).
    if (this->comm_.rank() == 0) {
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
        this->logger_.info(fmt::format(
            "Sent number of injection targets {} and production targets {} to slave process with name: {}",
            num_injection_targets, num_production_targets, this->getSlaveName(slave_idx)
        ));
    }
}

template <class Scalar>
void
ReservoirCouplingMaster<Scalar>::
sendProductionTargetsToSlave(std::size_t slave_idx,
                             const std::vector<ProductionGroupTarget>& production_targets) const
{
    // Only rank 0 sends data to slaves. Other ranks in the master's MPI communicator
    // do not participate in master-slave communication (no else branch needed).
    if (this->comm_.rank() == 0) {
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
        this->logger_.info(fmt::format(
            "Sent {} production targets to slave process with name: {}", num_production_targets, this->getSlaveName(slave_idx)
        ));
    }
}

template <class Scalar>
void
ReservoirCouplingMaster<Scalar>::
updateMasterGroupNameOrderMap(
    const std::string& slave_name, const std::map<std::string, std::size_t>& master_group_name_map
)
{
    this->master_group_name_order_[slave_name] = master_group_name_map;
}

// ------------------
// Private methods
// ------------------

template <class Scalar>
double
ReservoirCouplingMaster<Scalar>::
getMasterActivationDate_() const
{
    // Assume master mode is activated when the first SLAVES keyword is encountered in the schedule
    // NOTE: getStartTime() returns a std::time_t value, which is typically a long integer representing
    //     the number of seconds since the epoch (1970-01-01 00:00:00 UTC)
    //     The maximum integer that can be represented by a double is 2^53 - 1, which is approximately
    //     9e15. This corresponds to a date in the year 2.85e8 or 285 million years into the future.
    //     So we should be able to represent reasonable epoch values within a double.
    double start_date = this->schedule_.getStartTime();
    for (std::size_t report_step = 0; report_step < this->schedule_.size(); ++report_step) {
        auto rescoup = this->schedule_[report_step].rescoup();
        if (rescoup.slaveCount() > 0) {
            return start_date + this->schedule_.seconds(report_step);
        }
    }
    // NOTE: Consistency between SLAVES and GRUPMAST keywords has already been checked in
    //       init() in SimulatorFullyImplicitBlackoil.hpp
    OPM_THROW(std::runtime_error, "Reservoir coupling: Failed to find master activation time: "
              "No SLAVES keyword found in schedule");
}

template class ReservoirCouplingMaster<double>;

#if FLOW_INSTANTIATE_FLOAT
template class ReservoirCouplingMaster<float>;
#endif

} // namespace Opm

