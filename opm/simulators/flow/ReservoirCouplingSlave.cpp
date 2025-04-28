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
#include <opm/simulators/flow/ReservoirCoupling.hpp>
#include <opm/simulators/flow/ReservoirCouplingMpiTraits.hpp>
#include <opm/simulators/flow/ReservoirCouplingSlave.hpp>

#include <opm/input/eclipse/Schedule/ResCoup/ReservoirCouplingInfo.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/MasterGroup.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/Slaves.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <dune/common/parallel/mpitraits.hh>

#include <vector>
#include <fmt/format.h>

namespace Opm {
// NOTE: All slave-master communicators have set a custom error handler, which eventually
//   will call MPI_Abort() so there is no need to check the return value of any MPI_Recv()
//   or MPI_Send() calls below.

ReservoirCouplingSlave::
ReservoirCouplingSlave(
    const Parallel::Communication &comm,
    const Schedule &schedule,
    const SimulatorTimer &timer
) :
    comm_{comm},
    schedule_{schedule},
    timer_{timer}
{
    this->slave_master_comm_ = MPI_COMM_NULL;
    MPI_Comm_get_parent(&this->slave_master_comm_);
    if (this->slave_master_comm_ == MPI_COMM_NULL) {
        OPM_THROW(std::runtime_error, "Slave process is not spawned by a master process");
    }
    // NOTE: By installing a custom error handler for all slave-master communicators, which
    //   eventually will call MPI_Abort(), there is no need to check the return value of any
    //   MPI_Recv() or MPI_Send() calls as errors will be caught by the error handler.
    ReservoirCoupling::setErrhandler(this->slave_master_comm_, /*is_master=*/false);
}

void
ReservoirCouplingSlave::
sendAndReceiveInitialData() {
    this->sendActivationDateToMasterProcess_();
    this->sendSimulationStartDateToMasterProcess_();
    this->receiveSlaveNameFromMasterProcess_();
    this->receiveMasterGroupNamesFromMasterProcess_();
}

double
ReservoirCouplingSlave::
receiveNextTimeStepFromMaster() {
    double timestep;
    if (this->comm_.rank() == 0) {
        // NOTE: See comment about error handling at the top of this file.
        MPI_Recv(
            &timestep,
            /*count=*/1,
            /*datatype=*/MPI_DOUBLE,
            /*source_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::SlaveNextTimeStep),
            this->slave_master_comm_,
            MPI_STATUS_IGNORE
        );
        OpmLog::info(
            fmt::format("Slave rank 0 received next timestep {} from master.", timestep)
        );
    }
    this->comm_.broadcast(&timestep, /*count=*/1, /*emitter_rank=*/0);
    OpmLog::info("Broadcasted slave next time step to all ranks");
    return timestep;
}


void
ReservoirCouplingSlave::
sendNextReportDateToMasterProcess() const
{
    if (this->comm_.rank() == 0) {
        double elapsed_time = this->timer_.simulationTimeElapsed();
        double current_step_length = this->timer_.currentStepLength();
        // NOTE: This is an offset in seconds from the start date, so it will be 0 if the next report
        //      would be the start date. In general, it should be a positive number.
        double next_report_time_offset = elapsed_time + current_step_length;
        // NOTE: See comment about error handling at the top of this file.
        MPI_Send(
            &next_report_time_offset,
            /*count=*/1,
            /*datatype=*/MPI_DOUBLE,
            /*dest_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::SlaveNextReportDate),
            this->slave_master_comm_
        );
        OpmLog::info("Sent next report date to master process from rank 0");
   }
}


void
ReservoirCouplingSlave::
sendPotentialsToMaster(const std::vector<Potentials> &potentials) const
{
    // NOTE: The master can determine from the ordering of the potentials in the vector
    //   which slave group for a given slave name the given potentials belong to,
    //   so we do not need to send the slave group names also.
    if (this->comm_.rank() == 0) {
        //auto num_groups = potentials.size();
        auto MPI_POTENTIALS_TYPE = Dune::MPITraits<Potentials>::getType();
        MPI_Send(
            potentials.data(),
            /*count=*/potentials.size(),
            /*datatype=*/MPI_POTENTIALS_TYPE,
            /*dest_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::Potentials),
            this->slave_master_comm_
        );
        this->logger_.info(
            "Sent potentials to master process from rank 0, slave name: " + this->slave_name_
        );
    }
}



// ------------------
// Private methods
// ------------------

void
ReservoirCouplingSlave::
checkGrupSlavGroupNames_()
{
    // Validate that each slave group name has a corresponding master group name
    bool grup_slav_found = false;
    for (std::size_t report_step = 0; report_step < this->schedule_.size(); ++report_step) {
        auto rescoup = this->schedule_[report_step].rescoup();
        if (rescoup.grupSlavCount() > 0) {
            grup_slav_found = true;
            auto grup_slavs = rescoup.grupSlavs();
            for (const auto& [slave_group_name, grup_slav] : grup_slavs) {
                auto map_iter = this->slave_to_master_group_map_.find(slave_group_name);
                if (map_iter == this->slave_to_master_group_map_.end()) {
                    OPM_THROW(std::runtime_error,
                              "Reservoir coupling: Failed to find master group name for slave group: "
                              + slave_group_name);
                }
                else {
                    const auto& master_group_name = map_iter->second;
                    if (grup_slav.masterGroupName() != master_group_name) {
                        OPM_THROW(std::runtime_error,
                                  "Reservoir coupling: Inconsistent master group name for slave group: "
                                  + slave_group_name);
                    }
                }
            }
        }
    }
    if (!grup_slav_found) {
        OPM_THROW(std::runtime_error, "Reservoir coupling: Failed to find slave group names: "
                  "No GRUPSLAV keyword found in schedule");
    }
}

double
ReservoirCouplingSlave::
getGrupSlavActivationDate_() const
{
    double start_date = this->schedule_.getStartTime();
    for (std::size_t report_step = 0; report_step < this->schedule_.size(); ++report_step) {
        auto rescoup = this->schedule_[report_step].rescoup();
        if (rescoup.grupSlavCount() > 0) {
            return start_date + this->schedule_.seconds(report_step);
        }
    }
    OPM_THROW(std::runtime_error, "Reservoir coupling: Failed to find slave activation time: "
              "No GRUPSLAV keyword found in schedule");
}

void
ReservoirCouplingSlave::
maybeActivate(int report_step) {
    if (!this->activated()) {
        auto rescoup = this->schedule_[report_step].rescoup();
        if (rescoup.grupSlavCount() > 0) {
            this->activated_ = true;
        }
    }
}

void
ReservoirCouplingSlave::
receiveMasterGroupNamesFromMasterProcess_() {
    std::size_t size;
    std::vector<char> group_names;
    if (this->comm_.rank() == 0) {
        auto MPI_SIZE_T_TYPE = Dune::MPITraits<std::size_t>::getType();
        // NOTE: See comment about error handling at the top of this file.
        MPI_Recv(
            &size,
            /*count=*/1,
            /*datatype=*/MPI_SIZE_T_TYPE,
            /*source_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::MasterGroupNamesSize),
            this->slave_master_comm_,
            MPI_STATUS_IGNORE
        );
        OpmLog::info("Received master group names size from master process rank 0");
        group_names.resize(size);
        MPI_Recv(
            group_names.data(),
            /*count=*/size,
            /*datatype=*/MPI_CHAR,
            /*source_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::MasterGroupNames),
            this->slave_master_comm_,
            MPI_STATUS_IGNORE
        );
        OpmLog::info("Received master group names from master process rank 0");
    }
    this->comm_.broadcast(&size, /*count=*/1, /*emitter_rank=*/0);
    if (this->comm_.rank() != 0) {
        group_names.resize(size);
    }
    this->comm_.broadcast(group_names.data(), /*count=*/size, /*emitter_rank=*/0);
    this->saveMasterGroupNamesAsMap_(group_names);
    this->checkGrupSlavGroupNames_();
}

void
ReservoirCouplingSlave::
receiveSlaveNameFromMasterProcess_() {
    std::size_t size;
    std::string slave_name;
    if (this->comm_.rank() == 0) {
        auto MPI_SIZE_T_TYPE = Dune::MPITraits<std::size_t>::getType();
        // NOTE: See comment about error handling at the top of this file.
        MPI_Recv(
            &size,
            /*count=*/1,
            /*datatype=*/MPI_SIZE_T_TYPE,
            /*source_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::SlaveNameSize),
            this->slave_master_comm_,
            MPI_STATUS_IGNORE
        );
        OpmLog::info("Received slave name size from master process rank 0");
        slave_name.resize(size+1); // +1 for the null terminator
        MPI_Recv(
            slave_name.data(),
            /*count=*/size,
            /*datatype=*/MPI_CHAR,
            /*source_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::SlaveName),
            this->slave_master_comm_,
            MPI_STATUS_IGNORE
        );
        slave_name[size] = '\0';  // Add null terminator
        OpmLog::info("Received slave name from master process rank 0");
    }
    this->comm_.broadcast(&size, /*count=*/1, /*emitter_rank=*/0);
    if (this->comm_.rank() != 0) {
        slave_name.resize(size+1); // +1 for the null terminator
    }
    this->comm_.broadcast(slave_name.data(), /*count=*/size+1, /*emitter_rank=*/0);
    OpmLog::info(fmt::format("Received slave name: {}", slave_name));
    this->slave_name_ = slave_name;
}

void
ReservoirCouplingSlave::
saveMasterGroupNamesAsMap_(const std::vector<char>& group_names) {
    // Deserialize the group names vector into a map of slavegroup names -> mastergroup names
    auto total_size = group_names.size();
    std::size_t offset = 0;
    while (offset < total_size) {
        std::string master_group{group_names.data() + offset};
        offset += master_group.size() + 1;
        assert(offset < total_size);
        std::string slave_group{group_names.data() + offset};
        offset += slave_group.size() + 1;
        this->slave_to_master_group_map_[slave_group] = master_group;
    }
}

void
ReservoirCouplingSlave::
sendActivationDateToMasterProcess_() const
{
    if (this->comm_.rank() == 0) {
        // NOTE: The master process will use this date to check that no slave starts before the master
        double activation_date = this->getGrupSlavActivationDate_();
        // NOTE: See comment about error handling at the top of this file.
        MPI_Send(
            &activation_date,
            /*count=*/1,
            /*datatype=*/MPI_DOUBLE,
            /*dest_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::SlaveActivationDate),
            this->slave_master_comm_
        );
        OpmLog::info("Sent simulation activation date to master process from rank 0");
   }
}

void
ReservoirCouplingSlave::
sendSimulationStartDateToMasterProcess_() const
{
    if (this->comm_.rank() == 0) {
        // NOTE: The master process needs the s
        double start_date = this->schedule_.getStartTime();
        // NOTE: See comment about error handling at the top of this file.
        MPI_Send(
            &start_date,
            /*count=*/1,
            /*datatype=*/MPI_DOUBLE,
            /*dest_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::SlaveSimulationStartDate),
            this->slave_master_comm_
        );
        OpmLog::info("Sent simulation start date to master process from rank 0");
   }
}

} // namespace Opm
