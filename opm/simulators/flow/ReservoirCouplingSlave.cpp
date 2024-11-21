/*
  Copyright 2024 Equinor AS

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
#include <opm/simulators/flow/ReservoirCouplingSlave.hpp>

#include <opm/input/eclipse/Schedule/ResCoup/ReservoirCouplingInfo.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/MasterGroup.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/Slaves.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <vector>
#include <fmt/format.h>

namespace Opm {

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
    this->slave_master_comm_ = MPI_Comm_Ptr(new MPI_Comm(MPI_COMM_NULL));
    MPI_Comm_get_parent(this->slave_master_comm_.get());
    if (*(this->slave_master_comm_) == MPI_COMM_NULL) {
        OPM_THROW(std::runtime_error, "Slave process is not spawned by a master process");
    }
    ReservoirCoupling::setErrhandler(*this->slave_master_comm_, /*is_master=*/false);
}

double
ReservoirCouplingSlave::
receiveNextTimeStepFromMaster() {
    double timestep;
    if (this->comm_.rank() == 0) {
        int result = MPI_Recv(
            &timestep,
            /*count=*/1,
            /*datatype=*/MPI_DOUBLE,
            /*source_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::SlaveNextTimeStep),
            *this->slave_master_comm_,
            MPI_STATUS_IGNORE
        );
        if (result != MPI_SUCCESS) {
            OPM_THROW(std::runtime_error, "Failed to receive next time step from master");
        }
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
receiveMasterGroupNamesFromMasterProcess() {
    std::size_t size;
    std::vector<char> group_names;
    if (this->comm_.rank() == 0) {
        MPI_Aint asize = 0;
        int result = MPI_Recv(
            &asize,
            /*count=*/1,
            /*datatype=*/MPI_AINT,
            /*source_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::MasterGroupNamesSize),
            *this->slave_master_comm_,
            MPI_STATUS_IGNORE
        );
        OpmLog::info("Received master group names size from master process rank 0");
        if (result != MPI_SUCCESS) {
            OPM_THROW(std::runtime_error,
                      "Failed to receive master group names (size) from master process");
        }
        // NOTE: MPI_Aint and std::size_t should be compatible on most systems, but we will
        //     cast it to std::size_t to avoid any potential issues
        size = static_cast<std::size_t>(asize);
        group_names.resize(size);
        int result2 = MPI_Recv(
            group_names.data(),
            /*count=*/size,
            /*datatype=*/MPI_CHAR,
            /*source_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::MasterGroupNames),
            *this->slave_master_comm_,
            MPI_STATUS_IGNORE
        );
        if (result2 != MPI_SUCCESS) {
            OPM_THROW(std::runtime_error,
                      "Failed to receive master group names from master process");
        }
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
sendNextReportDateToMasterProcess() const
{
    if (this->comm_.rank() == 0) {
        double elapsed_time = this->timer_.simulationTimeElapsed();
        double current_step_length = this->timer_.currentStepLength();
        // NOTE: This is an offset in seconds from the start date, so it will be 0 if the next report
        //      would be the start date. In general, it should be a positive number.
        double next_report_time_offset = elapsed_time + current_step_length;
        MPI_Send(
            &next_report_time_offset,
            /*count=*/1,
            /*datatype=*/MPI_DOUBLE,
            /*dest_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::SlaveNextReportDate),
            *this->slave_master_comm_
        );
        OpmLog::info("Sent next report date to master process from rank 0");
   }
}

void
ReservoirCouplingSlave::
sendActivationDateToMasterProcess() const
{
    if (this->comm_.rank() == 0) {
        // NOTE: The master process needs the s
        double activation_date = this->getGrupSlavActivationDate_();
        MPI_Send(
            &activation_date,
            /*count=*/1,
            /*datatype=*/MPI_DOUBLE,
            /*dest_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::SlaveActivationDate),
            *this->slave_master_comm_
        );
        OpmLog::info("Sent simulation start date to master process from rank 0");
   }
}

void
ReservoirCouplingSlave::
sendSimulationStartDateToMasterProcess() const
{
    if (this->comm_.rank() == 0) {
        // NOTE: The master process needs the s
        double start_date = this->schedule_.getStartTime();
        MPI_Send(
            &start_date,
            /*count=*/1,
            /*datatype=*/MPI_DOUBLE,
            /*dest_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::SlaveSimulationStartDate),
            *this->slave_master_comm_
        );
        OpmLog::info("Sent simulation start date to master process from rank 0");
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
                auto master_group_name_it = this->master_group_names_.find(slave_group_name);
                if (master_group_name_it == this->master_group_names_.end()) {
                    OPM_THROW(std::runtime_error,
                              "Reservoir coupling: Failed to find master group name for slave group: "
                              + slave_group_name);
                }
                else {
                    auto master_group_name = master_group_name_it->second;
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
        this->master_group_names_[slave_group] = master_group;
    }
}


} // namespace Opm
