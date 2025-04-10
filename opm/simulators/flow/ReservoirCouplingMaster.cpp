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

ReservoirCouplingMaster::
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

void
ReservoirCouplingMaster::
maybeSpawnSlaveProcesses(int report_step)
{
    if (this->numSlavesStarted() > 0) {  // We have already spawned the slave processes
        return;
    }
    const auto& rescoup = this->schedule_[report_step].rescoup();
    auto slave_count = rescoup.slaveCount();
    auto master_group_count = rescoup.masterGroupCount();
    if (slave_count > 0 && master_group_count > 0) {
        ReservoirCouplingSpawnSlaves spawn_slaves{*this, rescoup};
        spawn_slaves.spawn();
    }
}


double
ReservoirCouplingMaster::
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
        double slave_next_report_date{this->slave_next_report_time_offsets_[i] + slave_start_date};
        if (Seconds::compare_gt_or_eq(slave_start_date, step_end_date)) {
            // The slave process has not started yet, and will not start during this timestep
            continue;
        }
        double slave_elapsed_time;
        if (Seconds::compare_lt_or_eq(slave_start_date,step_start_date)) {
            // The slave process has already started, and will continue during this timestep
            if (Seconds::compare_gt(slave_next_report_date, step_end_date)) {
                // The slave process will not report during this timestep
                continue;
            }
            // The slave process will report during this timestep
            slave_elapsed_time = slave_next_report_date - step_start_date;
        }
        else {
            // The slave process will start during the timestep, but not at the beginning
            slave_elapsed_time = slave_start_date - step_start_date;
        }
        suggested_timestep = slave_elapsed_time;
        step_end_date = step_start_date + suggested_timestep;
    }
    return suggested_timestep;
}

void
ReservoirCouplingMaster::
sendNextTimeStepToSlaves(double timestep)
{
    if (this->comm_.rank() == 0) {
        for (unsigned int i = 0; i < this->master_slave_comm_.size(); i++) {
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


void
ReservoirCouplingMaster::
receiveNextReportDateFromSlaves()
{
    auto num_slaves = this->numSlavesStarted();
    OpmLog::info("Receiving next report dates from slave processes");
    if (this->comm_.rank() == 0) {
        for (unsigned int i = 0; i < num_slaves; i++) {
            double slave_next_report_time_offset; // Elapsed time from the beginning of the simulation
            // NOTE: All slave-master communicators have set a custom error handler, which eventually
            //   will call MPI_Abort() so there is no need to check the return value of any MPI_Recv()
            //   or MPI_Send() calls.
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


std::size_t
ReservoirCouplingMaster::
numSlavesStarted() const
{
    return this->slave_names_.size();
}

// ------------------
// Private methods
// ------------------

double
ReservoirCouplingMaster::
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


} // namespace Opm

