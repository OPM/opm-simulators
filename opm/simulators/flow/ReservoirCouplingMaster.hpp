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

#ifndef OPM_RESERVOIR_COUPLING_MASTER_HPP
#define OPM_RESERVOIR_COUPLING_MASTER_HPP

#include <opm/simulators/utils/ParallelCommunication.hpp>
#include <opm/simulators/flow/ReservoirCoupling.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <mpi.h>

#include <filesystem>
#include <vector>

namespace Opm {

class ReservoirCouplingMaster {
public:
    using MessageTag = ReservoirCoupling::MessageTag;
    using Seconds = ReservoirCoupling::Seconds;
    ReservoirCouplingMaster(
        const Parallel::Communication &comm,
        const Schedule &schedule,
        int argc, char **argv
    );

    bool activated() { return this->numSlavesStarted() > 0; }
    void addSlaveCommunicator(MPI_Comm comm) {
         this->master_slave_comm_.push_back(comm);
    }
    void addSlaveName(const std::string &name) { this->slave_names_.push_back(name); }
    void addSlaveNextReportTimeOffset(double offset) {
         this->slave_next_report_time_offsets_.push_back(offset);
    }
    void addSlaveStartDate(std::time_t date) { this->slave_start_dates_.push_back(date); }
    double getActivationDate() const { return this->activation_date_; }
    int getArgc() const { return this->argc_; }
    char *getArgv(int index) const { return this->argv_[index]; }
    char **getArgv() const { return this->argv_; }
    const Parallel::Communication &getComm() const { return this->comm_; }
    double getSimulationStartDate() const { return this->schedule_.getStartTime(); }
    MPI_Comm getSlaveComm(int index) const { return this->master_slave_comm_[index]; }
    const std::string &getSlaveName(int index) const { return this->slave_names_[index]; }
    const double *getSlaveStartDates() { return this->slave_start_dates_.data(); }
    double maybeChopSubStep(double suggested_timestep, double current_time) const;
    void maybeSpawnSlaveProcesses(int report_step);
    std::size_t numSlavesStarted() const;
    void receiveNextReportDateFromSlaves();
    void resizeSlaveStartDates(int size) { this->slave_start_dates_.resize(size); }
    void resizeNextReportDates(int size) { this->slave_next_report_time_offsets_.resize(size); }
    void sendNextTimeStepToSlaves(double next_time_step);
    // These are currently only used for unit testing
    void setSlaveStartDate(int index, std::time_t date) { this->slave_start_dates_[index] = date; }
    void setSlaveNextReportTimeOffset(int index, double offset) {
         this->slave_next_report_time_offsets_[index] = offset;
    }

private:
    double getMasterActivationDate_() const;

    const Parallel::Communication &comm_;
    const Schedule& schedule_;
    int argc_;
    char **argv_;
    // NOTE: MPI_Comm is just an integer handle, so we can just copy it into the vector
    std::vector<MPI_Comm> master_slave_comm_; // MPI communicators for the slave processes
    std::vector<std::string> slave_names_;
    // The start dates are in whole seconds since the epoch. We use a double to store the value
    // since both schedule_.getStartTime() and schedule_.stepLength(report_step) returns
    // a double value representing whole seconds.
    // However, note that schedule_[report_step].start_time() returns a time_point
    // which can include milliseconds. The double values are also convenient when we need to
    // to add fractions of seconds for sub steps to the start date.
    std::vector<double> slave_start_dates_;
    // Elapsed time from the beginning of the simulation
    std::vector<double> slave_next_report_time_offsets_;
    double activation_date_{0.0};  // The date when SLAVES is encountered in the schedule
};

} // namespace Opm
#endif // OPM_RESERVOIR_COUPLING_MASTER_HPP
