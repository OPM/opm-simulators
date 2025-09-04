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
    using Potentials = ReservoirCoupling::Potentials;
    ReservoirCouplingMaster(
        const Parallel::Communication &comm,
        const Schedule &schedule,
        int argc, char **argv
    );

    bool activated() { return this->activated_; }
    void addSlaveCommunicator(MPI_Comm comm) {
         this->master_slave_comm_.push_back(comm);
    }
    void addSlaveName(const std::string &name) { this->slave_names_.push_back(name); }
    void addSlaveNextReportTimeOffset(double offset) {
         this->slave_next_report_time_offsets_.push_back(offset);
    }
    void addSlaveActivationDate(double date) { this->slave_activation_dates_.push_back(date); }
    void addSlaveStartDate(std::time_t date) { this->slave_start_dates_.push_back(date); }
    void clearDeferredLogger() { logger_.clearDeferredLogger(); }
    double getActivationDate() const { return this->activation_date_; }
    int getArgc() const { return this->argc_; }
    char *getArgv(int index) const { return this->argv_[index]; }
    char **getArgv() const { return this->argv_; }
    const Parallel::Communication &getComm() const { return this->comm_; }
    /// @brief Get the index of the master group potential for a given slave name and master group name.
    /// The index is used to map the slave group potentials to the master group potentials.
    /// @param slave_name The name of the slave reservoir.
    /// @param master_group_name The name of the master group.
    /// @return The index of the master group potential for the given slave name and master group name.
    std::size_t getMasterGroupPotIdx(
        const std::string &slave_name, const std::string &master_group_name) const;
    std::map<std::string, std::string>& getMasterGroupToSlaveNameMap() {
         return this->master_group_slave_names_;
    }
    double getSlaveActivationDate(int index) const { return this->slave_activation_dates_[index]; }
    const double *getSlaveActivationDates() const { return this->slave_activation_dates_.data(); }
    double getSimulationStartDate() const { return this->schedule_.getStartTime(); }
    MPI_Comm getSlaveComm(int index) const { return this->master_slave_comm_[index]; }
    const Potentials& getSlaveGroupPotentials(const std::string &master_group_name);
    const std::string &getSlaveName(int index) const { return this->slave_names_[index]; }
    const double *getSlaveStartDates() { return this->slave_start_dates_.data(); }
    bool isMasterGroup(const std::string &group_name) const;
    void maybeActivate(int report_step);
    void maybeReceiveActivationHandshakeFromSlaves(double current_time);
    double maybeChopSubStep(double suggested_timestep, double current_time) const;
    void maybeSpawnSlaveProcesses(int report_step);
    std::size_t numSlaveGroups(unsigned int index);
    std::size_t numSlavesStarted() const;
    void receiveNextReportDateFromSlaves();
    void receivePotentialsFromSlaves();
    void resizeSlaveActivationDates(int size) { this->slave_activation_dates_.resize(size); }
    void resizeSlaveStartDates(int size) { this->slave_start_dates_.resize(size); }
    void resizeNextReportDates(int size) { this->slave_next_report_time_offsets_.resize(size); }
    void sendNextTimeStepToSlaves(double next_time_step);
    void setDeferredLogger(DeferredLogger *deferred_logger) {
         this->logger_.setDeferredLogger(deferred_logger);
    }
    // These are currently only used for unit testing
    void setSlaveActivationDate(int index, double date) { this->slave_activation_dates_[index] = date; }
    void setSlaveStartDate(int index, std::time_t date) { this->slave_start_dates_[index] = date; }
    void setSlaveNextReportTimeOffset(int index, double offset) {
         this->slave_next_report_time_offsets_[index] = offset;
    }
    bool slaveIsActivated(int index) const { return this->slave_activation_status_[index] != 0; }
    void updateMasterGroupNameOrderMap(
        const std::string& slave_name, const std::map<std::string, std::size_t>& master_group_map);

private:
    double getMasterActivationDate_() const;

    const Parallel::Communication &comm_;
    const Schedule& schedule_;
    int argc_;
    char **argv_;
    // Whether the master process has activated the reservoir coupling
    bool activated_{false};
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
    // The activation dates are in whole seconds since the epoch.
    std::vector<double> slave_activation_dates_;
    // Elapsed time from the beginning of the simulation
    std::vector<double> slave_next_report_time_offsets_;
    double activation_date_{0.0};  // The date when SLAVES is encountered in the schedule
    // A mapping from a slave name to the master group name order used when slaves send
    //  potentials to the master process.
    std::map<std::string, std::map<std::string, std::size_t>> master_group_name_order_;
    ReservoirCoupling::Logger logger_;
    // Whether the slave has activated. Unfortunatley, we cannot use std::vector<bool> since
    // it is not supported to get a pointer to the underlying array of bools needed
    // with MPI broadcast().
    std::vector<std::uint8_t> slave_activation_status_;
    // Potentials for oil, gas, and water rates for each slave group
    std::map<std::string, std::vector<Potentials>> slave_group_potentials_;
    // A mapping from master group names to slave names
    std::map<std::string, std::string> master_group_slave_names_;
};

} // namespace Opm
#endif // OPM_RESERVOIR_COUPLING_MASTER_HPP
