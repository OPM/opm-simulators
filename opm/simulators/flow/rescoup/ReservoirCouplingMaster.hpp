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
#include <opm/simulators/flow/rescoup/ReservoirCoupling.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingMasterReportStep.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingTimeStepper.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <mpi.h>

#include <filesystem>
#include <vector>

namespace Opm {

template <class Scalar>
class ReservoirCouplingMaster {
public:
    using MessageTag = ReservoirCoupling::MessageTag;
    using Seconds = ReservoirCoupling::Seconds;
    using Potentials = ReservoirCoupling::Potentials<Scalar>;
    using SlaveGroupProductionData = ReservoirCoupling::SlaveGroupProductionData<Scalar>;
    using InjectionGroupTarget = ReservoirCoupling::InjectionGroupTarget<Scalar>;
    using ProductionGroupTarget = ReservoirCoupling::ProductionGroupTarget<Scalar>;

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
    void addSlaveActivationDate(double date) { this->slave_activation_dates_.push_back(date); }
    void addSlaveStartDate(std::time_t date) { this->slave_start_dates_.push_back(date); }
    void clearDeferredLogger() { logger_.clearDeferredLogger(); }
    double getActivationDate() const { return this->activation_date_; }
    int getArgc() const { return this->argc_; }
    char *getArgv(int index) const { return this->argv_[index]; }
    char **getArgv() const { return this->argv_; }
    const Parallel::Communication &getComm() const { return this->comm_; }
    /// @brief Get the master group names associated with a slave reservoir by index.
    ///
    /// This method retrieves the list of master group names that are associated with a
    /// specific slave reservoir identified by its index.
    ///
    /// @param slave_idx The zero-based index of the slave reservoir (must be < numSlaves())
    /// @return A const reference to a vector of master group names for the specified slave
    /// @throws std::runtime_error if slave_idx is out of bounds
    ///
    /// @note Performance: This method uses O(1) direct vector access when possible,
    ///       falling back to O(log n) map lookup for error handling.
    /// @see RescoupTargetCalculator::calculateAndSendTargets() for primary usage context
    const std::vector<std::string>& getMasterGroupNamesForSlave(std::size_t slave_idx) const;
    /// @brief Get the canonical index of the master group for a given slave name and master group name.
    /// The index is used to map slave group data sent from the slaves, like potentials to the corresponding
    /// master group.
    /// @param slave_name The name of the slave reservoir.
    /// @param master_group_name The name of the master group.
    /// @return The canonical index of the master group for the given slave name and master group name.
    std::size_t getMasterGroupCanonicalIdx(
        const std::string &slave_name, const std::string &master_group_name) const;
    Scalar getMasterGroupInjectionRate(const std::string &group_name, ReservoirCoupling::Phase phase, bool res_rates) const;
    Scalar getMasterGroupProductionRate(const std::string &group_name, ReservoirCoupling::Phase phase, bool res_rates) const;
    std::map<std::string, std::string>& getMasterGroupToSlaveNameMap() {
         return this->master_group_slave_names_;
    }
    double getSimulationStartDate() const { return this->schedule_.getStartTime(); }
    double getSlaveActivationDate(int index) const { return this->slave_activation_dates_[index]; }
    const double *getSlaveActivationDates() const { return this->slave_activation_dates_.data(); }
    MPI_Comm getSlaveComm(int index) const { return this->master_slave_comm_[index]; }
    std::map<std::string, std::vector<std::string>> &getSlaveNameToMasterGroupsMap() {
        return this->slave_name_to_master_groups_map_;
    }
    const Potentials& getSlaveGroupPotentials(const std::string &master_group_name);
    int getSlaveIdx(const std::string &slave_name) const;
    const std::string &getSlaveName(int index) const { return this->slave_names_[index]; }
    double getSlaveStartDate(int index) const { return this->slave_start_dates_[index]; }
    const double *getSlaveStartDates() { return this->slave_start_dates_.data(); }
    void initStartOfReportStep(int report_step_idx);
    void initTimeStepping();
    bool isFirstSubstepOfSyncTimestep() const;
    bool isMasterGroup(const std::string &group_name) const;
    ReservoirCoupling::Logger& logger() { return this->logger_; }
    ReservoirCoupling::Logger& logger() const { return this->logger_; }
    void maybeActivate(int report_step);
    void maybeReceiveActivationHandshakeFromSlaves(double current_time);
    double maybeChopSubStep(double suggested_timestep, double current_time) const;
    void maybeSpawnSlaveProcesses(int report_step);
    std::size_t numSlaveGroups(unsigned int index);
    std::size_t numSlaves() const { return this->numSlavesStarted(); }
    std::size_t numSlavesStarted() const;
    void rebuildSlaveIdxToMasterGroupsVector();
    void receiveNextReportDateFromSlaves();
    void receiveProductionDataFromSlaves();
    void receiveInjectionDataFromSlaves();
    void resizeNextReportDates(int size);
    void resizeSlaveActivationDates(int size) { this->slave_activation_dates_.resize(size); }
    void resizeSlaveStartDates(int size) { this->slave_start_dates_.resize(size); }
    const Schedule& schedule() const { return this->schedule_; }
    void sendNextTimeStepToSlaves(double next_time_step) {
        this->time_stepper_->sendNextTimeStepToSlaves(next_time_step);
    }
    void sendInjectionTargetsToSlave(
        std::size_t slave_idx,
        const std::vector<InjectionGroupTarget>& injection_targets
    ) const;
    void sendNumGroupTargetsToSlave(
        std::size_t slave_idx,
        std::size_t num_injection_targets,
        std::size_t num_production_targets
    ) const;
    void sendProductionTargetsToSlave(
        std::size_t slave_idx,
        const std::vector<ProductionGroupTarget>& production_targets
    ) const;
    void setDeferredLogger(DeferredLogger *deferred_logger) {
         this->logger_.setDeferredLogger(deferred_logger);
    }
    void setFirstSubstepOfSyncTimestep(bool value);
    // These are currently only used for unit testing
    void setSlaveActivationDate(int index, double date) { this->slave_activation_dates_[index] = date; }
    void setSlaveNextReportTimeOffset(int index, double offset);
    void setSlaveStartDate(int index, std::time_t date) { this->slave_start_dates_[index] = date; }
    bool slaveIsActivated(int index) const { return this->slave_activation_status_[index] != 0; }
    void updateMasterGroupNameOrderMap(
        const std::string& slave_name, const std::map<std::string, std::size_t>& master_group_map);

    /// @brief Send "don't terminate" signal (value=0) to all active slaves.
    ///
    /// This method is called at the start of each iteration in the master's substep loop,
    /// before receiving the next report date from slaves. The slave waits for this signal
    /// at the start of each iteration - if it receives 0, it continues; if it receives 1
    /// (from sendTerminateAndDisconnect()), it terminates.
    void sendDontTerminateSignalToSlaves();

    /// @brief Send terminate signal to all active slaves and disconnect intercommunicators.
    ///
    /// This method must be called at the end of the simulation to cleanly shut down
    /// the MPI intercommunicators created by MPI_Comm_spawn(). It performs two steps:
    /// 1. Sends a terminate signal to all active slaves (only from rank 0)
    /// 2. Disconnects the intercommunicators (collective operation)
    ///
    /// Both master and slaves must call their respective disconnect methods for
    /// MPI_Comm_disconnect() to complete - it is a collective operation.
    void sendTerminateAndDisconnect();

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
    double activation_date_{0.0};  // The date when SLAVES is encountered in the schedule
    // A mapping from a slave name to the master group name order used when slaves send
    //  potentials to the master process.
    std::map<std::string, std::map<std::string, std::size_t>> master_group_name_order_;
    mutable ReservoirCoupling::Logger logger_;
    // Whether the slave has activated. Unfortunatley, we cannot use std::vector<bool> since
    // it is not supported to get a pointer to the underlying array of bools needed
    // with MPI broadcast().
    std::vector<std::uint8_t> slave_activation_status_;
    // A mapping from master group names to slave names
    std::map<std::string, std::string> master_group_slave_names_;
    // A mapping from slave names to master group names
    // NOTE: The order of the master groups in the vector is important,
    //   as the slaves will communicate the indices of the master groups in
    //   the vector instead of the group names themselves.
    // NOTE: This map is created by ReservoirCouplingSpawnSlaves.cpp
    std::map<std::string, std::vector<std::string>> slave_name_to_master_groups_map_;
    // Direct index-based lookup for performance optimization (O(1) instead of O(log n))
    // This vector is populated in parallel with slave_name_to_master_groups_map_
    // and maintains the same ordering as slave_names_ vector for consistent indexing
    std::vector<std::vector<std::string>> slave_idx_to_master_groups_;
    // Stores data that changes for a single report step or for timesteps within a report step.
    std::unique_ptr<ReservoirCouplingMasterReportStep<Scalar>> report_step_data_{nullptr};
    // Handles time stepping for the master and slaves
    std::unique_ptr<ReservoirCouplingTimeStepper<Scalar>> time_stepper_{nullptr};
};

} // namespace Opm

#endif // OPM_RESERVOIR_COUPLING_MASTER_HPP
