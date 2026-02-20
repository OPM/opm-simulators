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

#ifndef OPM_RESERVOIR_COUPLING_SLAVE_HPP
#define OPM_RESERVOIR_COUPLING_SLAVE_HPP

#include <opm/simulators/flow/rescoup/ReservoirCoupling.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingSlaveReportStep.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>
#include <opm/simulators/timestepping/SimulatorTimer.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <mpi.h>

#include <vector>

namespace Opm {

template <class Scalar>
class ReservoirCouplingSlaveReportStep;

template <class Scalar>
class ReservoirCouplingSlave {
public:
    using MessageTag = ReservoirCoupling::MessageTag;
    using Potentials = ReservoirCoupling::Potentials<Scalar>;
    using SlaveGroupInjectionData = ReservoirCoupling::SlaveGroupInjectionData<Scalar>;
    using SlaveGroupProductionData = ReservoirCoupling::SlaveGroupProductionData<Scalar>;
    using InjectionGroupTarget = ReservoirCoupling::InjectionGroupTarget<Scalar>;
    using ProductionGroupConstraints = ReservoirCoupling::ProductionGroupConstraints<Scalar>;

    ReservoirCouplingSlave(
        const Parallel::Communication &comm, const Schedule &schedule, const SimulatorTimer &timer
    );
    bool activated() const { return activated_; }
    void clearDeferredLogger() { logger_.clearDeferredLogger(); }
    const Parallel::Communication& getComm() const { return comm_; }
    MPI_Comm getMasterComm() const { return slave_master_comm_; }
    const std::string& getSlaveName() const { return slave_name_; }
    const std::map<std::string, std::string>& getSlaveToMasterGroupNameMap() const {
        return slave_to_master_group_map_; }
    /// @brief Check if a master-imposed injection target exists for a group and phase
    /// @details Delegates to ReservoirCouplingSlaveReportStep
    bool hasMasterInjectionTarget(const std::string& gname, Phase phase) const;

    /// @brief Check if a master-imposed production target exists for a group
    /// @details Delegates to ReservoirCouplingSlaveReportStep
    bool hasMasterProductionTarget(const std::string& gname) const;
    void initTimeStepping();
    bool isFirstSubstepOfSyncTimestep() const;
    bool isSlaveGroup(const std::string& group_name) const;
    ReservoirCoupling::Logger& logger() { return this->logger_; }
    ReservoirCoupling::Logger& logger() const { return this->logger_; }
    /// @brief Get the master-imposed injection target and control mode for a group and phase
    /// @details Delegates to ReservoirCouplingSlaveReportStep
    std::pair<Scalar, Group::InjectionCMode> masterInjectionTarget(const std::string& gname, Phase phase) const;

    /// @brief Get the master-imposed production target and control mode for a group
    /// @details Delegates to ReservoirCouplingSlaveReportStep
    std::pair<Scalar, Group::ProductionCMode> masterProductionTarget(const std::string& gname) const;
    void maybeActivate(int report_step);
    std::size_t numSlaveGroups() const { return this->slave_group_order_.size(); }
    double receiveNextTimeStepFromMaster();
    std::pair<std::size_t, std::size_t> receiveNumGroupConstraintsFromMaster() const;
    /// @brief Receive injection group targets from master and store them locally
    /// @details Delegates to ReservoirCouplingSlaveReportStep
    void receiveInjectionGroupTargetsFromMaster(std::size_t num_targets);

    /// @brief Receive production group constraints from master and store them locally
    /// @details Delegates to ReservoirCouplingSlaveReportStep
    void receiveProductionGroupConstraintsFromMaster(std::size_t num_targets);

    void sendAndReceiveInitialData();
    void sendInjectionDataToMaster(const std::vector<SlaveGroupInjectionData> &injection_data) const;
    void sendNextReportDateToMasterProcess() const;
    void sendProductionDataToMaster(const std::vector<SlaveGroupProductionData> &production_data) const;
    void setDeferredLogger(DeferredLogger *deferred_logger) {
        this->logger_.setDeferredLogger(deferred_logger);
    }
    void setFirstSubstepOfSyncTimestep(bool value);
    const std::string& slaveGroupIdxToGroupName(std::size_t group_idx) const {
        return this->slave_group_order_.at(group_idx);
    }
    bool terminated() const { return this->terminated_; }

    /// @brief Blocking receive for terminate/continue signal from master.
    ///
    /// This method is called at the start of each timestep iteration in the slave's
    /// substep loop. Master sends 0 (continue) or 1 (terminate) at each iteration.
    /// If terminate signal (value != 0) is received, disconnects the intercommunicator,
    /// sets the terminated flag, and returns true.
    ///
    /// @return true if terminate signal received and disconnect completed, false to continue.
    bool maybeReceiveTerminateSignalFromMaster();

    /// @brief Receive terminate signal from master and disconnect the intercommunicator.
    ///
    /// This method must be called at the end of the simulation to cleanly shut down
    /// the MPI intercommunicator created when the slave was spawned. It performs two steps:
    /// 1. Receives a terminate signal from master (only on rank 0, then broadcast)
    /// 2. Disconnects the intercommunicator (collective operation)
    ///
    /// Both master and slaves must call their respective disconnect methods for
    /// MPI_Comm_disconnect() to complete - it is a collective operation.
    void receiveTerminateAndDisconnect();

private:
    void checkGrupSlavGroupNames_();
    std::pair<double, bool> getGrupSlavActivationDateAndCheckHistoryMatchingMode_() const;
    bool historyMatchingMode_() const { return this->history_matching_mode_; }
    std::size_t numMasterGroups_() const { return this->slave_to_master_group_map_.size(); }
    void receiveMasterGroupNamesFromMasterProcess_();
    void receiveSlaveNameFromMasterProcess_();
    void saveMasterGroupNamesAsMapAndEstablishOrder_(const std::vector<char>& group_names);
    void sendActivationDateToMasterProcess_();
    void sendActivationHandshakeToMasterProcess_() const;
    void sendSimulationStartDateToMasterProcess_() const;

    const Parallel::Communication &comm_;
    const Schedule& schedule_;
    const SimulatorTimer &timer_;
    // MPI parent communicator for a slave process
    MPI_Comm slave_master_comm_{MPI_COMM_NULL};
    std::map<std::string, std::string> slave_to_master_group_map_;
    bool activated_{false};
    // True if the slave was terminated by the master process
    bool terminated_{false};
    // True if no GRUPMAST keyword in the master schedule and no GRUPSLAV keyword in the slave schedule
    bool history_matching_mode_{false};
    std::string slave_name_;  // This is the slave name as defined in the master process
    mutable ReservoirCoupling::Logger logger_;
    // Order of the slave groups. A mapping from slave group index to slave group name.
    // The indices are determined by the order the master process sends us the group names, see
    // receiveMasterGroupNamesFromMasterProcess_()
    // Later, the master process will send us group name indices, and not the group names themselves,
    // so we use this mapping to recover the slave group names from the indices.
    std::map<std::size_t, std::string> slave_group_order_;
    // Stores data that changes for a single report step or for timesteps within a report step.
    std::unique_ptr<ReservoirCouplingSlaveReportStep<Scalar>> report_step_data_{nullptr};
};

} // namespace Opm

#endif // OPM_RESERVOIR_COUPLING_SLAVE_HPP
