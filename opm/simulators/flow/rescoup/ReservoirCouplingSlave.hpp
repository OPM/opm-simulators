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
    using ProductionGroupTarget = ReservoirCoupling::ProductionGroupTarget<Scalar>;

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
    void initTimeStepping();
    ReservoirCoupling::Logger& logger() { return this->logger_; }
    ReservoirCoupling::Logger& logger() const { return this->logger_; }
    void maybeActivate(int report_step);
    std::size_t numSlaveGroups() const { return this->slave_group_order_.size(); }
    double receiveNextTimeStepFromMaster();
    std::pair<std::size_t, std::size_t> receiveNumGroupTargetsFromMaster() const;
    void receiveInjectionGroupTargetsFromMaster(std::size_t num_targets) const;
    void receiveProductionGroupTargetsFromMaster(std::size_t num_targets) const;
    void sendAndReceiveInitialData();
    void sendInjectionDataToMaster(const std::vector<SlaveGroupInjectionData> &injection_data) const;
    void sendNextReportDateToMasterProcess() const;
    void sendProductionDataToMaster(const std::vector<SlaveGroupProductionData> &production_data) const;
    void setDeferredLogger(DeferredLogger *deferred_logger) {
        this->logger_.setDeferredLogger(deferred_logger);
    }
    const std::string& slaveGroupIdxToGroupName(std::size_t group_idx) const {
        return this->slave_group_order_.at(group_idx);
    }

private:
    void checkGrupSlavGroupNames_();
    double getGrupSlavActivationDate_() const;
    void receiveMasterGroupNamesFromMasterProcess_();
    void receiveSlaveNameFromMasterProcess_();
    void saveMasterGroupNamesAsMapAndEstablishOrder_(const std::vector<char>& group_names);
    void sendActivationDateToMasterProcess_() const;
    void sendActivationHandshakeToMasterProcess_() const;
    void sendSimulationStartDateToMasterProcess_() const;

    const Parallel::Communication &comm_;
    const Schedule& schedule_;
    const SimulatorTimer &timer_;
    // MPI parent communicator for a slave process
    MPI_Comm slave_master_comm_{MPI_COMM_NULL};
    std::map<std::string, std::string> slave_to_master_group_map_;
    bool activated_{false};
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
