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

#include <opm/simulators/flow/ReservoirCoupling.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>
#include <opm/simulators/timestepping/SimulatorTimer.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <mpi.h>

#include <vector>

namespace Opm {

class ReservoirCouplingSlave {
public:
    using MessageTag = ReservoirCoupling::MessageTag;
    using Potentials = ReservoirCoupling::Potentials;

    ReservoirCouplingSlave(
        const Parallel::Communication &comm, const Schedule &schedule, const SimulatorTimer &timer
    );
    bool activated() const { return activated_; }
    void clearDeferredLogger() { logger_.clearDeferredLogger(); }
    const Parallel::Communication& getComm() const { return comm_; }
    const std::map<std::string, std::string>& getSlaveToMasterGroupNameMap() const {
        return slave_to_master_group_map_; }
    void maybeActivate(int report_step);
    double receiveNextTimeStepFromMaster();
    void sendAndReceiveInitialData();
    void sendNextReportDateToMasterProcess() const;
    void sendPotentialsToMaster(const std::vector<Potentials> &potentials) const;
    void setDeferredLogger(DeferredLogger *deferred_logger) {
        this->logger_.setDeferredLogger(deferred_logger);
    }

private:
    void checkGrupSlavGroupNames_();
    double getGrupSlavActivationDate_() const;
    void receiveMasterGroupNamesFromMasterProcess_();
    void receiveSlaveNameFromMasterProcess_();
    void saveMasterGroupNamesAsMap_(const std::vector<char>& group_names);
    void sendActivationDateToMasterProcess_() const;
    void sendSimulationStartDateToMasterProcess_() const;

    const Parallel::Communication &comm_;
    const Schedule& schedule_;
    const SimulatorTimer &timer_;
    // MPI parent communicator for a slave process
    MPI_Comm slave_master_comm_{MPI_COMM_NULL};
    std::map<std::string, std::string> slave_to_master_group_map_;
    bool activated_{false};
    std::string slave_name_;  // This is the slave name as defined in the master process
    ReservoirCoupling::Logger logger_;
};

} // namespace Opm
#endif // OPM_RESERVOIR_COUPLING_SLAVE_HPP
