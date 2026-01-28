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
#include <opm/simulators/flow/rescoup/ReservoirCoupling.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingErrorMacros.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingMpiTraits.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingSlaveReportStep.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingSlave.hpp>

#include <opm/input/eclipse/Schedule/ResCoup/ReservoirCouplingInfo.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/MasterGroup.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/Slaves.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <dune/common/parallel/mpitraits.hh>

#include <vector>
#include <fmt/format.h>

namespace Opm {
// NOTE: All slave-master communicators have set a custom error handler, which eventually
//   will call MPI_Abort() so there is no need to check the return value of any MPI_Recv()
//   or MPI_Send() calls below.

template <class Scalar>
ReservoirCouplingSlave<Scalar>::
ReservoirCouplingSlave(
    const Parallel::Communication &comm,
    const Schedule &schedule,
    const SimulatorTimer &timer
) :
    comm_{comm},
    schedule_{schedule},
    timer_{timer},
    logger_{comm}
{
    this->slave_master_comm_ = MPI_COMM_NULL;
    MPI_Comm_get_parent(&this->slave_master_comm_);
    if (this->slave_master_comm_ == MPI_COMM_NULL) {
        RCOUP_LOG_THROW(std::runtime_error, "Slave process is not spawned by a master process");
    }
    // NOTE: By installing a custom error handler for all slave-master communicators, which
    //   eventually will call MPI_Abort(), there is no need to check the return value of any
    //   MPI_Recv() or MPI_Send() calls as errors will be caught by the error handler.
    ReservoirCoupling::setErrhandler(this->slave_master_comm_, /*is_master=*/false);
}

template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
initTimeStepping()
{
    assert(!this->report_step_data_);
    this->report_step_data_ = std::make_unique<ReservoirCouplingSlaveReportStep<Scalar>>(*this);
}

template <class Scalar>
bool
ReservoirCouplingSlave<Scalar>::
isFirstSubstepOfSyncTimestep() const
{
    assert(this->report_step_data_);
    return this->report_step_data_->isFirstSubstepOfSyncTimestep();
}

// NOTE: It is not legal for a slave to activate before the master has activated. This problem
//       will be caught by the master when it receives the slave activation date. See:
//       ReservoirCouplingSpawnSlaves::receiveActivationDateFromSlaves_()
template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
maybeActivate(int report_step) {
    if (!this->activated()) {
        bool do_activate = false;
        if (this->historyMatchingMode_()) {
            do_activate = true;
        }
        else {
            const auto& rescoup = this->schedule_[report_step].rescoup();
            if (rescoup.grupSlavCount() > 0) {
                do_activate = true;
            }
        }
        if (do_activate) {
            this->activated_ = true;
            this->sendActivationHandshakeToMasterProcess_();
        }
    }
}

template <class Scalar>
bool
ReservoirCouplingSlave<Scalar>::
maybeReceiveTerminateSignalFromMaster()
{
    // Blocking receive for terminate signal from master.
    // This is called at the start of each timestep iteration in the slave's substep loop.
    // Master sends 0 (continue) or 1 (terminate) at each iteration.

    int terminate_signal = 0;
    if (this->comm_.rank() == 0) {
        MPI_Recv(
            &terminate_signal,
            /*count=*/1,
            /*datatype=*/MPI_INT,
            /*source_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::SlaveProcessTermination),
            this->slave_master_comm_,
            MPI_STATUS_IGNORE
        );
    }
    this->comm_.broadcast(&terminate_signal, /*count=*/1, /*emitter_rank=*/0);

    if (terminate_signal != 0) {
        // Terminate signal received - disconnect and return true
        this->logger_.info("Received terminate signal from master process");
        MPI_Comm_disconnect(&this->slave_master_comm_);
        this->terminated_ = true;
        this->logger_.info("Disconnected intercommunicator with master process");
        return true;
    }
    return false;
}

template <class Scalar>
double
ReservoirCouplingSlave<Scalar>::
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
        this->logger_.info(
            fmt::format("Slave rank 0 received next timestep {} from master.", timestep)
        );
    }
    this->comm_.broadcast(&timestep, /*count=*/1, /*emitter_rank=*/0);
    this->logger_.info("Broadcasted slave next time step to all ranks");
    return timestep;
}

template <class Scalar>
std::pair<std::size_t, std::size_t>
ReservoirCouplingSlave<Scalar>::
receiveNumGroupTargetsFromMaster() const {
    assert(this->report_step_data_);
    return this->report_step_data_->receiveNumGroupTargetsFromMaster();
}

template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
receiveInjectionGroupTargetsFromMaster(std::size_t num_targets) const
{
    assert(this->report_step_data_);
    this->report_step_data_->receiveInjectionGroupTargetsFromMaster(num_targets);
}

template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
receiveProductionGroupTargetsFromMaster(std::size_t num_targets) const
{
    assert(this->report_step_data_);
    this->report_step_data_->receiveProductionGroupTargetsFromMaster(num_targets);
}

template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
receiveTerminateAndDisconnect()
{
    // Receive terminate signal from master (only on rank 0, then broadcast)
    int terminate_signal = 0;
    if (this->comm_.rank() == 0) {
        // NOTE: See comment about error handling at the top of this file.
        MPI_Recv(
            &terminate_signal,
            /*count=*/1,
            /*datatype=*/MPI_INT,
            /*source_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::SlaveProcessTermination),
            this->slave_master_comm_,
            MPI_STATUS_IGNORE
        );
        this->logger_.info("Received terminate signal from master process");
    }
    this->comm_.broadcast(&terminate_signal, /*count=*/1, /*emitter_rank=*/0);

    // Disconnect the intercommunicator (collective operation - all ranks must participate)
    MPI_Comm_disconnect(&this->slave_master_comm_);
    this->terminated_ = true;
    this->logger_.info("Disconnected intercommunicator with master process");
}

template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
sendAndReceiveInitialData() {
    // Communication order must match master's spawn() method.
    // We receive master group names before sending activation date so that
    // we can detect history matching mode (no master groups) and use start_date as activation.
    this->sendSimulationStartDateToMasterProcess_();
    this->receiveSlaveNameFromMasterProcess_();
    this->receiveMasterGroupNamesFromMasterProcess_();
    this->sendActivationDateToMasterProcess_();
}

template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
sendInjectionDataToMaster(
    const std::vector<ReservoirCoupling::SlaveGroupInjectionData<Scalar>> &injection_data
) const
{
    assert(this->report_step_data_);
    this->report_step_data_->sendInjectionDataToMaster(injection_data);
}

template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
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
        this->logger_.info("Sent next report date to master process from rank 0");
   }
}

template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
sendProductionDataToMaster(
    const std::vector<ReservoirCoupling::SlaveGroupProductionData<Scalar>> &production_data
) const
{
    assert(this->report_step_data_);
    this->report_step_data_->sendProductionDataToMaster(production_data);
}

template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
setFirstSubstepOfSyncTimestep(bool value)
{
    assert(this->report_step_data_);
    this->report_step_data_->setFirstSubstepOfSyncTimestep(value);
}

// ------------------
// Private methods
// ------------------

template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
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
                    RCOUP_LOG_THROW(std::runtime_error,
                              "Reservoir coupling: Failed to find master group name for slave group: "
                              + slave_group_name);
                }
                else {
                    const auto& master_group_name = map_iter->second;
                    if (grup_slav.masterGroupName() != master_group_name) {
                        RCOUP_LOG_THROW(std::runtime_error,
                                  "Reservoir coupling: Inconsistent master group name for slave group: "
                                  + slave_group_name);
                    }
                }
            }
        }
    }
    if (!grup_slav_found) {
        RCOUP_LOG_THROW(std::runtime_error, "Reservoir coupling: Failed to find slave group names: "
                  "No GRUPSLAV keyword found in schedule");
    }
}

template <class Scalar>
std::pair<double, bool>
ReservoirCouplingSlave<Scalar>::
getGrupSlavActivationDateAndCheckHistoryMatchingMode_() const
{
    double start_date = this->schedule_.getStartTime();
    for (std::size_t report_step = 0; report_step < this->schedule_.size(); ++report_step) {
        auto rescoup = this->schedule_[report_step].rescoup();
        if (rescoup.grupSlavCount() > 0) {
            return {start_date + this->schedule_.seconds(report_step), false};
        }
    }
    if (this->numMasterGroups_() == 0) {
        // History matching mode: No GRUPSLAV found, activate at simulation start
        this->logger_.warning("Reservoir coupling: No GRUPSLAV keyword found in schedule. "
        "Assuming slave activation at simulation start (history matching mode).");
        return {start_date, true};
    }
    RCOUP_LOG_THROW(std::runtime_error, "Reservoir coupling: Failed to find slave activation time: "
              "No GRUPSLAV keyword found in schedule");
    return {0.0, false};  // Should never happen
}

template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
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
        this->logger_.info(fmt::format(
            "Received master group names size from master process rank 0: {}", size));
        // size can be 0 for history matching mode (no GRUPMAST on master)
        if (size > 0) {
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
            this->logger_.info("Received master group names from master process rank 0");
        }
    }
    this->comm_.broadcast(&size, /*count=*/1, /*emitter_rank=*/0);
    if (size > 0) {
        if (this->comm_.rank() != 0) {
            group_names.resize(size);
        }
        this->comm_.broadcast(group_names.data(), /*count=*/size, /*emitter_rank=*/0);
        this->saveMasterGroupNamesAsMapAndEstablishOrder_(group_names);
        this->checkGrupSlavGroupNames_();
    }
    // If size == 0, master_group_names_ maps remain empty (history matching mode)
}

template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
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
        this->logger_.info("Received slave name size from master process rank 0");
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
        this->logger_.info("Received slave name from master process rank 0");
    }
    this->comm_.broadcast(&size, /*count=*/1, /*emitter_rank=*/0);
    if (this->comm_.rank() != 0) {
        slave_name.resize(size+1); // +1 for the null terminator
    }
    this->comm_.broadcast(slave_name.data(), /*count=*/size+1, /*emitter_rank=*/0);
    this->logger_.info(fmt::format("Received slave name: {}", slave_name));
    this->slave_name_ = slave_name;
}

template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
saveMasterGroupNamesAsMapAndEstablishOrder_(const std::vector<char>& group_names) {
    // Deserialize the group names vector into a map of slavegroup names -> mastergroup names
    // and establish the order the master process sends us the group data. This will enable the
    // master to send group data without the need to send the group names themselves.
    //
    // Call chain: sendAndReceiveInitialData() -> receiveMasterGroupNamesFromMasterProcess_()
    auto total_size = group_names.size();
    std::size_t offset = 0;
    std::size_t idx = 0;
    while (offset < total_size) {
        std::string master_group{group_names.data() + offset};
        offset += master_group.size() + 1;
        assert(offset < total_size);
        std::string slave_group{group_names.data() + offset};
        offset += slave_group.size() + 1;
        this->slave_to_master_group_map_[slave_group] = master_group;
        this->slave_group_order_[idx] = slave_group;
        idx++;
    }
}

template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
sendActivationDateToMasterProcess_()
{
    auto [activation_date, history_matching_mode] = this->getGrupSlavActivationDateAndCheckHistoryMatchingMode_();
    this->history_matching_mode_ = history_matching_mode;
    if (this->comm_.rank() == 0) {
        // NOTE: The master process will use this date to check that no slave starts before the master
        // NOTE: See comment about error handling at the top of this file.
        MPI_Send(
            &activation_date,
            /*count=*/1,
            /*datatype=*/MPI_DOUBLE,
            /*dest_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::SlaveActivationDate),
            this->slave_master_comm_
        );
        this->logger_.info("Sent simulation activation date to master process from rank 0");
   }
}

template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
sendActivationHandshakeToMasterProcess_() const
{
    if (this->comm_.rank() == 0) {
        std::uint8_t activation_handshake = 1u;
        auto MPI_UINT8_T_TYPE = Dune::MPITraits<std::uint8_t>::getType();
        // NOTE: See comment about error handling at the top of this file.
        MPI_Send(
            &activation_handshake,
            /*count=*/1,
            /*datatype=*/MPI_UINT8_T_TYPE,
            /*dest_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::SlaveActivationHandshake),
            this->slave_master_comm_
        );
        this->logger_.info("Sent simulation activation handshake to master process from rank 0");
    }
    this->comm_.barrier();
}

template <class Scalar>
void
ReservoirCouplingSlave<Scalar>::
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
        this->logger_.info("Sent simulation start date to master process from rank 0");
   }
}

// Explicit template instantiations
template class ReservoirCouplingSlave<double>;
#if FLOW_INSTANTIATE_FLOAT
template class ReservoirCouplingSlave<float>;
#endif

} // namespace Opm
