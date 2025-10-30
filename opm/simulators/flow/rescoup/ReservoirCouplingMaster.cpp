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

#include <opm/simulators/flow/rescoup/ReservoirCouplingMaster.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingMasterReportStep.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingMpiTraits.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingSpawnSlaves.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingTimeStepper.hpp>

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
// NOTE: All slave-master communicators have set a custom error handler, which eventually
//   will call MPI_Abort() so there is no need to check the return value of any MPI_Recv()
//   or MPI_Send() calls below.

template <class Scalar>
ReservoirCouplingMaster<Scalar>::
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

template <class Scalar>
const std::vector<std::string>&
ReservoirCouplingMaster<Scalar>::
getMasterGroupNamesForSlave(std::size_t slave_idx) const
{
    if (slave_idx >= this->slave_idx_to_master_groups_.size()) {
        OPM_THROW(
            std::runtime_error,
            fmt::format(
                "Slave index {} out of bounds. Valid range is [0, {})",
                slave_idx, this->slave_idx_to_master_groups_.size()
            )
        );
    }

    // Direct O(1) vector access for performance-critical usage
    return this->slave_idx_to_master_groups_[slave_idx];
}

template <class Scalar>
std::size_t
ReservoirCouplingMaster<Scalar>::
getMasterGroupCanonicalIdx(
    const std::string &slave_name, const std::string &master_group_name) const
{
    // NOTE: The master group name order is the canonical order of the master group names
    //       for a given slave name. This is the order in which the slave will communicate
    //       slave group data for its slave groups.
    return this->master_group_name_order_.at(slave_name).at(master_group_name);
}

template <class Scalar>
void
ReservoirCouplingMaster<Scalar>::
initTimeStepping()
{
    assert(!this->time_stepper_);
    this->time_stepper_ = std::make_unique<ReservoirCouplingTimeStepper<Scalar>>(*this);
    assert(!this->report_step_data_);
    this->report_step_data_ = std::make_unique<ReservoirCouplingMasterReportStep<Scalar>>(*this);
}

template <class Scalar>
void
ReservoirCouplingMaster<Scalar>::
initStartOfReportStep(int report_step_idx)
{
    assert(this->report_step_data_);
    this->report_step_data_->setReportStepIdx(report_step_idx);
    this->logger_.info("Initializing start of report step");
}

template <class Scalar>
bool
ReservoirCouplingMaster<Scalar>::
isMasterGroup(const std::string &group_name) const
{
    return this->master_group_slave_names_.find(group_name) !=
           this->master_group_slave_names_.end();
}

template <class Scalar>
void
ReservoirCouplingMaster<Scalar>::
maybeActivate(int report_step) {
    if (!this->activated()) {
        double start_date = this->schedule_.getStartTime();
        auto current_time = start_date + this->schedule_.seconds(report_step);
        if (Seconds::compare_gt_or_eq(current_time, this->activation_date_)) {
            this->activated_ = true;
        }
    }
}

template <class Scalar>
double
ReservoirCouplingMaster<Scalar>::
maybeChopSubStep(double suggested_timestep, double current_time) const
{
    assert(this->time_stepper_);
    return this->time_stepper_->maybeChopSubStep(suggested_timestep, current_time);
}

template <class Scalar>
void
ReservoirCouplingMaster<Scalar>::
maybeSpawnSlaveProcesses(int report_step)
{
    if (this->numSlavesStarted() > 0) {  // We have already spawned the slave processes
        return;
    }
    const auto& rescoup = this->schedule_[report_step].rescoup();
    auto slave_count = rescoup.slaveCount();
    auto master_group_count = rescoup.masterGroupCount();
    if (slave_count > 0 && master_group_count > 0) {
        ReservoirCouplingSpawnSlaves<Scalar> spawn_slaves{*this, rescoup};
        spawn_slaves.spawn();
    }
}

template <class Scalar>
void
ReservoirCouplingMaster<Scalar>::
maybeReceiveActivationHandshakeFromSlaves(double current_time)
{
    // Initialize on first call
    if (this->slave_activation_status_.empty()) {
        this->slave_activation_status_.resize(this->numSlavesStarted(), false);
    }

    if (this->comm_.rank() == 0) {
        auto current_date = this->schedule_.getStartTime() + current_time;
        for (unsigned int i = 0; i < this->numSlavesStarted(); i++) {
            // Skip if already activated
            if (this->slaveIsActivated(i)) {
                continue;
            }
            // Check if slave should activate during this timestep
            double slave_activation_date = this->slave_activation_dates_[i];
            // NOTE: The master will adjust its time stepping to ensure that its step will always
            //    conincide with the slave report steps (and hence the slave activation date)
            assert(Seconds::compare_gt_or_eq(slave_activation_date, current_date));
            if (Seconds::compare_eq(slave_activation_date, current_date)) {
                // Use non-blocking probe first to check for handshake
                int flag;
                MPI_Status status;
                MPI_Iprobe(0, static_cast<int>(MessageTag::SlaveActivationHandshake),
                          this->master_slave_comm_[i], &flag, &status);

                if (!flag) {
                    // Handshake not received yet, for debugging reasons it can be useful to know which slave
                    //   is waiting for the handshake
                    OpmLog::info(fmt::format("Waiting for activation handshake from slave {}",
                        this->slave_names_[i]));
                }
                std::uint8_t handshake;
                auto MPI_UINT8_T_TYPE = Dune::MPITraits<std::uint8_t>::getType();
                MPI_Recv(&handshake, 1, MPI_UINT8_T_TYPE, 0,
                        static_cast<int>(MessageTag::SlaveActivationHandshake),
                        this->master_slave_comm_[i], MPI_STATUS_IGNORE);
                this->slave_activation_status_[i] = handshake;
                OpmLog::info(fmt::format("Received activation handshake from slave {}",
                                        this->slave_names_[i]));

            }
        }
    }
    // Broadcast status to all ranks
    comm_.broadcast(
        this->slave_activation_status_.data(), /*count=*/this->numSlavesStarted(), /*emitter_rank=*/0
    );
}

template <class Scalar>
std::size_t
ReservoirCouplingMaster<Scalar>::
numSlaveGroups(unsigned int index)
{
    return this->master_group_name_order_[this->slave_names_[index]].size();
}

template <class Scalar>
std::size_t
ReservoirCouplingMaster<Scalar>::
numSlavesStarted() const
{
    return this->slave_names_.size();
}

template <class Scalar>
void
ReservoirCouplingMaster<Scalar>::
receiveNextReportDateFromSlaves()
{
    assert(this->time_stepper_);
    this->time_stepper_->receiveNextReportDateFromSlaves();
}

template <class Scalar>
void
ReservoirCouplingMaster<Scalar>::
rebuildSlaveIdxToMasterGroupsVector()
{
    // Rebuild the index-based vector from the name-based map
    // This ensures consistent ordering with slave_names_ vector
    this->slave_idx_to_master_groups_.clear();
    this->slave_idx_to_master_groups_.resize(this->slave_names_.size());

    for (std::size_t slave_idx = 0; slave_idx < this->slave_names_.size(); ++slave_idx) {
        const auto& slave_name = this->slave_names_[slave_idx];
        auto it = this->slave_name_to_master_groups_map_.find(slave_name);
        if (it != this->slave_name_to_master_groups_map_.end()) {
            this->slave_idx_to_master_groups_[slave_idx] = it->second;
        }
        // If slave_name not found in map, leave the vector empty for that index
    }
}

template <class Scalar>
void
ReservoirCouplingMaster<Scalar>::
receiveInjectionDataFromSlaves()
{
    assert(this->report_step_data_);
    this->report_step_data_->receiveInjectionDataFromSlaves();
}

template <class Scalar>
void
ReservoirCouplingMaster<Scalar>::
receiveProductionDataFromSlaves()
{
    assert(this->report_step_data_);
    this->report_step_data_->receiveProductionDataFromSlaves();
}

template <class Scalar>
const ReservoirCoupling::Potentials<Scalar>&
ReservoirCouplingMaster<Scalar>::
getSlaveGroupPotentials(const std::string &master_group_name)
{
    assert(this->report_step_data_);
    return this->report_step_data_->getSlaveGroupPotentials(master_group_name);
}

template <class Scalar>
void
ReservoirCouplingMaster<Scalar>::
resizeNextReportDates(int size)
{
    assert(this->time_stepper_);
    this->time_stepper_->resizeNextReportDates(size);
}

template <class Scalar>
void
ReservoirCouplingMaster<Scalar>::
setSlaveNextReportTimeOffset(int index, double offset)
{
    assert(this->time_stepper_);
    this->time_stepper_->setSlaveNextReportTimeOffset(index, offset);
}

template <class Scalar>
void
ReservoirCouplingMaster<Scalar>::
updateMasterGroupNameOrderMap(
    const std::string& slave_name, const std::map<std::string, std::size_t>& master_group_name_map
)
{
    this->master_group_name_order_[slave_name] = master_group_name_map;
}

// ------------------
// Private methods
// ------------------

template <class Scalar>
double
ReservoirCouplingMaster<Scalar>::
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


// Explicit template instantiations
template class ReservoirCouplingMaster<double>;
#if FLOW_INSTANTIATE_FLOAT
template class ReservoirCouplingMaster<float>;
#endif

} // namespace Opm
