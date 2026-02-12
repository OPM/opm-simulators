/*
  Copyright 2025 Equinor ASA

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
#include <opm/simulators/flow/rescoup/ReservoirCouplingMpiTraits.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingMaster.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingTimeStepper.hpp>

#include <opm/input/eclipse/Schedule/ResCoup/ReservoirCouplingInfo.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/MasterGroup.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/Slaves.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <dune/common/parallel/mpitraits.hh>

#include <vector>
#include <fmt/format.h>

namespace Opm {

template <class Scalar>
ReservoirCouplingTimeStepper<Scalar>::
ReservoirCouplingTimeStepper(
    ReservoirCouplingMaster<Scalar> &master
) :
    master_{master}
{
    this->slave_next_report_time_offsets_.resize(this->numSlaves());
}

// ------------------
// Public methods
// ------------------

template <class Scalar>
double
ReservoirCouplingTimeStepper<Scalar>::
maybeChopSubStep(double suggested_timestep_original, double elapsed_time) const
{
    // Check if the suggested timestep needs to be adjusted based on the slave processes'
    // next report step, or if the slave process has not started yet: the start of a slave process.
    // NOTE: getStartTime() returns a std::time_t value, which is typically a long integer. It should
    //     be possible to represent reasonable epoch values within a double. See comment for
    //     getMasterActivationDate_() for more information.
    double start_date = this->schedule().getStartTime();
    double step_start_date{start_date + elapsed_time};
    double step_end_date{step_start_date + suggested_timestep_original};
    double suggested_timestep{suggested_timestep_original};
    auto num_slaves = this->numSlaves();
    // Determine the minimum step_end_date and the corresponding suggested_timestep such that no
    // slave process will report or start during the timestep [step_start_date, step_end_date]
    // where suggested_timestep = step_end_date - step_start_date
    for (std::size_t i = 0; i < num_slaves; i++) {
        double slave_start_date = this->slaveStartDate(i);
        double slave_activation_date = this->slaveActivationDate(i);
        double slave_next_report_date{this->slave_next_report_time_offsets_[i] + slave_start_date};
        if (Seconds::compare_gt_or_eq(slave_activation_date, step_end_date)) {
            // The slave process has not activated yet, and will not activate during this timestep
            continue;
        }
        double slave_elapsed_time;
        if (Seconds::compare_lt_or_eq(slave_activation_date,step_start_date)) {
            // The slave process has already activated, and will continue during this timestep
            if (Seconds::compare_gt(slave_next_report_date, step_end_date)) {
                // The slave process will not report during this timestep
                continue;
            }
            // The slave process will report during this timestep
            slave_elapsed_time = slave_next_report_date - step_start_date;
        }
        else {
            // The slave process will activate during the timestep, but not at the beginning
            // Ensure that we synchronize the slave activation report date with the
            //   start of a master timestep by chopping the current master timestep.
            slave_elapsed_time = slave_activation_date - step_start_date;
        }
        suggested_timestep = slave_elapsed_time;
        step_end_date = step_start_date + suggested_timestep;
    }
    return suggested_timestep;
}


template <class Scalar>
void
ReservoirCouplingTimeStepper<Scalar>::
receiveNextReportDateFromSlaves()
{
    auto num_slaves = this->numSlaves();
    if (this->comm().rank() == 0) {
        this->logger().debug("Receiving next report dates from slave processes");
        for (unsigned int i = 0; i < num_slaves; i++) {
            if (!this->slaveIsActivated(i)) {
                // Set to zero to indicate that the slave has not activated yet
                this->slave_next_report_time_offsets_[i] = 0.0;
                this->logger().debug(fmt::format(
                    "Slave {} has not activated yet, setting next report date to 0.0",
                    this->slaveName(i)));
                continue;
            }
            double slave_next_report_time_offset; // Elapsed time from the beginning of the simulation
            // NOTE: See comment about error handling at the top of this file.
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
            this->logger().debug(fmt::format(
                "Received next report date from {}: {} (offset from slave start)",
                this->slaveName(i), ReservoirCoupling::formatDays(slave_next_report_time_offset)
            ));
        }
    }
    this->comm().broadcast(
        this->slave_next_report_time_offsets_.data(), /*count=*/num_slaves, /*emitter_rank=*/0
    );
    this->logger().debug("Broadcasted slave next report dates to all ranks");
}

template <class Scalar>
void
ReservoirCouplingTimeStepper<Scalar>::
sendNextTimeStepToSlaves(double timestep)
{
    if (this->comm().rank() == 0) {
        for (unsigned int slave_idx = 0; slave_idx < this->numSlaves(); slave_idx++) {
            if (!this->slaveIsActivated(slave_idx)) {
                this->logger().debug(fmt::format(
                    "Slave {} has not activated yet, skipping sending next time step",
                    this->slaveName(slave_idx)
                ));
                continue;
            }
            // Then send the next time step
            // NOTE: See comment about error handling at the top of this file.
            MPI_Send(
                &timestep,
                /*count=*/1,
                /*datatype=*/MPI_DOUBLE,
                /*dest_rank=*/0,
                /*tag=*/static_cast<int>(MessageTag::SlaveNextTimeStep),
                this->getSlaveComm(slave_idx)
            );
            this->logger().debug(fmt::format(
                "Sent next time step {} to {}",
                ReservoirCoupling::formatDays(timestep), this->slaveName(slave_idx))
            );
        }
   }
}

// ------------------
// Private methods
// ------------------


// Explicit instantiations
template class ReservoirCouplingTimeStepper<double>;
#if FLOW_INSTANTIATE_FLOAT
template class ReservoirCouplingTimeStepper<float>;
#endif

} // namespace Opm
