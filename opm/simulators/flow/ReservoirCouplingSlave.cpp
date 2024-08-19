/*
  Copyright 2024 Equinor AS

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
#include <opm/simulators/flow/ReservoirCoupling.hpp>
#include <opm/simulators/flow/ReservoirCouplingSlave.hpp>

#include <opm/input/eclipse/Schedule/ResCoup/ReservoirCouplingInfo.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/MasterGroup.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/Slaves.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <vector>
#include <fmt/format.h>

namespace Opm {

ReservoirCouplingSlave::ReservoirCouplingSlave(
    const Parallel::Communication &comm,
    const Schedule &schedule,
    const SimulatorTimer &timer
) :
    comm_{comm},
    schedule_{schedule},
    timer_{timer}
{
    this->slave_master_comm_ = MPI_Comm_Ptr(new MPI_Comm(MPI_COMM_NULL));
    MPI_Comm_get_parent(this->slave_master_comm_.get());
    if (*(this->slave_master_comm_) == MPI_COMM_NULL) {
        OPM_THROW(std::runtime_error, "Slave process is not spawned by a master process");
    }
}

double ReservoirCouplingSlave::receiveNextTimeStepFromMaster() {
    double timestep;
    if (this->comm_.rank() == 0) {
        int result = MPI_Recv(
            &timestep,
            /*count=*/1,
            /*datatype=*/MPI_DOUBLE,
            /*source_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::SlaveNextTimeStep),
            *this->slave_master_comm_,
            MPI_STATUS_IGNORE
        );
        if (result != MPI_SUCCESS) {
            OPM_THROW(std::runtime_error, "Failed to receive next time step from master");
        }
        OpmLog::info(
            fmt::format("Slave rank 0 received next timestep {} from master.", timestep)
        );
    }
    this->comm_.broadcast(&timestep, 1, /*emitter_rank=*/0);
    OpmLog::info("Broadcasted slave next time step to all ranks");
    return timestep;
}


void ReservoirCouplingSlave::sendNextReportDateToMasterProcess() {
    if (this->comm_.rank() == 0) {
        double elapsed_time = this->timer_.simulationTimeElapsed();
        double current_step_length = this->timer_.currentStepLength();
        // NOTE: This is an offset in seconds from the start date, so it will be 0 if the next report
        //      would be the start date. In general, it should be a positive number.
        double next_report_time_offset = elapsed_time + current_step_length;
        MPI_Send(
            &next_report_time_offset,
            /*count=*/1,
            /*datatype=*/MPI_DOUBLE,
            /*dest_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::SlaveNextReportDate),
            *this->slave_master_comm_
        );
        OpmLog::info("Sent next report date to master process from rank 0");
   }
}

void ReservoirCouplingSlave::sendSimulationStartDateToMasterProcess() {
    if (this->comm_.rank() == 0) {
        // Ensure that std::time_t is of type long since we are sending it over MPI with MPI_LONG
        static_assert(std::is_same<std::time_t, long>::value, "std::time_t is not of type long");
        std::time_t start_date = this->schedule_.getStartTime();
        MPI_Send(
            &start_date,
            /*count=*/1,
            /*datatype=*/MPI_LONG,
            /*dest_rank=*/0,
            /*tag=*/static_cast<int>(MessageTag::SlaveSimulationStartDate),
            *this->slave_master_comm_
        );
        OpmLog::info("Sent simulation start date to master process from rank 0");
   }
}

} // namespace Opm
