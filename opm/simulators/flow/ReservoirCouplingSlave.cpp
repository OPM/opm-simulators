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
#include <opm/simulators/flow/ReservoirCouplingSlave.hpp>

#include <opm/input/eclipse/Schedule/ResCoup/ReservoirCouplingInfo.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/MasterGroup.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/Slaves.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <vector>

namespace Opm {

ReservoirCouplingSlave::ReservoirCouplingSlave(
    const Parallel::Communication &comm,
    const Schedule &schedule
) :
    comm_{comm},
    schedule_{schedule}
{ }
void ReservoirCouplingSlave::sendSimulationStartDateToMasterProcess() {
    // TODO: Implement this function next

    //this->slave_master_comm_ = MPI_Comm_Ptr(new MPI_Comm(MPI_COMM_NULL));
    //MPI_Comm_get_parent(this->slave_master_comm_.get());
    //if (*(this->slave_master_comm_) == MPI_COMM_NULL) {
    //    OPM_THROW(std::runtime_error, "Slave process is not spawned by a master process");
    //}
    OpmLog::info("Sent simulation start date to master process");
}

} // namespace Opm
