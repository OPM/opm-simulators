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

#ifndef OPM_RESERVOIR_COUPLING_SLAVE_HPP
#define OPM_RESERVOIR_COUPLING_SLAVE_HPP

#include <opm/simulators/flow/ReservoirCouplingMaster.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <mpi.h>

#include <vector>

namespace Opm {

class ReservoirCouplingSlave {
public:
    using MPI_Comm_Ptr = ReservoirCouplingMaster::MPI_Comm_Ptr;

    ReservoirCouplingSlave(const Parallel::Communication &comm, const Schedule &schedule);
    void sendSimulationStartDateToMasterProcess();

private:
    const Parallel::Communication &comm_;
    const Schedule& schedule_;
    // MPI parent communicator for a slave process
    MPI_Comm_Ptr slave_master_comm_{nullptr};

};

} // namespace Opm
#endif // OPM_RESERVOIR_COUPLING_SLAVE_HPP
