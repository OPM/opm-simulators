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

#ifndef OPM_RESERVOIR_COUPLING_MASTER_HPP
#define OPM_RESERVOIR_COUPLING_MASTER_HPP

#include <opm/simulators/utils/ParallelCommunication.hpp>
#include <opm/simulators/flow/ReservoirCoupling.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <mpi.h>

#include <filesystem>
#include <vector>

namespace Opm {

class ReservoirCouplingMaster {
public:
    using MPI_Comm_Ptr = ReservoirCoupling::MPI_Comm_Ptr;
    using MessageTag = ReservoirCoupling::MessageTag;

    ReservoirCouplingMaster(const Parallel::Communication &comm, const Schedule &schedule);

    void spawnSlaveProcesses(int argc, char **argv);
    void receiveSimulationStartDateFromSlaves();

private:
    std::vector<char *> getSlaveArgv(
        int argc,
        char **argv,
        const std::filesystem::path &data_file,
        const std::string &slave_name,
        std::string &log_filename
    );

    const Parallel::Communication &comm_;
    const Schedule& schedule_;
    std::size_t num_slaves_ = 0;  // Initially zero, will be updated in spawnSlaveProcesses()
    std::vector<MPI_Comm_Ptr> master_slave_comm_; // MPI communicators for the slave processes
    std::vector<std::string> slave_names_;
    std::vector<std::time_t> slave_start_dates_;
};

} // namespace Opm
#endif // OPM_RESERVOIR_COUPLING_MASTER_HPP
