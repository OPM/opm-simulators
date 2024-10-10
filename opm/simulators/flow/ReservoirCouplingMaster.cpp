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

#include <opm/simulators/flow/ReservoirCouplingMaster.hpp>

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

ReservoirCouplingMaster::ReservoirCouplingMaster(
    const Parallel::Communication &comm,
    const Schedule &schedule
) :
    comm_{comm},
    schedule_{schedule}
{ }

// ------------------
// Public methods
// ------------------

void ReservoirCouplingMaster::receiveSimulationStartDateFromSlaves() {
    this->slave_start_dates_.resize(this->num_slaves_);
    if (this->comm_.rank() == 0) {
        // Ensure that std::time_t is of type long since we are sending it over MPI with MPI_LONG
        static_assert(std::is_same<std::time_t, long>::value, "std::time_t is not of type long");
        for (unsigned int i = 0; i < this->master_slave_comm_.size(); i++) {
            std::time_t start_date;
            int result = MPI_Recv(
                &start_date,
                /*count=*/1,
                /*datatype=*/MPI_LONG,
                /*source_rank=*/0,
                /*tag=*/static_cast<int>(MessageTag::SimulationStartDate),
                *this->master_slave_comm_[i].get(),
                MPI_STATUS_IGNORE
            );
            if (result != MPI_SUCCESS) {
                OPM_THROW(std::runtime_error, "Failed to receive simulation start date from slave process");
            }
            this->slave_start_dates_[i] = start_date;
            OpmLog::info(
                fmt::format(
                    "Received simulation start date from slave process with name: {}. "
                    "Start date: {}", this->slave_names_[i], start_date
                )
            );
        }
    }
    this->comm_.broadcast(this->slave_start_dates_.data(), this->num_slaves_, /*emitter_rank=*/0);
}

// NOTE: This functions is executed for all ranks, but only rank 0 will spawn
//   the slave processes
void ReservoirCouplingMaster::spawnSlaveProcesses(int argc, char **argv) {
    const auto& rescoup = this->schedule_[0].rescoup();
    char *flow_program_name = argv[0];
    for (const auto& [slave_name, slave] : rescoup.slaves()) {
        auto master_slave_comm = MPI_Comm_Ptr(new MPI_Comm(MPI_COMM_NULL));
        const auto& data_file_name = slave.dataFilename();
        const auto& directory_path = slave.directoryPath();
        // Concatenate the directory path and the data file name to get the full path
        std::filesystem::path dir_path(directory_path);
        std::filesystem::path data_file(data_file_name);
        std::filesystem::path full_path = dir_path / data_file;
        std::string log_filename; // the getSlaveArgv() function will set this
        std::vector<char *> slave_argv = getSlaveArgv(
            argc, argv, full_path, slave_name, log_filename
        );
        auto num_procs = slave.numprocs();
        std::vector<int> errcodes(num_procs);
        // TODO: We need to decide how to handle the output from the slave processes..
        //    As far as I can tell, open MPI does not support redirecting the output
        //    to a file, so we might need to implement a custom solution for this
        int spawn_result = MPI_Comm_spawn(
            flow_program_name,
            slave_argv.data(),
            /*maxprocs=*/num_procs,
            /*info=*/MPI_INFO_NULL,
            /*root=*/0,  // Rank 0 spawns the slave processes
            /*comm=*/this->comm_,
            /*intercomm=*/master_slave_comm.get(),
            /*array_of_errcodes=*/errcodes.data()
        );
        if (spawn_result != MPI_SUCCESS || (*master_slave_comm == MPI_COMM_NULL)) {
            for (unsigned int i = 0; i < num_procs; i++) {
                if (errcodes[i] != MPI_SUCCESS) {
                    char error_string[MPI_MAX_ERROR_STRING];
                    int length_of_error_string;
                    MPI_Error_string(errcodes[i], error_string, &length_of_error_string);
                    OpmLog::info(fmt::format("Error spawning process {}: {}", i, error_string));
                }
            }
            OPM_THROW(std::runtime_error, "Failed to spawn slave process");
        }
        OpmLog::info(
            fmt::format(
                "Spawned reservoir coupling slave simulation for slave with name: "
                "{}. Standard output logfile name: {}.log", slave_name, slave_name
            )
        );
        this->master_slave_comm_.push_back(std::move(master_slave_comm));
        this->slave_names_.push_back(slave_name);
        this->num_slaves_++;
    }
}

// ------------------
// Private methods
// ------------------

std::vector<char *> ReservoirCouplingMaster::getSlaveArgv(
    int argc,
    char **argv,
    const std::filesystem::path &data_file,
    const std::string &slave_name,
    std::string &log_filename
) {
    // Calculate the size of the slave_argv vector like this:
    // - We will not use the first argument in argv, as it is the program name
    // - Replace the data file name in argv with the data_file path
    // - Insert as first argument --slave-log-file=<slave_name>.log
    // - Also add the argument "--slave=true" to the argv
    // - Add a nullptr at the end of the argv
    // So the size of the slave_argv vector will be argc + 2
    //
    // Assume: All parameters will be on the form --parameter=value (as a string without spaces)
    //
    // Important: The returned vector will have pointers to argv pointers,
    //   data_file string buffer, and slave_name string buffer. So the caller
    //   must ensure that these buffers are not deallocated before the slave_argv has
    //   been used.
    std::vector<char *> slave_argv(argc + 2);
    log_filename = "--slave-log-file=" + slave_name;  // .log extension will be added by the slave process
    slave_argv[0] = const_cast<char*>(log_filename.c_str());
    for (int i = 1; i < argc; i++) {
        // Check if the argument starts with "--", if not, we will assume it is a positional argument
        //   and we will replace it with the data file path
        if (std::string(argv[i]).substr(0, 2) == "--") {
            slave_argv[i] = argv[i];
        } else {
            slave_argv[i] = const_cast<char*>(data_file.c_str());
        }
    }
    slave_argv[argc] = const_cast<char *>("--slave=true");
    slave_argv[argc+1] = nullptr;
    return slave_argv;
}

} // namespace Opm

