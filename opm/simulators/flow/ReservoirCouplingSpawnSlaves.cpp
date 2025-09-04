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

#include <opm/simulators/flow/ReservoirCouplingSpawnSlaves.hpp>

#include <opm/input/eclipse/Schedule/ResCoup/ReservoirCouplingInfo.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/MasterGroup.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/Slaves.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <dune/common/parallel/mpitraits.hh>

#include <algorithm>
#include <filesystem>
#include <vector>

#include <fmt/format.h>

namespace Opm {
// NOTE: All slave-master communicators have set a custom error handler, which eventually
//   will call MPI_Abort() so there is no need to check the return value of any MPI_Recv()
//   or MPI_Send() calls below.

// ------------------
// Public methods
// ------------------

ReservoirCouplingSpawnSlaves::
ReservoirCouplingSpawnSlaves(
    ReservoirCouplingMaster &master,
    const ReservoirCoupling::CouplingInfo &rescoup
) :
    master_{master},
    rescoup_{rescoup},
    comm_{master.getComm()}
{
}

void
ReservoirCouplingSpawnSlaves::
spawn()
{
    // spawn the slave processes and get the simulation start date from the slaves,
    // and finally send the master group names to the slaves
    this->spawnSlaveProcesses_();
    this->receiveActivationDateFromSlaves_();
    this->receiveSimulationStartDateFromSlaves_();
    this->sendSlaveNamesToSlaves_();
    this->createMasterGroupNameOrder_();
    this->createMasterGroupToSlaveNameMap_();
    this->sendMasterGroupNamesToSlaves_();
    this->prepareTimeStepping_();
    this->logger_.info("Reservoir coupling slave processes was spawned successfully");
}


// ------------------
// Private methods
// ------------------

void
ReservoirCouplingSpawnSlaves::
createMasterGroupNameOrder_()
{
    // When the slaves send master/slave group potentials to us, we need to know
    // which master group corresponds to which potentials. We will use the convention
    // that the slaves will send the potentials in the order of lexicographically sorted
    // slave group names
    auto num_slaves = this->master_.numSlavesStarted();
    const auto& master_groups = this->rescoup_.masterGroups();
    for (unsigned int i = 0; i < num_slaves; i++) {
        auto slave_name = this->master_.getSlaveName(i);
        std::vector<std::pair<std::string, std::string>> slave_group_names;
        for (const auto& [master_group_name, master_group] : master_groups) {
            if (master_group.slaveName() == slave_name) {
                slave_group_names.push_back({master_group_name, master_group.slaveGroupName()});
            }
        }
        // Sort the vector based on the slave group names in lexicographical order
        std::sort(
            slave_group_names.begin(),
            slave_group_names.end(),
            [](const auto &lhs, const auto &rhs) {
                return lhs.second < rhs.second;
            }
        );
        // Create a map from the master group name to the index in the vector
        // This is used to determine the order in which the slaves will send the potentials
        std::map<std::string, std::size_t> master_group_map;
        for (std::size_t j = 0; j < slave_group_names.size(); ++j) {
            master_group_map[slave_group_names[j].first] = j;
        }
        this->master_.updateMasterGroupNameOrderMap(slave_name, master_group_map);
    }
}
void
ReservoirCouplingSpawnSlaves::
createMasterGroupToSlaveNameMap_()
{
    // When the slaves send master/slave group potentials to us, we need to know
    // which master group corresponds to which potentials. The potentials are associated
    // with a slave name. So we make lookup table: master_group_name -> slave_name.
    const auto& master_groups = this->rescoup_.masterGroups();
    auto& master_group_slave_names = this->master_.getMasterGroupToSlaveNameMap();
    for (const auto& [master_group_name, master_group] : master_groups) {
        master_group_slave_names[master_group_name] = master_group.slaveName();
    }
}

std::vector<char *> ReservoirCouplingSpawnSlaves::
getSlaveArgv_(
    const std::filesystem::path &data_file,
    const std::string &slave_name,
    std::string &log_filename) const
{
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
    auto argc = this->master_.getArgc();
    auto argv = this->master_.getArgv();
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

std::pair<std::vector<char>, std::size_t>
ReservoirCouplingSpawnSlaves::
getMasterGroupNamesForSlave_(const std::string &slave_name) const
{
    // For the given slave name, get all pairs of master group names and slave group names
    // Serialize the data such that it can be sent over MPI in one chunk
    auto master_groups = this->rescoup_.masterGroups();
    std::vector<std::string> data;
    std::vector<std::string> master_group_names;
    for (const auto& [group_name, master_group] : master_groups) {
        if (master_group.slaveName() == slave_name) {
            data.push_back(group_name);
            data.push_back(master_group.slaveGroupName());
        }
    }
    assert(data.size() % 2 == 0);
    assert(data.size() > 0);
    return ReservoirCoupling::serializeStrings(data);
}

void
ReservoirCouplingSpawnSlaves::
prepareTimeStepping_()
{
    // Prepare the time stepping for the master process
    // This is done after the slave processes have been spawned
    // and the master group names have been sent to the slaves
    auto num_slaves = this->master_.numSlavesStarted();
    this->master_.resizeNextReportDates(num_slaves);
}

void
ReservoirCouplingSpawnSlaves::
receiveActivationDateFromSlaves_()
{
    // Currently, we only use the activation date to check that no slave process
    // starts before the master process.
    auto num_slaves = this->master_.numSlavesStarted();
    if (this->comm_.rank() == 0) {
        for (unsigned int i = 0; i < num_slaves; i++) {
            double slave_activation_date;
            // NOTE: See comment about error handling at the top of this file.
            MPI_Recv(
                &slave_activation_date,
                /*count=*/1,
                /*datatype=*/MPI_DOUBLE,
                /*source_rank=*/0,
                /*tag=*/static_cast<int>(MessageTag::SlaveActivationDate),
                this->master_.getSlaveComm(i),
                MPI_STATUS_IGNORE
            );
            if (slave_activation_date < this->master_.getActivationDate()) {
                OPM_THROW(std::runtime_error, "Slave process activation date is earlier than "
                                              "the master process' activation date");
            }
            this->master_.addSlaveActivationDate(slave_activation_date);
            this->logger_.info(
                fmt::format(
                    "Received activation date from slave process with name: {}. "
                    "Activation date: {}", this->master_.getSlaveName(i), slave_activation_date
                )
            );
        }
    }
    if (this->comm_.rank() != 0) {
        // Ensure that all ranks have the same number of slave activation dates
        this->master_.resizeSlaveActivationDates(num_slaves);
    }
    const double* data = this->master_.getSlaveActivationDates();
    this->comm_.broadcast(const_cast<double *>(data), /*count=*/num_slaves, /*emitter_rank=*/0);
    this->logger_.info("Broadcasted slave activation dates to all ranks");
}

void
ReservoirCouplingSpawnSlaves::
receiveSimulationStartDateFromSlaves_()
{
    auto num_slaves = this->master_.numSlavesStarted();
    if (this->comm_.rank() == 0) {
        for (unsigned int i = 0; i < num_slaves; i++) {
            double slave_start_date;
            // NOTE: See comment about error handling at the top of this file.
            MPI_Recv(
                &slave_start_date,
                /*count=*/1,
                /*datatype=*/MPI_DOUBLE,
                /*source_rank=*/0,
                /*tag=*/static_cast<int>(MessageTag::SlaveSimulationStartDate),
                this->master_.getSlaveComm(i),
                MPI_STATUS_IGNORE
            );
            this->master_.addSlaveStartDate(slave_start_date);
            this->logger_.info(
                fmt::format(
                    "Received start date from slave process with name: {}. "
                    "Start date: {}", this->master_.getSlaveName(i), slave_start_date
                )
            );
        }
    }
    if (this->comm_.rank() != 0) {
        // Ensure that all ranks have the same number of slave start dates
        this->master_.resizeSlaveStartDates(num_slaves);
    }
    const double* data = this->master_.getSlaveStartDates();
    this->comm_.broadcast(const_cast<double *>(data), /*count=*/num_slaves, /*emitter_rank=*/0);
    this->logger_.info("Broadcasted slave start dates to all ranks");
}

void
ReservoirCouplingSpawnSlaves::
sendMasterGroupNamesToSlaves_()
{
    if (this->comm_.rank() == 0) {
        auto num_slaves = this->master_.numSlavesStarted();
        for (unsigned int i = 0; i < num_slaves; i++) {
            auto slave_name = this->master_.getSlaveName(i);
            auto [group_names, size] = this->getMasterGroupNamesForSlave_(slave_name);
            static_assert(std::is_same_v<decltype(size), std::size_t>, "size must be of type std::size_t");
            auto MPI_SIZE_T_TYPE = Dune::MPITraits<std::size_t>::getType();
            // NOTE: See comment about error handling at the top of this file.
            MPI_Send(
                &size,
                /*count=*/1,
                /*datatype=*/MPI_SIZE_T_TYPE,
                /*dest_rank=*/0,
                /*tag=*/static_cast<int>(MessageTag::MasterGroupNamesSize),
                this->master_.getSlaveComm(i)
            );
            MPI_Send(
                group_names.data(),
                /*count=*/size,
                /*datatype=*/MPI_CHAR,
                /*dest_rank=*/0,
                /*tag=*/static_cast<int>(MessageTag::MasterGroupNames),
                this->master_.getSlaveComm(i)
            );
            this->logger_.info(fmt::format(
                "Sent master group names to slave process rank 0 with name: {}", slave_name)
            );
        }
   }
}

void
ReservoirCouplingSpawnSlaves::
sendSlaveNamesToSlaves_()
{
    if (this->comm_.rank() == 0) {
        auto num_slaves = this->master_.numSlavesStarted();
        for (unsigned int i = 0; i < num_slaves; i++) {
            const auto& slave_name = this->master_.getSlaveName(i);
            auto slave_name_size = slave_name.size();
            auto MPI_SIZE_T_TYPE = Dune::MPITraits<std::size_t>::getType();
            // NOTE: See comment about error handling at the top of this file.
            MPI_Send(
                &slave_name_size,
                /*count=*/1,
                /*datatype=*/MPI_SIZE_T_TYPE,
                /*dest_rank=*/0,
                /*tag=*/static_cast<int>(MessageTag::SlaveNameSize),
                this->master_.getSlaveComm(i)
            );
            MPI_Send(
                slave_name.c_str(),
                /*count=*/slave_name_size,
                /*datatype=*/MPI_CHAR,
                /*dest_rank=*/0,
                /*tag=*/static_cast<int>(MessageTag::SlaveName),
                this->master_.getSlaveComm(i)
            );
            this->logger_.info(fmt::format(
                "Sent slave name to slave process rank 0 with name: {}", slave_name)
            );
        }
    }
}


// NOTE: This functions is executed for all ranks, but only rank 0 will spawn
//   the slave processes
void
ReservoirCouplingSpawnSlaves::
spawnSlaveProcesses_()
{
    char *flow_program_name = this->master_.getArgv(0);
    for (const auto& [slave_name, slave] : this->rescoup_.slaves()) {
        MPI_Comm master_slave_comm = MPI_COMM_NULL;
        const auto& data_file_name = slave.dataFilename();
        const auto& directory_path = slave.directoryPath();
        // Concatenate the directory path and the data file name to get the full path
        std::filesystem::path dir_path{directory_path};
        std::filesystem::path data_file{data_file_name};
        std::filesystem::path full_path = dir_path / data_file;
        std::string log_filename; // the getSlaveArgv() function will set this
        std::vector<char *> slave_argv = this->getSlaveArgv_(
            full_path, slave_name, log_filename
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
            /*intercomm=*/&master_slave_comm,
            /*array_of_errcodes=*/errcodes.data()
        );
        if (spawn_result != MPI_SUCCESS || (master_slave_comm == MPI_COMM_NULL)) {
            for (unsigned int i = 0; i < num_procs; i++) {
                if (errcodes[i] != MPI_SUCCESS) {
                    char error_string[MPI_MAX_ERROR_STRING];
                    int length_of_error_string;
                    MPI_Error_string(errcodes[i], error_string, &length_of_error_string);
                    this->logger_.info(fmt::format("Error spawning process {}: {}", i, error_string));
                }
            }
            OPM_THROW(std::runtime_error, "Failed to spawn slave process");
        }
        // NOTE: By installing a custom error handler for all slave-master communicators, which
        //   eventually will call MPI_Abort(), there is no need to check the return value of any
        //   MPI_Recv() or MPI_Send() calls as errors will be caught by the error handler.
        ReservoirCoupling::setErrhandler(master_slave_comm, /*is_master=*/true);
        this->logger_.info(
            fmt::format(
                "Spawned reservoir coupling slave simulation for slave with name: "
                "{}. Standard output logfile name: {}.log", slave_name, slave_name
            )
        );
        this->master_.addSlaveCommunicator(master_slave_comm);
        this->master_.addSlaveName(slave_name);
    }
}

}  // namespace Opm
