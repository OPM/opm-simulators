/*
  Copyright 2013, 2014, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 IRIS AS
  Copyright 2014 STATOIL ASA.

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
#include <opm/simulators/flow/Main.hpp>

#include <opm/input/eclipse/Python/Python.hpp>
#include <opm/input/eclipse/Schedule/Action/State.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQState.hpp>

#include <opm/simulators/flow/Banners.hpp>
#include <opm/simulators/utils/readDeck.hpp>

#if HAVE_DAMARIS
#include <Damaris.h>
#include <opm/simulators/utils/DamarisOutputModule.hpp>
#endif

#if HAVE_HYPRE
#include <HYPRE_config.h>
#include <HYPRE_utilities.h>
#endif

#if HAVE_AMGX
#include <amgx_c.h>
#endif

#include <iostream>
// NOTE: There is no C++ header replacement for these C posix headers (as of C++17)
#include <fcntl.h>  // for open()
#include <unistd.h> // for dup2(), close()

#include <iostream>

namespace Opm {

Main::Main(int argc, char** argv, bool ownMPI)
    : argc_(argc), argv_(argv), ownMPI_(ownMPI)
{
#if HAVE_MPI
    maybeSaveReservoirCouplingSlaveLogFilename_();
#endif
    if (ownMPI_) {
        initMPI();
    }
}

Main::Main(const std::string& filename, bool mpi_init, bool mpi_finalize)
    : mpi_init_{mpi_init}
    , mpi_finalize_{mpi_finalize}
{
    setArgvArgc_(filename);
    initMPI();
}

Main::Main(const std::string& filename,
           std::shared_ptr<EclipseState> eclipseState,
           std::shared_ptr<Schedule> schedule,
           std::shared_ptr<SummaryConfig> summaryConfig,
           bool mpi_init,
           bool mpi_finalize)
    : eclipseState_{std::move(eclipseState)}
    , schedule_{std::move(schedule)}
    , summaryConfig_{std::move(summaryConfig)}
    , mpi_init_{mpi_init}
    , mpi_finalize_{mpi_finalize}
{
    setArgvArgc_(filename);
    initMPI();
}

Main::~Main()
{
#if HAVE_MPI
    if (test_split_comm_) {
        // Cannot use EclGenericVanguard::comm()
        // to get world size here, as it may be
        // a split communication at this point.
        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        if (world_size > 1) {
            MPI_Comm new_comm = FlowGenericVanguard::comm();
            int result;
            MPI_Comm_compare(MPI_COMM_WORLD, new_comm, &result);
            assert(result == MPI_UNEQUAL);
            MPI_Comm_free(&new_comm);
        }
    }
#endif // HAVE_MPI

#if HAVE_HYPRE
    HYPRE_Finalize();
#endif

#if HAVE_AMGX
    AMGX_SAFE_CALL(AMGX_finalize());
#endif

    if (ownMPI_) {
        FlowGenericVanguard::setCommunication(nullptr);
    }

#if HAVE_DAMARIS
    if (enableDamarisOutput_) {
        int err;
        if (isSimulationRank_) {
            err = damaris_stop();
            if (err != DAMARIS_OK) {
                std::cerr << "ERROR: Damaris library produced an error result for damaris_stop()" << std::endl;
            }
        }
        err = damaris_finalize();
        if (err != DAMARIS_OK) {
            std::cerr << "ERROR: Damaris library produced an error result for damaris_finalize()" << std::endl;
        }
    }
#endif // HAVE_DAMARIS

#if HAVE_MPI && !HAVE_DUNE_FEM
    if (ownMPI_ && this->mpi_finalize_) {
        MPI_Finalize();
    }
#endif
}

#if HAVE_MPI
void Main::maybeSaveReservoirCouplingSlaveLogFilename_()
{
    // If first command line argument is "--slave-log-file=<filename>",
    // then redirect stdout and stderr to the specified file.
    if (this->argc_ >= 2) {
        std::string_view arg = this->argv_[1];
        if (arg.substr(0, 17) == "--slave-log-file=") {
            // For now, just save the basename of the filename and we will append the rank
            // to the filename after having called MPI_Init() to avoid all ranks outputting
            // to the same file.
            this->reservoirCouplingSlaveOutputFilename_ = arg.substr(17);
            this->argc_ -= 1;
            char *program_name = this->argv_[0];
            this->argv_ += 1;
            // We assume the "argv" array pointers remain valid (not freed) for the lifetime
            //   of this program, so the following assignment is safe.
            // Overwrite the first argument with the program name, i.e. pretend the program
            //   was called with the same arguments, but without the "--slave-log-file" argument.
            this->argv_[0] = program_name;
        }
    }
}
#endif
#if HAVE_MPI
void Main::maybeRedirectReservoirCouplingSlaveOutput_() {
    if (!this->reservoirCouplingSlaveOutputFilename_.empty()) {
        std::string filename = this->reservoirCouplingSlaveOutputFilename_
                     + "." + std::to_string(FlowGenericVanguard::comm().rank()) + ".log";
        int fd = open(filename.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        if (fd == -1) {
            std::string error_msg = "Slave: Failed to open stdout+stderr file" + filename;
            perror(error_msg.c_str());
            MPI_Abort(MPI_COMM_WORLD, /*status=*/1);
        }
        // Redirect stdout and stderr to the file.
        if (dup2(fd, fileno(stdout)) == -1 || dup2(fileno(stdout), fileno(stderr)) == -1) {
            std::string error_msg = "Slave: Failed to redirect stdout+stderr to " + filename;
            perror(error_msg.c_str());
            close(fd);
            MPI_Abort(MPI_COMM_WORLD, /*status=*/1);
        }
        close(fd);
    }
}
#endif

void Main::setArgvArgc_(const std::string& filename)
{
    this->deckFilename_ = filename;
    this->flowProgName_ = "flow";

    this->argc_ = 2;
    this->saveArgs_[0] = const_cast<char *>(this->flowProgName_.c_str());
    this->saveArgs_[1] = const_cast<char *>(this->deckFilename_.c_str());

    // Note: argv[argc] must exist and be nullptr
    assert ((sizeof this->saveArgs_) > (this->argc_ * sizeof this->saveArgs_[0]));
    this->saveArgs_[this->argc_] = nullptr;

    this->argv_ = this->saveArgs_;
}

void Main::initMPI()
{
#if HAVE_DUNE_FEM
    // The instance() method already checks if MPI has been initialized so we may
    // not need to check mpi_init_ here.
    if (this->mpi_init_) {
        Dune::MPIHelper::instance(argc_, argv_);
    }
#elif HAVE_MPI
    if (this->mpi_init_) {
        MPI_Init(&argc_, &argv_);
    }
#endif
    FlowGenericVanguard::setCommunication(std::make_unique<Parallel::Communication>());

    handleTestSplitCommunicatorCmdLine_();

#if HAVE_MPI
    maybeRedirectReservoirCouplingSlaveOutput_();
    if (test_split_comm_ && FlowGenericVanguard::comm().size() > 1) {
        int world_rank = FlowGenericVanguard::comm().rank();
        int color = (world_rank == 0);
        MPI_Comm new_comm;
        MPI_Comm_split(FlowGenericVanguard::comm(), color, world_rank, &new_comm);
        isSimulationRank_ = (world_rank > 0);
        FlowGenericVanguard::setCommunication(std::make_unique<Parallel::Communication>(new_comm));
    }

#if HAVE_CUDA
    Opm::gpuistl::setDevice();
#endif

#endif // HAVE_MPI

#if HAVE_HYPRE
#if HYPRE_RELEASE_NUMBER >= 22900
    HYPRE_Initialize();
#else
    HYPRE_Init();
#endif
#endif

#if HAVE_AMGX
    AMGX_SAFE_CALL(AMGX_initialize());
#endif
}

void Main::handleVersionCmdLine_(int argc, char** argv,
                                 std::string_view moduleVersionName)
{
    auto pos = std::find_if(argv, argv + argc,
        [](const char* arg)
    {
        return std::strcmp(arg, "--version") == 0;
    });

    if (pos != argv + argc) {
        std::cout << "flow " << moduleVersionName << std::endl;
        std::exit(EXIT_SUCCESS);
    }
}

void Main::handleTestSplitCommunicatorCmdLine_()
{
    if (argc_ >= 2 && std::strcmp(argv_[1], "--test-split-communicator=true") == 0) {
        test_split_comm_ = true;
        --argc_;             // We have one less argument.
        argv_[1] = argv_[0]; // What used to be the first proper argument now becomes the command argument.
        ++argv_;             // Pretend this is what it always was.
    }
}

void Main::readDeck(const std::string& deckFilename,
                    const std::string& outputDir,
                    const std::string& outputMode,
                    const bool init_from_restart_file,
                    const bool allRanksDbgPrtLog,
                    const std::string& parsingStrictness,
                    const std::string& actionParsingStrictness,
                    const std::string& inputSkipMode,
                    const bool keepKeywords,
                    const std::size_t numThreads,
                    const int output_param,
                    const bool slaveMode,
                    const std::string& parameters,
                    std::string_view moduleVersion,
                    std::string_view compileTimestamp)
{
    auto omode = setupLogging(FlowGenericVanguard::comm(),
                              deckFilename,
                              outputDir,
                              outputMode,
                              outputCout_, "STDOUT_LOGGER", allRanksDbgPrtLog);

    if (outputCout_) {
        printPRTHeader(FlowGenericVanguard::comm().size(), numThreads,
                       parameters, moduleVersion, compileTimestamp);
        OpmLog::info("Reading deck file '" + deckFilename + "'");
    }

    std::optional<int> outputInterval;
    if (output_param >= 0)
        outputInterval = output_param;

    Opm::readDeck(FlowGenericVanguard::comm(),
                  deckFilename,
                  eclipseState_,
                  schedule_,
                  udqState_,
                  actionState_,
                  wtestState_,
                  summaryConfig_,
                  std::make_shared<Python>(),
                  parsingStrictness,
                  actionParsingStrictness,
                  inputSkipMode,
                  init_from_restart_file,
                  outputCout_,
                  keepKeywords,
                  outputInterval,
                  slaveMode);

    verifyValidCellGeometry(FlowGenericVanguard::comm(), *this->eclipseState_);

    outputFiles_ = (omode != FileOutputMode::OUTPUT_NONE);
}

void Main::setupVanguard()
{
    FlowGenericVanguard::modelParams_.setupTime_ = this->setupTime_;
    FlowGenericVanguard::modelParams_.actionState_ = std::move(this->actionState_);
    FlowGenericVanguard::modelParams_.eclSchedule_ = this->schedule_;
    FlowGenericVanguard::modelParams_.eclState_ = this->eclipseState_;
    FlowGenericVanguard::modelParams_.eclSummaryConfig_ = this->summaryConfig_;
    FlowGenericVanguard::modelParams_.udqState_ = std::move(udqState_);
    FlowGenericVanguard::modelParams_.wtestState_ = std::move(wtestState_);
}

#if HAVE_DAMARIS
void Main::setupDamaris(const std::string& outputDir )
{
    if (!outputDir.empty()) {
        ensureOutputDirExists(outputDir);
    }

    //const auto find_replace_map;
    //const auto find_replace_map = Opm::DamarisOutput::DamarisKeywords<PreTypeTag>(EclGenericVanguard::comm(), outputDir);
    std::map<std::string, std::string> find_replace_map;
    find_replace_map = DamarisOutput::getDamarisKeywords(FlowGenericVanguard::comm(), outputDir);

    // By default EnableDamarisOutputCollective is true so all simulation results will
    // be written into one single file for each iteration using Parallel HDF5.
    // If set to false, FilePerCore mode is used in Damaris, then simulation results in each
    // node are aggregated by dedicated Damaris cores and stored to separate files per Damaris core.
    // Irrespective of mode, output is written asynchronously at the end of each timestep.
    // Using the ModifyModel class to set the XML file for Damaris.
    DamarisOutput::initializeDamaris(FlowGenericVanguard::comm(),
                                     FlowGenericVanguard::comm().rank(),
                                     find_replace_map);
    int is_client;
    MPI_Comm new_comm;
    // damaris_start() is where the Damaris Server ranks will block, until damaris_stop()
    // is called from the client ranks
    int err = damaris_start(&is_client);
    isSimulationRank_ = (is_client > 0);
    if (isSimulationRank_ && err == DAMARIS_OK) {
        damaris_client_comm_get(&new_comm);
        FlowGenericVanguard::setCommunication(std::make_unique<Parallel::Communication>(new_comm));
    }

    if (err != DAMARIS_OK) {
        OPM_THROW(std::runtime_error, "Failed to configure Damaris: " + std::to_string(err));
    }
}
#endif

} // namespace Opm
