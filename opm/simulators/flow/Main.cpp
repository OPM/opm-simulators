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
#include <opm/simulators/utils/DamarisOutputModule.hpp>
#endif

#if HAVE_CUDA
#include <opm/simulators/linalg/cuistl/set_device.hpp>
#endif

namespace Opm {

Main::Main(int argc, char** argv)
    : argc_(argc), argv_(argv)
{
    initMPI();
}

Main::Main(const std::string& filename)
{
    setArgvArgc_(filename);
    initMPI();
}

Main::Main(const std::string& filename,
           std::shared_ptr<EclipseState> eclipseState,
           std::shared_ptr<Schedule> schedule,
           std::shared_ptr<SummaryConfig> summaryConfig)
    : eclipseState_{std::move(eclipseState)}
    , schedule_{std::move(schedule)}
    , summaryConfig_{std::move(summaryConfig)}
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
            MPI_Comm new_comm = EclGenericVanguard::comm();
            int result;
            MPI_Comm_compare(MPI_COMM_WORLD, new_comm, &result);
            assert(result == MPI_UNEQUAL);
            MPI_Comm_free(&new_comm);
        }
    }
#endif // HAVE_MPI

    EclGenericVanguard::setCommunication(nullptr);

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
    MPI_Finalize();
#endif
}

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
    Dune::Fem::MPIManager::initialize(argc_, argv_);
#elif HAVE_MPI
    MPI_Init(&argc_, &argv_);
#endif
    EclGenericVanguard::setCommunication(std::make_unique<Parallel::Communication>());

    handleTestSplitCommunicatorCmdLine_();

#if HAVE_MPI
    if (test_split_comm_ && EclGenericVanguard::comm().size() > 1) {
        int world_rank = EclGenericVanguard::comm().rank();
        int color = (world_rank == 0);
        MPI_Comm new_comm;
        MPI_Comm_split(EclGenericVanguard::comm(), color, world_rank, &new_comm);
        isSimulationRank_ = (world_rank > 0);
        EclGenericVanguard::setCommunication(std::make_unique<Parallel::Communication>(new_comm));
    }

#if HAVE_CUDA
    Opm::cuistl::setDevice(EclGenericVanguard::comm().rank(), EclGenericVanguard::comm().size());
#endif

#endif // HAVE_MPI
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
                    const int mpiRank,
                    const std::size_t numThreads,
                    const int output_param,
                    const std::string& parameters,
                    std::string_view moduleVersion,
                    std::string_view compileTimestamp)
{
    auto omode = setupLogging(mpiRank,
                              deckFilename,
                              outputDir,
                              outputMode,
                              outputCout_, "STDOUT_LOGGER", allRanksDbgPrtLog);

    if (outputCout_) {
        printPRTHeader(EclGenericVanguard::comm().size(), numThreads,
                       parameters, moduleVersion, compileTimestamp);
        OpmLog::info("Reading deck file '" + deckFilename + "'");
    }

    std::optional<int> outputInterval;
    if (output_param >= 0)
        outputInterval = output_param;

    Opm::readDeck(EclGenericVanguard::comm(),
                  deckFilename,
                  eclipseState_,
                  schedule_,
                  udqState_,
                  actionState_,
                  wtestState_,
                  summaryConfig_,
                  std::make_shared<Python>(),
                  parsingStrictness,
                  init_from_restart_file,
                  outputCout_,
                  outputInterval);

    verifyValidCellGeometry(EclGenericVanguard::comm(), *this->eclipseState_);

    outputFiles_ = (omode != FileOutputMode::OUTPUT_NONE);
}

void Main::setupVanguard()
{
    EclGenericVanguard::modelParams_.setupTime_ = this->setupTime_;
    EclGenericVanguard::modelParams_.actionState_ = std::move(this->actionState_);
    EclGenericVanguard::modelParams_.eclSchedule_ = this->schedule_;
    EclGenericVanguard::modelParams_.eclState_ = this->eclipseState_;
    EclGenericVanguard::modelParams_.eclSummaryConfig_ = this->summaryConfig_;
    EclGenericVanguard::modelParams_.udqState_ = std::move(udqState_);
    EclGenericVanguard::modelParams_.wtestState_ = std::move(wtestState_);
}

#if HAVE_DAMARIS
void Main::setupDamaris(const std::string& outputDir,
                        const bool enableDamarisOutputCollective)
{
    if (!outputDir.empty()) {
        ensureOutputDirExists(outputDir);
    }

    // By default EnableDamarisOutputCollective is true so all simulation results will
    // be written into one single file for each iteration using Parallel HDF5.
    // It set to false, FilePerCore mode is used in Damaris, then simulation results in each
    // node are aggregated by dedicated Damaris cores and stored to separate files per Damaris core.
    // Irrespective of mode, output is written asynchronously at the end of each timestep.
    // Using the ModifyModel class to set the XML file for Damaris.
    DamarisOutput::initializeDamaris(EclGenericVanguard::comm(), EclGenericVanguard::comm().rank(), outputDir, enableDamarisOutputCollective);
    int is_client;
    MPI_Comm new_comm;
    int err = damaris_start(&is_client);
    isSimulationRank_ = (is_client > 0);
    if (isSimulationRank_ && err == DAMARIS_OK) {
        damaris_client_comm_get(&new_comm);
        EclGenericVanguard::setCommunication(std::make_unique<Parallel::Communication>(new_comm));
    }

    if (err != DAMARIS_OK) {
        OPM_THROW(std::runtime_error, "Failed to configure Damaris: " + std::to_string(err));
    }
}
#endif

} // namespace Opm
