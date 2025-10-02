/*
  Copyright 2013, 2014, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 IRIS AS
  Copyright 2014 STATOIL ASA.
  Copyright 2023 Inria

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
#ifndef OPM_MAIN_HEADER_INCLUDED
#define OPM_MAIN_HEADER_INCLUDED

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hpp>

#include <opm/simulators/flow/Banners.hpp>
#include <opm/simulators/flow/FlowMain.hpp>
#include <opm/simulators/flow/FlowProblemBlackoil.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

#if HAVE_MPI
#include <opm/simulators/utils/ParallelEclipseState.hpp>
#endif

#if HAVE_CUDA
#include <opm/simulators/linalg/gpuistl/device_management.hpp>
#endif

#if HAVE_DAMARIS
#include <opm/simulators/utils/DamarisKeywords.hpp>
#endif

#include <cassert>
#include <charconv>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

namespace Opm::Properties {

// this is a dummy type tag that is used to setup the parameters before the actual
// simulator.
namespace TTag {
struct FlowEarlyBird {
    using InheritsFrom = std::tuple<FlowProblem>;
};
}

} // namespace Opm::Properties

namespace Opm {

namespace Action { class State; }
class UDQState;
class WellTestState;

// ----------------- Main program -----------------
template <class TypeTag>
int flowMain(int argc, char** argv, bool outputCout, bool outputFiles)
{
    // we always want to use the default locale, and thus spare us the trouble
    // with incorrect locale settings.
    resetLocale();

    FlowMain<TypeTag> mainfunc(argc, argv, outputCout, outputFiles);
    return mainfunc.execute();
}

// ----------------- Main class -----------------
//   For now, we will either be instantiated from main() in flow.cpp,
//   or from a Python pybind11 module..
// NOTE (March 2020): When used from a pybind11 module, we do not neccessarily
//   want to run the whole simulation by calling run(), it is also
//   useful to just run one report step at a time. According to these different
//   usage scenarios, we refactored the original run() in flow.cpp into this class.
class Main
{
public:
    Main(int argc, char** argv, bool ownMPI = true);

    // This constructor can be called from Python
    explicit Main(const std::string& filename, bool mpi_init = true, bool mpi_finalize = true);

    // This constructor can be called from Python when Python has
    // already parsed a deck
    Main(const std::string& filename,
         std::shared_ptr<EclipseState> eclipseState,
         std::shared_ptr<Schedule> schedule,
         std::shared_ptr<SummaryConfig> summaryConfig,
         bool mpi_init = true,
         bool mpi_finalize = true);

    ~Main();

    void setArgvArgc_(const std::string& filename);
    void maybeSaveReservoirCouplingSlaveLogFilename_();
    void maybeRedirectReservoirCouplingSlaveOutput_();
    void initMPI();

    /// Run simulation.
    ///
    /// Selects an appropriate simulator based on runtime information in the
    /// input deck.
    ///
    /// \return Simulation's status/exit code.
    int runDynamic()
    {
        int exitCode = EXIT_SUCCESS;
        if (initialize_<Properties::TTag::FlowEarlyBird>(exitCode)) {
            Parameters::reset();
            if (isSimulationRank_) {
                return this->dispatchDynamic_();
            }
        }

        return exitCode;
    }

    /// Run simulation.
    ///
    /// Uses staticially configured simulator defined at call site.
    ///
    /// \tparam TypeTag Simulation type's statically configured properties.
    ///
    /// \return Simulation's status/exit code.
    template <class TypeTag>
    int runStatic()
    {
        int exitCode = EXIT_SUCCESS;
        if (initialize_<TypeTag>(exitCode)) {
            if (isSimulationRank_) {
                return this->dispatchStatic_<TypeTag>();
            }
        }

        return exitCode;
    }

    //! \brief Used for test_outputdir.
    int justInitialize()
    {
        int exitCode = EXIT_SUCCESS;
        initialize_<Properties::TTag::FlowEarlyBird>(exitCode);
        return exitCode;
    }

protected:
    /// \brief Initialize
    /// \param exitCode The exitCode of the program.
    /// \param keepKeywords Keep Schedule keywords even if there are no actions
    ///
    /// \return Whether to actually run the simulator. I.e. true if
    /// parsing of command line was successful and no --help,
    /// --print-properties, or --print-parameters have been found.
    template <class TypeTagEarlyBird>
    bool initialize_(int& exitCode, bool keepKeywords = false)
    {
        Dune::Timer externalSetupTimer;
        externalSetupTimer.start();

        handleVersionCmdLine_(argc_, argv_, Opm::moduleVersionName());

        // we always want to use the default locale, and thus spare us the trouble
        // with incorrect locale settings.
        resetLocale();

        // this is a work-around for a catch 22: we do not know what code path to use without
        // parsing the deck, but we don't know the deck without having access to the
        // parameters and this requires to know the type tag to be used. To solve this, we
        // use a type tag just for parsing the parameters before we instantiate the actual
        // simulator object. (Which parses the parameters again, but since this is done in an
        // identical manner it does not matter.)
        typedef TypeTagEarlyBird PreTypeTag;
        using PreProblem = GetPropType<PreTypeTag, Properties::Problem>;

        PreProblem::setBriefDescription("Flow, an advanced reservoir simulator for ECL-decks provided by the Open Porous Media project.");
        int status = FlowMain<PreTypeTag>::setupParameters_(argc_, argv_, FlowGenericVanguard::comm());
        if (status != 0) {
            // if setupParameters_ returns a value smaller than 0, there was no error, but
            // the program should abort. This is the case e.g. for the --help and the
            // --print-properties parameters.
#if HAVE_MPI
            if (status >= 0)
                MPI_Abort(MPI_COMM_WORLD, status);
#endif
            exitCode = (status > 0) ? status : EXIT_SUCCESS;
            return false; //  Whether to run the simulator
        }

        OpmLog::setDebugVerbosityLevel(Parameters::Get<Parameters::DebugVerbosityLevel>());

        std::string deckFilename;
        std::string outputDir;
        if ( eclipseState_ ) {
            deckFilename = eclipseState_->getIOConfig().fullBasePath();
            outputDir = eclipseState_->getIOConfig().getOutputDir();
        }
        else {
            deckFilename = Parameters::Get<Parameters::EclDeckFileName>();
            outputDir = Parameters::Get<Parameters::OutputDir>();
        }

#if HAVE_DAMARIS
        enableDamarisOutput_ = Parameters::Get<Parameters::EnableDamarisOutput>();

        // Reset to false as we cannot use Damaris if there is only one rank.
        if ((enableDamarisOutput_ == true) && (FlowGenericVanguard::comm().size() == 1)) {
            std::string msg ;
            msg = "\nUse of Damaris (command line argument --enable-damaris-output=true) has been disabled for run with only one rank.\n" ;
            OpmLog::warning(msg);
            enableDamarisOutput_ = false ;
        }

        if (enableDamarisOutput_) {
            // Deal with empty (defaulted) output dir, should be deck dir
            auto damarisOutputDir = outputDir;
            if (outputDir.empty()) {
                auto odir = std::filesystem::path{deckFilename}.parent_path();
                if (odir.empty()) {
                    damarisOutputDir = ".";
                } else {
                    damarisOutputDir = odir.generic_string();
                }
            }
            // Damaris server ranks will block here until damaris_stop() is called by client ranks
            this->setupDamaris(damarisOutputDir);
        }
#endif // HAVE_DAMARIS

        // Guard for when the Damaris core(s) return from damaris_start()
        // which happens when damaris_stop() is called in main simulation
        if (!isSimulationRank_) {
            exitCode = EXIT_SUCCESS;
            return true;
        }

        int mpiRank = FlowGenericVanguard::comm().rank();
        outputCout_ = false;
        if (mpiRank == 0)
            outputCout_ = Parameters::Get<Parameters::EnableTerminalOutput>();

        if (deckFilename.empty()) {
            if (mpiRank == 0) {
                std::cerr << "No input case given. Try '--help' for a usage description.\n";
            }
            exitCode = EXIT_FAILURE;
            return false;
        }

        using PreVanguard = GetPropType<PreTypeTag, Properties::Vanguard>;
        try {
            deckFilename = PreVanguard::canonicalDeckPath(deckFilename);
        }
        catch (const std::exception& e) {
            if ( mpiRank == 0 ) {
                std::cerr << "Exception received: " << e.what() << ". Try '--help' for a usage description.\n";
            }
#if HAVE_MPI
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
            exitCode = EXIT_FAILURE;
            return false;
        }

        std::string cmdline_params;
        if (outputCout_) {
            printFlowBanner(FlowGenericVanguard::comm().size(),
                            getNumThreads(),
                            Opm::moduleVersionName());
            std::ostringstream str;
            Parameters::printValues(str);
            cmdline_params = str.str();
        }

        // Create Deck and EclipseState.
        try {
            this->readDeck(deckFilename,
                           outputDir,
                           Parameters::Get<Parameters::OutputMode>(),
                           !Parameters::Get<Parameters::SchedRestart>(),
                           Parameters::Get<Parameters::EnableLoggingFalloutWarning>(),
                           Parameters::Get<Parameters::ParsingStrictness>(),
                           Parameters::Get<Parameters::ActionParsingStrictness>(),
                           Parameters::Get<Parameters::InputSkipMode>(),
                           keepKeywords,
                           getNumThreads(),
                           Parameters::Get<Parameters::EclOutputInterval>(),
                           Parameters::Get<Parameters::Slave>(),
                           cmdline_params,
                           Opm::moduleVersion(),
                           Opm::compileTimestamp());
            setupTime_ = externalSetupTimer.elapsed();
        }
        catch (const std::invalid_argument& e)
        {
            if (outputCout_) {
                std::cerr << "Failed to create valid EclipseState object." << std::endl;
                std::cerr << "Exception caught: " << e.what() << std::endl;
            }
#if HAVE_MPI
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
            exitCode = EXIT_FAILURE;
            return false;
        }

#if HAVE_CUDA
    Opm::gpuistl::printDevice();
#endif

        exitCode = EXIT_SUCCESS;
        return true;
    }

    void setupVanguard();

private:
    // This function is an extreme special case, if the program has been invoked
    // *exactly* as:
    //
    //    flow   --version
    //
    // the call is intercepted by this function which will print "flow $version"
    // on stdout and exit(0).
    void handleVersionCmdLine_(int argc, char** argv,
                               std::string_view moduleVersionName);

    // This function is a special case, if the program has been invoked
    // with the argument "--test-split-communicator=true" as the FIRST
    // argument, it will be removed from the argument list and we set the
    // test_split_comm_ flag to true.
    // Note: initializing the parameter system before MPI could make this
    // use the parameter system instead.
    void handleTestSplitCommunicatorCmdLine_();

    /// Dispatch to actual simulation functions based on input deck's setup.
    ///
    /// Called from runDynamic()
    ///
    /// \return Simulation's status/exit code.
    int dispatchDynamic_();

    /// Dispatch to actual simulation function based on statically
    /// configured simulation type.
    ///
    /// Called from runStatic()
    ///
    /// \tparam TypeTag Simulation type's statically configured properties.
    ///
    /// \return Simulation's status/exit code.
    template <class TypeTag>
    int dispatchStatic_()
    {
        this->setupVanguard();
        return flowMain<TypeTag>(argc_, argv_, outputCout_, outputFiles_);
    }

    /// Run a simulation with MICP effects.
    ///
    /// Called from dispatchDynamic_()
    ///
    /// \param[in] phases Run's active phases.  Needed to determine whether
    /// or not the run's phase setup is supported.
    ///
    /// \return Simulation's status/exit code.
    int runMICP(const Phases& phases);

    /// Run a simulation with two active phases.
    ///
    /// Called from dispatchDynamic_()
    ///
    /// \param[in] phases Run's active phases.  Needed to determine whether
    /// or not the run's phase setup is supported.
    ///
    /// \return Simulation's status/exit code.
    int runTwoPhase(const Phases& phases);

    /// Run a simulation with Biofilm effects.
    ///
    /// Called from dispatchDynamic_()
    ///
    /// \param[in] phases Run's active phases.  Needed to determine whether
    /// or not the run's phase setup is supported.
    ///
    /// \return Simulation's status/exit code.
    int runBiofilm(const Phases& phases);

    /// Run a simulation with polymers.
    ///
    /// Called from dispatchDynamic_()
    ///
    /// \param[in] phases Run's active phases.  Needed to determine whether
    /// or not the run's phase setup is supported.
    ///
    /// \return Simulation's status/exit code.
    int runPolymer(const Phases& phases);

    /// Run a simulation with foam.
    ///
    /// Called from dispatchDynamic_()
    ///
    /// \return Simulation's status/exit code.
    int runFoam();

    /// Run a single phase, water-only simulation.
    ///
    /// Called from dispatchDynamic_()
    ///
    /// \param[in] phases Run's active phases.  Needed to determine whether
    /// or not the run's phase setup is supported.
    ///
    /// \return Simulation's status/exit code.
    int runWaterOnly(const Phases& phases);

    /// Run a single phase, water-only simulation with themal/energy effects.
    ///
    /// Called from dispatchDynamic_()
    ///
    /// \param[in] phases Run's active phases.  Needed to determine whether
    /// or not the run's phase setup is supported.
    ///
    /// \return Simulation's status/exit code.
    int runWaterOnlyEnergy(const Phases& phases);

    /// Run a simulation with brine, typically in a CCS workflow
    ///
    /// Called from dispatchDynamic_()
    ///
    /// \param[in] phases Run's active phases.  Needed to determine whether
    /// or not the run's phase setup is supported.
    ///
    /// \return Simulation's status/exit code.
    int runBrine(const Phases& phases);

    /// Run a simulation with solvents.
    ///
    /// Called from dispatchDynamic_()
    ///
    /// \param[in] phases Run's active phases.  Needed to determine whether
    /// or not the run's phase setup is supported.
    ///
    /// \return Simulation's status/exit code.
    int runSolvent(const Phases& phases);

    /// Run a simulation with the extended black-oil model.
    ///
    /// Called from dispatchDynamic_()
    ///
    /// \return Simulation's status/exit code.
    int runExtendedBlackOil();

    /// Run a three-phase simulation with thermal effects.
    ///
    /// Called from dispatchDynamic_()
    ///
    /// \param[in] phases Run's active phases.  Needed to determine whether
    /// or not the run's phase setup is supported.
    ///
    /// \return Simulation's status/exit code.
    int runThermal(const Phases& phases);

    /// Run a regular three-phase simulation without thermal effects.
    ///
    /// Called from dispatchDynamic_()
    ///
    /// \return Simulation's status/exit code.
    int runBlackOil();

    void readDeck(const std::string& deckFilename,
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
                  std::string_view compileTimestamp);

    static int getNumThreads()
    {
#ifdef _OPENMP
        return omp_get_max_threads();
#else
        return 1;
#endif
    }

#if HAVE_DAMARIS
    void setupDamaris(const std::string& outputDir);
#endif

protected:
    int argc_{0};
    char** argv_{nullptr};
    bool outputCout_{false};
    bool outputFiles_{false};

private:
    bool ownMPI_{true}; //!< True if we "own" MPI and should init / finalize
    double setupTime_{0.0};
    std::string deckFilename_{};
    std::string flowProgName_{};
    char *saveArgs_[3]{nullptr};
    std::unique_ptr<UDQState> udqState_{};
    std::unique_ptr<Action::State> actionState_{};
    std::unique_ptr<WellTestState> wtestState_{};

    // These variables may be owned by both Python and the simulator
    std::shared_ptr<EclipseState> eclipseState_{};
    std::shared_ptr<Schedule> schedule_{};
    std::shared_ptr<SummaryConfig> summaryConfig_{};
    bool mpi_init_{true}; //!< True if MPI_Init should be called
    bool mpi_finalize_{true}; //!< True if MPI_Finalize should be called

    // To demonstrate run with non_world_comm
    bool test_split_comm_ = false;
    bool isSimulationRank_ = true;
#if HAVE_MPI
    std::string reservoirCouplingSlaveOutputFilename_{};
#endif
#if HAVE_DAMARIS
    bool enableDamarisOutput_ = false;
#endif
};

} // namespace Opm

#endif // OPM_MAIN_HEADER_INCLUDED
