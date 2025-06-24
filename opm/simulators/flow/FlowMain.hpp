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
#ifndef OPM_FLOW_MAIN_HEADER_INCLUDED
#define OPM_FLOW_MAIN_HEADER_INCLUDED

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/IOConfig/IOConfig.hpp>
#include <opm/input/eclipse/EclipseState/InitConfig/InitConfig.hpp>

#include <opm/models/utils/start.hh>

#include <opm/simulators/flow/Banners.hpp>
#include <opm/simulators/flow/FlowUtils.hpp>
#include <opm/simulators/flow/NlddReporting.hpp>
#include <opm/simulators/flow/SimulatorFullyImplicitBlackoil.hpp>

#if HAVE_MPI
#define RESERVOIR_COUPLING_ENABLED
#endif
#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <charconv>
#include <cstddef>
#include <memory>

namespace Opm::Parameters {

// Do not merge parallel output files or warn about them
struct EnableLoggingFalloutWarning { static constexpr bool value = false; };
struct OutputInterval { static constexpr int value = 1; };
// Set global debug verbosity level
struct DebugVerbosityLevel { static constexpr int value = 1; };
} // namespace Opm::Parameters

namespace Opm {

    class Deck;

    // The FlowMain class is the black-oil simulator.
    template <class TypeTag>
    class FlowMain
    {
    public:
        using MaterialLawManager = typename GetProp<TypeTag, Properties::MaterialLaw>::EclMaterialLawManager;
        using ModelSimulator = GetPropType<TypeTag, Properties::Simulator>;
        using Grid = GetPropType<TypeTag, Properties::Grid>;
        using GridView = GetPropType<TypeTag, Properties::GridView>;
        using Problem = GetPropType<TypeTag, Properties::Problem>;
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

        using Simulator = SimulatorFullyImplicitBlackoil<TypeTag>;

        FlowMain(int argc, char **argv, bool output_cout, bool output_files )
            : argc_{argc}, argv_{argv},
              output_cout_{output_cout}, output_files_{output_files}
        {

        }

        // Read the command line parameters. Throws an exception if something goes wrong.
        static int setupParameters_(int argc, char** argv, Parallel::Communication comm)
        {
            if (!Parameters::IsRegistrationOpen()) {
                // We have already successfully run setupParameters_().
                // For the dynamically chosen runs (as from the main flow
                // executable) we must run this function again with the
                // real typetag to be used, as the first time was with the
                // "FlowEarlyBird" typetag. However, for the static ones (such
                // as 'flow_onephase_energy') it has already been run with the
                // correct typetag.
                return EXIT_SUCCESS;
            }
            // register the flow specific parameters
            Parameters::Register<Parameters::OutputInterval>
                ("Specify the number of report steps between two consecutive writes of restart data");
            Parameters::Register<Parameters::EnableLoggingFalloutWarning>
                ("Developer option to see whether logging was on non-root processors. "
                 "In that case it will be appended to the *.DBG or *.PRT files");
            Parameters::Register<Parameters::DebugVerbosityLevel>
                ("Set debug verbosity level globally. Default is 1, increasing values give additional output and 0 disables most messages to the .DBG file");

            // register the base parameters
            registerAllParameters_<TypeTag>(/*finalizeRegistration=*/false);

            Simulator::registerParameters();

            detail::hideUnusedParameters<Scalar>();

            Parameters::endRegistration();

            int mpiRank = comm.rank();

            // read in the command line parameters
            int status = ::Opm::setupParameters_<TypeTag>(argc,
                                                          const_cast<const char**>(argv),
                                                          /*doRegistration=*/false,
                                                          /*allowUnused=*/true,
                                                          /*handleHelp=*/(mpiRank==0),
                                                          mpiRank);
            if (status == 0) {

                // deal with unknown parameters.

                int unknownKeyWords = 0;
                if (mpiRank == 0) {
                    unknownKeyWords = Parameters::printUnused(std::cerr);
                }
                int globalUnknownKeyWords = comm.sum(unknownKeyWords);
                unknownKeyWords = globalUnknownKeyWords;
                if ( unknownKeyWords )
                {
                    if ( mpiRank == 0 )
                    {
                        std::string msg = "Aborting simulation due to unknown "
                            "parameters. Please query \"flow --help\" for "
                            "supported command line parameters.";
                        if (OpmLog::hasBackend("STREAMLOG"))
                        {
                            OpmLog::error(msg);
                        }
                        else {
                            std::cerr << msg << std::endl;
                        }
                    }
                    return EXIT_FAILURE;
                }

                // deal with --print-parameters and unknown parameters.
                if (Parameters::Get<Parameters::PrintParameters>() == 1) {
                    if (mpiRank == 0) {
                        Parameters::printValues(std::cout);
                    }
                    return -1;
                }
            }

            // set the maximum limit on OMP threads
            setMaxThreads();

            return status;
        }

        /// This is the main function of Flow.  It runs a complete simulation with the
        /// given grid and simulator classes, based on the user-specified command-line
        /// input.
        int execute()
        {
            return execute_(&FlowMain::runSimulator, /*cleanup=*/true);
        }

        int executeInitStep()
        {
            return execute_(&FlowMain::runSimulatorInit, /*cleanup=*/false);
        }

        // Returns true unless "EXIT" was encountered in the schedule
        //   section of the input datafile.
        int executeStep()
        {
            return simulator_->runStep(*simtimer_);
        }

        // Called from Python to cleanup after having executed the last
        // executeStep()
        int executeStepsCleanup()
        {
            SimulatorReport report = simulator_->finalize();
            runSimulatorAfterSim_(report);
            return report.success.exit_status;
        }

        ModelSimulator* getSimulatorPtr()
        {
            return modelSimulator_.get();
        }

        SimulatorTimer* getSimTimer()
        {
            return simtimer_.get();
        }

        /// Get the size of the previous report step
        double getPreviousReportStepSize()
        {
            return simtimer_->stepLengthTaken();
        }

    private:
        // called by execute() or executeInitStep()
        int execute_(int (FlowMain::* runOrInitFunc)(), bool cleanup)
        {
            auto logger = [this](const std::exception& e, const std::string& message_start) {
                std::ostringstream message;
                message  << message_start << e.what();

                if (this->output_cout_) {
                    // in some cases exceptions are thrown before the logging system is set
                    // up.
                    if (OpmLog::hasBackend("STREAMLOG")) {
                        OpmLog::error(message.str());
                    }
                    else {
                        std::cout << message.str() << "\n";
                    }
                }
                detail::checkAllMPIProcesses();
                return EXIT_FAILURE;
            };

            try {
                // deal with some administrative boilerplate

                Dune::Timer setupTimerAfterReadingDeck;
                setupTimerAfterReadingDeck.start();

                int status = setupParameters_(this->argc_, this->argv_, FlowGenericVanguard::comm());
                if (status) {
                    return status;
                }

                setupParallelism();
                setupModelSimulator();
                createSimulator();

                this->deck_read_time_ = modelSimulator_->vanguard().setupTime();
                this->total_setup_time_ = setupTimerAfterReadingDeck.elapsed() + this->deck_read_time_;

                // if run, do the actual work, else just initialize
                int exitCode = (this->*runOrInitFunc)();
                if (cleanup) {
                    executeCleanup_();
                }
                return exitCode;
            }
            catch (const TimeSteppingBreakdown& e) {
                auto exitCode = logger(e, "Simulation aborted: ");
                executeCleanup_();
                return exitCode;
            }
            catch (const std::exception& e) {
                auto exitCode =  logger(e, "Simulation aborted as program threw an unexpected exception: ");
                executeCleanup_();
                return exitCode;
            }
        }

        void executeCleanup_() {
            // clean up
            mergeParallelLogFiles();
        }

    protected:
        void setupParallelism()
        {
            // determine the rank of the current process and the number of processes
            // involved in the simulation. MPI must have already been initialized
            // here. (yes, the name of this method is misleading.)
            auto comm = FlowGenericVanguard::comm();
            mpi_rank_ = comm.rank();
            mpi_size_ = comm.size();

            setMaxThreads();
        }

        static void setMaxThreads()
        {
#if _OPENMP
            // If openMP is available, default to 2 threads per process unless
            // OMP_NUM_THREADS is set or command line --threads-per-process used.
            // Issue a warning if both OMP_NUM_THREADS and --threads-per-process are set,
            // but let the environment variable take precedence.
            constexpr int default_threads = 2;
            const bool isSet = Parameters::IsSet<Parameters::ThreadsPerProcess>();
            const int requested_threads = Parameters::Get<Parameters::ThreadsPerProcess>();
            int threads = requested_threads > 0 ? requested_threads : default_threads;

            const char* env_var = getenv("OMP_NUM_THREADS");
            if (env_var) {
                int omp_num_threads = -1;
                auto result = std::from_chars(env_var, env_var + std::strlen(env_var), omp_num_threads);
                if (result.ec == std::errc() && omp_num_threads > 0) {
                    // Set threads to omp_num_threads if it was successfully parsed and is positive
                    threads = omp_num_threads;
                    if (isSet) {
                        OpmLog::warning("Environment variable OMP_NUM_THREADS takes precedence over the --threads-per-process cmdline argument.");
                    }
                } else {
                    OpmLog::warning("Invalid value for OMP_NUM_THREADS environment variable.");
                }
            }

            // Requesting -1 thread will let OMP automatically deduce the number
            // but setting OMP_NUM_THREADS takes precedence.
            if (env_var || !(isSet && requested_threads == -1)) {
                omp_set_num_threads(threads);
            }
#endif

            using TM = GetPropType<TypeTag, Properties::ThreadManager>;
            TM::init(false);
        }

        void mergeParallelLogFiles()
        {
            // force closing of all log files.
            OpmLog::removeAllBackends();

            if (mpi_rank_ != 0 || mpi_size_ < 2 || !this->output_files_ || !modelSimulator_) {
                return;
            }

            detail::mergeParallelLogFiles(eclState().getIOConfig().getOutputDir(),
                                          Parameters::Get<Parameters::EclDeckFileName>(),
                                          Parameters::Get<Parameters::EnableLoggingFalloutWarning>());
        }

        void setupModelSimulator()
        {
            modelSimulator_ = std::make_unique<ModelSimulator>(FlowGenericVanguard::comm(), /*verbose=*/false);
            modelSimulator_->executionTimer().start();
            modelSimulator_->model().applyInitialSolution();
        }

        const EclipseState& eclState() const
        { return modelSimulator_->vanguard().eclState(); }

        EclipseState& eclState()
        { return modelSimulator_->vanguard().eclState(); }

        const Schedule& schedule() const
        { return modelSimulator_->vanguard().schedule(); }

        // Run the simulator.
        int runSimulator()
        {
            return runSimulatorInitOrRun_(&FlowMain::runSimulatorRunCallback_);
        }

        int runSimulatorInit()
        {
            return runSimulatorInitOrRun_(&FlowMain::runSimulatorInitCallback_);
        }

    private:
        // Callback that will be called from runSimulatorInitOrRun_().
        int runSimulatorRunCallback_()
        {
#ifdef RESERVOIR_COUPLING_ENABLED
            SimulatorReport report = simulator_->run(*simtimer_, this->argc_, this->argv_);
#else
            SimulatorReport report = simulator_->run(*simtimer_);
#endif
            runSimulatorAfterSim_(report);
            return report.success.exit_status;
        }

        // Callback that will be called from runSimulatorInitOrRun_().
        int runSimulatorInitCallback_()
        {
#ifdef RESERVOIR_COUPLING_ENABLED
            simulator_->init(*simtimer_, this->argc_, this->argv_);
#else
            simulator_->init(*simtimer_);
#endif
            return EXIT_SUCCESS;
        }

        // Output summary after simulation has completed
        void runSimulatorAfterSim_(SimulatorReport &report)
        {
            if (simulator_->model().hasNlddSolver()) {
                const auto& odir = eclState().getIOConfig().getOutputDir();
                // Write the number of nonlinear iterations per cell to a file in ResInsight compatible format
                simulator_->model().writeNonlinearIterationsPerCell(odir);
                // Write the NLDD statistics to the DBG file
                reportNlddStatistics(simulator_->model().domainAccumulatedReports(),
                                     simulator_->model().localAccumulatedReports(),
                                     this->output_cout_,
                                     FlowGenericVanguard::comm());
            }

            if (! this->output_cout_) {
                return;
            }

            const int threads
#if !defined(_OPENMP) || !_OPENMP
                = 1;
#else
                = omp_get_max_threads();
#endif

            printFlowTrailer(mpi_size_, threads, total_setup_time_, deck_read_time_, report);

            detail::handleExtraConvergenceOutput(report,
                                                 Parameters::Get<Parameters::OutputExtraConvergenceInfo>(),
                                                 R"(OutputExtraConvergenceInfo (--output-extra-convergence-info))",
                                                 eclState().getIOConfig().getOutputDir(),
                                                 eclState().getIOConfig().getBaseName());
        }

        // Run the simulator.
        int runSimulatorInitOrRun_(int (FlowMain::* initOrRunFunc)())
        {

            const auto& schedule = this->schedule();
            auto& ioConfig = eclState().getIOConfig();
            simtimer_ = std::make_unique<SimulatorTimer>();

            // initialize variables
            const auto& initConfig = eclState().getInitConfig();
            simtimer_->init(schedule, static_cast<std::size_t>(initConfig.getRestartStep()));

            if (this->output_cout_) {
                std::ostringstream oss;

                // This allows a user to catch typos and misunderstandings in the
                // use of simulator parameters.
                if (Parameters::printUnused(oss)) {
                    std::cout << "-----------------   Unrecognized parameters:   -----------------\n";
                    std::cout << oss.str();
                    std::cout << "----------------------------------------------------------------" << std::endl;
                }
            }

            if (!ioConfig.initOnly()) {
                if (this->output_cout_) {
                    std::string msg;
                    msg = "\n\n================ Starting main simulation loop ===============\n";
                    OpmLog::info(msg);
                }

                return (this->*initOrRunFunc)();
            }
            else {
                if (this->output_cout_) {
                    std::cout << "\n\n================ Simulation turned off ===============\n" << std::flush;
                }
                return EXIT_SUCCESS;
            }
        }

    protected:

        /// This is the main function of Flow.
        // Create simulator instance.
        // Writes to:
        //   simulator_
        void createSimulator()
        {
            // Create the simulator instance.
            simulator_ = std::make_unique<Simulator>(*modelSimulator_);
        }

        Grid& grid()
        { return modelSimulator_->vanguard().grid(); }

    private:
        std::unique_ptr<ModelSimulator> modelSimulator_;
        int  mpi_rank_ = 0;
        int  mpi_size_ = 1;
        std::any parallel_information_;
        std::unique_ptr<Simulator> simulator_;
        std::unique_ptr<SimulatorTimer> simtimer_;
        int argc_;
        char **argv_;
        bool output_cout_;
        bool output_files_;
        double total_setup_time_ = 0.0;
        double deck_read_time_ = 0.0;
    };

} // namespace Opm

#endif // OPM_FLOW_MAIN_HEADER_INCLUDED
