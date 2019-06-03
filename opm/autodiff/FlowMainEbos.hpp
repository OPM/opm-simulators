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
#ifndef OPM_FLOW_MAIN_EBOS_HEADER_INCLUDED
#define OPM_FLOW_MAIN_EBOS_HEADER_INCLUDED


#include <sys/utsname.h>

#include <opm/simulators/utils/ParallelFileMerger.hpp>

#include <opm/autodiff/BlackoilModelEbos.hpp>
#include <opm/autodiff/MissingFeatures.hpp>
#include <opm/simulators/utils/moduleVersion.hpp>
#include <opm/simulators/linalg/ExtractParallelGridInformationToISTL.hpp>
#include <opm/autodiff/SimulatorFullyImplicitBlackoilEbos.hpp>

#include <opm/core/props/satfunc/RelpermDiagnostics.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/OpmLog/EclipsePRTLog.hpp>
#include <opm/common/OpmLog/LogUtil.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/IOConfig/IOConfig.hpp>
#include <opm/parser/eclipse/EclipseState/InitConfig/InitConfig.hpp>
#include <opm/parser/eclipse/EclipseState/checkDeck.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

BEGIN_PROPERTIES

NEW_PROP_TAG(OutputMode);
NEW_PROP_TAG(EnableDryRun);
NEW_PROP_TAG(OutputInterval);
NEW_PROP_TAG(UseAmg);
NEW_PROP_TAG(EnableLoggingFalloutWarning);

SET_STRING_PROP(EclFlowProblem, OutputMode, "all");

// TODO: enumeration parameters. we use strings for now.
SET_STRING_PROP(EclFlowProblem, EnableDryRun, "auto");
// Do not merge parallel output files or warn about them
SET_BOOL_PROP(EclFlowProblem, EnableLoggingFalloutWarning, false);
SET_INT_PROP(EclFlowProblem, OutputInterval, 1);

END_PROPERTIES

namespace Opm
{
    // The FlowMain class is the ebos based black-oil simulator.
    template <class TypeTag>
    class FlowMainEbos
    {
        enum FileOutputMode
        {
            //! \brief No output to files.
            OUTPUT_NONE = 0,
            //! \brief Output only to log files, no eclipse output.
            OUTPUT_LOG_ONLY = 1,
            //! \brief Output to all files.
            OUTPUT_ALL = 3
        };

    public:
        typedef typename GET_PROP(TypeTag, MaterialLaw)::EclMaterialLawManager MaterialLawManager;
        typedef typename GET_PROP_TYPE(TypeTag, Simulator) EbosSimulator;
        typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
        typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
        typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

        typedef Opm::SimulatorFullyImplicitBlackoilEbos<TypeTag> Simulator;

        // Read the command line parameters. Throws an exception if something goes wrong.
        static int setupParameters_(int argc, char** argv)
        {
            // register the flow specific parameters
            EWOMS_REGISTER_PARAM(TypeTag, std::string, OutputMode,
                                 "Specify which messages are going to be printed. Valid values are: none, log, all (default)");
            EWOMS_REGISTER_PARAM(TypeTag, std::string, EnableDryRun,
                                 "Specify if the simulation ought to be actually run, or just pretended to be");
            EWOMS_REGISTER_PARAM(TypeTag, int, OutputInterval,
                                 "Specify the number of report steps between two consecutive writes of restart data");
            EWOMS_REGISTER_PARAM(TypeTag, bool, EnableLoggingFalloutWarning,
                                 "Developer option to see whether logging was on non-root processors. In that case it will be appended to the *.DBG or *.PRT files");

            Simulator::registerParameters();

            // register the parameters inherited from ebos
            Ewoms::registerAllParameters_<TypeTag>(/*finalizeRegistration=*/false);

            // hide the parameters unused by flow. TODO: this is a pain to maintain
            EWOMS_HIDE_PARAM(TypeTag, EnableGravity);
            EWOMS_HIDE_PARAM(TypeTag, EnableGridAdaptation);

            // this parameter is actually used in eWoms, but the flow well model
            // hard-codes the assumption that the intensive quantities cache is enabled,
            // so flow crashes. Let's hide the parameter for that reason.
            EWOMS_HIDE_PARAM(TypeTag, EnableIntensiveQuantityCache);

            // thermodynamic hints are not implemented/required by the eWoms blackoil
            // model
            EWOMS_HIDE_PARAM(TypeTag, EnableThermodynamicHints);

            // in flow only the deck file determines the end time of the simulation
            EWOMS_HIDE_PARAM(TypeTag, EndTime);

            // time stepping is not done by the eWoms code in flow
            EWOMS_HIDE_PARAM(TypeTag, InitialTimeStepSize);
            EWOMS_HIDE_PARAM(TypeTag, MaxTimeStepDivisions);
            EWOMS_HIDE_PARAM(TypeTag, MaxTimeStepSize);
            EWOMS_HIDE_PARAM(TypeTag, MinTimeStepSize);
            EWOMS_HIDE_PARAM(TypeTag, PredeterminedTimeStepsFile);

            EWOMS_HIDE_PARAM(TypeTag, EclMaxTimeStepSizeAfterWellEvent);
            EWOMS_HIDE_PARAM(TypeTag, EclRestartShrinkFactor);
            EWOMS_HIDE_PARAM(TypeTag, EclMaxFails);
            EWOMS_HIDE_PARAM(TypeTag, EclEnableTuning);

            // flow also does not use the eWoms Newton method
            EWOMS_HIDE_PARAM(TypeTag, NewtonMaxError);
            EWOMS_HIDE_PARAM(TypeTag, NewtonMaxIterations);
            EWOMS_HIDE_PARAM(TypeTag, NewtonTolerance);
            EWOMS_HIDE_PARAM(TypeTag, NewtonTargetIterations);
            EWOMS_HIDE_PARAM(TypeTag, NewtonVerbose);
            EWOMS_HIDE_PARAM(TypeTag, NewtonWriteConvergence);
            EWOMS_HIDE_PARAM(TypeTag, EclNewtonSumTolerance);
            EWOMS_HIDE_PARAM(TypeTag, EclNewtonSumToleranceExponent);
            EWOMS_HIDE_PARAM(TypeTag, EclNewtonStrictIterations);
            EWOMS_HIDE_PARAM(TypeTag, EclNewtonRelaxedVolumeFraction);
            EWOMS_HIDE_PARAM(TypeTag, EclNewtonRelaxedTolerance);

            // the default eWoms checkpoint/restart mechanism does not work with flow
            EWOMS_HIDE_PARAM(TypeTag, RestartTime);
            EWOMS_HIDE_PARAM(TypeTag, RestartWritingInterval);

            EWOMS_END_PARAM_REGISTRATION(TypeTag);

            // read in the command line parameters
            int status = Ewoms::setupParameters_<TypeTag>(argc, const_cast<const char**>(argv), /*doRegistration=*/false);
            if (status == 0) {
                // deal with --print-properties and --print-parameters

                bool doExit = false;

                int mpiRank = 0;
#if HAVE_MPI
                MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
#endif
                if (EWOMS_GET_PARAM(TypeTag, int, PrintProperties) == 1) {
                    doExit = true;
                    if (mpiRank == 0)
                        Ewoms::Properties::printValues<TypeTag>();
                }

                if (EWOMS_GET_PARAM(TypeTag, int, PrintParameters) == 1) {
                    doExit = true;
                    if (mpiRank == 0)
                        Ewoms::Parameters::printValues<TypeTag>();
                }

                if (doExit)
                    return -1;
            }

            return status;
        }

        static void printBanner()
        {
            const int lineLen = 70;
            const std::string version = moduleVersionName();
            const std::string banner = "This is flow "+version;
            const int bannerPreLen = (lineLen - 2 - banner.size())/2;
            const int bannerPostLen = bannerPreLen + (lineLen - 2 - banner.size())%2;
            std::cout << "**********************************************************************\n";
            std::cout << "*                                                                    *\n";
            std::cout << "*" << std::string(bannerPreLen, ' ') << banner << std::string(bannerPostLen, ' ') << "*\n";
            std::cout << "*                                                                    *\n";
            std::cout << "* Flow is a simulator for fully implicit three-phase black-oil flow, *\n";
            std::cout << "*             including solvent and polymer capabilities.            *\n";
            std::cout << "*          For more information, see https://opm-project.org         *\n";
            std::cout << "*                                                                    *\n";
            std::cout << "**********************************************************************\n\n";
        }

        /// This is the main function of Flow.  It runs a complete simulation with the
        /// given grid and simulator classes, based on the user-specified command-line
        /// input.
        int execute(int argc, char** argv)
        {
            try {
                // deal with some administrative boilerplate

                int status = setupParameters_(argc, argv);
                if (status)
                    return status;

                setupParallelism();
                setupOutput();
                setupEbosSimulator();
                setupLogging();
                int unknownKeyWords = printPRTHeader();
#if HAVE_MPI
                int globalUnknownKeyWords;
                MPI_Allreduce(&unknownKeyWords,  &globalUnknownKeyWords, 1, MPI_INT,  MPI_SUM, MPI_COMM_WORLD);
                unknownKeyWords = globalUnknownKeyWords;
#endif
                if ( unknownKeyWords )
                {
#if HAVE_MPI
                    MPI_Finalize();
#endif
                    exit(EXIT_FAILURE);
                }
                runDiagnostics();
                createSimulator();

                // do the actual work
                runSimulator();

                // clean up
                mergeParallelLogFiles();

                return EXIT_SUCCESS;
            }
            catch (const std::exception& e) {
                std::ostringstream message;
                message  << "Program threw an exception: " << e.what();

                if (output_cout_) {
                    // in some cases exceptions are thrown before the logging system is set
                    // up.
                    if (OpmLog::hasBackend("STREAMLOG")) {
                        OpmLog::error(message.str());
                    }
                    else {
                        std::cout << message.str() << "\n";
                    }
                }

                return EXIT_FAILURE;
            }
        }

    protected:
        void setupParallelism()
        {
            // determine the rank of the current process and the number of processes
            // involved in the simulation. MPI must have already been initialized
            // here. (yes, the name of this method is misleading.)
#if HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
            MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);
#else
            mpi_rank_ = 0;
            mpi_size_ = 1;
#endif

#if _OPENMP
            // if openMP is available, default to 2 threads per process.
            if (!getenv("OMP_NUM_THREADS"))
                omp_set_num_threads(std::min(2, omp_get_num_procs()));
#endif

            typedef typename GET_PROP_TYPE(TypeTag, ThreadManager) ThreadManager;
            ThreadManager::init();
        }

        // Extract the minimum priority and determines if log files ought to be created.
        // Writes to:
        //   output_to_files_
        //   output_
        void setupOutput()
        {
            const std::string outputModeString =
                EWOMS_GET_PARAM(TypeTag, std::string, OutputMode);
            static std::map<std::string, FileOutputMode> stringToOutputMode =
                { {"none", OUTPUT_NONE },
                  {"false", OUTPUT_LOG_ONLY },
                  {"log", OUTPUT_LOG_ONLY },
                  {"all" , OUTPUT_ALL },
                  {"true" , OUTPUT_ALL }};
            auto outputModeIt = stringToOutputMode.find(outputModeString);
            if (outputModeIt != stringToOutputMode.end()) {
                output_ = outputModeIt->second;
            }
            else {
                output_ = OUTPUT_ALL;
                std::cerr << "Value " << outputModeString <<
                    " is not a recognized output mode. Using \"all\" instead."
                          << std::endl;
            }

            output_cout_ = false;
            if (mpi_rank_ == 0) {
                output_cout_ = EWOMS_GET_PARAM(TypeTag, bool, EnableTerminalOutput);
                output_to_files_ = (output_ != OUTPUT_NONE);
            }
        }

        // Setup the OpmLog backends
        void setupLogging()
        {
            std::string deck_filename = EWOMS_GET_PARAM(TypeTag, std::string, EclDeckFileName);
            // create logFile
            using boost::filesystem::path;
            path fpath(deck_filename);
            std::string baseName;
            std::ostringstream debugFileStream;
            std::ostringstream logFileStream;

            // Strip extension "." or ".DATA"
            std::string extension = boost::to_upper_copy(fpath.extension().string());
            if ( extension == ".DATA" || extension == "." )
            {
                baseName = boost::to_upper_copy(fpath.stem().string());
            }
            else
            {
                baseName = boost::to_upper_copy(fpath.filename().string());
            }

            const std::string& output_dir = eclState().getIOConfig().getOutputDir();
            logFileStream << output_dir << "/" << baseName;
            debugFileStream << output_dir << "/" << baseName;

            if (mpi_rank_ != 0) {
                // Added rank to log file for non-zero ranks.
                // This prevents message loss.
                debugFileStream << "."<< mpi_rank_;
                // If the following file appears then there is a bug.
                logFileStream << "." << mpi_rank_;
            }
            logFileStream << ".PRT";
            debugFileStream << ".DBG";

            logFile_ = logFileStream.str();

            if (output_ > OUTPUT_NONE) {
                std::shared_ptr<EclipsePRTLog> prtLog = std::make_shared<EclipsePRTLog>(logFile_ , Log::NoDebugMessageTypes, false, output_cout_);
                OpmLog::addBackend( "ECLIPSEPRTLOG" , prtLog );
                prtLog->setMessageLimiter(std::make_shared<MessageLimiter>());
                prtLog->setMessageFormatter(std::make_shared<SimpleMessageFormatter>(false));
            }

            if (output_ >= OUTPUT_LOG_ONLY) {
                std::string debugFile = debugFileStream.str();
                std::shared_ptr<StreamLog> debugLog = std::make_shared<EclipsePRTLog>(debugFile, Log::DefaultMessageTypes, false, output_cout_);
                OpmLog::addBackend("DEBUGLOG",  debugLog);
            }

            std::shared_ptr<StreamLog> streamLog = std::make_shared<StreamLog>(std::cout, Log::StdoutMessageTypes);
            OpmLog::addBackend( "STREAMLOG", streamLog);
            const auto& msgLimits = schedule().getMessageLimits();
            const std::map<int64_t, int> limits = {{Log::MessageType::Note, msgLimits.getCommentPrintLimit(0)},
                                                   {Log::MessageType::Info, msgLimits.getMessagePrintLimit(0)},
                                                   {Log::MessageType::Warning, msgLimits.getWarningPrintLimit(0)},
                                                   {Log::MessageType::Error, msgLimits.getErrorPrintLimit(0)},
                                                   {Log::MessageType::Problem, msgLimits.getProblemPrintLimit(0)},
                                                   {Log::MessageType::Bug, msgLimits.getBugPrintLimit(0)}};
            streamLog->setMessageLimiter(std::make_shared<MessageLimiter>(10, limits));
            streamLog->setMessageFormatter(std::make_shared<SimpleMessageFormatter>(true));
        }

        // Print an ASCII-art header to the PRT and DEBUG files.
        // \return Whether unkown keywords were seen during parsing.
        bool printPRTHeader()
        {
          if (output_cout_) {
              const std::string version = moduleVersion();
              const double megabyte = 1024 * 1024;
              unsigned num_cpu = std::thread::hardware_concurrency();
              struct utsname arch;
              const char* user = getlogin();
              time_t now = std::time(0);
              struct tm  tstruct;
              char      tmstr[80];
              tstruct = *localtime(&now);
              strftime(tmstr, sizeof(tmstr), "%d-%m-%Y at %X", &tstruct);
              const double mem_size = getTotalSystemMemory() / megabyte;
              std::ostringstream ss;
              ss << "\n\n\n";
              ss << " ########  #          ######   #           #\n";
              ss << " #         #         #      #   #         # \n";
              ss << " #####     #         #      #    #   #   #  \n";
              ss << " #         #         #      #     # # # #   \n";
              ss << " #         #######    ######       #   #    \n\n";
              ss << "Flow is a simulator for fully implicit three-phase black-oil flow,";
              ss << " and is part of OPM.\nFor more information visit: https://opm-project.org \n\n";
              ss << "Flow Version  =  " + version + "\n";
              if (uname(&arch) == 0) {
                 ss << "System        =  " << arch.nodename << " (Number of logical cores: " << num_cpu;
                 ss << ", RAM: " << std::fixed << std::setprecision (2) << mem_size << " MB) \n";
                 ss << "Architecture  =  " << arch.sysname << " " << arch.machine << " (Release: " << arch.release;
                 ss << ", Version: " << arch.version << " )\n";
                 }
              if (user) {
                 ss << "User          =  " << user << std::endl;
                 }
              ss << "Simulation started on " << tmstr << " hrs\n";

              ss << "Parameters used by Flow:\n";
              Ewoms::Parameters::printValues<TypeTag>(ss);

              OpmLog::note(ss.str());
          }

              if ( mpi_rank_ == 0 )
              {
                  return Ewoms::Parameters::printUnused<TypeTag>(std::cerr);
              }
              else
              {
                  return false;
              }
        }

        void mergeParallelLogFiles()
        {
            // force closing of all log files.
            OpmLog::removeAllBackends();

            if (mpi_rank_ != 0 || mpi_size_ < 2 || !output_to_files_) {
                return;
            }

            namespace fs = boost::filesystem;
            const std::string& output_dir = eclState().getIOConfig().getOutputDir();
            fs::path output_path(output_dir);
            fs::path deck_filename(EWOMS_GET_PARAM(TypeTag, std::string, EclDeckFileName));
            std::string basename;
            // Strip extension "." and ".DATA"
            std::string extension = boost::to_upper_copy(deck_filename.extension().string());
            if ( extension == ".DATA" || extension == "." )
            {
                basename = boost::to_upper_copy(deck_filename.stem().string());
            }
            else
            {
                basename = boost::to_upper_copy(deck_filename.filename().string());
            }
            std::for_each(fs::directory_iterator(output_path),
                          fs::directory_iterator(),
                          detail::ParallelFileMerger(output_path, basename,
                                                     EWOMS_GET_PARAM(TypeTag, bool, EnableLoggingFalloutWarning)));
        }

        void setupEbosSimulator()
        {
            ebosSimulator_.reset(new EbosSimulator(/*verbose=*/false));
            ebosSimulator_->executionTimer().start();
            ebosSimulator_->model().applyInitialSolution();

            try {
                if (output_cout_) {
                    MissingFeatures::checkKeywords(deck());
                }

                // Possible to force initialization only behavior (NOSIM).
                const std::string& dryRunString = EWOMS_GET_PARAM(TypeTag, std::string, EnableDryRun);
                if (dryRunString != "" && dryRunString != "auto") {
                    bool yesno;
                    if (dryRunString == "true"
                        || dryRunString == "t"
                        || dryRunString == "1")
                        yesno = true;
                    else if (dryRunString == "false"
                             || dryRunString == "f"
                             || dryRunString == "0")
                        yesno = false;
                    else
                        throw std::invalid_argument("Invalid value for parameter EnableDryRun: '"
                                                    +dryRunString+"'");
                    auto& ioConfig = eclState().getIOConfig();
                    ioConfig.overrideNOSIM(yesno);
                }
            }
            catch (const std::invalid_argument& e) {
                std::cerr << "Failed to create valid EclipseState object. See logfile: " << logFile_ << std::endl;
                std::cerr << "Exception caught: " << e.what() << std::endl;
                throw;
            }
        }

        const Deck& deck() const
        { return ebosSimulator_->vanguard().deck(); }

        Deck& deck()
        { return ebosSimulator_->vanguard().deck(); }

        const EclipseState& eclState() const
        { return ebosSimulator_->vanguard().eclState(); }

        EclipseState& eclState()
        { return ebosSimulator_->vanguard().eclState(); }

        const Schedule& schedule() const
        { return ebosSimulator_->vanguard().schedule(); }


        // Run diagnostics.
        // Writes to:
        //   OpmLog singleton.
        void runDiagnostics()
        {
            if (!output_cout_) {
                return;
            }

            // Run relperm diagnostics
            RelpermDiagnostics diagnostic;
            diagnostic.diagnosis(eclState(), deck(), this->grid());
        }

        // Run the simulator.
        void runSimulator()
        {
            const auto& schedule = this->schedule();
            const auto& timeMap = schedule.getTimeMap();
            auto& ioConfig = eclState().getIOConfig();
            SimulatorTimer simtimer;

            // initialize variables
            const auto& initConfig = eclState().getInitConfig();
            simtimer.init(timeMap, (size_t)initConfig.getRestartStep());

            if (output_cout_) {
                std::ostringstream oss;

                // This allows a user to catch typos and misunderstandings in the
                // use of simulator parameters.
                if (Ewoms::Parameters::printUnused<TypeTag>(oss)) {
                    std::cout << "-----------------   Unrecognized parameters:   -----------------\n";
                    std::cout << oss.str();
                    std::cout << "----------------------------------------------------------------" << std::endl;
                }
            }

            if (!ioConfig.initOnly()) {
                if (output_cout_) {
                    std::string msg;
                    msg = "\n\n================ Starting main simulation loop ===============\n";
                    OpmLog::info(msg);
                }

                SimulatorReport successReport = simulator_->run(simtimer);
                SimulatorReport failureReport = simulator_->failureReport();

                if (output_cout_) {
                    std::ostringstream ss;
                    ss << "\n\n================    End of simulation     ===============\n\n";
                    successReport.reportFullyImplicit(ss, &failureReport);
                    OpmLog::info(ss.str());
                }

            } else {
                if (output_cout_) {
                    std::cout << "\n\n================ Simulation turned off ===============\n" << std::flush;
                }

            }
        }

        /// This is the main function of Flow.
        // Create simulator instance.
        // Writes to:
        //   simulator_
        void createSimulator()
        {
            // Create the simulator instance.
            simulator_.reset(new Simulator(*ebosSimulator_));
        }

        unsigned long long getTotalSystemMemory()
        {
            long pages = sysconf(_SC_PHYS_PAGES);
            long page_size = sysconf(_SC_PAGE_SIZE);
            return pages * page_size;
        }


        Grid& grid()
        { return ebosSimulator_->vanguard().grid(); }

    private:
        std::unique_ptr<EbosSimulator> ebosSimulator_;
        int  mpi_rank_ = 0;
        int  mpi_size_ = 1;
        bool output_cout_ = false;
        FileOutputMode output_ = OUTPUT_ALL;
        bool output_to_files_ = false;
        boost::any parallel_information_;
        std::unique_ptr<Simulator> simulator_;
        std::string logFile_;
    };
} // namespace Opm

#endif // OPM_FLOW_MAIN_EBOS_HEADER_INCLUDED
