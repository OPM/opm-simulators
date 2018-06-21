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

#include <opm/simulators/ParallelFileMerger.hpp>

#include <opm/autodiff/BlackoilModelEbos.hpp>
#include <opm/autodiff/NewtonIterationBlackoilInterface.hpp>
#include <opm/autodiff/MissingFeatures.hpp>
#include <opm/autodiff/moduleVersion.hpp>
#include <opm/autodiff/ExtractParallelGridInformationToISTL.hpp>
#include <opm/autodiff/RedistributeDataHandles.hpp>
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

BEGIN_PROPERTIES;

NEW_PROP_TAG(FlowOutputMode);
NEW_PROP_TAG(FlowEnableDryRun);
NEW_PROP_TAG(FlowOutputInterval);
NEW_PROP_TAG(FlowUseAmg);

SET_STRING_PROP(EclFlowProblem, FlowOutputMode, "all");

// TODO: enumeration parameters. we use strings for now.
SET_STRING_PROP(EclFlowProblem, FlowEnableDryRun, "auto");

SET_INT_PROP(EclFlowProblem, FlowOutputInterval, 1);

END_PROPERTIES;

namespace Opm
{
    // The FlowMain class is the ebos based black-oil simulator.
    template <class TypeTag>
    class FlowMainEbos
    {
        enum FileOutputMode
        {
            //! \brief No output to files.
            OutputModeNone = 0,
            //! \brief Output only to log files, no eclipse output.
            OutputModeLogOnly = 1,
            //! \brief Output to all files.
            OutputModeAll = 3
        };

    public:
        typedef typename GET_PROP(TypeTag, MaterialLaw)::EclMaterialLawManager MaterialLawManager;
        typedef typename GET_PROP_TYPE(TypeTag, Simulator) EbosSimulator;
        typedef typename GET_PROP_TYPE(TypeTag, ThreadManager) EbosThreadManager;
        typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
        typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
        typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

        typedef Opm::SimulatorFullyImplicitBlackoilEbos<TypeTag> Simulator;
        typedef typename Simulator::ReservoirState ReservoirState;

        // Read the command line parameters. Throws an exception if something goes wrong.
        static int setupParameters_(int argc, char** argv)
        {
            // register the flow specific parameters
            EWOMS_REGISTER_PARAM(TypeTag, std::string, FlowOutputMode,
                                 "Specify which messages are going to be printed. Valid values are: none, log, all (default)");
            EWOMS_REGISTER_PARAM(TypeTag, std::string, FlowEnableDryRun,
                                 "Specify if the simulation ought to be actually run, or just pretended to be");
            EWOMS_REGISTER_PARAM(TypeTag, int, FlowOutputInterval,
                                 "Specify the number of report steps between two consecutive writes of restart data");
            Simulator::registerParameters();

            typedef typename BlackoilModelEbos<TypeTag>::ISTLSolverType ISTLSolverType;
            ISTLSolverType::registerParameters();

            // register the parameters inherited from ebos
            Ewoms::registerAllParameters_<TypeTag>();

            // read in the command line parameters
            return Ewoms::setupParameters_<TypeTag>(argc, const_cast<const char**>(argv), /*doRegistration=*/false);
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

                setupParallelism_();
                setupOutput_();
                printStartupMessage_();
                setupEbosSimulator_();
                setupLogging_();
                printPRTHeader_();
                runDiagnostics_();
                setupLinearSolver_();
                createSimulator_();

                // do the actual work
                runSimulator_();

                // clean up
                mergeParallelLogFiles_();

                return EXIT_SUCCESS;
            }
            catch (const std::exception& e) {
                std::ostringstream message;
                message  << "Program threw an exception: " << e.what();

                if (outputCout_) {
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
        void setupParallelism_()
        {
            // determine the rank of the current process and the number of processes
            // involved in the simulation. MPI must have already been initialized
            // here. (yes, the name of this method is misleading.)
#if HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank_);
            MPI_Comm_size(MPI_COMM_WORLD, &mpiSize_);
#else
            mpiRank_ = 0;
            mpiSize_ = 1;
#endif
        }

        // Print startup message if on output rank.
        void printStartupMessage_()
        {

            if (outputCout_) {
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
                std::cout << "*          For more information, see https://opm-project.org          *\n";
                std::cout << "*                                                                    *\n";
                std::cout << "**********************************************************************\n\n";
            }
        }

        // Extract the minimum priority and determines if log files ought to be created.
        // Writes to:
        //   outputToFiles_
        //   outputMode_
        void setupOutput_()
        {
            const std::string outputModeString =
                EWOMS_GET_PARAM(TypeTag, std::string, FlowOutputMode);
            static std::map<std::string, FileOutputMode> stringToOutputMode =
                { {"none", OutputModeNone },
                  {"false", OutputModeLogOnly },
                  {"log", OutputModeLogOnly },
                  {"all" , OutputModeAll },
                  {"true" , OutputModeAll }};
            auto outputModeIt = stringToOutputMode.find(outputModeString);
            if (outputModeIt != stringToOutputMode.end()) {
                outputMode_ = outputModeIt->second;
            }
            else {
                outputMode_ = OutputModeAll;
                std::cerr << "Value " << outputModeString <<
                    " is not a recognized output mode. Using \"all\" instead."
                          << std::endl;
            }

            outputCout_ = false;
            if (mpiRank_ == 0) {
                outputCout_ = EWOMS_GET_PARAM(TypeTag, bool, FlowEnableTerminalOutput);
                outputToFiles_ = (outputMode_ != OutputModeNone);
            }
        }

        // Setup the OpmLog backends
        void setupLogging_()
        {
            std::string deckFilename = EWOMS_GET_PARAM(TypeTag, std::string, EclDeckFileName);
            // create logFile
            using boost::filesystem::path;
            path fpath(deckFilename);
            std::string baseName;
            std::ostringstream debugFileStream;
            std::ostringstream logFileStream;

            if (boost::to_upper_copy(path(fpath.extension()).string()) == ".DATA") {
                baseName = path(fpath.stem()).string();
            } else {
                baseName = path(fpath.filename()).string();
            }

            const std::string& outputDir = eclState_().getIOConfig().getOutputDir();
            logFileStream << outputDir << "/" << baseName;
            debugFileStream << outputDir << "/" << baseName;

            if (mpiRank_ != 0) {
                // Added rank to log file for non-zero ranks.
                // This prevents message loss.
                debugFileStream << "."<< mpiRank_;
                // If the following file appears then there is a bug.
                logFileStream << "." << mpiRank_;
            }
            logFileStream << ".PRT";
            debugFileStream << ".DBG";

            logFile_ = logFileStream.str();

            if (outputMode_ > OutputModeNone) {
                std::shared_ptr<EclipsePRTLog> prtLog = std::make_shared<EclipsePRTLog>(logFile_ , Log::NoDebugMessageTypes, false, outputCout_);
                OpmLog::addBackend( "ECLIPSEPRTLOG" , prtLog );
                prtLog->setMessageLimiter(std::make_shared<MessageLimiter>());
                prtLog->setMessageFormatter(std::make_shared<SimpleMessageFormatter>(false));
            }

            if (outputMode_ >= OutputModeLogOnly) {
                std::string debugFile = debugFileStream.str();
                std::shared_ptr<StreamLog> debugLog = std::make_shared<EclipsePRTLog>(debugFile, Log::DefaultMessageTypes, false, outputCout_);
                OpmLog::addBackend("DEBUGLOG",  debugLog);
            }

            std::shared_ptr<StreamLog> streamLog = std::make_shared<StreamLog>(std::cout, Log::StdoutMessageTypes);
            OpmLog::addBackend( "STREAMLOG", streamLog);
            const auto& msgLimits = schedule_().getMessageLimits();
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
        void printPRTHeader_()
        {
          if (outputCout_) {
              const std::string version = moduleVersionName();
              const double megabyte = 1024 * 1024;
              unsigned numCores = std::thread::hardware_concurrency();
              struct utsname arch;
              const char* user = getlogin();
              time_t now = std::time(0);
              struct tm  tstruct;
              char      tmstr[80];
              tstruct = *localtime(&now);
              strftime(tmstr, sizeof(tmstr), "%d-%m-%Y at %X", &tstruct);
              const double memSize = getTotalSystemMemory_() / megabyte;
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
                 ss << "System        =  " << arch.nodename << " (Number of logical cores: " << numCores;
                 ss << ", RAM: " << std::fixed << std::setprecision (2) << memSize << " MB) \n";
                 ss << "Architecture  =  " << arch.sysname << " " << arch.machine << " (Release: " << arch.release;
                 ss << ", Version: " << arch.version << " )\n";
                 }
              if (user) {
                 ss << "User          =  " << user << std::endl;
                 }
              ss << "Simulation started on " << tmstr << " hrs\n";
              OpmLog::note(ss.str());
            }
        }

        void mergeParallelLogFiles_()
        {
            // force closing of all log files.
            OpmLog::removeAllBackends();

            if (mpiRank_ != 0 || mpiSize_ < 2 || !outputToFiles_) {
                return;
            }

            namespace fs = boost::filesystem;
            fs::path outputPath(".");
            const std::string& outputDir = eclState_().getIOConfig().getOutputDir();

            fs::path deckFileName(EWOMS_GET_PARAM(TypeTag, std::string, EclDeckFileName));
            std::for_each(fs::directory_iterator(outputPath),
                          fs::directory_iterator(),
                          detail::ParallelFileMerger(outputPath, deckFileName.stem().string()));
        }

        void setupEbosSimulator_()
        {
            ebosSimulator_.reset(new EbosSimulator(/*verbose=*/false));
            ebosSimulator_->model().applyInitialSolution();

            try {
                if (outputCout_) {
                    MissingFeatures::checkKeywords(deck_());
                }

                // Possible to force initialization only behavior (NOSIM).
                const std::string& dryRunString = EWOMS_GET_PARAM(TypeTag, std::string, FlowEnableDryRun);
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
                        throw std::invalid_argument("Invalid value for parameter FlowEnableDryRun: '"
                                                    +dryRunString+"'");
                    auto& ioConfig = eclState_().getIOConfig();
                    ioConfig.overrideNOSIM(yesno);
                }
            }
            catch (const std::invalid_argument& e) {
                std::cerr << "Failed to create valid EclipseState object. See logfile: " << logFile_ << std::endl;
                std::cerr << "Exception caught: " << e.what() << std::endl;
                throw;
            }
        }

        const Deck& deck_() const
        { return ebosSimulator_->vanguard().deck(); }

        Deck& deck_()
        { return ebosSimulator_->vanguard().deck(); }

        const EclipseState& eclState_() const
        { return ebosSimulator_->vanguard().eclState(); }

        EclipseState& eclState_()
        { return ebosSimulator_->vanguard().eclState(); }

        const Schedule& schedule_() const
        { return ebosSimulator_->vanguard().schedule(); }


        // Run diagnostics.
        // Writes to:
        //   OpmLog singleton.
        void runDiagnostics_()
        {
            if (!outputCout_) {
                return;
            }

            // Run relperm diagnostics
            RelpermDiagnostics diagnostic;
            diagnostic.diagnosis(eclState_(), deck_(), grid_());
        }

        // Run the simulator.
        void runSimulator_()
        {
            const auto& schedule = schedule_();
            const auto& timeMap = schedule.getTimeMap();
            auto& ioConfig = eclState_().getIOConfig();
            SimulatorTimer simtimer;

            // initialize variables
            const auto& initConfig = eclState_().getInitConfig();
            simtimer.init(timeMap, (size_t)initConfig.getRestartStep());

            if (outputCout_) {
                // This allows a user to catch typos and misunderstandings in the
                // use of simulator parameters.
                std::cout << "-----------------   Unrecognized parameters:   -----------------\n";
                Ewoms::Parameters::printUnused<TypeTag>();
                std::cout << "----------------------------------------------------------------" << std::endl;
            }

            if (!ioConfig.initOnly()) {
                if (outputCout_) {
                    std::string msg;
                    msg = "\n\n================ Starting main simulation loop ===============\n";
                    OpmLog::info(msg);
                }

                SimulatorReport successReport = simulator_->run(simtimer);
                SimulatorReport failureReport = simulator_->failureReport();

                if (outputCout_) {
                    std::ostringstream ss;
                    ss << "\n\n================    End of simulation     ===============\n\n";
                    successReport.reportFullyImplicit(ss, &failureReport);
                    OpmLog::info(ss.str());
                }

            } else {
                if (outputCout_) {
                    std::cout << "\n\n================ Simulation turned off ===============\n" << std::flush;
                }

            }
        }

        // Setup linear solver.
        // Writes to:
        //   linearSolver_
        void setupLinearSolver_()
        {
            typedef typename BlackoilModelEbos<TypeTag>::ISTLSolverType ISTLSolverType;

            extractParallelGridInformationToISTL(grid_(), parallelInformation_);
            linearSolver_.reset(new ISTLSolverType(parallelInformation_));

            // Deactivate selection of CPR via eclipse keyword
            // as this preconditioner is still considered experimental
            // and fails miserably for some models.
            if (outputCout_) {
                std::ostringstream message;
                message << "Ignoring request for CPRPreconditioner "
                        << "via Eclipse keyword as it is considered "
                        <<" experimental. To activate use "
                        <<"\"--flow-use-cpr=true\" command "
                        <<"line parameter.";
                OpmLog::info(message.str());
            }
        }

        /// This is the main function of Flow.
        // Create simulator instance.
        // Writes to:
        //   simulator_
        void createSimulator_()
        {
            // Create the simulator instance.
            simulator_.reset(new Simulator(*ebosSimulator_, *linearSolver_));
        }

        unsigned long long getTotalSystemMemory_()
        {
            long pages = sysconf(_SC_PHYS_PAGES);
            long pageSize = sysconf(_SC_PAGE_SIZE);
            return pages * pageSize;
        }


        Grid& grid_()
        { return ebosSimulator_->vanguard().grid(); }

    private:
        std::unique_ptr<EbosSimulator> ebosSimulator_;
        int  mpiRank_ = 0;
        int  mpiSize_ = 1;
        bool outputCout_ = false;
        FileOutputMode outputMode_;
        bool outputToFiles_ = false;
        boost::any parallelInformation_;
        std::unique_ptr<NewtonIterationBlackoilInterface> linearSolver_;
        std::unique_ptr<Simulator> simulator_;
        std::string logFile_;
    };
} // namespace Opm

#endif // OPM_FLOW_MAIN_EBOS_HEADER_INCLUDED
