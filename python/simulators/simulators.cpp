#include "config.h"

#include <flow/flow_ebos_blackoil.hpp>

#include <opm/simulators/flow/FlowMainEbos.hpp>
#include <opm/simulators/utils/moduleVersion.hpp>
#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/material/common/ResetLocale.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/OpmLog/EclipsePRTLog.hpp>
#include <opm/common/OpmLog/LogUtil.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>

#include <opm/parser/eclipse/EclipseState/Schedule/ArrayDimChecker.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
#include <string>
#include <iostream>

namespace py = pybind11;

BEGIN_PROPERTIES

// this is a dummy type tag that is used to setup the parameters before the actual
// simulator.
NEW_TYPE_TAG(FlowEarlyBird, INHERITS_FROM(EclFlowProblem));

END_PROPERTIES

namespace detail
{
    boost::filesystem::path simulationCaseName( const std::string& casename ) {
        namespace fs = boost::filesystem;

        const auto exists = []( const fs::path& f ) -> bool {
            if( !fs::exists( f ) ) return false;

            if( fs::is_regular_file( f ) ) return true;

            return fs::is_symlink( f )
            && fs::is_regular_file( fs::read_symlink( f ) );
        };

        auto simcase = fs::path( casename );

        if( exists( simcase ) ) {
            return simcase;
        }

        for( const auto& ext : { std::string("data"), std::string("DATA") } ) {
            if( exists( simcase.replace_extension( ext ) ) ) {
                return simcase;
            }
        }

        throw std::invalid_argument( "Cannot find input case " + casename );
    }


    // This function is an extreme special case, if the program has been invoked
    // *exactly* as:
    //
    //    flow   --version
    //
    // the call is intercepted by this function which will print "flow $version"
    // on stdout and exit(0).
    void handleVersionCmdLine(int argc, char** argv) {
        for ( int i = 1; i < argc; ++i )
        {
            if (std::strcmp(argv[i], "--version") == 0) {
                std::cout << "flow " << Opm::moduleVersionName() << std::endl;
                std::exit(EXIT_SUCCESS);
            }
        }
    }
}

enum class FileOutputMode {
    //! \brief No output to files.
    OUTPUT_NONE = 0,
    //! \brief Output only to log files, no eclipse output.
    OUTPUT_LOG_ONLY = 1,
    //! \brief Output to all files.
    OUTPUT_ALL = 3
};

void ensureOutputDirExists(const std::string& cmdline_output_dir)
{
    if (!boost::filesystem::is_directory(cmdline_output_dir)) {
        try {
            boost::filesystem::create_directories(cmdline_output_dir);
        }
        catch (...) {
            throw std::runtime_error("Creation of output directory '" + cmdline_output_dir + "' failed\n");
        }
    }
}

// Setup the OpmLog backends
FileOutputMode setupLogging(int mpi_rank_, const std::string& deck_filename, const std::string& cmdline_output_dir, const std::string& cmdline_output, bool output_cout_, const std::string& stdout_log_id) {

    if (!cmdline_output_dir.empty()) {
        ensureOutputDirExists(cmdline_output_dir);
    }

    // create logFile
    using boost::filesystem::path;
    path fpath(deck_filename);
    std::string baseName;
    std::ostringstream debugFileStream;
    std::ostringstream logFileStream;

    // Strip extension "." or ".DATA"
    std::string extension = boost::to_upper_copy(fpath.extension().string());
    if (extension == ".DATA" || extension == ".") {
        baseName = boost::to_upper_copy(fpath.stem().string());
    } else {
        baseName = boost::to_upper_copy(fpath.filename().string());
    }

    std::string output_dir = cmdline_output_dir;
    if (output_dir.empty()) {
        output_dir = absolute(path(baseName).parent_path()).string();
    }

    logFileStream << output_dir << "/" << baseName;
    debugFileStream << output_dir << "/" << baseName;

    if (mpi_rank_ != 0) {
        // Added rank to log file for non-zero ranks.
        // This prevents message loss.
        debugFileStream << "." << mpi_rank_;
        // If the following file appears then there is a bug.
        logFileStream << "." << mpi_rank_;
    }
    logFileStream << ".PRT";
    debugFileStream << ".DBG";

    FileOutputMode output;
    {
        static std::map<std::string, FileOutputMode> stringToOutputMode =
            { {"none", FileOutputMode::OUTPUT_NONE },
              {"false", FileOutputMode::OUTPUT_LOG_ONLY },
              {"log", FileOutputMode::OUTPUT_LOG_ONLY },
              {"all" , FileOutputMode::OUTPUT_ALL },
              {"true" , FileOutputMode::OUTPUT_ALL }};
        auto outputModeIt = stringToOutputMode.find(cmdline_output);
        if (outputModeIt != stringToOutputMode.end()) {
            output = outputModeIt->second;
        }
        else {
            output = FileOutputMode::OUTPUT_ALL;
            std::cerr << "Value " << cmdline_output <<
                " is not a recognized output mode. Using \"all\" instead."
                      << std::endl;
        }
    }

    if (output > FileOutputMode::OUTPUT_NONE) {
        std::shared_ptr<Opm::EclipsePRTLog> prtLog = std::make_shared<Opm::EclipsePRTLog>(logFileStream.str(), Opm::Log::NoDebugMessageTypes, false, output_cout_);
        Opm::OpmLog::addBackend("ECLIPSEPRTLOG", prtLog);
        prtLog->setMessageLimiter(std::make_shared<Opm::MessageLimiter>());
        prtLog->setMessageFormatter(std::make_shared<Opm::SimpleMessageFormatter>(false));
    }

    if (output >= FileOutputMode::OUTPUT_LOG_ONLY) {
        std::string debugFile = debugFileStream.str();
        std::shared_ptr<Opm::StreamLog> debugLog = std::make_shared<Opm::EclipsePRTLog>(debugFileStream.str(), Opm::Log::DefaultMessageTypes, false, output_cout_);
        Opm::OpmLog::addBackend("DEBUGLOG", debugLog);
    }

    if (mpi_rank_ == 0) {
        std::shared_ptr<Opm::StreamLog> streamLog = std::make_shared<Opm::StreamLog>(std::cout, Opm::Log::StdoutMessageTypes);
        Opm::OpmLog::addBackend(stdout_log_id, streamLog);
        streamLog->setMessageFormatter(std::make_shared<Opm::SimpleMessageFormatter>(true));
    }
    return output;
}

void setupMessageLimiter(const Opm::MessageLimits msgLimits,  const std::string& stdout_log_id) {
    std::shared_ptr<Opm::StreamLog> stream_log = Opm::OpmLog::getBackend<Opm::StreamLog>(stdout_log_id);

    const std::map<int64_t, int> limits = {{Opm::Log::MessageType::Note,
                                            msgLimits.getCommentPrintLimit(0)},
                                           {Opm::Log::MessageType::Info,
                                            msgLimits.getMessagePrintLimit(0)},
                                           {Opm::Log::MessageType::Warning,
                                            msgLimits.getWarningPrintLimit(0)},
                                           {Opm::Log::MessageType::Error,
                                            msgLimits.getErrorPrintLimit(0)},
                                           {Opm::Log::MessageType::Problem,
                                            msgLimits.getProblemPrintLimit(0)},
                                           {Opm::Log::MessageType::Bug,
                                            msgLimits.getBugPrintLimit(0)}};
    stream_log->setMessageLimiter(std::make_shared<Opm::MessageLimiter>(10, limits));
}

class BlackOilSimulator
{
public:
    BlackOilSimulator( const std::string& filename )
    {
        setupCmdLineArgs();
        std::strcpy(argv_[1], filename.c_str());
    }

    BlackOilSimulator( const Opm::Deck& deck, 
                       const Opm::EclipseState& eclipseState, 
                       const Opm::Schedule& schedule,
                       const Opm::SummaryConfig& summaryConfig )
    {
        setupCmdLineArgs();
        setDeck(deck);
        setEclipseState(eclipseState);
        setSchedule(schedule);
        setSummaryConfig(summaryConfig);
    }

    ~BlackOilSimulator() 
    {
        delete[] argv_[0];
        delete[] argv_[1];
        delete[] argv_;
    }

    void setupCmdLineArgs()
    {
        argc_ = 2;
        argv_ = new char*[2];
        argv_[0] = new char[200];
        char argv0[] = "flow";
        std::strcpy(argv_[0], argv0);
        argv_[1] = new char[200];
    }

    void setDeck( const Opm::Deck& deck )
    {
        deck_ = std::make_shared< Opm::Deck >(deck);
    }

    void setEclipseState( const Opm::EclipseState& eclipseState )
    {
        eclipseState_ = std::make_shared< Opm::EclipseState >(eclipseState);
    }

    void setSchedule( const Opm::Schedule& schedule )
    {
        schedule_ = std::make_shared< Opm::Schedule >(schedule);
    }

    void setSummaryConfig( const Opm::SummaryConfig& summaryConfig )
    {
        summaryConfig_ = std::make_shared< Opm::SummaryConfig >(summaryConfig);
    }

    int run()
    {
        if (prepareRun_()) {
            return Opm::flowEbosBlackoilMain(argc_, argv_, outputCout_, outputFiles_);
        }
        else {
            return EXIT_FAILURE;
        }
    }

    int step_init()
    {
        if (hasRunInit_) {
            // Running step_init() multiple times is not implemented yet,
            // currently we just do nothing and return
            return EXIT_SUCCESS;
        }
        if (prepareRun_()) {
            mainfunc_ = Opm::flowEbosBlackoilMainInit(
                argc_, argv_, outputCout_, outputFiles_);
            int result = mainfunc_->executeInitStep(
                argc_, argv_, outputCout_, outputFiles_);
            hasRunInit_ = true;
            return result;
        }
        else {
            return EXIT_FAILURE;
        }
    }
private:

    bool prepareRun_()
    {
        int argc = argc_;
        char **argv = argv_;
        Dune::Timer externalSetupTimer;
        externalSetupTimer.start();

        detail::handleVersionCmdLine(argc, argv);
        // MPI setup.
#if HAVE_DUNE_FEM
        Dune::Fem::MPIManager::initialize(argc, argv);
        int mpiRank = Dune::Fem::MPIManager::rank();
#else
        // the design of the plain dune MPIHelper class is quite flawed: there is no way to
        // get the instance without having the argc and argv parameters available and it is
        // not possible to determine the MPI rank and size without an instance. (IOW: the
        // rank() and size() methods are supposed to be static.)
        const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);
        int mpiRank = mpiHelper.rank();
#endif

        // we always want to use the default locale, and thus spare us the trouble
        // with incorrect locale settings.
        Opm::resetLocale();

        // this is a work-around for a catch 22: we do not know what code path to use without
        // parsing the deck, but we don't know the deck without having access to the
        // parameters and this requires to know the type tag to be used. To solve this, we
        // use a type tag just for parsing the parameters before we instantiate the actual
        // simulator object. (Which parses the parameters again, but since this is done in an
        // identical manner it does not matter.)
        typedef TTAG(FlowEarlyBird) PreTypeTag;
        typedef GET_PROP_TYPE(PreTypeTag, Problem) PreProblem;

        PreProblem::setBriefDescription("Flow, an advanced reservoir simulator for ECL-decks provided by the Open Porous Media project.");

        int status = Opm::FlowMainEbos<PreTypeTag>::setupParameters_(argc, argv);
        if (status != 0) {
            // if setupParameters_ returns a value smaller than 0, there was no error, but
            // the program should abort. This is the case e.g. for the --help and the
            // --print-properties parameters.
#if HAVE_MPI
            MPI_Finalize();
#endif
            return (status >= 0)?status:0;
        }

        FileOutputMode outputMode = FileOutputMode::OUTPUT_NONE;
        outputCout_ = false;
        if (mpiRank == 0)
            outputCout_ = EWOMS_GET_PARAM(PreTypeTag, bool, EnableTerminalOutput);

        std::string deckFilename;
        std::string outputDir;
        if ( eclipseState_ )
        {
            deckFilename = eclipseState_->getIOConfig().fullBasePath();
            outputDir = eclipseState_->getIOConfig().getOutputDir();
        }
        else
        {
            deckFilename = EWOMS_GET_PARAM(PreTypeTag, std::string, EclDeckFileName);
            typedef typename GET_PROP_TYPE(PreTypeTag, Vanguard) PreVanguard;
            try {
                deckFilename = PreVanguard::canonicalDeckPath(deckFilename).string();
            }
            catch (const std::exception& e) {
                if ( mpiRank == 0 )
                    std::cerr << "Exception received: " << e.what() << ". Try '--help' for a usage description.\n";
#if HAVE_MPI
                MPI_Finalize();
#endif
                return 1;
            }
        }

        if (outputCout_) {
            Opm::FlowMainEbos<PreTypeTag>::printBanner();
        }

        // Create Deck, EclipseState, Schedule and SummaryConfig if they don't exist
        if ( deck_ && eclipseState_ && schedule_ && summaryConfig_ )
        {
            outputMode = setupLogging(mpiRank,
                                    deckFilename,
                                    outputDir,
                                    EWOMS_GET_PARAM(PreTypeTag, std::string, OutputMode),
                                    outputCout_, "STDOUT_LOGGER");

            if (mpiRank == 0) {
                setupMessageLimiter(schedule_->getMessageLimits(), "STDOUT_LOGGER");
            }
        }
        else
        {
            Opm::Parser parser;
            Opm::ParseContext parseContext({{Opm::ParseContext::PARSE_RANDOM_SLASH, Opm::InputError::IGNORE},
                                            {Opm::ParseContext::PARSE_MISSING_DIMS_KEYWORD, Opm::InputError::WARN},
                                            {Opm::ParseContext::SUMMARY_UNKNOWN_WELL, Opm::InputError::WARN},
                                            {Opm::ParseContext::SUMMARY_UNKNOWN_GROUP, Opm::InputError::WARN}});
            Opm::ErrorGuard errorGuard;
            outputMode = setupLogging(mpiRank,
                                      deckFilename,
                                      EWOMS_GET_PARAM(PreTypeTag, std::string, OutputDir),
                                      EWOMS_GET_PARAM(PreTypeTag, std::string, OutputMode),
                                      outputCout_, "STDOUT_LOGGER");

            if (EWOMS_GET_PARAM(PreTypeTag, bool, EclStrictParsing))
                parseContext.update( Opm::InputError::DELAYED_EXIT1);

            Opm::FlowMainEbos<PreTypeTag>::printPRTHeader(outputCout_);

            deck_.reset( new Opm::Deck( parser.parseFile(deckFilename , parseContext, errorGuard)));
            Opm::MissingFeatures::checkKeywords(*deck_, parseContext, errorGuard);
            if ( outputCout_ )
                Opm::checkDeck(*deck_, parser, parseContext, errorGuard);

            eclipseState_.reset( new Opm::EclipseState(*deck_, parseContext, errorGuard ));
            schedule_.reset(new Opm::Schedule(*deck_, *eclipseState_, parseContext, errorGuard));
            summaryConfig_.reset( new Opm::SummaryConfig(*deck_, *schedule_, eclipseState_->getTableManager(), parseContext, errorGuard));
            if (mpiRank == 0) {
                setupMessageLimiter(schedule_->getMessageLimits(), "STDOUT_LOGGER");
            }

            Opm::checkConsistentArrayDimensions(*eclipseState_, *schedule_, parseContext, errorGuard);

            if (errorGuard) {
                errorGuard.dump();
                errorGuard.clear();

                throw std::runtime_error("Unrecoverable errors were encountered while loading input.");
            }
        }

        const auto& phases = Opm::Runspec(*deck_).phases();
        outputFiles_ = (outputMode != FileOutputMode::OUTPUT_NONE);

        // Blackoil case
        if( phases.size() == 3 ) {
            Opm::flowEbosBlackoilSetDeck(externalSetupTimer.elapsed(), *deck_, *eclipseState_, *schedule_, *summaryConfig_);
            return true;
        }
        else
        {
            if (outputCout_)
                std::cerr << "Only blackoil configuration is supported!" << std::endl;
            return false;
        }
    }


    int argc_;
    char **argv_;
    bool outputCout_; // copy of EWOMS parameter "EnableTerminalOutput"
    bool outputFiles_; // output files?
    bool hasRunInit_{false};

    std::shared_ptr<Opm::Deck>          deck_;
    std::shared_ptr<Opm::EclipseState>  eclipseState_;
    std::shared_ptr<Opm::Schedule>      schedule_;
    std::shared_ptr<Opm::SummaryConfig> summaryConfig_;

    std::unique_ptr<Opm::FlowMainEbos<TTAG(EclFlowProblem)>> mainfunc_;
};

PYBIND11_MODULE(simulators, m)
{
    py::class_<BlackOilSimulator>(m, "BlackOilSimulator")
        .def(py::init< const std::string& >())
        .def(py::init< const Opm::Deck&, const Opm::EclipseState&, const Opm::Schedule&, const Opm::SummaryConfig& >())
        .def("run", &BlackOilSimulator::run)
        .def("set_deck", &BlackOilSimulator::setDeck)
        .def("set_eclipse_state", &BlackOilSimulator::setEclipseState)
        .def("set_schedule", &BlackOilSimulator::setSchedule)
        .def("set_summary_config", &BlackOilSimulator::setSummaryConfig)
        .def("step_init", &BlackOilSimulator::step_init);
}
