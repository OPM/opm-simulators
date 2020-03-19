/*
  Copyright 2013, 2014, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015, 2017 IRIS AS

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
#ifndef FLOW_TAG_HPP
#define FLOW_TAG_HPP

#include <opm/simulators/flow/FlowMainEbos.hpp>
#include <opm/simulators/flow/MissingFeatures.hpp>
#include <opm/simulators/flow/SimulatorFullyImplicitBlackoilEbos.hpp>
#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/common/utility/parameters/ParameterGroup.hpp>
#include <opm/material/common/ResetLocale.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/OpmLog/EclipsePRTLog.hpp>
#include <opm/common/OpmLog/LogUtil.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/checkDeck.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/ArrayDimChecker.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/UDQ/UDQAssign.hpp>

//#include <opm/material/fluidsystems/BlackOilFluidSystemSimple.hpp>
//#include <opm/material/fluidsystems/BlackOilFluidSystemSimple.hpp>
#include <opm/models/blackoil/blackoilintensivequantities.hh>
#include <opm/material/fluidstates/BlackOilFluidState.hpp>
//#include <opm/material/fluidstates/BlackOilFluidStateSimple.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

#if HAVE_MPI
#include <opm/simulators/utils/ParallelEclipseState.hpp>
#include <opm/simulators/utils/ParallelSerialization.hpp>
#endif


BEGIN_PROPERTIES

// this is a dummy type tag that is used to setup the parameters before the actual
// simulator.
NEW_TYPE_TAG(FlowEarlyBird, INHERITS_FROM(EclFlowProblem));

END_PROPERTIES




namespace Opm {
  template <class TypeTag>
  void flowEbosSetDeck(Deck *deck, EclipseState& eclState, Schedule& schedule, SummaryConfig& summaryConfig)
  {
    typedef typename GET_PROP_TYPE(TypeTag, Vanguard) Vanguard;
    Vanguard::setExternalDeck(deck);
    Vanguard::setExternalEclState(&eclState);
    Vanguard::setExternalSchedule(&schedule);
    Vanguard::setExternalSummaryConfig(&summaryConfig);
  }

// ----------------- Main program -----------------
  template <class TypeTag>
  int flowEbosMain(int argc, char** argv, bool outputCout, bool outputFiles)
  {
    // we always want to use the default locale, and thus spare us the trouble
    // with incorrect locale settings.
    Opm::resetLocale();

#if HAVE_DUNE_FEM
    Dune::Fem::MPIManager::initialize(argc, argv);
#else
    Dune::MPIHelper::instance(argc, argv);
#endif
    Opm::FlowMainEbos<TypeTag> mainfunc;
    return mainfunc.execute(argc, argv, outputCout, outputFiles);
  }

}





namespace detail
{
    Opm::filesystem::path simulationCaseName( const std::string& casename ) {
        namespace fs = Opm::filesystem;

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
    if (!Opm::filesystem::is_directory(cmdline_output_dir)) {
        try {
            Opm::filesystem::create_directories(cmdline_output_dir);
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
    using Opm::filesystem::path;
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

    std::shared_ptr<Opm::StreamLog> streamLog = std::make_shared<Opm::StreamLog>(std::cout, Opm::Log::StdoutMessageTypes);
    Opm::OpmLog::addBackend(stdout_log_id, streamLog);
    streamLog->setMessageFormatter(std::make_shared<Opm::SimpleMessageFormatter>(true));

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


// ----------------- Main program -----------------
template<class TypeTag>
int mainFlow(int argc, char** argv)
{
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
    bool outputCout = false;
    if (mpiRank == 0)
        outputCout = EWOMS_GET_PARAM(PreTypeTag, bool, EnableTerminalOutput);

    std::string deckFilename = EWOMS_GET_PARAM(PreTypeTag, std::string, EclDeckFileName);
    typedef typename GET_PROP_TYPE(PreTypeTag, Vanguard) PreVanguard;
    try {
        deckFilename = PreVanguard::canonicalDeckPath(deckFilename).string();
    }
    catch (const std::exception& e) {
        if (mpiRank == 0)
        {
            Opm::Parameters::printUsage<PreTypeTag>(PreProblem::helpPreamble(argc, const_cast<const char**>(argv)),
                                                    e.what());
        }
#if HAVE_MPI
        MPI_Finalize();
#endif
        return 1;
    }

    // Create Deck and EclipseState.
    try {
        if (outputCout) {
            std::cout << "Reading deck file '" << deckFilename << "'\n";
            std::cout.flush();
        }
        std::shared_ptr<Opm::Deck> deck;
        std::shared_ptr<Opm::EclipseState> eclipseState;
        std::shared_ptr<Opm::Schedule> schedule;
        std::shared_ptr<Opm::SummaryConfig> summaryConfig;
        {
            Opm::Parser parser;
            Opm::ParseContext parseContext;
            Opm::ErrorGuard errorGuard;
            outputMode = setupLogging(mpiRank,
                                      deckFilename,
                                      EWOMS_GET_PARAM(PreTypeTag, std::string, OutputDir),
                                      EWOMS_GET_PARAM(PreTypeTag, std::string, OutputMode),
                                      outputCout, "STDOUT_LOGGER");

            if (EWOMS_GET_PARAM(PreTypeTag, bool, EclStrictParsing))
                parseContext.update( Opm::InputError::DELAYED_EXIT1);
            else {
                parseContext.update(Opm::ParseContext::PARSE_RANDOM_SLASH, Opm::InputError::IGNORE);
                parseContext.update(Opm::ParseContext::PARSE_MISSING_DIMS_KEYWORD, Opm::InputError::WARN);
                parseContext.update(Opm::ParseContext::SUMMARY_UNKNOWN_WELL, Opm::InputError::WARN);
                parseContext.update(Opm::ParseContext::SUMMARY_UNKNOWN_GROUP, Opm::InputError::WARN);
            }

            Opm::FlowMainEbos<PreTypeTag>::printPRTHeader(outputCout);

            if (mpiRank == 0) {
                deck.reset( new Opm::Deck( parser.parseFile(deckFilename , parseContext, errorGuard)));
                Opm::MissingFeatures::checkKeywords(*deck, parseContext, errorGuard);
                if ( outputCout )
                    Opm::checkDeck(*deck, parser, parseContext, errorGuard);

#if HAVE_MPI
                eclipseState.reset(new Opm::ParallelEclipseState(*deck));
#else
                eclipseState.reset(new Opm::EclipseState(*deck));
#endif
                schedule.reset(new Opm::Schedule(*deck, *eclipseState, parseContext, errorGuard));
                setupMessageLimiter(schedule->getMessageLimits(), "STDOUT_LOGGER");
                summaryConfig.reset( new Opm::SummaryConfig(*deck, *schedule, eclipseState->getTableManager(), parseContext, errorGuard));
            }
#if HAVE_MPI
            else {
                summaryConfig.reset(new Opm::SummaryConfig);
                schedule.reset(new Opm::Schedule);
                eclipseState.reset(new Opm::ParallelEclipseState);
            }
            Opm::eclStateBroadcast(*eclipseState, *schedule, *summaryConfig);
#endif

            Opm::checkConsistentArrayDimensions(*eclipseState, *schedule, parseContext, errorGuard);

            if (errorGuard) {
                errorGuard.dump();
                errorGuard.clear();

                throw std::runtime_error("Unrecoverable errors were encountered while loading input.");
            }
        }
        bool outputFiles = (outputMode != FileOutputMode::OUTPUT_NONE);
        Opm::flowEbosSetDeck<TypeTag>(deck.get(), *eclipseState, *schedule, *summaryConfig);
        return Opm::flowEbosMain<TypeTag>(argc, argv, outputCout, outputFiles);
    }
  catch (const std::invalid_argument& e)
    {
      if (outputCout) {
        std::cerr << "Failed to create valid EclipseState object." << std::endl;
        std::cerr << "Exception caught: " << e.what() << std::endl;
      }
      throw;
    }

    return EXIT_SUCCESS;
}
#endif
