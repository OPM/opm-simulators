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
#ifndef OPM_MAIN_HEADER_INCLUDED
#define OPM_MAIN_HEADER_INCLUDED

#include <flow/flow_ebos_blackoil.hpp>

#ifdef OPM_FLOW_MAIN
# ifndef FLOW_BLACKOIL_ONLY
#  include <flow/flow_ebos_gasoil.hpp>
#  include <flow/flow_ebos_oilwater.hpp>
#  include <flow/flow_ebos_solvent.hpp>
#  include <flow/flow_ebos_polymer.hpp>
#  include <flow/flow_ebos_foam.hpp>
#  include <flow/flow_ebos_brine.hpp>
#  include <flow/flow_ebos_energy.hpp>
#  include <flow/flow_ebos_oilwater_polymer.hpp>
#  include <flow/flow_ebos_oilwater_polymer_injectivity.hpp>
# endif
# include <opm/simulators/utils/moduleVersion.hpp>
# include <opm/io/eclipse/rst/state.hpp>
# include <opm/io/eclipse/ERst.hpp>
# include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
# include <opm/parser/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#else  /* OPM_FLOW_MAIN */
# include <opm/simulators/flow/SimulatorFullyImplicitBlackoilEbos.hpp>
# include <opm/common/utility/parameters/ParameterGroup.hpp>
# include <opm/material/common/ResetLocale.hpp>
# include <opm/parser/eclipse/Python/Python.hpp>
# include <opm/parser/eclipse/Parser/ParseContext.hpp>
# include <opm/parser/eclipse/Parser/ErrorGuard.hpp>
# include <opm/parser/eclipse/EclipseState/Schedule/UDQ/UDQAssign.hpp>
//# include <opm/material/fluidsystems/BlackOilFluidSystemSimple.hpp>
//# include <opm/material/fluidsystems/BlackOilFluidSystemSimple.hpp>
# include <opm/models/blackoil/blackoilintensivequantities.hh>
# include <opm/material/fluidstates/BlackOilFluidState.hpp>
//# include <opm/material/fluidstates/BlackOilFluidStateSimple.hpp>
#endif

#include <opm/simulators/flow/FlowMainEbos.hpp>
#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/simulators/flow/MissingFeatures.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/OpmLog/EclipsePRTLog.hpp>
#include <opm/common/OpmLog/LogUtil.hpp>


#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/checkDeck.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/ArrayDimChecker.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

#if HAVE_MPI
#include <opm/simulators/utils/ParallelEclipseState.hpp>
#include <opm/simulators/utils/ParallelSerialization.hpp>
#endif

#include <type_traits>

BEGIN_PROPERTIES

// this is a dummy type tag that is used to setup the parameters before the actual
// simulator.
NEW_TYPE_TAG(FlowEarlyBird, INHERITS_FROM(EclFlowProblem));

END_PROPERTIES

#ifndef OPM_FLOW_MAIN
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

# if HAVE_DUNE_FEM
    Dune::Fem::MPIManager::initialize(argc, argv);
# else
    Dune::MPIHelper::instance(argc, argv);
# endif
    Opm::FlowMainEbos<TypeTag> mainfunc;
    return mainfunc.execute(argc, argv, outputCout, outputFiles);
  }

}

#endif  /* #ifndef OPM_FLOW_MAIN */

namespace Opm
{
    // ----------------- Main class -----------------
    //   For now, we will either be instantiated from main() in flow.cpp,
    //   or from a Python pybind11 module..
    // NOTE (March 2020): When used from a pybind11 module, we do not neccessarily
    //   want to run the whole simulation by calling run(), it is also
    //   useful to just run one report step at a time. According to these different
    //   usage scenarios, we refactored the original run() in flow.cpp into this class.
    template <class TypeTag>
    class Main
    {
    private:
        enum class FileOutputMode {
            //! \brief No output to files.
            OUTPUT_NONE = 0,
            //! \brief Output only to log files, no eclipse output.
            OUTPUT_LOG_ONLY = 1,
            //! \brief Output to all files.
            OUTPUT_ALL = 3
        };
    public:
        Main(int argc, char** argv) : argc_(argc), argv_(argv)  {  }

        int run() {
            return main_(argc_, argv_);
        }

    private:

        Opm::filesystem::path simulationCaseName_( const std::string& casename ) {
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
        void handleVersionCmdLine_(int argc, char** argv) {
            for ( int i = 1; i < argc; ++i )
            {
                if (std::strcmp(argv[i], "--version") == 0) {
                    std::cout << "flow " << Opm::moduleVersionName() << std::endl;
                    std::exit(EXIT_SUCCESS);
                }
            }
        }


        void ensureOutputDirExists_(const std::string& cmdline_output_dir)
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
        FileOutputMode setupLogging_(int mpi_rank_, const std::string& deck_filename, const std::string& cmdline_output_dir, const std::string& cmdline_output, bool output_cout_, const std::string& stdout_log_id) {

            if (!cmdline_output_dir.empty()) {
                ensureOutputDirExists_(cmdline_output_dir);
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
                output_dir = fpath.has_parent_path()
                    ? absolute(fpath.parent_path()).generic_string()
                    : Opm::filesystem::current_path().generic_string();
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



        void setupMessageLimiter_(const Opm::MessageLimits msgLimits,  const std::string& stdout_log_id) {
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
        int main_(int argc, char** argv)
        {
            Dune::Timer externalSetupTimer;
            externalSetupTimer.start();

            handleVersionCmdLine_(argc, argv);
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
                if ( mpiRank == 0 )
#ifdef OPM_FLOW_MAIN
                    std::cerr << "Exception received: " << e.what() << ". Try '--help' for a usage description.\n";
#else
                    Opm::Parameters::printUsage<PreTypeTag>(PreProblem::helpPreamble(
                            argc, const_cast<const char**>(argv)), e.what());
#endif
#if HAVE_MPI
                MPI_Finalize();
#endif
                return 1;
            }
#ifdef OPM_FLOW_MAIN
            if (outputCout) {
                Opm::FlowMainEbos<PreTypeTag>::printBanner();
            }
#endif
            // Create Deck and EclipseState.
            try {
                if (outputCout) {
                    std::cout << "Reading deck file '" << deckFilename << "'\n";
                    std::cout.flush();
                }
                auto python = std::make_shared<Opm::Python>();
                std::shared_ptr<Opm::Deck> deck;
                std::shared_ptr<Opm::EclipseState> eclipseState;
                std::shared_ptr<Opm::Schedule> schedule;
                std::shared_ptr<Opm::SummaryConfig> summaryConfig;
                {
                    Opm::Parser parser;
                    Opm::ParseContext parseContext({{Opm::ParseContext::PARSE_RANDOM_SLASH, Opm::InputError::IGNORE},
                                                    {Opm::ParseContext::PARSE_MISSING_DIMS_KEYWORD, Opm::InputError::WARN},
                                                    {Opm::ParseContext::SUMMARY_UNKNOWN_WELL, Opm::InputError::WARN},
                                                    {Opm::ParseContext::SUMMARY_UNKNOWN_GROUP, Opm::InputError::WARN}});
                    Opm::ErrorGuard errorGuard;
                    outputMode = setupLogging_(mpiRank,
                                      deckFilename,
                                      EWOMS_GET_PARAM(PreTypeTag, std::string, OutputDir),
                                      EWOMS_GET_PARAM(PreTypeTag, std::string, OutputMode),
                                      outputCout, "STDOUT_LOGGER");

                    if (EWOMS_GET_PARAM(PreTypeTag, bool, EclStrictParsing))
                        parseContext.update( Opm::InputError::DELAYED_EXIT1);

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
#ifdef OPM_FLOW_MAIN
                        /*
                          For the time being initializing wells and groups from the
                          restart file is not possible, but work is underways and it is
                          included here as a switch.
                        */
                        const bool init_from_restart_file = !EWOMS_GET_PARAM(PreTypeTag, bool, SchedRestart);
                        const auto& init_config = eclipseState->getInitConfig();
                        if (init_config.restartRequested() && init_from_restart_file) {
                            int report_step = init_config.getRestartStep();
                            const auto& rst_filename = eclipseState->getIOConfig().getRestartFileName( init_config.getRestartRootName(), report_step, false );
                            Opm::EclIO::ERst rst_file(rst_filename);
                            const auto& rst_state = Opm::RestartIO::RstState::load(rst_file, report_step);
                            schedule.reset(new Opm::Schedule(*deck, *eclipseState, parseContext, errorGuard, python, &rst_state) );
                        } else
                            schedule.reset(new Opm::Schedule(*deck, *eclipseState, parseContext, errorGuard, python));
#else
                        schedule.reset(new Opm::Schedule(*deck, *eclipseState, parseContext, errorGuard, python));
#endif /*OPM_FLOW_MAIN */
                        setupMessageLimiter_(schedule->getMessageLimits(), "STDOUT_LOGGER");
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
#ifdef OPM_FLOW_MAIN
                const auto& phases = eclipseState->runspec().phases();
                // run the actual simulator
                //
                // TODO: make sure that no illegal combinations like thermal and twophase are
                //       requested.

                if ( false ) {}
#ifndef FLOW_BLACKOIL_ONLY
                // Twophase cases
                else if( phases.size() == 2 ) {
                    // oil-gas
                    if (phases.active( Opm::Phase::GAS ))
                    {
                        Opm::flowEbosGasOilSetDeck(externalSetupTimer.elapsed(), deck.get(), *eclipseState, *schedule, *summaryConfig);
                        return Opm::flowEbosGasOilMain(argc, argv, outputCout, outputFiles);
                    }
                    // oil-water
                    else if ( phases.active( Opm::Phase::WATER ) )
                    {
                        Opm::flowEbosOilWaterSetDeck(externalSetupTimer.elapsed(), deck.get(), *eclipseState, *schedule, *summaryConfig);
                        return Opm::flowEbosOilWaterMain(argc, argv, outputCout, outputFiles);
                    }
                    else {
                        if (outputCout)
                            std::cerr << "No suitable configuration found, valid are Twophase (oilwater and oilgas), polymer, solvent, or blackoil" << std::endl;
                        return EXIT_FAILURE;
                    }
                }
                // Polymer case
                else if ( phases.active( Opm::Phase::POLYMER ) ) {

                    if ( !phases.active( Opm::Phase::WATER) ) {
                        if (outputCout)
                            std::cerr << "No valid configuration is found for polymer simulation, valid options include "
                              << "oilwater + polymer and blackoil + polymer" << std::endl;
                        return EXIT_FAILURE;
                    }

                    // Need to track the polymer molecular weight
                    // for the injectivity study
                    if ( phases.active( Opm::Phase::POLYMW ) ) {
                        // only oil water two phase for now
                        assert( phases.size() == 4);
                        return Opm::flowEbosOilWaterPolymerInjectivityMain(argc, argv, outputCout, outputFiles);
                    }

                    if ( phases.size() == 3 ) { // oil water polymer case
                        Opm::flowEbosOilWaterPolymerSetDeck(externalSetupTimer.elapsed(), deck.get(), *eclipseState, *schedule, *summaryConfig);
                        return Opm::flowEbosOilWaterPolymerMain(argc, argv, outputCout, outputFiles);
                    } else {
                        Opm::flowEbosPolymerSetDeck(externalSetupTimer.elapsed(), deck.get(), *eclipseState, *schedule, *summaryConfig);
                        return Opm::flowEbosPolymerMain(argc, argv, outputCout, outputFiles);
                    }
                }
                // Foam case
                else if ( phases.active( Opm::Phase::FOAM ) ) {
                    Opm::flowEbosFoamSetDeck(externalSetupTimer.elapsed(), deck.get(), *eclipseState, *schedule, *summaryConfig);
                    return Opm::flowEbosFoamMain(argc, argv, outputCout, outputFiles);
                }
                // Brine case
                else if ( phases.active( Opm::Phase::BRINE ) ) {
                    Opm::flowEbosBrineSetDeck(externalSetupTimer.elapsed(), deck.get(), *eclipseState, *schedule, *summaryConfig);
                    return Opm::flowEbosBrineMain(argc, argv, outputCout, outputFiles);
                }
                // Solvent case
                else if ( phases.active( Opm::Phase::SOLVENT ) ) {
                    Opm::flowEbosSolventSetDeck(externalSetupTimer.elapsed(), deck.get(), *eclipseState, *schedule, *summaryConfig);
                    return Opm::flowEbosSolventMain(argc, argv, outputCout, outputFiles);
                }
                // Energy case
                else if (eclipseState->getSimulationConfig().isThermal()) {
                    Opm::flowEbosEnergySetDeck(externalSetupTimer.elapsed(), deck.get(), *eclipseState, *schedule, *summaryConfig);
                    return Opm::flowEbosEnergyMain(argc, argv, outputCout, outputFiles);
                }
#endif // FLOW_BLACKOIL_ONLY
                // Blackoil case
                else if( phases.size() == 3 ) {
                    Opm::flowEbosBlackoilSetDeck(externalSetupTimer.elapsed(), deck.get(), *eclipseState, *schedule, *summaryConfig);
                    return Opm::flowEbosBlackoilMain(argc, argv, outputCout, outputFiles);
                }
                else
                {
                    if (outputCout)
                        std::cerr << "No suitable configuration found, valid are Twophase, polymer, foam, brine, solvent, energy, blackoil." << std::endl;
                    return EXIT_FAILURE;
                }
#else
            Opm::flowEbosSetDeck<TypeTag>(deck.get(), *eclipseState, *schedule, *summaryConfig);
            return Opm::flowEbosMain<TypeTag>(argc, argv, outputCout, outputFiles);
#endif /* OPM_FLOW_MAIN */
            }
            catch (const std::invalid_argument& e)
            {
                if (outputCout) {
                    std::cerr << "Failed to create valid EclipseState object." << std::endl;
                    std::cerr << "Exception caught: " << e.what() << std::endl;
                }
                return EXIT_FAILURE;
            }

            return EXIT_SUCCESS;
        }

        int argc_;
        char** argv_;
    };
    
} // namespace Opm

#endif // OPM_MAIN_HEADER_INCLUDED
