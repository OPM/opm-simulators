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

# ifndef FLOW_BLACKOIL_ONLY
#  include <flow/flow_ebos_gasoil.hpp>
#  include <flow/flow_ebos_oilwater.hpp>
#  include <flow/flow_ebos_solvent.hpp>
#  include <flow/flow_ebos_polymer.hpp>
#  include <flow/flow_ebos_foam.hpp>
#  include <flow/flow_ebos_brine.hpp>
#  include <flow/flow_ebos_oilwater_brine.hpp>
#  include <flow/flow_ebos_energy.hpp>
#  include <flow/flow_ebos_oilwater_polymer.hpp>
#  include <flow/flow_ebos_oilwater_polymer_injectivity.hpp>
# endif

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/OpmLog/EclipsePRTLog.hpp>
#include <opm/common/OpmLog/LogUtil.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/checkDeck.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/ArrayDimChecker.hpp>

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>

#include <opm/simulators/flow/FlowMainEbos.hpp>
#include <opm/simulators/flow/MissingFeatures.hpp>


#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

#if HAVE_MPI
#include <opm/simulators/utils/ParallelEclipseState.hpp>
#include <opm/simulators/utils/ParallelSerialization.hpp>
#endif

#include <string>
#include <type_traits>

namespace Opm::Properties {

// this is a dummy type tag that is used to setup the parameters before the actual
// simulator.
NEW_TYPE_TAG(FlowEarlyBird, INHERITS_FROM(EclFlowProblem));

} // namespace Opm::Properties

namespace Opm {
  template <class TypeTag>
  void flowEbosSetDeck(Deck *deck, EclipseState& eclState, Schedule& schedule, SummaryConfig& summaryConfig)
  {
    using Vanguard = typename GET_PROP_TYPE(TypeTag, Vanguard);
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


namespace Opm
{
    // ----------------- Main class -----------------
    //   For now, we will either be instantiated from main() in flow.cpp,
    //   or from a Python pybind11 module..
    // NOTE (March 2020): When used from a pybind11 module, we do not neccessarily
    //   want to run the whole simulation by calling run(), it is also
    //   useful to just run one report step at a time. According to these different
    //   usage scenarios, we refactored the original run() in flow.cpp into this class.
    class Main
    {
    private:
        using FlowMainEbosType = Opm::FlowMainEbos<TTAG(EclFlowProblem)>;
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

        Main(const std::string &filename)
        {
            deckFilename_.assign(filename);
            flowProgName_.assign("flow");
            argc_ = 2;
            saveArgs_[0] = const_cast<char *>(flowProgName_.c_str());
            saveArgs_[1] = const_cast<char *>(deckFilename_.c_str());
            argv_ = saveArgs_;
        }

        Main(int argc,
             char** argv,
             std::shared_ptr<Opm::Deck> deck,
             std::shared_ptr<Opm::EclipseState> eclipseState,
             std::shared_ptr<Opm::Schedule> schedule,
             std::shared_ptr<Opm::SummaryConfig> summaryConfig)
            : argc_(argc)
            , argv_(argv)
            , deck_(deck)
            , eclipseState_(eclipseState)
            , schedule_(schedule)
            , summaryConfig_(summaryConfig)
        {
        }

        int runDynamic()
        {
            int exitCode = EXIT_SUCCESS;
            if (initialize_<TTAG(FlowEarlyBird)>(exitCode)) {
                return dispatchDynamic_();
            } else {
                return exitCode;
            }
        }

        template <class TypeTag>
        int runStatic()
        {
            int exitCode = EXIT_SUCCESS;
            if (initialize_<TypeTag>(exitCode)) {
                return dispatchStatic_<TypeTag>();
            } else {
                return exitCode;
            }
        }

        // To be called from the Python interface code. Only do the
        // initialization and then return a pointer to the FlowEbosMain
        // object that can later be accessed directly from the Python interface
        // to e.g. advance the simulator one report step
        std::unique_ptr<FlowMainEbosType> initFlowEbosBlackoil(int& exitCode)
        {
            exitCode = EXIT_SUCCESS;
            if (initialize_<TTAG(FlowEarlyBird)>(exitCode)) {
                // TODO: check that this deck really represents a blackoil
                // case. E.g. check that number of phases == 3
                Opm::flowEbosBlackoilSetDeck(
                    setupTime_,
                    deck_.get(),
                    *eclipseState_,
                    *schedule_,
                    *summaryConfig_);
                return Opm::flowEbosBlackoilMainInit(
                    argc_, argv_, outputCout_, outputFiles_);
            } else {
                //NOTE: exitCode was set by initialize_() above;
                return std::unique_ptr<FlowMainEbosType>(); // nullptr
            }
        }

    private:
        int dispatchDynamic_()
        {
            const auto& phases = eclipseState_->runspec().phases();
            // run the actual simulator
            //
            // TODO: make sure that no illegal combinations like thermal and twophase are
            //       requested.

            if ( false ) {}
#ifndef FLOW_BLACKOIL_ONLY
            // Twophase cases
            else if( phases.size() == 2 ) {
                // oil-gas
                if (phases.active( Opm::Phase::GAS )) {
                    Opm::flowEbosGasOilSetDeck(setupTime_, deck_.get(), *eclipseState_, *schedule_, *summaryConfig_);
                    return Opm::flowEbosGasOilMain(argc_, argv_, outputCout_, outputFiles_);
                }
                // oil-water
                else if ( phases.active( Opm::Phase::WATER ) ) {
                    Opm::flowEbosOilWaterSetDeck(setupTime_, deck_.get(), *eclipseState_, *schedule_, *summaryConfig_);
                    return Opm::flowEbosOilWaterMain(argc_, argv_, outputCout_, outputFiles_);
                }
                else {
                    if (outputCout_)
                        std::cerr << "No suitable configuration found, valid are Twophase (oilwater and oilgas), polymer, solvent, or blackoil" << std::endl;
                    return EXIT_FAILURE;
                }
            }
            // Polymer case
            else if ( phases.active( Opm::Phase::POLYMER ) ) {
                if ( !phases.active( Opm::Phase::WATER) ) {
                    if (outputCout_)
                        std::cerr << "No valid configuration is found for polymer simulation, valid options include "
                                  << "oilwater + polymer and blackoil + polymer" << std::endl;
                    return EXIT_FAILURE;
                }

                // Need to track the polymer molecular weight
                // for the injectivity study
                if ( phases.active( Opm::Phase::POLYMW ) ) {
                    // only oil water two phase for now
                    assert( phases.size() == 4);
                    return Opm::flowEbosOilWaterPolymerInjectivityMain(argc_, argv_, outputCout_, outputFiles_);
                }

                if ( phases.size() == 3 ) { // oil water polymer case
                    Opm::flowEbosOilWaterPolymerSetDeck(setupTime_, deck_.get(), *eclipseState_, *schedule_, *summaryConfig_);
                    return Opm::flowEbosOilWaterPolymerMain(argc_, argv_, outputCout_, outputFiles_);
                } else {
                    Opm::flowEbosPolymerSetDeck(setupTime_, deck_.get(), *eclipseState_, *schedule_, *summaryConfig_);
                    return Opm::flowEbosPolymerMain(argc_, argv_, outputCout_, outputFiles_);
                }
            }
            // Foam case
            else if ( phases.active( Opm::Phase::FOAM ) ) {
                Opm::flowEbosFoamSetDeck(setupTime_, deck_.get(), *eclipseState_, *schedule_, *summaryConfig_);
                return Opm::flowEbosFoamMain(argc_, argv_, outputCout_, outputFiles_);
            }
            // Brine case
            else if ( phases.active( Opm::Phase::BRINE ) ) {
                if ( !phases.active( Opm::Phase::WATER) ) {
                    if (outputCout_)
                        std::cerr << "No valid configuration is found for brine simulation, valid options include "
                                  << "oilwater + brine and blackoil + brine" << std::endl;
                    return EXIT_FAILURE;
                }
                if ( phases.size() == 3 ) { // oil water brine case
                    Opm::flowEbosOilWaterBrineSetDeck(setupTime_, deck_.get(), *eclipseState_, *schedule_, *summaryConfig_);
                    return Opm::flowEbosOilWaterBrineMain(argc_, argv_, outputCout_, outputFiles_);
                } else {
                    Opm::flowEbosBrineSetDeck(setupTime_, deck_.get(), *eclipseState_, *schedule_, *summaryConfig_);
                    return Opm::flowEbosBrineMain(argc_, argv_, outputCout_, outputFiles_);
                }
            }
            // Solvent case
            else if ( phases.active( Opm::Phase::SOLVENT ) ) {
                Opm::flowEbosSolventSetDeck(setupTime_, deck_.get(), *eclipseState_, *schedule_, *summaryConfig_);
                return Opm::flowEbosSolventMain(argc_, argv_, outputCout_, outputFiles_);
            }
            // Energy case
            else if (eclipseState_->getSimulationConfig().isThermal()) {
                Opm::flowEbosEnergySetDeck(setupTime_, deck_.get(), *eclipseState_, *schedule_, *summaryConfig_);
                return Opm::flowEbosEnergyMain(argc_, argv_, outputCout_, outputFiles_);
            }
#endif // FLOW_BLACKOIL_ONLY
            // Blackoil case
            else if( phases.size() == 3 ) {
                Opm::flowEbosBlackoilSetDeck(setupTime_, deck_.get(), *eclipseState_, *schedule_, *summaryConfig_);
                return Opm::flowEbosBlackoilMain(argc_, argv_, outputCout_, outputFiles_);
            }
            else {
                if (outputCout_)
                    std::cerr << "No suitable configuration found, valid are Twophase, polymer, foam, brine, solvent, energy, blackoil." << std::endl;
                return EXIT_FAILURE;
            }
        }

        template <class TypeTag>
        int dispatchStatic_()
        {
            Opm::flowEbosSetDeck<TypeTag>(deck_.get(), *eclipseState_, *schedule_, *summaryConfig_);
            return Opm::flowEbosMain<TypeTag>(argc_, argv_, outputCout_, outputFiles_);
        }

        /// \brief Initialize
        /// \param exitCode The exitCode of the program.
        /// \return Whether to actually run the simulator. I.e. true if parsing of command line
        /// was successful and no --help, --print-properties, or --print-parameters have been found.
        template <class TypeTagEarlyBird>
        bool initialize_(int& exitCode)
        {
            Dune::Timer externalSetupTimer;
            externalSetupTimer.start();

            handleVersionCmdLine_(argc_, argv_);
            // MPI setup.
#if HAVE_DUNE_FEM
            Dune::Fem::MPIManager::initialize(argc_, argv_);
            int mpiRank = Dune::Fem::MPIManager::rank();
#else
            // the design of the plain dune MPIHelper class is quite flawed: there is no way to
            // get the instance without having the argc and argv parameters available and it is
            // not possible to determine the MPI rank and size without an instance. (IOW: the
            // rank() and size() methods are supposed to be static.)
            const auto& mpiHelper = Dune::MPIHelper::instance(argc_, argv_);
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
            typedef TypeTagEarlyBird PreTypeTag;
            typedef typename GET_PROP_TYPE(PreTypeTag, Problem) PreProblem;

            PreProblem::setBriefDescription("Flow, an advanced reservoir simulator for ECL-decks provided by the Open Porous Media project.");
            int status = Opm::FlowMainEbos<PreTypeTag>::setupParameters_(argc_, argv_);
            if (status != 0) {
                // if setupParameters_ returns a value smaller than 0, there was no error, but
                // the program should abort. This is the case e.g. for the --help and the
                // --print-properties parameters.
#if HAVE_MPI
                if (status < 0)
                    MPI_Finalize(); // graceful stop for --help or --print-properties command line.
                else
                    MPI_Abort(MPI_COMM_WORLD, status);
#endif
                exitCode = (status > 0) ? status : EXIT_SUCCESS;
                return false; //  Whether to run the simulator
            }

            FileOutputMode outputMode = FileOutputMode::OUTPUT_NONE;
            outputCout_ = false;
            if (mpiRank == 0)
                outputCout_ = EWOMS_GET_PARAM(PreTypeTag, bool, EnableTerminalOutput);

            std::string deckFilename;
            std::string outputDir;
            if ( eclipseState_ ) {
                deckFilename = eclipseState_->getIOConfig().fullBasePath();
                outputDir = eclipseState_->getIOConfig().getOutputDir();
            }
            else {
                deckFilename = EWOMS_GET_PARAM(PreTypeTag, std::string, EclDeckFileName);
            }

            typedef typename GET_PROP_TYPE(PreTypeTag, Vanguard) PreVanguard;
            try {
                deckFilename = PreVanguard::canonicalDeckPath(deckFilename).string();
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
            if (outputCout_) {
                Opm::FlowMainEbos<PreTypeTag>::printBanner();
            }
            // Create Deck and EclipseState.
            try {
                if (outputCout_) {
                    std::cout << "Reading deck file '" << deckFilename << "'\n";
                    std::cout.flush();
                }
                auto python = std::make_shared<Opm::Python>();
                {
                    Opm::Parser parser;
                    Opm::ParseContext parseContext({{Opm::ParseContext::PARSE_RANDOM_SLASH, Opm::InputError::IGNORE},
                                                    {Opm::ParseContext::PARSE_MISSING_DIMS_KEYWORD, Opm::InputError::WARN},
                                                    {Opm::ParseContext::SUMMARY_UNKNOWN_WELL, Opm::InputError::WARN},
                                                    {Opm::ParseContext::SUMMARY_UNKNOWN_GROUP, Opm::InputError::WARN}});
                    Opm::ErrorGuard errorGuard;
                    if (outputDir.empty())
                        outputDir = EWOMS_GET_PARAM(PreTypeTag, std::string, OutputDir);
                    outputMode = setupLogging_(mpiRank,
                                      deckFilename,
                                      outputDir,
                                      EWOMS_GET_PARAM(PreTypeTag, std::string, OutputMode),
                                      outputCout_, "STDOUT_LOGGER");

                    if (EWOMS_GET_PARAM(PreTypeTag, bool, EclStrictParsing))
                        parseContext.update( Opm::InputError::DELAYED_EXIT1);

                    Opm::FlowMainEbos<PreTypeTag>::printPRTHeader(outputCout_);

#if HAVE_MPI
                    int parseSuccess = 0;
#endif
                    std::string failureMessage;

                    if (mpiRank == 0) {
                        try
                        {
                            if (!deck_)
                                deck_.reset( new Opm::Deck( parser.parseFile(deckFilename , parseContext, errorGuard)));
                            Opm::MissingFeatures::checkKeywords(*deck_, parseContext, errorGuard);
                            if ( outputCout_ )
                                Opm::checkDeck(*deck_, parser, parseContext, errorGuard);

                            if (!eclipseState_) {
#if HAVE_MPI
                                eclipseState_.reset(new Opm::ParallelEclipseState(*deck_));
#else
                                eclipseState_.reset(new Opm::EclipseState(*deck_));
#endif
                            }
                            /*
                              For the time being initializing wells and groups from the
                              restart file is not possible, but work is underways and it is
                              included here as a switch.
                            */
                            const bool init_from_restart_file = !EWOMS_GET_PARAM(PreTypeTag, bool, SchedRestart);
                            const auto& init_config = eclipseState_->getInitConfig();
                            if (init_config.restartRequested() && init_from_restart_file) {
                                int report_step = init_config.getRestartStep();
                                const auto& rst_filename = eclipseState_->getIOConfig().getRestartFileName( init_config.getRestartRootName(), report_step, false );
                                Opm::EclIO::ERst rst_file(rst_filename);
                                const auto& rst_state = Opm::RestartIO::RstState::load(rst_file, report_step);
                                if (!schedule_)
                                    schedule_.reset(new Opm::Schedule(*deck_, *eclipseState_, parseContext, errorGuard, python, &rst_state) );
                            }
                            else {
                                if (!schedule_)
                                    schedule_.reset(new Opm::Schedule(*deck_, *eclipseState_, parseContext, errorGuard, python));
                            }
                            setupMessageLimiter_(schedule_->getMessageLimits(), "STDOUT_LOGGER");
                            if (!summaryConfig_)
                                summaryConfig_.reset( new Opm::SummaryConfig(*deck_, *schedule_, eclipseState_->getTableManager(), parseContext, errorGuard));
#if HAVE_MPI
                            parseSuccess = 1;
#endif
                        }
                        catch(const std::exception& e)
                        {
                            failureMessage = e.what();
                        }
                    }
#if HAVE_MPI
                    else {
                        if (!summaryConfig_)
                            summaryConfig_.reset(new Opm::SummaryConfig);
                        if (!schedule_)
                            schedule_.reset(new Opm::Schedule(python));
                        if (!eclipseState_)
                            eclipseState_.reset(new Opm::ParallelEclipseState);
                    }

                    auto comm = Dune::MPIHelper::getCollectiveCommunication();
                    parseSuccess = comm.max(parseSuccess);
                    if (!parseSuccess)
                    {
                        if (errorGuard) {
                            errorGuard.dump();
                            errorGuard.clear();
                        }
                        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
                    }

                    Opm::eclStateBroadcast(*eclipseState_, *schedule_, *summaryConfig_);
#endif

                    Opm::checkConsistentArrayDimensions(*eclipseState_, *schedule_, parseContext, errorGuard);

                    if (errorGuard) {
                        errorGuard.dump();
                        errorGuard.clear();

                        throw std::runtime_error("Unrecoverable errors were encountered while loading input.");
                    }
                }
                setupTime_ = externalSetupTimer.elapsed();
                outputFiles_ = (outputMode != FileOutputMode::OUTPUT_NONE);
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

            exitCode = EXIT_SUCCESS;
            return true;
        }

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

        int argc_;
        char** argv_;
        bool outputCout_;
        bool outputFiles_;
        double setupTime_;
        std::string deckFilename_;
        std::string flowProgName_;
        char *saveArgs_[2];
        std::shared_ptr<Opm::Deck> deck_;
        std::shared_ptr<Opm::EclipseState> eclipseState_;
        std::shared_ptr<Opm::Schedule> schedule_;
        std::shared_ptr<Opm::SummaryConfig> summaryConfig_;
    };

} // namespace Opm

#endif // OPM_MAIN_HEADER_INCLUDED
