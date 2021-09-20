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
#  include <flow/flow_ebos_gaswater.hpp>
#  include <flow/flow_ebos_solvent.hpp>
#  include <flow/flow_ebos_polymer.hpp>
#  include <flow/flow_ebos_extbo.hpp>
#  include <flow/flow_ebos_foam.hpp>
#  include <flow/flow_ebos_brine.hpp>
#  include <flow/flow_ebos_oilwater_brine.hpp>
#  include <flow/flow_ebos_energy.hpp>
#  include <flow/flow_ebos_oilwater_polymer.hpp>
#  include <flow/flow_ebos_oilwater_polymer_injectivity.hpp>
# endif

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Parser/ErrorGuard.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/checkDeck.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/ArrayDimChecker.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/UDQ/UDQState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Action/State.hpp>

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>

#include <opm/simulators/flow/FlowMainEbos.hpp>
#include <opm/simulators/utils/readDeck.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

#if HAVE_MPI
#include <opm/simulators/utils/ParallelEclipseState.hpp>
#endif

#include <string>
#include <type_traits>

namespace Opm::Properties {

// this is a dummy type tag that is used to setup the parameters before the actual
// simulator.
namespace TTag {
struct FlowEarlyBird {
    using InheritsFrom = std::tuple<EclFlowProblem>;
};
}

} // namespace Opm::Properties

namespace Opm {
  template <class TypeTag>
  void flowEbosSetDeck(std::unique_ptr<Deck> deck, std::unique_ptr<EclipseState> eclState, std::unique_ptr<Schedule> schedule, std::unique_ptr<SummaryConfig> summaryConfig)
  {
    using Vanguard = GetPropType<TypeTag, Properties::Vanguard>;
    Vanguard::setExternalDeck(std::move(deck));
    Vanguard::setExternalEclState(std::move(eclState));
    Vanguard::setExternalSchedule(std::move(schedule));
    Vanguard::setExternalSummaryConfig(std::move(summaryConfig));
  }

// ----------------- Main program -----------------
  template <class TypeTag>
  int flowEbosMain(int argc, char** argv, bool outputCout, bool outputFiles)
  {
    // we always want to use the default locale, and thus spare us the trouble
    // with incorrect locale settings.
    resetLocale();

    FlowMainEbos<TypeTag> mainfunc(argc, argv, outputCout, outputFiles);
    return mainfunc.execute();
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
        using FlowMainEbosType = FlowMainEbos<Properties::TTag::EclFlowProblem>;

    public:
        Main(int argc, char** argv) : argc_(argc), argv_(argv)  { initMPI();  }

        Main(const std::string &filename)
        {
            deckFilename_.assign(filename);
            flowProgName_.assign("flow");
            argc_ = 2;
            saveArgs_[0] = const_cast<char *>(flowProgName_.c_str());
            saveArgs_[1] = const_cast<char *>(deckFilename_.c_str());
            argv_ = saveArgs_;
            initMPI();
        }

        Main(int argc,
             char** argv,
             std::unique_ptr<Deck> deck,
             std::unique_ptr<EclipseState> eclipseState,
             std::unique_ptr<Schedule> schedule,
             std::unique_ptr<SummaryConfig> summaryConfig)
            : argc_(argc)
            , argv_(argv)
            , deck_(std::move(deck))
            , eclipseState_(std::move(eclipseState))
            , schedule_(std::move(schedule))
            , summaryConfig_(std::move(summaryConfig))
        {
          initMPI();
        }

        ~Main()
        {
            EclGenericVanguard::setCommunication(nullptr);

#if HAVE_MPI && !HAVE_DUNE_FEM
            MPI_Finalize();
#endif
        }

        void initMPI()
        {
#if HAVE_DUNE_FEM
            Dune::Fem::MPIManager::initialize(argc_, argv_);
#elif HAVE_MPI
            MPI_Init(&argc_, &argv_);
#endif
            EclGenericVanguard::setCommunication(std::make_unique<EclGenericVanguard::CommunicationType>());
        }

        int runDynamic()
        {
            int exitCode = EXIT_SUCCESS;
            if (initialize_<Properties::TTag::FlowEarlyBird>(exitCode)) {
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
            if (initialize_<Properties::TTag::FlowEarlyBird>(exitCode)) {
                // TODO: check that this deck really represents a blackoil
                // case. E.g. check that number of phases == 3
                flowEbosBlackoilSetDeck(
                    setupTime_,
                    std::move(deck_),
                    std::move(eclipseState_),
                    std::move(schedule_),
                    std::move(udqState_),
                    std::move(summaryConfig_));
                return flowEbosBlackoilMainInit(
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
                if (phases.active( Phase::OIL ) && phases.active( Phase::GAS )) {
                    flowEbosGasOilSetDeck(setupTime_, std::move(deck_), std::move(eclipseState_),
                                          std::move(schedule_), std::move(summaryConfig_));
                    return flowEbosGasOilMain(argc_, argv_, outputCout_, outputFiles_);
                }
                // oil-water
                else if ( phases.active( Phase::OIL ) && phases.active( Phase::WATER ) ) {
                    flowEbosOilWaterSetDeck(setupTime_, std::move(deck_), std::move(eclipseState_), std::move(schedule_), std::move(summaryConfig_));
                    return flowEbosOilWaterMain(argc_, argv_, outputCout_, outputFiles_);
                }
                // gas-water
                else if ( phases.active( Phase::GAS ) && phases.active( Phase::WATER ) ) {
                    flowEbosGasWaterSetDeck(setupTime_, std::move(deck_), std::move(eclipseState_), std::move(schedule_), std::move(summaryConfig_));
                    return flowEbosGasWaterMain(argc_, argv_, outputCout_, outputFiles_);
                }
                else {
                    if (outputCout_)
                        std::cerr << "No suitable configuration found, valid are Twophase (oilwater, oilgas and gaswater), polymer, solvent, or blackoil" << std::endl;
                    return EXIT_FAILURE;
                }
            }
            // Polymer case
            else if ( phases.active( Phase::POLYMER ) ) {
                if ( !phases.active( Phase::WATER) ) {
                    if (outputCout_)
                        std::cerr << "No valid configuration is found for polymer simulation, valid options include "
                                  << "oilwater + polymer and blackoil + polymer" << std::endl;
                    return EXIT_FAILURE;
                }

                // Need to track the polymer molecular weight
                // for the injectivity study
                if ( phases.active( Phase::POLYMW ) ) {
                    // only oil water two phase for now
                    assert( phases.size() == 4);
                    return flowEbosOilWaterPolymerInjectivityMain(argc_, argv_, outputCout_, outputFiles_);
                }

                if ( phases.size() == 3 ) { // oil water polymer case
                    flowEbosOilWaterPolymerSetDeck(setupTime_, std::move(deck_),
                                                   std::move(eclipseState_),
                                                   std::move(schedule_),
                                                   std::move(summaryConfig_));
                    return flowEbosOilWaterPolymerMain(argc_, argv_, outputCout_, outputFiles_);
                } else {
                    flowEbosPolymerSetDeck(setupTime_, std::move(deck_),
                                                std::move(eclipseState_),
                                                std::move(schedule_),
                                                std::move(summaryConfig_));
                    return flowEbosPolymerMain(argc_, argv_, outputCout_, outputFiles_);
                }
            }
            // Foam case
            else if ( phases.active( Phase::FOAM ) ) {
                flowEbosFoamSetDeck(setupTime_, std::move(deck_),
                                    std::move(eclipseState_),
                                    std::move(schedule_),
                                    std::move(summaryConfig_));
                return flowEbosFoamMain(argc_, argv_, outputCout_, outputFiles_);
            }
            // Brine case
            else if ( phases.active( Phase::BRINE ) ) {
                if ( !phases.active( Phase::WATER) ) {
                    if (outputCout_)
                        std::cerr << "No valid configuration is found for brine simulation, valid options include "
                                  << "oilwater + brine and blackoil + brine" << std::endl;
                    return EXIT_FAILURE;
                }
                if ( phases.size() == 3 ) { // oil water brine case
                    flowEbosOilWaterBrineSetDeck(setupTime_, std::move(deck_),
                                                 std::move(eclipseState_),
                                                 std::move(schedule_),
                                                 std::move(summaryConfig_));
                    return flowEbosOilWaterBrineMain(argc_, argv_, outputCout_, outputFiles_);
                } else {
                    flowEbosBrineSetDeck(setupTime_, std::move(deck_),
                                         std::move(eclipseState_),
                                         std::move(schedule_),
                                         std::move(summaryConfig_));
                    return flowEbosBrineMain(argc_, argv_, outputCout_, outputFiles_);
                }
            }
            // Solvent case
            else if ( phases.active( Phase::SOLVENT ) ) {
                flowEbosSolventSetDeck(setupTime_, std::move(deck_),
                                       std::move(eclipseState_),
                                       std::move(schedule_),
                                       std::move(summaryConfig_));
                return flowEbosSolventMain(argc_, argv_, outputCout_, outputFiles_);
            }
            // Extended BO case
            else if ( phases.active( Phase::ZFRACTION ) ) {
                flowEbosExtboSetDeck(setupTime_, std::move(deck_),
                                     std::move(eclipseState_),
                                     std::move(schedule_),
                                     std::move(summaryConfig_));
                return flowEbosExtboMain(argc_, argv_, outputCout_, outputFiles_);
            }
            // Energy case
            else if (eclipseState_->getSimulationConfig().isThermal()) {
                flowEbosEnergySetDeck(setupTime_, std::move(deck_),
                                      std::move(eclipseState_),
                                      std::move(schedule_),
                                      std::move(summaryConfig_));
                return flowEbosEnergyMain(argc_, argv_, outputCout_, outputFiles_);
            }
#endif // FLOW_BLACKOIL_ONLY
            // Blackoil case
            else if( phases.size() == 3 ) {
                flowEbosBlackoilSetDeck(setupTime_, std::move(deck_),
                                        std::move(eclipseState_),
                                        std::move(schedule_),
                                        std::move(udqState_),
                                        std::move(summaryConfig_));
                return flowEbosBlackoilMain(argc_, argv_, outputCout_, outputFiles_);
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
            flowEbosSetDeck<TypeTag>(std::move(deck_),
                                     std::move(eclipseState_),
                                     std::move(schedule_),
                                     std::move(summaryConfig_));
            return flowEbosMain<TypeTag>(argc_, argv_, outputCout_, outputFiles_);
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
#if HAVE_DUNE_FEM
            int mpiRank = Dune::Fem::MPIManager::rank();
#else
            int mpiRank = EclGenericVanguard::comm().rank();
#endif

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
            int status = FlowMainEbos<PreTypeTag>::setupParameters_(argc_, argv_);
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
            if (outputCout_) {
                FlowMainEbos<PreTypeTag>::printBanner();
            }
            // Create Deck and EclipseState.
            try {
                auto python = std::make_shared<Python>();
                const bool init_from_restart_file = !EWOMS_GET_PARAM(PreTypeTag, bool, SchedRestart);
                if (outputDir.empty())
                    outputDir = EWOMS_GET_PARAM(PreTypeTag, std::string, OutputDir);
                outputMode = setupLogging(mpiRank,
                                          deckFilename,
                                          outputDir,
                                          EWOMS_GET_PARAM(PreTypeTag, std::string, OutputMode),
                                          outputCout_, "STDOUT_LOGGER");
                auto parseContext =
                    std::make_unique<ParseContext>(std::vector<std::pair<std::string , InputError::Action>>
                                                   {{ParseContext::PARSE_RANDOM_SLASH, InputError::IGNORE},
                                                    {ParseContext::PARSE_MISSING_DIMS_KEYWORD, InputError::WARN},
                                                    {ParseContext::SUMMARY_UNKNOWN_WELL, InputError::WARN},
                                                    {ParseContext::SUMMARY_UNKNOWN_GROUP, InputError::WARN}});
                if (EWOMS_GET_PARAM(PreTypeTag, bool, EclStrictParsing))
                    parseContext->update(InputError::DELAYED_EXIT1);

                FlowMainEbos<PreTypeTag>::printPRTHeader(outputCout_);

                if (outputCout_) {
                    OpmLog::info("Reading deck file '" + deckFilename + "'");
                }

                std::optional<int> outputInterval;
                int output_param = EWOMS_GET_PARAM(PreTypeTag, int, EclOutputInterval);
                if (output_param >= 0)
                    outputInterval = output_param;

                readDeck(mpiRank, deckFilename, deck_, eclipseState_, schedule_, udqState_, actionState_,
                         summaryConfig_, nullptr, python, std::move(parseContext),
                         init_from_restart_file, outputCout_, outputInterval);

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

        filesystem::path simulationCaseName_( const std::string& casename ) {
            namespace fs = ::Opm::filesystem;

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
                    std::cout << "flow " << moduleVersionName() << std::endl;
                    std::exit(EXIT_SUCCESS);
                }
            }
        }


        int argc_;
        char** argv_;
        bool outputCout_;
        bool outputFiles_;
        double setupTime_;
        std::string deckFilename_;
        std::string flowProgName_;
        char *saveArgs_[2];
        std::unique_ptr<Deck> deck_;
        std::unique_ptr<EclipseState> eclipseState_;
        std::unique_ptr<Schedule> schedule_;
        std::unique_ptr<UDQState> udqState_;
        std::unique_ptr<Action::State> actionState_;
        std::unique_ptr<SummaryConfig> summaryConfig_;
    };

} // namespace Opm

#endif // OPM_MAIN_HEADER_INCLUDED
