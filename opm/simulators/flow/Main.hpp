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
#  include <flow/flow_ebos_micp.hpp>
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
#include <opm/parser/eclipse/EclipseState/Schedule/Well/WellTestState.hpp>

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

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

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
    void flowEbosSetDeck(std::shared_ptr<Deck> deck,
                         std::shared_ptr<EclipseState> eclState,
                         std::shared_ptr<Schedule> schedule,
                         std::shared_ptr<SummaryConfig> summaryConfig)
    {
        using Vanguard = GetPropType<TypeTag, Properties::Vanguard>;

        Vanguard::setExternalDeck(deck);
        Vanguard::setExternalEclState(eclState);
        Vanguard::setExternalSchedule(schedule);
        Vanguard::setExternalSummaryConfig(summaryConfig);
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

        // This constructor can be called from Python
        Main(const std::string& filename)
        {
            setArgvArgc_(filename);
            initMPI();
        }

        // This constructor can be called from Python when Python has
        // already parsed a deck
        Main(std::shared_ptr<Deck> deck,
             std::shared_ptr<EclipseState> eclipseState,
             std::shared_ptr<Schedule> schedule,
             std::shared_ptr<SummaryConfig> summaryConfig)
            : deck_{std::move(deck)}
            , eclipseState_{std::move(eclipseState)}
            , schedule_{std::move(schedule)}
            , summaryConfig_{std::move(summaryConfig)}
        {
            setArgvArgc_(deck_->getDataFile());
            initMPI();
        }

#define DEMONSTRATE_RUN_WITH_NONWORLD_COMM 0

        ~Main()
        {
#if DEMONSTRATE_RUN_WITH_NONWORLD_COMM
#if HAVE_MPI
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
#endif // HAVE_MPI
#endif // DEMONSTRATE_RUN_WITH_NONWORLD_COMM

            EclGenericVanguard::setCommunication(nullptr);

#if HAVE_MPI && !HAVE_DUNE_FEM
            MPI_Finalize();
#endif
        }

        void setArgvArgc_(const std::string& filename)
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

        void initMPI()
        {
#if HAVE_DUNE_FEM
            Dune::Fem::MPIManager::initialize(argc_, argv_);
#elif HAVE_MPI
            MPI_Init(&argc_, &argv_);
#endif
            EclGenericVanguard::setCommunication(std::make_unique<Parallel::Communication>());

#if DEMONSTRATE_RUN_WITH_NONWORLD_COMM
#if HAVE_MPI
            if (EclGenericVanguard::comm().size() > 1) {
                int world_rank = EclGenericVanguard::comm().rank();
                int color = (world_rank == 0);
                MPI_Comm new_comm;
                MPI_Comm_split(EclGenericVanguard::comm(), color, world_rank, &new_comm);
                isSimulationRank_ = (world_rank > 0);
                EclGenericVanguard::setCommunication(std::make_unique<Parallel::Communication>(new_comm));
            }
#endif // HAVE_MPI
#endif // DEMONSTRATE_RUN_WITH_NONWORLD_COMM
        }

        int runDynamic()
        {
            int exitCode = EXIT_SUCCESS;
            if (isSimulationRank_) {
                if (initialize_<Properties::TTag::FlowEarlyBird>(exitCode)) {
                    return this->dispatchDynamic_();
                }
            }

            return exitCode;
        }

        template <class TypeTag>
        int runStatic()
        {
            int exitCode = EXIT_SUCCESS;
            if (isSimulationRank_) {
                if (initialize_<TypeTag>(exitCode)) {
                    return this->dispatchStatic_<TypeTag>();
                }
            }

            return exitCode;
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
                    deck_,
                    eclipseState_,
                    schedule_,
                    std::move(udqState_),
                    std::move(this->actionState_),
                    std::move(this->wtestState_),
                    summaryConfig_);
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
            // TODO: make sure that no illegal combinations like thermal and
            //       twophase are requested.

            if (false) {}

#ifndef FLOW_BLACKOIL_ONLY
            // Single-phase case
            else if( eclipseState_->runspec().micp() ) {
              // micp
              if ( !phases.active( Phase::WATER) || phases.size() > 2) {
                if (outputCout_)
                std::cerr << "No valid configuration is found for MICP simulation, the only valid option is "
                      << "water + MICP" << std::endl;
                return EXIT_FAILURE;
            }
            flowEbosMICPSetDeck(
                        setupTime_, deck_, eclipseState_, schedule_, summaryConfig_);
                    return flowEbosMICPMain(argc_, argv_, outputCout_, outputFiles_);
            }
            // Twophase cases
            else if (phases.size() == 2) {
                return this->runTwoPhase(phases);
            }

            // Polymer case
            else if (phases.active(Phase::POLYMER)) {
                return this->runPolymer(phases);
            }

            // Foam case
            else if (phases.active(Phase::FOAM)) {
                return this->runFoam();
            }

            // Brine case
            else if (phases.active(Phase::BRINE)) {
                return this->runBrine(phases);
            }

            // Solvent case
            else if (phases.active(Phase::SOLVENT)) {
                return this->runSolvent();
            }

            // Extended BO case
            else if (phases.active(Phase::ZFRACTION)) {
                return this->runExtendedBlackOil();
            }

            // Energy case
            else if (eclipseState_->getSimulationConfig().isThermal()) {
                return this->runThermal();
            }
#endif // FLOW_BLACKOIL_ONLY

            // Blackoil case
            else if (phases.size() == 3) {
                return this->runBlackOil();
            }

            else {
                if (outputCout_) {
                    std::cerr << "No suitable configuration found, valid are "
                              << "Twophase, polymer, foam, brine, solvent, "
                              << "energy, and blackoil.\n";
                }

                return EXIT_FAILURE;
            }
        }

        template <class TypeTag>
        int dispatchStatic_()
        {
            flowEbosSetDeck<TypeTag>(
                deck_, eclipseState_, schedule_, summaryConfig_);
            return flowEbosMain<TypeTag>(argc_, argv_, outputCout_, outputFiles_);
        }

        /// \brief Initialize
        /// \param exitCode The exitCode of the program.
        ///
        /// \return Whether to actually run the simulator. I.e. true if
        /// parsing of command line was successful and no --help,
        /// --print-properties, or --print-parameters have been found.
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
            int status = FlowMainEbos<PreTypeTag>::setupParameters_(argc_, argv_, EclGenericVanguard::comm());
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
                FlowMainEbos<PreTypeTag>::printBanner(EclGenericVanguard::comm());
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

                readDeck(EclGenericVanguard::comm(), deckFilename, deck_, eclipseState_, schedule_, udqState_, actionState_, wtestState_,
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

        filesystem::path simulationCaseName_(const std::string& casename)
        {
            namespace fs = ::Opm::filesystem;

            auto exists = [](const fs::path& f)
            {
                return (fs::exists(f) && fs::is_regular_file(f))
                    || (fs::is_symlink(f) &&
                        fs::is_regular_file(fs::read_symlink(f)));
            };

            auto simcase = fs::path { casename };

            if (exists(simcase)) {
                return simcase;
            }

            for (const auto& ext : { std::string("DATA"), std::string("data") }) {
                if (exists(simcase.replace_extension(ext))) {
                    return simcase;
                }
            }

            throw std::invalid_argument {
                "Cannot find input case '" + casename + '\''
            };
        }

        // This function is an extreme special case, if the program has been invoked
        // *exactly* as:
        //
        //    flow   --version
        //
        // the call is intercepted by this function which will print "flow $version"
        // on stdout and exit(0).
        void handleVersionCmdLine_(int argc, char** argv)
        {
            auto pos = std::find_if(argv, argv + argc,
                [](const char* arg)
            {
                return std::strcmp(arg, "--version") == 0;
            });

            if (pos != argv + argc) {
                std::cout << "flow " << moduleVersionName() << std::endl;
                std::exit(EXIT_SUCCESS);
            }
        }

#ifndef FLOW_BLACKOIL_ONLY
        int runTwoPhase(const Phases& phases)
        {
            // oil-gas
            if (phases.active( Phase::OIL ) && phases.active( Phase::GAS )) {
                flowEbosGasOilSetDeck(
                    setupTime_, deck_, eclipseState_, schedule_, summaryConfig_);
                return flowEbosGasOilMain(argc_, argv_, outputCout_, outputFiles_);
            }

            // oil-water
            else if ( phases.active( Phase::OIL ) && phases.active( Phase::WATER ) ) {
                flowEbosOilWaterSetDeck(
                    setupTime_, deck_, eclipseState_, schedule_, summaryConfig_);
                return flowEbosOilWaterMain(argc_, argv_, outputCout_, outputFiles_);
            }

            // gas-water
            else if ( phases.active( Phase::GAS ) && phases.active( Phase::WATER ) ) {
                flowEbosGasWaterSetDeck(
                    setupTime_, deck_, eclipseState_, schedule_, summaryConfig_);
                return flowEbosGasWaterMain(argc_, argv_, outputCout_, outputFiles_);
            }
            else {
                if (outputCout_) {
                    std::cerr << "No suitable configuration found, valid are Twophase (oilwater, oilgas and gaswater), polymer, solvent, or blackoil" << std::endl;
                }

                return EXIT_FAILURE;
            }
        }

        int runPolymer(const Phases& phases)
        {
            if (! phases.active(Phase::WATER)) {
                if (outputCout_)
                    std::cerr << "No valid configuration is found for polymer simulation, valid options include "
                              << "oilwater + polymer and blackoil + polymer" << std::endl;

                return EXIT_FAILURE;
            }

            // Need to track the polymer molecular weight
            // for the injectivity study
            if (phases.active(Phase::POLYMW)) {
                // only oil water two phase for now
                assert (phases.size() == 4);
                return flowEbosOilWaterPolymerInjectivityMain(argc_, argv_, outputCout_, outputFiles_);
            }

            if (phases.size() == 3) { // oil water polymer case
                flowEbosOilWaterPolymerSetDeck(
                    setupTime_, deck_, eclipseState_, schedule_, summaryConfig_);
                return flowEbosOilWaterPolymerMain(argc_, argv_, outputCout_, outputFiles_);
            }
            else {
                flowEbosPolymerSetDeck(
                    setupTime_, deck_, eclipseState_, schedule_, summaryConfig_);
                return flowEbosPolymerMain(argc_, argv_, outputCout_, outputFiles_);
            }
        }

        int runFoam()
        {
            flowEbosFoamSetDeck(
                setupTime_, deck_, eclipseState_, schedule_, summaryConfig_);

            return flowEbosFoamMain(argc_, argv_, outputCout_, outputFiles_);
        }

        int runBrine(const Phases& phases)
        {
            if (! phases.active(Phase::WATER)) {
                if (outputCout_)
                    std::cerr << "No valid configuration is found for brine simulation, valid options include "
                              << "oilwater + brine and blackoil + brine" << std::endl;

                return EXIT_FAILURE;
            }

            if (phases.size() == 3) { // oil water brine case
                flowEbosOilWaterBrineSetDeck(
                    setupTime_, deck_, eclipseState_, schedule_, summaryConfig_);
                return flowEbosOilWaterBrineMain(argc_, argv_, outputCout_, outputFiles_);
            }
            else {
                flowEbosBrineSetDeck(
                    setupTime_, deck_, eclipseState_, schedule_, summaryConfig_);
                return flowEbosBrineMain(argc_, argv_, outputCout_, outputFiles_);
            }
        }

        int runSolvent()
        {
            flowEbosSolventSetDeck(
                setupTime_, deck_, eclipseState_, schedule_, summaryConfig_);

            return flowEbosSolventMain(argc_, argv_, outputCout_, outputFiles_);
        }

        int runExtendedBlackOil()
        {
            flowEbosExtboSetDeck(
                setupTime_, deck_, eclipseState_, schedule_, summaryConfig_);

            return flowEbosExtboMain(argc_, argv_, outputCout_, outputFiles_);
        }

        int runThermal()
        {
            flowEbosEnergySetDeck(
                setupTime_, deck_, eclipseState_, schedule_, summaryConfig_);

            return flowEbosEnergyMain(argc_, argv_, outputCout_, outputFiles_);
        }
#endif  // FLOW_BLACKOIL_ONLY

        int runBlackOil()
        {
            flowEbosBlackoilSetDeck(this->setupTime_,
                                    this->deck_,
                                    this->eclipseState_,
                                    this->schedule_,
                                    std::move(this->udqState_),
                                    std::move(this->actionState_),
                                    std::move(this->wtestState_),
                                    this->summaryConfig_);

            return flowEbosBlackoilMain(argc_, argv_, outputCout_, outputFiles_);
        }

        int argc_{0};
        char** argv_{nullptr};
        bool outputCout_{false};
        bool outputFiles_{false};
        double setupTime_{0.0};
        std::string deckFilename_{};
        std::string flowProgName_{};
        char *saveArgs_[3]{nullptr};
        std::unique_ptr<UDQState> udqState_{};
        std::unique_ptr<Action::State> actionState_{};
        std::unique_ptr<WellTestState> wtestState_{};

        // These variables may be owned by both Python and the simulator
        std::shared_ptr<Deck> deck_{};
        std::shared_ptr<EclipseState> eclipseState_{};
        std::shared_ptr<Schedule> schedule_{};
        std::shared_ptr<SummaryConfig> summaryConfig_{};

        // To demonstrate run with non_world_comm
        bool isSimulationRank_ = true;
    };

} // namespace Opm

#endif // OPM_MAIN_HEADER_INCLUDED
