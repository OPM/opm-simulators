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

#include <flow/flow_ebos_blackoil.hpp>
#include <flow/flow_ebos_blackoil_legacyassembly.hpp>

#include <flow/flow_ebos_gasoil.hpp>
#include <flow/flow_ebos_gasoildiffuse.hpp>
#include <flow/flow_ebos_gasoil_energy.hpp>
#include <flow/flow_ebos_oilwater.hpp>
#include <flow/flow_ebos_gaswater.hpp>
#include <flow/flow_ebos_gaswater_solvent.hpp>
#include <flow/flow_ebos_solvent.hpp>
#include <flow/flow_ebos_solvent_foam.hpp>
#include <flow/flow_ebos_polymer.hpp>
#include <flow/flow_ebos_extbo.hpp>
#include <flow/flow_ebos_foam.hpp>
#include <flow/flow_ebos_brine.hpp>
#include <flow/flow_ebos_brine_saltprecipitation.hpp>
#include <flow/flow_ebos_gaswater_saltprec_vapwat.hpp>
#include <flow/flow_ebos_gaswater_saltprec_energy.hpp>
#include <flow/flow_ebos_brine_precsalt_vapwat.hpp>
#include <flow/flow_ebos_onephase.hpp>
#include <flow/flow_ebos_onephase_energy.hpp>
#include <flow/flow_ebos_oilwater_brine.hpp>
#include <flow/flow_ebos_gaswater_brine.hpp>
#include <flow/flow_ebos_gaswater_energy.hpp>
#include <flow/flow_ebos_gaswater_dissolution.hpp>
#include <flow/flow_ebos_gaswater_dissolution_diffuse.hpp>
#include <flow/flow_ebos_energy.hpp>
#include <flow/flow_ebos_oilwater_polymer.hpp>
#include <flow/flow_ebos_oilwater_polymer_injectivity.hpp>
#include <flow/flow_ebos_micp.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>

#include <opm/simulators/flow/Banners.hpp>
#include <opm/simulators/flow/FlowMainEbos.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

#if HAVE_MPI
#include <opm/simulators/utils/ParallelEclipseState.hpp>
#endif

#if HAVE_DAMARIS
#include <opm/simulators/utils/DamarisKeywords.hpp>
#endif

#include <cassert>
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
    using InheritsFrom = std::tuple<EclFlowProblem>;
};
}

} // namespace Opm::Properties

namespace Opm {

namespace Action { class State; }
class UDQState;
class WellTestState;

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
    Main(int argc, char** argv);

    // This constructor can be called from Python
    Main(const std::string& filename);

    // This constructor can be called from Python when Python has
    // already parsed a deck
    Main(const std::string& filename,
         std::shared_ptr<EclipseState> eclipseState,
         std::shared_ptr<Schedule> schedule,
         std::shared_ptr<SummaryConfig> summaryConfig);

    ~Main();

    void setArgvArgc_(const std::string& filename);

    void initMPI();

    int runDynamic()
    {
        int exitCode = EXIT_SUCCESS;
        if (initialize_<Properties::TTag::FlowEarlyBird>(exitCode)) {
            if (isSimulationRank_) {
                return this->dispatchDynamic_();
            }
        }

        return exitCode;
    }

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

    using FlowMainEbosType = FlowMainEbos<Properties::TTag::EclFlowProblemTPFA>;
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
            this->setupVanguard();
            return flowEbosBlackoilTpfaMainInit(
                argc_, argv_, outputCout_, outputFiles_);
        } else {
            //NOTE: exitCode was set by initialize_() above;
            return std::unique_ptr<FlowMainEbosType>(); // nullptr
        }
    }

private:
    int dispatchDynamic_()
    {
        const auto& rspec = this->eclipseState_->runspec();
        const auto& phases = rspec.phases();

        this->setupVanguard();

        // run the actual simulator
        //
        // TODO: make sure that no illegal combinations like thermal and
        //       twophase are requested.
        const bool thermal = eclipseState_->getSimulationConfig().isThermal();

        // Single-phase case
        if (rspec.micp()) {
            return this->runMICP(phases);
        }

        // water-only case
        else if (phases.size() == 1 && phases.active(Phase::WATER) && !thermal) {
            return this->runWaterOnly(phases);
        }

        // water-only case with energy
        else if (phases.size() == 2 && phases.active(Phase::WATER) && thermal) {
            return this->runWaterOnlyEnergy(phases);
        }

        // Twophase cases
        else if (phases.size() == 2 && !thermal) {
            return this->runTwoPhase(phases);
        }

        // Polymer case
        else if (phases.active(Phase::POLYMER)) {
            return this->runPolymer(phases);
        }

        // Foam case
        else if (phases.active(Phase::FOAM) && !phases.active(Phase::SOLVENT)) {
            return this->runFoam();
        }

        // Solvent case
        else if (phases.active(Phase::SOLVENT)) {
            return this->runSolvent(phases);
        }

        // Brine case
        else if (phases.active(Phase::BRINE) && !thermal) {
            return this->runBrine(phases);
        }

        // Extended BO case
        else if (phases.active(Phase::ZFRACTION)) {
            return this->runExtendedBlackOil();
        }

        // Energy case
        else if (thermal) {
            return this->runThermal(phases);
        }

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
        this->setupVanguard();
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

        std::string deckFilename;
        std::string outputDir;
        if ( eclipseState_ ) {
            deckFilename = eclipseState_->getIOConfig().fullBasePath();
            outputDir = eclipseState_->getIOConfig().getOutputDir();
        }
        else {
            deckFilename = EWOMS_GET_PARAM(PreTypeTag, std::string, EclDeckFileName);
            outputDir = EWOMS_GET_PARAM(PreTypeTag, std::string, OutputDir);
        }

#if HAVE_DAMARIS
        enableDamarisOutput_ = EWOMS_GET_PARAM(PreTypeTag, bool, EnableDamarisOutput);
        
        // Reset to false as we cannot use Damaris if there is only one rank.
        if ((enableDamarisOutput_ == true) && (EclGenericVanguard::comm().size() == 1)) {
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
        
        int mpiRank = EclGenericVanguard::comm().rank();
        outputCout_ = false;
        if (mpiRank == 0)
            outputCout_ = EWOMS_GET_PARAM(PreTypeTag, bool, EnableTerminalOutput);

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
            printFlowBanner(EclGenericVanguard::comm().size(),
                            getNumThreads<PreTypeTag>(),
                            Opm::moduleVersionName());
            std::ostringstream str;
            Parameters::printValues<PreTypeTag>(str);
            cmdline_params = str.str();
        }

        // Create Deck and EclipseState.
        try {
            this->readDeck(deckFilename,
                           outputDir,
                           EWOMS_GET_PARAM(PreTypeTag, std::string, OutputMode),
                           !EWOMS_GET_PARAM(PreTypeTag, bool, SchedRestart),
                           EWOMS_GET_PARAM(PreTypeTag, bool,  EnableLoggingFalloutWarning),
                           EWOMS_GET_PARAM(PreTypeTag, std::string, ParsingStrictness),
                           mpiRank,
                           getNumThreads<PreTypeTag>(),
                           EWOMS_GET_PARAM(PreTypeTag, int, EclOutputInterval),
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

        exitCode = EXIT_SUCCESS;
        return true;
    }

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

    int runMICP(const Phases& phases)
    {
        if (!phases.active(Phase::WATER) || (phases.size() > 2)) {
            if (outputCout_) {
                std::cerr << "No valid configuration is found for MICP simulation, "
                          << "the only valid option is water + MICP\n";
            }

            return EXIT_FAILURE;
        }

        return flowEbosMICPMain(this->argc_,
                                this->argv_,
                                this->outputCout_,
                                this->outputFiles_);
    }

    int runTwoPhase(const Phases& phases)
    {
        const bool diffusive = eclipseState_->getSimulationConfig().isDiffusive();
        const bool disgasw = eclipseState_->getSimulationConfig().hasDISGASW();
        const bool vapwat = eclipseState_->getSimulationConfig().hasVAPWAT();

        // oil-gas
        if (phases.active( Phase::OIL ) && phases.active( Phase::GAS )) {
            if (diffusive) {
                return flowEbosGasOilDiffuseMain(argc_, argv_, outputCout_, outputFiles_);
            } else {
                return flowEbosGasOilMain(argc_, argv_, outputCout_, outputFiles_);
            }
        }

        // oil-water
        else if ( phases.active( Phase::OIL ) && phases.active( Phase::WATER ) ) {
            if (diffusive) {
                if (outputCout_) {
                    std::cerr << "The DIFFUSE option is not available for the two-phase water/oil model." << std::endl;
                }
                return EXIT_FAILURE;
            }
            return flowEbosOilWaterMain(argc_, argv_, outputCout_, outputFiles_);
        }

        // gas-water
        else if ( phases.active( Phase::GAS ) && phases.active( Phase::WATER ) ) {
            if (disgasw || vapwat) {
                if (diffusive) {
                    return flowEbosGasWaterDissolutionDiffuseMain(argc_, argv_, outputCout_, outputFiles_);
                }
                return flowEbosGasWaterDissolutionMain(argc_, argv_, outputCout_, outputFiles_);
            }
            if (diffusive) {
                if (outputCout_) {
                    std::cerr << "The DIFFUSE option is not available for the two-phase gas/water model without disgasw or vapwat." << std::endl;
                }
                return EXIT_FAILURE;
            }

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
            return flowEbosOilWaterPolymerMain(argc_, argv_, outputCout_, outputFiles_);
        }
        else {
            return flowEbosPolymerMain(argc_, argv_, outputCout_, outputFiles_);
        }
    }

    int runFoam()
    {
        return flowEbosFoamMain(argc_, argv_, outputCout_, outputFiles_);
    }

    int runWaterOnly(const Phases& phases)
    {
        if (!phases.active(Phase::WATER) || phases.size() != 1) {
            if (outputCout_)
                std::cerr << "No valid configuration is found for water-only simulation, valid options include "
                          << "water, water + thermal" << std::endl;

            return EXIT_FAILURE;
        }

        return flowEbosWaterOnlyMain(argc_, argv_, outputCout_, outputFiles_);
    }

    int runWaterOnlyEnergy(const Phases& phases)
    {
        if (!phases.active(Phase::WATER) || phases.size() != 2) {
            if (outputCout_)
                std::cerr << "No valid configuration is found for water-only simulation, valid options include "
                          << "water, water + thermal" << std::endl;

            return EXIT_FAILURE;
        }

        return flowEbosWaterOnlyEnergyMain(argc_, argv_, outputCout_, outputFiles_);
    }

    int runBrine(const Phases& phases)
    {
        if (! phases.active(Phase::WATER) || phases.size() == 2) {
            if (outputCout_)
                std::cerr << "No valid configuration is found for brine simulation, valid options include "
                          << "oilwater + brine, gaswater + brine and blackoil + brine" << std::endl;

            return EXIT_FAILURE;
        }

        if (phases.size() == 3) {

            if (phases.active(Phase::OIL)){ // oil water brine case
                return flowEbosOilWaterBrineMain(argc_, argv_, outputCout_, outputFiles_);
            }
            if (phases.active(Phase::GAS)){ // gas water brine case
                if (eclipseState_->getSimulationConfig().hasPRECSALT() &&
                    eclipseState_->getSimulationConfig().hasVAPWAT()) {
                    //case with water vaporization into gas phase and salt precipitation
                    return flowEbosGasWaterSaltprecVapwatMain(argc_, argv_, outputCout_, outputFiles_);
                }
                else {
                    return flowEbosGasWaterBrineMain(argc_, argv_, outputCout_, outputFiles_);
                }
            }
        }
        else if (eclipseState_->getSimulationConfig().hasPRECSALT()) {
            if (eclipseState_->getSimulationConfig().hasVAPWAT()) {
                    //case with water vaporization into gas phase and salt precipitation
                    return flowEbosBrinePrecsaltVapwatMain(argc_, argv_, outputCout_, outputFiles_);
            }
            else {
                return flowEbosBrineSaltPrecipitationMain(argc_, argv_, outputCout_, outputFiles_);
            }
        }
        else {
            return flowEbosBrineMain(argc_, argv_, outputCout_, outputFiles_);
        }

        return EXIT_FAILURE;
    }

    int runSolvent(const Phases& phases)
    {
        if (phases.active(Phase::FOAM)) {
            return flowEbosSolventFoamMain(argc_, argv_, outputCout_, outputFiles_);
        }
        // solvent + gas + water
        if (!phases.active( Phase::OIL ) && phases.active( Phase::WATER ) && phases.active( Phase::GAS )) {
            return flowEbosGasWaterSolventMain(argc_, argv_, outputCout_, outputFiles_);
        }

        // solvent + gas + water + oil
        if (phases.active( Phase::OIL ) && phases.active( Phase::WATER ) && phases.active( Phase::GAS )) {
            return flowEbosSolventMain(argc_, argv_, outputCout_, outputFiles_);
        }

        if (outputCout_)
            std::cerr << "No valid configuration is found for solvent simulation, valid options include "
                      << "gas + water + solvent and gas + oil + water + solvent" << std::endl;

        return EXIT_FAILURE;
    }

    int runExtendedBlackOil()
    {
        return flowEbosExtboMain(argc_, argv_, outputCout_, outputFiles_);
    }

    int runThermal(const Phases& phases)
    {
        // oil-gas-thermal
        if (!phases.active( Phase::WATER ) && phases.active( Phase::OIL ) && phases.active( Phase::GAS )) {
            return flowEbosGasOilEnergyMain(argc_, argv_, outputCout_, outputFiles_);
        }

        // water-gas-thermal
        if (!phases.active( Phase::OIL ) && phases.active( Phase::WATER ) && phases.active( Phase::GAS )) {

            if (phases.active(Phase::BRINE)){
                return flowEbosGasWaterSaltprecEnergyMain(argc_, argv_, outputCout_, outputFiles_);
            }
            return flowEbosGasWaterEnergyMain(argc_, argv_, outputCout_, outputFiles_);
        }

        return flowEbosEnergyMain(argc_, argv_, outputCout_, outputFiles_);
    }

    int runBlackOil()
    {
        const bool diffusive = eclipseState_->getSimulationConfig().isDiffusive();
        if (diffusive) {
            // Use the traditional linearizer, as the TpfaLinearizer does not
            // support the diffusion module yet.
            return flowEbosBlackoilMain(argc_, argv_, outputCout_, outputFiles_);
        } else {
            return flowEbosBlackoilTpfaMain(argc_, argv_, outputCout_, outputFiles_);
        }
    }

    void readDeck(const std::string& deckFilename,
                  const std::string& outputDir,
                  const std::string& outputMode,
                  const bool init_from_restart_file,
                  const bool allRanksDbgPrtLog,
                  const std::string& parsingStrictness,
                  const int mpiRank,
                  const std::size_t numThreads,
                  const int output_param,
                  const std::string& parameters,
                  std::string_view moduleVersion,
                  std::string_view compileTimestamp);

    void setupVanguard();

    template<class TypeTag>
    static int getNumThreads()
    {

        int threads;

#ifdef _OPENMP
        // This function is called before the parallel OpenMP stuff gets initialized.
        // That initialization happends after the deck is read and we want this message.
        // Hence we duplicate the code of setupParallelism to get the number of threads.
        if (std::getenv("OMP_NUM_THREADS")) {
            threads =  omp_get_max_threads();
        }
        else {
            threads = 2;

            const int input_threads = EWOMS_GET_PARAM(TypeTag, int, ThreadsPerProcess);

            if (input_threads > 0)
                threads = input_threads;
        }
#else
        threads = 1;
#endif

        return threads;
    }

#if HAVE_DAMARIS
    void setupDamaris(const std::string& outputDir);
#endif

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
    std::shared_ptr<EclipseState> eclipseState_{};
    std::shared_ptr<Schedule> schedule_{};
    std::shared_ptr<SummaryConfig> summaryConfig_{};

    // To demonstrate run with non_world_comm
    bool test_split_comm_ = false;
    bool isSimulationRank_ = true;
#if HAVE_DAMARIS
    bool enableDamarisOutput_ = false;
#endif
};

} // namespace Opm

#endif // OPM_MAIN_HEADER_INCLUDED
