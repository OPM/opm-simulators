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
#include "config.h"

#include <flow/flow_ebos_blackoil.hpp>
#include <flow/flow_ebos_gasoil.hpp>
#include <flow/flow_ebos_oilwater.hpp>
#include <flow/flow_ebos_solvent.hpp>
#include <flow/flow_ebos_polymer.hpp>
#include <flow/flow_ebos_energy.hpp>
#include <flow/flow_ebos_oilwater_polymer.hpp>
#include <flow/flow_ebos_oilwater_polymer_injectivity.hpp>

#include <opm/simulators/flow/SimulatorFullyImplicitBlackoilEbos.hpp>
#include <opm/simulators/flow/FlowMainEbos.hpp>
#include <opm/simulators/utils/moduleVersion.hpp>
#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>
#include <opm/simulators/flow/MissingFeatures.hpp>
#include <opm/material/common/ResetLocale.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/Parser/ErrorGuard.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/checkDeck.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>


#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

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

// ----------------- Main program -----------------
int main(int argc, char** argv)
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
    if (status != 0)
        // if setupParameters_ returns a value smaller than 0, there was no error, but
        // the program should abort. This is the case e.g. for the --help and the
        // --print-properties parameters.
        return (status >= 0)?status:0;

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
            std::cerr << "Exception received: " << e.what() << ". Try '--help' for a usage description.\n";
#if HAVE_MPI
        MPI_Finalize();
#endif
        return 1;
    }

    if (outputCout) {
        Opm::FlowMainEbos<PreTypeTag>::printBanner();
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

            if (EWOMS_GET_PARAM(PreTypeTag, bool, EclStrictParsing))
                parseContext.update( Opm::InputError::DELAYED_EXIT1);
            else {
                parseContext.update(Opm::ParseContext::PARSE_RANDOM_SLASH, Opm::InputError::IGNORE);
                parseContext.update(Opm::ParseContext::PARSE_MISSING_DIMS_KEYWORD, Opm::InputError::WARN);
                parseContext.update(Opm::ParseContext::SUMMARY_UNKNOWN_WELL, Opm::InputError::WARN);
                parseContext.update(Opm::ParseContext::SUMMARY_UNKNOWN_GROUP, Opm::InputError::WARN);
            }

            deck.reset( new Opm::Deck( parser.parseFile(deckFilename , parseContext, errorGuard)));
            Opm::MissingFeatures::checkKeywords(*deck, parseContext, errorGuard);

            if ( outputCout )
                Opm::checkDeck(*deck, parser, parseContext, errorGuard);

            eclipseState.reset( new Opm::EclipseState(*deck, parseContext, errorGuard ));
            schedule.reset(new Opm::Schedule(*deck, *eclipseState, parseContext, errorGuard));
            summaryConfig.reset( new Opm::SummaryConfig(*deck, *schedule, eclipseState->getTableManager(), parseContext, errorGuard));

            if (errorGuard) {
                errorGuard.dump();
                errorGuard.clear();

                throw std::runtime_error("Unrecoverable errors were encountered while loading input.");
            }
        }
        const auto& phases = Opm::Runspec(*deck).phases();

        // run the actual simulator
        //
        // TODO: make sure that no illegal combinations like thermal and twophase are
        //       requested.

        // Twophase cases
        if( phases.size() == 2 ) {
            // oil-gas
            if (phases.active( Opm::Phase::GAS ))
            {
                Opm::flowEbosGasOilSetDeck(externalSetupTimer.elapsed(), *deck, *eclipseState, *schedule, *summaryConfig);
                return Opm::flowEbosGasOilMain(argc, argv);
            }
            // oil-water
            else if ( phases.active( Opm::Phase::WATER ) )
            {
                Opm::flowEbosOilWaterSetDeck(externalSetupTimer.elapsed(), *deck, *eclipseState, *schedule, *summaryConfig);
                return Opm::flowEbosOilWaterMain(argc, argv);
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
                return Opm::flowEbosOilWaterPolymerInjectivityMain(argc, argv);
            }

            if ( phases.size() == 3 ) { // oil water polymer case
                Opm::flowEbosOilWaterPolymerSetDeck(externalSetupTimer.elapsed(), *deck, *eclipseState, *schedule, *summaryConfig);
                return Opm::flowEbosOilWaterPolymerMain(argc, argv);
            } else {
                Opm::flowEbosPolymerSetDeck(externalSetupTimer.elapsed(), *deck, *eclipseState, *schedule, *summaryConfig);
                return Opm::flowEbosPolymerMain(argc, argv);
            }
        }
        // Solvent case
        else if ( phases.active( Opm::Phase::SOLVENT ) ) {
            Opm::flowEbosSolventSetDeck(externalSetupTimer.elapsed(), *deck, *eclipseState, *schedule, *summaryConfig);
            return Opm::flowEbosSolventMain(argc, argv);
        }
        // Energy case
        else if (eclipseState->getSimulationConfig().isThermal()) {
            Opm::flowEbosEnergySetDeck(externalSetupTimer.elapsed(), *deck, *eclipseState, *schedule, *summaryConfig);
            return Opm::flowEbosEnergyMain(argc, argv);
        }
        // Blackoil case
        else if( phases.size() == 3 ) {
            Opm::flowEbosBlackoilSetDeck(externalSetupTimer.elapsed(), *deck, *eclipseState, *schedule, *summaryConfig);
            return Opm::flowEbosBlackoilMain(argc, argv);
        }
        else
        {
            if (outputCout)
                std::cerr << "No suitable configuration found, valid are Twophase, polymer, solvent, energy, or blackoil" << std::endl;
            return EXIT_FAILURE;
        }
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
