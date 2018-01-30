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

#include <opm/simulators/flow_ebos_blackoil.hpp>
#include <opm/simulators/flow_ebos_gasoil.hpp>
#include <opm/simulators/flow_ebos_oilwater.hpp>
#include <opm/simulators/flow_ebos_solvent.hpp>
#include <opm/simulators/flow_ebos_polymer.hpp>

#include <opm/autodiff/MissingFeatures.hpp>
#include <opm/common/utility/parameters/ParameterGroup.hpp>
#include <opm/common/utility/ResetLocale.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/checkDeck.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

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
}


// ----------------- Main program -----------------
int main(int argc, char** argv)
{
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

    const bool outputCout = (mpiRank == 0);

    // we always want to use the default locale, and thus spare us the trouble
    // with incorrect locale settings.
    Opm::resetLocale();

    Opm::ParameterGroup param(argc, argv, false, outputCout);

    // See if a deck was specified on the command line.
    if (!param.unhandledArguments().empty()) {
        if (param.unhandledArguments().size() != 1) {
            if (outputCout)
                std::cerr << "You can only specify a single input deck on the command line.\n";
            return EXIT_FAILURE;
        } else {
            const auto casename = detail::simulationCaseName( param.unhandledArguments()[ 0 ] );
            param.insertParameter("deck_filename", casename.string() );
        }
    }

    // We must have an input deck. Grid and props will be read from that.
    if (!param.has("deck_filename")) {
        if (outputCout)
            std::cerr << "This program must be run with an input deck.\n"
                "Specify the deck filename either\n"
                "    a) as a command line argument by itself\n"
                "    b) as a command line parameter with the syntax deck_filename=<path to your deck>, or\n"
                "    c) as a parameter in a parameter file (.param or .xml) passed to the program.\n";
        return EXIT_FAILURE;
    }

    std::string deckFilename = param.get<std::string>("deck_filename");

    // Create Deck and EclipseState.
    try {
        Opm::Parser parser;
        typedef std::pair<std::string, Opm::InputError::Action> ParseModePair;
        typedef std::vector<ParseModePair> ParseModePairs;
        ParseModePairs tmp;
        tmp.push_back(ParseModePair(Opm::ParseContext::PARSE_RANDOM_SLASH, Opm::InputError::IGNORE));
        tmp.push_back(ParseModePair(Opm::ParseContext::PARSE_MISSING_DIMS_KEYWORD, Opm::InputError::WARN));
        tmp.push_back(ParseModePair(Opm::ParseContext::SUMMARY_UNKNOWN_WELL, Opm::InputError::WARN));
        tmp.push_back(ParseModePair(Opm::ParseContext::SUMMARY_UNKNOWN_GROUP, Opm::InputError::WARN));
        Opm::ParseContext parseContext(tmp);

        std::shared_ptr<Opm::Deck> deck = std::make_shared< Opm::Deck >( parser.parseFile(deckFilename , parseContext) );
        if ( outputCout ) {
            Opm::checkDeck(*deck, parser);
            Opm::MissingFeatures::checkKeywords(*deck);
        }
        Opm::Runspec runspec( *deck );
        const auto& phases = runspec.phases();

        std::shared_ptr<Opm::EclipseState> eclipseState = std::make_shared< Opm::EclipseState > ( *deck, parseContext );
        std::shared_ptr<Opm::Schedule> schedule = std::make_shared<Opm::Schedule>(*deck, eclipseState->getInputGrid(), eclipseState->get3DProperties(), phases, parseContext);
        std::shared_ptr<Opm::SummaryConfig> summary_config = std::make_shared<Opm::SummaryConfig>(*deck, *schedule, eclipseState->getTableManager(), parseContext);
                // Twophase cases
        if( phases.size() == 2 ) {
            // oil-gas
            if (phases.active( Opm::Phase::GAS ))
            {
                Opm::flowEbosGasOilSetDeck(*deck, *eclipseState, *schedule, *summary_config);
                return Opm::flowEbosGasOilMain(argc, argv);
            }
            // oil-water
            else if ( phases.active( Opm::Phase::WATER ) )
            {
                Opm::flowEbosOilWaterSetDeck(*deck, *eclipseState, *schedule, *summary_config);
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
            Opm::flowEbosPolymerSetDeck(*deck, *eclipseState, *schedule, *summary_config);
            return Opm::flowEbosPolymerMain(argc, argv);
        }
        // Solvent case
        else if ( phases.active( Opm::Phase::SOLVENT ) ) {
            Opm::flowEbosSolventSetDeck(*deck, *eclipseState, *schedule, *summary_config);
            return Opm::flowEbosSolventMain(argc, argv);
        }
        // Blackoil case
        else if( phases.size() == 3 ) {
            Opm::flowEbosBlackoilSetDeck(*deck, *eclipseState, *schedule, *summary_config);
            return Opm::flowEbosBlackoilMain(argc, argv);
        }
        else
        {
            if (outputCout)
                std::cerr << "No suitable configuration found, valid are Twophase, polymer, solvent, or blackoil" << std::endl;
            return EXIT_FAILURE;
        }
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
