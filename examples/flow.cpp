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

#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <memory>

#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/IOConfig/IOConfig.hpp>
#include <opm/parser/eclipse/EclipseState/InitConfig/InitConfig.hpp>
#include <opm/parser/eclipse/EclipseState/checkDeck.hpp>

// Define making clear that the simulator supports AMG
#define FLOW_SUPPORT_AMG 1

#include <opm/material/densead/Evaluation.hpp>
#include <ewoms/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/autodiff/DuneMatrix.hpp>
#include <opm/autodiff/SimulatorFullyImplicitBlackoilEbos.hpp>
#include <opm/autodiff/FlowMainEbos.hpp>

namespace Ewoms {
namespace Properties {

    ///////////////////////////////////
    //   Twophase case
    ///////////////////////////////////

    NEW_TYPE_TAG(EclFlowOilWaterProblem, INHERITS_FROM(EclFlowProblem));
    //! The indices required by the model
    SET_TYPE_PROP(EclFlowOilWaterProblem, Indices,
      Ewoms::BlackOilTwoPhaseIndices<GET_PROP_VALUE(TypeTag, EnableSolvent)?1:0, GET_PROP_VALUE(TypeTag, EnablePolymer)?1:0,  /*PVOffset=*/0, /*disabledCompIdx=*/2>);


    NEW_TYPE_TAG(EclFlowGasOilProblem, INHERITS_FROM(EclFlowProblem));
    //! The indices required by the model
    SET_TYPE_PROP(EclFlowGasOilProblem, Indices,
      Ewoms::BlackOilTwoPhaseIndices<GET_PROP_VALUE(TypeTag, EnableSolvent)?1:0, GET_PROP_VALUE(TypeTag, EnablePolymer)?1:0,  /*PVOffset=*/0, /*disabledCompIdx=*/1>);

    ///////////////////////////////////
    //   Polymer case
    ///////////////////////////////////

    NEW_TYPE_TAG(EclFlowPolymerProblem, INHERITS_FROM(EclFlowProblem));
    SET_BOOL_PROP(EclFlowPolymerProblem, EnablePolymer, true);


    ///////////////////////////////////
    //   Solvent case
    ///////////////////////////////////

    NEW_TYPE_TAG(EclFlowSolventProblem, INHERITS_FROM(EclFlowProblem));
    SET_BOOL_PROP(EclFlowSolventProblem, EnableSolvent, true);

}} // end namespaces


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
    // Must ensure an instance of the helper is created to initialise MPI.
    // For a build without MPI the Dune::FakeMPIHelper is used, so rank will
    // be 0 and size 1.
    const Dune::MPIHelper& mpi_helper = Dune::MPIHelper::instance(argc, argv);
    const bool outputCout = mpi_helper.rank() == 0;

    Opm::ParameterGroup param(argc, argv, false, outputCout);

    // See if a deck was specified on the command line.
    if (!param.unhandledArguments().empty()) {
        if (param.unhandledArguments().size() != 1) {
            std::cerr << "You can only specify a single input deck on the command line.\n";
            return EXIT_FAILURE;
        } else {
            const auto casename = detail::simulationCaseName( param.unhandledArguments()[ 0 ] );
            param.insertParameter("deck_filename", casename.string() );
        }
    }

    // We must have an input deck. Grid and props will be read from that.
    if (!param.has("deck_filename")) {
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
        Opm::checkDeck(*deck, parser);
        if ( outputCout ) {
            Opm::MissingFeatures::checkKeywords(*deck);
        }

        std::shared_ptr<Opm::EclipseState> eclipseState =
            std::make_shared< Opm::EclipseState > ( *deck, parseContext );

        Opm::Runspec runspec( *deck );
        const auto& phases = runspec.phases();

        // Twophase case
        if( phases.size() == 2 ) {
            // oil-gas
            if (phases.active( Opm::Phase::GAS ))
            {
                Opm::FlowMainEbos<TTAG(EclFlowGasOilProblem)> mainfunc;
                return mainfunc.execute(argc, argv, deck, eclipseState );
            }
            // oil-water
            else if ( phases.active( Opm::Phase::WATER ) )
            {
                Opm::FlowMainEbos<TTAG(EclFlowOilWaterProblem)> mainfunc;
                return mainfunc.execute(argc, argv, deck, eclipseState );
            }
            else {
                std::cerr << "No suitable configuration found, valid are Twophase (oilwater and oilgas), polymer, solvent, or blackoil" << std::endl;
                return EXIT_FAILURE;
            }
        }
        // Polymer case
        else if ( phases.active( Opm::Phase::POLYMER ) ) {
            Opm::FlowMainEbos<TTAG(EclFlowPolymerProblem)> mainfunc;
            return mainfunc.execute(argc, argv, deck, eclipseState );

        }
        // Solvent case
        else if ( phases.active( Opm::Phase::SOLVENT ) ) {
            Opm::FlowMainEbos<TTAG(EclFlowSolventProblem)> mainfunc;
            return mainfunc.execute(argc, argv, deck, eclipseState );

        }
        // Blackoil case
        else if( phases.size() == 3 ) {
            Opm::FlowMainEbos<TTAG(EclFlowProblem)> mainfunc;
            return mainfunc.execute(argc, argv, deck, eclipseState );
        }
        else
        {
            std::cerr << "No suitable configuration found, valid are Twophase, polymer, solvent, or blackoil" << std::endl;
            return EXIT_FAILURE;
        }
    }
    catch (const std::invalid_argument& e)
    {
        std::cerr << "Failed to create valid EclipseState object." << std::endl;
        std::cerr << "Exception caught: " << e.what() << std::endl;
        throw;
    }

    return EXIT_SUCCESS;
}
