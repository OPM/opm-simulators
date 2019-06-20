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
#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>
#include <opm/common/utility/parameters/ParameterGroup.hpp>
#include <opm/material/common/ResetLocale.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/checkDeck.hpp>

//#include <opm/material/fluidsystems/BlackOilFluidSystemSimple.hpp>
//#include <opm/material/fluidsystems/BlackOilFluidSystemSimple.hpp> 
#include <ewoms/models/blackoil/blackoilintensivequantities.hh>
#include <opm/material/fluidstates/BlackOilFluidState.hpp>
//#include <opm/material/fluidstates/BlackOilFluidStateSimple.hpp>

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




namespace Opm {
  template <class TypeTag>
  void flowEbosSetDeck(Deck &deck, EclipseState& eclState)
  {
    typedef typename GET_PROP_TYPE(TypeTag, Vanguard) Vanguard;
    Vanguard::setExternalDeck(&deck);
    Vanguard::setExternalEclState(&eclState);
  }
  
// ----------------- Main program -----------------
  template <class TypeTag>
  int flowEbosMain(int argc, char** argv)
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
    return mainfunc.execute(argc, argv);
  }

}





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
template<class TypeTag>
int mainFlow(int argc, char** argv)
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

    PreProblem::setBriefDescription("Simple Flow, an advanced reservoir simulator for ECL-decks provided by the Open Porous Media project.");
    
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
        Ewoms::Parameters::printUsage<PreTypeTag>(PreProblem::helpPreamble(argc, const_cast<const char**>(argv)),
                                                  e.what());
        return 1;
    }

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
      Opm::ErrorGuard errorGuard;

      std::shared_ptr<Opm::Deck> deck = std::make_shared< Opm::Deck >( parser.parseFile(deckFilename , parseContext, errorGuard) );
      if ( outputCout ) {
          Opm::checkDeck(*deck, parser, parseContext, errorGuard);
	Opm::MissingFeatures::checkKeywords(*deck);
      }
      //Opm::Runspec runspec( *deck );
      //const auto& phases = runspec.phases();

      std::shared_ptr<Opm::EclipseState> eclipseState = std::make_shared< Opm::EclipseState > ( *deck, parseContext, errorGuard );
      Opm::flowEbosSetDeck<TypeTag>(*deck, *eclipseState);
      return Opm::flowEbosMain<TypeTag>(argc, argv);
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
