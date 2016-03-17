/*
  Copyright 2015 Statoil ASA.

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

/* --- Boost.Test boilerplate --- */
#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif

#define NVERBOSE  // Suppress own messages when throw()ing

#define BOOST_TEST_MODULE RelpermDiagnostics


#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <opm/core/grid.h>
#include <opm/core/grid/cart_grid.h>
#include <opm/core/grid/GridManager.hpp>

#include <opm/core/props/satfunc/RelpermDiagnostics.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

BOOST_AUTO_TEST_SUITE ()

BOOST_AUTO_TEST_CASE(diagnosis)
{
    using namespace Opm;
    EclipseStateConstPtr eclState; 
    ParserPtr parser(new Opm::Parser);
    Opm::ParseContext parseContext({{ ParseContext::PARSE_RANDOM_SLASH , InputError::IGNORE }, 
                              { ParseContext::PARSE_UNKNOWN_KEYWORD, InputError::IGNORE},
                              { ParseContext::PARSE_RANDOM_TEXT, InputError::IGNORE}
                             });
    Opm::DeckConstPtr deck(parser->parseFile("../tests/relpermDiagnostics.DATA", parseContext));
    eclState.reset(new EclipseState(deck, parseContext));

    GridManager gm(deck);
    const UnstructuredGrid& grid = *gm.c_grid();
    std::string logFile = "LOGFILE.txt";
    RelpermDiagnostics diagnostics(logFile);
    diagnostics.diagnosis(eclState, deck, grid);
    auto msg = diagnostics.getMessages();
    BOOST_CHECK(!msg.empty());
    BOOST_CHECK_EQUAL(msg.size(), 1);
}
BOOST_AUTO_TEST_SUITE_END()
