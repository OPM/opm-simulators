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

#define NVERBOSE  // Suppress own messages when throw()ing

#define BOOST_TEST_MODULE RelpermDiagnostics


#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif
#include <opm/common/utility/platform_dependent/reenable_warnings.h>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/OpmLog/CounterLog.hpp>

#include <opm/grid/UnstructuredGrid.h>
#include <opm/grid/cart_grid.h>
#include <opm/grid/GridManager.hpp>

#include <opm/core/props/satfunc/RelpermDiagnostics.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

BOOST_AUTO_TEST_SUITE ()

BOOST_AUTO_TEST_CASE(diagnosis)
{
    using namespace Opm;
    Parser parser;

    Opm::Deck deck = parser.parseFile("../tests/relpermDiagnostics.DATA");
    EclipseState eclState(deck);
    GridManager gm(eclState.getInputGrid());
    const UnstructuredGrid& grid = *gm.c_grid();
    std::shared_ptr<CounterLog> counterLog = std::make_shared<CounterLog>(Log::DefaultMessageTypes);
    OpmLog::addBackend( "COUNTERLOG" , counterLog );
    RelpermDiagnostics diagnostics;
    diagnostics.diagnosis(eclState, grid);
    BOOST_CHECK_EQUAL(1, counterLog->numMessages(Log::MessageType::Warning));
}
BOOST_AUTO_TEST_SUITE_END()
