/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#include <config.h>

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif

#define NVERBOSE  // Suppress own messages when throw()ing

#define BOOST_TEST_MODULE WellCollectionTest
#include <boost/test/unit_test.hpp>
#include <opm/core/wells/WellCollection.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/GroupTree.hpp>

using namespace Opm;

BOOST_AUTO_TEST_CASE(AddWellsAndGroupToCollection) {
    Parser parser;
    std::string scheduleFile("wells_group.data");
    ParseContext parseContext;
    Deck deck = parser.parseFile(scheduleFile, parseContext);
    EclipseState eclipseState(deck, parseContext);
    PhaseUsage pu = phaseUsageFromDeck(eclipseState);

    WellCollection collection;

    // Add groups to WellCollection
    const auto& fieldGroup =  eclipseState.getSchedule().getGroup("FIELD");
    collection.addField(fieldGroup, 2, pu);

    collection.addGroup( eclipseState.getSchedule().getGroup( "G1" ), fieldGroup.name(), 2, pu);
    collection.addGroup( eclipseState.getSchedule().getGroup( "G2" ), fieldGroup.name(), 2, pu);

    BOOST_CHECK_EQUAL("FIELD", collection.findNode("FIELD")->name());
    BOOST_CHECK_EQUAL("FIELD", collection.findNode("G1")->getParent()->name());
    BOOST_CHECK_EQUAL("FIELD", collection.findNode("G2")->getParent()->name());

    // Add wells to WellCollection
    WellCollection wellCollection;
    auto wells = eclipseState.getSchedule().getWells();
    for (size_t i=0; i<wells.size(); i++) {
        collection.addWell(wells[i], 2, pu);
    }

    BOOST_CHECK_EQUAL("G1", collection.findNode("INJ1")->getParent()->name());
    BOOST_CHECK_EQUAL("G1", collection.findNode("INJ2")->getParent()->name());
    BOOST_CHECK_EQUAL("G2", collection.findNode("PROD1")->getParent()->name());
    BOOST_CHECK_EQUAL("G2", collection.findNode("PROD2")->getParent()->name());
}

