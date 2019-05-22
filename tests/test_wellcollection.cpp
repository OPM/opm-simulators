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

#define NVERBOSE  // Suppress own messages when throw()ing

#define BOOST_TEST_MODULE WellCollectionTest
#include <boost/test/unit_test.hpp>
#include <opm/core/wells/WellCollection.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well2.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/GroupTree.hpp>

using namespace Opm;

BOOST_AUTO_TEST_CASE(AddWellsAndGroupToCollection) {
    Parser parser;
    std::string scheduleFile("wells_group.data");
    Deck deck = parser.parseFile(scheduleFile);
    EclipseState eclipseState(deck);
    PhaseUsage pu = phaseUsageFromDeck(eclipseState);
    const auto& grid = eclipseState.getInputGrid();
    const TableManager table ( deck );
    const Eclipse3DProperties eclipseProperties ( deck , table, grid);
    const Runspec runspec(deck);
    const Schedule sched(deck, grid, eclipseProperties, runspec);


    WellCollection collection;

    // Add groups to WellCollection
    const auto& fieldGroup =  sched.getGroup("FIELD");
    collection.addField(fieldGroup, 2, pu);

    collection.addGroup( sched.getGroup( "G1" ), fieldGroup.name(), 2, pu);
    collection.addGroup( sched.getGroup( "G2" ), fieldGroup.name(), 2, pu);

    BOOST_CHECK_EQUAL("FIELD", collection.findNode("FIELD")->name());
    BOOST_CHECK_EQUAL("FIELD", collection.findNode("G1")->getParent()->name());
    BOOST_CHECK_EQUAL("FIELD", collection.findNode("G2")->getParent()->name());

    // Add wells to WellCollection
    WellCollection wellCollection;
    const auto wells = sched.getWells2atEnd();
    for (size_t i=0; i<wells.size(); i++) {
        collection.addWell(wells[i], 2, pu);
    }

    BOOST_CHECK_EQUAL("G1", collection.findNode("INJ1")->getParent()->name());
    BOOST_CHECK_EQUAL("G1", collection.findNode("INJ2")->getParent()->name());
    BOOST_CHECK_EQUAL("G2", collection.findNode("PROD1")->getParent()->name());
    BOOST_CHECK_EQUAL("G2", collection.findNode("PROD2")->getParent()->name());
}

BOOST_AUTO_TEST_CASE(EfficiencyFactor) {
    Parser parser;
    std::string scheduleFile("wells_group.data");
    Deck deck = parser.parseFile(scheduleFile);
    EclipseState eclipseState(deck);
    PhaseUsage pu = phaseUsageFromDeck(eclipseState);
    const auto& grid = eclipseState.getInputGrid();
    const TableManager table ( deck );
    const Eclipse3DProperties eclipseProperties ( deck , table, grid);
    const Runspec runspec(deck);
    const Schedule sched(deck, grid, eclipseProperties, runspec);

    size_t timestep = 2;
    WellCollection collection;
    // Add groups to WellCollection
    const auto& fieldGroup =  sched.getGroup("FIELD");
    collection.addField(fieldGroup, timestep, pu);
    collection.addGroup( sched.getGroup( "G1" ), fieldGroup.name(), timestep, pu);
    collection.addGroup( sched.getGroup( "G2" ), fieldGroup.name(), timestep, pu);

    BOOST_CHECK_EQUAL(1.0, collection.findNode("FIELD")->efficiencyFactor());
    BOOST_CHECK_EQUAL(1.0, collection.findNode("G1")->getParent()->efficiencyFactor());
    BOOST_CHECK_EQUAL(1.0, collection.findNode("G2")->getParent()->efficiencyFactor());

    // Add wells to WellCollection
    const auto wells1 = sched.getWells2(timestep);
    for (size_t i=0; i<wells1.size(); i++) {
        collection.addWell(wells1[i], timestep, pu);
    }

    // 0.5(inj1) * 0.8(G1)
    BOOST_CHECK_CLOSE(0.4, collection.findWellNode("INJ1").getAccumulativeEfficiencyFactor(), 1e-10);
    // 0.8(inj2) * 0.8(G1)
    BOOST_CHECK_CLOSE(0.64, collection.findWellNode("INJ2").getAccumulativeEfficiencyFactor(), 1e-10);
    // 0.5 (prod1) * 1.0 (G2)
    BOOST_CHECK_CLOSE(0.5, collection.findWellNode("PROD1").getAccumulativeEfficiencyFactor(), 1e-10);
    // 1.0 (prod2) * 1.0 (G2)
    BOOST_CHECK_CLOSE(1.0, collection.findWellNode("PROD2").getAccumulativeEfficiencyFactor(), 1e-10);
}
