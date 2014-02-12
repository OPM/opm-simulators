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
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/GroupTreeNode.hpp>

using namespace Opm;

BOOST_AUTO_TEST_CASE(AddWellsAndGroupToCollection) {
    ParserPtr parser(new Parser());
    boost::filesystem::path scheduleFile("wells_group.data");
    DeckConstPtr deck =  parser->parseFile(scheduleFile.string());
    EclipseStateConstPtr eclipseState(new EclipseState(deck));
    PhaseUsage pu = phaseUsageFromDeck(eclipseState);

    GroupTreeNodePtr field=eclipseState->getSchedule()->getGroupTree(2)->getNode("FIELD");
    GroupTreeNodePtr g1=eclipseState->getSchedule()->getGroupTree(2)->getNode("G1");
    GroupTreeNodePtr g2=eclipseState->getSchedule()->getGroupTree(2)->getNode("G2");

    WellCollection collection;

    // Add groups to WellCollection
    GroupConstPtr fieldGroup =  eclipseState->getSchedule()->getGroup(field->name());
    for (auto iter = field->begin(); iter != field->end(); ++iter) {
        GroupConstPtr childGroupNode = eclipseState->getSchedule()->getGroup((*iter).second->name());
        collection.addChild(childGroupNode, fieldGroup, 2, pu);
    }

    GroupConstPtr g1Group =  eclipseState->getSchedule()->getGroup(g1->name());
    for (auto iter = g1->begin(); iter != g1->end(); ++iter) {
        GroupConstPtr childGroupNode = eclipseState->getSchedule()->getGroup((*iter).second->name());
        collection.addChild(childGroupNode, g1Group, 2, pu);
    }


    GroupConstPtr g2Group =  eclipseState->getSchedule()->getGroup(g2->name());
    for (auto iter = g2->begin(); iter != g2->end(); ++iter) {
        GroupConstPtr childGroupNode = eclipseState->getSchedule()->getGroup((*iter).second->name());
        collection.addChild(childGroupNode, g2Group, 2, pu);
    }

    BOOST_CHECK_EQUAL("FIELD", collection.findNode("FIELD")->name());
    BOOST_CHECK_EQUAL("FIELD", collection.findNode("G1")->getParent()->name());
    BOOST_CHECK_EQUAL("FIELD", collection.findNode("G2")->getParent()->name());

    // Add wells to WellCollection
    WellCollection wellCollection;
    std::vector<WellConstPtr> wells = eclipseState->getSchedule()->getWells();
    for (size_t i=0; i<wells.size(); i++) {
        GroupConstPtr parentGroup = eclipseState->getSchedule()->getGroup(wells[i]->getGroupName(2));
        collection.addChild(wells[i], parentGroup, 2, pu);
    }

    BOOST_CHECK_EQUAL("G1", collection.findNode("INJ1")->getParent()->name());
    BOOST_CHECK_EQUAL("G1", collection.findNode("INJ2")->getParent()->name());
    BOOST_CHECK_EQUAL("G2", collection.findNode("PROD1")->getParent()->name());
    BOOST_CHECK_EQUAL("G2", collection.findNode("PROD2")->getParent()->name());
}

