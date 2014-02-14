/*
  Copyright 2014 Statoil.

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

#define BOOST_TEST_MODULE WellsGroupTest

#include <memory>
#include <vector>

#include <boost/test/unit_test.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include <opm/core/props/phaseUsageFromDeck.hpp>
#include <opm/core/wells/WellsGroup.hpp>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/GroupTreeNode.hpp>

#include <opm/parser/eclipse/EclipseState/Schedule/Well.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

using namespace Opm;

BOOST_AUTO_TEST_CASE(ConstructGroupFromWell) {
    ParserPtr parser(new Parser());
    boost::filesystem::path scheduleFile("wells_group.data");
    DeckConstPtr deck =  parser->parseFile(scheduleFile.string());
    EclipseStateConstPtr eclipseState(new EclipseState(deck));
    PhaseUsage pu = phaseUsageFromDeck(eclipseState);

    std::vector<WellConstPtr> wells = eclipseState->getSchedule()->getWells();

    for (size_t i=0; i<wells.size(); i++) {
        WellConstPtr well = wells[i];
        std::shared_ptr<WellsGroupInterface> wellsGroup = createWellWellsGroup(well, 2, pu);
        BOOST_CHECK_EQUAL(well->name(), wellsGroup->name());
        if (well->isInjector(2)) {
            BOOST_CHECK_EQUAL(well->getSurfaceInjectionRate(2), wellsGroup->injSpec().surface_flow_max_rate_);
            BOOST_CHECK_EQUAL(well->getBHPLimit(2), wellsGroup->injSpec().BHP_limit_);
            BOOST_CHECK_EQUAL(well->getReservoirInjectionRate(2), wellsGroup->injSpec().reservoir_flow_max_rate_);
            BOOST_CHECK_EQUAL(0.0, wellsGroup->prodSpec().guide_rate_);
        }
        if (well->isProducer(2)) {
            BOOST_CHECK_EQUAL(well->getResVRate(2), wellsGroup->prodSpec().reservoir_flow_max_rate_);
            BOOST_CHECK_EQUAL(well->getBHPLimit(2), wellsGroup->prodSpec().BHP_limit_);
            BOOST_CHECK_EQUAL(well->getOilRate(2), wellsGroup->prodSpec().oil_max_rate_);
            BOOST_CHECK_EQUAL(well->getWaterRate(2), wellsGroup->prodSpec().water_max_rate_);
            BOOST_CHECK_EQUAL(0.0, wellsGroup->injSpec().guide_rate_);
        }
    }
}


BOOST_AUTO_TEST_CASE(ConstructGroupFromGroup) {
    ParserPtr parser(new Parser());
    boost::filesystem::path scheduleFile("wells_group.data");
    DeckConstPtr deck =  parser->parseFile(scheduleFile.string());
    EclipseStateConstPtr eclipseState(new EclipseState(deck));
    PhaseUsage pu = phaseUsageFromDeck(eclipseState);

    std::vector<GroupTreeNodeConstPtr> nodes = eclipseState->getSchedule()->getGroupTree(2)->getNodes();

    for (size_t i=0; i<nodes.size(); i++) {
        GroupConstPtr group = eclipseState->getSchedule()->getGroup(nodes[i]->name());
        std::shared_ptr<WellsGroupInterface> wellsGroup = createGroupWellsGroup(group, 2, pu);
        BOOST_CHECK_EQUAL(group->name(), wellsGroup->name());
        if (group->isInjectionGroup(2)) {
            BOOST_CHECK_EQUAL(group->getSurfaceMaxRate(2), wellsGroup->injSpec().surface_flow_max_rate_);
            BOOST_CHECK_EQUAL(group->getReservoirMaxRate(2), wellsGroup->injSpec().reservoir_flow_max_rate_);
            BOOST_CHECK_EQUAL(group->getTargetReinjectFraction(2), wellsGroup->injSpec().reinjection_fraction_target_);
            BOOST_CHECK_EQUAL(group->getTargetVoidReplacementFraction(2), wellsGroup->injSpec().voidage_replacment_fraction_);
        }
        if (group->isProductionGroup(2)) {
            BOOST_CHECK_EQUAL(group->getReservoirMaxRate(2), wellsGroup->prodSpec().reservoir_flow_max_rate_);
            BOOST_CHECK_EQUAL(group->getGasTargetRate(2), wellsGroup->prodSpec().gas_max_rate_);
            BOOST_CHECK_EQUAL(group->getOilTargetRate(2), wellsGroup->prodSpec().oil_max_rate_);
            BOOST_CHECK_EQUAL(group->getWaterTargetRate(2), wellsGroup->prodSpec().water_max_rate_);
            BOOST_CHECK_EQUAL(group->getLiquidTargetRate(2), wellsGroup->prodSpec().liquid_max_rate_);
        }
    }
}



