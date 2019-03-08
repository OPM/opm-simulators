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
#include <opm/parser/eclipse/EclipseState/Schedule/GroupTree.hpp>

#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

using namespace Opm;

BOOST_AUTO_TEST_CASE(ConstructGroupFromWell) {
    std::string scheduleFile("wells_group.data");
    Parser parser;
    Deck deck =  parser.parseFile(scheduleFile);
    EclipseState eclipseState(deck);
    const auto& grid = eclipseState.getInputGrid();
    const TableManager table ( deck );
    const Eclipse3DProperties eclipseProperties ( deck , table, grid);
    const Opm::Runspec runspec (deck);
    const Schedule sched(deck, grid, eclipseProperties, runspec);

   PhaseUsage pu = phaseUsageFromDeck(eclipseState);

    auto wells = sched.getWells();

    for (size_t i=0; i<wells.size(); i++) {
        const auto* well = wells[i];
        std::shared_ptr<WellsGroupInterface> wellsGroup = createWellWellsGroup(well, 2, pu);
        BOOST_CHECK_EQUAL(well->name(), wellsGroup->name());
        if (well->isInjector(2)) {
            const WellInjectionProperties& properties = well->getInjectionProperties(2);
            BOOST_CHECK_EQUAL(properties.surfaceInjectionRate, wellsGroup->injSpec().surface_flow_max_rate_);
            BOOST_CHECK_EQUAL(properties.BHPLimit, wellsGroup->injSpec().BHP_limit_);
            BOOST_CHECK_EQUAL(properties.reservoirInjectionRate, wellsGroup->injSpec().reservoir_flow_max_rate_);
            BOOST_CHECK_EQUAL(0.0, wellsGroup->prodSpec().guide_rate_);
        }
        if (well->isProducer(2)) {
            const WellProductionProperties& properties = well->getProductionProperties(2);
            BOOST_CHECK_EQUAL(properties.ResVRate, wellsGroup->prodSpec().reservoir_flow_max_rate_);
            BOOST_CHECK_EQUAL(properties.BHPLimit, wellsGroup->prodSpec().BHP_limit_);
            BOOST_CHECK_EQUAL(properties.OilRate, wellsGroup->prodSpec().oil_max_rate_);
            BOOST_CHECK_EQUAL(properties.WaterRate, wellsGroup->prodSpec().water_max_rate_);
            BOOST_CHECK_EQUAL(0.0, wellsGroup->injSpec().guide_rate_);
        }
    }
}


BOOST_AUTO_TEST_CASE(ConstructGroupFromGroup) {
    Parser parser;
    std::string scheduleFile("wells_group.data");
    Deck deck =  parser.parseFile(scheduleFile);
    EclipseState eclipseState(deck);
    PhaseUsage pu = phaseUsageFromDeck(eclipseState);
    const auto& grid = eclipseState.getInputGrid();
    const TableManager table ( deck );
    const Eclipse3DProperties eclipseProperties ( deck , table, grid);
    const Opm::Runspec runspec (deck);
    const Schedule sched(deck, grid, eclipseProperties, runspec);


    const auto& nodes = sched.getGroupTree(2);

    for( const auto& grp : sched.getGroups() ) {
        if( !nodes.exists( grp->name() ) ) continue;
        const auto& group = *grp;

        std::shared_ptr<WellsGroupInterface> wellsGroup = createGroupWellsGroup(group, 2, pu);
        BOOST_CHECK_EQUAL(group.name(), wellsGroup->name());
        if (group.isInjectionGroup(2)) {
            BOOST_CHECK_EQUAL(group.getSurfaceMaxRate(2), wellsGroup->injSpec().surface_flow_max_rate_);
            BOOST_CHECK_EQUAL(group.getReservoirMaxRate(2), wellsGroup->injSpec().reservoir_flow_max_rate_);
            BOOST_CHECK_EQUAL(group.getTargetReinjectFraction(2), wellsGroup->injSpec().reinjection_fraction_target_);
            BOOST_CHECK_EQUAL(group.getTargetVoidReplacementFraction(2), wellsGroup->injSpec().voidage_replacment_fraction_);
        }
        if (group.isProductionGroup(2)) {
            BOOST_CHECK_EQUAL(group.getReservoirVolumeTargetRate(2), wellsGroup->prodSpec().reservoir_flow_max_rate_);
            BOOST_CHECK_EQUAL(group.getGasTargetRate(2), wellsGroup->prodSpec().gas_max_rate_);
            BOOST_CHECK_EQUAL(group.getOilTargetRate(2), wellsGroup->prodSpec().oil_max_rate_);
            BOOST_CHECK_EQUAL(group.getWaterTargetRate(2), wellsGroup->prodSpec().water_max_rate_);
            BOOST_CHECK_EQUAL(group.getLiquidTargetRate(2), wellsGroup->prodSpec().liquid_max_rate_);
        }
    }
}

BOOST_AUTO_TEST_CASE(EfficiencyFactor) {
    Parser parser;
    std::string scheduleFile("wells_group.data");
    Deck deck =  parser.parseFile(scheduleFile);
    EclipseState eclipseState(deck);
    PhaseUsage pu = phaseUsageFromDeck(eclipseState);
    const auto& grid = eclipseState.getInputGrid();
    const TableManager table ( deck );
    const Eclipse3DProperties eclipseProperties ( deck , table, grid);
    const Opm::Runspec runspec (deck);
    const Schedule sched(deck, grid, eclipseProperties, runspec);


    const auto& nodes = sched.getGroupTree(2);

    for( const auto& grp : sched.getGroups() ) {
        if( !nodes.exists( grp->name() ) ) continue;
        const auto& group = *grp;

        std::shared_ptr<WellsGroupInterface> wellsGroup = createGroupWellsGroup(group, 2, pu);
        BOOST_CHECK_EQUAL(group.name(), wellsGroup->name());
        BOOST_CHECK_EQUAL(group.getGroupEfficiencyFactor(2), wellsGroup->efficiencyFactor());
        BOOST_CHECK_EQUAL(group.getGroupEfficiencyFactor(2), wellsGroup->efficiencyFactor());

    }
}

