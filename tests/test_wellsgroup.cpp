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
#include <opm/parser/eclipse/EclipseState/Schedule/Group/Group2.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group/GroupTree.hpp>

#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well2.hpp>
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
    SummaryState summaryState;
   PhaseUsage pu = phaseUsageFromDeck(eclipseState);

   auto wells = sched.getWells2atEnd();

    for (size_t i=0; i<wells.size(); i++) {
        const auto& well = wells[i];
        std::shared_ptr<WellsGroupInterface> wellsGroup = createWellWellsGroup(well, summaryState, pu);
        BOOST_CHECK_EQUAL(well.name(), wellsGroup->name());
        if (well.isInjector()) {
            const auto controls = well.injectionControls(summaryState);
            BOOST_CHECK_EQUAL(controls.surface_rate, wellsGroup->injSpec().surface_flow_max_rate_);
            BOOST_CHECK_EQUAL(controls.bhp_limit, wellsGroup->injSpec().BHP_limit_);
            BOOST_CHECK_EQUAL(controls.reservoir_rate, wellsGroup->injSpec().reservoir_flow_max_rate_);
            BOOST_CHECK_EQUAL(0.0, wellsGroup->prodSpec().guide_rate_);
        }
        if (well.isProducer()) {
            const auto controls = well.productionControls(summaryState);
            BOOST_CHECK_EQUAL(controls.resv_rate, wellsGroup->prodSpec().reservoir_flow_max_rate_);
            BOOST_CHECK_EQUAL(controls.bhp_limit, wellsGroup->prodSpec().BHP_limit_);
            BOOST_CHECK_EQUAL(controls.oil_rate, wellsGroup->prodSpec().oil_max_rate_);
            BOOST_CHECK_EQUAL(controls.water_rate, wellsGroup->prodSpec().water_max_rate_);
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

    for( const auto& grp_name : sched.groupNames() ) {
        if( !nodes.exists( grp_name ) ) continue;
        const auto& group = sched.getGroup2(grp_name, 2);

        std::shared_ptr<WellsGroupInterface> wellsGroup = createGroupWellsGroup(group, 2, pu);
        BOOST_CHECK_EQUAL(group.name(), wellsGroup->name());
        if (group.isInjectionGroup()) {
            const auto& injection = group.injectionProperties();

            BOOST_CHECK_EQUAL(injection.surface_max_rate, wellsGroup->injSpec().surface_flow_max_rate_);
            BOOST_CHECK_EQUAL(injection.resv_max_rate, wellsGroup->injSpec().reservoir_flow_max_rate_);
            BOOST_CHECK_EQUAL(injection.target_reinj_fraction, wellsGroup->injSpec().reinjection_fraction_target_);
            BOOST_CHECK_EQUAL(injection.target_void_fraction, wellsGroup->injSpec().voidage_replacment_fraction_);
        }

        if (group.isProductionGroup()) {
            const auto& production = group.productionProperties();
            BOOST_CHECK_EQUAL(production.resv_target, wellsGroup->prodSpec().reservoir_flow_max_rate_);
            BOOST_CHECK_EQUAL(production.gas_target, wellsGroup->prodSpec().gas_max_rate_);
            BOOST_CHECK_EQUAL(production.oil_target, wellsGroup->prodSpec().oil_max_rate_);
            BOOST_CHECK_EQUAL(production.water_target, wellsGroup->prodSpec().water_max_rate_);
            BOOST_CHECK_EQUAL(production.liquid_target, wellsGroup->prodSpec().liquid_max_rate_);
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
    for( const auto& grp_name : sched.groupNames() ) {
        if( !nodes.exists( grp_name ) ) continue;
        const auto& group = sched.getGroup2(grp_name, 2);

        std::shared_ptr<WellsGroupInterface> wellsGroup = createGroupWellsGroup(group, 2, pu);
        BOOST_CHECK_EQUAL(group.name(), wellsGroup->name());
        BOOST_CHECK_EQUAL(group.getGroupEfficiencyFactor(), wellsGroup->efficiencyFactor());
        BOOST_CHECK_EQUAL(group.getGroupEfficiencyFactor(), wellsGroup->efficiencyFactor());

    }
}

