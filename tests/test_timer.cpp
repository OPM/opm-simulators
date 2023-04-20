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

#define BOOST_TEST_MODULE OPM-TimerTest
#include <boost/test/unit_test.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/simulators/timestepping/SimulatorTimer.hpp>
#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/input/eclipse/Python/Python.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FieldPropsManager.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <string>
#include <iostream>
#include <vector>
#include <memory>

BOOST_AUTO_TEST_CASE(CreateTimer)
{
    const std::string filename1 = "TESTTIMER.DATA";
    Opm::Parser parser;
    Opm::Deck parserDeck = parser.parseFile( filename1);
    auto python = std::make_shared<Opm::Python>();
    Opm::EclipseGrid grid(10,10,10);
    Opm::FieldPropsManager fp(parserDeck, Opm::Phases{true, true, true}, grid, Opm::TableManager());
    Opm::Runspec runspec(parserDeck);
    Opm::Schedule schedule( parserDeck, grid, fp, runspec, python );
    Opm::SimulatorTimer simtimer;

    boost::gregorian::date defaultStartDate( 2012, 1, 1);
    BOOST_CHECK_EQUAL(  boost::posix_time::ptime(defaultStartDate), simtimer.currentDateTime() );

    simtimer.init(schedule);
    boost::gregorian::date startDate( 2014, 3, 26);
    BOOST_CHECK_EQUAL(  boost::posix_time::ptime(startDate), simtimer.currentDateTime() );

    BOOST_CHECK_EQUAL( 0, simtimer.currentStepNum() );
    BOOST_CHECK_EQUAL( 0., simtimer.simulationTimeElapsed() );
    BOOST_CHECK_EQUAL( 125, simtimer.numSteps() );
    // 1200 + 1000 * 365 * 5
    BOOST_CHECK_EQUAL( 1826200, Opm::unit::convert::to(simtimer.totalTime(), Opm::unit::day) );
    BOOST_CHECK_EQUAL( 0., Opm::unit::convert::to(simtimer.simulationTimeElapsed(), Opm::unit::day) );

    double testCurrentTime = 0.;
    BOOST_CHECK_EQUAL( Opm::unit::convert::to(testCurrentTime, Opm::unit::day),
                       Opm::unit::convert::to(simtimer.simulationTimeElapsed(), Opm::unit::day) );

    for ( int i = 0; i < simtimer.numSteps(); ++i ) {
        BOOST_CHECK_EQUAL( i, simtimer.currentStepNum() );
        BOOST_CHECK_EQUAL( Opm::unit::convert::to(testCurrentTime, Opm::unit::minute),
                           Opm::unit::convert::to(simtimer.simulationTimeElapsed(), Opm::unit::minute) );
        testCurrentTime += simtimer.currentStepLength();
        ++simtimer;
    }

    for ( int i = 0; i <= simtimer.numSteps(); ++i ) {
        simtimer.setCurrentStepNum(i);
        BOOST_CHECK_EQUAL( i, simtimer.currentStepNum() );
    }

    BOOST_CHECK_EQUAL( true, simtimer.done() );
    simtimer.setCurrentStepNum(0);
    BOOST_CHECK_EQUAL( false, simtimer.done() );
    BOOST_CHECK_EQUAL( 0., Opm::unit::convert::to(simtimer.simulationTimeElapsed(), Opm::unit::day) );

    simtimer.setCurrentStepNum(125);
    BOOST_CHECK_EQUAL( Opm::unit::convert::to(simtimer.simulationTimeElapsed(), Opm::unit::day),
                       Opm::unit::convert::to(simtimer.totalTime(), Opm::unit::day) );

    boost::gregorian::date endDate( 7014, 3, 14);
    BOOST_CHECK_EQUAL ( simtimer.currentDateTime(), boost::posix_time::ptime(endDate));

    int i = 0;
    double testCurrentTime1 = 0.;
    double testCurrentTime2 = 0.;
    simtimer.setCurrentStepNum(0);

    while (!simtimer.done()) {
        testCurrentTime1 += simtimer.currentStepLength();
        BOOST_CHECK_EQUAL( i, simtimer.currentStepNum() );
        ++i;
        ++simtimer;
        testCurrentTime2 += simtimer.stepLengthTaken();
        BOOST_CHECK_EQUAL( Opm::unit::convert::to(testCurrentTime1, Opm::unit::minute),
                           Opm::unit::convert::to(simtimer.simulationTimeElapsed(), Opm::unit::minute) );
        BOOST_CHECK_EQUAL( Opm::unit::convert::to(testCurrentTime2, Opm::unit::minute),
                           Opm::unit::convert::to(simtimer.simulationTimeElapsed(), Opm::unit::minute) );
    }

    BOOST_CHECK_EQUAL( true, simtimer.done() );
    BOOST_CHECK_EQUAL( Opm::unit::convert::to(testCurrentTime1, Opm::unit::minute),
                       Opm::unit::convert::to(simtimer.totalTime(), Opm::unit::minute) );
    BOOST_CHECK_EQUAL( Opm::unit::convert::to(testCurrentTime2, Opm::unit::minute),
                       Opm::unit::convert::to(simtimer.totalTime(), Opm::unit::minute) );



}
