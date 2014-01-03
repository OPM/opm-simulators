/*
  Copyright 2013 Statoil ASA.

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

#define BOOST_TEST_MODULE WellsManagerTests
#include <boost/test/unit_test.hpp>
#include <opm/core/wells/WellsManager.hpp>
#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/core/grid/GridManager.hpp>


void wells_static_check(const Wells* wells) {
    BOOST_CHECK_EQUAL(2 , wells->number_of_wells);
    BOOST_CHECK_EQUAL(3 , wells->number_of_phases);
    
    BOOST_CHECK_EQUAL("INJ1" , wells->name[0]);
    BOOST_CHECK_EQUAL("PROD1" , wells->name[1]);

    /* The mapping from well number into the wells->WI and wells->well_cells arrays. */
    BOOST_CHECK_EQUAL(0 , wells->well_connpos[0]);
    BOOST_CHECK_EQUAL(1 , wells->well_connpos[1]);
    BOOST_CHECK_EQUAL(2 , wells->well_connpos[2]);
    
    /* Connection factor */
    BOOST_CHECK_CLOSE(1.2279166666666664e-12 , wells->WI[0] , 0.001);
    BOOST_CHECK_CLOSE(1.2279166666666664e-12 , wells->WI[1] , 0.001);
    
    /* Completed cells */
    BOOST_CHECK_EQUAL(0 , wells->well_cells[0]);
    BOOST_CHECK_EQUAL(9 + 2*10 + 2*10*10 , wells->well_cells[1]);
}


/* 
   The number of controls is determined by looking at which elements
   have been given explicit - non-default - values in the WCONxxxx
   keyword. Is that at all interesting?
*/


void check_controls_epoch0( struct WellControls ** ctrls) {
    // The injector
    {
        const struct WellControls * ctrls0 = ctrls[0];        
        BOOST_CHECK_EQUAL( 3 , ctrls0->num);   // The number of controls for the injector == 3??

        BOOST_CHECK_EQUAL( SURFACE_RATE   , ctrls0->type[0] );
        BOOST_CHECK_EQUAL( RESERVOIR_RATE , ctrls0->type[1] );
        BOOST_CHECK_EQUAL( BHP            , ctrls0->type[2] );

        // The different targets
        BOOST_CHECK_EQUAL( 100.0 / 86400 , ctrls0->target[0]);
        BOOST_CHECK_EQUAL( 200.0 / 86400 , ctrls0->target[1]);
        BOOST_CHECK_EQUAL( 400 * 100000  , ctrls0->target[2]);

        // Which control is active
        BOOST_CHECK_EQUAL( 0 , ctrls0->current );
        
        // The phase distribution in the active target
        BOOST_CHECK_EQUAL( 0 , ctrls0->distr[0] );  // Water
        BOOST_CHECK_EQUAL( 0 , ctrls0->distr[1] );  // Oil
        BOOST_CHECK_EQUAL( 1 , ctrls0->distr[2] );  // Gas
    }
    
    // The producer
    {
        const struct WellControls * ctrls1 = ctrls[1];
        BOOST_CHECK_EQUAL( 2 , ctrls1->num);   // The number of controls for the producer == 2??
        BOOST_CHECK_EQUAL( SURFACE_RATE   , ctrls1->type[0] );
        BOOST_CHECK_EQUAL( BHP            , ctrls1->type[1] );

        // The different targets
        BOOST_CHECK_EQUAL( -20000.0 / 86400 , ctrls1->target[0]);
        BOOST_CHECK_EQUAL(  1000 * 100000  , ctrls1->target[1]);

        // Which control is active
        BOOST_CHECK_EQUAL( 0 , ctrls1->current );

        // The phase distribution in the active target
        BOOST_CHECK_EQUAL( 0 , ctrls1->distr[0] );  // Water
        BOOST_CHECK_EQUAL( 1 , ctrls1->distr[1] );  // Oil
        BOOST_CHECK_EQUAL( 0 , ctrls1->distr[2] );  // Gas
    }
}




void check_controls_epoch1( struct WellControls ** ctrls) {
    // The injector
    {
        const struct WellControls * ctrls0 = ctrls[0];        
        BOOST_CHECK_EQUAL( 3 , ctrls0->num);   // The number of controls for the injector == 3??

        BOOST_CHECK_EQUAL( SURFACE_RATE   , ctrls0->type[0] );
        BOOST_CHECK_EQUAL( RESERVOIR_RATE , ctrls0->type[1] );
        BOOST_CHECK_EQUAL( BHP            , ctrls0->type[2] );

        // The different targets
        BOOST_CHECK_CLOSE( 10.0 / 86400 , ctrls0->target[0] , 0.001);
        BOOST_CHECK_CLOSE( 20.0 / 86400 , ctrls0->target[1] , 0.001);
        BOOST_CHECK_CLOSE( 40 * 100000  , ctrls0->target[2] , 0.001);

        // Which control is active
        BOOST_CHECK_EQUAL( 1 , ctrls0->current );
        
        // The phase distribution in the active target
        BOOST_CHECK_EQUAL( 1 , ctrls0->distr[3] );  // Water
        BOOST_CHECK_EQUAL( 0 , ctrls0->distr[4] );  // Oil
        BOOST_CHECK_EQUAL( 0 , ctrls0->distr[5] );  // Gas
    }
    
    // The producer
    {
        const struct WellControls * ctrls1 = ctrls[1];
        BOOST_CHECK_EQUAL( 3 , ctrls1->num);   // The number of controls for the producer - now 3.
        BOOST_CHECK_EQUAL( SURFACE_RATE   , ctrls1->type[0] );
        BOOST_CHECK_EQUAL( RESERVOIR_RATE , ctrls1->type[1] );
        BOOST_CHECK_EQUAL( BHP            , ctrls1->type[2] );

        // The different targets
        BOOST_CHECK_CLOSE( -999.0 / 86400 , ctrls1->target[0], 0.001);
        BOOST_CHECK_CLOSE( -123.0 / 86400 , ctrls1->target[1], 0.001);
        BOOST_CHECK_CLOSE(  100 * 100000  , ctrls1->target[2], 0.001);

        // Which control is active
        BOOST_CHECK_EQUAL( 1 , ctrls1->current );

        // The phase distribution in the active target
        BOOST_CHECK_EQUAL( 1 , ctrls1->distr[3] );  // Water
        BOOST_CHECK_EQUAL( 1 , ctrls1->distr[4] );  // Oil
        BOOST_CHECK_EQUAL( 1 , ctrls1->distr[5] );  // Gas
    }
}




BOOST_AUTO_TEST_CASE(Constructor_Works) {
    Opm::EclipseGridParser Deck("wells_manager_data.data");
    Opm::GridManager gridManager(Deck);

    Deck.setCurrentEpoch(0);
    {
        Opm::WellsManager wellsManager(Deck, *gridManager.c_grid(), NULL);
        const Wells* wells = wellsManager.c_wells();
        wells_static_check( wells );
        check_controls_epoch0( wells->ctrls );
    }


    Deck.setCurrentEpoch(1);
    {
        Opm::WellsManager wellsManager(Deck, *gridManager.c_grid(), NULL);
        const Wells* wells = wellsManager.c_wells();
        
        wells_static_check( wells );    
        check_controls_epoch1( wells->ctrls );
    }
}



BOOST_AUTO_TEST_CASE(WellsEqual) {
    Opm::EclipseGridParser Deck("wells_manager_data.data");
    Opm::GridManager gridManager(Deck);

    Deck.setCurrentEpoch(0);
    Opm::WellsManager wellsManager0(Deck, *gridManager.c_grid(), NULL);

    Deck.setCurrentEpoch(1);
    Opm::WellsManager wellsManager1(Deck, *gridManager.c_grid(), NULL);

    BOOST_CHECK(  wells_equal( wellsManager0.c_wells() , wellsManager0.c_wells()) ); 
    BOOST_CHECK( !wells_equal( wellsManager0.c_wells() , wellsManager1.c_wells()) ); 
}


BOOST_AUTO_TEST_CASE(ControlsEqual) {
    Opm::EclipseGridParser Deck("wells_manager_data.data");
    Opm::GridManager gridManager(Deck);

    Deck.setCurrentEpoch(0);
    Opm::WellsManager wellsManager0(Deck, *gridManager.c_grid(), NULL);

    Deck.setCurrentEpoch(1);
    Opm::WellsManager wellsManager1(Deck, *gridManager.c_grid(), NULL);

    BOOST_CHECK(  well_controls_equal( wellsManager0.c_wells()->ctrls[0] , wellsManager0.c_wells()->ctrls[0]));
    BOOST_CHECK(  well_controls_equal( wellsManager0.c_wells()->ctrls[1] , wellsManager0.c_wells()->ctrls[1]));
    BOOST_CHECK(  well_controls_equal( wellsManager1.c_wells()->ctrls[0] , wellsManager1.c_wells()->ctrls[0]));
    BOOST_CHECK(  well_controls_equal( wellsManager1.c_wells()->ctrls[1] , wellsManager1.c_wells()->ctrls[1]));

    BOOST_CHECK(  !well_controls_equal( wellsManager0.c_wells()->ctrls[0] , wellsManager0.c_wells()->ctrls[1]));
    BOOST_CHECK(  !well_controls_equal( wellsManager0.c_wells()->ctrls[1] , wellsManager0.c_wells()->ctrls[0]));
    BOOST_CHECK(  !well_controls_equal( wellsManager1.c_wells()->ctrls[0] , wellsManager0.c_wells()->ctrls[0]));
    BOOST_CHECK(  !well_controls_equal( wellsManager1.c_wells()->ctrls[1] , wellsManager0.c_wells()->ctrls[1]));
}


