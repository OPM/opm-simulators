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

#define BOOST_TEST_MODULE WellsModuleTest
#include <boost/test/unit_test.hpp>

#include <opm/core/wells.h>
#define HAVE_WELLCONTROLS
#include <opm/core/well_controls.h>

#include <iostream>
#include <vector>
#include <memory>


BOOST_AUTO_TEST_CASE(Construction)
{
    struct WellControls * ctrls = well_controls_create();

    well_controls_set_current( ctrls , 1 );
    BOOST_CHECK_EQUAL( 1 , well_controls_get_current( ctrls ));
    well_controls_set_current( ctrls , 2 );
    BOOST_CHECK_EQUAL( 2 , well_controls_get_current( ctrls ));
    
    well_controls_invert_current( ctrls );
    BOOST_CHECK( well_controls_get_current( ctrls ) < 0 );
    well_controls_invert_current( ctrls );
    BOOST_CHECK_EQUAL( 2 , well_controls_get_current( ctrls ));

    {
        enum WellControlType type1 = BHP;
        enum WellControlType type2 = SURFACE_RATE;
        int num_phases = 3;
        double dist1[3] = {0 , 1  ,  2};
        double dist2[3] = {10, 11 , 12};
        double target = 77;

        well_controls_assert_number_of_phases( ctrls , num_phases );
        well_controls_add_new( type1 ,   target , dist1 , ctrls );
        well_controls_add_new( type2 , 2*target , dist2 , ctrls );

        BOOST_CHECK_EQUAL( target , well_controls_iget_target(ctrls , 0 ));
        BOOST_CHECK_EQUAL( type1 , well_controls_iget_type(ctrls , 0 ));

        BOOST_CHECK_EQUAL( 2*target , well_controls_iget_target(ctrls , 1 ));
        BOOST_CHECK_EQUAL( type2 , well_controls_iget_type(ctrls , 1 ));
        
        {
            const double * d1 = well_controls_iget_distr( ctrls , 0 );
            const double * d2 = well_controls_iget_distr( ctrls , 1 );
            BOOST_CHECK( memcmp(d1 , dist1 , num_phases * sizeof * d1 ) == 0);
            BOOST_CHECK( memcmp(d2 , dist2 , num_phases * sizeof * d2 ) == 0);
        }
    }
    well_controls_iset_target( ctrls , 0 , 123);
    BOOST_CHECK_EQUAL( 123 , well_controls_iget_target( ctrls , 0 ));
    well_controls_iset_target( ctrls , 1 , 456);
    BOOST_CHECK_EQUAL( 456 , well_controls_iget_target( ctrls , 1 ));

    well_controls_iset_type( ctrls , 0 , SURFACE_RATE);
    BOOST_CHECK_EQUAL( SURFACE_RATE , well_controls_iget_type( ctrls , 0 ));
    well_controls_iset_type( ctrls , 1 , BHP);
    BOOST_CHECK_EQUAL( BHP, well_controls_iget_type( ctrls , 1 ));
    

    well_controls_destroy( ctrls );
}


