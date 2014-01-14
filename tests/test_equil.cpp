/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
*/

#include "config.h"

/* --- Boost.Test boilerplate --- */
#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif

#define NVERBOSE  // Suppress own messages when throw()ing

#define BOOST_TEST_MODULE UnitsTest
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

/* --- our own headers --- */

#include <opm/core/simulator/initStateEquil.hpp>

#include <opm/core/grid.h>
#include <opm/core/grid/cart_grid.h>

#include <opm/core/props/BlackoilPropertiesBasic.hpp>
#include <opm/core/props/BlackoilPhases.hpp>

#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/Units.hpp>

#include <array>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

BOOST_AUTO_TEST_SUITE ()

BOOST_AUTO_TEST_CASE (PhasePressure)
{
    typedef std::vector<double> PVal;
    typedef std::vector<PVal> PPress;

    std::shared_ptr<UnstructuredGrid>
        G(create_grid_cart3d(10, 1, 10), destroy_grid);

    Opm::parameter::ParameterGroup param;
    {
        using Opm::unit::kilogram;
        using Opm::unit::meter;
        using Opm::unit::cubic;

        std::stringstream dens; dens << 700*kilogram/cubic(meter);
        param.insertParameter("rho2", dens.str());
    }

    typedef Opm::BlackoilPropertiesBasic Props;
    Props props(param, G->dimensions, G->number_of_cells);

    typedef Opm::equil::DensityCalculator<Opm::BlackoilPropertiesInterface> RhoCalc;
    RhoCalc calc(props, 0);

    Opm::equil::EquilRecord record =
        {
            { 0 , 1e5 } , // Datum depth, pressure
            { 5 , 0   } , // Zwoc       , Pcow_woc
            { 0 , 0   }   // Zgoc       , Pcgo_goc
        };

    Opm::equil::EquilReg<RhoCalc>
        region(record, calc,
               Opm::equil::miscibility::NoMixing(),
               Opm::equil::miscibility::NoMixing(),
               props.phaseUsage());

    const double grav   = 10;
    const PPress ppress = Opm::equil::phasePressures(*G, region, grav);

    const int first = 0, last = G->number_of_cells - 1;
    const double reltol = 1.0e-8;
    BOOST_CHECK_CLOSE(ppress[0][first] ,  90e3   , reltol);
    BOOST_CHECK_CLOSE(ppress[0][last ] , 180e3   , reltol);
    BOOST_CHECK_CLOSE(ppress[1][first] , 103.5e3 , reltol);
    BOOST_CHECK_CLOSE(ppress[1][last ] , 166.5e3 , reltol);
}

BOOST_AUTO_TEST_SUITE_END()
