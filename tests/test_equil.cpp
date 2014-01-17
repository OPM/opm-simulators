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

    std::vector<int> cells(G->number_of_cells);
    std::iota(cells.begin(), cells.end(), 0);

    const double grav   = 10;
    const PPress ppress = Opm::equil::phasePressures(*G, region, cells, grav);

    const int first = 0, last = G->number_of_cells - 1;
    const double reltol = 1.0e-8;
    BOOST_CHECK_CLOSE(ppress[0][first] ,  90e3   , reltol);
    BOOST_CHECK_CLOSE(ppress[0][last ] , 180e3   , reltol);
    BOOST_CHECK_CLOSE(ppress[1][first] , 103.5e3 , reltol);
    BOOST_CHECK_CLOSE(ppress[1][last ] , 166.5e3 , reltol);
}

BOOST_AUTO_TEST_CASE (CellSubset)
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

    Opm::equil::EquilRecord record[] =
        {
            {
                { 0   ,  1e5     } , // Datum depth, pressure
                { 2.5 , -0.075e5 } , // Zwoc       , Pcow_woc
                { 0   ,  0       }   // Zgoc       , Pcgo_goc
            }
            ,
            {
                { 5   ,  1.35e5  } , // Datum depth, pressure
                { 7.5 , -0.225e5 } , // Zwoc       , Pcow_woc
                { 5   ,  0       }   // Zgoc       , Pcgo_goc
            }
        };

    Opm::equil::EquilReg<RhoCalc> region[] =
        {
            Opm::equil::EquilReg<RhoCalc>(record[0], calc,
                                          Opm::equil::miscibility::NoMixing(),
                                          Opm::equil::miscibility::NoMixing(),
                                          props.phaseUsage())
            ,
            Opm::equil::EquilReg<RhoCalc>(record[0], calc,
                                          Opm::equil::miscibility::NoMixing(),
                                          Opm::equil::miscibility::NoMixing(),
                                          props.phaseUsage())
            ,
            Opm::equil::EquilReg<RhoCalc>(record[1], calc,
                                          Opm::equil::miscibility::NoMixing(),
                                          Opm::equil::miscibility::NoMixing(),
                                          props.phaseUsage())
            ,
            Opm::equil::EquilReg<RhoCalc>(record[1], calc,
                                          Opm::equil::miscibility::NoMixing(),
                                          Opm::equil::miscibility::NoMixing(),
                                          props.phaseUsage())
        };

    const int cdim[] = { 2, 1, 2 };
    int ncoarse = cdim[0];
    for (std::size_t d = 1; d < 3; ++d) { ncoarse *= cdim[d]; }

    std::vector< std::vector<int> > cells(ncoarse);
    for (int c = 0; c < G->number_of_cells; ++c) {
        int ci = c;
        const int i = ci % G->cartdims[0];  ci /= G->cartdims[0];
        const int j = ci % G->cartdims[1];
        const int k = ci / G->cartdims[1];

        const int ic = (i / (G->cartdims[0] / cdim[0]));
        const int jc = (j / (G->cartdims[1] / cdim[1]));
        const int kc = (k / (G->cartdims[2] / cdim[2]));
        const int ix = ic + cdim[0]*(jc + cdim[1]*kc);

        assert ((0 <= ix) && (ix < ncoarse));
        cells[ix].push_back(c);
    }

    PPress ppress(2, PVal(G->number_of_cells, 0));
    for (std::vector< std::vector<int> >::const_iterator
             r = cells.begin(), e = cells.end();
         r != e; ++r)
    {
        const int    rno  = int(r - cells.begin());
        const double grav = 10;
        const PPress p    =
            Opm::equil::phasePressures(*G, region[rno], *r, grav);

        PVal::size_type i = 0;
        for (std::vector<int>::const_iterator
                 c = r->begin(), ce = r->end();
             c != ce; ++c, ++i)
        {
            assert (i < p[0].size());

            ppress[0][*c] = p[0][i];
            ppress[1][*c] = p[1][i];
        }
    }

    const int first = 0, last = G->number_of_cells - 1;
    const double reltol = 1.0e-8;
    BOOST_CHECK_CLOSE(ppress[0][first] , 105e3   , reltol);
    BOOST_CHECK_CLOSE(ppress[0][last ] , 195e3   , reltol);
    BOOST_CHECK_CLOSE(ppress[1][first] , 103.5e3 , reltol);
    BOOST_CHECK_CLOSE(ppress[1][last ] , 166.5e3 , reltol);
}

BOOST_AUTO_TEST_SUITE_END()
