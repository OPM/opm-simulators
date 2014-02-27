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
#include <opm/core/grid/GridManager.hpp>

#include <opm/core/props/BlackoilPropertiesBasic.hpp>
#include <opm/core/props/BlackoilPropertiesFromDeck.hpp>
#include <opm/core/props/BlackoilPhases.hpp>

#include <opm/core/pressure/msmfem/partition.h>

#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/Units.hpp>

#include <array>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
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

    typedef Opm::Equil::DensityCalculator<Opm::BlackoilPropertiesInterface> RhoCalc;
    RhoCalc calc(props, 0);

    Opm::Equil::EquilRecord record =
        {
            { 0 , 1e5 } , // Datum depth, pressure
            { 5 , 0   } , // Zwoc       , Pcow_woc
            { 0 , 0   }   // Zgoc       , Pcgo_goc
        };

    Opm::Equil::EquilReg<RhoCalc>
        region(record, calc,
               std::make_shared<Opm::Equil::Miscibility::NoMixing>(),
               std::make_shared<Opm::Equil::Miscibility::NoMixing>(),
               props.phaseUsage());

    std::vector<int> cells(G->number_of_cells);
    std::iota(cells.begin(), cells.end(), 0);

    const double grav   = 10;
    const PPress ppress = Opm::Equil::phasePressures(*G, region, cells, grav);

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

    typedef Opm::Equil::DensityCalculator<Opm::BlackoilPropertiesInterface> RhoCalc;
    RhoCalc calc(props, 0);

    Opm::Equil::EquilRecord record[] =
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

    Opm::Equil::EquilReg<RhoCalc> region[] =
        {
            Opm::Equil::EquilReg<RhoCalc>(record[0], calc,
                                          std::make_shared<Opm::Equil::Miscibility::NoMixing>(),
                                          std::make_shared<Opm::Equil::Miscibility::NoMixing>(),
                                          props.phaseUsage())
            ,
            Opm::Equil::EquilReg<RhoCalc>(record[0], calc,
                                          std::make_shared<Opm::Equil::Miscibility::NoMixing>(),
                                          std::make_shared<Opm::Equil::Miscibility::NoMixing>(),
                                          props.phaseUsage())
            ,
            Opm::Equil::EquilReg<RhoCalc>(record[1], calc,
                                          std::make_shared<Opm::Equil::Miscibility::NoMixing>(),
                                          std::make_shared<Opm::Equil::Miscibility::NoMixing>(),
                                          props.phaseUsage())
            ,
            Opm::Equil::EquilReg<RhoCalc>(record[1], calc,
                                          std::make_shared<Opm::Equil::Miscibility::NoMixing>(),
                                          std::make_shared<Opm::Equil::Miscibility::NoMixing>(),
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
            Opm::Equil::phasePressures(*G, region[rno], *r, grav);

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




BOOST_AUTO_TEST_CASE (RegMapping)
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

    typedef Opm::Equil::DensityCalculator<Opm::BlackoilPropertiesInterface> RhoCalc;
    RhoCalc calc(props, 0);

    Opm::Equil::EquilRecord record[] =
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

    Opm::Equil::EquilReg<RhoCalc> region[] =
        {
            Opm::Equil::EquilReg<RhoCalc>(record[0], calc,
                                          std::make_shared<Opm::Equil::Miscibility::NoMixing>(),
                                          std::make_shared<Opm::Equil::Miscibility::NoMixing>(),
                                          props.phaseUsage())
            ,
            Opm::Equil::EquilReg<RhoCalc>(record[0], calc,
                                          std::make_shared<Opm::Equil::Miscibility::NoMixing>(),
                                          std::make_shared<Opm::Equil::Miscibility::NoMixing>(),
                                          props.phaseUsage())
            ,
            Opm::Equil::EquilReg<RhoCalc>(record[1], calc,
                                          std::make_shared<Opm::Equil::Miscibility::NoMixing>(),
                                          std::make_shared<Opm::Equil::Miscibility::NoMixing>(),
                                          props.phaseUsage())
            ,
            Opm::Equil::EquilReg<RhoCalc>(record[1], calc,
                                          std::make_shared<Opm::Equil::Miscibility::NoMixing>(),
                                          std::make_shared<Opm::Equil::Miscibility::NoMixing>(),
                                          props.phaseUsage())
        };

    std::vector<int> eqlnum(G->number_of_cells);
    {
        std::vector<int> cells(G->number_of_cells);
        std::iota(cells.begin(), cells.end(), 0);

        const int cdim[] = { 2, 1, 2 };
        int ncoarse = cdim[0];
        for (std::size_t d = 1; d < 3; ++d) { ncoarse *= cdim[d]; }

        partition_unif_idx(G->dimensions, G->number_of_cells,
                           G->cartdims, cdim,
                           &cells[0], &eqlnum[0]);
    }
    Opm::RegionMapping<> eqlmap(eqlnum);

    PPress ppress(2, PVal(G->number_of_cells, 0));
    for (int r = 0, e = eqlmap.numRegions(); r != e; ++r)
    {
        const Opm::RegionMapping<>::CellRange&
            rng = eqlmap.cells(r);

        const int    rno  = r;
        const double grav = 10;
        const PPress p    =
            Opm::Equil::phasePressures(*G, region[rno], rng, grav);

        PVal::size_type i = 0;
        for (Opm::RegionMapping<>::CellRange::const_iterator
                 c = rng.begin(), ce = rng.end();
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



BOOST_AUTO_TEST_CASE (DeckAllDead)
{
    std::shared_ptr<UnstructuredGrid>
        grid(create_grid_cart3d(1, 1, 10), destroy_grid);
    Opm::EclipseGridParser deck("deadfluids.DATA");
    Opm::BlackoilPropertiesFromDeck props(deck, *grid, false);
    Opm::Equil::DeckDependent::InitialStateComputer<Opm::EclipseGridParser> comp(props, deck, *grid, 10.0);
    const auto& pressures = comp.press();
    BOOST_REQUIRE(pressures.size() == 3);
    BOOST_REQUIRE(int(pressures[0].size()) == grid->number_of_cells);

    const int first = 0, last = grid->number_of_cells - 1;
    // The relative tolerance is too loose to be very useful,
    // but the answer we are checking is the result of an ODE
    // solver, and it is unclear if we should check it against
    // the true answer or something else.
    const double reltol = 1.0e-3;
    BOOST_CHECK_CLOSE(pressures[0][first] , 14955e3   , reltol);
    BOOST_CHECK_CLOSE(pressures[0][last ] , 15045e3   , reltol);
    BOOST_CHECK_CLOSE(pressures[1][last] , 1.50473e7   , reltol);
}



BOOST_AUTO_TEST_CASE (CapillaryInversion)
{
    // Test setup.
    Opm::GridManager gm(1, 1, 40, 1.0, 1.0, 2.5);
    const UnstructuredGrid& grid = *(gm.c_grid());
    Opm::EclipseGridParser deck("capillary.DATA");
    Opm::BlackoilPropertiesFromDeck props(deck, grid, false);

    // Test the capillary inversion for oil-water.
    const int cell = 0;
    const double reltol = 1.0e-7;
    {
        const int phase = 0;
        const bool increasing = false;
        const std::vector<double> pc = { 10.0e5, 0.5e5, 0.4e5, 0.3e5, 0.2e5, 0.1e5, 0.099e5, 0.0e5, -10.0e5 };
        const std::vector<double> s = { 0.2, 0.2, 0.2, 0.466666666666, 0.733333333333, 1.0, 1.0, 1.0, 1.0 };
        BOOST_REQUIRE(pc.size() == s.size());
        for (size_t i = 0; i < pc.size(); ++i) {
            const double s_computed = Opm::Equil::satFromPc(props, phase, cell, pc[i], increasing);
            BOOST_CHECK_CLOSE(s_computed, s[i], reltol);
        }
    }

    // Test the capillary inversion for gas-oil.
    {
        const int phase = 2;
        const bool increasing = true;
        const std::vector<double> pc = { 10.0e5, 0.6e5, 0.5e5, 0.4e5, 0.3e5, 0.2e5, 0.1e5, 0.0e5, -10.0e5 };
        const std::vector<double> s = { 0.8, 0.8, 0.8, 0.533333333333, 0.266666666666, 0.0, 0.0, 0.0, 0.0 };
        BOOST_REQUIRE(pc.size() == s.size());
        for (size_t i = 0; i < pc.size(); ++i) {
            const double s_computed = Opm::Equil::satFromPc(props, phase, cell, pc[i], increasing);
            BOOST_CHECK_CLOSE(s_computed, s[i], reltol);
        }
    }

    // Test the capillary inversion for gas-water.
    {
        const int water = 0;
        const int gas = 2;
        const std::vector<double> pc = { 0.9e5, 0.8e5, 0.6e5, 0.4e5, 0.3e5 };
        const std::vector<double> s = { 0.2, 0.333333333333, 0.6, 0.866666666666, 1.0 };
        BOOST_REQUIRE(pc.size() == s.size());
        for (size_t i = 0; i < pc.size(); ++i) {
            const double s_computed = Opm::Equil::satFromSumOfPcs(props, water, gas, cell, pc[i]);
            BOOST_CHECK_CLOSE(s_computed, s[i], reltol);
        }
    }
}



BOOST_AUTO_TEST_CASE (DeckWithCapillary)
{
    Opm::GridManager gm(1, 1, 20, 1.0, 1.0, 5.0);
    const UnstructuredGrid& grid = *(gm.c_grid());
    Opm::EclipseGridParser deck("capillary.DATA");
    Opm::BlackoilPropertiesFromDeck props(deck, grid, false);

    Opm::Equil::DeckDependent::InitialStateComputer<Opm::EclipseGridParser> comp(props, deck, grid, 10.0);
    const auto& pressures = comp.press();
    BOOST_REQUIRE(pressures.size() == 3);
    BOOST_REQUIRE(int(pressures[0].size()) == grid.number_of_cells);

    const int first = 0, last = grid.number_of_cells - 1;
    // The relative tolerance is too loose to be very useful,
    // but the answer we are checking is the result of an ODE
    // solver, and it is unclear if we should check it against
    // the true answer or something else.
    const double reltol = 1.0e-6;
    BOOST_CHECK_CLOSE(pressures[0][first] , 1.45e7   , reltol);
    BOOST_CHECK_CLOSE(pressures[0][last ] , 1.545e7   , reltol);
    BOOST_CHECK_CLOSE(pressures[1][last] , 1.5351621345e7   , reltol);

    const auto& sats = comp.saturation();
    const std::vector<double> s[3]{
        { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.425893333333, 0.774026666666, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 0, 0, 0, 0.00736, 0.792746666666, 0.8, 0.8, 0.8, 0.8, 0.574106666666, 0.225973333333, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0.8, 0.8, 0.8, 0.79264, 0.007253333333, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }
    };
    for (int phase = 0; phase < 3; ++phase) {
        BOOST_REQUIRE(sats[phase].size() == s[phase].size());
        for (size_t i = 0; i < s[phase].size(); ++i) {
            BOOST_CHECK_CLOSE(sats[phase][i], s[phase][i], reltol);
        }
    }
}



BOOST_AUTO_TEST_CASE (DeckWithCapillaryOverlap)
{
    Opm::GridManager gm(1, 1, 20, 1.0, 1.0, 5.0);
    const UnstructuredGrid& grid = *(gm.c_grid());
    Opm::EclipseGridParser deck("capillary_overlap.DATA");
    Opm::BlackoilPropertiesFromDeck props(deck, grid, false);

    Opm::Equil::DeckDependent::InitialStateComputer<Opm::EclipseGridParser> comp(props, deck, grid, 10.0);
    const auto& pressures = comp.press();
    BOOST_REQUIRE(pressures.size() == 3);
    BOOST_REQUIRE(int(pressures[0].size()) == grid.number_of_cells);

    const int first = 0, last = grid.number_of_cells - 1;
    // The relative tolerance is too loose to be very useful,
    // but the answer we are checking is the result of an ODE
    // solver, and it is unclear if we should check it against
    // the true answer or something else.
    const double reltol = 1.0e-6;
    BOOST_CHECK_CLOSE(pressures[0][first] , 1.45e7   , reltol);
    BOOST_CHECK_CLOSE(pressures[0][last ] , 1.545e7   , reltol);
    BOOST_CHECK_CLOSE(pressures[1][last] , 1.5351621345e7   , reltol);

    const auto& sats = comp.saturation();
    // std::cout << "Saturations:\n";
    // for (const auto& sat : sats) {
    //     for (const double s : sat) {
    //         std::cout << s << ' ';
    //     }
    //     std::cout << std::endl;
    // }
    const std::vector<double> s[3]{
        { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.223141818182, 0.532269090909, 0.78471, 0.91526, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.207743333333, 0.08474, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.776858181818, 0.467730909091, 0.0075466666666, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }
    };
    for (int phase = 0; phase < 3; ++phase) {
        BOOST_REQUIRE(sats[phase].size() == s[phase].size());
        for (size_t i = 0; i < s[phase].size(); ++i) {
            BOOST_CHECK_CLOSE(sats[phase][i], s[phase][i], reltol);
        }
    }
}



BOOST_AUTO_TEST_CASE (DeckWithLiveOil)
{
    Opm::GridManager gm(1, 1, 20, 1.0, 1.0, 5.0);
    const UnstructuredGrid& grid = *(gm.c_grid());
    Opm::EclipseGridParser deck("equil_liveoil.DATA");
    Opm::BlackoilPropertiesFromDeck props(deck, grid, false);

    Opm::Equil::DeckDependent::InitialStateComputer<Opm::EclipseGridParser> comp(props, deck, grid, 10.0);
    const auto& pressures = comp.press();
    BOOST_REQUIRE(pressures.size() == 3);
    BOOST_REQUIRE(int(pressures[0].size()) == grid.number_of_cells);

    const int first = 0, last = grid.number_of_cells - 1;
    // The relative tolerance is too loose to be very useful,
    // but the answer we are checking is the result of an ODE
    // solver, and it is unclear if we should check it against
    // the true answer or something else.
    const double reltol = 1.0e-6;
    BOOST_CHECK_CLOSE(pressures[0][first] , 1.45e7   , reltol);
    BOOST_CHECK_CLOSE(pressures[0][last ] , 1.545e7   , reltol);
    BOOST_CHECK_CLOSE(pressures[1][last] , 1.5489764605846416e7   , reltol);

    const auto& sats = comp.saturation();
    // std::cout << "Saturations:\n";
    // for (const auto& sat : sats) {
    //     for (const double s : sat) {
    //         std::cout << s << ' ';
    //     }
    //     std::cout << std::endl;
    // }
    const std::vector<double> s[3]{
        { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.223141818182, 0.532269090909, 0.78471, 0.91526, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.207743333333, 0.08474, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.776858181818, 0.467730909091, 0.0075466666666, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }
    };
    for (int phase = 0; phase < 3; ++phase) {
        BOOST_REQUIRE(sats[phase].size() == s[phase].size());
        for (size_t i = 0; i < s[phase].size(); ++i) {
            std::cout << sats[phase][i] << '\n';
            //BOOST_CHECK_CLOSE(sats[phase][i], s[phase][i], reltol);
        }
        std::cout << std::endl;
    }
}


BOOST_AUTO_TEST_SUITE_END()
