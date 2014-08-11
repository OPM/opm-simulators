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

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

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
            { 0 , 0   } , // Zgoc       , Pcgo_goc
            0, 0, 0
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
                { 0   ,  0       } , // Zgoc       , Pcgo_goc
                0, 0, 0
            }
            ,
            {
                { 5   ,  1.35e5  } , // Datum depth, pressure
                { 7.5 , -0.225e5 } , // Zwoc       , Pcow_woc
                { 5   ,  0       } , // Zgoc       , Pcgo_goc
                0, 0, 0
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
                { 0   ,  0       } , // Zgoc       , Pcgo_goc
                0, 0, 0
            }
            ,
            {
                { 5   ,  1.35e5  } , // Datum depth, pressure
                { 7.5 , -0.225e5 } , // Zwoc       , Pcow_woc
                { 5   ,  0       } , // Zgoc       , Pcgo_goc
                0, 0, 0
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
    Opm::ParserPtr parser(new Opm::Parser() );
    Opm::DeckConstPtr deck = parser->parseFile("deadfluids.DATA");
    Opm::EclipseStateConstPtr eclipseState(new Opm::EclipseState(deck));
    Opm::BlackoilPropertiesFromDeck props(deck, eclipseState, *grid, false);
    Opm::Equil::DeckDependent::InitialStateComputer comp(props, deck, eclipseState, *grid, 10.0);
    const auto& pressures = comp.press();
    BOOST_REQUIRE(pressures.size() == 3);
    BOOST_REQUIRE(int(pressures[0].size()) == grid->number_of_cells);

    const int first = 0, last = grid->number_of_cells - 1;
    // The relative tolerance is too loose to be very useful,
    // but the answer we are checking is the result of an ODE
    // solver, and it is unclear if we should check it against
    // the true answer or something else.
    const double reltol = 1.0e-3;
    BOOST_CHECK_CLOSE(pressures[0][first] , 1.496329839e7   , reltol);
    BOOST_CHECK_CLOSE(pressures[0][last ] , 1.50473245e7   , reltol);
    BOOST_CHECK_CLOSE(pressures[1][last] , 1.50473245e7   , reltol);
}



BOOST_AUTO_TEST_CASE (CapillaryInversion)
{
    // Test setup.
    Opm::GridManager gm(1, 1, 40, 1.0, 1.0, 2.5);
    const UnstructuredGrid& grid = *(gm.c_grid());
    Opm::ParserPtr parser(new Opm::Parser() );
    Opm::DeckConstPtr deck = parser->parseFile("capillary.DATA");
    Opm::EclipseStateConstPtr eclipseState(new Opm::EclipseState(deck));
    Opm::BlackoilPropertiesFromDeck props(deck, eclipseState, grid, false);

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
    Opm::ParserPtr parser(new Opm::Parser() );
    Opm::DeckConstPtr deck = parser->parseFile("capillary.DATA");
    Opm::EclipseStateConstPtr eclipseState(new Opm::EclipseState(deck));
    Opm::BlackoilPropertiesFromDeck props(deck, eclipseState, grid, false);

    Opm::Equil::DeckDependent::InitialStateComputer comp(props, deck, eclipseState, grid, 10.0);
    const auto& pressures = comp.press();
    BOOST_REQUIRE(pressures.size() == 3);
    BOOST_REQUIRE(int(pressures[0].size()) == grid.number_of_cells);

    const int first = 0, last = grid.number_of_cells - 1;
    // The relative tolerance is too loose to be very useful,
    // but the answer we are checking is the result of an ODE
    // solver, and it is unclear if we should check it against
    // the true answer or something else.
    const double reltol = 1.0e-6;
    BOOST_CHECK_CLOSE(pressures[0][first] , 1.469769063e7   , reltol);
    BOOST_CHECK_CLOSE(pressures[0][last ] , 15452880.328284413   , reltol);
    BOOST_CHECK_CLOSE(pressures[1][last] , 15462880.328284413   , reltol);

    const auto& sats = comp.saturation();
    const std::vector<double> s[3]{
        { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.42192000000000002, 0.77802666666666664, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 0, 0, 0, 0.00736, 0.792746666666, 0.8, 0.8, 0.8, 0.8, 0.57807999999999993, 0.22197333333333336, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
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
    Opm::ParserPtr parser(new Opm::Parser() );
    Opm::DeckConstPtr deck = parser->parseFile("capillary_overlap.DATA");
    Opm::EclipseStateConstPtr eclipseState(new Opm::EclipseState(deck));
    Opm::BlackoilPropertiesFromDeck props(deck, eclipseState, grid, false);

    Opm::Equil::DeckDependent::InitialStateComputer comp(props, deck, eclipseState, grid, 9.80665);
    const auto& pressures = comp.press();
    BOOST_REQUIRE(pressures.size() == 3);
    BOOST_REQUIRE(int(pressures[0].size()) == grid.number_of_cells);

    const int first = 0, last = grid.number_of_cells - 1;
    // The relative tolerance is too loose to be very useful,
    // but the answer we are checking is the result of an ODE
    // solver, and it is unclear if we should check it against
    // the true answer or something else.
    const double reltol = 1.0e-6;
    const double reltol_ecl = 1.0;
    BOOST_CHECK_CLOSE(pressures[0][first], 1.48324e+07, reltol_ecl);  // eclipse 
    BOOST_CHECK_CLOSE(pressures[0][last],  1.54801e+07, reltol_ecl);
    BOOST_CHECK_CLOSE(pressures[1][first], 1.49224e+07, reltol_ecl);
    BOOST_CHECK_CLOSE(pressures[1][last],  1.54901e+07, reltol_ecl);
    
    BOOST_CHECK_CLOSE(pressures[0][first] , 14832467.14, reltol); // opm
    BOOST_CHECK_CLOSE(pressures[0][last ] , 15479883.47, reltol);
    BOOST_CHECK_CLOSE(pressures[1][last ] , 15489883.47, reltol);

    const auto& sats = comp.saturation();
    // std::cout << "Saturations:\n";
    // for (const auto& sat : sats) {
    //     for (const double s : sat) {
    //         std::cout << s << ' ';
    //     }
    //     std::cout << std::endl;
    // }
   
    const std::vector<double> s_ecl[3]{// eclipse
        { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.22874042, 0.53397995, 0.78454906,  0.91542006, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 0,   0,   0,   0,   0,   0,   0,   0,          0,          0.20039,     0.08458,    0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.77125955, 0.46602005, 0.015063271, 0,          0, 0, 0, 0, 0, 0, 0, 0, 0 }
    };

    const std::vector<double> s_opm[3]{ // opm 
        { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2289309090909091,  0.53406545454545451, 0.78458,              0.9154, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 0,   0,   0,   0,   0,   0,   0,   0,                   0,                   0.2002466666666666,   0.0846, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.77106909090909093, 0.46593454545454549, 0.015173333333333336, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0 }
    };
    for (int phase = 0; phase < 3; ++phase) {
        BOOST_REQUIRE(sats[phase].size() == s_opm[phase].size());
        for (size_t i = 0; i < s_opm[phase].size(); ++i) {
            //std::cout << std::setprecision(10) << sats[phase][i] << '\n';
            BOOST_CHECK_CLOSE(sats[phase][i], s_ecl[phase][i], reltol_ecl);
            BOOST_CHECK_CLOSE(sats[phase][i], s_opm[phase][i], reltol);
        }
    }
}



BOOST_AUTO_TEST_CASE (DeckWithLiveOil)
{
    Opm::GridManager gm(1, 1, 20, 1.0, 1.0, 5.0);
    const UnstructuredGrid& grid = *(gm.c_grid());
    Opm::ParserPtr parser(new Opm::Parser() );
    Opm::DeckConstPtr deck = parser->parseFile("equil_liveoil.DATA");
    Opm::EclipseStateConstPtr eclipseState(new Opm::EclipseState(deck));
    Opm::BlackoilPropertiesFromDeck props(deck, eclipseState, grid, false);

    Opm::Equil::DeckDependent::InitialStateComputer comp(props, deck, eclipseState, grid, 9.80665);
    const auto& pressures = comp.press();
    BOOST_REQUIRE(pressures.size() == 3);
    BOOST_REQUIRE(int(pressures[0].size()) == grid.number_of_cells);

    const int first = 0, last = grid.number_of_cells - 1;
    // The relative tolerance is too loose to be very useful,
    // but the answer we are checking is the result of an ODE
    // solver, and it is unclear if we should check it against
    // the true answer or something else.
    const double reltol = 1.0e-6;
    const double reltol_ecl = 1.0;
    BOOST_CHECK_CLOSE(pressures[0][first], 1.48324e+07, reltol_ecl);  // eclipse
    BOOST_CHECK_CLOSE(pressures[0][last],  1.54801e+07, reltol_ecl);
    BOOST_CHECK_CLOSE(pressures[1][first], 1.49224e+07, reltol_ecl);
    BOOST_CHECK_CLOSE(pressures[1][last],  1.54901e+07, reltol_ecl);
    
    BOOST_CHECK_CLOSE(pressures[0][first], 1.483246714e7, reltol);  // opm
    BOOST_CHECK_CLOSE(pressures[0][last],  1.547991652e7, reltol);
    BOOST_CHECK_CLOSE(pressures[1][first], 1.492246714e7, reltol);
    BOOST_CHECK_CLOSE(pressures[1][last],  1.548991652e7, reltol);

    const auto& sats = comp.saturation();
    // std::cout << "Saturations:\n";
    // for (const auto& sat : sats) {
    //     for (const double s : sat) {
    //         std::cout << s << ' ';
    //     }
    //     std::cout << std::endl;
    // }
    const std::vector<double> s_ecl[3]{ // eclipse
        { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.22898, 0.53422, 0.78470, 0.91531, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 0,   0,   0,   0,   0,   0,   0,   0,       0,       0.20073, 0.08469, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.77102, 0.46578, 0.01458, 0,       0, 0, 0, 0, 0, 0, 0, 0, 0 }
    };
    const std::vector<double> s_opm[3]{ // opm 
        { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2291709091, 0.5343054545, 0.78472,        0.91529, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 0,   0,   0,   0,   0,   0,   0,   0,            0,            0.2005866667,   0.08471, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.7708290909, 0.4656945455, 0.01469333333,  0,       0, 0, 0, 0, 0, 0, 0, 0, 0 }
    };
    for (int phase = 0; phase < 3; ++phase) {
        BOOST_REQUIRE(sats[phase].size() == s_opm[phase].size());
        for (size_t i = 0; i < s_opm[phase].size(); ++i) {
            //std::cout << std::setprecision(10) << sats[phase][i] << '\n';
            BOOST_CHECK_CLOSE(sats[phase][i], s_opm[phase][i], reltol);
            BOOST_CHECK_CLOSE(sats[phase][i], s_ecl[phase][i], reltol_ecl);
        }
        std::cout << std::endl;
    }
    
    const auto& rs = comp.rs();
    const std::vector<double> rs_opm {74.61233568, 74.64905212, 74.68578656, 74.72253902, // opm
                                      74.75930951, 74.79609803, 74.83290459, 74.87519876,
                                      74.96925416, 75.09067512, 75.0,        75.0, 
                                      75.0,        75.0,        75.0,        75.0, 
                                      75.0,        75.0,        75.0,        75.0};
    const std::vector<double> rs_ecl {74.612228, 74.648956, 74.685707, 74.722473,  // eclipse
                                      74.759254, 74.796051, 74.832870, 74.875145,
                                      74.969231, 75.090706, 75.000000, 75.000000,
                                      75.000000, 75.000000, 75.000000, 75.000000,
                                      75.000000, 75.000000, 75.000000, 75.000000}; 
    for (size_t i = 0; i < rs_opm.size(); ++i) {
        //std::cout << std::setprecision(10) << rs[i] << '\n';
        BOOST_CHECK_CLOSE(rs[i], rs_opm[i], reltol);
        BOOST_CHECK_CLOSE(rs[i], rs_ecl[i], reltol_ecl);
    }
}



BOOST_AUTO_TEST_CASE (DeckWithLiveGas)
{
    Opm::GridManager gm(1, 1, 20, 1.0, 1.0, 5.0);
    const UnstructuredGrid& grid = *(gm.c_grid());
    Opm::ParserPtr parser(new Opm::Parser() );
    Opm::DeckConstPtr deck = parser->parseFile("equil_livegas.DATA");
    Opm::EclipseStateConstPtr eclipseState(new Opm::EclipseState(deck));
    Opm::BlackoilPropertiesFromDeck props(deck, eclipseState, grid, false);

    Opm::Equil::DeckDependent::InitialStateComputer comp(props, deck, eclipseState, grid, 9.80665);
    const auto& pressures = comp.press();
    BOOST_REQUIRE(pressures.size() == 3);
    BOOST_REQUIRE(int(pressures[0].size()) == grid.number_of_cells);

    const int first = 0, last = grid.number_of_cells - 1;
    // The relative tolerance is too loose to be very useful,
    // but the answer we are checking is the result of an ODE
    // solver, and it is unclear if we should check it against
    // the true answer or something else.
    const double reltol = 1.0e-6;
    const double reltol_ecl = 1.0;
    BOOST_CHECK_CLOSE(pressures[0][first], 1.48215e+07, reltol_ecl);  // eclipse
    BOOST_CHECK_CLOSE(pressures[0][last],  1.54801e+07, reltol_ecl);
    BOOST_CHECK_CLOSE(pressures[1][first], 1.49115e+07, reltol_ecl);
    BOOST_CHECK_CLOSE(pressures[1][last],  1.54901e+07, reltol_ecl);
    
    BOOST_CHECK_CLOSE(pressures[0][first], 1.482150311e7, reltol);  // opm
    BOOST_CHECK_CLOSE(pressures[0][last],  1.547988347e7, reltol);
    BOOST_CHECK_CLOSE(pressures[1][first], 1.491150311e7, reltol);
    BOOST_CHECK_CLOSE(pressures[1][last],  1.548988347e7, reltol);

    const auto& sats = comp.saturation();
    // std::cout << "Saturations:\n";
    // for (const auto& sat : sats) {
    //     for (const double s : sat) {
    //         std::cout << s << ' ';
    //     }
    //     std::cout << std::endl;
    // }
    const std::vector<double> s_ecl[3]{ // eclipse
        { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.24285614, 0.53869015, 0.78454906,  0.91542006, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 0,   0,   0,   0,   0,   0,   0,   0,          0,          0.18311,     0.08458,    0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.75714386, 0.46130988, 0.032345835, 0,          0, 0, 0, 0, 0, 0, 0, 0, 0 }
    };

    const std::vector<double> s_opm[3]{ // opm 
        { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.24310545, 0.5388, 0.78458,    0.91540, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 0,   0,   0,   0,   0,   0,   0,   0,          0,      0.18288667, 0.0846,  0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.75689455, 0.4612, 0.03253333, 0,       0, 0, 0, 0, 0, 0, 0, 0, 0 }
    };
    for (int phase = 0; phase < 3; ++phase) {
        BOOST_REQUIRE(sats[phase].size() == s_opm[phase].size());
        for (size_t i = 0; i < s_opm[phase].size(); ++i) {
            //std::cout << std::setprecision(10) << sats[phase][i] << '\n';
            BOOST_CHECK_CLOSE(sats[phase][i], s_opm[phase][i], 100.*reltol);
            BOOST_CHECK_CLOSE(sats[phase][i], s_ecl[phase][i], reltol_ecl);
        }
        std::cout << std::endl;
    }
    
    const auto& rv = comp.rv();
    const std::vector<double> rv_opm { // opm
        2.4884509e-4, 2.4910378e-4, 2.4936267e-4, 2.4962174e-4,
        2.4988100e-4, 2.5014044e-4, 2.5040008e-4, 2.5065990e-4, 
        2.5091992e-4, 2.5118012e-4, 2.5223082e-4, 2.5105e-4, 
        2.5105e-4,    2.5105e-4,    2.5105e-4,    2.5105e-4, 
        2.5105e-4,    2.5105e-4,    2.5105e-4,    2.5105e-4};

    const std::vector<double> rv_ecl {  // eclipse
        0.24884584E-03,   0.24910446E-03,   0.24936325E-03,   0.24962222E-03,
        0.24988138E-03,   0.25014076E-03,   0.25040031E-03,   0.25066003E-03,
        0.25091995E-03,   0.25118008E-03,   0.25223137E-03,   0.25104999E-03,
        0.25104999E-03,   0.25104999E-03,   0.25104999E-03,   0.25104999E-03,
        0.25104999E-03,   0.25104999E-03,   0.25104999E-03,   0.25104999E-03};
         
    for (size_t i = 0; i < rv_opm.size(); ++i) {
        //std::cout << std::setprecision(10) << sats[phase][i] << '\n';
        BOOST_CHECK_CLOSE(rv[i], rv_opm[i], 100.*reltol);
        BOOST_CHECK_CLOSE(rv[i], rv_ecl[i], reltol_ecl);
    }
}

BOOST_AUTO_TEST_CASE (DeckWithRSVDAndRVVD)
{
    Opm::GridManager gm(1, 1, 20, 1.0, 1.0, 5.0);
    const UnstructuredGrid& grid = *(gm.c_grid());
    Opm::ParserPtr parser(new Opm::Parser() );
    Opm::DeckConstPtr deck = parser->parseFile("equil_rsvd_and_rvvd.DATA");
    Opm::EclipseStateConstPtr eclipseState(new Opm::EclipseState(deck));
    Opm::BlackoilPropertiesFromDeck props(deck, eclipseState, grid, false);

    Opm::Equil::DeckDependent::InitialStateComputer comp(props, deck, eclipseState, grid, 9.80665);
    const auto& pressures = comp.press();
    BOOST_REQUIRE(pressures.size() == 3);
    BOOST_REQUIRE(int(pressures[0].size()) == grid.number_of_cells);

    const int first = 0, last = grid.number_of_cells - 1;
    // The relative tolerance is too loose to be very useful,
    // but the answer we are checking is the result of an ODE
    // solver, and it is unclear if we should check it against
    // the true answer or something else.
    const double reltol = 1.0e-6;
    const double reltol_ecl = 1.0;
    BOOST_CHECK_CLOSE(pressures[0][first], 1.48350e+07, reltol_ecl);  // eclipse
    BOOST_CHECK_CLOSE(pressures[0][last],  1.54794e+07, reltol_ecl);
    BOOST_CHECK_CLOSE(pressures[1][first], 1.49250e+07, reltol_ecl);
    BOOST_CHECK_CLOSE(pressures[1][last],  1.54894e+07, reltol_ecl);
    
    BOOST_CHECK_CLOSE(pressures[0][first], 1.483499660e7, reltol);  // opm
    BOOST_CHECK_CLOSE(pressures[0][last],  1.547924516e7, reltol);
    BOOST_CHECK_CLOSE(pressures[1][first], 1.492499660e7, reltol);
    BOOST_CHECK_CLOSE(pressures[1][last],  1.548924516e7, reltol);

    const auto& sats = comp.saturation();
    // std::cout << "Saturations:\n";
    // for (const auto& sat : sats) {
    //     for (const double s : sat) {
    //         std::cout << s << ' ';
    //     }
    //     std::cout << std::endl;
    // }
    const std::vector<double> s_ecl[3]{ // eclipse
        { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.22206347, 0.52871972, 0.78150368,  0.91819441,  1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 0,   0,   0,   0,   0,   0,   0,   0,          0,          0.19656529,  0.081805572, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.77793652, 0.47128031, 0.021931054, 0,           0, 0, 0, 0, 0, 0, 0, 0, 0 }
    };

    const std::vector<double> s_opm[3]{ // opm 
        { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.22232000, 0.52882909, 0.78153000,    0.91817000, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 0,   0,   0,   0,   0,   0,   0,   0,          0,          0.19636333,    0.08183000, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.77768000, 0.47117091, 0.02210667,    0,          0, 0, 0, 0, 0, 0, 0, 0, 0 }
    };
    
    for (int phase = 0; phase < 3; ++phase) {
        BOOST_REQUIRE(sats[phase].size() == s_opm[phase].size());
        for (size_t i = 0; i < s_opm[phase].size(); ++i) {
            //std::cout << std::setprecision(10) << sats[phase][i] << '\n';
            BOOST_CHECK_CLOSE(sats[phase][i], s_opm[phase][i], 100.*reltol);
            BOOST_CHECK_CLOSE(sats[phase][i], s_ecl[phase][i], reltol_ecl);
        }
        std::cout << std::endl;
    }
    
    const auto& rs = comp.rs();
    const std::vector<double> rs_opm { // opm
        74.62498302, 74.65959041, 74.69438035, 74.72935336,
        74.76450995, 74.79985061, 74.83537588, 74.87527125,
        74.96863769, 75.08891765, 52.5,        57.5,
        62.5,        67.5,        72.5,        76.45954841,
        76.70621045, 76.95287736, 77.19954913, 77.44622578};

    const std::vector<double> rs_ecl {  // eclipse
        74.625114, 74.659706, 74.694481, 74.729439,
        74.764580, 74.799904, 74.835419, 74.875252,
        74.968628, 75.088951, 52.500000, 57.500000,
        62.500000, 67.500000, 72.500000, 76.168388,
        76.349953, 76.531532, 76.713142, 76.894775,};
    
    const auto& rv = comp.rv();
    const std::vector<double> rv_opm { // opm
        2.50e-6, 7.50e-6,       1.25e-5,       1.75e-5,
        2.25e-5, 2.75e-5,       3.25e-5,       3.75e-5, 
        4.25e-5, 2.51158386e-4, 2.52203372e-4, 5.75e-5, 
        6.25e-5, 6.75e-5,       7.25e-5,       7.75e-5, 
        8.25e-5, 8.75e-5,       9.25e-5,       9.75e-5};

    const std::vector<double> rv_ecl {  // eclipse
        0.24999999E-05, 0.74999998E-05, 0.12500000E-04, 0.17500000E-04,
        0.22500000E-04, 0.27500000E-04, 0.32500000E-04, 0.37500002E-04,
        0.42500000E-04, 0.25115837E-03, 0.25220393E-03, 0.57500001E-04,
        0.62500003E-04, 0.67499997E-04, 0.72499999E-04, 0.77500001E-04,
        0.82500002E-04, 0.87499997E-04, 0.92499999E-04, 0.97500000E-04};
         
    for (size_t i = 0; i < rv_opm.size(); ++i) {
        //std::cout << std::setprecision(10) << rs[i] << '\n';
        BOOST_CHECK_CLOSE(rs[i], rs_opm[i], 100*reltol);
        BOOST_CHECK_CLOSE(rs[i], rs_ecl[i], reltol_ecl);
        BOOST_CHECK_CLOSE(rv[i], rv_opm[i], 100.*reltol);
        BOOST_CHECK_CLOSE(rv[i], rv_ecl[i], reltol_ecl);
    }
}

BOOST_AUTO_TEST_SUITE_END()
