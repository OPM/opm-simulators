/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.
  Copyright 2013 Statoil ASA.

  This file is part of the Open Porous Media Project (OPM).

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

#define BOOST_TEST_MODULE FluidPropertiesTest

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <boost/test/unit_test.hpp>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <opm/autodiff/BlackoilPropsAdFromDeck.hpp>

#include <opm/core/grid/GridManager.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Parser/ParseMode.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

#include <fstream>
#include <iostream>


struct SetupSimple {
    SetupSimple()
    {
        Opm::ParseMode parseMode;
        Opm::ParserPtr parser(new Opm::Parser());
        deck = parser->parseFile("fluid.data", parseMode);
        eclState.reset(new Opm::EclipseState(deck , parseMode));

        param.disableOutput();
        param.insertParameter("init_rock"       , "false" );
        param.insertParameter("threephase_model", "simple");
        param.insertParameter("pvt_tab_size"    , "0"     );
        param.insertParameter("sat_tab_size"    , "0"     );
    }

    Opm::parameter::ParameterGroup  param;
    Opm::DeckConstPtr               deck;
    Opm::EclipseStateConstPtr       eclState;
};


template <class Setup>
struct TestFixture : public Setup
{
    TestFixture()
        : Setup()
        , grid (deck)
        , boprops_ad(deck, eclState, *grid.c_grid(), param.getDefault("init_rock", false))
    {
    }

    using Setup::param;
    using Setup::deck;
    using Setup::eclState;

    Opm::GridManager             grid;
    Opm::BlackoilPropsAdFromDeck boprops_ad;
};

template <class Setup>
struct TestFixtureAd : public Setup
{
    TestFixtureAd()
        : Setup()
        , grid (deck)
        , props(deck, eclState, *grid.c_grid(),
                param.getDefault("init_rock", false))
    {
    }

    using Setup::param;
    using Setup::deck;
    using Setup::eclState;

    Opm::GridManager             grid;
    Opm::BlackoilPropsAdFromDeck props;
};


BOOST_FIXTURE_TEST_CASE(Construction, TestFixture<SetupSimple>)
{
}

BOOST_FIXTURE_TEST_CASE(SubgridConstruction, TestFixtureAd<SetupSimple>)
{
    Opm::BlackoilPropsAdFromDeck subgrid_props(props);
}

BOOST_FIXTURE_TEST_CASE(SurfaceDensity, TestFixture<SetupSimple>)
{
    const double* rho0AD = boprops_ad.surfaceDensity();

    enum { Water = Opm::BlackoilPropsAdFromDeck::Water };
    BOOST_CHECK_EQUAL(rho0AD[ Water ], 1000.0);

    enum { Oil = Opm::BlackoilPropsAdFromDeck::Oil };
    BOOST_CHECK_EQUAL(rho0AD[ Oil ], 800.0);

    enum { Gas = Opm::BlackoilPropsAdFromDeck::Gas };
    BOOST_CHECK_EQUAL(rho0AD[ Gas ], 1.0);
}


BOOST_FIXTURE_TEST_CASE(ViscosityValue, TestFixture<SetupSimple>)
{
    const Opm::BlackoilPropsAdFromDeck::Cells cells(5, 0);

    typedef Opm::BlackoilPropsAdFromDeck::V V;
    typedef Opm::BlackoilPropsAdFromDeck::ADB ADB;

    V Vpw;
    Vpw.resize(cells.size());
    Vpw[0] =  1*Opm::unit::barsa;
    Vpw[1] =  2*Opm::unit::barsa;
    Vpw[2] =  4*Opm::unit::barsa;
    Vpw[3] =  8*Opm::unit::barsa;
    Vpw[4] = 16*Opm::unit::barsa;

    // standard temperature
    V T = V::Constant(cells.size(), 273.15+20);

    BOOST_REQUIRE_EQUAL(Vpw.size(), cells.size());

    const V VmuWat = boprops_ad.muWat(ADB::constant(Vpw), ADB::constant(T), cells).value();

    BOOST_REQUIRE_EQUAL(Vpw.size(), cells.size());

    // Zero pressure dependence in water viscosity
    for (V::Index i = 0, n = VmuWat.size(); i < n; ++i) {
        BOOST_CHECK_EQUAL(VmuWat[i], VmuWat[0]);
    }
}


BOOST_FIXTURE_TEST_CASE(ViscosityAD, TestFixture<SetupSimple>)
{
    const Opm::BlackoilPropsAdFromDeck::Cells cells(5, 0);

    typedef Opm::BlackoilPropsAdFromDeck::V V;
    typedef Opm::BlackoilPropsAdFromDeck::ADB ADB;

    V Vpw;
    Vpw.resize(cells.size());
    Vpw[0] =  1*Opm::unit::barsa;
    Vpw[1] =  2*Opm::unit::barsa;
    Vpw[2] =  4*Opm::unit::barsa;
    Vpw[3] =  8*Opm::unit::barsa;
    Vpw[4] = 16*Opm::unit::barsa;

    // standard temperature
    V T = V::Constant(cells.size(), 273.15+20);

    typedef Opm::BlackoilPropsAdFromDeck::ADB ADB;

    const V VmuWat = boprops_ad.muWat(ADB::constant(Vpw), ADB::constant(T), cells).value();
    for (V::Index i = 0, n = Vpw.size(); i < n; ++i) {
        const std::vector<int> bp(1, grid.c_grid()->number_of_cells);

        const Opm::BlackoilPropsAdFromDeck::Cells c(1, 0);
        const V   pw     = V(1, 1) * Vpw[i];
        const ADB Apw    = ADB::variable(0, pw, bp);
        const ADB AT     = ADB::constant(T);
        const ADB AmuWat = boprops_ad.muWat(Apw, AT, c);

        BOOST_CHECK_EQUAL(AmuWat.value()[0], VmuWat[i]);
    }
}
