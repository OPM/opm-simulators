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
#include <opm/parser/eclipse/Units/Units.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

#include <fstream>
#include <iostream>


struct SetupSimple {
    SetupSimple() :
        deck( Opm::Parser{}.parseFile("fluid.data") ),
        eclState( deck, Opm::ParseContext() )
    {
        param.disableOutput();
        param.insertParameter("init_rock"       , "false" );
        param.insertParameter("threephase_model", "simple");
        param.insertParameter("pvt_tab_size"    , "0"     );
        param.insertParameter("sat_tab_size"    , "0"     );
    }

    Opm::parameter::ParameterGroup  param;
    Opm::Deck deck;
    Opm::EclipseState eclState;
};


template <class Setup>
struct TestFixture : public Setup
{
    TestFixture()
        : Setup()
        , grid (eclState.getInputGrid())
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
        , grid (eclState.getInputGrid())
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
    const Opm::BlackoilPropsAdFromDeck::Cells cells(1, 0);

    typedef Opm::BlackoilPropsAdFromDeck::V V;

    enum { Water = Opm::Water };
    V rho0AD_Water = boprops_ad.surfaceDensity(Water, cells);
    BOOST_REQUIRE_EQUAL(rho0AD_Water.size(), cells.size());
    BOOST_CHECK_EQUAL(rho0AD_Water[0], 1000.0);

    enum { Oil = Opm::Oil };
    V rho0AD_Oil = boprops_ad.surfaceDensity(Oil, cells);
    BOOST_REQUIRE_EQUAL(rho0AD_Oil.size(), cells.size());
    BOOST_CHECK_EQUAL(rho0AD_Oil[0], 800.0);

    enum { Gas = Opm::Gas };
    V rho0AD_Gas = boprops_ad.surfaceDensity(Gas, cells);
    BOOST_REQUIRE_EQUAL(rho0AD_Gas.size(), cells.size());
    BOOST_CHECK_EQUAL(rho0AD_Gas[0], 1.0);
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

BOOST_FIXTURE_TEST_CASE(criticalSaturations, TestFixture<SetupSimple>)
{
   const Opm::BlackoilPropsAdFromDeck::Cells cells(10, 0);

    typedef Opm::BlackoilPropsAdFromDeck::V V;

    V sgcr = boprops_ad.scaledCriticalGasSaturations(cells);
    V sogcr = boprops_ad.scaledCriticalOilinGasSaturations(cells);
    BOOST_CHECK_EQUAL(sgcr[0], 0.02);
    BOOST_CHECK_EQUAL(sogcr[0], 0.13);

}
