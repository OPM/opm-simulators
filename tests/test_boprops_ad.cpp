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

#include <opm/autodiff/BlackoilPropsAd.hpp>

#include <boost/test/unit_test.hpp>

#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/props/BlackoilPropertiesFromDeck.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <fstream>

struct SetupSimple {
    SetupSimple()
        : param()
        , deck()
    {
        std::ifstream str("fluid.data");
        deck.read(str);

        param.disableOutput();
        param.insertParameter("init_rock"       , "false" );
        param.insertParameter("threephase_model", "simple");
        param.insertParameter("pvt_tab_size"    , "0"     );
        param.insertParameter("sat_tab_size"    , "0"     );
    }

    Opm::parameter::ParameterGroup  param;
    Opm::EclipseGridParser          deck;
};


template <class Setup>
struct TestFixture : public Setup
{
    TestFixture()
        : Setup()
        , grid (deck)
        , props(deck, *grid.c_grid(), param,
                param.getDefault("init_rock", false))
    {
    }

    using Setup::param;
    using Setup::deck;

    Opm::GridManager                grid;
    Opm::BlackoilPropertiesFromDeck props;
};


BOOST_FIXTURE_TEST_CASE(Construction, TestFixture<SetupSimple>)
{
    Opm::BlackoilPropsAd boprops_ad(props);
}



BOOST_FIXTURE_TEST_CASE(SurfaceDensity, TestFixture<SetupSimple>)
{
    Opm::BlackoilPropsAd boprops_ad(props);

    const double* rho0   = props     .surfaceDensity();
    const double* rho0AD = boprops_ad.surfaceDensity();

    enum { Water = Opm::BlackoilPropsAd::Water };
    BOOST_CHECK_EQUAL(rho0AD[ Water ], rho0[ Water ]);

    enum { Oil = Opm::BlackoilPropsAd::Oil };
    BOOST_CHECK_EQUAL(rho0AD[ Oil ], rho0[ Oil ]);

    enum { Gas = Opm::BlackoilPropsAd::Gas };
    BOOST_CHECK_EQUAL(rho0AD[ Gas ], rho0[ Gas ]);
}


BOOST_FIXTURE_TEST_CASE(ViscosityValue, TestFixture<SetupSimple>)
{
    Opm::BlackoilPropsAd boprops_ad(props);

    const Opm::BlackoilPropsAd::Cells cells(5, 0);

    typedef Opm::BlackoilPropsAd::V V;

    V Vpw;
    Vpw.resize(cells.size());
    Vpw[0] =  1*Opm::unit::barsa;
    Vpw[1] =  2*Opm::unit::barsa;
    Vpw[2] =  4*Opm::unit::barsa;
    Vpw[3] =  8*Opm::unit::barsa;
    Vpw[4] = 16*Opm::unit::barsa;

    const Opm::BlackoilPropsAd::V VmuWat = boprops_ad.muWat(Vpw, cells);

    // Zero pressure dependence in water viscosity
    for (V::Index i = 0, n = VmuWat.size(); i < n; ++i) {
        BOOST_CHECK_EQUAL(VmuWat[i], VmuWat[0]);
    }
}


BOOST_FIXTURE_TEST_CASE(ViscosityAD, TestFixture<SetupSimple>)
{
    Opm::BlackoilPropsAd boprops_ad(props);

    const Opm::BlackoilPropsAd::Cells cells(5, 0);

    typedef Opm::BlackoilPropsAd::V V;

    V Vpw;
    Vpw.resize(cells.size());
    Vpw[0] =  1*Opm::unit::barsa;
    Vpw[1] =  2*Opm::unit::barsa;
    Vpw[2] =  4*Opm::unit::barsa;
    Vpw[3] =  8*Opm::unit::barsa;
    Vpw[4] = 16*Opm::unit::barsa;

    typedef Opm::BlackoilPropsAd::ADB ADB;

    const V VmuWat = boprops_ad.muWat(Vpw, cells);
    for (V::Index i = 0, n = Vpw.size(); i < n; ++i) {
        const std::vector<int> bp(1, grid.c_grid()->number_of_cells);

        const Opm::BlackoilPropsAd::Cells c(1, 0);
        const V   pw     = V(1, 1) * Vpw[i];
        const ADB Apw    = ADB::variable(0, pw, bp);
        const ADB AmuWat = boprops_ad.muWat(Apw, c);

        BOOST_CHECK_EQUAL(AmuWat.value()[0], VmuWat[i]);
    }
}
