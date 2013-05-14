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

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif

#define BOOST_TEST_MODULE FluidPropertiesTest

#include "BlackoilPropsAd.hpp"

#include <boost/test/unit_test.hpp>

#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/props/BlackoilPropertiesFromDeck.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <fstream>

struct SetupSimple {
    SetupSimple()
        : param()
        , deck()
    {
        std::ifstream str("fluid.data");
        deck.read(str);

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
