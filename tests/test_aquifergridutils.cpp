/*
  Copyright 2018 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2018 Equinor.

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
#define BOOST_TEST_MODULE AquiferGridUtils
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>

#include <dune/grid/common/rangegenerators.hh>

#include <opm/grid/polyhedralgrid.hh>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <opm/simulators/aquifers/AquiferGridUtils.hpp>

namespace {

Opm::Deck createNumericalAquiferDeck()
{
    const char *deckData = R"(
DIMENS
 8 15 3 /
AQUDIMS
    3      2      1*       1*     1*       50      1*      1*  /
GRID

DX
  360*10./
DY
   360*10./
DZ
   360*1./
TOPS
   360*100./

PORO
   0. 0.25 0. 357*0.25/
PERMX
    360*1000./
PERMY
    360*1000./
PERMZ
    360*10./

BOX
1	8 15 15 3 3 /

MULTY
 1e9 1e-9 1.0 2.0 3.0 4.0 5.0 6.0/

ENDBOX

-- setting the three cells for numerical aquifer to be inactive
ACTNUM
0 1 0 0 356*1 /

AQUNUM
--AQnr.  I  J  K     Area      Length PHI      K     Depth  Initial.Pr	PVTNUM   SATNUM
   1     1  1  1   1000000.0   10000   0.25   400    2585.00   285.00	 2   2  /
   1     3  1  1   1500000.0   20000   0.24   600    2585.00   285.00	 3   *  /
   1     4  1  1   2000000.0   30000   *   700    2585.00   285.00	 *   3  /
/
AQUCON
--  Connect numerical aquifer to the reservoir
--  Id.nr  I1	I2     J1  J2	 K1  K2    Face    Trans.mult.  Trans.opt.
     1     1	8      15    15	  3   3   'J+'      1.00      1  /
/
    )";

    Opm::Parser parser;
    return parser.parseString(deckData);
}

} // Anonymous namespace

struct Fixture {
    Fixture()
        : numaquifer_deck(createNumericalAquiferDeck())
        , ecl_state(numaquifer_deck)
        , ecl_grid(ecl_state.getInputGrid())
    {
    }

    Opm::Deck numaquifer_deck;
    Opm::EclipseState ecl_state;
    const Opm::EclipseGrid& ecl_grid;
};

BOOST_FIXTURE_TEST_CASE(NumericalAquiferCellUnsupported, Fixture)
{
    Dune::PolyhedralGrid<3,3,double> grid(ecl_grid);

    Opm::IsNumericalAquiferCell isNumericalAquiferCell(grid);

    BOOST_CHECK_EQUAL(isNumericalAquiferCell(1), false);
}

BOOST_FIXTURE_TEST_CASE(NumericalAquiferCpGrid, Fixture)
{
    Dune::CpGrid grid;
    grid.processEclipseFormat(&ecl_grid, &ecl_state, false, false, false);

    Opm::IsNumericalAquiferCell isNumericalAquiferCell(grid);

    for (const auto& elem : elements(grid.leafGridView())) {
        BOOST_CHECK_EQUAL(isNumericalAquiferCell(elem), elem.index() == 0 ||
                                                        elem.index() == 2 ||
                                                        elem.index() == 3);
    }
}

bool init_unit_test_func()
{
    return true;
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
