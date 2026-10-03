/*
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

#define BOOST_TEST_MODULE LevelZeroRanks
#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <opm/grid/CpGrid.hpp>

#include <opm/simulators/flow/CollectDataOnIORank.hpp>

#include <cmath>
#include <vector>

// A 3x1x1 grid of unit cells whose middle cell is refined into 2x1x1
// children.  Each leaf cell gets 10 plus the I index of the level-0 cell that
// contains it, found from the cell centre, so each level-0 cell, refined or
// not, must end up with its own value.
BOOST_AUTO_TEST_CASE(EachLevelZeroCellTakesTheRankOfItsLeafCells)
{
    Dune::CpGrid grid;
    grid.createCartesian({3, 1, 1}, {1.0, 1.0, 1.0});
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2, 1, 1}},
                               /* startIJK_vec = */ {{1, 0, 0}},
                               /* endIJK_vec = */ {{2, 1, 1}},
                               /* lgr_name_vec = */ {"LGR1"});

    const auto leafView = grid.leafGridView();
    BOOST_REQUIRE_EQUAL(leafView.size(0), 4);

    std::vector<int> leafRanks(leafView.size(0), -1);
    for (const auto& elem : elements(leafView)) {
        const auto i = static_cast<int>(std::floor(elem.geometry().center()[0]));
        leafRanks[leafView.indexSet().index(elem)] = 10 + i;
    }

    const auto ranks = Opm::levelZeroRanks(grid, leafRanks);

    const auto expected = std::vector<int>{10, 11, 12};
    BOOST_CHECK_EQUAL_COLLECTIONS(ranks.begin(), ranks.end(),
                                  expected.begin(), expected.end());
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
