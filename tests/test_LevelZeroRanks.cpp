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
#include <opm/grid/cpgrid/CartesianIndexMapper.hpp>

#include <opm/simulators/flow/CollectDataOnIORank.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <string>
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

// A 4x3x3 grid of unit cells distributed by hand: the process of a cell
// follows its I index.  LGR1 refines cells on two processes, LGR2 cells on the
// last one.  The grids are set up as in a parallel run: the I/O rank keeps a
// copy of the grid made before the distribution, which shares its global view,
// and the LGRs are added to the distributed view, then to the global view, and
// the ids of the refined cells are synchronised.  Each leaf cell collected on
// the I/O rank must carry the process of its level-zero origin, and the ranks
// folded onto level zero must give back the chosen distribution.
BOOST_AUTO_TEST_CASE(CollectedRanksGiveBackTheDistributionOfAGridWithLgrs)
{
    Dune::CpGrid grid;
    grid.createCartesian({4, 3, 3}, {1.0, 1.0, 1.0});

    const int numProcs = grid.comm().size();
    if (numProcs == 1) {
        return;
    }

    std::vector<int> parts(grid.size(0));
    for (std::size_t cell = 0; cell < parts.size(); ++cell) {
        parts[cell] = static_cast<int>(cell % 4) * numProcs / 4;
    }

    std::unique_ptr<Dune::CpGrid> equilGrid;
    if (grid.comm().rank() == 0) {
        equilGrid = std::make_unique<Dune::CpGrid>(grid);
    }

    grid.loadBalance(parts, /* ownersFirst = */ true,
                     /* addCornerCells = */ false, /* overlapLayers = */ 1);

    const std::vector<std::array<int,3>> cellsPerDim = {{2, 2, 2}, {3, 3, 3}};
    const std::vector<std::array<int,3>> startIJK = {{1, 0, 0}, {3, 2, 1}};
    const std::vector<std::array<int,3>> endIJK = {{3, 1, 1}, {4, 3, 3}};
    const std::vector<std::string> names = {"LGR1", "LGR2"};
    grid.addLgrsUpdateLeafView(cellsPerDim, startIJK, endIJK, names);
    grid.switchToGlobalView();
    grid.addLgrsUpdateLeafView(cellsPerDim, startIJK, endIJK, names);
    grid.switchToDistributedView();
    grid.syncDistributedGlobalCellIds();

    const Dune::CartesianIndexMapper<Dune::CpGrid> cartMapper(grid);
    std::unique_ptr<Dune::CartesianIndexMapper<Dune::CpGrid>> equilCartMapper;
    if (equilGrid) {
        equilCartMapper = std::make_unique<Dune::CartesianIndexMapper<Dune::CpGrid>>(*equilGrid);
    }

    const Opm::CollectDataOnIORank<Dune::CpGrid, Dune::CpGrid, Dune::CpGrid::LeafGridView>
        collect(grid, equilGrid.get(), grid.leafGridView(), cartMapper, equilCartMapper.get());

    if (collect.isIORank()) {
        BOOST_REQUIRE_EQUAL(equilGrid->maxLevel(), 2);

        const auto leafView = equilGrid->leafGridView();
        const auto& leafRanks = collect.globalRanks();
        BOOST_REQUIRE_EQUAL(leafRanks.size(), static_cast<std::size_t>(leafView.size(0)));
        for (const auto& elem : elements(leafView)) {
            BOOST_CHECK_EQUAL(leafRanks[leafView.indexSet().index(elem)],
                              parts[elem.getOrigin().index()]);
        }

        const auto ranks = Opm::levelZeroRanks(*equilGrid, leafRanks);
        BOOST_CHECK_EQUAL_COLLECTIONS(ranks.begin(), ranks.end(),
                                      parts.begin(), parts.end());
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
