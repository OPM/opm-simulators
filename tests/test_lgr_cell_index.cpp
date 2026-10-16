// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \brief Test the cell index of unrefined cells in a grid with an LGR
 *
 * LGR_CELL_INDEX.DATA refines cell (1,1,1) of a 3x1x1 grid into 2x1x1.  The
 * refined cells come first in the leaf grid, so cells (2,1,1) and (3,1,1) have
 * other leaf indices than level-0 indices.  compressedIndexForInterior(), which
 * places well connections and aquifer cells, must return their leaf cell.
 */

#define BOOST_TEST_MODULE LgrCellIndex

#include "SimulatorFixture.hpp"

#include <opm/simulators/flow/FlowProblemBlackoil.hpp>

#include <boost/test/unit_test.hpp>

using SimulatorFixture = Opm::SimulatorFixture;
BOOST_GLOBAL_FIXTURE(SimulatorFixture);

BOOST_AUTO_TEST_CASE(UnrefinedCellsMapToTheirLeafCell)
{
    using TypeTag = Opm::Properties::TTag::TestTypeTag;

    // CARFIN is parsed at low strictness, as in the other runs with local grids.
    auto simulator = Opm::initSimulator<TypeTag>("LGR_CELL_INDEX.DATA", "test_lgr_cell_index",
                                                 /*threads_per_process=*/1, "low");
    const auto& vanguard = simulator->vanguard();
    const auto& gridView = vanguard.gridView();

    BOOST_REQUIRE_EQUAL(vanguard.grid().maxLevel(), 1);
    BOOST_REQUIRE_EQUAL(gridView.size(0), 4);

    for (const int cartesianIdx : {1, 2}) {
        const int leafIdx = vanguard.compressedIndexForInterior(cartesianIdx);
        BOOST_REQUIRE_GE(leafIdx, 0);
        BOOST_CHECK_EQUAL(vanguard.cartesianIndex(leafIdx), static_cast<unsigned>(cartesianIdx));

        for (const auto& elem : elements(gridView)) {
            if (gridView.indexSet().index(elem) == leafIdx) {
                BOOST_CHECK_EQUAL(elem.level(), 0);
            }
        }
    }
}
