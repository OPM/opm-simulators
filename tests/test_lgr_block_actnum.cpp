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
 * \brief Test that the ACTNUM of a CARFIN block removes refined cells from
 *        the simulation grid.
 *
 * The conforming refinement leaves only the active refined cells on the leaf
 * grid, so only they connect to the host cells.  The trilinear refinement
 * cannot remove refined cells and stops.
 */

#define BOOST_TEST_MODULE LgrBlockActnum

#include "SimulatorFixture.hpp"

#include <opm/simulators/flow/FlowProblemBlackoil.hpp>

#include <boost/test/unit_test.hpp>

#include <dune/grid/common/rangegenerators.hh>

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

using SimulatorFixture = Opm::SimulatorFixture;
BOOST_GLOBAL_FIXTURE(SimulatorFixture);

namespace {

// A 4x4x1 grid.  LGR1 refines cell (2,2,1) into 3x3x3 cells, and its ACTNUM
// keeps the 9 cells of the middle layer, between depths 2003.33 and 2006.67.
const char* deck = R"(RUNSPEC
DIMENS
 4 4 1 /
OIL
WATER
METRIC
START
 1 'JAN' 2020 /
EQLDIMS
/
TABDIMS
/
GRID
CARFIN
'LGR1' 2 2 2 2 1 1 3 3 3 /
ACTNUM
 9*0 9*1 9*0 /
ENDFIN
DX
 16*100 /
DY
 16*100 /
DZ
 16*10 /
TOPS
 16*2000 /
PORO
 16*0.25 /
PERMX
 16*100 /
PERMY
 16*100 /
PERMZ
 16*10 /
PROPS
SWOF
 0.2 0.0 1.0 0
 1.0 1.0 0.0 0 /
PVTW
 200 1.0 4e-5 0.5 0 /
PVDO
 100 1.05 1.0
 300 1.03 1.0 /
DENSITY
 800 1000 1 /
ROCK
 200 1e-5 /
SOLUTION
EQUIL
 2000 200 2100 0 1000 0 /
SCHEDULE
TSTEP
 1 /
END
)";

std::string writeDeck(const std::filesystem::path& dir)
{
    std::filesystem::create_directories(dir);
    const auto deck_file = (dir / "LGR_BLOCK_ACTNUM.DATA").string();
    std::ofstream os(deck_file);
    os << deck;
    return deck_file;
}

} // Anonymous namespace

BOOST_AUTO_TEST_CASE(BlockActnumRemovesRefinedCells)
{
    using TypeTag = Opm::Properties::TTag::TestTypeTag;

    const auto dir = std::filesystem::temp_directory_path() / "test_lgr_block_actnum";
    const auto deck_file = writeDeck(dir);

    auto simulator = Opm::initSimulator<TypeTag>(deck_file.c_str(), "test_lgr_block_actnum",
                                                 /*threads_per_process=*/1, "low",
                                                 { "--output-dir=" + dir.string() });
    const auto& gridView = simulator->vanguard().grid().leafGridView();

    // The 15 cells outside the box and the 9 active refined cells.
    BOOST_CHECK_EQUAL(gridView.size(0), 24);

    int refined = 0;
    int toHost = 0;
    for (const auto& element : elements(gridView)) {
        if (element.level() == 0) {
            continue;
        }
        ++refined;
        BOOST_CHECK_CLOSE(element.geometry().center()[2], 2005.0, 1.0e-8);
        for (const auto& intersection : intersections(gridView, element)) {
            if (intersection.neighbor() && (intersection.outside().level() == 0)) {
                ++toHost;
            }
        }
    }
    BOOST_CHECK_EQUAL(refined, 9);

    // The 3 cells along each side of the middle layer touch a host cell.
    BOOST_CHECK_EQUAL(toHost, 12);

    std::filesystem::remove_all(dir);
}

BOOST_AUTO_TEST_CASE(TrilinearRefinementStopsOnBlockActnum)
{
    using TypeTag = Opm::Properties::TTag::TestTypeTag;

    const auto dir = std::filesystem::temp_directory_path() / "test_lgr_block_actnum_trilinear";
    const auto deck_file = writeDeck(dir);

    BOOST_CHECK_THROW(Opm::initSimulator<TypeTag>(deck_file.c_str(), "test_lgr_block_actnum",
                                                  /*threads_per_process=*/1, "low",
                                                  { "--output-dir=" + dir.string(),
                                                    "--lgr-backend=trilinear" }),
                      std::invalid_argument);

    std::filesystem::remove_all(dir);
}
