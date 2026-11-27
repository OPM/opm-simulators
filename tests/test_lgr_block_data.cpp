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
 * \brief Test where the block summary vectors of a run with an LGR are kept.
 *
 * A B* vector names a cell of the global grid and an LB* vector a cell of an
 * LGR: the B* vector goes to the block data, the LB* vector to the LGR block
 * data only.
 */

#define BOOST_TEST_MODULE LgrBlockData

#include "SimulatorFixture.hpp"

#include <opm/simulators/flow/FlowProblemBlackoil.hpp>

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <string>

using SimulatorFixture = Opm::SimulatorFixture;
BOOST_GLOBAL_FIXTURE(SimulatorFixture);

namespace {

// A 3x3x1 grid.  LGR1 refines cell (2,2,1) into 3x3x1 cells.  The summary asks
// for the pressure of cell (1,1,1) and of cell (2,2,1) of LGR1.
const char* deck = R"(RUNSPEC
DIMENS
 3 3 1 /
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
'LGR1' 2 2 2 2 1 1 3 3 1 /
ENDFIN
DX
 9*100 /
DY
 9*100 /
DZ
 9*10 /
TOPS
 9*2000 /
PORO
 9*0.25 /
PERMX
 9*100 /
PERMY
 9*100 /
PERMZ
 9*10 /
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
SUMMARY
BPR
 1 1 1 /
/
LBPR
 'LGR1' 2 2 1 /
/
SCHEDULE
TSTEP
 1 /
END
)";

} // Anonymous namespace

BOOST_AUTO_TEST_CASE(LgrBlockVectorsAreKeptAsLgrBlockData)
{
    using TypeTag = Opm::Properties::TTag::TestTypeTag;

    const auto dir = std::filesystem::temp_directory_path() / "test_lgr_block_data";
    std::filesystem::create_directories(dir);
    const auto deck_file = (dir / "LGR_BLOCK_DATA.DATA").string();
    {
        std::ofstream os(deck_file);
        os << deck;
    }

    auto simulator = Opm::initSimulator<TypeTag>(deck_file.c_str(), "test_lgr_block_data",
                                                 /*threads_per_process=*/1, "low",
                                                 { "--output-dir=" + dir.string() });
    auto& outputModule = simulator->problem().eclWriter().mutableOutputModule();

    const auto& blockData = outputModule.getBlockData();
    BOOST_CHECK_EQUAL(blockData.count({ "BPR", 1 }), 1U);
    BOOST_CHECK(std::none_of(blockData.begin(), blockData.end(),
                             [](const auto& entry) { return entry.first.first == "LBPR"; }));

    const auto& lgrBlockData = outputModule.getLgrBlockData();
    BOOST_CHECK(std::any_of(lgrBlockData.begin(), lgrBlockData.end(),
                            [](const auto& entry) { return std::get<0>(entry.first) == "LBPR"; }));

    std::filesystem::remove_all(dir);
}
