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
 * \brief Test the ranks of wells in LGRs that never flow
 *
 * A well that never flows is left out of the load balancing and goes to the
 * ranks that hold its cells.  The cells of a well in an LGR are named by
 * their (i,j,k) inside the LGR.  Every connection cell must be on exactly one
 * rank, that rank must have the well, and the well model must set up on all
 * ranks.
 */

#define BOOST_TEST_MODULE LgrWellRanks

#include "SimulatorFixture.hpp"

#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/simulators/flow/FlowProblemBlackoil.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <string>

using SimulatorFixture = Opm::SimulatorFixture;
BOOST_GLOBAL_FIXTURE(SimulatorFixture);

namespace {

// A dipping 12x8x2 grid.  LGR1 refines cells (5..6,4,1) into 6x3x3 cells and LGR2
// cell (9,6,2) into 3x3x3 cells.  Wells PROD and PRL2 are completed in the
// LGRs, NEXT and FAR in the global grid, and all of them are shut.
const char* deck = R"(RUNSPEC
DIMENS
 12 8 2 /
OIL
WATER
METRIC
START
 1 'JAN' 2020 /
WELLDIMS
 4 4 1 4 /
TABDIMS
/
EQLDIMS
/
GRID
CARFIN
'LGR1' 5 6 4 4 1 1 6 3 3 /
ENDFIN
CARFIN
'LGR2' 9 9 6 6 2 2 3 3 3 /
ENDFIN
DXV
 12*100.0 /
DYV
 8*100.0 /
DZV
 2*30.0 /
DEPTHZ
 2000 2004 2008 2012 2016 2020 2024 2028 2032 2036 2040 2044 2048
 2003 2007 2011 2015 2019 2023 2027 2031 2035 2039 2043 2047 2051
 2006 2010 2014 2018 2022 2026 2030 2034 2038 2042 2046 2050 2054
 2009 2013 2017 2021 2025 2029 2033 2037 2041 2045 2049 2053 2057
 2012 2016 2020 2024 2028 2032 2036 2040 2044 2048 2052 2056 2060
 2015 2019 2023 2027 2031 2035 2039 2043 2047 2051 2055 2059 2063
 2018 2022 2026 2030 2034 2038 2042 2046 2050 2054 2058 2062 2066
 2021 2025 2029 2033 2037 2041 2045 2049 2053 2057 2061 2065 2069
 2024 2028 2032 2036 2040 2044 2048 2052 2056 2060 2064 2068 2072 /
PORO
 192*0.25 /
PERMX
 192*500 /
COPY
 PERMX PERMY /
 PERMX PERMZ /
/
PROPS
PVTW
 250.0 1.02 4.0E-05 0.5 0.0 /
ROCK
 250.0 1.0E-03 /
DENSITY
 850 1014 1.0 /
PVDO
 100.0 1.05 2.0
 250.0 1.02 2.5
 400.0 1.00 3.0 /
SWOF
 0.2 0.0  1.0 0.0
 0.3 0.05 0.6 0.0
 0.5 0.15 0.3 0.0
 0.7 0.35 0.1 0.0
 0.8 0.5  0.0 0.0 /
SOLUTION
EQUIL
 2000 250.0 2100 0.0 /
SCHEDULE
WELSPECS
 'NEXT' 'G' 7 4 2020 'OIL' /
 'FAR' 'G' 11 2 2020 'OIL' /
/
WELSPECL
 'PROD' 'G' 'LGR1' 2 2 2020 'OIL' /
 'PRL2' 'G' 'LGR2' 2 2 2050 'OIL' /
/
COMPDAT
 'NEXT' 7 4 1 2 'OPEN' 2* 0.2 /
 'FAR' 11 2 1 2 'OPEN' 2* 0.2 /
/
COMPDATL
 'PROD' 'LGR1' 2 2 1 3 'OPEN' 2* 0.2 /
 'PRL2' 'LGR2' 2 2 2 2 'OPEN' 2* 0.2 /
/
WCONPROD
 'PROD' 'SHUT' 'BHP' 5* 150 /
 'NEXT' 'SHUT' 'BHP' 5* 150 /
 'FAR' 'SHUT' 'BHP' 5* 150 /
 'PRL2' 'SHUT' 'BHP' 5* 150 /
/
TSTEP
 1 /
END
)";

} // Anonymous namespace

BOOST_AUTO_TEST_CASE(ShutWellsInLgrsAreOnTheRanksOfTheirCells)
{
    using TypeTag = Opm::Properties::TTag::TestTypeTag;

    const auto& comm = Opm::FlowGenericVanguard::comm();
    const auto dir = std::filesystem::temp_directory_path() / ("test_lgr_well_ranks_np" + std::to_string(comm.size()));
    const auto deck_file = (dir / "LGR_WELL_RANKS.DATA").string();
    if (comm.rank() == 0) {
        std::filesystem::create_directories(dir);
        std::ofstream os(deck_file);
        os << deck;
    }
    comm.barrier();

    // CARFIN is parsed at low strictness, as in the other runs with local grids.
    // The conforming refinement refines the grid before distributing it.
    auto simulator = Opm::initSimulator<TypeTag>(deck_file.c_str(), "test_lgr_well_ranks",
                                                 /*threads_per_process=*/1, "low",
                                                 { "--lgr-backend=conforming" });
    const auto& vanguard = simulator->vanguard();
    const auto& wells = vanguard.parallelWells();

    for (const std::string name : { "PROD", "PRL2", "NEXT", "FAR" }) {
        BOOST_TEST_CONTEXT("Well " << name) {
            // A serial run distributes no wells: the one rank has them all.
            const auto well = std::ranges::find(wells, name, &std::pair<std::string, bool>::first);
            const bool onRank = (comm.size() == 1) || ((well != wells.end()) && well->second);

            for (const auto& conn : vanguard.schedule().getWell(name, 0).getConnections()) {
                const int cell = vanguard.compressedIndexForConnection(conn);
                BOOST_CHECK_EQUAL(comm.sum(cell >= 0 ? 1 : 0), 1);
                if (cell >= 0) {
                    BOOST_CHECK(onRank);
                }
            }
        }
    }

    simulator->model().applyInitialSolution();
    simulator->model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);
    simulator->setEpisodeIndex(-1);
    simulator->setEpisodeLength(0.0);
    simulator->startNextEpisode(/*episodeStartTime=*/0.0, /*episodeLength=*/1e30);
    simulator->setTimeStepSize(Opm::unit::day);
    simulator->problem().resetIterationForNewTimestep();

    BOOST_CHECK_NO_THROW(simulator->problem().wellModel().beginReportStep(0));

    comm.barrier();
    if (comm.rank() == 0) {
        std::filesystem::remove_all(dir);
    }
}
