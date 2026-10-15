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
 * \brief Test the INIT output of a parallel run with LGRs along a sealed fault.
 *
 * The partitioner leaves a cell out of a process' overlap layer when the
 * transmissibility of the face to it is zero, so the refined cells on the two
 * sides of a sealed fault are never neighbours on one process.  The INIT output
 * must still be extracted, with their connection sealed.
 */

#define BOOST_TEST_MODULE LgrSealedFaultRanks

#include "SimulatorFixture.hpp"

#include <opm/simulators/flow/FlowProblemBlackoil.hpp>

#include <boost/test/unit_test.hpp>

#include <filesystem>
#include <fstream>
#include <string>

using SimulatorFixture = Opm::SimulatorFixture;
BOOST_GLOBAL_FIXTURE(SimulatorFixture);

namespace {

// A 10x10x3 grid with a sealed fault on the X face of column (4,6..7), and
// four LGRs along it.
const char* deck = R"(RUNSPEC
DIMENS
 10 10 3 /
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
DX
 300*100 /
DY
 300*100 /
DZ
 100*10 100*10 100*10 /
TOPS
 100*2000 /
PORO
 300*0.25 /
PERMX
 300*100 /
PERMY
 300*100 /
PERMZ
 300*10 /
FAULTS
 'F1' 4 4 6 7 1 3 'X' /
/
MULTFLT
 'F1' 0 /
/
CARFIN
'LGR0' 4 4 6 7 1 3 3 6 3 /
ENDFIN
CARFIN
'LGR1' 5 5 5 6 1 3 3 6 3 /
ENDFIN
CARFIN
'LGR2' 6 6 4 5 1 3 3 6 3 /
ENDFIN
CARFIN
'LGR3' 7 7 3 4 1 3 3 6 3 /
ENDFIN
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

} // Anonymous namespace

BOOST_AUTO_TEST_CASE(InitOutputAlongASealedFault)
{
    using TypeTag = Opm::Properties::TTag::TestTypeTag;

    const auto& comm = Opm::FlowGenericVanguard::comm();
    const auto dir = std::filesystem::temp_directory_path() / ("test_lgr_sealed_fault_ranks_np" + std::to_string(comm.size()));
    const auto deck_file = (dir / "LGR_SEALED_FAULT.DATA").string();
    if (comm.rank() == 0) {
        std::filesystem::create_directories(dir);
        std::ofstream os(deck_file);
        os << deck;
    }
    comm.barrier();

    // The trilinear refinement is flow's default.  The INIT output is
    // extracted while the simulator is set up.
    BOOST_CHECK_NO_THROW(Opm::initSimulator<TypeTag>(deck_file.c_str(),
                                                     "test_lgr_sealed_fault_ranks",
                                                     /*threads_per_process=*/1, "low",
                                                     { "--output-dir=" + dir.string() }));
    comm.barrier();
    if (comm.rank() == 0) {
        // Setting up the simulator extracted the INIT output without a throw.
        std::filesystem::remove_all(dir);
    }
}
