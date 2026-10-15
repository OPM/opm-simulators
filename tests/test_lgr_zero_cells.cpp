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
 * \brief Test that a load balancing which leaves a process without cells
 *        stops every process with its error.
 *
 * The grid is refined before it is distributed, so the properties are sent
 * from the refined grid.  When the load balancing throws because a process
 * got no cells, every process must stop with that error; the handle that
 * sends the properties must not read a grid that was never distributed.
 */

#define BOOST_TEST_MODULE LgrZeroCells

#include "SimulatorFixture.hpp"

#include <opm/simulators/flow/FlowProblemBlackoil.hpp>

#include <boost/test/unit_test.hpp>

#include <filesystem>
#include <fstream>
#include <memory>
#include <string>

using SimulatorFixture = Opm::SimulatorFixture;
BOOST_GLOBAL_FIXTURE(SimulatorFixture);

namespace {

// Two cells; LGR1 refines the first one into two, so the refined grid has three
// cells: fewer than four processes.
const char* deck = R"(RUNSPEC
DIMENS
 1 1 2 /
OIL
WATER
METRIC
START
 1 'JAN' 2020 /
TABDIMS
/
EQLDIMS
/
GRID
CARFIN
'LGR1' 1 1 1 1 1 1 1 1 2 /
ENDFIN
DX
 2*100 /
DY
 2*100 /
DZ
 2*10 /
TOPS
 2000 /
PORO
 2*0.2 /
PERMX
 2*100 /
PERMY
 2*100 /
PERMZ
 2*10 /
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

// As initSimulator, with the conforming refinement.
template <class TypeTag>
void setupConforming(const std::string& filename)
{
    const std::string deckArg = "--ecl-deck-file-name=" + filename;
    const char* argv[] = { "test_lgr_zero_cells", deckArg.c_str(), "--lgr-backend=conforming" };

    Opm::Parameters::reset();
    Opm::registerAllParameters_<TypeTag>(false);
    Opm::registerEclTimeSteppingParameters<double>();
    Opm::BlackoilModelParameters<double>::registerParameters();
    Opm::Parameters::Register<Opm::Parameters::EnableTerminalOutput>("Do *NOT* use!");
    Opm::Parameters::SetDefault<Opm::Parameters::ThreadsPerProcess>(1);
    Opm::Parameters::endRegistration();
    Opm::setupParameters_<TypeTag>(sizeof(argv) / sizeof(argv[0]), argv,
                                   /*registerParams=*/false, /*allowUnused=*/false,
                                   /*handleHelp=*/true, /*myRank=*/0);
    Opm::FlowGenericVanguard::readDeck(filename, "low");
}

} // Anonymous namespace

BOOST_AUTO_TEST_CASE(ProcessWithoutCellsStopsEveryProcess)
{
    using TypeTag = Opm::Properties::TTag::TestTypeTag;
    using Simulator = Opm::GetPropType<TypeTag, Opm::Properties::Simulator>;

    const auto& comm = Opm::FlowGenericVanguard::comm();
    const auto dir = std::filesystem::temp_directory_path() / ("test_lgr_zero_cells_np" + std::to_string(comm.size()));
    const auto deck_file = (dir / "LGR_ZERO_CELLS.DATA").string();
    if (comm.rank() == 0) {
        std::filesystem::create_directories(dir);
        std::ofstream os(deck_file);
        os << deck;
    }
    comm.barrier();

    setupConforming<TypeTag>(deck_file);
    if (comm.size() > 3) {
        // More processes than cells: every process stops with the error.
        BOOST_CHECK_THROW(std::make_unique<Simulator>(), std::exception);
    }
    else {
        std::unique_ptr<Simulator> simulator;
        BOOST_CHECK_NO_THROW(simulator = std::make_unique<Simulator>());
        BOOST_CHECK_EQUAL(comm.sum(simulator->vanguard().grid().leafGridView().size(0) > 0 ? 1 : 0),
                          comm.size());
    }

    comm.barrier();
    if (comm.rank() == 0) {
        std::filesystem::remove_all(dir);
    }
}
