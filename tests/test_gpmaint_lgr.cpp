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
#define BOOST_TEST_MODULE GpmaintLgr
#include "SimulatorFixture.hpp"

#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/simulators/flow/FlowProblemBlackoil.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>

#include <boost/test/unit_test.hpp>

using SimulatorFixture = Opm::SimulatorFixture;
BOOST_GLOBAL_FIXTURE(SimulatorFixture);

// GPMAINT averages the pressure of a FIPNUM region over the refined grid; the
// region array must cover every refined cell.
BOOST_AUTO_TEST_CASE(RegionPressureWithLgr)
{
    using TypeTag = Opm::Properties::TTag::TestTypeTag;

    // CARFIN is parsed at low strictness, as in the other runs with local grids.
    auto simulator = Opm::initSimulator<TypeTag>("GPMAINT_LGR.DATA", "test_gpmaint_lgr",
                                                 /*threads_per_process=*/1, "low");
    simulator->model().applyInitialSolution();
    simulator->setEpisodeIndex(-1);
    simulator->setEpisodeLength(0.0);
    simulator->startNextEpisode(/*episodeStartTime=*/0.0, /*episodeLength=*/1e30);
    simulator->setTimeStepSize(Opm::unit::day);
    simulator->problem().resetIterationForNewTimestep();

    auto& wellModel = simulator->problem().wellModel();
    wellModel.beginReportStep(0);
    BOOST_CHECK_NO_THROW(wellModel.beginTimeStep());
}
