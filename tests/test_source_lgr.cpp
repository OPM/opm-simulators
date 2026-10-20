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
#define BOOST_TEST_MODULE SourceLgr
#include "SimulatorFixture.hpp"

#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/simulators/flow/FlowProblemBlackoil.hpp>

#include <boost/test/unit_test.hpp>

using SimulatorFixture = Opm::SimulatorFixture;
BOOST_GLOBAL_FIXTURE(SimulatorFixture);

// A SOURCE in a refined cell is shared by the cells of its local grid: summed
// over all cells, the injected rate is the rate in the deck.
BOOST_AUTO_TEST_CASE(SourceInRefinedCell)
{
    using TypeTag = Opm::Properties::TTag::TestTypeTag;
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;
    using RateVector = Opm::GetPropType<TypeTag, Opm::Properties::RateVector>;

    // CARFIN is parsed at low strictness, as in the other runs with local grids.
    auto simulator = Opm::initSimulator<TypeTag>("SOURCE_LGR.DATA", "test_source_lgr",
                                                 /*threads_per_process=*/1, "low");
    simulator->model().applyInitialSolution();
    simulator->setEpisodeIndex(-1);
    simulator->setEpisodeLength(0.0);
    simulator->startNextEpisode(/*episodeStartTime=*/0.0, /*episodeLength=*/1e30);
    simulator->setTimeStepSize(Opm::unit::day);

    const auto& problem = simulator->problem();
    const auto& model = simulator->model();
    const unsigned waterIdx = FluidSystem::canonicalToActiveCompIdx(FluidSystem::waterCompIdx);

    double injected = 0.0;
    for (unsigned dofIdx = 0; dofIdx < model.numGridDof(); ++dofIdx) {
        RateVector rate(0.0);
        problem.addToSourceDense(rate, dofIdx, /*timeIdx=*/0);
        injected += Opm::getValue(rate[waterIdx]) * model.dofTotalVolume(dofIdx);
    }

    // SOURCE 2 2 1 WATER 5000 kg/day, as a surface volume rate.
    const double expected = 5000.0 * Opm::unit::kilogram / Opm::unit::day
        / FluidSystem::referenceDensity(FluidSystem::waterPhaseIdx, /*regionIdx=*/0);
    BOOST_CHECK_CLOSE(injected, expected, 1.0e-8);
}
