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
#define BOOST_TEST_MODULE BcconLgr
#include "SimulatorFixture.hpp"

#include <opm/input/eclipse/Schedule/BCState.hpp>

#include <opm/simulators/flow/FlowProblemBlackoil.hpp>

#include <boost/test/unit_test.hpp>

using SimulatorFixture = Opm::SimulatorFixture;
BOOST_GLOBAL_FIXTURE(SimulatorFixture);

// A boundary condition on a face of a refined cell applies to every refined
// cell on that face: the faces that carry it cover the whole host face.
BOOST_AUTO_TEST_CASE(BoundaryConditionOnRefinedCell)
{
    using TypeTag = Opm::Properties::TTag::TestTypeTag;

    // CARFIN is parsed at low strictness, as in the other runs with local grids.
    auto simulator = Opm::initSimulator<TypeTag>("BCCON_LGR.DATA", "test_bccon_lgr",
                                                 /*threads_per_process=*/1, "low");
    simulator->setEpisodeIndex(-1);
    simulator->setEpisodeLength(0.0);
    simulator->startNextEpisode(/*episodeStartTime=*/0.0, /*episodeLength=*/1e30);

    const auto& problem = simulator->problem();
    const auto& gridView = simulator->vanguard().gridView();
    const auto& elemMapper = simulator->model().elementMapper();

    constexpr int yMinus = 2; // intersection index of the Y- face
    int faces = 0;
    double area = 0.0;
    for (const auto& elem : elements(gridView)) {
        const unsigned elemIdx = elemMapper.index(elem);
        for (const auto& intersection : intersections(gridView, elem)) {
            if (intersection.boundary() && (intersection.indexInInside() == yMinus) &&
                (problem.boundaryCondition(elemIdx, yMinus).first == Opm::BCType::DIRICHLET))
            {
                ++faces;
                area += intersection.geometry().volume();
            }
        }
    }

    // The Y- face of cell (2,1,1) is 100 m x 10 m, split into 3 x 3 faces.
    BOOST_CHECK_EQUAL(faces, 9);
    BOOST_CHECK_CLOSE(area, 100.0 * 10.0, 1.0e-8);
}
