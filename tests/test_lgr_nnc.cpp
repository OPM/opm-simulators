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
 * \brief Test the output of an NNC keyword in a grid with an LGR
 *
 * LGR_NNC.DATA refines cell (1,1,1) of a 4x1x1 grid into 2x1x1 and connects
 * cells (2,1,1) and (4,1,1) with the NNC keyword.  The refined cells come
 * first in the leaf grid, so the two connected cells have other leaf indices
 * than level-0 indices.  The connection must be written with the
 * transmissibility given in the deck.
 */

#define BOOST_TEST_MODULE LgrNnc

#include "SimulatorFixture.hpp"

#include <opm/simulators/flow/FlowProblemBlackoil.hpp>

#include <boost/test/unit_test.hpp>

#include <algorithm>

using SimulatorFixture = Opm::SimulatorFixture;
BOOST_GLOBAL_FIXTURE(SimulatorFixture);

BOOST_AUTO_TEST_CASE(DeckNncBetweenUnrefinedCells)
{
    using TypeTag = Opm::Properties::TTag::TestTypeTag;

    // CARFIN is parsed at low strictness, as in the other runs with local grids.
    auto simulator = Opm::initSimulator<TypeTag>("LGR_NNC.DATA", "test_lgr_nnc",
                                                 /*threads_per_process=*/1, "low");

    const auto& deckNnc = simulator->vanguard().eclState().getInputNNC().input();
    BOOST_REQUIRE_EQUAL(deckNnc.size(), 1U);

    const auto& outputNnc = simulator->problem().eclWriter().getOutputNnc();
    BOOST_REQUIRE(!outputNnc.empty());

    const auto written = std::find_if(outputNnc.front().begin(), outputNnc.front().end(),
                                      [&deckNnc](const auto& nnc)
                                      {
                                          return (nnc.cell1 == deckNnc.front().cell1)
                                              && (nnc.cell2 == deckNnc.front().cell2);
                                      });
    BOOST_REQUIRE(written != outputNnc.front().end());
    BOOST_CHECK_CLOSE(written->trans, deckNnc.front().trans, 1.0e-8);
}
