/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.

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

#include <config.h>

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif

#define NVERBOSE  // Suppress own messages when throw()ing

#define BOOST_TEST_MODULE OPM-multiPhaseUpwind
#include <boost/test/unit_test.hpp>

#include <opm/autodiff/multiPhaseUpwind.hpp>

// For all the following cases we test a setup with two cells,
// forming a gravity column:
//
//     -------------------
//     |                 |
//     |     Cell 1      |
//     |                 |
//     |                 |  |
//     -------------------  | flux
//     |                 |  V
//     |     Cell 2      |
//     |                 |
//     |                 |
//     -------------------
//
// The gravity-related head differences are (4, -1, -2) for (w, o, g).
// The mobilities are all 1 and the transmissibility is 1.
// The total flux from cell 1 to 2 will vary from case to case.

BOOST_AUTO_TEST_CASE(GravityColumnLowFlux)
{
    // Case 1: a gravity column with two cells and low total flux.
    // The total flux from cell 1 to 2 is 1.0.

    const std::array<double, 3> head_diff = {{ 4.0, -1.0, -2.0 }};
    const std::array<double, 3> mob1 = {{ 1.0, 1.0, 1.0 }};
    const std::array<double, 3> mob2 = {{ 1.0, 1.0, 1.0 }};
    const double transmissibility = 1.0;
    const double flux = 1.0;

    const std::array<double, 3> expected_upw = {{ 1.0, -1.0, -1.0 }};
    const std::array<double, 3> upw = Opm::connectionMultiPhaseUpwind(head_diff, mob1, mob2, transmissibility, flux);
    BOOST_CHECK_EQUAL(upw[0], expected_upw[0]);
    BOOST_CHECK_EQUAL(upw[1], expected_upw[1]);
    BOOST_CHECK_EQUAL(upw[2], expected_upw[2]);
}


BOOST_AUTO_TEST_CASE(GravityColumnMediumFlux)
{
    // Case 2: a gravity column with two cells and medium-sized total flux.
    // The total flux from cell 1 to 2 is 5.0.

    const std::array<double, 3> head_diff = {{ 4.0, -1.0, -2.0 }};
    const std::array<double, 3> mob1 = {{ 1.0, 1.0, 1.0 }};
    const std::array<double, 3> mob2 = {{ 1.0, 1.0, 1.0 }};
    const double transmissibility = 1.0;
    const double flux = 5.0;

    const std::array<double, 3> expected_upw = {{ 1.0, 1.0, -1.0 }};
    const std::array<double, 3> upw = Opm::connectionMultiPhaseUpwind(head_diff, mob1, mob2, transmissibility, flux);
    BOOST_CHECK_EQUAL(upw[0], expected_upw[0]);
    BOOST_CHECK_EQUAL(upw[1], expected_upw[1]);
    BOOST_CHECK_EQUAL(upw[2], expected_upw[2]);
}


BOOST_AUTO_TEST_CASE(GravityColumnHighFlux)
{
    // Case 3: a gravity column with two cell and high total flux.
    // The total flux from cell 1 to 2 is 10.0.

    const std::array<double, 3> head_diff = {{ 4.0, -1.0, -2.0 }};
    const std::array<double, 3> mob1 = {{ 1.0, 1.0, 1.0 }};
    const std::array<double, 3> mob2 = {{ 1.0, 1.0, 1.0 }};
    const double transmissibility = 1.0;
    const double flux = 10.0;

    const std::array<double, 3> expected_upw = {{ 1.0, 1.0, 1.0 }};
    const std::array<double, 3> upw = Opm::connectionMultiPhaseUpwind(head_diff, mob1, mob2, transmissibility, flux);
    BOOST_CHECK_EQUAL(upw[0], expected_upw[0]);
    BOOST_CHECK_EQUAL(upw[1], expected_upw[1]);
    BOOST_CHECK_EQUAL(upw[2], expected_upw[2]);
}
