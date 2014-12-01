/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.

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
#define NVERBOSE // to suppress our messages when throwing

#define BOOST_TEST_MODULE AnisotropicEikonalTest
#include <boost/test/unit_test.hpp>

#include <opm/core/tof/AnisotropicEikonal.hpp>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/grid.h>
#include <cmath>

using namespace Opm;

BOOST_AUTO_TEST_CASE(cartesian_2d_a)
{
    const GridManager gm(2, 2);
    const UnstructuredGrid& grid = *gm.c_grid();
    AnisotropicEikonal2d ae(grid);

    const std::vector<double> metric = {
        1, 0, 0, 1,
        1, 0, 0, 1,
        1, 0, 0, 1,
        1, 0, 0, 1
    };
    BOOST_REQUIRE_EQUAL(metric.size(), grid.number_of_cells*grid.dimensions*grid.dimensions);
    const std::vector<int> start = { 0 };
    std::vector<double> sol;
    ae.solve(metric.data(), start, sol);
    BOOST_REQUIRE(!sol.empty());
    BOOST_CHECK_EQUAL(sol.size(), grid.number_of_cells);
    std::vector<double> truth = { 0, 1, 1, std::sqrt(2) };
    BOOST_CHECK_EQUAL_COLLECTIONS(sol.begin(), sol.end(), truth.begin(), truth.end());
}


BOOST_AUTO_TEST_CASE(cartesian_2d_b)
{
    const GridManager gm(3, 2, 1.0, 2.0);
    const UnstructuredGrid& grid = *gm.c_grid();
    AnisotropicEikonal2d ae(grid);

    const std::vector<double> metric = {
        1, 0, 0, 1,
        1, 0, 0, 1,
        1, 0, 0, 1,
        1, 0, 0, 1,
        1, 0, 0, 1,
        1, 0, 0, 1
    };
    BOOST_REQUIRE_EQUAL(metric.size(), grid.number_of_cells*grid.dimensions*grid.dimensions);
    const std::vector<int> start = { 0 };
    std::vector<double> sol;
    ae.solve(metric.data(), start, sol);
    BOOST_REQUIRE(!sol.empty());
    BOOST_CHECK_EQUAL(sol.size(), grid.number_of_cells);
    std::vector<double> truth = { 0, 1, 2, 2, std::sqrt(5), std::sqrt(8) };
    BOOST_CHECK_EQUAL_COLLECTIONS(sol.begin(), sol.end(), truth.begin(), truth.end());
}
