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

#include "config.h"

/* --- Boost.Test boilerplate --- */
#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif

#define NVERBOSE  // Suppress own messages when throw()ing

#define BOOST_TEST_MODULE UnitsTest
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

/* --- our own headers --- */

#include <opm/core/utility/RegionMapping.hpp>


BOOST_AUTO_TEST_SUITE ()


BOOST_AUTO_TEST_CASE (RegionMapping)
{
    //                           0  1  2  3  4  5  6  7  8
    std::vector<int> regions = { 2, 5, 2, 4, 2, 7, 6, 3, 6 };
    Opm::RegionMapping<> rm(regions);
    for (size_t i = 0; i < regions.size(); ++i) {
        BOOST_CHECK_EQUAL(rm.region(i), regions[i]);
    }
    std::vector<int> region_ids = { 2, 3, 4, 5, 6, 7 };
    std::vector< std::vector<int> > region_cells = {  { 0, 2, 4 },  { 7 },  { 3 },  { 1 },  { 6, 8 },  { 5 }  };
    BOOST_REQUIRE_EQUAL(rm.numRegions(), region_ids.size());
    for (size_t i = 0; i < region_ids.size(); ++i) {
        auto cells = rm.cells(region_ids[i]);
        BOOST_REQUIRE_EQUAL(std::distance(cells.begin(), cells.end()), region_cells[i].size());
        size_t count = 0;
        for (auto iter = cells.begin(); iter != cells.end(); ++iter, ++count) {
            BOOST_CHECK_EQUAL(*iter, region_cells[i][count]);
        }
    }
}



BOOST_AUTO_TEST_SUITE_END()
