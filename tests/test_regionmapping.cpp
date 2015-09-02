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

#define BOOST_TEST_MODULE RegionMapping

#include <opm/core/utility/platform_dependent/disable_warnings.h>
#include <boost/test/unit_test.hpp>
#include <opm/core/utility/platform_dependent/reenable_warnings.h>

/* --- our own headers --- */

#include <opm/core/utility/RegionMapping.hpp>

#include <algorithm>
#include <map>

BOOST_AUTO_TEST_SUITE ()


BOOST_AUTO_TEST_CASE (Forward)
{
    std::vector<int> regions = { 2, 5, 2, 4, 2, 7, 6, 3, 6 };

    Opm::RegionMapping<> rm(regions);

    for (decltype(regions.size()) i = 0, n = regions.size(); i < n; ++i) {
        BOOST_CHECK_EQUAL(rm.region(i), regions[i]);
    }
}


BOOST_AUTO_TEST_CASE (ActiveRegions)
{
    //                           0  1  2  3  4  5  6  7  8
    std::vector<int> regions = { 2, 5, 2, 4, 2, 7, 6, 3, 6 };

    Opm::RegionMapping<> rm(regions);

    std::vector<int> region_ids = { 2, 3, 4, 5, 6, 7 };

    auto active = [&region_ids](const int reg)
        {
            auto b = region_ids.begin();
            auto e = region_ids.end();

            return std::find(b, e, reg) != e;
        };

    BOOST_CHECK_EQUAL(rm.activeRegions().size(), region_ids.size());

    for (const auto& reg : rm.activeRegions()) {
        BOOST_CHECK(active(reg));
    }
}


BOOST_AUTO_TEST_CASE (Consecutive)
{
    using RegionCells = std::map<int, std::vector<int>>;

    //                           0  1  2  3  4  5  6  7  8
    std::vector<int> regions = { 2, 5, 2, 4, 2, 7, 6, 3, 6 };

    Opm::RegionMapping<> rm(regions);

    std::vector<int> region_ids = { 2, 3, 4, 5, 6, 7 };
    RegionCells      region_cells;
    {
        using VT = RegionCells::value_type;

        region_cells.insert(VT(2, { 0, 2, 4 }));
        region_cells.insert(VT(3, { 7 }));
        region_cells.insert(VT(4, { 3 }));
        region_cells.insert(VT(5, { 1 }));
        region_cells.insert(VT(6, { 6, 8 }));
        region_cells.insert(VT(7, { 5 }));
    }

    for (const auto& reg : region_ids) {
        const auto& cells  = rm.cells(reg);
        const auto& expect = region_cells[reg];

        BOOST_CHECK_EQUAL_COLLECTIONS(cells .begin(), cells .end(),
                                      expect.begin(), expect.end());
    }

    // Verify that there are no cells in unused regions 0 and 1.
    for (const auto& r : { 0, 1 }) {
        BOOST_CHECK(rm.cells(r).empty());
    }
}


BOOST_AUTO_TEST_CASE (NonConsecutive)
{
    using RegionCells = std::map<int, std::vector<int>>;

    //                           0  1  2  3  4  5  6  7  8
    std::vector<int> regions = { 2, 4, 2, 4, 2, 7, 6, 3, 6 };

    Opm::RegionMapping<> rm(regions);

    std::vector<int> region_ids = { 2, 3, 4, 6, 7 };
    RegionCells      region_cells;
    {
        using VT = RegionCells::value_type;

        region_cells.insert(VT(2, { 0, 2, 4 }));
        region_cells.insert(VT(3, { 7 }));
        region_cells.insert(VT(4, { 1, 3 }));
        region_cells.insert(VT(6, { 6, 8 }));
        region_cells.insert(VT(7, { 5 }));
    }

    for (const auto& reg : region_ids) {
        const auto& cells  = rm.cells(reg);
        const auto& expect = region_cells[reg];

        BOOST_CHECK_EQUAL_COLLECTIONS(cells .begin(), cells .end(),
                                      expect.begin(), expect.end());
    }

    // Verify that there are no cells in unused regions 0, 1, and 5.
    for (const auto& r : { 0, 1, 5 }) {
        BOOST_CHECK(rm.cells(r).empty());
    }
}


BOOST_AUTO_TEST_SUITE_END()
