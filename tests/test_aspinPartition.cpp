/*
  Copyright 2021 SINTEF Digital, Mathematics and Cybernetics.

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

#define BOOST_TEST_MODULE OPM_test_aspinPartition
#include <boost/test/unit_test.hpp>

#include <opm/simulators/flow/aspinPartition.hpp>


BOOST_AUTO_TEST_CASE(fileBased)
{
    auto [part, num_part]  = Opm::partitionCells(10);
    BOOST_CHECK_EQUAL(num_part, 3);
    std::vector<int> expected = { 0, 0, 1, 1, 2, 2, 1, 1, 0, 0 };
    BOOST_CHECK_EQUAL_COLLECTIONS(expected.begin(), expected.end(), part.begin(), part.end());
}

