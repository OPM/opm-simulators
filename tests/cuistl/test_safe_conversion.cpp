/*
  Copyright 2023 SINTEF AS

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

#define BOOST_TEST_MODULE TestSafeConversion

#include <boost/test/unit_test.hpp>
#include <opm/simulators/linalg/cuistl/detail/safe_conversion.hpp>

BOOST_AUTO_TEST_CASE(TestThrowsOutofRange)
{
    BOOST_CHECK_THROW(Opm::cuistl::detail::convert(size_t(std::numeric_limits<int>::max()) + size_t(1));
                      , std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(TestConvertInRange)
{
    // This might seem slow, but it is really fast:
    for (size_t i = 0; i <= size_t(1024 * 1024); ++i) {
        BOOST_CHECK_EQUAL(int(i), Opm::cuistl::detail::convert(i));
    }

    BOOST_CHECK_EQUAL(std::numeric_limits<int>::max(),
                      Opm::cuistl::detail::convert(size_t(std::numeric_limits<int>::max())));
}
