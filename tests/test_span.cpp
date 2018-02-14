/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.

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

#define BOOST_TEST_MODULE SpanTest

#include <opm/autodiff/AutoDiffHelpers.hpp>

#include <boost/test/unit_test.hpp>

#include <iostream>

using namespace Opm;

BOOST_AUTO_TEST_CASE(OneArgConstr)
{
    const int num = 4;
    const Span s(num);
    BOOST_CHECK_EQUAL(s.size(), num);
    for (int i = 0; i < num; ++i) {
        BOOST_CHECK_EQUAL(s[i], i);
    }
    int count = 0;
    for (Span::const_iterator it = s.begin(); it != s.end(); ++it) {
        BOOST_CHECK_EQUAL(*it, count);
        ++count;
    }
    BOOST_CHECK_EQUAL(count, num);
}

BOOST_AUTO_TEST_CASE(ThreeArgConstr)
{
    const int num = 3;
    const int stride = 7;
    const int start = 5;
    const int seq[num] = { start, start + 1*stride, start + 2*stride };

    const Span s(num, stride, start);
    BOOST_CHECK_EQUAL(s.size(), num);
    for (int i = 0; i < num; ++i) {
        BOOST_CHECK_EQUAL(s[i], seq[i]);
    }
    int count = 0;
    for (Span::const_iterator it = s.begin(); it != s.end(); ++it) {
        BOOST_CHECK_EQUAL(*it, seq[count]);
        ++count;
    }
    BOOST_CHECK_EQUAL(count, num);
}
