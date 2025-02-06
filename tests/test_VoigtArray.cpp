/*
  Copyright 2025 Equinor ASA.

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

#include <opm/simulators/utils/VoigtArray.hpp>

#define BOOST_TEST_MODULE VoigArrayTest
#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

#include <stdexcept>

namespace {
    #if FLOW_INSTANTIATE_FLOAT
    using Types = boost::mpl::list<float,double>;
    #else
    using Types = boost::mpl::list<double>;
    #endif
}

BOOST_AUTO_TEST_CASE_TEMPLATE(DefaultConstructed, Scalar, Types)
{
    Opm::VoigtArray<Scalar> array;
    static constexpr auto& indices = Opm::VoigtArray<Scalar>::indices;
    BOOST_CHECK_EQUAL(array.size(), 6);
    for (const auto i : indices) {
        BOOST_CHECK(array[i].empty());
    }

    for (const auto i : indices) {
        BOOST_CHECK(array[i].empty());
        BOOST_CHECK_THROW(array(i, 1), std::out_of_range);
    }

    array.resize(5);
    for (const auto i : indices) {
        BOOST_CHECK(!array[i].empty());
        BOOST_CHECK_NO_THROW(array(i, 1));
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(SizeConstructed, Scalar, Types)
{
    Opm::VoigtArray<Scalar> array(4);
    BOOST_CHECK_EQUAL(array.size(), 6);
    static constexpr auto& indices = Opm::VoigtArray<Scalar>::indices;
    static constexpr auto& uindices = Opm::VoigtArray<Scalar>::unique_indices;
    for (const auto i : indices) {
        BOOST_CHECK_EQUAL(array[i].size(), 4);
    }

    for (const auto i : uindices) {
        std::fill(array[i].begin(), array[i].end(), static_cast<Scalar>(i));
    }

    for (const auto i: indices) {
        BOOST_CHECK_EQUAL(array[i].size(), 4);
        BOOST_CHECK_EQUAL(array[i][1], static_cast<Scalar>(i));
        BOOST_CHECK_EQUAL(array(i,1), static_cast<Scalar>(i));
    }

    BOOST_CHECK_EQUAL(array(Opm::VoigtIndex::XY, 1),
                      array(Opm::VoigtIndex::YX, 1));

    BOOST_CHECK_EQUAL(array(Opm::VoigtIndex::XZ, 2),
                      array(Opm::VoigtIndex::ZX, 2));

    BOOST_CHECK_EQUAL(array(Opm::VoigtIndex::YZ, 3),
                      array(Opm::VoigtIndex::ZY, 3));
}
