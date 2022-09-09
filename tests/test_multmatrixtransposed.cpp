/*
  Copyright 2018 Statoil ASA

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

#define BOOST_TEST_MODULE MultMatrixTransposed
#include <boost/test/unit_test.hpp>
#include <opm/simulators/linalg/SmallDenseMatrixUtils.hpp>

#include <dune/common/fmatrix.hh>

using namespace Dune;
using namespace Opm::detail;

BOOST_AUTO_TEST_CASE(testmultmatrixtrans)
{
    using Scalar = FieldMatrix<double,1,1>;
    Scalar a1 = {{ 4 }}, b1 = {{ 2 }}, res1 = {{}}, resExpect1 = {{ 8 }};
    multMatrixTransposed(a1, b1, res1);
    BOOST_CHECK_EQUAL(res1, resExpect1);

    using Mat2 =  FieldMatrix<double,2,2>;


    Mat2 a2 = {{ 1, 2 }, { 3, 4} }, b2 = {{ 3, 4 }, { 5, 6} };
    // resulting matrix filled with temporay values to avoid initialization error
    // from dune 2.4.1
    Mat2 res2 = {{0, 0}, {0, 0}};
    Mat2 resExpect2 = {{ 18, 22 }, {26, 32} };
    multMatrixTransposed(a2, b2, res2);
    BOOST_CHECK_EQUAL(res2, resExpect2);

    using Mat3 =  FieldMatrix<double,3,3>;


    Mat3 a3 = {{ 1, 2, 3 }, { 3, 4, 5}, {6, 7, 8} };

    Mat3 b3 = {{ 3, 4, 5 }, { 5, 6, 7}, {7, 8, 9} };
    // resulting matrix filled with temporay values to avoid initialization error
    // from dune 2.4.1
    Mat3 res3 = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    Mat3 resExpect3 = {{ 60, 70, 80 }, {75, 88, 101}, { 90, 106, 122 } };
    multMatrixTransposed(a3, b3, res3);
    BOOST_CHECK_EQUAL(res3, resExpect3);
}
