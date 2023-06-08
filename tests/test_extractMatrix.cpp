/*
  Copyright 2023 SINTEF Digital, Mathematics and Cybernetics.

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

#define BOOST_TEST_MODULE OPM_test_extractMatrix
#include <boost/test/unit_test.hpp>

#include <opm/simulators/linalg/extractMatrix.hpp>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>

using B = Dune::FieldMatrix<double, 2, 2>;
using M = Dune::BCRSMatrix<B>;
using V = Dune::BlockVector<Dune::FieldVector<double, 2>>;

M build3x3BlockMatrix()
{
    // Build matrix with this pattern:
    //   ( x     ) 
    //   (   x   )
    //   (   x x )
    M res(3, 3, 4, M::row_wise);
    for (auto row = res.createbegin(); row != res.createend(); ++row) {
        row.insert(row.index());
        if (row.index() == 2) {
            row.insert(1);
        }
    }
    // Set all entries to the block (1, 2; 3, x), where x starts at 4
    // but gets incremented for each entry.
    B b;
    b[0][0] = 1.0;
    b[0][1] = 2.0;
    b[1][0] = 3.0;
    b[1][1] = 4.0;
    res[0][0] = b;
    b[1][1] = b[1][1] + 1.0;
    res[1][1] = b;
    b[1][1] = b[1][1] + 1.0;
    res[2][2] = b;
    b[1][1] = b[1][1] + 1.0;
    res[2][1] = b;

    return res;
}

BOOST_AUTO_TEST_CASE(extractMatrix)
{
    auto m1 = build3x3BlockMatrix();
    BOOST_CHECK_EQUAL(m1[2][1][1][1], 7.0);
    std::vector<int> indices = {1, 2};
    auto m2 = Opm::Details::extractMatrix(m1, indices);
    BOOST_CHECK_EQUAL(m2[0][0], m1[1][1]);
    BOOST_CHECK_EQUAL(m2[1][1], m1[2][2]);
    BOOST_CHECK_EQUAL(m2[1][0], m1[2][1]);
}

BOOST_AUTO_TEST_CASE(copySubMatrixAndMatrixEqual)
{
    auto m1 = build3x3BlockMatrix();
    std::vector<int> indices = {1, 2};
    auto m2 = Opm::Details::extractMatrix(m1, indices);
    m1[2][1][0][0] = 0.1234;
    auto m3 = Opm::Details::extractMatrix(m1, indices);
    BOOST_CHECK(m2[1][0] != m1[2][1]);
    BOOST_CHECK(m3[1][0] == m1[2][1]);
    BOOST_CHECK(!Opm::Details::matrixEqual(m2, m3));
    Opm::Details::copySubMatrix(m1, indices, m2);
    BOOST_CHECK(Opm::Details::matrixEqual(m2, m3));
}

BOOST_AUTO_TEST_CASE(vectorOps)
{
    V v1(4);
    v1[0] = { 0.1, 0.2 };
    v1[1] = { 0.3, 0.4 };
    v1[2] = { 0.5, 0.6 };
    v1[3] = { 0.6, 0.8 };
    std::vector<int> indices = { 0, 2, 3 };

    V vref(3);
    vref[0] = { 0.1, 0.2 };
    vref[1] = { 0.5, 0.6 };
    vref[2] = { 0.6, 0.8 };

    auto v2 = Opm::Details::extractVector(v1, indices);
    BOOST_CHECK_EQUAL_COLLECTIONS(vref.begin(), vref.end(), v2.begin(), v2.end());

    V v3 = v1;
    v3[2] = { 1.1, 1.2 };
    v2[1] = { 1.1, 1.2 };
    Opm::Details::setGlobal(v2, indices, v1);
    BOOST_CHECK_EQUAL_COLLECTIONS(v1.begin(), v1.end(), v3.begin(), v3.end());
}
