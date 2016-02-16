/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.

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

#define BOOST_TEST_MODULE fastSparseSumTest

#include <opm/autodiff/fastSparseSum.hpp>

#include <boost/test/unit_test.hpp>

typedef Eigen::SparseMatrix<double> Sp;
using namespace Opm;


bool
operator ==(const Eigen::SparseMatrix<double>& A,
            const Eigen::SparseMatrix<double>& B)
{
    // Two SparseMatrices are equal if
    //   0) They have the same ordering (enforced by equal types)
    //   1) They have the same outer and inner dimensions
    //   2) They have the same number of non-zero elements
    //   3) They have the same sparsity structure
    //   4) The non-zero elements are equal

    // 1) Outer and inner dimensions
    bool eq =       (A.outerSize() == B.outerSize());
    eq      = eq && (A.innerSize() == B.innerSize());

    // 2) Equal number of non-zero elements
    eq      = eq && (A.nonZeros() == B.nonZeros());

    for (typename Eigen::SparseMatrix<double>::Index
             k0 = 0, kend = A.outerSize(); eq && (k0 < kend); ++k0) {
        for (typename Eigen::SparseMatrix<double>::InnerIterator
                 iA(A, k0), iB(B, k0); eq && (iA && iB); ++iA, ++iB) {
            // 3) Sparsity structure
            eq = (iA.row() == iB.row()) && (iA.col() == iB.col());

            // 4) Equal non-zero elements
            eq = eq && (iA.value() == iB.value());
        }
    }

    return eq;

    // Note: Investigate implementing this operator as
    // return A.cwiseNotEqual(B).count() == 0;
}


BOOST_AUTO_TEST_CASE(function_fastSparseSum)
{
    // s1 has last entry.
    {
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> s1(6, 1);
        s1 << 1, 2, 0, 4, 0, 6;
        Sp ss1(s1.sparseView());
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> s2(6, 1);
        s2 << 0, 2, 3, 4, 0, 0;
        Sp ss2(s2.sparseView());

        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> r(6, 1);
        r << 1, 4, 3, 8, 0, 6;
        Sp reference(r.sparseView());

        Sp sum;
        Opm::fastSparseSum(ss1, ss2, sum);
        for (int row = 0; row < sum.rows(); ++row) {
            BOOST_CHECK_EQUAL(sum.coeff(row, 0), reference.coeff(row, 0));
        }
    }
    // s2 has last entry.
    {
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> s1(6, 1);
        s1 << 1, 2, 0, 4, 0, 0;
        Sp ss1(s1.sparseView());
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> s2(6, 1);
        s2 << 0, 2, 3, 4, 0, 6;
        Sp ss2(s2.sparseView());

        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> r(6, 1);
        r << 1, 4, 3, 8, 0, 6;
        Sp reference(r.sparseView());

        Sp sum;
        Opm::fastSparseSum(ss1, ss2, sum);
        for (int row = 0; row < sum.rows(); ++row) {
            BOOST_CHECK_EQUAL(sum.coeff(row, 0), reference.coeff(row, 0));
        }
    }
    // Both have last entry.
    {
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> s1(6, 1);
        s1 << 1, 2, 0, 4, 5, 0;
        Sp ss1(s1.sparseView());
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> s2(6, 1);
        s2 << 0, 2, 3, 4, 5, 0;
        Sp ss2(s2.sparseView());

        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> r(6, 1);
        r << 1, 4, 3, 8, 10, 0;
        Sp reference(r.sparseView());

        Sp sum;
        Opm::fastSparseSum(ss1, ss2, sum);
        for (int row = 0; row < sum.rows(); ++row) {
            BOOST_CHECK_EQUAL(sum.coeff(row, 0), reference.coeff(row, 0));
        }
    }
    // Random structures.
    {
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> s1(2, 3);
        s1 << 1, 2, 0, 4, 5, 0;
        Sp ss1(s1.sparseView());
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> s2(2, 3);
        s2 << 0, 2, 3, 4, 5, 0;
        Sp ss2(s2.sparseView());

        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> r(2, 3);
        r << 1, 4, 3, 8, 10, 0;
        Sp reference(r.sparseView());

        Sp sum;
        Opm::fastSparseSum(ss1, ss2, sum);
        for (int row = 0; row < sum.rows(); ++row) {
            BOOST_CHECK(sum == ss1 + ss2);
        }
    }
    // With a zero column.
    {
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> s1(2, 3);
        s1 << 1, 2, 0, 4, 5, 0;
        Sp ss1(s1.sparseView());
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> s2(2, 3);
        s2 << 0, 2, 3, 4, 0, 0;
        Sp ss2(s2.sparseView());

        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> r(2, 3);
        r << 1, 4, 3, 8, 5, 0;
        Sp reference(r.sparseView());

        Sp sum;
        Opm::fastSparseSum(ss1, ss2, sum);
        for (int row = 0; row < sum.rows(); ++row) {
            BOOST_CHECK(sum == ss1 + ss2);
        }
    }
}

