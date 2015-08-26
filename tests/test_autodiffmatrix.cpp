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

#define BOOST_TEST_MODULE AutoDiffMatrixTest

#include <opm/autodiff/AutoDiffMatrix.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>

#include <boost/test/unit_test.hpp>

typedef Eigen::SparseMatrix<double> Sp;
typedef Opm::AutoDiffMatrix Mat;
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


BOOST_AUTO_TEST_CASE(Initialization)
{
    // Setup.
    Mat z = Mat(AutoDiffMatrix::ZeroMatrix, 3);

    Mat i = Mat(AutoDiffMatrix::IdentityMatrix, 3);

    Eigen::Array<double, Eigen::Dynamic, 1> d1(3);
    d1 << 0.2, 1.2, 13.4;
    Mat d = Mat(d1.matrix().asDiagonal());

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> s1(3,2);
    s1 <<
        1.0, 0.0, 2.0,
        0.0, 1.0, 0.0;
    Sp s2(s1.sparseView());
    Mat s = Mat(s2);
}

BOOST_AUTO_TEST_CASE(EigenConversion)
{
    // Setup.
    Mat z = Mat(AutoDiffMatrix::ZeroMatrix, 3);

    Mat i = Mat(AutoDiffMatrix::IdentityMatrix, 3);

    Eigen::Array<double, Eigen::Dynamic, 1> d1(3);
    d1 << 0.2, 1.2, 13.4;
    Mat d = Mat(d1.matrix().asDiagonal());

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> s1(3,2);
    s1 <<
        1.0, 0.0, 2.0,
        0.0, 1.0, 0.0;
    Sp s2(s1.sparseView());
    Mat s = Mat(s2);

    // Convert to Eigen::SparseMatrix
    Sp x;
    z.toSparse(x);
    Sp z1(3,3);
    BOOST_CHECK(x == z1);
    i.toSparse(x);
    Sp i1(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Identity(3,3).sparseView());
    BOOST_CHECK(x == i1);
    d.toSparse(x);
    Sp d2 = spdiag(d1);
    BOOST_CHECK(x == d2);
    s.toSparse(x);
    BOOST_CHECK(x == s2);
}



BOOST_AUTO_TEST_CASE(AdditionOps)
{
    // Setup.
    Mat z = Mat(AutoDiffMatrix::ZeroMatrix, 3);
    Sp zs(3,3);

    Mat i = Mat(AutoDiffMatrix::IdentityMatrix, 3);
    Sp is(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Identity(3,3).sparseView());

    Eigen::Array<double, Eigen::Dynamic, 1> d1(3);
    d1 << 0.2, 1.2, 13.4;
    Mat d = Mat(d1.matrix().asDiagonal());
    Sp ds = spdiag(d1);

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> s1(3,3);
    s1 <<
        1.0, 0.0, 2.0,
        0.0, 1.0, 0.0,
	0.0, 0.0, 2.0;
    Sp ss(s1.sparseView());
    Mat s = Mat(ss);

    // Convert to Eigen::SparseMatrix
    Sp x;
    z.toSparse(x);
    BOOST_CHECK(x == zs);
    i.toSparse(x);
    BOOST_CHECK(x == is);
    d.toSparse(x);
    BOOST_CHECK(x == ds);
    s.toSparse(x);
    BOOST_CHECK(x == ss);

    // Adding zero.
    auto zpz = z + z;
    zpz.toSparse(x);
    BOOST_CHECK(x == zs);
    auto ipz = i + z;
    ipz.toSparse(x);
    BOOST_CHECK(x == is);
    auto dpz = d + z;
    dpz.toSparse(x);
    BOOST_CHECK(x == ds);
    auto spz = s + z;
    spz.toSparse(x);
    BOOST_CHECK(x == ss);
}

BOOST_AUTO_TEST_CASE(MultOps)
{
    // Setup.
    Mat z = Mat(AutoDiffMatrix::ZeroMatrix, 3);
    Sp zs(3,3);

    Mat i = Mat(AutoDiffMatrix::IdentityMatrix, 3);
    Sp is(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Identity(3,3).sparseView());

    Eigen::Array<double, Eigen::Dynamic, 1> d1(3);
    d1 << 0.2, 1.2, 13.4;
    Mat d = Mat(d1.matrix().asDiagonal());
    Sp ds = spdiag(d1);

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> s1(3,3);
    s1 <<
        1.0, 0.0, 2.0,
        0.0, 1.0, 0.0,
	0.0, 0.0, 2.0;
    Sp ss(s1.sparseView());
    Mat s = Mat(ss);

    // Convert to Eigen::SparseMatrix
    Sp x;
    z.toSparse(x);
    BOOST_CHECK(x == zs);
    i.toSparse(x);
    BOOST_CHECK(x == is);
    d.toSparse(x);
    BOOST_CHECK(x == ds);
    s.toSparse(x);
    BOOST_CHECK(x == ss);

    // Multiply by diagonal matrix.
    auto ztd = z * d;
    ztd.toSparse(x);
    BOOST_CHECK(x == zs*ds);
    auto itd = i * d;
    itd.toSparse(x);
    BOOST_CHECK(x == is*ds);
    auto dtd = d * d;
    dtd.toSparse(x);
    BOOST_CHECK(x == ds*ds);
    auto std = s * d;
    std.toSparse(x);
    BOOST_CHECK(x == ss*ds);
}

