/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.
  Copyright 2016 IRIS AS

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

#define BOOST_TEST_MODULE AutoDiffBlockTest

#include <opm/autodiff/AutoDiffBlock.hpp>

#include <boost/test/unit_test.hpp>

#include <Eigen/Eigen>
#include <Eigen/Sparse>

using namespace Opm;

namespace {
    template <typename Scalar>
    bool
    operator ==(const Eigen::SparseMatrix<Scalar>& A,
                const Eigen::SparseMatrix<Scalar>& B)
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

        for (typename Eigen::SparseMatrix<Scalar>::Index
                 k0 = 0, kend = A.outerSize(); eq && (k0 < kend); ++k0) {
            for (typename Eigen::SparseMatrix<Scalar>::InnerIterator
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


    bool operator==(const AutoDiffMatrix& lhs,
                    const AutoDiffMatrix& rhs)
    {
        Eigen::SparseMatrix<double> lhs_s, rhs_s;

        lhs.toSparse(lhs_s);
        rhs.toSparse(rhs_s);

        return lhs_s == rhs_s;
    }

    void checkClose(const AutoDiffBlock<double>& lhs, const AutoDiffBlock<double>& rhs, double tolerance) {

        BOOST_CHECK(lhs.value().isApprox(rhs.value(), tolerance));

        auto lhs_d = lhs.derivative();
        auto rhs_d = rhs.derivative();

        Eigen::SparseMatrix<double> lhs_s, rhs_s;

        //If lhs has no derivatives, make sure all rhs derivatives are zero
        if (lhs_d.size() == 0) {
            for (size_t i=0; i<rhs_d.size(); ++i) {
                rhs.derivative()[i].toSparse(rhs_s);
                BOOST_CHECK_EQUAL(rhs_s.nonZeros(), 0);
            }
        }
        //If rhs has no derivatives, make sure all lhs derivatives are zero
        else if (rhs_d.size() == 0) {
            for (size_t i=0; i<lhs_d.size(); ++i) {
                lhs.derivative()[i].toSparse(lhs_s);
                BOOST_CHECK_EQUAL(lhs_s.nonZeros(), 0);
            }
        }
        //Else, check that derivatives are close to each other
        else {
            BOOST_CHECK_EQUAL(lhs_d.size(), rhs_d.size());
            for (size_t i=0; i<lhs_d.size(); ++i) {
                lhs.derivative()[i].toSparse(lhs_s);
                rhs.derivative()[i].toSparse(rhs_s);
                BOOST_CHECK(lhs_s.isApprox(rhs_s, tolerance));
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(ConstantInitialisation)
{
    typedef AutoDiffBlock<double> ADB;

    ADB::V v(3);
    v << 0.2, 1.2, 13.4;

    ADB a = ADB::constant(v);
    BOOST_REQUIRE(a.value().matrix() == v.matrix());

    const std::vector<ADB::M>& J = a.derivative();
    for (std::vector<ADB::M>::const_iterator
             b = J.begin(), e = J.end(); b != e; ++b) {
        BOOST_REQUIRE(b->nonZeros() == 0);
    }
}

BOOST_AUTO_TEST_CASE(VariableInitialisation)
{
    typedef AutoDiffBlock<double> ADB;

    std::vector<int> blocksizes = { 3, 1, 2 };

    ADB::V v(3);
    v << 1.0, 2.2, 3.4;

    enum { FirstVar = 0, SecondVar = 1, ThirdVar = 2 };

    ADB x = ADB::variable(FirstVar, v, blocksizes);

    BOOST_REQUIRE(x.value().matrix() == v.matrix());

    const std::vector<ADB::M>& J = x.derivative();
    BOOST_REQUIRE(J[0].nonZeros() == v.size());

    const ADB::M& d = J[0];
    for (int i=0; i<d.cols(); ++i) {
        for (int j=0; j<d.rows(); ++j) {
            if (i==j) {
                BOOST_REQUIRE_EQUAL(d.coeff(j, i), 1.0);
            }
            else {
                BOOST_REQUIRE_EQUAL(d.coeff(j, i), 0.0);
            }
        }
    }

    for (std::vector<ADB::M>::const_iterator
             b = J.begin() + 1, e = J.end(); b != e; ++b) {
        BOOST_REQUIRE(b->nonZeros() == 0);
    }
}

BOOST_AUTO_TEST_CASE(FunctionInitialisation)
{
    typedef AutoDiffBlock<double> ADB;

    std::vector<int>            blocksizes = { 3, 1, 2 };
    std::vector<int>::size_type num_blocks = blocksizes.size();

    enum { FirstVar = 0, SecondVar = 1, ThirdVar = 2 };

    ADB::V v(3);
    v << 1.0, 2.2, 3.4;

    std::vector<ADB::M> jacs(num_blocks);
    for (std::vector<int>::size_type j = 0; j < num_blocks; ++j) {
        Eigen::SparseMatrix<double> sm(blocksizes[FirstVar], blocksizes[j]);
        sm.insert(0,0) = -1.0;
        jacs[j] = ADB::M(sm);
    }

    ADB::V v_copy(v);
    std::vector<ADB::M> jacs_copy(jacs);
    ADB f = ADB::function(std::move(v_copy), std::move(jacs_copy));

    BOOST_REQUIRE(f.value().matrix() == v.matrix());

    const std::vector<ADB::M>& J = f.derivative();
    for (std::vector<ADB::M>::const_iterator
             bf = J.begin(), ef = J.end(), bj = jacs.begin();
         bf != ef; ++bf, ++bj) {

        BOOST_CHECK(*bf == *bj);
    }
}

BOOST_AUTO_TEST_CASE(Addition)
{
    typedef AutoDiffBlock<double> ADB;
    std::vector<int> blocksizes = { 3, 1, 2 };

    ADB::V va(3);
    va << 0.2, 1.2, 13.4;

    ADB::V vx(3);
    vx << 1.0, 2.2, 3.4;

    enum { FirstVar = 0, SecondVar = 1, ThirdVar = 2 };

    ADB a = ADB::constant(va, blocksizes);
    ADB x = ADB::variable(FirstVar, vx, blocksizes);

    ADB xpx = x + x;

    BOOST_CHECK_EQUAL(xpx.value().cwiseNotEqual(2 * x.value()).count(), 0);

    const std::vector<ADB::M>& J1x = x  .derivative();
    const std::vector<ADB::M>& J2x = xpx.derivative();
    BOOST_CHECK_EQUAL(J1x.size(), J2x.size());
    for (std::vector<ADB::M>::const_iterator
             j1b = J1x.begin(), j1e = J1x.end(), j2b = J2x.begin();
         j1b != j1e; ++j1b, ++j2b) {
        BOOST_CHECK(*j2b == ADB::M((*j1b) * 2));
    }

    ADB::V  r = 2*x.value() + a.value();
    ADB xpxpa = x + x + a;
    BOOST_CHECK_EQUAL(xpxpa.value().cwiseNotEqual(r).count(), 0);

    const std::vector<ADB::M>& J3 = xpxpa.derivative();
    for (std::vector<ADB::M>::const_iterator
             j1b = J1x.begin(), j1e = J1x.end(), j3b = J3.begin();
         j1b != j1e; ++j1b, ++j3b) {
        BOOST_CHECK(*j3b == ADB::M((*j1b) * 2));
    }
}

BOOST_AUTO_TEST_CASE(AssignAddSubtractOperators)
{
    typedef AutoDiffBlock<double> ADB;

    // Basic testing of += and -=.
    ADB::V vx(3);
    vx << 0.2, 1.2, 13.4;

    ADB::V vy(3);
    vy << 1.0, 2.2, 3.4;

    std::vector<ADB::V> vals{ vx, vy };
    std::vector<ADB> vars = ADB::variables(vals);

    const ADB x = vars[0];
    const ADB y = vars[1];

    ADB z = x;
    z += y;
    ADB sum = x + y;
    const double tolerance = 1e-14;
    checkClose(z, sum, tolerance);

    z -= y;
    checkClose(z, x, tolerance);

    // Testing the case when the left hand side has empty() jacobian.
    ADB yconst = ADB::constant(vy);
    z = yconst;
    z -= x;
    ADB diff = yconst - x;
    checkClose(z, diff, tolerance);

    z += x;
    checkClose(z, yconst, tolerance);
}

BOOST_AUTO_TEST_CASE(Pow)
{
    typedef AutoDiffBlock<double> ADB;

    // Basic testing of derivatives
    ADB::V vx(3);
    vx << 0.2, 1.2, 13.4;

    ADB::V vy(3);
    vy << 2.0, 3.0, 0.5;

    std::vector<ADB::V> vals{ vx, vy };
    std::vector<ADB> vars = ADB::variables(vals);

    const ADB x = vars[0];
    const ADB y = vars[1];

    const double tolerance = 1e-14;

    // test exp = double
    const ADB xx = x * x;
    ADB xxpow2 = Opm::pow(x,2.0);
    checkClose(xxpow2, xx, tolerance);

    const ADB xy = x * y;
    const ADB xxyy = xy * xy;
    ADB xypow2 = Opm::pow(xy,2.0);
    checkClose(xypow2, xxyy, tolerance);

    const ADB xxx = x * x * x;
    ADB xpow3 = Opm::pow(x,3.0);
    checkClose(xpow3, xxx, tolerance);

    ADB xpowhalf = Opm::pow(x,0.5);
    ADB::V x_sqrt(3);
    x_sqrt << 0.447214 , 1.095445 , 3.6606;
    for (int i = 0 ; i < 3; ++i){
        BOOST_CHECK_CLOSE(xpowhalf.value()[i], x_sqrt[i], 1e-4);
    }

    // test exp = ADB::V
    ADB xpowyval = Opm::pow(x,y.value());

    // each of the component of y is tested in the test above
    // we compare with the results from the above tests.
    ADB::V pick1(3);
    pick1 << 1,0,0;
    ADB::V pick2(3);
    pick2 << 0,1,0;
    ADB::V pick3(3);
    pick3 << 0,0,1;

    ADB compare = pick1 * xx + pick2 * xxx + pick3 * xpowhalf;
    checkClose(xpowyval, compare, tolerance);

    // test exponent = ADB::V and base = ADB
    ADB xvalpowy = Opm::pow(x.value(),y);

    // the value should be equal to xpowyval
    // the first jacobian should be trivial
    // the second jacobian is hand calculated
    // log(0.2)*0.2^2.0, log(1.2) * 1.2^3.0, log(13.4) * 13.4^0.5
    ADB::V jac2(3);
    jac2 << -0.0643775165 , 0.315051650 , 9.50019208855;
    for (int i = 0 ; i < 3; ++i){
        BOOST_CHECK_CLOSE(xvalpowy.value()[i], xpowyval.value()[i], tolerance);
        BOOST_CHECK_CLOSE(xvalpowy.derivative()[0].coeff(i,i), 0.0, tolerance);
        BOOST_CHECK_CLOSE(xvalpowy.derivative()[1].coeff(i,i), jac2[i], 1e-4);
    }

    // test exp = ADB
    ADB xpowy = Opm::pow(x,y);

    // the first jacobian should be equal to the xpowyval
    // the second jacobian should be equal to the xvalpowy
    for (int i = 0 ; i < 3; ++i){
        BOOST_CHECK_CLOSE(xpowy.value()[i], xpowyval.value()[i], tolerance);
        BOOST_CHECK_CLOSE(xpowy.derivative()[0].coeff(i,i), xpowyval.derivative()[0].coeff(i,i), tolerance);
        BOOST_CHECK_CLOSE(xpowy.derivative()[1].coeff(i,i), xvalpowy.derivative()[1].coeff(i,i), tolerance);
    }

}


