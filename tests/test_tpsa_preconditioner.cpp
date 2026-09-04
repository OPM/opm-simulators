// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2025 NORCE AS

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

#define BOOST_TEST_MODULE TpsaPreconditionerTest

#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/tpsa/TpsaMatrix.hpp>
#include <opm/simulators/linalg/tpsa/TpsaPreconditioner.hpp>
#include <opm/simulators/linalg/tpsa/TpsaVector.hpp>

#include <array>
#include <set>
#include <string>
#include <vector>

namespace
{

//! \brief Scalar types every test case is instantiated for.
#if FLOW_INSTANTIATE_FLOAT
using ScalarTypes = boost::mpl::list<double, float>;
#else
using ScalarTypes = boost::mpl::list<double>;
#endif // FLOW_INSTANTIATE_FLOAT

template <class Scalar>
using Matrix = Opm::Linear::TpsaMatrix<Scalar>;
template <class Scalar>
using Vector = Opm::Linear::TpsaVector<Scalar>;
template <class Scalar>
using MatrixBlock = typename Matrix<Scalar>::MatrixBlock;
template <class Scalar>
using EqVector = typename Vector<Scalar>::EqVector;
template <class Scalar>
using Preconditioner = Opm::TpsaPreconditioner<Scalar,
                                               Opm::SeqDispDispOperatorT<Scalar>,
                                               Opm::SeqRotRotOperatorT<Scalar>,
                                               Opm::SeqSPresSPresOperatorT<Scalar> >;

constexpr int numEq = Opm::Linear::numTpsaEq;

/*!
 * \brief Relative tolerance, in percent, of the inexact comparisons.
 *
 * \tparam Scalar Field type the comparison is made in.
 */
template <class Scalar>
struct Tolerance;

template <>
struct Tolerance<double>
{
    static constexpr double percent = 1e-8;
    static constexpr double abs = 1e-8;
};
#if FLOW_INSTANTIATE_FLOAT
template <>
struct Tolerance<float>
{
    static constexpr float percent = 1e-2;
    static constexpr float abs = 1e-3f;
};
#endif

/*!
 * \brief A single-cell TPSA matrix, small enough that the block-triangular sweep
 *        in apply() can be checked against a hand-derived reference.
 *
 * \tparam Scalar Field type of the matrix entries.
 *
 * \param[in,out] m Matrix the block is written into.
 */
template <class Scalar>
void
fillSingleCellMatrix(Matrix<Scalar>& m)
{
    const std::vector<std::set<unsigned> > pattern {{0}};
    m.reserve(pattern);

    MatrixBlock<Scalar> b(Scalar(0));
    b[0][0] = Scalar(4);
    b[1][1] = Scalar(5);
    b[2][2] = Scalar(6);

    b[3][3] = Scalar(10);
    b[3][4] = Scalar(1);
    b[3][5] = Scalar(1);
    b[4][3] = Scalar(1);
    b[4][4] = Scalar(10);
    b[4][5] = Scalar(1);
    b[5][3] = Scalar(1);
    b[5][4] = Scalar(1);
    b[5][5] = Scalar(10);

    b[6][6] = Scalar(7);

    // Rotation-displacement coupling, read via S_[_3][_0..2] in apply().
    b[3][0] = Scalar(1);
    b[4][1] = Scalar(1);
    b[5][2] = Scalar(1);

    // Solid pressure-displacement coupling, read via S_[_4][_0..2] in apply().
    b[6][0] = Scalar(0.5);
    b[6][1] = Scalar(0.5);
    b[6][2] = Scalar(0.5);

    m.setBlock(0, 0, b);
}

/*!
 * \brief Simple configuration of the three diagonal block solvers.
 *
 * \param[in] preconditionerType Preconditioner used by each of the three block
 *                               solvers.
 */
template <class Scalar>
Opm::PropertyTree
makeBlockSolverPrm(const std::string& preconditionerType)
{
    using namespace std::string_literals;
    Opm::PropertyTree prm;
    for (const auto& root : {"disp_disp_solver"s, "rot_rot_solver"s, "spres_spres_solver"s}) {
        prm.put(root + ".solver", "loopsolver"s);
        prm.put(root + ".maxiter", 1);
        prm.put(root + ".tol", 0.0);
        prm.put(root + ".verbosity", 0);
        prm.put(root + ".preconditioner.type", preconditionerType);
    }
    return prm;
}

/*!
 * \brief A two-cell TPSA matrix with cell-to-cell coupling.
 *
 * \tparam Scalar Field type of the matrix entries.
 *
 * \param[in,out] m Matrix the blocks are written into.
 */
template <class Scalar>
void
fillTwoCellMatrix(Matrix<Scalar>& m)
{
    const std::vector<std::set<unsigned> > pattern {{0, 1}, {0, 1}};
    m.reserve(pattern);

    MatrixBlock<Scalar> diag(Scalar(0));
    diag[0][0] = Scalar(8);
    diag[1][1] = Scalar(9);
    diag[2][2] = Scalar(10);

    diag[3][3] = Scalar(20);
    diag[3][4] = Scalar(1);
    diag[3][5] = Scalar(1);
    diag[4][3] = Scalar(1);
    diag[4][4] = Scalar(20);
    diag[4][5] = Scalar(1);
    diag[5][3] = Scalar(1);
    diag[5][4] = Scalar(1);
    diag[5][5] = Scalar(20);

    diag[6][6] = Scalar(15);

    diag[3][0] = Scalar(1);
    diag[4][1] = Scalar(1);
    diag[5][2] = Scalar(1);

    diag[6][0] = Scalar(0.5);
    diag[6][1] = Scalar(0.5);
    diag[6][2] = Scalar(0.5);

    MatrixBlock<Scalar> offDiag(Scalar(0));
    offDiag[0][0] = Scalar(2);
    offDiag[1][1] = Scalar(2);
    offDiag[2][2] = Scalar(2);

    offDiag[3][3] = Scalar(2);
    offDiag[3][4] = Scalar(0.3);
    offDiag[3][5] = Scalar(0.3);
    offDiag[4][3] = Scalar(0.3);
    offDiag[4][4] = Scalar(2);
    offDiag[4][5] = Scalar(0.3);
    offDiag[5][3] = Scalar(0.3);
    offDiag[5][4] = Scalar(0.3);
    offDiag[5][5] = Scalar(2);

    offDiag[6][6] = Scalar(1.5);

    m.setBlock(0, 0, diag);
    m.setBlock(1, 1, diag);
    m.setBlock(0, 1, offDiag);
    m.setBlock(1, 0, offDiag);
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE_TEMPLATE(ApplySolvesSingleCellBlockSystem, Scalar, ScalarTypes)
{
    Matrix<Scalar> m(1, 1);
    fillSingleCellMatrix(m);
    const auto prm = makeBlockSolverPrm<Scalar>("ilu0");
    Preconditioner<Scalar> precond(m.istlMatrix(), prm);

    Vector<Scalar> d(1);
    d = Scalar(0);
    EqVector<Scalar> dContribution;
    dContribution[0] = Scalar(8);
    dContribution[1] = Scalar(10);
    dContribution[2] = Scalar(12);
    dContribution[3] = Scalar(3);
    dContribution[4] = Scalar(5);
    dContribution[5] = Scalar(8);
    dContribution[6] = Scalar(5);
    d[0] = dContribution;

    Vector<Scalar> v(1);
    v = Scalar(0);
    precond.apply(v.istlVector(), d.istlVector());

    // Hand-derived reference: with a single cell there is no fill-in, so the
    // block-triangular sweep in apply() is exact.
    //   v0 = d0 / 4 = 2, v1 = d1 / 5 = 2, v2 = d2 / 6 = 2
    //   d3' = d3 - [v0, v1, v2] = [1, 3, 6]
    //   v3 = RotRot^-1 * d3', RotRot = [[10,1,1],[1,10,1],[1,1,10]]
    //   d4' = d4 - 0.5 * (v0 + v1 + v2) = 2
    //   v4 = d4' / 7
    const std::array<Scalar, numEq> expected {
        Scalar(2), Scalar(2), Scalar(2),
        Scalar(1) / Scalar(54), Scalar(13) / Scalar(54), Scalar(31) / Scalar(54),
        Scalar(2) / Scalar(7)};

    const EqVector<Scalar> result = v[0];
    for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
        BOOST_CHECK_CLOSE(result[eqIdx], expected[eqIdx], Tolerance<Scalar>::percent);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(ApplyIsIndependentOfTheIncomingUpdate, Scalar, ScalarTypes)
{
    // A "jac" block solver only ever looks at the diagonal block of whatever it
    // is given, so with cell-to-cell coupling (see fillTwoCellMatrix()) a
    // single loopsolver sweep is an *approximate* solve. That makes its result
    // depend on the initial guess it starts from, hence can be used to test v = 0.0
    // inside apply()
    Matrix<Scalar> m(2, 2);
    fillTwoCellMatrix(m);
    const auto prm = makeBlockSolverPrm<Scalar>("jac");
    Preconditioner<Scalar> precond(m.istlMatrix(), prm);

    Vector<Scalar> d(2);
    d = Scalar(0);
    EqVector<Scalar> d0;
    d0[0] = Scalar(8);
    d0[1] = Scalar(10);
    d0[2] = Scalar(12);
    d0[3] = Scalar(3);
    d0[4] = Scalar(5);
    d0[5] = Scalar(8);
    d0[6] = Scalar(5);
    d[0] = d0;
    EqVector<Scalar> d1;
    d1[0] = Scalar(-4);
    d1[1] = Scalar(6);
    d1[2] = Scalar(-9);
    d1[3] = Scalar(-2);
    d1[4] = Scalar(4);
    d1[5] = Scalar(-6);
    d1[6] = Scalar(-3);
    d[1] = d1;

    const auto applyFrom = [&](Vector<Scalar>& v) {
        precond.apply(v.istlVector(), d.istlVector());
        return std::array<EqVector<Scalar>, 2> {v[0], v[1]};
    };

    // A zero initial guess
    Vector<Scalar> vZero(2);
    vZero = Scalar(0);
    const auto resultZero = applyFrom(vZero);

    // A garbage initial guess, different in each cell
    Vector<Scalar> vGarbage(2);
    EqVector<Scalar> garbage0;
    EqVector<Scalar> garbage1;
    for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
        garbage0[eqIdx] = Scalar(1000 + eqIdx);
        garbage1[eqIdx] = Scalar(-2000 - eqIdx);
    }
    vGarbage[0] = garbage0;
    vGarbage[1] = garbage1;
    const auto resultGarbage = applyFrom(vGarbage);

    for (std::size_t cell = 0; cell < 2; ++cell) {
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            BOOST_CHECK_SMALL(resultZero[cell][eqIdx] - resultGarbage[cell][eqIdx],
                              Tolerance<Scalar>::abs);
        }
    }
}
