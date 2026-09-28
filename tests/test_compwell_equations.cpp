/*
  Copyright 2026, SINTEF Digital

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
/*!
 * \file
 *
 * \brief Unit tests for the CompWellEquations Schur-complement machinery.
 *
 * CompWellEquations stores the well system in block form
 *
 *     [ A    C^T ] [ x      ]   [ res      ]
 *     [ B    D   ] [ x_well ] = [ res_well ]
 *
 * and eliminates the well unknowns via a Schur complement. This test fills the
 * blocks B, C, D and the well residual with known values and checks the four
 * operations the model relies on against independent dense reference
 * computations (using plain Dune::FieldMatrix algebra):
 *
 *   apply(x, Ax)           ->  Ax -= C^T D^-1 B x
 *   solve(dx)              ->  dx  = D^-1 res_well
 *   recoverSolutionWell    ->  x_w = D^-1 (res_well - B x)
 *   apply(r)               ->  r  -= C^T D^-1 res_well
 *
 * It checks extract(jacobian), A -= C^T D^-1 B, against apply(x, Ax).
 *
 * It also checks that singular matrices are rejected for both the 4x4 block
 * size (which signals singularity by throwing) and the 3x3 block size (which
 * silently yields non-finite values, inf/NaN).
 */
#include "config.h"

#define BOOST_TEST_MODULE CompWellEquations
#include <boost/test/unit_test.hpp>

#include <flowexperimental/comp/wells/CompWellEquations.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace {

using Scalar = double;

// A well system with numWellEq = 4 (matches a three-component well) and a
// reservoir block size numEq = 3, using two connections so the off-diagonal
// blocks B and C carry more than one column.
constexpr int nw = 4;
constexpr int ne = 3;
constexpr int num_conn = 2;

using Eqns = Opm::CompWellEquations<Scalar, nw, ne>;

using DiagMat = Dune::FieldMatrix<Scalar, nw, nw>;
using OffDiagMat = Dune::FieldMatrix<Scalar, nw, ne>;
using WellVec = Dune::FieldVector<Scalar, nw>;
using ResVec = Dune::FieldVector<Scalar, ne>;

// Well-conditioned (diagonally dominant) diagonal block.
DiagMat makeD()
{
    DiagMat d(0.0);
    for (int i = 0; i < nw; ++i) {
        for (int j = 0; j < nw; ++j) {
            d[i][j] = 1.0 / (1.0 + i + j) + (i == j ? 5.0 : 0.0);
        }
    }
    return d;
}

OffDiagMat makeOffDiag(Scalar a, Scalar b)
{
    OffDiagMat m(0.0);
    for (int i = 0; i < nw; ++i) {
        for (int j = 0; j < ne; ++j) {
            m[i][j] = a * (i + 1) - b * j + 0.05 * (i * j);
        }
    }
    return m;
}

void checkClose(const WellVec& a, const WellVec& b, const std::string& what)
{
    for (int i = 0; i < nw; ++i) {
        BOOST_CHECK_MESSAGE(std::abs(a[i] - b[i]) < 1e-10,
                            what << "[" << i << "]: got " << a[i] << " expected " << b[i]);
    }
}

void checkClose(const ResVec& a, const ResVec& b, const std::string& what)
{
    for (int i = 0; i < ne; ++i) {
        BOOST_CHECK_MESSAGE(std::abs(a[i] - b[i]) < 1e-10,
                            what << "[" << i << "]: got " << a[i] << " expected " << b[i]);
    }
}

// Stands in for the reservoir matrix: sums the blocks added per (row, column) cell.
struct BlockCollector {
    using MatrixBlock = Dune::FieldMatrix<Scalar, ne, ne>;

    void addToBlock(std::size_t row, std::size_t col, const MatrixBlock& value)
    {
        blocks.try_emplace({row, col}, 0.0).first->second += value;
    }

    std::map<std::pair<std::size_t, std::size_t>, MatrixBlock> blocks;
};

} // anonymous namespace

BOOST_AUTO_TEST_CASE(SchurComplementOperations)
{
    // --- Reference data (plain dense) -------------------------------------
    const DiagMat Dref = makeD();
    const std::array<OffDiagMat, num_conn> Bref{makeOffDiag(0.10, 0.05), makeOffDiag(0.20, 0.03)};
    const std::array<OffDiagMat, num_conn> Cref{makeOffDiag(0.30, 0.04), makeOffDiag(-0.10, 0.07)};

    WellVec resWell;
    for (int i = 0; i < nw; ++i) {
        resWell[i] = 1.0 + i;
    }

    // Reservoir solution x (one ne-vector per connection).
    std::array<ResVec, num_conn> xres;
    for (int c = 0; c < num_conn; ++c) {
        for (int j = 0; j < ne; ++j) {
            xres[c][j] = 1.0 - 0.2 * j + 0.5 * c;
        }
    }

    // Independent inverse of D for the references.
    DiagMat Dinv = Dref;
    Dinv.invert();

    // --- Fill the CompWellEquations blocks --------------------------------
    Eqns eqns;
    eqns.init(num_conn, std::vector<std::size_t>{0, 1});
    eqns.clear();

    for (int i = 0; i < nw; ++i) {
        eqns.residual()[0][i] = resWell[i];
        for (int j = 0; j < nw; ++j) {
            eqns.D()[0][0][i][j] = Dref[i][j];
        }
        for (int c = 0; c < num_conn; ++c) {
            for (int j = 0; j < ne; ++j) {
                eqns.B()[0][c][i][j] = Bref[c][i][j];
                eqns.C()[0][c][i][j] = Cref[c][i][j];
            }
        }
    }

    eqns.invert();

    // --- apply: Ax -= C^T D^-1 B x -----------------------------------------
    // It reuses the storage of D^-1 res_well as scratch, so it runs first and
    // the checks below also cover the operations that follow it in a solve.
    {
        Eqns::BVector x(num_conn);
        Eqns::BVector Ax(num_conn);
        std::array<ResVec, num_conn> Ax_in;
        for (int c = 0; c < num_conn; ++c) {
            x[c] = xres[c];
            for (int j = 0; j < ne; ++j) {
                Ax_in[c][j] = 0.7 - c + 0.1 * j;
            }
            Ax[c] = Ax_in[c];
        }

        eqns.apply(x, Ax);

        // invDBx = D^-1 sum_c B_c x_c ; Ax_c = Ax_in_c - C_c^T invDBx
        WellVec bx(0.0);
        for (int c = 0; c < num_conn; ++c) {
            Bref[c].umv(xres[c], bx);
        }
        WellVec invDBx;
        Dinv.mv(bx, invDBx);
        for (int c = 0; c < num_conn; ++c) {
            ResVec ctx;
            Cref[c].mtv(invDBx, ctx); // C_c^T invDBx
            ResVec expected = Ax_in[c];
            expected -= ctx;
            checkClose(Ax[c], expected, "apply(x, Ax) conn " + std::to_string(c));
        }
    }

    // --- solve: dx = D^-1 res_well ----------------------------------------
    {
        Eqns::BVectorWell dx(1);
        eqns.solve(dx);

        WellVec expected;
        Dinv.mv(resWell, expected);
        checkClose(dx[0], expected, "solve");
    }

    // --- recoverSolutionWell: x_w = D^-1 (res_well - B x) ------------------
    {
        Eqns::BVector x(num_conn);
        for (int c = 0; c < num_conn; ++c) {
            x[c] = xres[c];
        }

        Eqns::BVectorWell xw(1);
        eqns.recoverSolutionWell(x, xw);

        // res_well - sum_c B_c x_c
        WellVec rhs = resWell;
        for (int c = 0; c < num_conn; ++c) {
            WellVec bx;
            Bref[c].mv(xres[c], bx);
            rhs -= bx;
        }
        WellVec expected;
        Dinv.mv(rhs, expected);
        checkClose(xw[0], expected, "recoverSolutionWell");
    }

    // --- apply: r -= C^T D^-1 res_well -------------------------------------
    {
        Eqns::BVector r(num_conn);
        std::array<ResVec, num_conn> r_in;
        for (int c = 0; c < num_conn; ++c) {
            for (int j = 0; j < ne; ++j) {
                r_in[c][j] = 2.0 + c - 0.3 * j;
            }
            r[c] = r_in[c];
        }

        eqns.apply(r);

        // invDrw = D^-1 res_well ; r_c = r_in_c - C_c^T invDrw
        WellVec invDrw;
        Dinv.mv(resWell, invDrw);
        for (int c = 0; c < num_conn; ++c) {
            ResVec ctx;
            Cref[c].mtv(invDrw, ctx); // C_c^T invDrw
            ResVec expected = r_in[c];
            expected -= ctx;
            checkClose(r[c], expected, "apply conn " + std::to_string(c));
        }
    }
}

// The blocks added to the matrix must act as the well operator, so summing them
// times x gives apply(x, Ax) from Ax = 0. The connection cells are not 0, 1, so
// this also checks that the blocks land in the cells' rows and columns.
BOOST_AUTO_TEST_CASE(MatrixContributionsMatchOperator)
{
    const DiagMat D = makeD();
    const std::array<OffDiagMat, num_conn> B{makeOffDiag(0.10, 0.05), makeOffDiag(0.20, 0.03)};
    const std::array<OffDiagMat, num_conn> C{makeOffDiag(0.30, 0.04), makeOffDiag(-0.10, 0.07)};
    const std::vector<std::size_t> cells{5, 2};

    Eqns eqns;
    eqns.init(num_conn, cells);
    eqns.clear();
    for (int i = 0; i < nw; ++i) {
        for (int j = 0; j < nw; ++j) {
            eqns.D()[0][0][i][j] = D[i][j];
        }
        for (int c = 0; c < num_conn; ++c) {
            for (int j = 0; j < ne; ++j) {
                eqns.B()[0][c][i][j] = B[c][i][j];
                eqns.C()[0][c][i][j] = C[c][i][j];
            }
        }
    }
    eqns.invert();

    BlockCollector jacobian;
    eqns.extract(jacobian);
    BOOST_CHECK_EQUAL(jacobian.blocks.size(), std::size_t{num_conn * num_conn});

    Eqns::BVector x(num_conn);
    for (int c = 0; c < num_conn; ++c) {
        for (int j = 0; j < ne; ++j) {
            x[c][j] = 1.0 - 0.2 * j + 0.5 * c;
        }
    }
    Eqns::BVector Ax(num_conn);
    Ax = 0.0;
    eqns.apply(x, Ax);

    for (int r = 0; r < num_conn; ++r) {
        ResVec y(0.0);
        for (int c = 0; c < num_conn; ++c) {
            const auto it = jacobian.blocks.find({cells[r], cells[c]});
            BOOST_REQUIRE(it != jacobian.blocks.end());
            it->second.umv(x[c], y);
        }
        checkClose(y, Ax[r], "extract conn " + std::to_string(r));
    }
}

// Singular 4x4 well matrix: detail::invertMatrix throws Dune::MatrixBlockError,
// which invert() reports as NumericalProblem.
BOOST_AUTO_TEST_CASE(SingularMatrix4x4)
{
    Eqns eqns;
    eqns.init(num_conn, std::vector<std::size_t>{0, 1});
    eqns.clear(); // D is now all zeros -> singular

    BOOST_CHECK_THROW(eqns.invert(), Opm::NumericalProblem);
}

// Singular 3x3 matrix: detail::invertMatrix silently produces inf/NaN, which
// invert() reports as NumericalProblem too.
BOOST_AUTO_TEST_CASE(SingularMatrix3x3)
{
    constexpr int nw3 = 3;
    using Eqns3 = Opm::CompWellEquations<Scalar, nw3, ne>;

    Eqns3 eqns;
    eqns.init(1, std::vector<std::size_t>{0});
    eqns.clear(); // singular D

    BOOST_CHECK_THROW(eqns.invert(), Opm::NumericalProblem);
}
