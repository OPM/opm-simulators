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
 *   solve(dx)              ->  dx  = D^-1 res_well
 *   recoverSolutionWell    ->  x_w = D^-1 (res_well - B x)
 *   apply(r)               ->  r  -= C^T D^-1 res_well
 *   extract(A)             ->  A  -= C^T D^-1 B
 *
 * It also checks that singular matrices are rejected for both the 4x4 block
 * size (which signals singularity by throwing) and the 3x3 block size (which
 * silently yields non-finite values, inf/NaN).
 *
 * Finally it checks sumAndPinRows(), which keeps the system of a wellbore
 * holding water alone solvable.
 */
#include "config.h"

#define BOOST_TEST_MODULE CompWellEquations
#include <boost/test/unit_test.hpp>

#include <flowexperimental/comp/wells/CompWellEquations.hpp>

#include <opm/common/Exceptions.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
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
using ResMat = Dune::FieldMatrix<Scalar, ne, ne>;

struct TestMatrixAdapter {
    using MatrixBlock = ResMat;
    std::array<std::array<MatrixBlock, 8>, 8> blocks;

    TestMatrixAdapter()
    {
        for (auto& row : blocks) {
            for (auto& block : row) {
                block = 0.0;
            }
        }
    }

    void addToBlock(std::size_t row, std::size_t col, const MatrixBlock& block)
    {
        blocks.at(row).at(col) += block;
    }
};

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

void
checkClose(const ResMat& a, const ResMat& b, const std::string& what)
{
    for (int i = 0; i < ne; ++i) {
        for (int j = 0; j < ne; ++j) {
            BOOST_CHECK_MESSAGE(std::abs(a[i][j] - b[i][j]) < 1e-10,
                                what << "[" << i << "][" << j << "]: got " << a[i][j]
                                     << " expected " << b[i][j]);
        }
    }
}

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
    // The cell indices differ from the connection indices: the vectors handed
    // to apply() and recoverSolutionWell() are local to the well.
    Eqns eqns;
    eqns.init(num_conn, std::vector<std::size_t> {3, 7});
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

    // --- extract: A_rc -= C_r^T D^-1 B_c ----------------------------------
    {
        TestMatrixAdapter jacobian;
        eqns.extract(jacobian);

        constexpr std::array<std::size_t, num_conn> cells {3, 7};
        for (int row = 0; row < num_conn; ++row) {
            for (int col = 0; col < num_conn; ++col) {
                ResMat expected(0.0);
                for (int i = 0; i < ne; ++i) {
                    for (int j = 0; j < ne; ++j) {
                        for (int k = 0; k < nw; ++k) {
                            for (int l = 0; l < nw; ++l) {
                                expected[i][j] -= Cref[row][k][i] * Dinv[k][l] * Bref[col][l][j];
                            }
                        }
                    }
                }
                checkClose(jacobian.blocks[cells[row]][cells[col]],
                           expected,
                           "extract conn " + std::to_string(row) + "," + std::to_string(col));
            }
        }
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

// A wellbore holding water alone: the component rows have lost everything but
// their dependence on the water fraction, so D is singular. sumAndPinRows()
// leaves a system that keeps the mole fractions and determines the remaining
// unknowns from the summed component row, the water row and the control row.
BOOST_AUTO_TEST_CASE(SumAndPinRowsOfWaterFilledWellbore)
{
    // three components and water: total rate, two mole fractions, water
    // fraction, bhp
    constexpr int nw5 = 5;
    constexpr int num_comp = 3;
    constexpr int first_mole_fraction = 1;
    constexpr std::array<int, 3> kept {0, 3, 4}; // unknowns that are still solved for
    using Eqns5 = Opm::CompWellEquations<Scalar, nw5, ne>;
    using Vec5 = Dune::FieldVector<Scalar, nw5>;

    Dune::FieldMatrix<Scalar, nw5, nw5> D(0.0);
    Dune::FieldMatrix<Scalar, nw5, ne> B(0.0);
    Vec5 res(0.0);
    for (int row = 0; row < num_comp; ++row) {
        D[row][3] = -2.0 - row; // only the water fraction is left in a component row
        res[row] = 0.1 * (row + 1);
    }
    for (int row = num_comp; row < nw5; ++row) {
        for (int col = 0; col < nw5; ++col) {
            D[row][col] = 1.0 / (1.0 + row + col) + (row == col ? 3.0 : 0.0);
        }
        res[row] = 1.0 + row;
    }
    for (int row = 0; row < nw5; ++row) {
        for (int col = 0; col < ne; ++col) {
            B[row][col] = 0.1 * (row + 1) - 0.05 * col;
        }
    }

    Eqns5 eqns;
    eqns.init(1, std::vector<std::size_t> {0});
    eqns.clear();
    eqns.D()[0][0] = D;
    eqns.B()[0][0] = B;
    eqns.residual()[0] = res;

    eqns.sumAndPinRows(num_comp - 1, first_mole_fraction);
    eqns.invert();

    // Reference: the summed component row, the water row and the control row
    // in the unknowns that are kept.
    const auto reducedRow = [](const auto& full_row, const int row) {
        auto sum = full_row(row);
        if (row == num_comp - 1) {
            for (int comp = 0; comp < num_comp - 1; ++comp) {
                sum += full_row(comp);
            }
        }
        return sum;
    };
    constexpr std::array<int, 3> rows {num_comp - 1, 3, 4};
    Dune::FieldMatrix<Scalar, 3, 3> Dred(0.0);
    Dune::FieldVector<Scalar, 3> res_red(0.0), rhs_red(0.0);
    ResVec x;
    for (int j = 0; j < ne; ++j) {
        x[j] = 1.0 - 0.3 * j;
    }
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            const int col = kept[j];
            Dred[i][j] = reducedRow([&D, col](const int r) { return D[r][col]; }, rows[i]);
        }
        res_red[i] = reducedRow([&res](const int r) { return res[r]; }, rows[i]);
        rhs_red[i] = res_red[i] - reducedRow([&B, &x](const int r) { return B[r] * x; }, rows[i]);
    }
    Dred.invert();

    const auto check = [&kept](const Vec5& got,
                               const Dune::FieldVector<Scalar, 3>& expected,
                               const std::string& what) {
        for (int comp = 0; comp < num_comp - 1; ++comp) {
            BOOST_CHECK_MESSAGE(got[first_mole_fraction + comp] == 0.0,
                                what << ": mole fraction " << comp << " moved by "
                                     << got[first_mole_fraction + comp]);
        }
        for (int i = 0; i < 3; ++i) {
            BOOST_CHECK_MESSAGE(std::abs(got[kept[i]] - expected[i]) < 1e-10,
                                what << "[" << kept[i] << "]: got " << got[kept[i]] << " expected "
                                     << expected[i]);
        }
    };

    {
        Eqns5::BVectorWell dx(1);
        eqns.solve(dx);
        Dune::FieldVector<Scalar, 3> expected;
        Dred.mv(res_red, expected);
        check(dx[0], expected, "solve");
    }
    {
        Eqns5::BVector xres(1);
        xres[0] = x;
        Eqns5::BVectorWell xw(1);
        eqns.recoverSolutionWell(xres, xw);
        Dune::FieldVector<Scalar, 3> expected;
        Dred.mv(rhs_red, expected);
        check(xw[0], expected, "recoverSolutionWell");
    }
}
