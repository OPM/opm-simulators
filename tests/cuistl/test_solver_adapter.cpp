/*
  Copyright 2022-2023 SINTEF AS

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

#define BOOST_TEST_MODULE TestSolverAdapter

#include <boost/test/unit_test.hpp>
#include <dune/istl/solvers.hh>
#include <opm/simulators/linalg/PreconditionerFactory.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/cuistl/SolverAdapter.hpp>

static const constexpr int dim = 3;
using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, dim, dim>>;
using Vector = Dune::BlockVector<Dune::FieldVector<double, dim>>;
using Moperator = Dune::MatrixAdapter<Matrix, Vector, Vector>;
using PrecondFactory = Opm::PreconditionerFactory<Moperator, Dune::Amg::SequentialInformation>;
using SolverAdapter = Opm::cuistl::SolverAdapter<Moperator, Dune::BiCGSTABSolver, Vector>;

namespace
{
auto
createSolverAdapterWithMatrix(const size_t N = 10)
{


    const int nonZeroes = N * 3 - 2;

    // We need to hold the matrix in memory somehow, but we don't want to deference
    // a pointer all the time (quality of life...):
    auto matrixPtr = std::make_shared<Matrix>(N, N, nonZeroes, Matrix::row_wise);
    auto& matrix = *matrixPtr;
    for (auto row = matrix.createbegin(); row != matrix.createend(); ++row) {
        // Add nonzeros for left neighbour, diagonal and right neighbour
        if (row.index() > 0) {
            row.insert(row.index() - 1);
        }
        row.insert(row.index());
        if (row.index() < matrix.N() - 1) {
            row.insert(row.index() + 1);
        }
    }
    // This might not be the most elegant way of filling in a Dune sparse matrix, but it works.
    for (size_t i = 0; i < N; ++i) {
        for (int k = 0; k < dim; ++k) {
            matrix[i][i][k][k] = -2;
        }
        if (i < N - 1) {
            for (int k = 0; k < dim; ++k) {
                matrix[i][i + 1][k][k] = 1;
            }
        }

        if (i > 0) {
            for (int k = 0; k < dim; ++k) {
                matrix[i][i - 1][k][k] = 1;
            }
        }
    }
    auto op = std::make_shared<Moperator>(matrix);
    auto sp = std::make_shared<Dune::ScalarProduct<Vector>>();
    auto prm = Opm::PropertyTree();
    prm.put<double>("relaxation", 1.0);
    prm.put<std::string>("type", "CUILU0");
    auto prec = PrecondFactory::create(*op, prm);
    auto solverAdapter = std::make_shared<SolverAdapter>(*op, *sp, prec, 1.0, 10, 0);

    return std::make_tuple(matrixPtr, solverAdapter, op, sp);
}
} // namespace

BOOST_AUTO_TEST_CASE(TestCreation)
{
    BOOST_CHECK_NO_THROW(createSolverAdapterWithMatrix(););
}

BOOST_AUTO_TEST_CASE(TestSolve)
{
    const size_t N = 10;
    auto [matrix, solverAdapter, op, sp] = createSolverAdapterWithMatrix(N);

    Vector xActual(N), xInitial(N), b(N);
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            xActual[i][j] = 1.0;
            xInitial[i][j] = 0.1 * i;
        }
    }

    matrix->mv(xActual, b);
    Dune::InverseOperatorResult res;
    solverAdapter->apply(xInitial, b, res);

    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            // This should actually be up to rounding exact since ILU is just the inverse
            // for this matrix.
            BOOST_CHECK_CLOSE(xActual[i][j], xInitial[i][j], 1e-13);
        }
    }
}
