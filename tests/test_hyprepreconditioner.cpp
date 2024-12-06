/*
  Copyright 2024 SINTEF AS
  Copyright 2024 Equinor ASA

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

#define BOOST_TEST_MODULE HyprePreconditionerTest
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>

#include <opm/simulators/linalg/HyprePreconditioner.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>

#include "MpiFixture.hpp"

#include <memory>

// Create a 2D Laplace problem matrix
template<class Matrix>
void setupLaplace2d(int N, Matrix& mat)
{
    const int nonZeroes = N*N * 5; // max 5 entries per row (diagonal + 4 neighbors)
    mat.setBuildMode(Matrix::row_wise);
    mat.setSize(N*N, N*N, nonZeroes);

    // Set up sparsity pattern
    for (auto row = mat.createbegin(); row != mat.createend(); ++row) {
        const int i = row.index();
        int x = i % N;
        int y = i / N;

        row.insert(i);  // diagonal
        if (x > 0)      // left neighbor
            row.insert(i-1);
        if (x < N-1)    // right neighbor
            row.insert(i+1);
        if (y > 0)      // upper neighbor
            row.insert(i-N);
        if (y < N-1)    // lower neighbor
            row.insert(i+N);
    }

    // Fill the matrix with values
    for (auto row = mat.begin(); row != mat.end(); ++row) {
        const int i = row.index();
        int x = i % N;
        int y = i / N;

        // Set diagonal
        (*row)[i] = 4.0;

        // Set off-diagonal entries
        if (x > 0)         // left neighbor
            (*row)[i-1] = -1.0;
        if (x < N-1)       // right neighbor
            (*row)[i+1] = -1.0;
        if (y > 0)         // upper neighbor
            (*row)[i-N] = -1.0;
        if (y < N-1)       // lower neighbor
            (*row)[i+N] = -1.0;
    }
}

void testHyprePreconditioner(bool use_gpu)
{
    constexpr int N = 100; // 100x100 grid
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;
    using Vector = Dune::BlockVector<Dune::FieldVector<double, 1>>;

    // Create matrix
    Matrix matrix;
    setupLaplace2d(N, matrix);

    // Create vectors
    Vector x(N*N), b(N*N);
    x = 100.0;  // Initial guess
    b = 0.0;    // RHS

    // Create operator
    using Operator = Dune::MatrixAdapter<Matrix, Vector, Vector>;
    Operator op(matrix);

    // Set up HYPRE parameters
    Opm::PropertyTree prm;
    prm.put("use_gpu", use_gpu ? 1 : 0);

    // Create preconditioner
    auto prec = std::make_shared<Hypre::HyprePreconditioner<Matrix, Vector, Vector>>(matrix, prm);

    // Create solver
    double reduction = 1e-8;
    int maxit = 300;
    int verbosity = 0;
    Dune::LoopSolver<Vector> solver(op, *prec, reduction, maxit, verbosity);

    // Solve
    Dune::InverseOperatorResult res;
    solver.apply(x, b, res);

    // Check convergence
    BOOST_CHECK(res.converged);
    BOOST_CHECK_LT(res.reduction, 1e-8);
}

BOOST_GLOBAL_FIXTURE(MPIFixture);

BOOST_AUTO_TEST_CASE(TestHyprePreconditionerCPU)
{
    testHyprePreconditioner(false);
}

#if HYPRE_USING_CUDA || HYPRE_USING_HIP
BOOST_AUTO_TEST_CASE(TestHyprePreconditionerGPU)
{
    testHyprePreconditioner(true);
}
#endif

bool init_unit_test_func()
{
    return true;
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

#if HYPRE_RELEASE_NUMBER >= 22900
    HYPRE_Initialize();
#else
    HYPRE_Init();
#endif

    int result = boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);

    HYPRE_Finalize();

    return result;
}
