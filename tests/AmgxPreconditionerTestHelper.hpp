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

#ifndef TEST_AMGXPRECONDITIONER_HELPER_HPP
#define TEST_AMGXPRECONDITIONER_HELPER_HPP

#include <boost/test/unit_test.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>

#include <opm/simulators/linalg/AmgxPreconditioner.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>

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

//==============================================================================
// Test preconditioner solve CPU input
//==============================================================================
template<typename MatrixScalar, typename VectorScalar>
void testAmgxPreconditionerCpuInput()
{
    constexpr int N = 100; // 100x100 grid
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<MatrixScalar, 1, 1>>;
    using Vector = Dune::BlockVector<Dune::FieldVector<VectorScalar, 1>>;

    // Create matrix
    Matrix matrix;
    setupLaplace2d(N, matrix);

    // Create vectors
    Vector x(N*N), b(N*N);
    x = 0.0;    // Initial guess
    b = 1.0;    // RHS

    // Create operator
    using Operator = Dune::MatrixAdapter<Matrix, Vector, Vector>;
    Operator op(matrix);

    // Create preconditioner
    Opm::PropertyTree prm;
    auto prec = std::make_shared<Amgx::AmgxPreconditioner<Matrix, Vector, Vector>>(matrix, prm);

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

    // Check solution has expected behavior (maximum at the center, decreasing outward)
    double center_value = x[N/2 + (N/2)*N];
    BOOST_CHECK_GT(center_value, 0.0);

    // Check a few points away from center (should have smaller values)
    double near_center = x[(N/2+1) + (N/2)*N];
    double far_from_center = x[0]; // Corner

    BOOST_CHECK_LT(near_center, center_value);
    BOOST_CHECK_LT(far_from_center, near_center);
}

//==============================================================================
// Test preconditioner update CPU input
//==============================================================================
template<typename MatrixScalar, typename VectorScalar>
void testAmgxPreconditionerUpdateCpuInput()
{
    using namespace Opm::gpuistl;

    constexpr int N = 20;
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<MatrixScalar, 1, 1>>;
    using Vector = Dune::BlockVector<Dune::FieldVector<VectorScalar, 1>>;

    // Create matrix
    Matrix matrix;
    setupLaplace2d(N, matrix);

    // Create vectors
    Vector x(N*N);
    Vector b(N*N);
    x = 0.0;
    b = 1.0;

    // Create preconditioner and apply it with initial values
    Opm::PropertyTree prm;
    auto prec = std::make_shared<Amgx::AmgxPreconditioner<Matrix, Vector, Vector>>(matrix, prm);
    prec->apply(x, b);

    // Update matrix with the new values and update the preconditioner
    matrix *= 2.0;
    prec->update();

    // Apply again after update but with same initial guess
    Vector x2(N*N);
    x2 = 0.0;
    prec->apply(x2, b);

    // Verify that the vectors are different
    bool all_same = true;
    for (size_t i = 0; i < x.size(); ++i) {
        if (std::abs(x[i] - x2[i]) > 1e-10) {
            all_same = false;
            break;
        }
    }

    // The vectors should be different after applying preconditioner with updated matrix
    BOOST_CHECK(!all_same);
}


//==============================================================================
// Test preconditioner solve GPU input
//==============================================================================
template<typename MatrixScalar, typename VectorScalar>
void testAmgxPreconditionerGpuInput()
{
    using namespace Opm::gpuistl;

    constexpr int N = 100; // 100x100 grid
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<MatrixScalar, 1, 1>>;
    using Vector = Dune::BlockVector<Dune::FieldVector<VectorScalar, 1>>;
    using GpuMatrixType = GpuSparseMatrix<MatrixScalar>;
    using GpuVectorType = GpuVector<VectorScalar>;

    // Create matrix on CPU first
    Matrix matrix;
    setupLaplace2d(N, matrix);

    // Convert to GPU matrix
    GpuMatrixType gpu_matrix = GpuMatrixType::fromMatrix(matrix);

    // Create cpu vectors
    Vector x(N*N);
    Vector b(N*N);
    x = 0.0;
    b = 1.0;

    // Convert to GPU vectors
    GpuVectorType gpu_x(x);
    GpuVectorType gpu_b(b);

    // Create GPU operator
    using GpuOperator = Dune::MatrixAdapter<GpuMatrixType, GpuVectorType, GpuVectorType>;
    auto op = std::make_shared<GpuOperator>(gpu_matrix);

    // Set up AMGX parameters
    Opm::PropertyTree prm;

    // Create preconditioner
    auto prec = std::make_shared<Amgx::AmgxPreconditioner<GpuMatrixType, GpuVectorType, GpuVectorType>>(gpu_matrix, prm);

    // Create solver
    double reduction = 1e-8;
    int maxit = 300;
    int verbosity = 0;
    Dune::LoopSolver<GpuVectorType> solver(*op, *prec, reduction, maxit, verbosity);

    // Solve the system
    Dune::InverseOperatorResult res;
    solver.apply(gpu_x, gpu_b, res);

    // Check convergence
    BOOST_CHECK(res.converged);
    BOOST_CHECK_LT(res.reduction, reduction);

    // Retrieve solution from GPU
    auto gpu_x_host = gpu_x.asStdVector();

    // Check solution has expected behavior (maximum at the center, decreasing outward)
    double center_value = gpu_x_host[N/2 + (N/2)*N];
    BOOST_CHECK_GT(center_value, 0.0);

    // Check a few points away from center (should have smaller values)
    double near_center = gpu_x_host[(N/2+1) + (N/2)*N];
    double far_from_center = gpu_x_host[0]; // Corner

    BOOST_CHECK_LT(near_center, center_value);
    BOOST_CHECK_LT(far_from_center, near_center);
}


//==============================================================================
// Test preconditioner update GPU input
//==============================================================================
template<typename MatrixScalar, typename VectorScalar>
void testAmgxPreconditionerUpdateGpuInput()
{
    using namespace Opm::gpuistl;

    constexpr int N = 20;
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<MatrixScalar, 1, 1>>;
    using GpuMatrixType = GpuSparseMatrix<MatrixScalar>;
    using GpuVectorType = GpuVector<VectorScalar>;

    // Create matrix
    Matrix matrix;
    setupLaplace2d(N, matrix);
    GpuMatrixType gpu_matrix = GpuMatrixType::fromMatrix(matrix);

    // Create vectors
    GpuVectorType gpu_x(N*N);
    GpuVectorType gpu_b(N*N);
    gpu_x = 0.0;
    gpu_b = 1.0;

    // Create preconditioner and apply it with initial values
    Opm::PropertyTree prm;
    auto prec = std::make_shared<Amgx::AmgxPreconditioner<GpuMatrixType, GpuVectorType, GpuVectorType>>(gpu_matrix, prm);
    prec->apply(gpu_x, gpu_b);

    // Update GPU matrix with the new values and update the preconditioner
    gpu_matrix.getNonZeroValues() *= 2.0;
    prec->update();

    // Apply again after update but with same initial guess
    GpuVectorType gpu_x2(N*N);
    gpu_x2 = 0.0;
    prec->apply(gpu_x2, gpu_b);

    // Verify that the vectors are different
    auto x1 = gpu_x.asStdVector();
    auto x2 = gpu_x2.asStdVector();
    bool all_same = true;
    for (size_t i = 0; i < x1.size(); ++i) {
        if (std::abs(x1[i] - x2[i]) > 1e-10) {
            all_same = false;
            break;
        }
    }

    // The vectors should be different after applying preconditioner with updated matrix
    BOOST_CHECK(!all_same);
}

#endif // TEST_AMGXPRECONDITIONER_HELPER_HPP
