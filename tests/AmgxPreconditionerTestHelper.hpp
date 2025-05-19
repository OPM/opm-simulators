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
#include <opm/simulators/linalg/gpuistl/GpuSparseMatrix.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_type_detection.hpp>


template <class Matrix>
void
setupLaplace2d(int N, Matrix& mat)
{
    const int nonZeroes = N * N * 5; // max 5 entries per row (diagonal + 4 neighbors)
    mat.setBuildMode(Matrix::row_wise);
    mat.setSize(N * N, N * N, nonZeroes);

    // Set up sparsity pattern
    for (auto row = mat.createbegin(); row != mat.createend(); ++row) {
        const int i = row.index();
        int x = i % N;
        int y = i / N;

        row.insert(i); // diagonal
        if (x > 0) // left neighbor
            row.insert(i - 1);
        if (x < N - 1) // right neighbor
            row.insert(i + 1);
        if (y > 0) // upper neighbor
            row.insert(i - N);
        if (y < N - 1) // lower neighbor
            row.insert(i + N);
    }

    // Fill the matrix with values
    for (auto row = mat.begin(); row != mat.end(); ++row) {
        const int i = row.index();
        int x = i % N;
        int y = i / N;

        // Set diagonal
        (*row)[i] = 4.0;

        // Set off-diagonal entries
        if (x > 0) // left neighbor
            (*row)[i - 1] = -1.0;
        if (x < N - 1) // right neighbor
            (*row)[i + 1] = -1.0;
        if (y > 0) // upper neighbor
            (*row)[i - N] = -1.0;
        if (y < N - 1) // lower neighbor
            (*row)[i + N] = -1.0;
    }
}

template <typename MatrixType>
MatrixType
createMatrix(int N)
{
    if constexpr (!Opm::gpuistl::is_gpu_type_v<MatrixType>) {
        MatrixType matrix;
        setupLaplace2d(N, matrix);
        return matrix;
    } else {
        // For GPU matrix, create CPU matrix first then convert
        using ScalarType = typename MatrixType::field_type;
        using CpuMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<ScalarType, 1, 1>>;
        CpuMatrix cpu_matrix;
        setupLaplace2d(N, cpu_matrix);
        return MatrixType::fromMatrix(cpu_matrix);
    }
}


template <typename MatrixType, typename VectorType>
void
testAmgxPreconditionerGeneric()
{
    constexpr int N = 100; // 100x100 grid

    // Create matrix
    auto matrix = createMatrix<MatrixType>(N);

    // Create vectors
    VectorType x(N * N), b(N * N);
    x = 0.0;
    b = 1.0;

    // Create operator
    using Operator = Dune::MatrixAdapter<MatrixType, VectorType, VectorType>;
    auto op = std::make_shared<Operator>(matrix);

    // Create preconditioner
    Opm::PropertyTree prm;
    auto prec = std::make_shared<Amgx::AmgxPreconditioner<MatrixType, VectorType, VectorType>>(matrix, prm);

    // Create solver
    double reduction = 1e-8;
    int maxit = 300;
    int verbosity = 0;
    Dune::LoopSolver<VectorType> solver(*op, *prec, reduction, maxit, verbosity);

    // Solve
    Dune::InverseOperatorResult res;
    solver.apply(x, b, res);

    // Check convergence
    BOOST_CHECK(res.converged);
    BOOST_CHECK_LT(res.reduction, reduction);

    std::vector<typename VectorType::field_type> x_values;
    if constexpr (Opm::gpuistl::is_gpu_type_v<VectorType>) {
        x_values = x.asStdVector();
    } else {
        x_values.assign(&x[0][0], &x[0][0] + x.dim());
    }

    // Validate solution - check that solution has expected spatial distribution
    double center_value = x_values[N / 2 + (N / 2) * N];
    double near_center = x_values[(N / 2 + 1) + (N / 2) * N];
    double far_from_center = x_values[0]; // Corner

    BOOST_CHECK_GT(center_value, 0.0);
    BOOST_CHECK_LT(near_center, center_value);
    BOOST_CHECK_LT(far_from_center, near_center);
}

template <typename MatrixType, typename VectorType>
void
testAmgxPreconditionerUpdateGeneric()
{
    constexpr int N = 20;

    // Create matrix
    auto matrix = createMatrix<MatrixType>(N);

    // Create vectors
    VectorType x(N * N), b(N * N);
    x = 0.0;
    b = 1.0;

    // Create preconditioner and apply it with initial values
    Opm::PropertyTree prm;
    auto prec = std::make_shared<Amgx::AmgxPreconditioner<MatrixType, VectorType, VectorType>>(matrix, prm);
    prec->apply(x, b);

    // Update matrix with new values and update the preconditioner
    if constexpr (!Opm::gpuistl::is_gpu_type_v<MatrixType>) {
        matrix *= 2.0;
    } else {
        matrix.getNonZeroValues() *= 2.0;
    }
    prec->update();

    // Apply again after update but with same initial guess
    VectorType x2(N * N);
    x2 = 0.0;
    prec->apply(x2, b);

    std::vector<typename VectorType::field_type> x_values, x2_values;
    if constexpr (Opm::gpuistl::is_gpu_type_v<VectorType>) {
        x_values = x.asStdVector();
        x2_values = x2.asStdVector();
    } else {
        x_values.assign(&x[0][0], &x[0][0] + x.dim());
        x2_values.assign(&x2[0][0], &x2[0][0] + x2.dim());
    }
    // The vectors should be different after applying preconditioner with updated matrix
    bool all_vectors_different = true;
    for (size_t i = 0; i < x_values.size(); ++i) {
        if (std::abs(x_values[i] - x2_values[i]) <= 1e-10) {
            all_vectors_different = false;
            break;
        }
    }
    BOOST_CHECK(all_vectors_different);
}

template <typename Scalar, bool GpuInput = false>
void
testAmgxPreconditioner()
{
    if constexpr (GpuInput) {
        using Matrix = Opm::gpuistl::GpuSparseMatrix<Scalar>;
        using Vector = Opm::gpuistl::GpuVector<Scalar>;
        testAmgxPreconditionerGeneric<Matrix, Vector>();
    } else {
        using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<Scalar, 1, 1>>;
        using Vector = Dune::BlockVector<Dune::FieldVector<Scalar, 1>>;
        testAmgxPreconditionerGeneric<Matrix, Vector>();
    }
}

template <typename Scalar, bool GpuInput = false>
void
testAmgxPreconditionerUpdate()
{
    if constexpr (GpuInput) {
        using Matrix = Opm::gpuistl::GpuSparseMatrix<Scalar>;
        using Vector = Opm::gpuistl::GpuVector<Scalar>;
        testAmgxPreconditionerUpdateGeneric<Matrix, Vector>();
    } else {
        using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<Scalar, 1, 1>>;
        using Vector = Dune::BlockVector<Dune::FieldVector<Scalar, 1>>;
        testAmgxPreconditionerUpdateGeneric<Matrix, Vector>();
    }
}

#endif // TEST_AMGXPRECONDITIONER_HELPER_HPP
