/*
  Copyright 2024 SINTEF AS
  Copyright 2024-2025 Equinor ASA

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

#ifndef HYPRE_TEST_HELPER_HPP
#define HYPRE_TEST_HELPER_HPP

#include <boost/test/unit_test.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/owneroverlapcopy.hh>

#include <opm/simulators/linalg/HyprePreconditioner.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/gpuistl/HypreInterface.hpp>

#if HAVE_CUDA
#if USE_HIP
#include <opm/simulators/linalg/gpuistl_hip/GpuSparseMatrixWrapper.hpp>
#include <opm/simulators/linalg/gpuistl_hip/GpuVector.hpp>
#else
#include <opm/simulators/linalg/gpuistl/GpuSparseMatrixWrapper.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>
#endif
#endif

#include <HYPRE.h>
#include <_hypre_utilities.h>

#include <numeric>
#include <vector>

namespace HypreTestHelpers {

/**
 * @brief Per-test fixture for HYPRE state isolation
 *
 * Ensures each test starts with a clean HYPRE state by reinitializing
 * the library before each test and cleaning up afterwards.
 */
#if HYPRE_RELEASE_NUMBER >= 22900
struct HypreTestFixture {
    HypreTestFixture()
    {
        // Reset HYPRE state for each test
        if (HYPRE_Initialized()) {
            HYPRE_Finalize();
        }
        HYPRE_Initialize();
    }
    ~HypreTestFixture()
    {
        // Clean state after test
        if (HYPRE_Initialized()) {
            HYPRE_Finalize();
        }
    }
};
#else
struct HypreTestFixture {
    static bool hypre_is_initialized;

    HypreTestFixture()
    {
        // For older versions, we need to track initialization state manually
        if (hypre_is_initialized) {
            HYPRE_Finalize();
            hypre_is_initialized = false;
        }
        HYPRE_Init();
        hypre_is_initialized = true;
    }
    ~HypreTestFixture()
    {
        // Clean state after test
        if (hypre_is_initialized) {
            HYPRE_Finalize();
            hypre_is_initialized = false;
        }
    }
};
bool HypreTestFixture::hypre_is_initialized = false;
#endif

/**
 * @brief Setup a 2D Laplace operator matrix for testing
 *
 * Creates an N×N grid discretization of the 2D Laplace operator with
 * standard 5-point stencil (diagonal = 4, neighbors = -1).
 *
 * @param N Grid size (N×N total points)
 * @param mat Output matrix to be filled
 */
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

/**
 * @brief Setup matrix helper arrays from CPU matrix for HypreInterface
 *
 * Extracts sparsity pattern information needed for HYPRE matrix operations.
 *
 * @param matrix Input CPU matrix
 * @param ncols Number of columns per row
 * @param rows Row indices
 * @param cols Column indices
 */
template <typename Matrix>
void setupMatrixHelperArrays(const Matrix& matrix,
                            std::vector<HYPRE_Int>& ncols,
                            std::vector<HYPRE_BigInt>& rows,
                            std::vector<HYPRE_BigInt>& cols)
{
    const int N = matrix.N();
    const int nnz = matrix.nonzeroes();

    ncols.resize(N);
    rows.resize(N);
    cols.resize(nnz);

    // Setup sparsity pattern from CPU matrix
    int pos = 0;
    for (auto row = matrix.begin(); row != matrix.end(); ++row) {
        const int rowIdx = row.index();
        rows[rowIdx] = rowIdx;
        ncols[rowIdx] = row->size();

        for (auto col = row->begin(); col != row->end(); ++col) {
            cols[pos++] = col.index();
        }
    }
}

#if HAVE_CUDA || HAVE_HIP
/**
 * @brief Setup matrix helper arrays from GPU matrix for HypreInterface
 *
 * Extracts sparsity pattern information from GPU sparse matrix.
 *
 * @param gpu_matrix Input GPU matrix
 * @param ncols Number of columns per row
 * @param rows Row indices
 * @param cols Column indices
 */
template <typename T>
void setupGpuMatrixHelperArrays(const Opm::gpuistl::GpuSparseMatrixWrapper<T>& gpu_matrix,
                               std::vector<HYPRE_Int>& ncols,
                               std::vector<HYPRE_BigInt>& rows,
                               std::vector<HYPRE_BigInt>& cols)
{
    const int N = gpu_matrix.N();
    const int nnz = gpu_matrix.nonzeroes();

    ncols.resize(N);
    rows.resize(N);
    cols.resize(nnz);

    // Get row pointers and column indices from GPU matrix
    auto host_row_ptrs = gpu_matrix.getRowIndices().asStdVector();
    auto host_col_indices = gpu_matrix.getColumnIndices().asStdVector();

    // Setup row information and ncols
    for (int i = 0; i < N; ++i) {
        rows[i] = i;
        ncols[i] = host_row_ptrs[i + 1] - host_row_ptrs[i];
    }

    // Convert column indices to HYPRE_BigInt format
    for (int i = 0; i < nnz; ++i) {
        cols[i] = host_col_indices[i];
    }
}
#endif

/**
 * @brief Generic test for HyprePreconditioner with any matrix/vector type
 *
 * Tests the preconditioner by solving a Laplace equation with iterative solver.
 * Works with both CPU and GPU matrix/vector types.
 *
 * @param matrix Input matrix (CPU or GPU)
 * @param use_gpu Whether to use GPU backend in HYPRE
 */
template<class MatrixType, class VectorType>
void testHyprePreconditioner(const MatrixType& matrix, bool use_gpu)
{
    constexpr int N = 100; // 100x100 grid

    // Create vectors
    VectorType x(N*N), b(N*N);
    x = 100.0;  // Initial guess
    b = 0.0;    // RHS

    // Create operator
    using Operator = Dune::MatrixAdapter<MatrixType, VectorType, VectorType>;
    Operator op(matrix);

    // Set up HYPRE parameters
    Opm::PropertyTree prm;
    prm.put("use_gpu", use_gpu ? 1 : 0);

    // Create preconditioner
    auto prec = std::make_shared<Hypre::HyprePreconditioner<MatrixType, VectorType, VectorType, Dune::Amg::SequentialInformation>>(matrix, prm, Dune::Amg::SequentialInformation());

    // Create solver
    double reduction = 1e-8;
    int maxit = 300;
    int verbosity = 0;
    Dune::LoopSolver<VectorType> solver(op, *prec, reduction, maxit, verbosity);

    // Solve
    Dune::InverseOperatorResult res;
    solver.apply(x, b, res);

    // Check convergence
    BOOST_CHECK(res.converged);
    BOOST_CHECK_LT(res.reduction, 1e-8);
}

/**
 * @brief Generic vector transfer test for HypreInterface
 *
 * Tests bidirectional vector transfer between Dune vectors and HYPRE vectors.
 *
 * @param use_gpu_backend Whether to use GPU backend
 */
template<typename VectorType>
void testVectorTransfer(bool use_gpu_backend)
{
    using namespace Opm::gpuistl;

    // Initialize HypreInterface
    HypreInterface::initialize(use_gpu_backend);

    // Create test vector with values {1,2,3,4,5} - start with CPU vector
    using CpuVector = Dune::BlockVector<Dune::FieldVector<double, 1>>;
    CpuVector cpu_vec(5);
    for (int i = 0; i < 5; ++i) {
        cpu_vec[i][0] = i + 1.0;
    }

    // Convert to target vector type
    VectorType input_vec = [&]() {
        if constexpr (std::is_same_v<VectorType, Dune::BlockVector<Dune::FieldVector<double, 1>>>) {
            return std::move(cpu_vec);
        }
#if HAVE_CUDA || HAVE_HIP
        else if constexpr (std::is_same_v<VectorType, Opm::gpuistl::GpuVector<double>>) {
            return Opm::gpuistl::GpuVector<double>(cpu_vec);
        }
#endif
    }();

    // Create HYPRE vector
    Dune::Amg::SequentialInformation comm;
    auto hypre_vec = HypreInterface::createVector(5, 0, comm);

    // Setup helper arrays
    HypreInterface::HostArrays host_arrays;
    host_arrays.indices.resize(5);
    std::iota(host_arrays.indices.begin(), host_arrays.indices.end(), 0);

    HypreInterface::DeviceArrays device_arrays;

    // Allocate device arrays if using GPU backend
    if (use_gpu_backend) {
#if HYPRE_USING_CUDA || HYPRE_USING_HIP
        device_arrays.indices_device = hypre_CTAlloc(HYPRE_BigInt, 5, HYPRE_MEMORY_DEVICE);
        device_arrays.vector_buffer_device = hypre_CTAlloc(HYPRE_Real, 5, HYPRE_MEMORY_DEVICE);

        // Copy indices to device
        hypre_TMemcpy(device_arrays.indices_device, host_arrays.indices.data(),
                     HYPRE_BigInt, 5, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
#endif
    }

    // Setup parallel info
    HypreInterface::ParallelInfo parallel_info;
    parallel_info.owner_first = true;
    parallel_info.local_hypre_to_local_dune.clear(); // Empty for owner_first = true

    // Test transfer to HYPRE
    HypreInterface::transferVectorToHypre(input_vec, hypre_vec, host_arrays, device_arrays, parallel_info, use_gpu_backend);

    // Test transfer back from HYPRE
    VectorType result_vec = [&]() {
        if constexpr (std::is_same_v<VectorType, Dune::BlockVector<Dune::FieldVector<double, 1>>>) {
            VectorType vec(5);
            return vec;
        }
#if HAVE_CUDA || HAVE_HIP
        else if constexpr (std::is_same_v<VectorType, Opm::gpuistl::GpuVector<double>>) {
            return Opm::gpuistl::GpuVector<double>(5);
        }
#endif
    }();

    HypreInterface::transferVectorFromHypre(hypre_vec, result_vec, host_arrays, device_arrays, parallel_info, use_gpu_backend);

    // Verify values
    if constexpr (std::is_same_v<VectorType, Dune::BlockVector<Dune::FieldVector<double, 1>>>) {
        for (int i = 0; i < 5; ++i) {
            BOOST_CHECK_CLOSE(input_vec[i][0], result_vec[i][0], 1e-12);
        }
    }
#if HAVE_CUDA || HAVE_HIP
    else if constexpr (std::is_same_v<VectorType, Opm::gpuistl::GpuVector<double>>) {
        auto input_data = input_vec.asStdVector();
        auto result_data = result_vec.asStdVector();
        for (int i = 0; i < 5; ++i) {
            BOOST_CHECK_CLOSE(input_data[i], result_data[i], 1e-12);
        }
    }
#endif

    // Clean up
    if (use_gpu_backend) {
#if HYPRE_USING_CUDA || HYPRE_USING_HIP
        if (device_arrays.indices_device) {
            hypre_TFree(device_arrays.indices_device, HYPRE_MEMORY_DEVICE);
        }
        if (device_arrays.vector_buffer_device) {
            hypre_TFree(device_arrays.vector_buffer_device, HYPRE_MEMORY_DEVICE);
        }
#endif
    }

    HypreInterface::destroyVector(hypre_vec);
}

/**
 * @brief Generic matrix transfer test for HypreInterface
 *
 * Tests matrix value updates in HYPRE matrices from Dune matrices.
 *
 * @param use_gpu_backend Whether to use GPU backend
 */
template<typename MatrixType>
void testMatrixTransfer(bool use_gpu_backend)
{
    using namespace Opm::gpuistl;

    // Initialize HypreInterface
    HypreInterface::initialize(use_gpu_backend);

    // Create test matrix - start with CPU matrix
    using CpuMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;
    CpuMatrix cpu_matrix;
    cpu_matrix.setBuildMode(CpuMatrix::row_wise);
    cpu_matrix.setSize(3, 3, 3);
    for (auto row = cpu_matrix.createbegin(); row != cpu_matrix.createend(); ++row) {
        row.insert(row.index());
    }
    for (auto row = cpu_matrix.begin(); row != cpu_matrix.end(); ++row) {
        for (auto col = row->begin(); col != row->end(); ++col) {
            if (row.index() == col.index()) {
                (*col)[0][0] = 1.0;
            }
        }
    }

    // Convert to target matrix type
    MatrixType matrix = [&]() {
        if constexpr (std::is_same_v<MatrixType, Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>>) {
            return std::move(cpu_matrix);
        }
#if HAVE_CUDA || HAVE_HIP
        else if constexpr (std::is_same_v<MatrixType, Opm::gpuistl::GpuSparseMatrixWrapper<double>>) {
            return Opm::gpuistl::GpuSparseMatrixWrapper<double>::fromMatrix(cpu_matrix, true);
        }
#endif
    }();

    // Setup helper arrays
    std::vector<HYPRE_Int> ncols;
    std::vector<HYPRE_BigInt> rows;
    std::vector<HYPRE_BigInt> cols;

    if constexpr (std::is_same_v<MatrixType, Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>>) {
        setupMatrixHelperArrays(matrix, ncols, rows, cols);
    }
#if HAVE_CUDA || HAVE_HIP
    else if constexpr (std::is_same_v<MatrixType, Opm::gpuistl::GpuSparseMatrixWrapper<double>>) {
        setupGpuMatrixHelperArrays(matrix, ncols, rows, cols);
    }
#endif

    // Create HYPRE matrix
    Dune::Amg::SequentialInformation comm;
    auto hypre_matrix = HypreInterface::createMatrix(3, 0, comm);

    // Setup arrays for HypreInterface
    HypreInterface::SparsityPattern sparsity_pattern;
    sparsity_pattern.ncols = ncols;
    sparsity_pattern.rows = rows;
    sparsity_pattern.cols = cols;
    sparsity_pattern.nnz = cols.size();

    HypreInterface::HostArrays host_arrays;
    std::vector<int> local_dune_to_local_hypre(3);
    for (int i = 0; i < 3; ++i) {
        local_dune_to_local_hypre[i] = i;
    }
    bool owner_first = true;
    host_arrays.row_indexes = HypreInterface::computeRowIndexes(matrix, ncols, local_dune_to_local_hypre, owner_first);

    HypreInterface::DeviceArrays device_arrays;

    // Allocate device arrays if using GPU backend
    if (use_gpu_backend) {
#if HYPRE_USING_CUDA || HYPRE_USING_HIP
        const int N = 3;
        const int nnz = cols.size();

        device_arrays.ncols_device = hypre_CTAlloc(HYPRE_Int, N, HYPRE_MEMORY_DEVICE);
        device_arrays.rows_device = hypre_CTAlloc(HYPRE_BigInt, N, HYPRE_MEMORY_DEVICE);
        device_arrays.cols_device = hypre_CTAlloc(HYPRE_BigInt, nnz, HYPRE_MEMORY_DEVICE);
        device_arrays.row_indexes_device = hypre_CTAlloc(HYPRE_Int, N, HYPRE_MEMORY_DEVICE);
        device_arrays.matrix_buffer_device = hypre_CTAlloc(HYPRE_Real, nnz, HYPRE_MEMORY_DEVICE);

        // Copy data to device
        hypre_TMemcpy(device_arrays.ncols_device, ncols.data(), HYPRE_Int, N, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
        hypre_TMemcpy(device_arrays.rows_device, rows.data(), HYPRE_BigInt, N, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
        hypre_TMemcpy(device_arrays.cols_device, cols.data(), HYPRE_BigInt, nnz, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
        hypre_TMemcpy(device_arrays.row_indexes_device, host_arrays.row_indexes.data(), HYPRE_Int, N, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
#endif
    }

    // Test matrix update
    HypreInterface::updateMatrixValues(matrix, hypre_matrix, sparsity_pattern, host_arrays, device_arrays, use_gpu_backend);

    // Verify values
    auto values_hypre = HypreInterface::getMatrixValues(hypre_matrix, ncols, rows, cols, use_gpu_backend);

    // Check that we got identity matrix values (all 1.0)
    BOOST_REQUIRE_EQUAL(values_hypre.size(), 3);
    for (size_t i = 0; i < values_hypre.size(); ++i) {
        BOOST_CHECK_CLOSE(values_hypre[i], 1.0, 1e-12);
    }

    // Clean up
    if (use_gpu_backend) {
#if HYPRE_USING_CUDA || HYPRE_USING_HIP
        if (device_arrays.ncols_device) hypre_TFree(device_arrays.ncols_device, HYPRE_MEMORY_DEVICE);
        if (device_arrays.rows_device) hypre_TFree(device_arrays.rows_device, HYPRE_MEMORY_DEVICE);
        if (device_arrays.cols_device) hypre_TFree(device_arrays.cols_device, HYPRE_MEMORY_DEVICE);
        if (device_arrays.row_indexes_device) hypre_TFree(device_arrays.row_indexes_device, HYPRE_MEMORY_DEVICE);
        if (device_arrays.matrix_buffer_device) hypre_TFree(device_arrays.matrix_buffer_device, HYPRE_MEMORY_DEVICE);
#endif
    }

    HypreInterface::destroyMatrix(hypre_matrix);
}

/**
 * @brief Test basic HypreInterface resource management
 */
inline void testResourceManagement(bool use_gpu_backend)
{
    using namespace Opm::gpuistl;

    // Initialize HypreInterface
    HypreInterface::initialize(use_gpu_backend);

    // Create and destroy resources
    auto solver = HypreInterface::createAMGSolver();
    Dune::Amg::SequentialInformation comm;
    auto matrix = HypreInterface::createMatrix(10, 0, comm);
    auto vector = HypreInterface::createVector(10, 0, comm);

    // Clean up
    HypreInterface::destroySolver(solver);
    HypreInterface::destroyMatrix(matrix);
    HypreInterface::destroyVector(vector);
}

/**
 * @brief Test error handling in HypreInterface
 */
inline void testErrorHandling(bool use_gpu_backend)
{
    using namespace Opm::gpuistl;

    // Initialize HypreInterface
    HypreInterface::initialize(use_gpu_backend);

    // Test that errors are properly thrown
    Dune::Amg::SequentialInformation comm;
    BOOST_CHECK_THROW(HypreInterface::createMatrix(-1, 0, comm), Opm::gpuistl::HypreError);
    BOOST_CHECK_THROW(HypreInterface::createVector(-1, 0, comm), Opm::gpuistl::HypreError);

    // Test null pointer handling (should not throw)
    HypreInterface::destroySolver(nullptr);
    HypreInterface::destroyMatrix(nullptr);
    HypreInterface::destroyVector(nullptr);
}

} // namespace HypreTestHelpers

#endif // HYPRE_TEST_HELPER_HPP
