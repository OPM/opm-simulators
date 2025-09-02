/*
  Copyright 2025 Equinor ASA

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

#define BOOST_TEST_MODULE TestHypreInterfaceGPU
#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>

#include "../MpiFixture.hpp"

#include <opm/simulators/linalg/gpuistl/HypreInterface.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>
#include <opm/simulators/linalg/gpuistl/GpuSparseMatrix.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <numeric>

using namespace Opm::gpuistl;

BOOST_GLOBAL_FIXTURE(MPIFixture);

// Per-test fixture for HYPRE state isolation
struct HypreTestFixture {
    HypreTestFixture() {
        // Reset HYPRE state for each test
        if (HYPRE_Initialized()) {
            HYPRE_Finalize();
        }

        // Re-initialize HYPRE for this test
#if HYPRE_RELEASE_NUMBER >= 22900
        HYPRE_Initialize();
#else
        HYPRE_Init();
#endif
    }

    ~HypreTestFixture() {
        // Clean state after test
        if (HYPRE_Initialized()) {
            HYPRE_Finalize();
        }
    }
};

// Helper function to setup matrix helper arrays from CPU matrix
template<typename Matrix>
void setupMatrixHelperArrays(const Matrix& matrix, 
                            std::vector<HYPRE_Int>& ncols,
                            std::vector<HYPRE_BigInt>& rows,
                            std::vector<HYPRE_BigInt>& cols) {
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

// Helper function to setup matrix helper arrays from GPU matrix
template<typename T>
void setupGpuMatrixHelperArrays(const GpuSparseMatrix<T>& gpu_matrix,
                               std::vector<HYPRE_Int>& ncols,
                               std::vector<HYPRE_BigInt>& rows,
                               std::vector<HYPRE_BigInt>& cols) {
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


// Test GPU functionality if available
BOOST_FIXTURE_TEST_CASE(TestGpuBackend, HypreTestFixture)
{
    // Initialize with GPU backend
    HypreInterface::initialize(true);

    // Create resources
    auto solver = HypreInterface::createSolver();
    Dune::Amg::SequentialInformation comm;
    auto matrix = HypreInterface::createMatrix(5, 0, comm);
    auto vector = HypreInterface::createVector(5, 0, comm);

    // Clean up
    HypreInterface::destroySolver(solver);
    HypreInterface::destroyMatrix(matrix);
    HypreInterface::destroyVector(vector);
}

#if HAVE_CUDA || HAVE_HIP
BOOST_FIXTURE_TEST_CASE(TestGpuVectorTransfer, HypreTestFixture)
{
    // Initialize with GPU backend
    HypreInterface::initialize(true);

    // Create a simple GPU vector
    std::vector<double> data = {1.0, 2.0, 3.0, 4.0, 5.0};
    GpuVector<double> gpu_vec(data);

    // Create Hypre vector
    Dune::Amg::SequentialInformation comm;
    auto hypre_vec = HypreInterface::createVector(5, 0, comm);

    // Setup helper arrays for vector operations
    std::vector<HYPRE_BigInt> indices(5);
    std::iota(indices.begin(), indices.end(), 0);

    // Allocate device indices for GPU backend
    HYPRE_BigInt* indices_device = hypre_CTAlloc(HYPRE_BigInt, 5, HYPRE_MEMORY_DEVICE);
    hypre_TMemcpy(indices_device, indices.data(), HYPRE_BigInt, 5, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
    HYPRE_Real* vector_buffer_device = hypre_CTAlloc(HYPRE_Real, 5, HYPRE_MEMORY_DEVICE);

    // Setup MPI parallel helper arrays
    std::vector<HYPRE_Real> continuous_vector_values;
    std::vector<int> local_hypre_to_local_dune; // Empty for owner_first = true
    bool owner_first = true;

    // Transfer to Hypre (GPU backend)
    HypreInterface::transferVectorToHypre(gpu_vec, hypre_vec, true, indices, indices_device, vector_buffer_device,
                                        continuous_vector_values, local_hypre_to_local_dune, owner_first);

    // Transfer back
    GpuVector<double> result_vec(5);
    HypreInterface::transferVectorFromHypre(hypre_vec, result_vec, true, indices, indices_device, vector_buffer_device,
                                          continuous_vector_values, local_hypre_to_local_dune, owner_first);

    // Check values
    auto result_data = result_vec.asStdVector();
    for (int i = 0; i < 5; ++i) {
        BOOST_CHECK_CLOSE(data[i], result_data[i], 1e-12);
    }

    // Clean up
    hypre_TFree(indices_device, HYPRE_MEMORY_DEVICE);
    hypre_TFree(vector_buffer_device, HYPRE_MEMORY_DEVICE);
    HypreInterface::destroyVector(hypre_vec);
}

BOOST_FIXTURE_TEST_CASE(TestGpuMatrixTransfer, HypreTestFixture)
{
    // Initialize with GPU backend
    HypreInterface::initialize(true);

    // Create a simple CPU matrix first
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;
    Matrix cpu_matrix;
    cpu_matrix.setBuildMode(Matrix::row_wise);
    cpu_matrix.setSize(3, 3, 3); // 3x3 matrix with 3 non-zeros (diagonal)

    // Fill matrix structure
    for (auto row = cpu_matrix.createbegin(); row != cpu_matrix.createend(); ++row) {
        row.insert(row.index());
    }

    // Set diagonal values
    for (auto row = cpu_matrix.begin(); row != cpu_matrix.end(); ++row) {
        for (auto col = row->begin(); col != row->end(); ++col) {
            if (row.index() == col.index()) {
                (*col)[0][0] = 1.0;
            }
        }
    }

    // Convert to GPU matrix
    auto gpu_matrix = GpuSparseMatrix<double>::fromMatrix(cpu_matrix, true);

    // Setup helper arrays for GPU matrix
    std::vector<HYPRE_Int> ncols;
    std::vector<HYPRE_BigInt> rows;
    std::vector<HYPRE_BigInt> cols;
    setupGpuMatrixHelperArrays(gpu_matrix, ncols, rows, cols);

    // Pre-compute row indexes for optimal performance
    std::vector<int> local_dune_to_local_hypre(3);
    for (int i = 0; i < 3; ++i) {
        local_dune_to_local_hypre[i] = i;
    }
    bool owner_first = true;
    auto row_indexes = HypreInterface::computeRowIndexes(gpu_matrix, ncols, local_dune_to_local_hypre, owner_first);

    // Allocate device arrays for GPU backend
    const int N = gpu_matrix.N();
    const int nnz = gpu_matrix.nonzeroes();
    HYPRE_Int* ncols_device = hypre_CTAlloc(HYPRE_Int, N, HYPRE_MEMORY_DEVICE);
    HYPRE_BigInt* rows_device = hypre_CTAlloc(HYPRE_BigInt, N, HYPRE_MEMORY_DEVICE);
    HYPRE_BigInt* cols_device = hypre_CTAlloc(HYPRE_BigInt, nnz, HYPRE_MEMORY_DEVICE);
    HYPRE_Real* matrix_buffer_device = hypre_CTAlloc(HYPRE_Real, nnz, HYPRE_MEMORY_DEVICE);

    // Copy data to device
    hypre_TMemcpy(ncols_device, ncols.data(), HYPRE_Int, N, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
    hypre_TMemcpy(rows_device, rows.data(), HYPRE_BigInt, N, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
    hypre_TMemcpy(cols_device, cols.data(), HYPRE_BigInt, nnz, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);

    // Create Hypre matrix and test initialization
    Dune::Amg::SequentialInformation comm;
    auto hypre_matrix = HypreInterface::createMatrix(3, 0, comm);

    // Setup MPI parallel helper arrays

    // This should not throw (GPU backend with GPU input)
    HypreInterface::initializeMatrix(gpu_matrix, hypre_matrix, true,
                                    ncols, rows, cols, row_indexes,
                                    ncols_device, rows_device, cols_device, nullptr,
                                    nullptr);

    // Test update
    HypreInterface::updateMatrixValues(gpu_matrix, hypre_matrix, true,
                                        ncols, rows, cols, row_indexes,
                                        ncols_device, rows_device, cols_device, nullptr,
                                        nullptr);

    // Verify matrix values after update
    auto values_hypre = HypreInterface::getMatrixValues(hypre_matrix, ncols, rows, cols, true);
    auto original_values = gpu_matrix.getNonZeroValues().asStdVector();
    for (size_t i = 0; i < values_hypre.size(); ++i) {
        BOOST_CHECK_CLOSE(values_hypre[i], original_values[i], 1e-12);
    }

    // Clean up device arrays
    hypre_TFree(ncols_device, HYPRE_MEMORY_DEVICE);
    hypre_TFree(rows_device, HYPRE_MEMORY_DEVICE);
    hypre_TFree(cols_device, HYPRE_MEMORY_DEVICE);
    hypre_TFree(matrix_buffer_device, HYPRE_MEMORY_DEVICE);
    HypreInterface::destroyMatrix(hypre_matrix);
}

BOOST_FIXTURE_TEST_CASE(TestCpuVectorWithGpuBackend, HypreTestFixture)
{
    // Initialize with GPU backend
    HypreInterface::initialize(true);

    // Test CPU vector transfer with GPU backend
    using Vector = Dune::BlockVector<Dune::FieldVector<double, 1>>;
    Vector cpu_vec(5);
    for (int i = 0; i < 5; ++i) {
        cpu_vec[i][0] = i + 1.0;
    }

    // Create Hypre vector
    Dune::Amg::SequentialInformation comm;
    auto hypre_vec = HypreInterface::createVector(5, 0, comm);

    // Setup helper arrays for vector operations
    std::vector<HYPRE_BigInt> indices(5);
    std::iota(indices.begin(), indices.end(), 0);

    // Allocate device indices for GPU backend
    HYPRE_BigInt* indices_device = hypre_CTAlloc(HYPRE_BigInt, 5, HYPRE_MEMORY_DEVICE);
    hypre_TMemcpy(indices_device, indices.data(), HYPRE_BigInt, 5, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
    HYPRE_Real* vector_buffer_device = hypre_CTAlloc(HYPRE_Real, 5, HYPRE_MEMORY_DEVICE);

    // Setup MPI parallel helper arrays
    std::vector<HYPRE_Real> continuous_vector_values;
    std::vector<int> local_hypre_to_local_dune; // Empty for owner_first = true
    bool owner_first = true;

    // Transfer to Hypre (GPU backend with CPU input)
    HypreInterface::transferVectorToHypre(cpu_vec, hypre_vec, true, indices, indices_device, vector_buffer_device,
                                        continuous_vector_values, local_hypre_to_local_dune, owner_first);

    // Transfer back
    Vector result_vec(5);
    HypreInterface::transferVectorFromHypre(hypre_vec, result_vec, true, indices, indices_device, vector_buffer_device,
                                          continuous_vector_values, local_hypre_to_local_dune, owner_first);

    // Check values
    for (int i = 0; i < 5; ++i) {
        BOOST_CHECK_CLOSE(cpu_vec[i][0], result_vec[i][0], 1e-12);
    }

    // Clean up
    hypre_TFree(indices_device, HYPRE_MEMORY_DEVICE);
    hypre_TFree(vector_buffer_device, HYPRE_MEMORY_DEVICE);
    HypreInterface::destroyVector(hypre_vec);
}

BOOST_FIXTURE_TEST_CASE(TestCpuMatrixWithGpuBackend, HypreTestFixture)
{
    // Initialize with GPU backend
    HypreInterface::initialize(true);

    // Test CPU matrix transfer with GPU backend
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;
    Matrix cpu_matrix;
    cpu_matrix.setBuildMode(Matrix::row_wise);
    cpu_matrix.setSize(3, 3, 3); // 3x3 matrix with 3 non-zeros (diagonal)

    // Fill matrix structure - make it a simple 3x3 identity matrix
    for (auto row = cpu_matrix.createbegin(); row != cpu_matrix.createend(); ++row) {
        row.insert(row.index());
    }

    // Set diagonal values
    for (auto row = cpu_matrix.begin(); row != cpu_matrix.end(); ++row) {
        for (auto col = row->begin(); col != row->end(); ++col) {
            if (row.index() == col.index()) {
                (*col)[0][0] = 1.0;
            }
        }
    }

    // Setup helper arrays for CPU matrix
    std::vector<HYPRE_Int> ncols;
    std::vector<HYPRE_BigInt> rows;
    std::vector<HYPRE_BigInt> cols;
    setupMatrixHelperArrays(cpu_matrix, ncols, rows, cols);

    // Setup MPI parallel helper arrays
    std::vector<int> local_dune_to_local_hypre(3);  // Sequential case: identity mapping
    for (int i = 0; i < 3; ++i) {
        local_dune_to_local_hypre[i] = i;
    }
    bool owner_first = true;

    // Pre-compute row indexes for optimal performance
    auto row_indexes = HypreInterface::computeRowIndexes(cpu_matrix, ncols, local_dune_to_local_hypre, owner_first);

    // Allocate device arrays for GPU backend
    const int N = cpu_matrix.N();
    const int nnz = cpu_matrix.nonzeroes();
    HYPRE_Int* ncols_device = hypre_CTAlloc(HYPRE_Int, N, HYPRE_MEMORY_DEVICE);
    HYPRE_BigInt* rows_device = hypre_CTAlloc(HYPRE_BigInt, N, HYPRE_MEMORY_DEVICE);
    HYPRE_BigInt* cols_device = hypre_CTAlloc(HYPRE_BigInt, nnz, HYPRE_MEMORY_DEVICE);
    HYPRE_Int* row_indexes_device = hypre_CTAlloc(HYPRE_Int, N, HYPRE_MEMORY_DEVICE);
    HYPRE_Real* matrix_buffer_device = hypre_CTAlloc(HYPRE_Real, nnz, HYPRE_MEMORY_DEVICE);

    // Copy data to device
    hypre_TMemcpy(ncols_device, ncols.data(), HYPRE_Int, N, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
    hypre_TMemcpy(rows_device, rows.data(), HYPRE_BigInt, N, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
    hypre_TMemcpy(cols_device, cols.data(), HYPRE_BigInt, nnz, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
    hypre_TMemcpy(row_indexes_device, row_indexes.data(), HYPRE_Int, N, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);

    // Create Hypre matrix and test initialization
    Dune::Amg::SequentialInformation comm;
    auto hypre_matrix = HypreInterface::createMatrix(3, 0, comm);

    // Setup MPI parallel helper arrays
    std::vector<HYPRE_Real> continuous_matrix_entries;

    // This should not throw (GPU backend with CPU input)
    HypreInterface::initializeMatrix(cpu_matrix, hypre_matrix, true,
                                    ncols, rows, cols, row_indexes,
                                    ncols_device, rows_device, cols_device, row_indexes_device,
                                    matrix_buffer_device);

    // Test update
    HypreInterface::updateMatrixValues(cpu_matrix, hypre_matrix, true,
                                        ncols, rows, cols, row_indexes,
                                        ncols_device, rows_device, cols_device, row_indexes_device,
                                        matrix_buffer_device);

    // Verify matrix values after update
    auto values_hypre = HypreInterface::getMatrixValues(hypre_matrix, ncols, rows, cols, true);

    // Create expected values from CPU matrix for comparison
    std::vector<double> expected_values;
    for (auto row = cpu_matrix.begin(); row != cpu_matrix.end(); ++row) {
        for (auto col = row->begin(); col != row->end(); ++col) {
            expected_values.push_back((*col)[0][0]);
        }
    }

    BOOST_REQUIRE_EQUAL(values_hypre.size(), expected_values.size());
    for (size_t i = 0; i < values_hypre.size(); ++i) {
        BOOST_CHECK_CLOSE(values_hypre[i], expected_values[i], 1e-12);
    }

    // Clean up device arrays
    hypre_TFree(ncols_device, HYPRE_MEMORY_DEVICE);
    hypre_TFree(rows_device, HYPRE_MEMORY_DEVICE);
    hypre_TFree(cols_device, HYPRE_MEMORY_DEVICE);
    hypre_TFree(row_indexes_device, HYPRE_MEMORY_DEVICE);
    hypre_TFree(matrix_buffer_device, HYPRE_MEMORY_DEVICE);
    HypreInterface::destroyMatrix(hypre_matrix);
}

BOOST_FIXTURE_TEST_CASE(TestGpuVectorWithCpuBackend, HypreTestFixture)
{
    // Initialize with CPU backend
    HypreInterface::initialize(false);

    // Test GPU vector transfer with CPU backend
    std::vector<double> data = {1.0, 2.0, 3.0, 4.0, 5.0};
    GpuVector<double> gpu_vec(data);

    // Create Hypre vector
    Dune::Amg::SequentialInformation comm;
    auto hypre_vec = HypreInterface::createVector(5, 0, comm);

    // Setup helper arrays for vector operations
    std::vector<HYPRE_BigInt> indices(5);
    std::iota(indices.begin(), indices.end(), 0);

    // For CPU backend, we don't need device arrays, but we still pass nullptr for consistency
    HYPRE_BigInt* indices_device = nullptr;
    HYPRE_Real* vector_buffer_device = nullptr;

    // Setup MPI parallel helper arrays
    std::vector<HYPRE_Real> continuous_vector_values;
    std::vector<int> local_hypre_to_local_dune; // Empty for owner_first = true
    bool owner_first = true;

    // Transfer to Hypre (CPU backend with GPU input)
    HypreInterface::transferVectorToHypre(gpu_vec, hypre_vec, false, indices, indices_device, vector_buffer_device,
                                        continuous_vector_values, local_hypre_to_local_dune, owner_first);

    // Transfer back
    GpuVector<double> result_vec(5);
    HypreInterface::transferVectorFromHypre(hypre_vec, result_vec, false, indices, indices_device, vector_buffer_device,
                                          continuous_vector_values, local_hypre_to_local_dune, owner_first);

    // Check values
    auto result_data = result_vec.asStdVector();
    for (int i = 0; i < 5; ++i) {
        BOOST_CHECK_CLOSE(data[i], result_data[i], 1e-12);
    }

    // Clean up
    HypreInterface::destroyVector(hypre_vec);
}

BOOST_FIXTURE_TEST_CASE(TestGpuMatrixWithCpuBackend, HypreTestFixture)
{
    // Initialize with CPU backend
    HypreInterface::initialize(false);

    // Test GPU matrix transfer with CPU backend
    // Create a simple CPU matrix first
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;
    Matrix cpu_matrix;
    cpu_matrix.setBuildMode(Matrix::row_wise);
    cpu_matrix.setSize(3, 3, 3); // 3x3 matrix with 3 non-zeros (diagonal)

    // Fill matrix structure
    for (auto row = cpu_matrix.createbegin(); row != cpu_matrix.createend(); ++row) {
        row.insert(row.index());
    }

    // Set diagonal values
    for (auto row = cpu_matrix.begin(); row != cpu_matrix.end(); ++row) {
        for (auto col = row->begin(); col != row->end(); ++col) {
            if (row.index() == col.index()) {
                (*col)[0][0] = 1.0;
            }
        }
    }

    // Convert to GPU matrix
    auto gpu_matrix = GpuSparseMatrix<double>::fromMatrix(cpu_matrix, true);

    // Setup helper arrays for GPU matrix
    std::vector<HYPRE_Int> ncols;
    std::vector<HYPRE_BigInt> rows;
    std::vector<HYPRE_BigInt> cols;
    setupGpuMatrixHelperArrays(gpu_matrix, ncols, rows, cols);

    // For CPU backend, we don't need device arrays, but we still pass nullptr for consistency
    HYPRE_Int* ncols_device = nullptr;
    HYPRE_BigInt* rows_device = nullptr;
    HYPRE_BigInt* cols_device = nullptr;

    // Create Hypre matrix and test initialization
    Dune::Amg::SequentialInformation comm;
    auto hypre_matrix = HypreInterface::createMatrix(3, 0, comm);

    // Setup MPI parallel helper arrays
    std::vector<HYPRE_Real> continuous_matrix_entries;
    std::vector<int> local_dune_to_local_hypre(3);  // empty for owner_first = true
    bool owner_first = true;

    // Pre-compute row indexes for optimal performance
    auto row_indexes = HypreInterface::computeRowIndexes(gpu_matrix, ncols, local_dune_to_local_hypre, owner_first);

    // This should not throw (CPU backend with GPU input)
    HypreInterface::initializeMatrix(gpu_matrix, hypre_matrix, false,
                                    ncols, rows, cols, row_indexes,
                                    ncols_device, rows_device, cols_device, nullptr,
                                    nullptr);

    // Test update
    HypreInterface::updateMatrixValues(gpu_matrix, hypre_matrix, false,
                                        ncols, rows, cols, row_indexes,
                                        ncols_device, rows_device, cols_device, nullptr,
                                        nullptr);

    // get the matrix values using the new interface function
    auto values_hypre = HypreInterface::getMatrixValues(hypre_matrix, ncols, rows, cols, false);

    // check the values (diagonal matrix should have 1.0 values)
    auto original_values = gpu_matrix.getNonZeroValues().asStdVector();
    for (size_t i = 0; i < values_hypre.size(); ++i) {
        BOOST_CHECK_CLOSE(values_hypre[i], original_values[i], 1e-12);
    }

    // Clean up
    HypreInterface::destroyMatrix(hypre_matrix);
}
#endif

BOOST_FIXTURE_TEST_CASE(TestErrorHandling, HypreTestFixture)
{
    // Initialize with GPU backend
    HypreInterface::initialize(true);

    // Test that errors are properly thrown
    Dune::Amg::SequentialInformation comm;
    BOOST_CHECK_THROW(HypreInterface::createMatrix(-1, 0, comm), HypreError);
    BOOST_CHECK_THROW(HypreInterface::createVector(-1, 0, comm), HypreError);

    // Test null pointer handling
    HypreInterface::destroySolver(nullptr); // Should not throw
    HypreInterface::destroyMatrix(nullptr); // Should not throw
    HypreInterface::destroyVector(nullptr); // Should not throw
}

bool init_unit_test_func()
{
    return true;
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    int result = boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);

    return result;
}