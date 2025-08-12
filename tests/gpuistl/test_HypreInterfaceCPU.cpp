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

#define BOOST_TEST_MODULE TestHypreInterfaceCPU
#define BOOST_TEST_NO_MAIN

#include <dune/common/parallel/mpihelper.hh>
#include "../MpiFixture.hpp"

#include <boost/test/unit_test.hpp>

#include <opm/simulators/linalg/gpuistl/HypreInterface.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>

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

        // Initialize HypreInterface with CPU backend
        HypreInterface::initialize(false);
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

// Global fixture for HypreInterface initialization and finalization
struct HypreGlobalFixture {
    HypreGlobalFixture() {
        // Initialize with CPU backend
        HypreInterface::initialize(false);
    }
    ~HypreGlobalFixture() {
        // No explicit finalize needed as it's handled by main()
    }
};

// Test that type traits work correctly
BOOST_FIXTURE_TEST_CASE(TestTypeTraits, HypreTestFixture)
{
    BOOST_CHECK(!is_gpu_type<std::vector<double>>::value);

    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;
    using Vector = Dune::BlockVector<Dune::FieldVector<double, 1>>;
    BOOST_CHECK(!is_gpu_type<Matrix>::value);
    BOOST_CHECK(!is_gpu_type<Vector>::value);
}

// Test that HypreInterface resource management works
BOOST_FIXTURE_TEST_CASE(TestResourceManagement, HypreTestFixture)
{
    // Create and destroy resources
    auto solver = HypreInterface::createSolver();
    Dune::Amg::SequentialInformation comm;
    auto matrix = HypreInterface::createMatrix(10, 0, comm);
    auto vector = HypreInterface::createVector(10, 0, comm);

    // Clean up
    HypreInterface::destroySolver(solver);
    HypreInterface::destroyMatrix(matrix);
    HypreInterface::destroyVector(vector);

    // Test passes if no exceptions are thrown
    BOOST_CHECK(true);
}

BOOST_FIXTURE_TEST_CASE(TestCpuVectorTransfer, HypreTestFixture)
{
    // Create a simple CPU vector
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

    // Transfer to Hypre (CPU backend)
    std::vector<HYPRE_Real> continuous_vector_values;
    std::vector<int> local_hypre_to_local_dune; // Empty for owner_first = true
    bool owner_first = true;
    bool use_gpu_backend = false;
    HypreInterface::transferVectorToHypre(cpu_vec, hypre_vec, use_gpu_backend, indices, nullptr, nullptr,
                                        continuous_vector_values, local_hypre_to_local_dune, owner_first);

    // Transfer back
    Vector result_vec(5);
    HypreInterface::transferVectorFromHypre(hypre_vec, result_vec, use_gpu_backend, indices, nullptr, nullptr,
                                          continuous_vector_values, local_hypre_to_local_dune, owner_first);

    // Check values
    for (int i = 0; i < 5; ++i) {
        BOOST_CHECK_CLOSE(cpu_vec[i][0], result_vec[i][0], 1e-12);
    }

    HypreInterface::destroyVector(hypre_vec);
}

BOOST_FIXTURE_TEST_CASE(TestCpuMatrixTransfer, HypreTestFixture)
{
    // Create a simple CPU matrix (3x3 identity)
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;
    Matrix cpu_matrix;
    cpu_matrix.setBuildMode(Matrix::row_wise);
    cpu_matrix.setSize(3, 3, 3); // 3x3 matrix with 3 non-zeros (diagonal)

    // Fill matrix structure - make it a simple 3x3 identity matrix
    for (auto row = cpu_matrix.createbegin(); row != cpu_matrix.createend(); ++row) {
        // Insert diagonal element for identity matrix
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

    // Setup helper arrays
    std::vector<HYPRE_Int> ncols;
    std::vector<HYPRE_BigInt> rows;
    std::vector<HYPRE_BigInt> cols;
    setupMatrixHelperArrays(cpu_matrix, ncols, rows, cols);

    // Create Hypre matrix and test initialization
    Dune::Amg::SequentialInformation comm;
    auto hypre_matrix = HypreInterface::createMatrix(3, 0, comm);

    // Setup additional helper arrays for matrix operations
    std::vector<HYPRE_Real> continuous_matrix_entries;
    std::vector<int> local_dune_to_local_hypre(3);
    // Sequential case: identity mapping
    for (int i = 0; i < 3; ++i) {
        local_dune_to_local_hypre[i] = i;
    }
    bool owner_first = true;

    // This should not throw (CPU backend with CPU input)
    HypreInterface::initializeMatrix(cpu_matrix, hypre_matrix, false,
                                   ncols, rows, cols,
                                   nullptr, nullptr, nullptr,
                                   nullptr, continuous_matrix_entries, local_dune_to_local_hypre, owner_first);

    // Test update
    HypreInterface::updateMatrixValues(cpu_matrix, hypre_matrix, false,
                                     ncols, rows, cols,
                                     nullptr, nullptr, nullptr,
                                     nullptr,
                                     continuous_matrix_entries,
                                     local_dune_to_local_hypre,
                                     owner_first);

    // Verify matrix values after update
    auto values_hypre = HypreInterface::getMatrixValues(hypre_matrix, ncols, rows, cols, false);

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

    HypreInterface::destroyMatrix(hypre_matrix);
}

BOOST_FIXTURE_TEST_CASE(TestErrorHandling, HypreTestFixture)
{
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