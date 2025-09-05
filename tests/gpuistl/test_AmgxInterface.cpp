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

#define BOOST_TEST_MODULE TestAmgxInterface
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>

#include <opm/simulators/linalg/gpuistl/AmgxInterface.hpp>
#include <opm/simulators/linalg/gpuistl/GpuSparseMatrixWrapper.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <cstddef>
#include <random>
#include <vector>

using namespace Opm::gpuistl;

// RAII class for managing AMGX resources during tests
class AmgxResourceGuard
{
public:
    AmgxResourceGuard()
    {
        AmgxInterface::initialize();
        cfg_ = AmgxInterface::createConfig("");
        rsrc_ = AmgxInterface::createResources(cfg_);
    }

    ~AmgxResourceGuard()
    {
        AmgxInterface::destroyResources(rsrc_);
        AmgxInterface::destroyConfig(cfg_);
        AmgxInterface::finalize();
    }

    AMGX_config_handle getConfig() const
    {
        return cfg_;
    }
    AMGX_resources_handle getResources() const
    {
        return rsrc_;
    }

private:
    AMGX_config_handle cfg_;
    AMGX_resources_handle rsrc_;
};

// Test vector transfer functionality
BOOST_AUTO_TEST_CASE(TestVectorTransfer)
{
    // Setup AMGX resources
    AmgxResourceGuard resources;

    // Create test data
    const int n = 100;
    std::vector<double> host_data(n);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-100.0, 100.0);

    // Initialize with random data
    for (size_t i = 0; i < n; ++i) {
        host_data[i] = dist(gen);
    }

    // Create GPU vector with test data
    GpuVector<double> gpu_vec(host_data);

    // Create AMGX vector
    AMGX_vector_handle amgx_vec = AmgxInterface::createVector(resources.getResources(), AMGX_mode_dDDI);

    // Set dimensions for the AMGX vector (needed before first transfer)
    OPM_AMGX_SAFE_CALL(AMGX_vector_set_zero(amgx_vec, n, 1));

    // Test transfer GPU -> AMGX
    AmgxInterface::updateAmgxFromGpuVector(gpu_vec, amgx_vec);

    // Create new GPU vector to receive AMGX data
    GpuVector<double> result_vec(n);

    // Test transfer AMGX -> GPU
    AmgxInterface::updateGpuVectorFromAmgx(amgx_vec, result_vec);

    // Verify data integrity
    std::vector<double> result_data = result_vec.asStdVector();
    BOOST_CHECK_EQUAL(result_data.size(), host_data.size());
    for (size_t i = 0; i < n; ++i) {
        BOOST_CHECK_CLOSE(result_data[i], host_data[i], 1e-10);
    }

    // Clean up
    AmgxInterface::destroyVector(amgx_vec);
}

// Test matrix transfer functionality
BOOST_AUTO_TEST_CASE(TestMatrixTransfer)
{
    // Setup AMGX resources
    AmgxResourceGuard resources;

    // Create a simple tridiagonal matrix
    const int N = 10;
    const int blockSize = 1;

    using FieldMatrix = Dune::FieldMatrix<double, blockSize, blockSize>;
    using BCRSMatrix = Dune::BCRSMatrix<FieldMatrix>;

    BCRSMatrix dune_matrix(N, N, BCRSMatrix::row_wise);

    // Set up matrix structure (tridiagonal)
    for (auto row = dune_matrix.createbegin(); row != dune_matrix.createend(); ++row) {
        const int i = row.index();
        if (i > 0) {
            row.insert(i - 1); // Lower diagonal
        }
        row.insert(i); // Diagonal
        if (i < N - 1) {
            row.insert(i + 1); // Upper diagonal
        }
    }

    // Fill with values
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-10.0, 10.0);

    for (int i = 0; i < N; ++i) {
        dune_matrix[i][i] = 10.0 + std::abs(dist(gen)); // Diagonally dominant

        if (i > 0) {
            dune_matrix[i][i - 1] = dist(gen);
        }

        if (i < N - 1) {
            dune_matrix[i][i + 1] = dist(gen);
        }
    }

    // Convert to GpuSparseMatrix
    auto gpu_matrix = GpuSparseMatrixWrapper<double>::fromMatrix(dune_matrix);

    // Create AMGX matrix and upload data
    AMGX_matrix_handle amgx_matrix = AmgxInterface::createMatrix(resources.getResources(), AMGX_mode_dDDI);
    AmgxInterface::updateAmgxMatrixFromGpuSparseMatrix(gpu_matrix, amgx_matrix);

    // Verify matrix dimensions after upload
    int n, block_dimx, block_dimy;
    OPM_AMGX_SAFE_CALL(AMGX_matrix_get_size(amgx_matrix, &n, &block_dimx, &block_dimy));
    BOOST_CHECK_EQUAL(n, N);
    BOOST_CHECK_EQUAL(block_dimx, blockSize);
    BOOST_CHECK_EQUAL(block_dimy, blockSize);

    // Test round-trip: upload to AMGX, then download back
    auto gpu_matrix_copy = GpuSparseMatrixWrapper<double>::fromMatrix(dune_matrix);
    AmgxInterface::updateGpuSparseMatrixFromAmgxMatrix(amgx_matrix, gpu_matrix_copy);

    // Verify round-trip preserves original values
    const auto original_values = gpu_matrix.getNonZeroValues().asStdVector();
    const auto downloaded_values = gpu_matrix_copy.getNonZeroValues().asStdVector();
    BOOST_CHECK_EQUAL(downloaded_values.size(), original_values.size());
    for (std::size_t i = 0; i < original_values.size(); ++i) {
        BOOST_CHECK_CLOSE(downloaded_values[i], original_values[i], 1e-10);
    }

    // Test coefficient update: modify GPU matrix and update AMGX coefficients only
    gpu_matrix.getNonZeroValues() *= 2.0;
    AmgxInterface::updateAmgxMatrixCoefficientsFromGpuSparseMatrix(gpu_matrix, amgx_matrix);

    // Verify coefficient update worked
    AmgxInterface::updateGpuSparseMatrixFromAmgxMatrix(amgx_matrix, gpu_matrix_copy);
    const auto updated_values = gpu_matrix_copy.getNonZeroValues().asStdVector();
    for (std::size_t i = 0; i < original_values.size(); ++i) {
        BOOST_CHECK_CLOSE(updated_values[i], 2.0 * original_values[i], 1e-10);
    }

    // Clean up
    AmgxInterface::destroyMatrix(amgx_matrix, gpu_matrix);
}

// Test templated utility functions
BOOST_AUTO_TEST_CASE(TestTemplatedUtilities)
{
    // Test is_gpu_type trait
    BOOST_CHECK(is_gpu_type<GpuVector<double>>::value);
    BOOST_CHECK(is_gpu_type<GpuSparseMatrixWrapper<double>>::value);
    BOOST_CHECK(!is_gpu_type<std::vector<double>>::value);

    // Test AMGX mode determination
    AMGX_Mode mode_dd = AmgxInterface::determineAmgxMode<double, double>();
    BOOST_CHECK_EQUAL(mode_dd, AMGX_mode_dDDI);

    AMGX_Mode mode_ff = AmgxInterface::determineAmgxMode<float, float>();
    BOOST_CHECK_EQUAL(mode_ff, AMGX_mode_dFFI);
}

// Test error handling
BOOST_AUTO_TEST_CASE(TestErrorHandling)
{
    // Setup AMGX resources
    AmgxResourceGuard resources;

    // Create vectors of different sizes
    GpuVector<double> vec1(10);
    GpuVector<double> vec2(5);

    // Create AMGX vector
    AMGX_vector_handle amgx_vec
        = AmgxInterface::createVector(resources.getResources(), AmgxInterface::determineAmgxMode<double, double>());

    // Set dimensions for the AMGX vector
    OPM_AMGX_SAFE_CALL(AMGX_vector_set_zero(amgx_vec, 10, 1));

    // This should work fine
    AmgxInterface::updateAmgxFromGpuVector(vec1, amgx_vec);

    // This should throw an exception due to size mismatch
    BOOST_CHECK_THROW(AmgxInterface::updateGpuVectorFromAmgx(amgx_vec, vec2), AmgxError);

    // Clean up
    AmgxInterface::destroyVector(amgx_vec);
}

// Custom main function to properly initialize/finalize AMGX
int
main(int argc, char** argv)
{
    return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
