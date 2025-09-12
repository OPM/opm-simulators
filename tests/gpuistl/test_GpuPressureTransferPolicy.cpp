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

#define BOOST_TEST_MODULE TestGpuPressureTransferPolicy

#include <boost/test/unit_test.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>

#include <opm/simulators/linalg/PressureTransferPolicy.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/getQuasiImpesWeights.hpp>

#include <opm/simulators/linalg/gpuistl/GpuPressureTransferPolicy.hpp>
#include <opm/simulators/linalg/gpuistl/GpuSparseMatrix.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>
#include <opm/simulators/linalg/gpuistl/detail/cpr_amg_operations.hpp>

#include <random>
#include <vector>

/**
 * @brief Test fixture that sets up matrices and common test data
 */
template <int blockSize>
struct TestFixture {
    // Define types
    using Scalar = double;
    static constexpr int N = 10;

    using BlockType = Dune::FieldMatrix<Scalar, blockSize, blockSize>;
    using MatrixType = Dune::BCRSMatrix<BlockType>;
    using VectorType = Dune::BlockVector<Dune::FieldVector<Scalar, blockSize>>;
    using CommInfo = Dune::Amg::SequentialInformation;

    // Define GPU types
    using GpuMatrixType = Opm::gpuistl::GpuSparseMatrix<Scalar>;
    using GpuVectorType = Opm::gpuistl::GpuVector<Scalar>;
    using GpuVectorIntType = Opm::gpuistl::GpuVector<int>;

    // Define operator types for the policies
    using CpuOperatorType = Dune::MatrixAdapter<MatrixType, VectorType, VectorType>;
    using GpuOperatorType = Dune::MatrixAdapter<GpuMatrixType, GpuVectorType, GpuVectorType>;

    MatrixType matrix;
    std::mt19937 generator;
    std::uniform_real_distribution<double> distribution;

    TestFixture()
        : matrix(N, N, 3 * N - 2, MatrixType::row_wise)
        , generator()
        , distribution(-10.0, 10.0)
    {
        setupMatrix();
    }

    // Helper method to create GPU matrix when needed
    GpuMatrixType createGpuMatrix() const
    {
        return GpuMatrixType::fromMatrix(matrix);
    }

private:
    void setupMatrix()
    {
        // Setup matrix sparsity pattern (tridiagonal block structure)
        for (auto row = matrix.createbegin(); row != matrix.createend(); ++row) {
            if (row.index() > 0) {
                row.insert(row.index() - 1);
            }
            row.insert(row.index());
            if (row.index() < N - 1) {
                row.insert(row.index() + 1);
            }
        }
        // Fill matrix with random values for testing
        for (int i = 0; i < N; ++i) {
            if (i > 0) {
                for (int row = 0; row < blockSize; ++row) {
                    for (int col = 0; col < blockSize; ++col) {
                        matrix[i][i - 1][row][col] = distribution(generator);
                    }
                }
            }
            for (int row = 0; row < blockSize; ++row) {
                for (int col = 0; col < blockSize; ++col) {
                    matrix[i][i][row][col] = distribution(generator);
                }
            }
            if (i < N - 1) {
                for (int row = 0; row < blockSize; ++row) {
                    for (int col = 0; col < blockSize; ++col) {
                        matrix[i][i + 1][row][col] = distribution(generator);
                    }
                }
            }
        }
    }
};

/**
 * @brief Test quasi-impes weight calculation between CPU and GPU implementations
 */
template <int blockSize, bool transpose>
void
testQuasiImpesWeights(int pressureVarIndex)
{
    TestFixture<blockSize> fixture;

    // Calculate quasiimpes weights for CPU
    auto cpuWeights = Opm::Amg::getQuasiImpesWeights<typename TestFixture<blockSize>::MatrixType,
                                                     typename TestFixture<blockSize>::VectorType>(
        fixture.matrix, pressureVarIndex, transpose);

    // Calculate quasiimpes weights for GPU
    auto gpuMatrix = fixture.createGpuMatrix();
    typename TestFixture<blockSize>::GpuVectorType gpuWeights(fixture.N * blockSize);
    auto diagonalIndices = Opm::Amg::precomputeDiagonalIndices(gpuMatrix);
    typename TestFixture<blockSize>::GpuVectorIntType gpuDiagonalIndices(diagonalIndices);
    Opm::gpuistl::detail::getQuasiImpesWeights<typename TestFixture<blockSize>::Scalar, transpose>(
        gpuMatrix, pressureVarIndex, gpuWeights, gpuDiagonalIndices);

    // Check that both implementations yield the same weights
    std::vector<double> gpuWeightsData = gpuWeights.asStdVector();
    BOOST_REQUIRE_EQUAL(fixture.N * blockSize, gpuWeightsData.size());

    for (int i = 0; i < fixture.N; ++i) {
        for (int j = 0; j < blockSize; ++j) {
            BOOST_CHECK_CLOSE(cpuWeights[i][j], gpuWeightsData[i * blockSize + j], 1e-10);
        }
    }
}

/**
 * @brief Test coarse matrix creation between CPU and GPU implementations
 */
template <int blockSize, bool transpose>
void
testCoarseMatrixCreation(int pressureVarIndex)
{
    TestFixture<blockSize> fixture;

    // Calculate weights for CPU
    auto cpuWeights = Opm::Amg::getQuasiImpesWeights<typename TestFixture<blockSize>::MatrixType,
                                                     typename TestFixture<blockSize>::VectorType>(
        fixture.matrix, pressureVarIndex, transpose);

    // Calculate quasiimpes weights for GPU
    auto gpuMatrix = fixture.createGpuMatrix();
    typename TestFixture<blockSize>::GpuVectorType gpuWeights(fixture.N * blockSize);
    auto diagonalIndices = Opm::Amg::precomputeDiagonalIndices(gpuMatrix);
    typename TestFixture<blockSize>::GpuVectorIntType gpuDiagonalIndices(diagonalIndices);
    Opm::gpuistl::detail::getQuasiImpesWeights<typename TestFixture<blockSize>::Scalar, transpose>(
        gpuMatrix, pressureVarIndex, gpuWeights, gpuDiagonalIndices);

    // Create empty property tree
    Opm::PropertyTree prm;

    // Create CPU transfer policy
    auto cpuPolicy = Opm::PressureTransferPolicy<typename TestFixture<blockSize>::CpuOperatorType,
                                                 typename TestFixture<blockSize>::CommInfo,
                                                 typename TestFixture<blockSize>::Scalar,
                                                 transpose>(
        typename TestFixture<blockSize>::CommInfo(), cpuWeights, prm, pressureVarIndex);

    // Create GPU transfer policy
    auto gpuPolicy = Opm::gpuistl::GpuPressureTransferPolicy<typename TestFixture<blockSize>::GpuOperatorType,
                                                             typename TestFixture<blockSize>::CommInfo,
                                                             typename TestFixture<blockSize>::Scalar,
                                                             transpose>(
        typename TestFixture<blockSize>::CommInfo(), gpuWeights, prm, pressureVarIndex);

    // Create operators
    typename TestFixture<blockSize>::CpuOperatorType cpuOperator(fixture.matrix);
    typename TestFixture<blockSize>::GpuOperatorType gpuOperator(gpuMatrix);

    // Create coarse level systems
    cpuPolicy.createCoarseLevelSystem(cpuOperator);
    gpuPolicy.createCoarseLevelSystem(gpuOperator);

    // Test coarse matrix entries
    const auto& cpuCoarseMatrix = cpuPolicy.getCoarseLevelOperator()->getmat();
    const auto& gpuCoarseMatrix = gpuPolicy.getCoarseLevelOperator()->getmat();

    // Get CPU coarse matrix values
    std::vector<double> cpuCoarseMatrixData;
    for (auto row = cpuCoarseMatrix.begin(); row != cpuCoarseMatrix.end(); ++row) {
        for (auto col = row->begin(); col != row->end(); ++col) {
            cpuCoarseMatrixData.push_back((*col)[0][0]);
        }
    }

    // Get GPU coarse matrix values
    std::vector<double> gpuCoarseMatrixData = gpuCoarseMatrix.getNonZeroValues().asStdVector();

    // Check that coarse matrix entries match
    BOOST_REQUIRE_EQUAL(cpuCoarseMatrixData.size(), gpuCoarseMatrixData.size());
    for (size_t i = 0; i < cpuCoarseMatrixData.size(); ++i) {
        BOOST_CHECK_CLOSE(cpuCoarseMatrixData[i], gpuCoarseMatrixData[i], 1e-10);
    }
}

/**
 * @brief Test restriction operation between CPU and GPU implementations
 */
template <int blockSize, bool transpose>
void
testRestriction(int pressureVarIndex)
{
    TestFixture<blockSize> fixture;

    // Calculate weights for CPU
    auto cpuWeights = Opm::Amg::getQuasiImpesWeights<typename TestFixture<blockSize>::MatrixType,
                                                     typename TestFixture<blockSize>::VectorType>(
        fixture.matrix, pressureVarIndex, transpose);

    // Calculate quasiimpes weights for GPU
    auto gpuMatrix = fixture.createGpuMatrix();
    typename TestFixture<blockSize>::GpuVectorType gpuWeights(fixture.N * blockSize);
    auto diagonalIndices = Opm::Amg::precomputeDiagonalIndices(gpuMatrix);
    typename TestFixture<blockSize>::GpuVectorIntType gpuDiagonalIndices(diagonalIndices);
    Opm::gpuistl::detail::getQuasiImpesWeights<typename TestFixture<blockSize>::Scalar, transpose>(
        gpuMatrix, pressureVarIndex, gpuWeights, gpuDiagonalIndices);

    // Create empty property tree
    Opm::PropertyTree prm;

    // Create CPU transfer policy
    auto cpuPolicy = Opm::PressureTransferPolicy<typename TestFixture<blockSize>::CpuOperatorType,
                                                 typename TestFixture<blockSize>::CommInfo,
                                                 typename TestFixture<blockSize>::Scalar,
                                                 transpose>(
        typename TestFixture<blockSize>::CommInfo(), cpuWeights, prm, pressureVarIndex);

    // Create GPU transfer policy
    auto gpuPolicy = Opm::gpuistl::GpuPressureTransferPolicy<typename TestFixture<blockSize>::GpuOperatorType,
                                                             typename TestFixture<blockSize>::CommInfo,
                                                             typename TestFixture<blockSize>::Scalar,
                                                             transpose>(
        typename TestFixture<blockSize>::CommInfo(), gpuWeights, prm, pressureVarIndex);

    // Create operators
    typename TestFixture<blockSize>::CpuOperatorType cpuOperator(fixture.matrix);
    typename TestFixture<blockSize>::GpuOperatorType gpuOperator(gpuMatrix);

    // Create coarse level systems
    cpuPolicy.createCoarseLevelSystem(cpuOperator);
    gpuPolicy.createCoarseLevelSystem(gpuOperator);

    // Create fine vector
    typename TestFixture<blockSize>::VectorType fineVector(fixture.N);
    for (int i = 0; i < fixture.N; ++i) {
        for (int j = 0; j < blockSize; ++j) {
            fineVector[i][j] = fixture.distribution(fixture.generator);
        }
    }

    // Copy fine vector to GPU
    typename TestFixture<blockSize>::GpuVectorType gpuFineVector(fineVector);

    // Perform restriction
    cpuPolicy.moveToCoarseLevel(fineVector);
    gpuPolicy.moveToCoarseLevel(gpuFineVector);

    // Compare restriction results
    auto& cpuCoarseRhs = cpuPolicy.getCoarseLevelRhs();
    auto& gpuCoarseRhs = gpuPolicy.getCoarseLevelRhs();

    // Check that restriction results match
    std::vector<double> gpuCoarseData = gpuCoarseRhs.asStdVector();
    for (size_t i = 0; i < fixture.N; ++i) {
        BOOST_CHECK_CLOSE(cpuCoarseRhs[i][pressureVarIndex], gpuCoarseData[i], 1e-10);
    }
}

/**
 * @brief Test prolongation operation between CPU and GPU implementations
 */
template <int blockSize, bool transpose>
void
testProlongation(int pressureVarIndex)
{
    TestFixture<blockSize> fixture;

    // Calculate weights and create policies
    auto cpuWeights = Opm::Amg::getQuasiImpesWeights<typename TestFixture<blockSize>::MatrixType,
                                                     typename TestFixture<blockSize>::VectorType>(
        fixture.matrix, pressureVarIndex, transpose);

    // Calculate quasiimpes weights for GPU
    auto gpuMatrix = fixture.createGpuMatrix();
    typename TestFixture<blockSize>::GpuVectorType gpuWeights(fixture.N * blockSize);
    auto diagonalIndices = Opm::Amg::precomputeDiagonalIndices(gpuMatrix);
    typename TestFixture<blockSize>::GpuVectorIntType gpuDiagonalIndices(diagonalIndices);
    Opm::gpuistl::detail::getQuasiImpesWeights<typename TestFixture<blockSize>::Scalar, transpose>(
        gpuMatrix, pressureVarIndex, gpuWeights, gpuDiagonalIndices);

    // Create empty property tree
    Opm::PropertyTree prm;

    // Create CPU transfer policy
    auto cpuPolicy = Opm::PressureTransferPolicy<typename TestFixture<blockSize>::CpuOperatorType,
                                                 typename TestFixture<blockSize>::CommInfo,
                                                 typename TestFixture<blockSize>::Scalar,
                                                 transpose>(
        typename TestFixture<blockSize>::CommInfo(), cpuWeights, prm, pressureVarIndex);

    // Create GPU transfer policy
    auto gpuPolicy = Opm::gpuistl::GpuPressureTransferPolicy<typename TestFixture<blockSize>::GpuOperatorType,
                                                             typename TestFixture<blockSize>::CommInfo,
                                                             typename TestFixture<blockSize>::Scalar,
                                                             transpose>(
        typename TestFixture<blockSize>::CommInfo(), gpuWeights, prm, pressureVarIndex);

    // Create operators
    typename TestFixture<blockSize>::CpuOperatorType cpuOperator(fixture.matrix);
    typename TestFixture<blockSize>::GpuOperatorType gpuOperator(gpuMatrix);

    // Create coarse level systems
    cpuPolicy.createCoarseLevelSystem(cpuOperator);
    gpuPolicy.createCoarseLevelSystem(gpuOperator);

    // Set some values in the coarse lhs for prolongation test
    auto& cpuCoarseLhs = cpuPolicy.getCoarseLevelLhs();
    auto& gpuCoarseLhs = gpuPolicy.getCoarseLevelLhs();
    for (int i = 0; i < fixture.N; ++i) {
        cpuCoarseLhs[i][pressureVarIndex] = fixture.distribution(fixture.generator);
    }

    // Copy coarse lhs to GPU
    gpuCoarseLhs.copyFromHost(cpuCoarseLhs);

    // Create fine result
    typename TestFixture<blockSize>::VectorType cpuFineResult(fixture.N);
    typename TestFixture<blockSize>::GpuVectorType gpuFineResult(fixture.N * blockSize);

    // Initialize fine result to zero, important when not using transpose
    // as we only update the pressure variable
    cpuFineResult = 0.0;
    gpuFineResult = 0.0;

    // Perform prolongation
    cpuPolicy.moveToFineLevel(cpuFineResult);
    gpuPolicy.moveToFineLevel(gpuFineResult);

    // Check that prolongation results match
    std::vector<double> cpuFineResultData(fixture.N * blockSize);
    for (int i = 0; i < fixture.N; ++i) {
        for (int j = 0; j < blockSize; ++j) {
            cpuFineResultData[i * blockSize + j] = cpuFineResult[i][j];
        }
    }
    std::vector<double> gpuFineResultData = gpuFineResult.asStdVector();
    BOOST_REQUIRE_EQUAL(cpuFineResultData.size(), gpuFineResultData.size());
    for (size_t i = 0; i < cpuFineResultData.size(); ++i) {
        BOOST_CHECK_CLOSE(cpuFineResultData[i], gpuFineResultData[i], 1e-10);
    }
}

/**
 * @brief Helper function to run a test function for all valid parameter combinations
 */
template <template <int, bool> class TestFunc>
void
runTestForAllCombinations()
{
    // Test all combinations of block sizes (1-3), pressure variable indices, and transpose flags
    for (int blockSize = 1; blockSize <= 3; ++blockSize) {
        for (int pressureVarIndex = 0; pressureVarIndex < blockSize; ++pressureVarIndex) {
            for (bool transpose : {false, true}) {
                // Dispatch to the appropriate template instantiation
                if (blockSize == 1) {
                    if (transpose) {
                        TestFunc<1, true>::run(pressureVarIndex);
                    } else {
                        TestFunc<1, false>::run(pressureVarIndex);
                    }
                } else if (blockSize == 2) {
                    if (transpose) {
                        TestFunc<2, true>::run(pressureVarIndex);
                    } else {
                        TestFunc<2, false>::run(pressureVarIndex);
                    }
                } else if (blockSize == 3) {
                    if (transpose) {
                        TestFunc<3, true>::run(pressureVarIndex);
                    } else {
                        TestFunc<3, false>::run(pressureVarIndex);
                    }
                }
            }
        }
    }
}

template <int blockSize, bool transpose>
struct QuasiImpesWeightsTestRunner {
    static void run(int pressureVarIndex)
    {
        testQuasiImpesWeights<blockSize, transpose>(pressureVarIndex);
    }
};

template <int blockSize, bool transpose>
struct CoarseMatrixCreationTestRunner {
    static void run(int pressureVarIndex)
    {
        testCoarseMatrixCreation<blockSize, transpose>(pressureVarIndex);
    }
};

template <int blockSize, bool transpose>
struct RestrictionTestRunner {
    static void run(int pressureVarIndex)
    {
        testRestriction<blockSize, transpose>(pressureVarIndex);
    }
};

template <int blockSize, bool transpose>
struct ProlongationTestRunner {
    static void run(int pressureVarIndex)
    {
        testProlongation<blockSize, transpose>(pressureVarIndex);
    }
};

BOOST_AUTO_TEST_CASE(TestQuasiImpesWeights)
{
    runTestForAllCombinations<QuasiImpesWeightsTestRunner>();
}

BOOST_AUTO_TEST_CASE(TestCoarseMatrixCreation)
{
    runTestForAllCombinations<CoarseMatrixCreationTestRunner>();
}

BOOST_AUTO_TEST_CASE(TestRestriction)
{
    runTestForAllCombinations<RestrictionTestRunner>();
}

BOOST_AUTO_TEST_CASE(TestProlongation)
{
    runTestForAllCombinations<ProlongationTestRunner>();
}
