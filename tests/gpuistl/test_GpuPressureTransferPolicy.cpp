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
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>

#include <opm/simulators/linalg/PressureTransferPolicy.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>

#include <opm/simulators/linalg/gpuistl/GpuPressureTransferPolicy.hpp>
#include <opm/simulators/linalg/gpuistl/GpuSparseMatrix.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>

#include <random>
#include <vector>

// Helper function to run tests for a specific transpose mode and pressure variable index
template <bool transpose>
void runPressureTransferPolicyTest(int pressureVarIndex)
{
    // Define types
    using Scalar = double;
    constexpr int blockSize = 3;
    constexpr int N = 10; // Size of the system

    using BlockType = Dune::FieldMatrix<Scalar, blockSize, blockSize>;
    using MatrixType = Dune::BCRSMatrix<BlockType>;
    using VectorType = Dune::BlockVector<Dune::FieldVector<Scalar, blockSize>>;
    using CommInfo = Dune::Amg::SequentialInformation;

    // Define GPU types
    using GpuMatrixType = Opm::gpuistl::GpuSparseMatrix<Scalar>;
    using GpuVectorType = Opm::gpuistl::GpuVector<Scalar>;

    // Define operator types for the policies
    using CpuOperatorType = Dune::MatrixAdapter<MatrixType, VectorType, VectorType>;
    using GpuOperatorType = Dune::MatrixAdapter<GpuMatrixType, GpuVectorType, GpuVectorType>;

    // Create a simple test matrix
    MatrixType matrix(N, N, 3*N-2, MatrixType::row_wise);

    // Setup matrix sparsity pattern
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
    std::mt19937 generator(123); // Fixed seed for reproducibility
    std::uniform_real_distribution<double> distribution(-10.0, 10.0);

    for (int i = 0; i < N; ++i) {
        if (i > 0) {
            for (int row = 0; row < blockSize; ++row) {
                for (int col = 0; col < blockSize; ++col) {
                    matrix[i][i-1][row][col] = distribution(generator);
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
                    matrix[i][i+1][row][col] = distribution(generator);
                }
            }
        }
    }

    // Convert to GPU matrix
    GpuMatrixType gpuMatrix = GpuMatrixType::fromMatrix(matrix);

    // Create operators
    CpuOperatorType cpuOperator(matrix);
    GpuOperatorType gpuOperator(gpuMatrix);

    // Create weights vector (typically pressure weights)
    VectorType weights(N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < blockSize; ++j) {
            weights[i][j] = (j == pressureVarIndex) ? 1.0 : 0.0; // Only pressure component has weight 1
        }
    }

    // Create GPU weights
    GpuVectorType gpuWeights(weights);

    // create empty property tree
    Opm::PropertyTree prm;

    // Create policies
    auto cpuPolicy = Opm::PressureTransferPolicy<CpuOperatorType, CommInfo, Scalar, transpose>(
        CommInfo(), weights, prm, pressureVarIndex);

    auto gpuPolicy = Opm::gpuistl::GpuPressureTransferPolicy<GpuOperatorType, CommInfo, Scalar, transpose>(
        CommInfo(), gpuWeights, prm, pressureVarIndex);

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
    BOOST_CHECK_EQUAL(cpuCoarseMatrixData.size(), gpuCoarseMatrixData.size());
    for (size_t i = 0; i < cpuCoarseMatrixData.size(); ++i) {
        BOOST_CHECK_CLOSE(cpuCoarseMatrixData[i], gpuCoarseMatrixData[i], 1e-10);
    }

    // Test restriction (moveToCoarseLevel)
    VectorType fineVector(N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < blockSize; ++j) {
            fineVector[i][j] = distribution(generator);
        }
    }

    // Copy fine vector to GPU
    GpuVectorType gpuFineVector(fineVector);

    // Perform restriction
    cpuPolicy.moveToCoarseLevel(fineVector);
    gpuPolicy.moveToCoarseLevel(gpuFineVector);

    // Compare restriction results
    auto& cpuCoarseRhs = cpuPolicy.getCoarseLevelRhs();
    auto& gpuCoarseRhs = gpuPolicy.getCoarseLevelRhs();

    // Check that restriction results match
    std::vector<double> gpuCoarseData = gpuCoarseRhs.asStdVector();
    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_CLOSE(cpuCoarseRhs[i][pressureVarIndex], gpuCoarseData[i], 1e-10);
    }

    // Set some values in the coarse lhs for prolongation test
    auto& cpuCoarseLhs = cpuPolicy.getCoarseLevelLhs();
    auto& gpuCoarseLhs = gpuPolicy.getCoarseLevelLhs();

    for (int i = 0; i < N; ++i) {
        cpuCoarseLhs[i][pressureVarIndex] = distribution(generator);
    }

    // Copy coarse lhs to GPU
    gpuCoarseLhs.copyFromHost(cpuCoarseLhs);

    // Test prolongation (moveToFineLevel)
    VectorType cpuFineResult(N);
    GpuVectorType gpuFineResult(N * blockSize);

    // Initialize fine result to zero, important when not using transpose
    // as we only update the pressure variable
    cpuFineResult = 0.0;
    gpuFineResult = 0.0;

    cpuPolicy.moveToFineLevel(cpuFineResult);
    gpuPolicy.moveToFineLevel(gpuFineResult);

    // Check that prolongation results match
    std::vector<double> cpuFineResultData(N * blockSize);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < blockSize; ++j) {
            cpuFineResultData[i * blockSize + j] = cpuFineResult[i][j];
        }
    }

    std::vector<double> gpuFineResultData = gpuFineResult.asStdVector();

    // Check that prolongation results match
    BOOST_CHECK_EQUAL(cpuFineResultData.size(), gpuFineResultData.size());
    for (size_t i = 0; i < cpuFineResultData.size(); ++i) {
        BOOST_CHECK_CLOSE(cpuFineResultData[i], gpuFineResultData[i], 1e-10);
    }
}

BOOST_DATA_TEST_CASE(TestPressureTransferPolicyTranspose,
                    boost::unit_test::data::make({0, 1, 2}), 
                    pressureVarIndex)
{
    runPressureTransferPolicyTest<true>(pressureVarIndex);
}

BOOST_DATA_TEST_CASE(TestPressureTransferPolicyStandard,
                    boost::unit_test::data::make({0, 1, 2}), 
                    pressureVarIndex)
{
    runPressureTransferPolicyTest<false>(pressureVarIndex);
}