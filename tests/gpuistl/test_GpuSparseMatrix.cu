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

#define BOOST_TEST_MODULE TestGpuSparseMatrix

#include <opm/simulators/linalg/gpuistl/GpuSparseMatrixWrapper.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpusparse_matrix_operations.hpp>
#include <opm/simulators/linalg/gpuistl/gpu_smart_pointer.hpp>
#include <opm/simulators/linalg/istlsparsematrixadapter.hh>

#include <opm/common/utility/pointerArithmetic.hpp>
#include <opm/simulators/linalg/gpuistl/detail/preconditionerKernels/ILU_variants_helper_kernels.hpp>

#include <dune/istl/bcrsmatrix.hh>

#include <boost/mpl/range_c.hpp>
#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include <cstddef>
#include <random>
#include <type_traits>

BOOST_AUTO_TEST_CASE(TestConstruction1D)
{
    // Here we will test a simple 1D finite difference scheme for
    // the Laplace equation:
    //
    //    -\Delta u = f on [0,1]
    //
    // Using a central difference approximation of \Delta u, this can
    // be approximated by
    //
    //    -(u_{i+1}-2u_i+u_{i-1})/Dx^2 = f(x_i)
    //
    // giving rise to the matrix
    //
    //     -2  1  0  0 ... 0  0
    //      1 -2  1  0  0 ... 0
    //      ....
    //      0  0  0  ...1 -2  1
    //      0  0  0  ...   1 -2

    const int N = 5;
    const int nonZeroes = N * 3 - 2;
    using M = Dune::FieldMatrix<double, 1, 1>;
    using SpMatrix = Dune::BCRSMatrix<M>;

    SpMatrix B(N, N, nonZeroes, SpMatrix::row_wise);
    for (auto row = B.createbegin(); row != B.createend(); ++row) {
        // Add nonzeros for left neighbour, diagonal and right neighbour
        if (row.index() > 0) {
            row.insert(row.index() - 1);
        }
        row.insert(row.index());
        if (row.index() < B.N() - 1) {
            row.insert(row.index() + 1);
        }
    }
    // This might not be the most elegant way of filling in a Dune sparse matrix, but it works.
    for (int i = 0; i < N; ++i) {
        B[i][i] = -2;
        if (i < N - 1) {
            B[i][i + 1] = 1;
        }

        if (i > 0) {
            B[i][i - 1] = 1;
        }
    }

    auto gpuSparseMatrix = Opm::gpuistl::GpuSparseMatrixWrapper<double>::fromMatrix(B);

    const auto& nonZeroValuesCuda = gpuSparseMatrix.getNonZeroValues();
    std::vector<double> buffer(gpuSparseMatrix.nonzeroes(), 0.0);
    nonZeroValuesCuda.copyToHost(buffer.data(), buffer.size());
    const double* nonZeroElements = static_cast<const double*>(&((B[0][0][0][0])));
    BOOST_CHECK_EQUAL_COLLECTIONS(buffer.begin(), buffer.end(), nonZeroElements, nonZeroElements + B.nonzeroes());
    BOOST_CHECK_EQUAL(N * 3 - 2, gpuSparseMatrix.nonzeroes());

    std::vector<int> rowIndicesFromCUDA(N + 1);
    gpuSparseMatrix.getRowIndices().copyToHost(rowIndicesFromCUDA.data(), rowIndicesFromCUDA.size());
    BOOST_CHECK_EQUAL(rowIndicesFromCUDA[0], 0);
    BOOST_CHECK_EQUAL(rowIndicesFromCUDA[1], 2);
    for (int i = 2; i < N; ++i) {
        BOOST_CHECK_EQUAL(rowIndicesFromCUDA[i], rowIndicesFromCUDA[i - 1] + 3);
    }


    std::vector<int> columnIndicesFromCUDA(B.nonzeroes(), 0);
    gpuSparseMatrix.getColumnIndices().copyToHost(columnIndicesFromCUDA.data(), columnIndicesFromCUDA.size());

    BOOST_CHECK_EQUAL(columnIndicesFromCUDA[0], 0);
    BOOST_CHECK_EQUAL(columnIndicesFromCUDA[1], 1);
    // TODO: Check rest
}

__global__ void checkDiagonalsKernel(double** diagPtrs, std::array<double, 3>* values)
{
    const int row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row == 0)
    {
        (*values)[0] = *diagPtrs[0];
        (*values)[1] = *diagPtrs[1];
        (*values)[2] = *diagPtrs[2];
    }
}

// This test ensures that getDiagPtrs works correctly for scalar matrices
BOOST_AUTO_TEST_CASE(TestGetDiagPtrs)
{
    const int N = 3;
    const int nonZeroes = 6;
    using M = Dune::FieldMatrix<double, 1, 1>;
    using SpMatrix = Dune::BCRSMatrix<M>;

    /*
    Fill in a matrix with this sparsity pattern:
        | * * * |
        |   *   |
        | *   * |
    */
    SpMatrix B(N, N, nonZeroes, SpMatrix::row_wise);
    for (auto row = B.createbegin(); row != B.createend(); ++row)
    {
        if (row.index() == 0)
        {
            row.insert(0);
            row.insert(1);
            row.insert(2);
        }
        else if (row.index() == 1)
        {
            row.insert(1);
        }
        else
        {
            row.insert(0);
            row.insert(2);
        }
    }

    // This might not be the most elegant way of filling in a Dune sparse matrix, but it works.
    B[0][0] = 1.0;
    B[0][1] = 2.0;
    B[0][2] = 3.0;
    B[1][1] = 4.0;
    B[2][0] = 5.0;
    B[2][2] = 6.0;

    auto gpuSparseMatrix = Opm::gpuistl::GpuSparseMatrixWrapper<double>::fromMatrix(B);

    auto cpuValues = gpuSparseMatrix.getNonZeroValues().asStdVector();
    BOOST_CHECK_EQUAL(cpuValues[0], 1.0);
    BOOST_CHECK_EQUAL(cpuValues[1], 2.0);
    BOOST_CHECK_EQUAL(cpuValues[2], 3.0);
    BOOST_CHECK_EQUAL(cpuValues[3], 4.0);
    BOOST_CHECK_EQUAL(cpuValues[4], 5.0);
    BOOST_CHECK_EQUAL(cpuValues[5], 6.0);

    auto diagPtrs = Opm::gpuistl::detail::getDiagPtrs(gpuSparseMatrix);
    auto diagPtrsView = Opm::gpuistl::GpuView<double*>(diagPtrs.data(), diagPtrs.size());

    auto d_result = Opm::gpuistl::make_gpu_unique_ptr<std::array<double, 3>>({-1});

    checkDiagonalsKernel<<<1, 1>>>(diagPtrsView.data(), d_result.get());

    std::array<double, 3> h_result = Opm::gpuistl::copyFromGPU(d_result);
    BOOST_CHECK_EQUAL(h_result[0], 1.0);
    BOOST_CHECK_EQUAL(h_result[1], 4.0);
    BOOST_CHECK_EQUAL(h_result[2], 6.0);
}

// Kernel to check diagonal block elements
__global__ void checkDiagBlocksKernel (double** diagPtrs, std::array<double, 12>* values) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx == 0) {
        // Extract diagonal block elements (2x2 blocks stored row-major)
        (*values)[0] = diagPtrs[0][0];  // B[0][0][0][0]
        (*values)[1] = diagPtrs[0][1];  // B[0][0][0][1]
        (*values)[2] = diagPtrs[0][2];  // B[0][0][1][0]
        (*values)[3] = diagPtrs[0][3];  // B[0][0][1][1]
        
        (*values)[4] = diagPtrs[1][0];  // B[1][1][0][0]
        (*values)[5] = diagPtrs[1][1];  // B[1][1][0][1]
        (*values)[6] = diagPtrs[1][2];  // B[1][1][1][0]
        (*values)[7] = diagPtrs[1][3];  // B[1][1][1][1]
        
        (*values)[8] = diagPtrs[2][0];  // B[2][2][0][0]
        (*values)[9] = diagPtrs[2][1];  // B[2][2][0][1]
        (*values)[10] = diagPtrs[2][2]; // B[2][2][1][0]
        (*values)[11] = diagPtrs[2][3]; // B[2][2][1][1]
    }
};

// Extra getDiagPtrs test to verify functionality on block-matrices
BOOST_AUTO_TEST_CASE(TestGetDiagPtrsBlockedNonScalar)
{
    const int N = 3;
    const int nonZeroes = 6;
    using M = Dune::FieldMatrix<double, 2, 2>;
    using SpMatrix = Dune::BCRSMatrix<M>;

    /*
    Fill in a blockmatrix with this sparsity pattern:
        | * * * |
        |   *   |
        | *   * |
    */
    SpMatrix B(N, N, nonZeroes, SpMatrix::row_wise);
    for (auto row = B.createbegin(); row != B.createend(); ++row)
    {
        if (row.index() == 0)
        {
            row.insert(0);
            row.insert(1);
            row.insert(2);
        }
        else if (row.index() == 1)
        {
            row.insert(1);
        }
        else
        {
            row.insert(0);
            row.insert(2);
        }
    }

    // Fill in the block matrix with 2x2 blocks
    B[0][0][0][0] = 1.0;
    B[0][0][0][1] = 1.1;
    B[0][0][1][0] = 1.2;
    B[0][0][1][1] = 1.3;
    
    B[0][1][0][0] = 2.0;
    B[0][1][0][1] = 2.1;
    B[0][1][1][0] = 2.2;
    B[0][1][1][1] = 2.3;
    
    B[0][2][0][0] = 3.0;
    B[0][2][0][1] = 3.1;
    B[0][2][1][0] = 3.2;
    B[0][2][1][1] = 3.3;
    
    B[1][1][0][0] = 4.0;
    B[1][1][0][1] = 4.1;
    B[1][1][1][0] = 4.2;
    B[1][1][1][1] = 4.3;
    
    B[2][0][0][0] = 5.0;
    B[2][0][0][1] = 5.1;
    B[2][0][1][0] = 5.2;
    B[2][0][1][1] = 5.3;
    
    B[2][2][0][0] = 6.0;
    B[2][2][0][1] = 6.1;
    B[2][2][1][0] = 6.2;
    B[2][2][1][1] = 6.3;

    auto gpuSparseMatrix = Opm::gpuistl::GpuSparseMatrixWrapper<double>::fromMatrix(B);

    auto diagPtrs = Opm::gpuistl::detail::getDiagPtrs(gpuSparseMatrix);
    
    // Check that we have the correct number of diagonal pointers
    BOOST_CHECK_EQUAL(diagPtrs.size(), N);

    // Create kernel to verify diagonal block pointers
    auto diagPtrsView = Opm::gpuistl::GpuView<double*>(diagPtrs.data(), diagPtrs.size());
    auto d_result = Opm::gpuistl::make_gpu_unique_ptr<std::array<double, 12>>({-1});

    checkDiagBlocksKernel<<<1, 1>>>(diagPtrsView.data(), d_result.get());

    std::array<double, 12> h_result = Opm::gpuistl::copyFromGPU(d_result);
    
    // Verify diagonal block elements
    BOOST_CHECK_EQUAL(h_result[0], 1.0);
    BOOST_CHECK_EQUAL(h_result[1], 1.1);
    BOOST_CHECK_EQUAL(h_result[2], 1.2);
    BOOST_CHECK_EQUAL(h_result[3], 1.3);
    
    BOOST_CHECK_EQUAL(h_result[4], 4.0);
    BOOST_CHECK_EQUAL(h_result[5], 4.1);
    BOOST_CHECK_EQUAL(h_result[6], 4.2);
    BOOST_CHECK_EQUAL(h_result[7], 4.3);
    
    BOOST_CHECK_EQUAL(h_result[8], 6.0);
    BOOST_CHECK_EQUAL(h_result[9], 6.1);
    BOOST_CHECK_EQUAL(h_result[10], 6.2);
    BOOST_CHECK_EQUAL(h_result[11], 6.3);
}

// Template function to run the random sparsity matrix test with a given block size
template<class GpuMatrix, std::size_t dim>
void runRandomSparsityMatrixTest()
{
    std::srand(0);
    double nonzeroPercent = 0.2;
    std::mt19937 generator;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    const std::size_t N = 300;
    using M = Dune::FieldMatrix<double, dim, dim>;
    using SpMatrix = Dune::BCRSMatrix<M>;
    using Vector = Dune::BlockVector<Dune::FieldVector<double, dim>>;

    std::vector<std::vector<std::size_t>> nonzerocols(N);
    int nonZeroes = 0;
    for (std::size_t row = 0; row < N; ++row) {
        // Always include the diagonal element to ensure each row has at least one entry
        nonzerocols.at(row).push_back(row);
        nonZeroes++;

        // Add other elements based on random sparsity
        for (std::size_t col = 0; col < N; ++col) {
            if (col != row && distribution(generator) < nonzeroPercent) {
                nonzerocols.at(row).push_back(col);
                nonZeroes++;
            }
        }
    }
    SpMatrix B(N, N, nonZeroes, SpMatrix::row_wise);
    for (auto row = B.createbegin(); row != B.createend(); ++row) {
        for (std::size_t j = 0; j < nonzerocols[row.index()].size(); ++j) {
            row.insert(nonzerocols[row.index()][j]);
        }
    }
    // This might not be the most elegant way of filling in a Dune sparse matrix, but it works.
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < nonzerocols[i].size(); ++j) {
            for (std::size_t c1 = 0; c1 < dim; ++c1) {
                for (std::size_t c2 = 0; c2 < dim; ++c2) {
                    B[i][nonzerocols[i][j]][c1][c2] = distribution(generator);
                }
            }
        }
    }

/*
    HIP currently does not properly support the generic sparse API
    CUDA 13 fully supports it, whereas lower CUDA versions only support it for block size 1
*/
#if USE_HIP || (!USE_HIP && CUDA_VERSION < 13000)
    constexpr bool isGeneric     = std::is_same_v<GpuMatrix, Opm::gpuistl::GpuSparseMatrixGeneric<double>>;
    constexpr bool expectFailure = isGeneric && (dim != 1);
#else
    constexpr bool expectFailure = false;
#endif

    if constexpr (expectFailure) {
        bool threw = false;
        try {
            auto shouldFail = GpuMatrix::fromMatrix(B);
            (void)shouldFail;
        } catch (...) {
            threw = true;
        }
        BOOST_CHECK(threw);
        return;
    }

    auto gpuSparseMatrix = GpuMatrix::fromMatrix(B);
    // check each column
    for (std::size_t component = 0; component < N * dim; component += N) {
        std::vector<double> inputDataX(N * dim, 0.0);
        inputDataX[component] = 1.0;
        std::vector<double> inputDataY(N * dim, .25);
        auto inputVectorX = Opm::gpuistl::GpuVector<double>(inputDataX.data(), inputDataX.size());
        auto inputVectorY = Opm::gpuistl::GpuVector<double>(inputDataY.data(), inputDataY.size());
        Vector xHost(N), yHost(N);
        yHost = inputDataY[0];
        inputVectorX.copyToHost(xHost);
        const double alpha = 1.42;
        gpuSparseMatrix.usmv(alpha, inputVectorX, inputVectorY);

        inputVectorY.copyToHost(inputDataY);

        B.usmv(alpha, xHost, yHost);
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t c = 0; c < dim; ++c) {
                BOOST_CHECK_CLOSE(inputDataY[i * dim + c], yHost[i][c], 1e-7);
            }
        }
        inputVectorX.copyToHost(xHost);

        gpuSparseMatrix.mv(inputVectorX, inputVectorY);

        inputVectorY.copyToHost(inputDataY);

        B.mv(xHost, yHost);
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t c = 0; c < dim; ++c) {
                BOOST_CHECK_CLOSE(inputDataY[i * dim + c], yHost[i][c], 1e-7);
            }
        }
    }
}

// Define test cases for block sizes 1-6 using boost::mpl::range_c
typedef boost::mpl::range_c<int, 1, 7> block_sizes;

BOOST_AUTO_TEST_CASE_TEMPLATE(RandomSparsityMatrix, T, block_sizes)
{
    runRandomSparsityMatrixTest<Opm::gpuistl::GpuSparseMatrixWrapper<double>, T::value>();
}

typedef boost::mpl::range_c<int, 1, 7> block_sizes;

BOOST_AUTO_TEST_CASE_TEMPLATE(RandomSparsityMatrixGeneric, T, block_sizes)
{
    runRandomSparsityMatrixTest<Opm::gpuistl::GpuSparseMatrixGeneric<double>, T::value>();
}

__global__ void checkPointers(double* data, std::array<double*, 6>* ptrs, std::array<bool, 6>* correctPtrs)
{
    const int row = blockIdx.x * blockDim.x + threadIdx.x;

    if (row == 0)
    {
        const std::array<double, 6> expectedValues = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

        for (int i = 0; i < 6; ++i) {
            double actualValue = *((*ptrs)[i]); // Dereference pointer to array, index the array, dereference the double* in the array
            (*correctPtrs)[i] = (actualValue == expectedValues[i]);
        }
    }
}

BOOST_AUTO_TEST_CASE(TestFindingPointersToElements)
{
    /*
        In this test I want to check that constructing a Dune Matrix with the istl wrapper will
        produce a matrix stored with in the same way as the GpuSparseMatrix

        That should be the case with the way we produce the jacobian matrix in the tpfalinearizer (using reserve and then filling in values)
    */
    using CPUBlockType = Dune::FieldMatrix<double, 1, 1>;

    Opm::Linear::IstlSparseMatrixAdapter<CPUBlockType> cpuMatrix(3,3);

    std::vector<std::set<int>> sparsityPattern =
        {
            {0,1,2},
            {1},
            {0,2}
        };
    cpuMatrix.reserve(sparsityPattern);

    // This way of filling in the matrix is chosen to mirror what happens in the tpfalinaerizer
    *(cpuMatrix.blockAddress(0, 0)) = CPUBlockType(1.0);
    *(cpuMatrix.blockAddress(0, 1)) = CPUBlockType(2.0);
    *(cpuMatrix.blockAddress(0, 2)) = CPUBlockType(3.0);
    *(cpuMatrix.blockAddress(1, 1)) = CPUBlockType(4.0);
    *(cpuMatrix.blockAddress(2, 0)) = CPUBlockType(5.0);
    *(cpuMatrix.blockAddress(2, 2)) = CPUBlockType(6.0);

    auto gpuMatrix = Opm::gpuistl::GpuSparseMatrixWrapper<double>::fromMatrix(cpuMatrix.istlMatrix());

    // Is getting access to the pure underlying data supposed to be this hard?
    // block
    std::array<double*, 6> ptrsCPU = {
        &(cpuMatrix.istlMatrix()[0][0][0][0]),
        &(cpuMatrix.istlMatrix()[0][1][0][0]),
        &(cpuMatrix.istlMatrix()[0][2][0][0]),
        &(cpuMatrix.istlMatrix()[1][1][0][0]),
        &(cpuMatrix.istlMatrix()[2][0][0][0]),
        &(cpuMatrix.istlMatrix()[2][2][0][0])
    };

    std::array<double*, 6> ptrsGpuComputed;
    double* gpuDataBuffer = gpuMatrix.getNonZeroValues().data();
    double* cpuDataBuffer = &((*cpuMatrix.blockAddress(0, 0))[0][0]);

    for (size_t i = 0; i < 6; ++i) {
        ptrsGpuComputed[i] = Opm::ComputePtrBasedOnOffsetInOtherBuffer(gpuDataBuffer, gpuMatrix.nonzeroes(), cpuDataBuffer, cpuMatrix.istlMatrix().nonzeroes(), ptrsCPU[i]);
    }

    auto ptrsGpu = Opm::gpuistl::make_gpu_unique_ptr<std::array<double*, 6>>(ptrsGpuComputed);
    auto correctValues = Opm::gpuistl::make_gpu_unique_ptr<std::array<bool, 6>>({false});

    checkPointers<<<1,1>>>(static_cast<double*>(gpuDataBuffer), ptrsGpu.get(), correctValues.get());

    auto correctValuesCPU = Opm::gpuistl::copyFromGPU(correctValues);

    BOOST_CHECK(correctValuesCPU[0]);
    BOOST_CHECK(correctValuesCPU[1]);
    BOOST_CHECK(correctValuesCPU[2]);
    BOOST_CHECK(correctValuesCPU[3]);
    BOOST_CHECK(correctValuesCPU[4]);
    BOOST_CHECK(correctValuesCPU[5]);
}
