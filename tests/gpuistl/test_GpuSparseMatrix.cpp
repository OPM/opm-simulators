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

#include <opm/simulators/linalg/gpuistl/GpuSparseMatrix.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>

#include <dune/istl/bcrsmatrix.hh>

#include <boost/mpl/range_c.hpp>
#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include <cstddef>
#include <random>

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

    auto gpuSparseMatrix = Opm::gpuistl::GpuSparseMatrix<double>::fromMatrix(B);

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

// Template function to run the random sparsity matrix test with a given block size
template<std::size_t dim>
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

    auto gpuSparseMatrix = Opm::gpuistl::GpuSparseMatrix<double>::fromMatrix(B);
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
    runRandomSparsityMatrixTest<T::value>();
}
