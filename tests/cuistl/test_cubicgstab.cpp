/*
  Copyright SINTEF AS 2022

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

#define BOOST_TEST_MODULE TestCuBICGSTAB
#define BOOST_TEST_NO_MAIN


#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixmarket.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <opm/simulators/linalg/ParallelOverlappingILU0.hpp>
#include <opm/simulators/linalg/cuistl/CuSeqILU0.hpp>
#include <opm/simulators/linalg/cuistl/CuVector.hpp>
#include <opm/simulators/linalg/cuistl/PreconditionerAdapter.hpp>
#include <opm/simulators/linalg/cuistl/impl/cuda_safe_call.hpp>

#include <limits>
#include <memory>
#include <random>

using NumericTypes = boost::mpl::list<double>;

BOOST_AUTO_TEST_CASE_TEMPLATE(TestFiniteDifference1D, T, NumericTypes)
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
    using M = Dune::FieldMatrix<T, 2, 2>;
    using SpMatrix = Dune::BCRSMatrix<M>;
    using Vector = Dune::BlockVector<Dune::FieldVector<T, 1>>;
    using CuILU0 = Opm::cuistl::CuSeqILU0<SpMatrix, Opm::cuistl::CuVector<T>, Opm::cuistl::CuVector<T>>;

    SpMatrix B(N, N, nonZeroes, SpMatrix::row_wise);
    for (auto row = B.createbegin(); row != B.createend(); ++row) {
        row.insert(row.index());
        if (row.index() < N - 1) {
            row.insert(row.index() + 1);
        }
        if (row.index() > 0) {
            row.insert(row.index() - 1);
        }
    }
    // This might not be the most elegant way of filling in a Dune sparse matrix, but it works.
    for (int i = 0; i < N; ++i) {
        B[i][i][0][0] = -2;
        B[i][i][1][1] = -2;
        B[i][i][0][1] = 1;
        B[i][i][1][0] = 1;
    }

    auto BonGPU = std::make_shared<Opm::cuistl::CuSparseMatrix<T>>(Opm::cuistl::CuSparseMatrix<T>::fromMatrix(B));
    auto BOperator = std::make_shared<
        Dune::MatrixAdapter<Opm::cuistl::CuSparseMatrix<T>, Opm::cuistl::CuVector<T>, Opm::cuistl::CuVector<T>>>(
        BonGPU);
    auto cuILU = std::make_shared<CuILU0>(B, 1.0);
    auto scalarProduct = std::make_shared<Dune::SeqScalarProduct<Opm::cuistl::CuVector<T>>>();

    auto solver = Dune::BiCGSTABSolver<Opm::cuistl::CuVector<T>>(BOperator, scalarProduct, cuILU, .001, 100, 0);
    std::vector<T> correct(N * 2, 2.0);
    std::vector<T> initialGuess(N * 2, 0.0);
    Opm::cuistl::CuVector<T> x(N * 2);
    Opm::cuistl::CuVector<T> y(N * 2);
    x.copyFromHost(correct.data(), correct.size());
    BonGPU->mv(x, y);
    x.copyFromHost(initialGuess.data(), initialGuess.size());

    Dune::InverseOperatorResult result;
    Opm::cuistl::CuVector<T> tmp(N * 2);
    tmp.copyFromHost(correct.data(), correct.size());
    tmp -= x;
    auto normBefore = tmp.two_norm();
    BOOST_CHECK_GT(normBefore, 0.5);
    solver.apply(x, y, result);
    tmp.copyFromHost(correct.data(), correct.size());
    tmp -= x;
    auto normAfter = tmp.two_norm();
    BOOST_CHECK_CLOSE(normAfter, 0.0, 1e-7);
}

BOOST_AUTO_TEST_CASE(TestLoadFromFile, *boost::unit_test::tolerance(1e-7))
{
    using T = double;
    static constexpr size_t dim = 3;
    using M = Dune::FieldMatrix<T, dim, dim>;
    using SpMatrix = Dune::BCRSMatrix<M>;
    using Vector = Dune::BlockVector<Dune::FieldVector<T, dim>>;
    using CuILU0 = Opm::cuistl::CuSeqILU0<SpMatrix, Opm::cuistl::CuVector<T>, Opm::cuistl::CuVector<T>>;
    using GPUVector = Opm::cuistl::CuVector<T>;
    std::srand(0);
    std::mt19937 generator;
    std::uniform_real_distribution<T> distribution(-100.0, 100.0);
    SpMatrix B;

    // TODO: Find a way to include this
    try {
        Dune::loadMatrixMarket(B, "../matrix.mm");
    } catch (Dune::IOError&) {
        // Try to load a smaller one
        Dune::loadMatrixMarket(B, "../tests/matr33.txt");
    }

    // B.compress();
    const size_t N = B.N();
    auto BonGPU = std::make_shared<Opm::cuistl::CuSparseMatrix<T>>(Opm::cuistl::CuSparseMatrix<T>::fromMatrix(B));
    auto BOperator = std::make_shared<
        Dune::MatrixAdapter<Opm::cuistl::CuSparseMatrix<T>, Opm::cuistl::CuVector<T>, Opm::cuistl::CuVector<T>>>(
        BonGPU);
    auto cuILU = std::make_shared<CuILU0>(B, 1.0);
    auto scalarProduct = std::make_shared<Dune::SeqScalarProduct<Opm::cuistl::CuVector<T>>>();
    const auto tolerance = 1e-7;
    auto solver = Dune::BiCGSTABSolver<Opm::cuistl::CuVector<T>>(BOperator, scalarProduct, cuILU, tolerance, 10000, 0);


    auto BCPUOperator = std::make_shared<Dune::MatrixAdapter<SpMatrix, Vector, Vector>>(B);
    auto scalarProductCPU = std::make_shared<Dune::SeqScalarProduct<Vector>>();
    auto ILUCPU
        = std::make_shared<Opm::ParallelOverlappingILU0<SpMatrix, Vector, Vector, Dune::Amg::SequentialInformation>>(
            B, 0, 1.0, Opm::MILU_VARIANT::ILU);
    auto solverCPU = Dune::BiCGSTABSolver<Vector>(BCPUOperator, scalarProductCPU, ILUCPU, tolerance, 10000, 0);

    auto checkCPUGPUEqual = [](const Vector& cpuVector, const GPUVector& gpuVector) {
        std::vector<T> gpuVectorOnHost(gpuVector.dim());
        gpuVector.copyToHost(gpuVectorOnHost);

        for (size_t i = 0; i < cpuVector.N(); ++i) {
            for (size_t j = 0; j < dim; ++j) {
                const auto normalizer = std::max(std::abs(cpuVector[i][j]), 0.1);
                BOOST_CHECK_LT(std::abs(cpuVector[i][j] - gpuVectorOnHost[i * dim + j]) / normalizer, 5e-3);
            }
        }
    };
    std::vector<T> correct(N * dim, 1.0);
    std::vector<T> initialGuess(N * dim);
    for (size_t i = 0; i < N; ++i) {
        initialGuess[i] = distribution(generator);
    }

    Opm::cuistl::CuVector<T> x(correct);
    Opm::cuistl::CuVector<T> y(initialGuess);

    Vector xHost(N), yHost(N);
    x.copyToHost(xHost);
    BonGPU->mv(x, y);
    B.mv(xHost, yHost);


    checkCPUGPUEqual(yHost, y);

    const size_t numberOfRandomVectorsToTry = 1000;


    for (size_t i = 0; i < numberOfRandomVectorsToTry; ++i) {
        std::vector<T> randomVector(N * dim);
        for (size_t j = 0; j < N * dim; ++j) {
            randomVector[j] = distribution(generator);
        }
        Vector xRandomCPU(N);
        Vector yRandomCPU(N);
        GPUVector xRandomGPU(randomVector);
        GPUVector yRandomGPU(N * dim);
        yRandomGPU = 0.0;
        yRandomCPU = 0.0;

        xRandomGPU.copyToHost(xRandomCPU);
        const double randomAlpha = distribution(generator);
        B.usmv(randomAlpha, xRandomCPU, yRandomCPU);
        BonGPU->usmv(randomAlpha, xRandomGPU, yRandomGPU);
        checkCPUGPUEqual(yRandomCPU, yRandomGPU);

        yRandomCPU = 0.0;
        yRandomGPU = 0.0;

        B.mv(xRandomCPU, yRandomCPU);
        BonGPU->mv(xRandomGPU, yRandomGPU);
        checkCPUGPUEqual(yRandomCPU, yRandomGPU);

        cuILU->apply(xRandomGPU, yRandomGPU);
        ILUCPU->apply(xRandomCPU, yRandomCPU);
        checkCPUGPUEqual(xRandomCPU, xRandomGPU);


        std::vector<T> randomVector2(N * dim);
        for (size_t j = 0; j < N * dim; ++j) {
            randomVector[j] = distribution(generator);
        }

        GPUVector xRandomGPU2(randomVector);
        Vector xRandomCPU2(N);

        Dune::InverseOperatorResult resultRandomGPU, resultRandomCPU;
        solver.apply(xRandomGPU2, yRandomGPU, resultRandomGPU);
        solverCPU.apply(xRandomCPU2, yRandomCPU, resultRandomCPU);

        checkCPUGPUEqual(xRandomCPU2, xRandomGPU2);
        BOOST_CHECK_LT(std::abs(resultRandomCPU.iterations - resultRandomGPU.iterations), 10);
        BOOST_CHECK_EQUAL(resultRandomCPU.converged, resultRandomGPU.converged);
        BOOST_CHECK_CLOSE(resultRandomCPU.conv_rate, resultRandomGPU.conv_rate, 50);
    }

    x.copyFromHost(initialGuess);
    x.copyToHost(xHost);

    Dune::InverseOperatorResult result, resultCPU;
    Opm::cuistl::CuVector<T> tmp(N * dim);
    tmp.copyFromHost(correct.data(), correct.size());
    tmp -= x;
    auto normBefore = tmp.two_norm();
    std::cout << normBefore << std::endl;
    BOOST_CHECK_GT(normBefore, 0.2);

    solver.apply(x, y, result);
    solverCPU.apply(xHost, yHost, resultCPU);
    tmp.copyFromHost(correct.data(), correct.size());
    tmp -= x;
    auto normAfter = tmp.two_norm();
    BOOST_CHECK_EQUAL(resultCPU.converged, result.converged);
    std::cout << resultCPU.iterations << std::endl;
    std::cout << result.iterations << std::endl;
    BOOST_CHECK_EQUAL(resultCPU.iterations, result.iterations);
    BOOST_CHECK_CLOSE(resultCPU.conv_rate, result.conv_rate, 50);
    BOOST_CHECK_LT(normAfter, tolerance);
}


bool
init_unit_test_func()
{
    return true;
}

int
main(int argc, char** argv)
{
    [[maybe_unused]] const auto& helper = Dune::MPIHelper::instance(argc, argv);
    boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
