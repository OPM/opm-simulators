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

#define BOOST_TEST_MODULE TestCuSeqILU0
#define BOOST_TEST_NO_MAIN


#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/preconditioners.hh>
#include <opm/simulators/linalg/cuistl/CuSeqILU0.hpp>
#include <opm/simulators/linalg/cuistl/CuVector.hpp>
#include <opm/simulators/linalg/cuistl/PreconditionerAdapter.hpp>
#include <opm/simulators/linalg/cuistl/detail/cuda_safe_call.hpp>

#include <limits>
#include <memory>


using NumericTypes = boost::mpl::list<double, float>;

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
    using M = Dune::FieldMatrix<T, 1, 1>;
    using SpMatrix = Dune::BCRSMatrix<M>;
    using Vector = Dune::BlockVector<Dune::FieldVector<T, 1>>;
    using CuILU0 = Opm::cuistl::CuSeqILU0<SpMatrix, Opm::cuistl::CuVector<T>, Opm::cuistl::CuVector<T>>;

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


    auto duneILU = Dune::SeqILU<SpMatrix, Vector, Vector>(B, 1.0);

    auto cuILU = Opm::cuistl::PreconditionerAdapter<Vector, Vector, CuILU0>(std::make_shared<CuILU0>(B, 1.0));

    // check for the standard basis {e_i}
    // (e_i=(0,...,0, 1 (i-th place), 0, ..., 0))
    for (int i = 0; i < N; ++i) {
        Vector inputVector(N);
        inputVector[i][0] = 1.0;
        Vector outputVectorDune(N);
        Vector outputVectorCuistl(N);
        duneILU.apply(outputVectorDune, inputVector);
        cuILU.apply(outputVectorCuistl, inputVector);

        for (int component = 0; component < N; ++component) {
            BOOST_CHECK_CLOSE(outputVectorDune[component][0],
                              outputVectorCuistl[component][0],
                              std::numeric_limits<T>::epsilon() * 1000);
        }
    }

    // Now we check that we can update the matrix. We basically just negate B
    B *= -1.0;
    auto duneILUNew = Dune::SeqILU<SpMatrix, Vector, Vector>(B, 1.0);
    cuILU.update();
    // check for the standard basis {e_i}
    // (e_i=(0,...,0, 1 (i-th place), 0, ..., 0))
    for (int i = 0; i < N; ++i) {
        Vector inputVector(N);
        inputVector[i][0] = 1.0;
        Vector outputVectorDune(N);
        Vector outputVectorCuistl(N);
        duneILUNew.apply(outputVectorDune, inputVector);
        cuILU.apply(outputVectorCuistl, inputVector);

        for (int component = 0; component < N; ++component) {
            BOOST_CHECK_CLOSE(outputVectorDune[component][0],
                              outputVectorCuistl[component][0],
                              std::numeric_limits<T>::epsilon() * 1000);
        }
    }
}


BOOST_AUTO_TEST_CASE_TEMPLATE(TestFiniteDifferenceBlock2, T, NumericTypes)
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
    using Vector = Dune::BlockVector<Dune::FieldVector<T, 2>>;
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


    auto duneILU = Dune::SeqILU<SpMatrix, Vector, Vector>(B, 1.0);

    auto cuILU = Opm::cuistl::PreconditionerAdapter<Vector, Vector, CuILU0>(std::make_shared<CuILU0>(B, 1.0));

    // check for the standard basis {e_i}
    // (e_i=(0,...,0, 1 (i-th place), 0, ..., 0))
    for (int i = 0; i < N; ++i) {
        Vector inputVector(N);
        inputVector[i][0] = 1.0;
        Vector outputVectorDune(N);
        Vector outputVectorCuistl(N);
        duneILU.apply(outputVectorDune, inputVector);
        cuILU.apply(outputVectorCuistl, inputVector);

        for (int component = 0; component < N; ++component) {
            BOOST_CHECK_CLOSE(outputVectorDune[component][0],
                              outputVectorCuistl[component][0],
                              std::numeric_limits<T>::epsilon() * 1000);
        }
    }

    // Now we check that we can update the matrix. We basically just negate B
    B *= -1.0;
    auto duneILUNew = Dune::SeqILU<SpMatrix, Vector, Vector>(B, 1.0);
    cuILU.update();
    // check for the standard basis {e_i}
    // (e_i=(0,...,0, 1 (i-th place), 0, ..., 0))
    for (int i = 0; i < N; ++i) {
        Vector inputVector(N);
        inputVector[i][0] = 1.0;
        Vector outputVectorDune(N);
        Vector outputVectorCuistl(N);
        duneILUNew.apply(outputVectorDune, inputVector);
        cuILU.apply(outputVectorCuistl, inputVector);

        for (int component = 0; component < N; ++component) {
            BOOST_CHECK_CLOSE(outputVectorDune[component][0],
                              outputVectorCuistl[component][0],
                              std::numeric_limits<T>::epsilon() * 1000);
        }
    }
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
