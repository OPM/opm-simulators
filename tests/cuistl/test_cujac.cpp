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

#define BOOST_TEST_MODULE TestCuJac
#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>
#include <cuda_runtime.h>
#include <dune/istl/bcrsmatrix.hh>
#include <opm/simulators/linalg/cuistl/CuJac.hpp>
#include <opm/simulators/linalg/cuistl/CuSparseMatrix.hpp>
#include <opm/simulators/linalg/cuistl/CuVector.hpp>
#include <opm/simulators/linalg/cuistl/PreconditionerAdapter.hpp>
#include <opm/simulators/linalg/cuistl/detail/cusparse_matrix_operations.hpp>
#include <opm/simulators/linalg/cuistl/detail/fix_zero_diagonal.hpp>
#include <opm/simulators/linalg/cuistl/detail/vector_operations.hpp>

using NumericTypes = boost::mpl::list<double, float>;

BOOST_AUTO_TEST_CASE_TEMPLATE(CUJACApplyBlocksize2, T, NumericTypes)
{
    /*
       Test data to validate jacobi preconditioner, expected result is x_1, and relaxation factor is 0.5
           | |3 1|  | 1  0| |     | |2| |       | | 0.5| |
           | |2 1|  | 0  1| |     | |1| |       | |-0.5| |
       A = |                | d = |     |   v = |        |
           | |0 0|  |-1  0| |     | |3| |       | |-1.5| |
           | |0 0|  | 0 -1| |     | |4| |       | |-2.0| |
   */
    const int N = 2;
    constexpr int blocksize = 2;
    const int nonZeroes = 3;
    using M = Dune::FieldMatrix<T, blocksize, blocksize>;
    using SpMatrix = Dune::BCRSMatrix<M>;
    using Vector = Dune::BlockVector<Dune::FieldVector<T, blocksize>>;
    using CuJac = Opm::cuistl::CuJac<SpMatrix, Opm::cuistl::CuVector<T>, Opm::cuistl::CuVector<T>>;

    SpMatrix B(N, N, nonZeroes, SpMatrix::row_wise);
    for (auto row = B.createbegin(); row != B.createend(); ++row) {
        row.insert(row.index());
        if (row.index() == 0) {
            row.insert(row.index() + 1);
        }
    }

    B[0][0][0][0] = 3.0;
    B[0][0][0][1] = 1.0;
    B[0][0][1][0] = 2.0;
    B[0][0][1][1] = 1.0;

    B[0][1][0][0] = 1.0;
    B[0][1][1][1] = 1.0;

    B[1][1][0][0] = -1.0;
    B[1][1][1][1] = -1.0;

    auto cujac = Opm::cuistl::PreconditionerAdapter<Vector, Vector, CuJac>(std::make_shared<CuJac>(B, 0.5));

    Vector vVector(2);
    Vector dVector(2);
    dVector[0][0] = 2.0;
    dVector[0][1] = 1.0;
    dVector[1][0] = 3.0;
    dVector[1][1] = 4.0;

    const T expectedAns[2][2] = {{1.0 / 2.0, -1.0 / 2.0}, {-3.0 / 2.0, -2.0}};

    cujac.apply(vVector, dVector);
    BOOST_CHECK_CLOSE(vVector[0][0], expectedAns[0][0], 1e-7);
    BOOST_CHECK_CLOSE(vVector[0][1], expectedAns[0][1], 1e-7);
    BOOST_CHECK_CLOSE(vVector[1][0], expectedAns[1][0], 1e-7);
    BOOST_CHECK_CLOSE(vVector[1][1], expectedAns[1][1], 1e-7);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(CUJACApplyBlocksize1, T, NumericTypes)
{
    /*
       Test data to validate jacobi preconditioner, expected result is x_1, relaxation factor is 0.5
           | 3  1  1  0|       |1|       | 1/3|
           | 2  1  0  1|       |2|       | 1/2|
       A = | 0  0 -1  0|   d = |1|   v = |-3/2|
           | 0  0  0 -1|       |1|       |  -2|
   */
    const int N = 4;
    constexpr int blocksize = 1;
    const int nonZeroes = 8;
    using M = Dune::FieldMatrix<T, blocksize, blocksize>;
    using SpMatrix = Dune::BCRSMatrix<M>;
    using Vector = Dune::BlockVector<Dune::FieldVector<T, blocksize>>;
    using CuJac = Opm::cuistl::CuJac<SpMatrix, Opm::cuistl::CuVector<T>, Opm::cuistl::CuVector<T>>;

    SpMatrix B(N, N, nonZeroes, SpMatrix::row_wise);
    for (auto row = B.createbegin(); row != B.createend(); ++row) {
        row.insert(row.index());
        if (row.index() == 0) {
            row.insert(row.index() + 1);
            row.insert(row.index() + 2);
        }
        if (row.index() == 1) {
            row.insert(row.index() - 1);
            row.insert(row.index() + 2);
        }
    }

    B[0][0][0][0] = 3.0;
    B[0][1][0][0] = 1.0;
    B[0][2][0][0] = 1.0;

    B[1][0][0][0] = 2.0;
    B[1][1][0][0] = 1.0;
    B[1][3][0][0] = 1.0;

    B[2][2][0][0] = -1.0;
    B[3][3][0][0] = -1.0;

    auto cujac = Opm::cuistl::PreconditionerAdapter<Vector, Vector, CuJac>(std::make_shared<CuJac>(B, 0.5));

    Vector vVector(4);
    Vector dVector(4);
    dVector[0] = 2.0;
    dVector[1] = 1.0;
    dVector[2] = 3.0;
    dVector[3] = 4.0;

    const T expectedAns[4] = {1.0 / 3.0, 1.0 / 2.0, -3.0 / 2.0, -2.0};

    cujac.apply(vVector, dVector);
    BOOST_CHECK_CLOSE(vVector[0], expectedAns[0], 1e-7);
    BOOST_CHECK_CLOSE(vVector[1], expectedAns[1], 1e-7);
    BOOST_CHECK_CLOSE(vVector[2], expectedAns[2], 1e-7);
    BOOST_CHECK_CLOSE(vVector[3], expectedAns[3], 1e-7);
}
