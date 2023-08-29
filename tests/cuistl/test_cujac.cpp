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
#include <opm/simulators/linalg/cuistl/detail/cusparse_matrix_operations.hpp>
#include <opm/simulators/linalg/cuistl/detail/vector_operations.hpp>
#include <opm/simulators/linalg/cuistl/detail/fix_zero_diagonal.hpp>
#include <opm/simulators/linalg/cuistl/PreconditionerAdapter.hpp>

#include <iostream>

using NumericTypes = boost::mpl::list<double, float>;

BOOST_AUTO_TEST_CASE_TEMPLATE(CUJACApplyBlocksize2, T, NumericTypes)
{
     /*
        Test data to validate jacobi preconditioner, expected result is x_1
            | |3 1|  | 1  0|       | |1| |     | |2| |       | |   1| |
            | |2 1|  | 0  1|       | |2| |     | |1| |       | |   0| |
        A = |              | x_0 = |     | b = |     | x_1 = |        |
            | |0 0|  |-1  0|       | |1| |     | |3| |       | |  -1| |
            | |0 0|  | 0 -1|       | |1| |     | |4| |       | |-1.5| |
    */
    const int N = 2;
    const int blocksize = 2;
    const int nonZeroes = 3;
    using M = Dune::FieldMatrix<T, blocksize, blocksize>;
    using SpMatrix = Dune::BCRSMatrix<M>;
    using Vector = Dune::BlockVector<Dune::FieldVector<T, blocksize>>;
    using cujac = Opm::cuistl::CuJac<SpMatrix, Opm::cuistl::CuVector<T>, Opm::cuistl::CuVector<T>>;
    
    SpMatrix B(N, N, nonZeroes, SpMatrix::row_wise);
    for (auto row = B.createbegin(); row != B.createend(); ++row) {
        row.insert(row.index());
        if (row.index() == 0) {
            row.insert(row.index() + 1);
        }
    }

    B[0][0][0][0]=3.0;
    B[0][0][0][1]=1.0;
    B[0][0][1][0]=2.0;
    B[0][0][1][1]=1.0;

    B[0][1][0][0]=1.0;
    B[0][1][1][1]=1.0;

    B[1][1][0][0]=-1.0;
    B[1][1][1][1]=-1.0;

    auto CuJac = Opm::cuistl::PreconditionerAdapter<Vector, Vector, cujac>(std::make_shared<cujac>(B, 0.5));
    
    Vector h_dune_x(2);
    h_dune_x[0][0] = 1.0;
    h_dune_x[0][1] = 2.0;
    h_dune_x[1][0] = 1.0;
    h_dune_x[1][1] = 1.0;

    Vector h_dune_b(2);
    h_dune_b[0][0] = 2.0;
    h_dune_b[0][1] = 1.0;
    h_dune_b[1][0] = 3.0;
    h_dune_b[1][1] = 4.0;
    
    const T expected_ans[2][2] = {{1.0, 0.0},{-1.0,-3.0/2.0}};

    CuJac.apply(h_dune_x, h_dune_b);
    BOOST_CHECK_CLOSE(h_dune_x[0][0], expected_ans[0][0], 1e-7);
    BOOST_CHECK_CLOSE(h_dune_x[0][1], expected_ans[0][1], 1e-7);
    BOOST_CHECK_CLOSE(h_dune_x[1][0], expected_ans[1][0], 1e-7);
    BOOST_CHECK_CLOSE(h_dune_x[1][1], expected_ans[1][1], 1e-7);
    BOOST_CHECK(true);
}