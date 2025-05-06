/*
  Copyright 2024 SINTEF AS

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

#define BOOST_TEST_MODULE TestGpuDILU

#include <boost/test/unit_test.hpp>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <memory>
#include <opm/simulators/linalg/DILU.hpp>
#include <opm/simulators/linalg/gpuistl/GpuDILU.hpp>
#include <opm/simulators/linalg/gpuistl/GpuSparseMatrix.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpusparse_matrix_operations.hpp>
#include <opm/simulators/linalg/gpuistl/PreconditionerCPUMatrixToGPUMatrix.hpp>
#include <random>
#include <vector>


using T = double;
using FM1x1 = Dune::FieldMatrix<T, 1, 1>;
using FM2x2 = Dune::FieldMatrix<T, 2, 2>;
using B1x1Vec = Dune::BlockVector<Dune::FieldVector<double, 1>>;
using B2x2Vec = Dune::BlockVector<Dune::FieldVector<double, 2>>;
using Sp1x1BlockMatrix = Dune::BCRSMatrix<FM1x1>;
using Sp2x2BlockMatrix = Dune::BCRSMatrix<FM2x2>;
using CuMatrix = Opm::gpuistl::GpuSparseMatrix<T>;
using CuIntVec = Opm::gpuistl::GpuVector<int>;
using CuFloatingPointVec = Opm::gpuistl::GpuVector<T>;
using GpuDilu1x1 = Opm::gpuistl::GpuDILU<Sp1x1BlockMatrix, CuFloatingPointVec, CuFloatingPointVec>;
using Wrapped1x1 = Opm::gpuistl::PreconditionerCPUMatrixToGPUMatrix<CuFloatingPointVec, CuFloatingPointVec, GpuDilu1x1, Sp1x1BlockMatrix>;
using GpuDilu2x2 = Opm::gpuistl::GpuDILU<Sp2x2BlockMatrix, CuFloatingPointVec, CuFloatingPointVec>;
using Wrapped2x2 = Opm::gpuistl::PreconditionerCPUMatrixToGPUMatrix<CuFloatingPointVec, CuFloatingPointVec, GpuDilu2x2, Sp2x2BlockMatrix>;

Sp1x1BlockMatrix
get1x1BlockTestMatrix()
{
    /*
        matA:
        1  2  0  3  0  0
        4  5  0  6  0  7
        0  0  8  0  0  0
        9 10  0 11 12  0
        0  0  0 13 14  0
        0 15  0  0  0 16

        Expected reordering:
        1  2  0  3  0  0
        0  0  8  0  0  0
        4  5  0  6  0  7
        9 10  0 11 12  0
        0 15  0  0  0 16
        0  0  0 13 14  0

        Expected lowerTriangularReorderedMatrix:
        0  0  0  0  0  0
        0  0  0  0  0  0
        4  0  0  0  0  0
        9 10  0  0  0  0
        0 15  0  0  0  0
        0  0  0 13  0  0

        Expected lowerTriangularReorderedMatrix:
        0  2  0  3  0  0
        0  0  0  0  0  0
        0  0  0  6  0  7
        0  0  0  0  12 0
        0  0  0  0  0  0
    */

    const int N = 6;
    const int nonZeroes = 16;

    // Create the Dune A matrix
    Sp1x1BlockMatrix matA(N, N, nonZeroes, Sp1x1BlockMatrix::row_wise);
    for (auto row = matA.createbegin(); row != matA.createend(); ++row) {
        row.insert(row.index());
        if (row.index() == 0) {
            row.insert(row.index() + 1);
            row.insert(row.index() + 3);
        }
        if (row.index() == 1) {
            row.insert(row.index() - 1);
            row.insert(row.index() + 2);
            row.insert(row.index() + 4);
        }
        if (row.index() == 2) {
        }
        if (row.index() == 3) {
            row.insert(row.index() - 3);
            row.insert(row.index() - 2);
            row.insert(row.index() + 1);
        }
        if (row.index() == 4) {
            row.insert(row.index() - 1);
        }
        if (row.index() == 5) {
            row.insert(row.index() - 4);
        }
    }

    matA[0][0][0][0] = 1.0;
    matA[0][1][0][0] = 2.0;
    matA[0][3][0][0] = 3.0;
    matA[1][0][0][0] = 4.0;
    matA[1][1][0][0] = 5.0;
    matA[1][3][0][0] = 6.0;
    matA[1][5][0][0] = 7.0;
    matA[2][2][0][0] = 8.0;
    matA[3][0][0][0] = 9.0;
    matA[3][1][0][0] = 10.0;
    matA[3][3][0][0] = 11.0;
    matA[3][4][0][0] = 12.0;
    matA[4][3][0][0] = 13.0;
    matA[4][4][0][0] = 14.0;
    matA[5][1][0][0] = 15.0;
    matA[5][5][0][0] = 16.0;

    return matA;
}

Sp2x2BlockMatrix
get2x2BlockTestMatrix()
{
    /*
    matA:
    1  2    0  3    0  0
    4  5    0  6    0  7

    0  0    1  0    0  0
    9 10    0  1   12  0

    0  0    0 13   14  0
    0 15    0  0    0 16

    */
    const int N = 3;
    const int nonZeroes = 9;

    // Create the Dune A matrix
    Sp2x2BlockMatrix matA(N, N, nonZeroes, Sp2x2BlockMatrix::row_wise);
    for (auto row = matA.createbegin(); row != matA.createend(); ++row) {
        row.insert(row.index());
        if (row.index() == 0) {
            row.insert(row.index() + 1);
            row.insert(row.index() + 2);
        }
        if (row.index() == 1) {
            row.insert(row.index() - 1);
            row.insert(row.index() + 1);
        }
        if (row.index() == 2) {
            row.insert(row.index() - 1);
            row.insert(row.index() - 2);
        }
    }

    matA[0][0][0][0] = 1.0;
    matA[0][0][0][1] = 2.0;
    matA[0][0][1][0] = 4.0;
    matA[0][0][1][1] = 5.0;
    matA[0][1][0][1] = 3.0;
    matA[0][1][1][1] = 6.0;
    matA[0][2][1][1] = 7.0;
    matA[1][0][1][0] = 9.0;
    matA[1][0][1][1] = 10.0;
    matA[1][1][0][0] = 1.0;
    matA[1][1][1][1] = 1.0;
    matA[1][2][1][0] = 12.0;
    matA[2][0][1][1] = 15.0;
    matA[2][1][0][1] = 13.0;
    matA[2][2][0][0] = 14.0;
    matA[2][2][1][1] = 16.0;

    return matA;
}

BOOST_AUTO_TEST_CASE(TestDiluApply)
{
    Sp1x1BlockMatrix matA = get1x1BlockTestMatrix();

    std::vector<double> input = {1.1, 1.2, 1.3, 1.4, 1.5, 1.6};
    std::vector<double> output(6);

    CuFloatingPointVec d_input(input);
    CuFloatingPointVec d_output(output);

    B1x1Vec h_input(6);
    h_input[0] = 1.1;
    h_input[1] = 1.2;
    h_input[2] = 1.3;
    h_input[3] = 1.4;
    h_input[4] = 1.5;
    h_input[5] = 1.6;
    B1x1Vec h_output(6);

    // Initialize preconditioner objects
    Dune::MultithreadDILU<Sp1x1BlockMatrix, B1x1Vec, B1x1Vec> cpudilu(matA);
    auto gpudilu = Wrapped1x1(matA, matA, true, true, false);

    // Use the apply
    gpudilu.apply(d_output, d_input);
    cpudilu.apply(h_output, h_input);

    // put results in std::vector
    std::vector<T> cpudilures;
    for (auto e : h_output) {
        cpudilures.push_back(e);
    }
    auto cudilures = d_output.asStdVector();

    // check that GpuDilu results matches that of CPU dilu
    for (size_t i = 0; i < cudilures.size(); ++i) {
        BOOST_CHECK_CLOSE(cudilures[i], cpudilures[i], 1e-7);
    }
}

BOOST_AUTO_TEST_CASE(TestDiluApplyBlocked)
{

    // init matrix with 2x2 blocks
    Sp2x2BlockMatrix matA = get2x2BlockTestMatrix();
    auto gpudilu = Wrapped2x2(matA, matA, true, true, false);
    Dune::MultithreadDILU<Sp2x2BlockMatrix, B2x2Vec, B2x2Vec> cpudilu(matA);

    // create input/output buffers for the apply
    std::vector<double> input = {1.1, 1.2, 1.3, 1.4, 1.5, 1.6};
    std::vector<double> output(6);
    CuFloatingPointVec d_input(input);
    CuFloatingPointVec d_output(output);

    B2x2Vec h_input(3);
    h_input[0][0] = 1.1;
    h_input[0][1] = 1.2;
    h_input[1][0] = 1.3;
    h_input[1][1] = 1.4;
    h_input[2][0] = 1.5;
    h_input[2][1] = 1.6;
    B2x2Vec h_output(3);

    // call apply with cpu and gpu dilu
    cpudilu.apply(h_output, h_input);
    gpudilu.apply(d_output, d_input);

    auto cudilures = d_output.asStdVector();
    std::vector<T> cpudilures;
    for (auto v : h_output) {
        for (auto e : v) {
            cpudilures.push_back(e);
        }
    }

    // check that the values are close
    for (size_t i = 0; i < cudilures.size(); ++i) {
        BOOST_CHECK_CLOSE(cudilures[i], cpudilures[i], 1e-7);
    }
}

BOOST_AUTO_TEST_CASE(TestDiluInitAndUpdateLarge)
{
    // create gpu dilu preconditioner
    Sp1x1BlockMatrix matA = get1x1BlockTestMatrix();
    auto gpudilu = Wrapped1x1(matA, matA, true, true, false);

    matA[0][0][0][0] = 11.0;
    matA[0][1][0][0] = 12.0;
    matA[0][3][0][0] = 13.0;
    matA[1][0][0][0] = 14.0;
    matA[1][1][0][0] = 15.0;
    matA[1][3][0][0] = 16.0;
    matA[1][5][0][0] = 17.0;
    matA[2][2][0][0] = 18.0;
    matA[3][0][0][0] = 19.0;
    matA[3][1][0][0] = 110.0;
    matA[3][3][0][0] = 111.0;
    matA[3][4][0][0] = 112.0;
    matA[4][3][0][0] = 113.0;
    matA[4][4][0][0] = 114.0;
    matA[5][1][0][0] = 115.0;
    matA[5][5][0][0] = 116.0;

    // make sure the function is updated
    gpudilu.update();
    // create a cpu dilu preconditioner on the matrix that is definitely updated
    Dune::MultithreadDILU<Sp1x1BlockMatrix, B1x1Vec, B1x1Vec> cpudilu(matA);

    std::vector<double> input = {1.1, 1.2, 1.3, 1.4, 1.5, 1.6};
    std::vector<double> output(6);

    CuFloatingPointVec d_input(input);
    CuFloatingPointVec d_output(output);

    B1x1Vec h_input(6);
    h_input[0] = 1.1;
    h_input[1] = 1.2;
    h_input[2] = 1.3;
    h_input[3] = 1.4;
    h_input[4] = 1.5;
    h_input[5] = 1.6;
    B1x1Vec h_output(6);

    // run an apply to see effect of update
    gpudilu.apply(d_output, d_input);
    cpudilu.apply(h_output, h_input);

    // put results in std::vector
    std::vector<T> cpudilures;
    for (auto e : h_output) {
        cpudilures.push_back(e);
    }
    auto cudilures = d_output.asStdVector();

    // check that GpuDilu results matches that of CPU dilu
    for (size_t i = 0; i < cudilures.size(); ++i) {
        BOOST_CHECK_CLOSE(cudilures[i], cpudilures[i], 1e-7);
    }
}
