/*
  Copyright 2025 EQUINOR ASA

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

#define BOOST_TEST_MODULE TestMiniMatrix

#include <boost/test/unit_test.hpp>
#include <cuda_runtime.h>
#include <opm/simulators/linalg/gpuistl/MiniMatrix.hpp>
#include <utility> // for std::ignore

using MatType = Opm::gpuistl::MiniMatrix<double, 3>;

__global__ void doNothingKernel(MatType m)
{
    auto idx = threadIdx.x;
    return;
}

BOOST_AUTO_TEST_CASE(TestPassingMatrixToKernel)
{
    MatType m;
    doNothingKernel<<<1, 1>>>(m);
    std::ignore = cudaDeviceSynchronize();
    cudaError_t err = cudaGetLastError();
    BOOST_CHECK(err == cudaSuccess);
}

__global__ void MiniMatrixOperationsInKernel(MatType m1, MatType m2)
{
    {
      auto tmp = m1 * m2;
    }
    {
      auto tmp = m1 + m2;
    }
    {
      auto tmp = m1 - m2;
    }
    m1 += m2;
    m1 -= m2;

    // Check equality and iterator usage
    for (auto it1 = m1.begin(), it2 = m2.begin(); it1 != m1.end() && it2 != m2.end(); ++it1, ++it2)
    {
        assert(*it1 == *it2);
    }

    return;
}

BOOST_AUTO_TEST_CASE(TestMiniMatrixOperationsInKernel)
{
    MatType m1 = 1.0;
    MatType m2 = 1.0;
    MiniMatrixOperationsInKernel<<<1, 1>>>(m1, m2);
    std::ignore = cudaDeviceSynchronize();
    cudaError_t err = cudaGetLastError();
    BOOST_CHECK(err == cudaSuccess);
}

__global__ void WriteToMatrixInKernel(MatType* m)
{
    (*m) = MatType(1.0);
    (*m)[1][1] = 3.14;
    return;
}

BOOST_AUTO_TEST_CASE(TestWritingToMatrixInKernel)
{
    MatType* d_m = nullptr;
    std::ignore = cudaMalloc(&d_m, sizeof(MatType));

    WriteToMatrixInKernel<<<1, 1>>>(d_m);
    std::ignore = cudaDeviceSynchronize();
    cudaError_t err = cudaGetLastError();
    BOOST_CHECK(err == cudaSuccess);

    MatType h_m;
    std::ignore = cudaMemcpy(&h_m, d_m, sizeof(MatType), cudaMemcpyDeviceToHost);
    std::ignore = cudaFree(d_m);

    BOOST_CHECK(h_m[1][1] == 3.14);
}
