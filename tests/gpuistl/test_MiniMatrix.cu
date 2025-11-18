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
#include <opm/simulators/linalg/gpuistl/MiniVector.hpp>
#include <opm/simulators/linalg/gpuistl/gpu_smart_pointer.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>

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

__global__ void MiniMatrixOperationsInKernel(MatType m1, MatType m2, bool* result)
{
    *result = true;
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
        *result &= (*it1 == *it2);
    }

    return;
}

BOOST_AUTO_TEST_CASE(TestMiniMatrixOperationsDefined)
{
    MatType m1 = 1.0;
    MatType m2 = 1.0;
    auto d_result = Opm::gpuistl::make_gpu_unique_ptr<bool>(true);

    MiniMatrixOperationsInKernel<<<1, 1>>>(m1, m2, d_result.get());
    std::ignore = cudaDeviceSynchronize();
    cudaError_t err = cudaGetLastError();
    BOOST_CHECK(err == cudaSuccess);

    bool h_result = Opm::gpuistl::copyFromGPU(d_result.get());
    BOOST_CHECK(h_result);
}

__global__ void WriteToMatrixInKernel(MatType* m)
{
    (*m) = MatType(1.0);
    (*m)[1][1] = 3.14;
    return;
}

BOOST_AUTO_TEST_CASE(TestWritingToMatrixInKernel)
{
    auto d_m = Opm::gpuistl::make_gpu_unique_ptr<MatType>();

    WriteToMatrixInKernel<<<1, 1>>>(d_m.get());
    std::ignore = cudaDeviceSynchronize();
    cudaError_t err = cudaGetLastError();
    BOOST_CHECK(err == cudaSuccess);

    MatType h_m = Opm::gpuistl::copyFromGPU(d_m.get());

    BOOST_CHECK(h_m[1][1] == 3.14);
}

BOOST_AUTO_TEST_CASE(TestMiniMatrixOperators)
{
    // A bit verbose but at least everything that is tested is clear
    MatType m1 = {
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
        7.0, 8.0, 9.0
    };

    MatType m2 = {
        1.0, 3.0, 2.0,
        2.0, 2.0, 3.0,
        3.0, 1.0, 1.0
    };

    // test +, - and * matrix operators
    MatType m3 = m1 + m2;
    BOOST_CHECK(m3[0][0] == 2.0);
    BOOST_CHECK(m3[0][1] == 5.0);
    BOOST_CHECK(m3[0][2] == 5.0);
    BOOST_CHECK(m3[1][0] == 6.0);
    BOOST_CHECK(m3[1][1] == 7.0);
    BOOST_CHECK(m3[1][2] == 9.0);
    BOOST_CHECK(m3[2][0] == 10.0);
    BOOST_CHECK(m3[2][1] == 9.0);
    BOOST_CHECK(m3[2][2] == 10.0);

    MatType m4 = m1 - m2;
    BOOST_CHECK(m4[0][0] == 0.0);
    BOOST_CHECK(m4[0][1] == -1.0);
    BOOST_CHECK(m4[0][2] == 1.0);
    BOOST_CHECK(m4[1][0] == 2.0);
    BOOST_CHECK(m4[1][1] == 3.0);
    BOOST_CHECK(m4[1][2] == 3.0);
    BOOST_CHECK(m4[2][0] == 4.0);
    BOOST_CHECK(m4[2][1] == 7.0);
    BOOST_CHECK(m4[2][2] == 8.0);

    MatType m5 = m1 * m2;
    BOOST_CHECK(m5[0][0] == 14.0);
    BOOST_CHECK(m5[0][1] == 10.0);
    BOOST_CHECK(m5[0][2] == 11.0);
    BOOST_CHECK(m5[1][0] == 32.0);
    BOOST_CHECK(m5[1][1] == 28.0);
    BOOST_CHECK(m5[1][2] == 29.0);
    BOOST_CHECK(m5[2][0] == 50.0);
    BOOST_CHECK(m5[2][1] == 46.0);
    BOOST_CHECK(m5[2][2] == 47.0);

    // test +=, -=, *= matrix operators
    MatType m6 = m1;
    m6 += m2;
    BOOST_CHECK(m6[0][0] == 2.0);
    BOOST_CHECK(m6[0][1] == 5.0);
    BOOST_CHECK(m6[0][2] == 5.0);
    BOOST_CHECK(m6[1][0] == 6.0);
    BOOST_CHECK(m6[1][1] == 7.0);
    BOOST_CHECK(m6[1][2] == 9.0);
    BOOST_CHECK(m6[2][0] == 10.0);
    BOOST_CHECK(m6[2][1] == 9.0);
    BOOST_CHECK(m6[2][2] == 10.0);

    MatType m7 = m1;;
    m7 -= m2;
    BOOST_CHECK(m7[0][0] == 0.0);
    BOOST_CHECK(m7[0][1] == -1.0);
    BOOST_CHECK(m7[0][2] == 1.0);
    BOOST_CHECK(m7[1][0] == 2.0);
    BOOST_CHECK(m7[1][1] == 3.0);
    BOOST_CHECK(m7[1][2] == 3.0);
    BOOST_CHECK(m7[2][0] == 4.0);
    BOOST_CHECK(m7[2][1] == 7.0);
    BOOST_CHECK(m7[2][2] == 8.0);

    MatType m8 = m1;;
    m8 *= m2;
    BOOST_CHECK(m8[0][0] == 14.0);
    BOOST_CHECK(m8[0][1] == 10.0);
    BOOST_CHECK(m8[0][2] == 11.0);
    BOOST_CHECK(m8[1][0] == 32.0);
    BOOST_CHECK(m8[1][1] == 28.0);
    BOOST_CHECK(m8[1][2] == 29.0);
    BOOST_CHECK(m8[2][0] == 50.0);
    BOOST_CHECK(m8[2][1] == 46.0);
    BOOST_CHECK(m8[2][2] == 47.0);
}

BOOST_AUTO_TEST_CASE(TestMiniMatrixWithMiniVector)
{
    MatType m = {
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
        7.0, 8.0, 9.0
    };

    Opm::gpuistl::MiniVector<double, 3> v = {1.0, 2.0, 3.0};

    auto result = m * v;

    BOOST_CHECK(result[0] == 14.0);
    BOOST_CHECK(result[1] == 32.0);
    BOOST_CHECK(result[2] == 50.0);
}
