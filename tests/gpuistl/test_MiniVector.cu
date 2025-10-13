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

#define BOOST_TEST_MODULE TestMiniVector

#include <boost/test/unit_test.hpp>
#include <cuda_runtime.h>
#include <opm/simulators/linalg/gpuistl/MiniVector.hpp>

template<typename VecType>
__global__ void doNothingKernel(VecType v)
{
    auto idx = threadIdx.x;
    return;
}

// Check that we can create a minivector and pass it to a GPU kernel without errors
BOOST_AUTO_TEST_CASE(TestPassingToKernel)
{
    Opm::gpuistl::MiniVector<double, 3> v = {1.0, 2.0, 3.0};
    doNothingKernel<<<1, 1>>>(v);
    cudaDeviceSynchronize();
    cudaError_t err = cudaGetLastError();
    BOOST_CHECK(err == cudaSuccess);
}

template <typename VecType>
__global__ void kernelWithVectorOperations(VecType v1, VecType v2)
{
    assert(v1.size() == v2.size());
    assert(v1 != v2);
    assert(v1.at(0) == v2.at(0));
    assert(v1[0] == v2[0]);
}

BOOST_AUTO_TEST_CASE(TestVectorOperationsOnDevice)
{
    Opm::gpuistl::MiniVector<double, 3> v1 = {1.0, 2.0, 3.0};
    Opm::gpuistl::MiniVector<double, 3> v2(1.0);
    kernelWithVectorOperations<<<1, 1>>>(v1, v2);
    cudaDeviceSynchronize();
    cudaError_t err = cudaGetLastError();
    BOOST_CHECK(err == cudaSuccess);
}
