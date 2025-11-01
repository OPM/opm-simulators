/*
  Copyright 2025 Equinor ASA
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

#include "config.h"
#include "opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp"

#define BOOST_TEST_MODULE OpmVectorTest
#include <boost/test/included/unit_test.hpp>

#include <opm/simulators/linalg/gpuistl/MiniVector.hpp>
#include <opm/simulators/linalg/gpuistl/gpu_smart_pointer.hpp>
#include <cuda_runtime.h>

BOOST_AUTO_TEST_CASE(DefaultConstructor)
{
    ::Opm::gpuistl::MiniVector<int, 3> v; // Default constructor initializes all elements to zero
    BOOST_TEST(v.size() == 3u);
    for (std::size_t i = 0; i < v.size(); ++i) {
        BOOST_TEST(v[i] == 0);
    }
}

BOOST_AUTO_TEST_CASE(UniformValueConstructor)
{
    constexpr float kVal = 2.5f;
    ::Opm::gpuistl::MiniVector<float, 4> v(kVal);
    for (auto e : v) {
        BOOST_TEST(e == kVal, boost::test_tools::tolerance(1e-6f));
    }
}

BOOST_AUTO_TEST_CASE(InitializerListConstructor)
{
    ::Opm::gpuistl::MiniVector<double, 2> v {1.1, 2.2};
    BOOST_TEST(v[0] == 1.1, boost::test_tools::tolerance(1e-12));
    BOOST_TEST(v[1] == 2.2, boost::test_tools::tolerance(1e-12));
}

BOOST_AUTO_TEST_CASE(AtBoundsChecking)
{
    ::Opm::gpuistl::MiniVector<int, 2> v {5, 6};
    BOOST_CHECK_THROW(v.at(2), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(IteratorTraversal)
{
    ::Opm::gpuistl::MiniVector<int, 3> v {1, 2, 3};
    int sum = 0;
    for (int e : v)
        sum += e;
    BOOST_TEST(sum == 6);
}

template <class VecT>
__global__ void
fillKernel(VecT* v, typename VecT::value_type val)
{
    (*v)[threadIdx.x] = val;
}

template <class VecT, class IntT>
__global__ void
sumKernel(const VecT* v, IntT* out)
{
    if (threadIdx.x == 0) {
        IntT s = 0;
        for (std::size_t i = 0; i < VecT::size(); ++i) {
            s += (*v)[i];
        }
        *out = s;
    }
}

BOOST_AUTO_TEST_CASE(GPUFill)
{
    // make_gpu_shared_ptr<T>() should:
    //   1. allocate device memory for a single T
    //   2. return a std::shared_ptr<T> whose custom deleter frees GPU memory
    //   3. optionally take a host value to copy to device (overload)

    auto dvec = ::Opm::gpuistl::make_gpu_shared_ptr<::Opm::gpuistl::MiniVector<int, 3>>();

    // Fill the vector from the device
    fillKernel<<<1, 3>>>(dvec.get(), 7);
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());

    // Copy back and verify
    ::Opm::gpuistl::MiniVector<int, 3> hvec = ::Opm::gpuistl::copyFromGPU(dvec);
    for (int i = 0; i < 3; ++i) {
        BOOST_TEST(hvec[i] == 7);
    }
}

BOOST_AUTO_TEST_CASE(GPUSum)
{
    ::Opm::gpuistl::MiniVector<int, 4> hsrc(1); // [1, 1, 1, 1]

    auto dsrc = ::Opm::gpuistl::make_gpu_shared_ptr<::Opm::gpuistl::MiniVector<int, 4>>(hsrc);
    auto dout = ::Opm::gpuistl::make_gpu_shared_ptr<int>();

    sumKernel<<<1, 1>>>(dsrc.get(), dout.get());
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());

    int sum = ::Opm::gpuistl::copyFromGPU(dout);
    BOOST_TEST(sum == 4);
}
