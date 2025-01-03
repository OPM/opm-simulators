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

#define BOOST_TEST_MODULE TestGpuView

#include <boost/test/unit_test.hpp>
#include <cuda_runtime.h>
#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>
#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>
#include <random>
#include <array>
#include <algorithm>
#include <type_traits>

using GpuViewDouble = ::Opm::gpuistl::GpuView<double>;
using GpuBufferDouble = ::Opm::gpuistl::GpuBuffer<double>;

__global__ void useGpuViewOnGPU(GpuViewDouble a, GpuViewDouble b){
    b[0] = a.front();
    b[1] = a.back();
    b[2] = *a.begin();
    b[3] = *(--a.end());

    a[0] = a[2];
}

BOOST_AUTO_TEST_CASE(TestCreationAndIndexing)
{
    // A simple test to check that we can move data to and from the GPU
    auto cpubuffer = std::vector<double>({1.0, 2.0, 42.0, 59.9451743, 10.7132692});
    auto cubuffer = GpuBufferDouble(cpubuffer);
    auto gpuview = GpuViewDouble(cubuffer.data(), cubuffer.size());
    const auto const_gpuview = GpuViewDouble(cubuffer.data(), cubuffer.size());

    auto stdVecOfGpuView = gpuview.asStdVector();
    auto const_stdVecOfGpuView = gpuview.asStdVector();

    BOOST_CHECK_EQUAL_COLLECTIONS(
        stdVecOfGpuView.begin(), stdVecOfGpuView.end(), cpubuffer.begin(), cpubuffer.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(
        stdVecOfGpuView.begin(), stdVecOfGpuView.end(), const_stdVecOfGpuView.begin(), const_stdVecOfGpuView.end());
}

BOOST_AUTO_TEST_CASE(TestGpuViewOnCPUTypes)
{
    auto buf = std::vector<double>({1.0, 2.0, 42.0, 59.9451743, 10.7132692});
    auto cpuview = GpuViewDouble(buf.data(), buf.size());
    const auto const_cpuview = GpuViewDouble(buf.data(), buf.size());

    // check that indexing a const view produces a value
    bool correct_type_of_const_cpu_front = std::is_same_v<double, decltype(const_cpuview.front())>;
    bool correct_type_of_const_cpu_back = std::is_same_v<double, decltype(const_cpuview.back())>;

    BOOST_CHECK(correct_type_of_const_cpu_front);
    BOOST_CHECK(correct_type_of_const_cpu_back);

    // check that the values are correct
    BOOST_CHECK(const_cpuview.front() == buf.front());
    BOOST_CHECK(const_cpuview.back() == buf.back());
}

BOOST_AUTO_TEST_CASE(TestGpuViewOnCPUWithSTLIteratorAlgorithm)
{
    auto buf = std::vector<double>({1.0, 2.0, 42.0, 59.9451743, 10.7132692});
    auto cpuview = GpuViewDouble(buf.data(), buf.size());
    std::sort(buf.begin(), buf.end());
    BOOST_CHECK(42.0 == cpuview[3]);
}

BOOST_AUTO_TEST_CASE(TestGpuViewOnGPU)
{
    auto buf = std::vector<double>({1.0, 2.0, 42.0, 59.9451743, 10.7132692});
    auto cubufA = GpuBufferDouble(buf);
    auto gpuviewA = GpuViewDouble(cubufA.data(), cubufA.size());
    auto cubufB = GpuBufferDouble(4);
    auto gpuviewB = GpuViewDouble(cubufB.data(), cubufB.size());

    useGpuViewOnGPU<<<1,1>>>(gpuviewA, gpuviewB);

    auto vecA = gpuviewA.asStdVector();
    auto vecB = gpuviewB.asStdVector();

    // checks that front/back/begin/end works
    BOOST_CHECK(vecB[0] == buf[0]);
    BOOST_CHECK(vecB[1] == buf[4]);
    BOOST_CHECK(vecB[2] == buf[0]);
    BOOST_CHECK(vecB[3] == buf[4]);

    // checks that view[0] = view[2] works
    BOOST_CHECK(buf[2] == vecA[0]);
}
