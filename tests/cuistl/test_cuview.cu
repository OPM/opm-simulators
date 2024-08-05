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

#define BOOST_TEST_MODULE TestCuView

#include <boost/test/unit_test.hpp>
#include <cuda_runtime.h>
#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>
#include <opm/simulators/linalg/cuistl/CuView.hpp>
#include <opm/simulators/linalg/cuistl/CuBuffer.hpp>
#include <opm/simulators/linalg/cuistl/detail/cuda_safe_call.hpp>
#include <random>
#include <array>
#include <algorithm>
#include <type_traits>

using CuViewDouble = ::Opm::cuistl::CuView<double>;
using CuBufferDouble = ::Opm::cuistl::CuBuffer<double>;

__global__ void useCuViewOnGPU(CuViewDouble a, CuViewDouble b){
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
    auto cubuffer = CuBufferDouble(cpubuffer);
    auto cuview = CuViewDouble(cubuffer.data(), cubuffer.size());
    const auto const_cuview = CuViewDouble(cubuffer.data(), cubuffer.size());

    auto stdVecOfCuView = cuview.asStdVector();
    auto const_stdVecOfCuView = cuview.asStdVector();

    BOOST_CHECK_EQUAL_COLLECTIONS(
        stdVecOfCuView.begin(), stdVecOfCuView.end(), cpubuffer.begin(), cpubuffer.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(
        stdVecOfCuView.begin(), stdVecOfCuView.end(), const_stdVecOfCuView.begin(), const_stdVecOfCuView.end());
}

BOOST_AUTO_TEST_CASE(TestCuViewOnCPUTypes)
{
    auto buf = std::vector<double>({1.0, 2.0, 42.0, 59.9451743, 10.7132692});
    auto cpuview = CuViewDouble(buf.data(), buf.size());
    const auto const_cpuview = CuViewDouble(buf.data(), buf.size());

    // check that indexing a mutable view gives references when indexing it
    bool correct_type_of_cpu_front = std::is_same_v<double&, decltype(cpuview.front())>;
    bool correct_type_of_cpu_back = std::is_same_v<double&, decltype(cpuview.back())>;
    bool correct_type_of_const_cpu_front = std::is_same_v<double, decltype(const_cpuview.front())>;
    bool correct_type_of_const_cpu_back = std::is_same_v<double, decltype(const_cpuview.back())>;

    BOOST_CHECK(correct_type_of_cpu_front);
    BOOST_CHECK(correct_type_of_cpu_back);
    BOOST_CHECK(correct_type_of_const_cpu_front);
    BOOST_CHECK(correct_type_of_const_cpu_back);

    // check that the values are correct
    BOOST_CHECK(cpuview.front() == buf.front());
    BOOST_CHECK(cpuview.back() == buf.back());
}

BOOST_AUTO_TEST_CASE(TestCuViewOnCPUWithSTLIteratorAlgorithm)
{
    auto buf = std::vector<double>({1.0, 2.0, 42.0, 59.9451743, 10.7132692});
    auto cpuview = CuViewDouble(buf.data(), buf.size());
    std::sort(buf.begin(), buf.end());
    BOOST_CHECK(42.0 == cpuview[3]);
}

BOOST_AUTO_TEST_CASE(TestCuViewOnGPU)
{
    auto buf = std::vector<double>({1.0, 2.0, 42.0, 59.9451743, 10.7132692});
    auto cubufA = CuBufferDouble(buf);
    auto cuviewA = CuViewDouble(cubufA.data(), cubufA.size());
    auto cubufB = CuBufferDouble(4);
    auto cuviewB = CuViewDouble(cubufB.data(), cubufB.size());

    useCuViewOnGPU<<<1,1>>>(cuviewA, cuviewB);

    auto vecA = cuviewA.asStdVector();
    auto vecB = cuviewB.asStdVector();

    // checks that front/back/begin/end works
    BOOST_CHECK(vecB[0] == buf[0]);
    BOOST_CHECK(vecB[1] == buf[4]);
    BOOST_CHECK(vecB[2] == buf[0]);
    BOOST_CHECK(vecB[3] == buf[4]);

    // checks that view[0] = view[2] works
    BOOST_CHECK(buf[2] == vecA[0]);
}
