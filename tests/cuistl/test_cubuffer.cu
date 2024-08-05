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

#define BOOST_TEST_MODULE TestCuBuffer

#include <boost/test/unit_test.hpp>
#include <cuda_runtime.h>

#include <opm/simulators/linalg/cuistl/CuBuffer.hpp>
#include <opm/simulators/linalg/cuistl/CuView.hpp>
#include <opm/simulators/linalg/cuistl/detail/cuda_safe_call.hpp>

#include <array>
#include <algorithm>
#include <type_traits>

BOOST_AUTO_TEST_CASE(TestMakeView)
{
    // test that we can create buffers and make views of the buffers using the pointer constructor
    auto buf = std::vector<int>({1, 2, 3, 4, 5, 6});
    const auto gpubuf = ::Opm::cuistl::CuBuffer<int>(buf);
    auto gpuview = ::Opm::cuistl::CuView<int>(buf.data(), buf.size());
    bool gpuBufCreatedView = std::is_same<::Opm::cuistl::CuView<int>, decltype(gpuview)>::value;

    BOOST_CHECK(gpuBufCreatedView);

    // test that we can make views of buffers by using the cubuffer constructor
    auto gpuview2 = ::Opm::cuistl::make_view(gpubuf);
    bool gpuBufCreatedView2 = std::is_same<::Opm::cuistl::CuView<const int>, decltype(gpuview2)>::value;

    BOOST_CHECK(gpuBufCreatedView2);

    // check that we retrieve the same values when pulling the data back to the cpu as a vector
    auto gpuBufOnCpu = gpubuf.asStdVector();
    BOOST_CHECK_EQUAL_COLLECTIONS(gpuBufOnCpu.begin(), gpuBufOnCpu.end(), buf.begin(), buf.end());
}
