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

#define BOOST_TEST_MODULE TestGpuBuffer

#include <boost/test/unit_test.hpp>
#include <cuda_runtime.h>

#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>

#include <array>
#include <algorithm>
#include <type_traits>

BOOST_AUTO_TEST_CASE(TestMakeView)
{
    // check creation of buffers and views for mutable buffers
    auto buf = std::vector<int>({1, 2, 3, 4, 5, 6});
    auto gpubuf = ::Opm::gpuistl::GpuBuffer<int>(buf);
    auto gpuview = ::Opm::gpuistl::make_view(gpubuf);
    bool gpuBufCreatedView = std::is_same_v<::Opm::gpuistl::GpuView<int>, decltype(gpuview)>;
    BOOST_CHECK(gpuBufCreatedView);

    auto gpubufOnCpu = gpubuf.asStdVector();
    BOOST_CHECK_EQUAL_COLLECTIONS(gpubufOnCpu.begin(), gpubufOnCpu.end(), buf.begin(), buf.end());

    // check creation of buffers and views for const buffers
    const auto buf2 = std::vector<int>({2, 3, 4, 5, 6});
    const auto gpubuf2 = ::Opm::gpuistl::GpuBuffer<int>(buf2);
    auto gpuview2 = ::Opm::gpuistl::make_view(gpubuf2);
    bool gpuBufCreatedView2 = std::is_same_v<::Opm::gpuistl::GpuView<const int>, decltype(gpuview2)>;
    BOOST_CHECK(gpuBufCreatedView2);

    auto gpubufOnCpu2 = gpubuf2.asStdVector();
    BOOST_CHECK_EQUAL_COLLECTIONS(gpubufOnCpu2.begin(), gpubufOnCpu2.end(), buf2.begin(), buf2.end());
}
