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
#include <config.h>

#define BOOST_TEST_MODULE TestIsGPUPointer

#include <boost/test/unit_test.hpp>
#include <opm/simulators/linalg/gpuistl/detail/is_gpu_pointer.hpp>
#include <opm/simulators/linalg/gpuistl/gpu_smart_pointer.hpp>
BOOST_AUTO_TEST_CASE(TestIsGPUPointer)
{
    using namespace Opm::gpuistl::detail;

    int* hostPtr = nullptr;
    auto hostSmartPtr = std::make_unique<int>(1);
    auto devicePtr = Opm::gpuistl::make_gpu_unique_ptr<int>(1);
    auto devicePtrShared = Opm::gpuistl::make_gpu_shared_ptr<double>(23.0);

    BOOST_CHECK(isGPUPointer(devicePtr));
    BOOST_CHECK(isGPUPointer(devicePtr.get()));
    BOOST_CHECK(!isGPUPointer(hostPtr));
    BOOST_CHECK(!isGPUPointer(hostSmartPtr.get()));
    BOOST_CHECK(!isGPUPointer(hostSmartPtr));
    BOOST_CHECK(isGPUPointer(devicePtrShared));
}
