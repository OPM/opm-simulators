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
#include <opm/common/utility/gpuDecorators.hpp>

#define BOOST_TEST_MODULE TestGPUMacros


#include <boost/test/unit_test.hpp>
#include <opm/simulators/linalg/gpuistl/gpu_smart_pointer.hpp>
namespace
{

__global__ void testKernelSetIfGPU(int* ptr) {
    #if OPM_IS_INSIDE_DEVICE_FUNCTION
    *ptr = 42;
    #else
    *ptr = 123;
    #endif
}

__global__ void testKernelSetIfNotCPU(int* ptr) {
    #if OPM_IS_INSIDE_HOST_FUNCTION
    *ptr = 42;
    #else
    *ptr = 123;
    #endif
}


} // namespace


BOOST_AUTO_TEST_CASE(TestInsideDevice)
{
    auto sharedPtr = Opm::gpuistl::make_gpu_shared_ptr<int>(1);

    testKernelSetIfGPU<<<1, 1>>>(sharedPtr.get());
    auto valueFromDevice = Opm::gpuistl::copyFromGPU(sharedPtr);
    BOOST_CHECK_EQUAL(valueFromDevice, 42);
}

BOOST_AUTO_TEST_CASE(TestOutsideHost)
{
    auto sharedPtr = Opm::gpuistl::make_gpu_shared_ptr<int>(1);

    testKernelSetIfNotCPU<<<1, 1>>>(sharedPtr.get());
    auto valueFromDevice = Opm::gpuistl::copyFromGPU(sharedPtr);
    BOOST_CHECK_EQUAL(valueFromDevice, 123);
}

BOOST_AUTO_TEST_CASE(TestMacrosOnHost)
{
    #if OPM_IS_INSIDE_HOST_FUNCTION
    BOOST_CHECK(true);
    #else
    BOOST_CHECK(false);
    #endif

    #if OPM_IS_INSIDE_DEVICE_FUNCTION
    BOOST_CHECK(false);
    #else
    BOOST_CHECK(true);
    #endif
}