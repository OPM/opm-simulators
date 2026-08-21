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

#define BOOST_TEST_MODULE TestGPUPointerAttributes

#include <boost/test/unit_test.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_pointer_attributes.hpp>
#include <opm/simulators/linalg/gpuistl/gpu_smart_pointer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>

# include <cuda.h>
BOOST_AUTO_TEST_CASE(TestGPUPointerAttributes)
{
    using namespace Opm::gpuistl::detail;

    int* hostPtr = nullptr;
    auto hostSmartPtr = std::make_unique<int>(1);
    auto devicePtr = Opm::gpuistl::make_gpu_unique_ptr<int>(1);
    auto devicePtrShared = Opm::gpuistl::make_gpu_shared_ptr<double>(23.0);

    std::vector<double> vectorOnCPU {{1, 2, 3, 4, 5, 6, 7}};
    auto vectorOnGPU = Opm::gpuistl::GpuVector<double>(vectorOnCPU.data(), vectorOnCPU.size());

    // Checking GPU pointers
    BOOST_CHECK(isGPUPointer(devicePtr));
    BOOST_CHECK(isGPUPointer(devicePtr.get()));
    BOOST_CHECK(!isGPUPointer(hostPtr));
    BOOST_CHECK(!isGPUPointer(hostSmartPtr.get()));
    BOOST_CHECK(!isGPUPointer(hostSmartPtr));
    BOOST_CHECK(isGPUPointer(devicePtrShared));
    BOOST_CHECK(isGPUPointer(vectorOnGPU.data()));
    BOOST_CHECK(isGPUPointer(vectorOnGPU.data() + 4));
    BOOST_CHECK(!isGPUPointer(vectorOnCPU.data()));
    BOOST_CHECK(!isGPUPointer(vectorOnCPU.data() + 4));

    // Checking CPU pointers
    BOOST_CHECK(!isCPUPointer(devicePtr));
    BOOST_CHECK(!isCPUPointer(devicePtr.get()));
    BOOST_CHECK(!isCPUPointer(hostPtr));
    BOOST_CHECK(isCPUPointer(hostSmartPtr.get()));
    BOOST_CHECK(isCPUPointer(hostSmartPtr));
    BOOST_CHECK(!isCPUPointer(devicePtrShared));
    BOOST_CHECK(!isCPUPointer(vectorOnGPU.data()));
    BOOST_CHECK(!isCPUPointer(vectorOnGPU.data() + 4));
    BOOST_CHECK(isCPUPointer(vectorOnCPU.data()));
    BOOST_CHECK(isCPUPointer(vectorOnCPU.data() + 4));
}
