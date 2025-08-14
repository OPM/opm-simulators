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

#define BOOST_TEST_MODULE TestConditionalStorageGPU

#include <boost/test/unit_test.hpp>
#include <cuda.h>
#include <cuda_runtime.h>
#include <opm/material/common/ConditionalStorage.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>
#include <opm/simulators/linalg/gpuistl/gpu_smart_pointer.hpp>

namespace
{

template <bool enabled, class T>
__global__ void
createCondtionalStorage()
{
    // just make sure we can run the constructor
    Opm::ConditionalStorage<enabled, T> someStorage;


    Opm::ConditionalStorage<enabled, T> other = someStorage;

    other = someStorage;
}
template <class T>
__global__ void
testEnabledStorage(Opm::ConditionalStorage<true, T> storage, T* output)
{
    output[0] = *storage;
}

template <class T, class S>
__global__ void
testEnabledStorageArrow(Opm::ConditionalStorage<true, T> storage, S* output)
{
    output[0] = storage->someFunc();
}

template <class T>
__global__ void
testEnabledStorage(Opm::ConditionalStorage<true, T>* storage, T* output)
{
    output[0] = **storage;
}

template <class T, class S>
__global__ void
testEnabledStorageArrow(Opm::ConditionalStorage<true, T>* storage, S* output)
{
    output[0] = (*storage)->someFunc();
}

struct SomeStruct {
    OPM_HOST_DEVICE int someFunc()
    {
        return 123;
    }
};

} // namespace

BOOST_AUTO_TEST_CASE(TestRunConstructor)
{
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
    OPM_GPU_SAFE_CALL(cudaGetLastError());
    createCondtionalStorage<true, double><<<1, 1>>>();
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
    OPM_GPU_SAFE_CALL(cudaGetLastError());
    createCondtionalStorage<false, double><<<1, 1>>>();
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
    OPM_GPU_SAFE_CALL(cudaGetLastError());
}

BOOST_AUTO_TEST_CASE(TestEnabledStoragePointer)
{
    using namespace Opm;
    using CS = ConditionalStorage<true, double>;
    auto storage = Opm::gpuistl::make_gpu_unique_ptr<CS>(CS(32.2));
    auto numberFromGPU = Opm::gpuistl::make_gpu_unique_ptr<double>(0.0);
    OPM_GPU_SAFE_CALL(cudaGetLastError());
    testEnabledStorage<<<1, 1>>>(storage.get(), numberFromGPU.get());
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
    OPM_GPU_SAFE_CALL(cudaGetLastError());
    auto number = Opm::gpuistl::copyFromGPU(numberFromGPU);
    BOOST_CHECK_EQUAL(32.2, number);

    auto numberFromGPUFromCall = Opm::gpuistl::make_gpu_unique_ptr<int>(0);

    auto storageSomeStruct = Opm::gpuistl::make_gpu_unique_ptr<ConditionalStorage<true, SomeStruct>>(
        ConditionalStorage<true, SomeStruct>());
    testEnabledStorageArrow<<<1, 1>>>(storageSomeStruct.get(), numberFromGPUFromCall.get());
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
    OPM_GPU_SAFE_CALL(cudaGetLastError());
    auto numberFromCall = Opm::gpuistl::copyFromGPU(numberFromGPUFromCall);
    BOOST_CHECK_EQUAL(123, numberFromCall);
}
BOOST_AUTO_TEST_CASE(TestEnabledStorageCopy)
{
    using namespace Opm;
    using CS = ConditionalStorage<true, double>;
    auto storage = CS(32.2);
    auto numberFromGPU = Opm::gpuistl::make_gpu_unique_ptr<double>(0.0);
    testEnabledStorage<<<1, 1>>>(storage, numberFromGPU.get());
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
    OPM_GPU_SAFE_CALL(cudaGetLastError());
    auto number = Opm::gpuistl::copyFromGPU(numberFromGPU);
    BOOST_CHECK_EQUAL(32.2, number);

    auto numberFromGPUFromCall = Opm::gpuistl::make_gpu_unique_ptr<int>(0);

    auto storageSomeStruct = ConditionalStorage<true, SomeStruct>();
    testEnabledStorageArrow<<<1, 1>>>(storageSomeStruct, numberFromGPUFromCall.get());
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
    OPM_GPU_SAFE_CALL(cudaGetLastError());
    auto numberFromCall = Opm::gpuistl::copyFromGPU(numberFromGPUFromCall);
    BOOST_CHECK_EQUAL(123, numberFromCall);
}
