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

#define BOOST_TEST_MODULE TestSmartPointers

#include <array>

#include <boost/test/unit_test.hpp>
#include <opm/simulators/linalg/gpuistl/gpu_smart_pointer.hpp>
namespace
{

struct SomeStruct {
    __device__ void someFunction()
    {
        this->isCalled = true;
    }

    bool isCalled = false;
};

template <class T>
__global__ void
setValue(Opm::gpuistl::PointerView<T> ptrIn, Opm::gpuistl::PointerView<T> ptrOut)
{
    *ptrOut = *ptrIn;
}

template <class T>
__global__ void
setValueGet(Opm::gpuistl::PointerView<T> ptrIn, Opm::gpuistl::PointerView<T> ptrOut)
{
    *ptrOut.get() = *ptrIn.get();
}

template <class T>
__global__ void
callFunction(Opm::gpuistl::PointerView<T> ptrIn)
{
    ptrIn->someFunction();
}


} // namespace


BOOST_AUTO_TEST_CASE(TestSharedPointer)
{
    auto sharedPtr = Opm::gpuistl::make_gpu_shared_ptr<int>(1);

    int valueFromDevice = 0;
    OPM_GPU_SAFE_CALL(cudaMemcpy(&valueFromDevice, sharedPtr.get(), sizeof(int), cudaMemcpyDeviceToHost));
    BOOST_CHECK_EQUAL(valueFromDevice, 1);
}

BOOST_AUTO_TEST_CASE(TestUniquePointer)
{
    auto uniquePtr = Opm::gpuistl::make_gpu_unique_ptr<int>(1);

    int valueFromDevice = 0;
    OPM_GPU_SAFE_CALL(cudaMemcpy(&valueFromDevice, uniquePtr.get(), sizeof(int), cudaMemcpyDeviceToHost));
    BOOST_CHECK_EQUAL(valueFromDevice, 1);
}

BOOST_AUTO_TEST_CASE(TestPointerView)
{
    auto pointerDestination = Opm::gpuistl::make_gpu_shared_ptr<double>(92);
    auto pointerSource = Opm::gpuistl::make_gpu_shared_ptr<double>(128.5);

    setValue<<<1, 1>>>(Opm::gpuistl::make_view(pointerSource), Opm::gpuistl::make_view(pointerDestination));

    double valueFromDevice = 0;
    OPM_GPU_SAFE_CALL(cudaMemcpy(&valueFromDevice, pointerDestination.get(), sizeof(double), cudaMemcpyDeviceToHost));
    BOOST_CHECK_EQUAL(valueFromDevice, 128.5);

    auto newSource = Opm::gpuistl::make_gpu_shared_ptr<double>(-1.0);
    setValueGet<<<1, 1>>>(Opm::gpuistl::make_view(newSource), Opm::gpuistl::make_view(pointerDestination));
    OPM_GPU_SAFE_CALL(cudaMemcpy(&valueFromDevice, pointerDestination.get(), sizeof(double), cudaMemcpyDeviceToHost));
    BOOST_CHECK_EQUAL(valueFromDevice, -1.0);

    auto structPtr = Opm::gpuistl::make_gpu_shared_ptr<SomeStruct>();
    callFunction<<<1, 1>>>(Opm::gpuistl::make_view(structPtr));
    bool isCalled = false;
    OPM_GPU_SAFE_CALL(cudaMemcpy(&isCalled, structPtr.get(), sizeof(bool), cudaMemcpyDeviceToHost));
    BOOST_CHECK_EQUAL(isCalled, true);


    auto uniquePtr = Opm::gpuistl::make_gpu_unique_ptr<double>(1.0);
    auto uniqueView = Opm::gpuistl::make_view(uniquePtr);

    double valueFromDeviceUnique = 0;
    OPM_GPU_SAFE_CALL(cudaMemcpy(&valueFromDeviceUnique, uniqueView.get(), sizeof(double), cudaMemcpyDeviceToHost));
    BOOST_CHECK_EQUAL(valueFromDeviceUnique, 1.0);
}

BOOST_AUTO_TEST_CASE(TestCopyFromGPU)
{
    auto sharedPtr = Opm::gpuistl::make_gpu_shared_ptr<int>(42);
    auto fromGPU = Opm::gpuistl::copyFromGPU(sharedPtr);
    BOOST_CHECK_EQUAL(fromGPU, 42);

    auto uniquePtr = Opm::gpuistl::make_gpu_unique_ptr<int>(128);
    auto fromGPUUnique = Opm::gpuistl::copyFromGPU(uniquePtr);
    BOOST_CHECK_EQUAL(fromGPUUnique, 128);

    auto sharedPtrArray
        = Opm::gpuistl::make_gpu_shared_ptr<std::array<double, 4>>(std::array<double, 4> {1.0, 2.0, 3.0, 4.0});
    auto fromGPUArray = Opm::gpuistl::copyFromGPU(sharedPtrArray);
    BOOST_CHECK_EQUAL(fromGPUArray[0], 1.0);
    BOOST_CHECK_EQUAL(fromGPUArray[1], 2.0);
    BOOST_CHECK_EQUAL(fromGPUArray[2], 3.0);
    BOOST_CHECK_EQUAL(fromGPUArray[3], 4.0);
}

BOOST_AUTO_TEST_CASE(TestCopyToGPU) {
    auto sharedPtr = Opm::gpuistl::make_gpu_shared_ptr<double>();
    Opm::gpuistl::copyToGPU(42.0, sharedPtr);
    auto valueFromDevice = Opm::gpuistl::copyFromGPU(sharedPtr);
    BOOST_CHECK_EQUAL(valueFromDevice, 42.0);

    auto uniquePtr = Opm::gpuistl::make_gpu_unique_ptr<double>();
    Opm::gpuistl::copyToGPU(128.0, uniquePtr);
    auto valueFromDeviceUnique = Opm::gpuistl::copyFromGPU(uniquePtr);
    BOOST_CHECK_EQUAL(valueFromDeviceUnique, 128.0);

    auto sharedPtrArray = Opm::gpuistl::make_gpu_shared_ptr<std::array<double, 4>>();
    Opm::gpuistl::copyToGPU(std::array<double, 4> {1.0, 2.0, 3.0, 4.0}, sharedPtrArray);
    auto valueFromDeviceArray = Opm::gpuistl::copyFromGPU(sharedPtrArray);
    BOOST_CHECK_EQUAL(valueFromDeviceArray[0], 1.0);
    BOOST_CHECK_EQUAL(valueFromDeviceArray[1], 2.0);
    BOOST_CHECK_EQUAL(valueFromDeviceArray[2], 3.0);
    BOOST_CHECK_EQUAL(valueFromDeviceArray[3], 4.0);
}
