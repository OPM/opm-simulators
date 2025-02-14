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
#include <boost/test/tools/old/interface.hpp>
#include <config.h>
#include <stdexcept>

#define BOOST_TEST_MODULE TestThrowMacrosOnGPU

#include <cuda.h>
#include <cuda_runtime.h>
#include <boost/test/unit_test.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>

namespace {
__global__ void codeThatContainsMacros(bool call) {
    if (call) {
        OPM_THROW(std::logic_error, "Something went wrong");
        OPM_THROW_NOLOG(std::logic_error, "Something went wrong");
        OPM_THROW_PROBLEM(std::logic_error, "Something went wrong");
    }
    OPM_ERROR_IF(!call, "Something went horribly wrong");
}

// TODO: Check if this is better on HIP
#if 0 // I am leaving this here to show that this is not possible due to limitations in CUDA
      // the assert will indeed cause an error, but the CUDA context will be broken for
      // the rest of the lifetime of the process, see 
      // https://forums.developer.nvidia.com/t/how-to-clear-cuda-errors/296393/5
__global__ void checkThrow() {
    OPM_THROW(std::logic_error, "Something went wrong");
}

__global__ void checkThrowNoLog() {
    OPM_THROW_NOLOG(std::logic_error, "Something went wrong");
}

__global__ void checkThrowProblem() {
    OPM_THROW_PROBLEM(std::logic_error, "Something went wrong");
}

__global__ void checkErrorIf() {
    OPM_ERROR_IF(true, "Something went horribly wrong");
}
#endif
}

BOOST_AUTO_TEST_CASE(TestKernel)
{
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
    OPM_GPU_SAFE_CALL(cudaGetLastError());
    codeThatContainsMacros<<<1, 1>>>(false);
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
    OPM_GPU_SAFE_CALL(cudaGetLastError());

    #if 0 // I am leaving this here to show that this is not possible due to limitations in CUDA
          // the assert will indeed cause an error, but the CUDA context will be broken for
          // the rest of the lifetime of the process, see 
          // https://forums.developer.nvidia.com/t/how-to-clear-cuda-errors/296393/5
    codeThatContainsMacros<<<1, 1>>>(true);
    // Make sure this actually throws
    BOOST_CHECK_THROW(OPM_GPU_SAFE_CALL(cudaDeviceSynchronize()), std::runtime_error);
    OPM_GPU_SAFE_CALL(cudaDeviceReset());
    OPM_GPU_SAFE_CALL(cudaGetLastError());

    checkThrow<<<1, 1>>>();
    BOOST_CHECK_THROW(OPM_GPU_SAFE_CALL(cudaDeviceSynchronize()), std::runtime_error);
    OPM_GPU_SAFE_CALL(cudaDeviceReset());
    OPM_GPU_SAFE_CALL(cudaGetLastError());

    checkThrowNoLog<<<1, 1>>>();
    BOOST_CHECK_THROW(OPM_GPU_SAFE_CALL(cudaDeviceSynchronize()), std::runtime_error);
    OPM_GPU_SAFE_CALL(cudaDeviceReset());
    OPM_GPU_SAFE_CALL(cudaGetLastError());

    checkThrowProblem<<<1, 1>>>();
    BOOST_CHECK_THROW(OPM_GPU_SAFE_CALL(cudaDeviceSynchronize()), std::runtime_error);
    OPM_GPU_SAFE_CALL(cudaDeviceReset());
    OPM_GPU_SAFE_CALL(cudaGetLastError());

    checkErrorIf<<<1, 1>>>();
    BOOST_CHECK_THROW(OPM_GPU_SAFE_CALL(cudaDeviceSynchronize()), std::runtime_error);
    OPM_GPU_SAFE_CALL(cudaDeviceReset());
    OPM_GPU_SAFE_CALL(cudaGetLastError());
    #endif
}

BOOST_AUTO_TEST_CASE(TestOutsideKernel) 
{
    // This is to make sure that the macros work outside of kernels but inside a .cu file
    // ie. inside a file compiled by nvcc/hipcc.
    BOOST_CHECK_THROW(OPM_THROW(std::runtime_error, "THROW"), std::runtime_error);
    BOOST_CHECK_THROW(OPM_THROW_NOLOG(std::runtime_error, "THROW_NOLOG"), std::runtime_error);
    BOOST_CHECK_THROW(OPM_THROW_PROBLEM(std::runtime_error, "THROW_PROBLEM"), std::runtime_error);
    BOOST_CHECK_THROW(OPM_ERROR_IF(true, "ERROR_IF"), std::logic_error);
}
