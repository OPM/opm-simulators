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
#include <csignal>
#include <cstdlib>
#include <stdexcept>

#include <sys/wait.h>
#include <unistd.h>

#define BOOST_TEST_MODULE TestThrowMacrosOnGPU

#include <boost/test/unit_test.hpp>
#include <cuda.h>
#include <cuda_runtime.h>
#include <opm/common/ErrorMacros.hpp>

namespace {
// NOTE: We have to split this into a separate function due
// to some weirdness of hipcc. Note however that this is
// the realistic use case of the macro.
__device__ __host__ void functionThatContainsMacros(bool call) {
    if (call) {
        OPM_THROW(std::logic_error, "Something went wrong");
        OPM_THROW_NOLOG(std::logic_error, "Something went wrong");
        OPM_THROW_PROBLEM(std::logic_error, "Something went wrong");
    }
    OPM_ERROR_IF(!call, "Something went horribly wrong");
}
__global__ void codeThatContainsMacros(bool call) {
    functionThatContainsMacros(call);
}
}

BOOST_AUTO_TEST_CASE(TestKernel)
{
    // A device-side abort leaves the GPU context unusable, and ROCm also aborts
    // the process. Run the kernel in a child so the remaining tests can proceed.
    constexpr auto macroFailureExitCode = EXIT_FAILURE;
    constexpr auto testSetupFailureExitCode = 2;
    const auto childProcessId = fork();
    BOOST_REQUIRE_NE(childProcessId, -1);

    if (childProcessId == 0) {
        std::signal(SIGABRT, SIG_DFL);

        if (cudaDeviceSynchronize() != cudaSuccess) {
            _exit(testSetupFailureExitCode);
        }

        codeThatContainsMacros<<<1, 1>>>(false);
        if (cudaGetLastError() != cudaSuccess) {
            _exit(testSetupFailureExitCode);
        }

        const auto synchronizeResult = cudaDeviceSynchronize();
        _exit(synchronizeResult == cudaSuccess ? testSetupFailureExitCode : macroFailureExitCode);
    }

    auto childStatus = 0;
    BOOST_REQUIRE_EQUAL(waitpid(childProcessId, &childStatus, 0), childProcessId);

    const auto childAborted = WIFSIGNALED(childStatus) && (WTERMSIG(childStatus) == SIGABRT);
    const auto reportedKernelFailure
        = WIFEXITED(childStatus) && (WEXITSTATUS(childStatus) == macroFailureExitCode);
    BOOST_CHECK(childAborted || reportedKernelFailure);
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
