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

#define BOOST_TEST_MODULE TestDeviceBlockOperations

#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>
#include <opm/simulators/linalg/gpuistl/detail/deviceBlockOperations.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>

#include <cuda_runtime.h>

#include <boost/test/unit_test.hpp>

// Simple GPU kernel to test transposeBlock function
__global__ void
testTransposeKernel(const double* input, double* output)
{
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        transposeBlock<double, 3>(input, output);
    }
}

// Simple GPU kernel to test solveBlock function
__global__ void
testSolveKernel(const double* matrix, const double* rhs, double* solution)
{
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        solveBlock<double, 3>(matrix, rhs, solution);
    }
}

BOOST_AUTO_TEST_CASE(TestTransposeBlock)
{
    std::vector<double> hostMatrix = {
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
        7.0, 8.0, 9.0
    };

    std::vector<double> expectedTranspose = {
        1.0, 4.0, 7.0,
        2.0, 5.0, 8.0,
        3.0, 6.0, 9.0
    };

    Opm::gpuistl::GpuBuffer<double> gpuInput(hostMatrix);
    Opm::gpuistl::GpuBuffer<double> gpuOutput(9);

    testTransposeKernel<<<1, 1>>>(gpuInput.data(), gpuOutput.data());
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
    OPM_GPU_SAFE_CALL(cudaGetLastError());

    std::vector<double> hostResult = gpuOutput.asStdVector();

    BOOST_REQUIRE_EQUAL(hostResult.size(), expectedTranspose.size());
    for (size_t i = 0; i < hostResult.size(); ++i) {
        BOOST_CHECK_CLOSE(hostResult[i], expectedTranspose[i], 1e-12);
    }
}

BOOST_AUTO_TEST_CASE(TestSolveBlock)
{
    std::vector<double> hostMatrix = {
        4.0, 1.0, 2.0,
        1.0, 5.0, 1.0,
        2.0, 1.0, 6.0
    };

    std::vector<double> hostRhs = {1.0, 0.0, 1.0};

    Opm::gpuistl::GpuBuffer<double> gpuMatrix(hostMatrix);
    Opm::gpuistl::GpuVector<double> gpuRhs(hostRhs);
    Opm::gpuistl::GpuVector<double> gpuSolution(3);

    testSolveKernel<<<1, 1>>>(gpuMatrix.data(), gpuRhs.data(), gpuSolution.data());
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
    OPM_GPU_SAFE_CALL(cudaGetLastError());

    std::vector<double> hostSolution = gpuSolution.asStdVector();
    std::vector<double> expectedSolution = {0.2127659574468085, -0.06382978723404256, 0.10638297872340424};
    // Verify that the solution matches the expected solution
    BOOST_REQUIRE_EQUAL(hostSolution.size(), expectedSolution.size());
    for (size_t i = 0; i < hostSolution.size(); ++i) {
        BOOST_CHECK_CLOSE(hostSolution[i], expectedSolution[i], 1e-10);
    }

    // Also verify by substituting back into original equation: A * x = rhs
    std::vector<double> lhs(3);
    lhs[0] = hostMatrix[0] * hostSolution[0] + hostMatrix[1] * hostSolution[1] + hostMatrix[2] * hostSolution[2];
    lhs[1] = hostMatrix[3] * hostSolution[0] + hostMatrix[4] * hostSolution[1] + hostMatrix[5] * hostSolution[2];
    lhs[2] = hostMatrix[6] * hostSolution[0] + hostMatrix[7] * hostSolution[1] + hostMatrix[8] * hostSolution[2];

    for (size_t i = 0; i < 3; ++i) {
        BOOST_CHECK_SMALL(std::abs(lhs[i] - hostRhs[i]), 1e-14);
    }
}
