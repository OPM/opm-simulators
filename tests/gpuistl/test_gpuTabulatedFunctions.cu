/*
  Copyright 2026 Equinor ASA

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
#include <config.h>

#define BOOST_TEST_MODULE TestGpuTabulatedFunctions

#include <boost/test/unit_test.hpp>

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/common/UniformTabulated2DFunction.hpp>
#include <opm/material/common/UniformXTabulated2DFunction.hpp>
#include <opm/material/common/UniformXTabulated2DFunctionBuilder.hpp>

#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>
#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>

#include <cuda_runtime.h>
#include <cmath>
#include <utility>
#include <vector>

/*
    This file contains unit tests exercising the UniformTabulated2DFunction and
    UniformXTabulated2DFunction classes on the GPU. Each test builds a CPU table,
    copies it to the GPU via copy_to_gpu()/make_view(), evaluates it on the device
    and checks that the result matches the CPU evaluation of the very same table.
*/

using Evaluation = Opm::DenseAd::Evaluation<double, 3>;

using GpuViewUniformTab = Opm::UniformTabulated2DFunction<double, Opm::gpuistl::GpuView>;

using XTabulatedFunction = Opm::UniformXTabulated2DFunction<double>;
using GpuViewXTab = Opm::UniformXTabulated2DFunction<double, Opm::gpuistl::GpuView>;
using XTabBuilder = Opm::UniformXTabulated2DFunctionBuilder<double>;
using InterpolationPolicy = XTabulatedFunction::InterpolationPolicy;

namespace {

const double ABS_TOL = 1e-6;

// Kernel to evaluate a UniformTabulated2DFunction on the GPU
__global__ void gpuEvaluateUniformTabulated2DFunction(GpuViewUniformTab gpuTab, Evaluation* inputX, Evaluation* inputY, double* result) {
    *result = gpuTab.eval(*inputX, *inputY, true).value();
}

// Kernel to evaluate a UniformXTabulated2DFunction on the GPU
__global__ void gpuEvaluateUniformXTabulated2DFunction(GpuViewXTab gpuTab, double x, double y, double* result) {
    *result = gpuTab.eval(x, y, /*extrapolate=*/false);
}

// Helper function to launch a kernel and retrieve the result on the CPU to reduce code duplication
template <typename KernelFunc, typename... Args>
double launchKernelAndRetrieveResult(KernelFunc kernel, Args... args) {
    double* resultOnGpu;
    double gpuComputedResultOnCpu;

    // Allocate memory for the result on the GPU
    OPM_GPU_SAFE_CALL(cudaMalloc(&resultOnGpu, sizeof(double)));

    // Launch the kernel
    kernel<<<1, 1>>>(args..., resultOnGpu);

    // Check for any errors in kernel launch
    OPM_GPU_SAFE_CALL(cudaPeekAtLastError());
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());

    // Retrieve the result from the GPU to the CPU
    OPM_GPU_SAFE_CALL(cudaMemcpy(&gpuComputedResultOnCpu, resultOnGpu, sizeof(double), cudaMemcpyDeviceToHost));

    // Free allocated GPU memory
    OPM_GPU_SAFE_CALL(cudaFree(resultOnGpu));

    return gpuComputedResultOnCpu;
}

// Regular grid: same X spacing and the same number/position of Y samples in every
// column. This is the direct GPU analogue of UniformTabulated2DFunction.
XTabulatedFunction buildRegularXTab(InterpolationPolicy policy = InterpolationPolicy::Vertical)
{
    using Scalar = double;

    const Scalar xMin = -2.0;
    const Scalar xMax = 3.0;
    const unsigned m = 5;

    const Scalar yMin = -0.5;
    const Scalar yMax = 1.0 / 3.0;
    const unsigned n = 4;

    XTabBuilder builder(policy);

    auto f = [](Scalar x, Scalar y) { return x * y; };

    for (unsigned i = 0; i < m; ++i) {
        const Scalar x = xMin + Scalar(i) / (m - 1) * (xMax - xMin);
        builder.appendXPos(x);
        for (unsigned j = 0; j < n; ++j) {
            const Scalar y = yMin + Scalar(j) / (n - 1) * (yMax - yMin);
            builder.appendSamplePoint(i, y, f(x, y));
        }
    }

    return std::move(builder).build();
}

// Ragged Y layout: every X column has a different number of Y sample points, so the
// SparseTable backing samples_ is jagged rather than rectangular. This is the layout
// that stresses the GPU port of the class (valueAt/ySegmentIndex/yToBeta indexing into
// rows of differing length on the device).
XTabulatedFunction buildRaggedXTab()
{
    using Scalar = double;

    const Scalar xMin = -2.0;
    const Scalar xMax = 3.0;
    const unsigned m = 6;

    const Scalar yMin = -4.0;
    const Scalar yMax = 5.0;

    XTabBuilder builder(InterpolationPolicy::Vertical);

    auto f = [](Scalar x, Scalar y) { return x * y; };

    for (unsigned i = 0; i < m; ++i) {
        const Scalar x = xMin + Scalar(i) / (m - 1) * (xMax - xMin);
        builder.appendXPos(x);

        const unsigned n = i + 4; // deliberately different number of Y points per column
        for (unsigned j = 0; j < n; ++j) {
            const Scalar y = yMin + Scalar(j) / (n - 1) * (yMax - yMin);
            builder.appendSamplePoint(i, y, f(x, y));
        }
    }

    return std::move(builder).build();
}

// Sloped guide curve: the upper Y bound grows with X (e.g. a bubble-point line), which
// gives interpolationGuide_ != Vertical a non-trivial, X-dependent shift to apply in
// findPoints(). With a flat (regular) guide the shift collapses to zero and the guided
// code path would not actually be exercised.
XTabulatedFunction buildSlopedGuideXTab(InterpolationPolicy policy)
{
    using Scalar = double;

    const Scalar xMin = 0.0;
    const Scalar xMax = 4.0;
    const unsigned m = 5;
    const unsigned n = 6;

    XTabBuilder builder(policy);

    auto f = [](Scalar x, Scalar y) { return x + y; };

    for (unsigned i = 0; i < m; ++i) {
        const Scalar x = xMin + Scalar(i) / (m - 1) * (xMax - xMin);
        builder.appendXPos(x);

        const Scalar yMin = 0.0;
        const Scalar yMax = 1.0 + 0.5 * x;
        for (unsigned j = 0; j < n; ++j) {
            const Scalar y = yMin + Scalar(j) / (n - 1) * (yMax - yMin);
            builder.appendSamplePoint(i, y, f(x, y));
        }
    }

    return std::move(builder).build();
}

void checkGpuMatchesCpu(const XTabulatedFunction& cpuTab, const std::vector<std::pair<double, double>>& points)
{
    auto gpuBufTab = Opm::gpuistl::copy_to_gpu(cpuTab);
    GpuViewXTab gpuViewTab = Opm::gpuistl::make_view(gpuBufTab);

    for (const auto& [x, y] : points) {
        const double cpuResult = cpuTab.eval(x, y, false);
        const double gpuResult = launchKernelAndRetrieveResult(gpuEvaluateUniformXTabulated2DFunction, gpuViewTab, x, y);
        BOOST_CHECK_MESSAGE(std::fabs(gpuResult - cpuResult) < ABS_TOL,
                            "eval(" << x << ", " << y << "): gpu=" << gpuResult << " cpu=" << cpuResult);
    }
}

} // END EMPTY NAMESPACE

// Test case for evaluating a UniformTabulated2DFunction on both CPU and GPU
BOOST_AUTO_TEST_CASE(TestEvaluateUniformTabulated2DFunctionOnGpu) {
    // Example tabulated data (2D)
    std::vector<std::vector<double>> tabData = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};

    // CPU-side function definition
    Opm::UniformTabulated2DFunction<double> cpuTab(1.0, 6.0, 3, 1.0, 6.0, 2, tabData);

    // Move data to GPU buffer and create a view for GPU operations
    Opm::UniformTabulated2DFunction<double, Opm::gpuistl::GpuBuffer> gpuBufTab = Opm::gpuistl::copy_to_gpu(cpuTab);
    GpuViewUniformTab gpuViewTab = Opm::gpuistl::make_view(gpuBufTab);

    // Evaluation points on the CPU
    Evaluation a(2.3);
    Evaluation b(4.5);

    // Allocate GPU memory for the Evaluation inputs
    Evaluation* gpuA = nullptr;
    Evaluation* gpuB = nullptr;
    OPM_GPU_SAFE_CALL(cudaMalloc(&gpuA, sizeof(Evaluation)));
    OPM_GPU_SAFE_CALL(cudaMemcpy(gpuA, &a, sizeof(Evaluation), cudaMemcpyHostToDevice));
    OPM_GPU_SAFE_CALL(cudaMalloc(&gpuB, sizeof(Evaluation)));
    OPM_GPU_SAFE_CALL(cudaMemcpy(gpuB, &b, sizeof(Evaluation), cudaMemcpyHostToDevice));

    const double gpuComputedResultOnCpu = launchKernelAndRetrieveResult(gpuEvaluateUniformTabulated2DFunction, gpuViewTab, gpuA, gpuB);

    // Free allocated GPU memory
    OPM_GPU_SAFE_CALL(cudaFree(gpuA));
    OPM_GPU_SAFE_CALL(cudaFree(gpuB));

    // Verify that the CPU and GPU results match within a reasonable tolerance
    const double cpuComputedResult = cpuTab.eval(a, b, true).value();
    BOOST_CHECK(std::fabs(gpuComputedResultOnCpu - cpuComputedResult) < ABS_TOL);
}

// Regular grid: baseline round trip for UniformXTabulated2DFunction, analogous to the
// UniformTabulated2DFunction test above.
BOOST_AUTO_TEST_CASE(TestEvaluateUniformXTabulated2DFunctionOnGpu) {
    const auto cpuTab = buildRegularXTab();

    const std::vector<std::pair<double, double>> points = {
        {-1.5, -0.4},
        {0.0, 0.0},
        {1.2, 0.2},
        {2.8, 0.3},
    };

    checkGpuMatchesCpu(cpuTab, points);
}

// Ragged Y layout: each X column has a different number of Y sample points, exercising
// the jagged SparseTable storage on the device. The chosen (x, y) points land in the
// first and last columns (fast-path segment lookup) as well as interior columns that
// force the bisection branch of xSegmentIndex/ySegmentIndex.
BOOST_AUTO_TEST_CASE(TestEvaluateUniformXTabulated2DFunctionRaggedYOnGpu) {
    const auto cpuTab = buildRaggedXTab();

    const std::vector<std::pair<double, double>> points = {
        {-1.9, -3.9}, // first column, fewest Y samples
        {-0.5, 0.0},  // interior column, forces X bisection
        {0.6, 2.0},   // interior column, deeper X bisection
        {2.9, 4.9},   // last column, most Y samples
    };

    checkGpuMatchesCpu(cpuTab, points);
}

// Non-default interpolation guide: a sloped guide curve gives RightExtreme a
// non-trivial, position-dependent shift to compute in findPoints(). yPos_ is a
// separate Storage<Scalar> member from samples_, so this checks that it also
// survived copy_to_gpu()/make_view() correctly.
BOOST_AUTO_TEST_CASE(TestEvaluateUniformXTabulated2DFunctionRightExtremeOnGpu) {
    const auto cpuTab = buildSlopedGuideXTab(InterpolationPolicy::RightExtreme);

    const std::vector<std::pair<double, double>> points = {
        {0.4, 0.3},
        {1.3, 0.7},
        {2.5, 1.0},
        {3.6, 1.5},
    };

    checkGpuMatchesCpu(cpuTab, points);
}
