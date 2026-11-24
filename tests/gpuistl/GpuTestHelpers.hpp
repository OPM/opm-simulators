/*
  Copyright 2026 Equinor ASA

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

/*!
 * \file
 *
 * \brief Helpers shared by the gpuistl unit tests.
 *
 * Contains kernels, so it can only be included from .cu files. It uses Boost.Test, so it
 * must be included after BOOST_TEST_MODULE is defined and <boost/test/unit_test.hpp> is
 * included; otherwise Boost.Test does not generate the test module's main().
 */
#ifndef OPM_GPUISTL_TESTS_GPU_TEST_HELPERS_HPP
#define OPM_GPUISTL_TESTS_GPU_TEST_HELPERS_HPP

#include <boost/test/unit_test.hpp>

#include <opm/common/utility/gpuDecorators.hpp>
#include <opm/material/densead/Evaluation.hpp>

#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpuThreadUtils.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>
#include <opm/simulators/linalg/gpuistl/detail/safe_conversion.hpp>
#include <opm/simulators/linalg/gpuistl/gpu_smart_pointer.hpp>

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <source_location>
#include <string>
#include <type_traits>
#include <vector>

namespace Opm::gpuistl::test
{

// ---------------------------------------------------------------------------------------------
// Kernel launches
// ---------------------------------------------------------------------------------------------

/*!
 * \brief Throws if the last kernel launch failed, or if the kernel failed while running.
 *
 * Call directly after a kernel launch. The error message names \p location, by default the
 * location of the call. A launch error is cleared, so that it is not reported again by the
 * next test case. An error from the running kernel cannot be cleared: it corrupts the CUDA
 * context, so the following test cases in the same executable fail as well.
 *
 * This is a function rather than OPM_CUDA_CHECK_LAST_ERROR followed by
 * OPM_CUDA_CHECK_DEVICE_SYNCHRONIZE (detail/cuda_check_last_error.hpp): a macro reports the
 * location where it is expanded, which inside a helper is the helper, while
 * std::source_location lets the helpers here pass on the location of the test.
 */
inline void
checkLaunchAndSynchronize(const std::source_location location = std::source_location::current())
{
    detail::cudaSafeCall(cudaGetLastError(),
                         "cudaGetLastError()",
                         location.file_name(),
                         location.function_name(),
                         location.line());
    detail::cudaSafeCall(cudaDeviceSynchronize(),
                         "cudaDeviceSynchronize()",
                         location.file_name(),
                         location.function_name(),
                         location.line());
}

constexpr unsigned defaultThreadsPerBlock = 128;

//! Number of blocks of \p threadsPerBlock threads needed for \p numThreads threads.
inline int
numBlocks(std::size_t numThreads, unsigned threadsPerBlock = defaultThreadsPerBlock)
{
    return detail::getNumberOfBlocks(detail::to_int(numThreads), static_cast<int>(threadsPerBlock));
}

/*!
 * \brief Launches \p kernel with one thread and returns the value it writes.
 *
 * \p kernel is called as kernel(args..., result), where result is a \c Result* in device
 * memory. The device memory is freed also when the launch fails. A launch error is reported
 * with the location in this function, not the caller's, since a defaulted
 * std::source_location parameter cannot follow \p args; the failing test case is still
 * named by Boost.Test.
 */
template <class Result = double, class Kernel, class... Args>
Result
runSingleThreadKernel(Kernel kernel, Args... args)
{
    auto result = make_gpu_unique_ptr<Result>();
    kernel<<<1, 1>>>(args..., result.get());
    checkLaunchAndSynchronize();
    return copyFromGPU(result);
}

// ---------------------------------------------------------------------------------------------
// Pointwise evaluation on the CPU and the GPU
//
// A functor f(object, point) evaluates the object (a PVT object, a material law, ...) at one
// point. The same functor is used for the CPU reference and in the kernel, so it must be
// OPM_HOST_DEVICE. Point and the result of f must be trivially copyable.
// ---------------------------------------------------------------------------------------------

template <class Functor, class Object, class Point, class Result>
__global__ void
evaluateKernel(Functor f, Object object, const Point* points, Result* results, std::size_t n)
{
    const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        results[i] = f(object, points[i]);
    }
}

//! As evaluateKernel(), but with the object in device memory, accessed through a pointer to
//! const (as when it is a member of another object in device memory).
template <class Functor, class Object, class Point, class Result>
__global__ void
evaluateThroughConstPointerKernel(
    Functor f, const Object* object, const Point* points, Result* results, std::size_t n)
{
    const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        results[i] = f(*object, points[i]);
    }
}

//! f(object, point) for each point, on the CPU.
template <class Functor, class Object, class Point>
auto
evaluateOnCpu(Functor f, const Object& object, const std::vector<Point>& points)
{
    using Result = decltype(f(object, points.front()));
    std::vector<Result> results;
    results.reserve(points.size());
    for (const auto& pt : points) {
        results.push_back(f(object, pt));
    }
    return results;
}

//! f(object, point) for each point, on the GPU. \p object is passed to the kernel by value,
//! so it must be trivially copyable (typically a view). A launch error is reported with
//! \p location, by default the location of the call.
template <class Functor, class Object, class Point>
auto
evaluateOnGpu(Functor f,
              Object object,
              const std::vector<Point>& points,
              const std::source_location location = std::source_location::current())
{
    using Result = decltype(f(object, points.front()));
    const GpuBuffer<Point> gpuPoints(points);
    GpuBuffer<Result> gpuResults(points.size());

    evaluateKernel<<<numBlocks(points.size()), defaultThreadsPerBlock>>>(
        f, object, gpuPoints.data(), gpuResults.data(), points.size());
    checkLaunchAndSynchronize(location);

    return gpuResults.asStdVector();
}

//! As evaluateOnGpu(), but with \p object copied to device memory and used in the kernel
//! through a pointer to const.
template <class Functor, class Object, class Point>
auto
evaluateOnGpuThroughConstPointer(Functor f,
                                 Object object,
                                 const std::vector<Point>& points,
                                 const std::source_location location
                                 = std::source_location::current())
{
    using Result = decltype(f(object, points.front()));
    const GpuBuffer<Object> gpuObject(std::vector<Object> {object});
    const GpuBuffer<Point> gpuPoints(points);
    GpuBuffer<Result> gpuResults(points.size());

    evaluateThroughConstPointerKernel<<<numBlocks(points.size()), defaultThreadsPerBlock>>>(
        f, gpuObject.data(), gpuPoints.data(), gpuResults.data(), points.size());
    checkLaunchAndSynchronize(location);

    return gpuResults.asStdVector();
}

// ---------------------------------------------------------------------------------------------
// Automatic differentiation
// ---------------------------------------------------------------------------------------------

//! \p value as an \c Eval: the plain value if \c Eval is a scalar type, otherwise an AD
//! variable with index \p varIdx.
template <class Eval, class Scalar>
OPM_HOST_DEVICE Eval
makeVariable(Scalar value, int varIdx)
{
    if constexpr (std::is_same_v<Eval, Scalar>) {
        return value;
    } else {
        return Eval::createVariable(value, varIdx);
    }
}

// ---------------------------------------------------------------------------------------------
// Copy semantics of GPU-ported classes
// ---------------------------------------------------------------------------------------------

/*!
 * \brief Checks the copy semantics of a class with a CPU, an owning GPU and a GPU view variant.
 *
 * The views in the owning GPU object point into device memory it owns, so a copy would have
 * views into the memory of the original: the owning GPU object must be move-only. The CPU
 * object must stay copyable, and the view is passed to kernels by value, so it must be
 * trivially copyable. Use as static_assert(hasGpuPortCopySemantics<Cpu, GpuBuf, GpuView>()).
 */
template <class Cpu, class GpuBuf, class GpuView>
constexpr bool
hasGpuPortCopySemantics()
{
    static_assert(!std::is_copy_constructible_v<GpuBuf>, "owning GPU object must be move-only");
    static_assert(!std::is_copy_assignable_v<GpuBuf>, "owning GPU object must be move-only");
    static_assert(std::is_move_constructible_v<GpuBuf>);
    static_assert(std::is_move_assignable_v<GpuBuf>);
    static_assert(std::is_copy_constructible_v<Cpu>);
    static_assert(std::is_trivially_copyable_v<GpuView>, "view is passed to kernels by value");
    return true;
}

// ---------------------------------------------------------------------------------------------
// Comparing CPU and GPU results
// ---------------------------------------------------------------------------------------------

//! Whether \p a and \p b agree to the relative tolerance \p relTol.
template <std::floating_point Scalar>
bool
isRelClose(Scalar a, Scalar b, std::type_identity_t<Scalar> relTol)
{
    if (a == b) {
        return true;
    }
    return std::abs(a - b) <= relTol * std::max(std::abs(a), std::abs(b));
}

/*!
 * \brief Whether the values and all derivatives of \p a and \p b agree to the relative
 *        tolerance \p relTol.
 *
 * A derivative that is zero up to rounding would be compared on rounding noise. For such
 * derivatives, \p derivativeFloor sets the scale the tolerance is relative to.
 */
template <class Scalar, int numDerivs, unsigned staticSize>
bool
isRelClose(const DenseAd::Evaluation<Scalar, numDerivs, staticSize>& a,
           const DenseAd::Evaluation<Scalar, numDerivs, staticSize>& b,
           std::type_identity_t<Scalar> relTol,
           std::type_identity_t<Scalar> derivativeFloor = 0)
{
    if (a.size() != b.size()) {
        return false;
    }
    bool close = isRelClose(a.value(), b.value(), relTol);
    for (int varIdx = 0; varIdx < a.size(); ++varIdx) {
        const Scalar da = a.derivative(varIdx);
        const Scalar db = b.derivative(varIdx);
        close = close
            && (da == db
                || std::abs(da - db)
                    <= relTol * std::max({std::abs(da), std::abs(db), derivativeFloor}));
    }
    return close;
}

/*!
 * \brief Checks that each GPU result is close to its CPU result.
 *
 * The comparison is passed in as \p isClose(cpu, gpu), so that it can be built from
 * overloads in the test's own namespace. \p describePoint(i) returns something printable
 * that identifies point i in the error message.
 */
template <class Result, class IsClose, class DescribePoint>
void
checkResults(const std::vector<Result>& cpu,
             const std::vector<Result>& gpu,
             IsClose isClose,
             DescribePoint describePoint,
             const std::string& name)
{
    BOOST_REQUIRE_EQUAL(cpu.size(), gpu.size());
    for (std::size_t i = 0; i < cpu.size(); ++i) {
        BOOST_CHECK_MESSAGE(isClose(cpu[i], gpu[i]),
                            name << " at point " << i << " (" << describePoint(i)
                                 << "): CPU = " << cpu[i] << ", GPU = " << gpu[i]);
    }
}

} // namespace Opm::gpuistl::test

#endif // OPM_GPUISTL_TESTS_GPU_TEST_HELPERS_HPP
