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

#ifndef OPM_GPUISTL_DETAIL_GPU_STREAM_HPP
#define OPM_GPUISTL_DETAIL_GPU_STREAM_HPP

#include <opm/common/ErrorMacros.hpp>

#include <cuda_runtime.h>

#include <format>
#include <source_location>
#include <stdexcept>
#include <string_view>

namespace Opm::gpuistl::detail
{

/**
 * @brief Debug-only check that @p stream is a valid GPU stream handle.
 * @note Uses cudaStreamQuery; cudaSuccess and cudaErrorNotReady both mean the handle is valid.
 * @note Prefer OPM_GPUISTL_DETAIL_ASSERT_GPU_STREAM so the stream name is captured.
 *       Call-site file/line/function come from @p location (default: current()).
 */
inline void
assertGPUStream([[maybe_unused]] cudaStream_t stream,
                 [[maybe_unused]] std::string_view name,
                 [[maybe_unused]] const std::source_location location
                 = std::source_location::current())
{
#ifndef NDEBUG
    const cudaError_t err = cudaStreamQuery(stream);
    if (err == cudaSuccess || err == cudaErrorNotReady) {
        return;
    }

    if (err != cudaErrorInvalidResourceHandle) {
        OPM_THROW(std::runtime_error,
                  std::format("GPU stream query for {} failed\n"
                              "  file: {}\n"
                              "  line: {}\n"
                              "  column: {}\n"
                              "  function: {}\n"
                              "  GPU stream error: {}",
                              name,
                              location.file_name(),
                              location.line(),
                              location.column(),
                              location.function_name(),
                              cudaGetErrorString(err)));
    }

    OPM_THROW(std::invalid_argument,
              std::format("{} is not a valid GPU stream\n"
                          "  file: {}\n"
                          "  line: {}\n"
                          "  column: {}\n"
                          "  function: {}\n"
                          "  GPU stream error: {} (invalid stream handle)",
                          name,
                          location.file_name(),
                          location.line(),
                          location.column(),
                          location.function_name(),
                          cudaGetErrorString(err)));
#endif
}

} // namespace Opm::gpuistl::detail

/**
 * @brief Captures the argument name for assertCudaStream (call site via source_location).
 */
#define OPM_GPUISTL_DETAIL_ASSERT_GPU_STREAM(x) ::Opm::gpuistl::detail::assertGPUStream((x), #x)

#endif // OPM_GPUISTL_DETAIL_GPU_STREAM_HPP
