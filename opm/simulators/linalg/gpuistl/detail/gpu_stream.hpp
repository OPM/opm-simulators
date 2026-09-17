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

#include <source_location>
#include <sstream>
#include <stdexcept>
#include <string_view>

namespace Opm::gpuistl::detail
{

/**
 * @brief Debug-only check that @p stream is a valid CUDA stream handle.
 * @note Uses cudaStreamQuery; cudaSuccess and cudaErrorNotReady both mean the handle is valid.
 * @note Prefer OPM_GPUISTL_DETAIL_ASSERT_CUDA_STREAM so the stream name is captured.
 *       Call-site file/line/function come from @p location (default: current()).
 */
inline void
assertCudaStream([[maybe_unused]] cudaStream_t stream,
                 [[maybe_unused]] std::string_view name,
                 [[maybe_unused]] const std::source_location location
                 = std::source_location::current())
{
#ifndef NDEBUG
    const cudaError_t err = cudaStreamQuery(stream);
    if (err == cudaSuccess || err == cudaErrorNotReady) {
        return;
    }

    std::ostringstream str;
    str << name << " is not a valid CUDA stream\n"
        << "  file: " << location.file_name() << '\n'
        << "  line: " << location.line() << '\n'
        << "  column: " << location.column() << '\n'
        << "  function: " << location.function_name() << '\n'
        << "  cudaError: " << cudaGetErrorString(err);
    if (err == cudaErrorInvalidResourceHandle) {
        str << " (invalid stream handle)";
    }
    OPM_THROW(std::invalid_argument, str.str());
#endif
}

} // namespace Opm::gpuistl::detail

/**
 * @brief Captures the argument name for assertCudaStream (call site via source_location).
 */
#define OPM_GPUISTL_DETAIL_ASSERT_CUDA_STREAM(x) ::Opm::gpuistl::detail::assertCudaStream((x), #x)

#endif // OPM_GPUISTL_DETAIL_GPU_STREAM_HPP
