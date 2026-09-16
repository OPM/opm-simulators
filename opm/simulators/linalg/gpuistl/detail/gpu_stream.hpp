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

#include <cstddef>
#include <sstream>
#include <stdexcept>
#include <string_view>

namespace Opm::gpuistl::detail
{

/**
 * @brief Debug-only check that @p stream is a valid CUDA stream handle.
 * @note Uses cudaStreamQuery; cudaSuccess and cudaErrorNotReady both mean the handle is valid.
 * @note Prefer OPM_GPUISTL_DETAIL_ASSERT_CUDA_STREAM so the stream name and call site are captured.
 *
 * @todo Refactor to use std::source_location once we shift to C++20
 */
inline void
assertCudaStream([[maybe_unused]] cudaStream_t stream,
                 [[maybe_unused]] std::string_view name,
                 [[maybe_unused]] std::string_view filename,
                 [[maybe_unused]] std::string_view functionName,
                 [[maybe_unused]] std::size_t lineNumber)
{
#ifndef NDEBUG
    const cudaError_t err = cudaStreamQuery(stream);
    if (err == cudaSuccess || err == cudaErrorNotReady) {
        return;
    }

    std::ostringstream str;
    str << name << " is not a valid CUDA stream\n"
        << "  file: " << filename << '\n'
        << "  line: " << lineNumber << '\n'
        << "  function: " << functionName << '\n'
        << "  cudaError: " << cudaGetErrorString(err);
    if (err == cudaErrorInvalidResourceHandle) {
        str << " (invalid stream handle)";
    }
    OPM_THROW(std::invalid_argument, str.str());
#endif
}

} // namespace Opm::gpuistl::detail

/**
 * @brief Captures the argument name and call site for assertCudaStream.
 */
#define OPM_GPUISTL_DETAIL_ASSERT_CUDA_STREAM(x)                                                   \
    ::Opm::gpuistl::detail::assertCudaStream((x), #x, __FILE__, __func__, __LINE__)

#endif // OPM_GPUISTL_DETAIL_GPU_STREAM_HPP
