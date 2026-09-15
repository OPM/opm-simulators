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

#ifndef OPM_GPUISTL_DETAIL_GPU_MEMCPY_HPP
#define OPM_GPUISTL_DETAIL_GPU_MEMCPY_HPP

#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_pointer_attributes.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>

#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>

#include <cuda_runtime.h>

namespace Opm::gpuistl::detail {

/**
 * @brief Debug-only check that @p ptr is host/CPU memory.
 * @note Uses cudaPointerGetAttributes; intended for debug builds only.
 */
inline void
assertHostPointer([[maybe_unused]] const void* ptr, [[maybe_unused]] std::string_view name)
{
#ifndef NDEBUG
    if (!isCPUPointer(ptr)) {
        OPM_THROW(std::invalid_argument, std::string(name) + " is not a CPU pointer");
    }
#endif
}

/**
 * @brief Debug-only check that @p ptr is device/GPU memory.
 * @note Uses cudaPointerGetAttributes; intended for debug builds only.
 */
inline void
assertDevicePointer([[maybe_unused]] const void* ptr, [[maybe_unused]] std::string_view name)
{
#ifndef NDEBUG
    if (!isGPUPointer(ptr)) {
        OPM_THROW(std::invalid_argument, std::string(name) + " is not a device/GPU pointer");
    }
#endif
}


/**
 * @brief gpuMemcpyHostToDevice copies count elements of type T from host to device.
 * @param dstDevice raw pointer to GPU memory (destination)
 * @param srcHost raw pointer to CPU memory (source)
 * @param count number of elements to copy
 *
 * @tparam T element type
 *
 * @note This does synchronous transfer
 * @note In debug builds, checks that @p srcHost is a CPU pointer and @p dstDevice is a GPU pointer.
 */
template <typename T>
inline void
gpuMemcpyHostToDevice(T* dstDevice, const T* srcHost, std::size_t count)
{
    if (count == 0) {
        return;
    }
    assertDevicePointer(dstDevice, "dstDevice");
    assertHostPointer(srcHost, "srcHost");
    OPM_GPU_SAFE_CALL(cudaMemcpy(dstDevice, srcHost, count * sizeof(T), cudaMemcpyHostToDevice));
}

/**
 * @brief gpuMemcpyDeviceToHost copies count elements of type T from device to host.
 * @param dstHost raw pointer to CPU memory (destination)
 * @param srcDevice raw pointer to GPU memory (source)
 * @param count number of elements to copy
 *
 * @tparam T element type
 *
 * @note This does synchronous transfer
 * @note In debug builds, checks that @p dstHost is a CPU pointer and @p srcDevice is a GPU pointer.
 */
template <typename T>
inline void
gpuMemcpyDeviceToHost(T* dstHost, const T* srcDevice, std::size_t count)
{
    if (count == 0) {
        return;
    }
    assertHostPointer(dstHost, "dstHost");
    assertDevicePointer(srcDevice, "srcDevice");
    OPM_GPU_SAFE_CALL(cudaMemcpy(dstHost, srcDevice, count * sizeof(T), cudaMemcpyDeviceToHost));
}

/**
 * @brief gpuMemcpyDeviceToDevice copies count elements of type T from device to device.
 * @param dstDevice raw pointer to GPU memory (destination)
 * @param srcDevice raw pointer to GPU memory (source)
 * @param count number of elements to copy
 *
 * @tparam T element type
 *
 * @note This does synchronous transfer
 * @note In debug builds, checks that @p dstDevice and @p srcDevice are GPU pointers.
 */
template <typename T>
inline void
gpuMemcpyDeviceToDevice(T* dstDevice, const T* srcDevice, std::size_t count)
{
    if (count == 0) {
        return;
    }
    assertDevicePointer(dstDevice, "dstDevice");
    assertDevicePointer(srcDevice, "srcDevice");
    OPM_GPU_SAFE_CALL(cudaMemcpy(dstDevice, srcDevice, count * sizeof(T), cudaMemcpyDeviceToDevice));
}

/**
 * @brief gpuMemcpyHostToDeviceAsync copies count elements of type T from host to device asynchronously.
 * @param dstDevice raw pointer to GPU memory (destination)
 * @param srcHost raw pointer to CPU memory (source)
 * @param count number of elements to copy
 * @param stream CUDA stream to use for the asynchronous copy
 *
 * @tparam T element type
 *
 * @note This does asynchronous transfer. If the memory region pointed to by @p srcHost
 *       has been previously registered (e.g., using cudaHostRegister by an external mechanism
 *       like PinnedMemoryHolder), the transfer may be faster.
 * @note In debug builds, checks that @p srcHost is a CPU pointer and @p dstDevice is a GPU pointer.
 * @note Does not synchronize the stream; the caller is responsible for stream completion.
 * @note Expects caller to specify the @p stream (no default is provided).
 */
template <typename T>
inline void
gpuMemcpyHostToDeviceAsync(T* dstDevice, const T* srcHost, std::size_t count, cudaStream_t stream)
{
    if (count == 0) {
        return;
    }
    assertDevicePointer(dstDevice, "dstDevice");
    assertHostPointer(srcHost, "srcHost");
    OPM_GPU_SAFE_CALL(cudaMemcpyAsync(dstDevice, srcHost, count * sizeof(T), cudaMemcpyHostToDevice, stream));
}

/**
 * @brief gpuMemcpyDeviceToHostAsync copies count elements of type T from device to host asynchronously.
 * @param dstHost raw pointer to CPU memory (destination)
 * @param srcDevice raw pointer to GPU memory (source)
 * @param count number of elements to copy
 * @param stream CUDA stream to use for the asynchronous copy
 *
 * @tparam T element type
 *
 * @note This does asynchronous transfer. If the memory region pointed to by @p dstHost
 *       has been previously registered (e.g., using cudaHostRegister by an external mechanism
 *       like PinnedMemoryHolder), the transfer may be faster.
 * @note In debug builds, checks that @p dstHost is a CPU pointer and @p srcDevice is a GPU pointer.
 * @note Does not synchronize the stream; the caller is responsible for stream completion.
 * @note Expects caller to specify the @p stream (no default is provided).
 */
template <typename T>
inline void
gpuMemcpyDeviceToHostAsync(T* dstHost, const T* srcDevice, std::size_t count, cudaStream_t stream)
{
    if (count == 0) {
        return;
    }
    assertHostPointer(dstHost, "dstHost");
    assertDevicePointer(srcDevice, "srcDevice");
    OPM_GPU_SAFE_CALL(cudaMemcpyAsync(dstHost, srcDevice, count * sizeof(T), cudaMemcpyDeviceToHost, stream));
}

} // namespace Opm::gpuistl::detail

#endif // OPM_GPUISTL_DETAIL_GPU_MEMCPY_HPP
