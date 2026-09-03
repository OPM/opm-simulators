/*
  Copyright 2025 Equinor ASA
  Copyright 2026 SINTEF AS
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

#ifndef OPM_SIMULATORS_LINALG_GPUISTL_DETAIL_GPU_POINTER_ATTRIBUTES_HPP
#define OPM_SIMULATORS_LINALG_GPUISTL_DETAIL_GPU_POINTER_ATTRIBUTES_HPP

#include <cuda.h>
#include <cuda_runtime.h>
#include <memory>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>

namespace Opm::gpuistl::detail
{

/**
 * @brief Checks whether the given pointer is associated with GPU device memory.
 *
 * This function retrieves CUDA pointer attributes for the provided pointer and
 * determines whether it references device memory. It returns true if the pointer
 * corresponds to GPU memory; otherwise, it returns false.
 *
 * @tparam T Type of elements pointed to by the input pointer.
 * @param ptr Pointer to the memory that needs to be checked.
 * @return True if the pointer represents GPU memory, false otherwise.
 */
template <class T>
inline bool
isGPUPointer(const T* ptr)
{
    if (ptr == nullptr) {
        return false;
    }
    cudaPointerAttributes attributes;
    OPM_GPU_SAFE_CALL(cudaPointerGetAttributes(&attributes, ptr));
    return attributes.type == cudaMemoryTypeDevice;
}


/**
 * @brief Checks if the given std::unique_ptr with custom deleter refers to GPU memory.
 *
 * @tparam T The type stored within the pointer.
 * @tparam D The custom deleter type.
 * @param ptr The std::unique_ptr object to inspect.
 * @return true if the pointer addresses GPU memory; otherwise false.
 */
template <class T, class D>
inline bool
isGPUPointer(const std::unique_ptr<T, D>& ptr)
{
    return isGPUPointer(ptr.get());
}

/**
 * @brief Checks if the given std::shared_ptr refers to GPU memory.
 *
 * @tparam T The type stored within the pointer.
 * @param ptr The std::shared_ptr object to inspect.
 * @return true if the pointer addresses GPU memory; otherwise false.
 */
template <class T>
inline bool
isGPUPointer(const std::shared_ptr<T>& ptr)
{
    return isGPUPointer(ptr.get());
}


// Check if the pointer is a CPU pointer -----------------------------------------

/**
 * @brief Checks whether the given pointer is associated with CPU device memory.
 *
 * This function retrieves CUDA pointer attributes for the provided pointer and
 * determines whether it references device memory. It returns true if the pointer
 * corresponds to CPU memory; otherwise, it returns false.
 *
 * @note This does not prove that @p ptr is a valid CPU allocation. CUDA reports
 *       cudaMemoryTypeUnregistered for any address it does not know about,
 *       such as stack allocations, malloc, etc, but also including out-of-range
 *       offsets into a device allocation. But for most use cases, this is sufficient.
 *
 * @tparam T Type of elements pointed to by the input pointer.
 * @param ptr Pointer to the memory that needs to be checked.
 * @return True if the pointer represents CPU memory, false otherwise.
 */
template <class T>
inline bool
isCPUPointer(const T* ptr)
{
    if (ptr == nullptr) {
        return false;
    }
    cudaPointerAttributes attributes;
    OPM_GPU_SAFE_CALL(cudaPointerGetAttributes(&attributes, ptr));
    return attributes.type == cudaMemoryTypeHost || attributes.type == cudaMemoryTypeUnregistered;
}


/**
 * @brief Checks if the given std::unique_ptr with custom deleter refers to CPU memory.
 *
 * @tparam T The type stored within the pointer.
 * @tparam D The custom deleter type.
 * @param ptr The std::unique_ptr object to inspect.
 * @return true if the pointer addresses CPU memory; otherwise false.
 */
template <class T, class D>
inline bool
isCPUPointer(const std::unique_ptr<T, D>& ptr)
{
    return isCPUPointer(ptr.get());
}

/**
 * @brief Checks if the given std::shared_ptr refers to CPU memory.
 *
 * @tparam T The type stored within the pointer.
 * @param ptr The std::shared_ptr object to inspect.
 * @return true if the pointer addresses CPU memory; otherwise false.
 */
template <class T>
inline bool
isCPUPointer(const std::shared_ptr<T>& ptr)
{
    return isCPUPointer(ptr.get());
}
} // namespace Opm::gpuistl::detail
#endif
