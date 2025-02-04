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

#ifndef OPM_SIMULATORS_LINALG_GPUISTL_GPU_SMART_POINTER_HPP
#define OPM_SIMULATORS_LINALG_GPUISTL_GPU_SMART_POINTER_HPP

#include <cuda_runtime.h>

#include <memory>

#include <opm/common/utility/gpuDecorators.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>

/**
 * @file gpu_smart_pointer.hpp defines convenience classes and functions for using std::shared_ptr and std::unique_ptr
 * with GPU allocated memory.
 */

namespace Opm::gpuistl
{


/**
 * @brief Creates a shared pointer managing GPU-allocated memory of the specified element type.
 *
 * This function allocates memory on the GPU for the type \c T, using \c cudaMalloc.
 * It returns a \c std::shared_ptr that automatically handles the release of
 * GPU memory with cudaFree when no longer in use.
 *
 * @tparam T The element type to allocate on the GPU.
 * @return A std::shared_ptr to the GPU-allocated memory.
 */
template <typename T>
std::shared_ptr<T>
make_gpu_shared_ptr()
{
    T* ptr = nullptr;
    OPM_GPU_SAFE_CALL(cudaMalloc(&ptr, sizeof(T)));
    auto deleter = [](T* ptr) { OPM_GPU_WARN_IF_ERROR(cudaFree(ptr)); };
    return std::shared_ptr<T>(ptr, deleter);
}

/**
 * @brief Creates a shared pointer managing GPU-allocated memory of the specified element type.
 *
 * This function allocates memory on the GPU for the type \c T, using \c cudaMalloc.
 * It returns a std::shared_ptr that automatically handles the release of
 * GPU memory with cudaFree when no longer in use.
 *
 * @tparam T The element type to allocate on the GPU.
 * @param value The value to copy to the GPU-allocated memory.
 * @return A std::shared_ptr to the GPU-allocated memory.
 */
template <typename T>
std::shared_ptr<T>
make_gpu_shared_ptr(const T& value)
{
    auto ptr = make_gpu_shared_ptr<T>();
    OPM_GPU_SAFE_CALL(cudaMemcpy(ptr.get(), &value, sizeof(T), cudaMemcpyHostToDevice));
    return ptr;
}


/**
 * @brief Creates a unique pointer managing GPU-allocated memory of the specified element type.
 *
 * This function allocates memory on the GPU for the type \c T, using \c cudaMalloc .
 * It returns a std::unique_ptr that automatically handles the release of
 * GPU memory with cudaFree when no longer in use.
 *
 * @tparam T The element type to allocate on the GPU.
 * @return A std::unique_ptr to the GPU-allocated memory.
 */
template <typename T>
auto
make_gpu_unique_ptr()
{
    T* ptr = nullptr;
    OPM_GPU_SAFE_CALL(cudaMalloc(&ptr, sizeof(T)));

    auto deleter = [](T* ptr) { OPM_GPU_WARN_IF_ERROR(cudaFree(ptr)); };
    return std::unique_ptr<T, decltype(deleter)>(ptr, deleter);
}

/**
 * @brief Creates a unique pointer managing GPU-allocated memory of the specified element type.
 *
 * This function allocates memory on the GPU for the type \c T, using \c cudaMalloc.
 * It returns a std::unique_ptr that automatically handles the release of
 * GPU memory with cudaFree when no longer in use.
 *
 * @tparam T The element type to allocate on the GPU.
 * @param value The value to copy to the GPU-allocated memory.
 * @return A std::unique_ptr to the GPU-allocated memory.
 */
template <typename T>
auto
make_gpu_unique_ptr(const T& value)
{
    auto ptr = make_gpu_unique_ptr<T>();
    OPM_GPU_SAFE_CALL(cudaMemcpy(ptr.get(), &value, sizeof(T), cudaMemcpyHostToDevice));
    return ptr;
}

/**
 * @brief A view towards a smart pointer to GPU-allocated memory.
 *
 * This will emulate a smart pointer to GPU-allocated memory, but without ownership semantics, and
 * being compatible with the requirements of the GPU kernels. This is useful when we want to pass
 * a smart pointer to a GPU kernel, but we do not want to transfer the ownership of the memory.
 */
template <class T>
class PointerView
{
public:
    PointerView(const PointerView& other) = default;

    PointerView(const std::shared_ptr<T>& ptr)
        : ptr_(ptr.get())
    {
    }

    template <class Deleter>
    PointerView(const std::unique_ptr<T, Deleter>& ptr)
        : ptr_(ptr.get())
    {
    }

    PointerView(T* ptr)
        : ptr_(ptr)
    {
    }

    OPM_HOST_DEVICE T* get() const
    {
        return ptr_;
    }

    OPM_HOST_DEVICE T& operator*() const
    {
        return *ptr_;
    }

    OPM_HOST_DEVICE T* operator->() const
    {
        return ptr_;
    }

private:
    T* ptr_;
};

template <class T>
PointerView<T>
make_view(const std::shared_ptr<T>& ptr)
{
    return PointerView<T>(ptr);
}

template <class T, class Deleter>
PointerView<T>
make_view(const std::unique_ptr<T, Deleter>& ptr)
{
    return PointerView<T>(ptr);
}
} // namespace Opm::gpuistl
#endif
