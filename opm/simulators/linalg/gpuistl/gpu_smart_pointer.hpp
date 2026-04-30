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
#include <opm/simulators/linalg/gpuistl/detail/is_gpu_pointer.hpp>

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
    auto deleter = [](T* ptrToDelete) { OPM_GPU_WARN_IF_ERROR(cudaFree(ptrToDelete)); };
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

    auto deleter = [](T* ptrToDelete) { OPM_GPU_WARN_IF_ERROR(cudaFree(ptrToDelete)); };
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
 * @brief Deleter that releases a GPU array allocation made with \c cudaMalloc.
 *
 * Used by \c make_gpu_unique_ptr_array. Having a named deleter type (rather than a lambda)
 * makes the resulting \c std::unique_ptr type stable across translation units, which is
 * important when storing it as a class member.
 */
template <class T>
struct GpuArrayDeleter {
    void operator()(T* ptr) const noexcept
    {
        OPM_GPU_WARN_IF_ERROR(cudaFree(ptr));
    }
};

/**
 * @brief Deleter for objects living in unified (managed) GPU memory.
 *
 * Calls the destructor of the contained object before releasing the allocation with
 * \c cudaFree. Used by \c make_gpu_managed_unique_ptr.
 */
template <class T>
struct GpuManagedDeleter {
    void operator()(T* ptr) const noexcept
    {
        if (ptr != nullptr) {
            ptr->~T();
            OPM_GPU_WARN_IF_ERROR(cudaFree(ptr));
        }
    }
};

/**
 * @brief Creates a unique pointer managing a GPU-allocated array of \p numElements elements.
 *
 * This function allocates memory on the GPU for \p numElements elements of type \c T using
 * \c cudaMalloc. It returns a \c std::unique_ptr that automatically releases the GPU memory
 * with \c cudaFree when destroyed.
 *
 * The returned smart pointer uses the \c T[] specialization, so element access via
 * \c operator[] is well defined; however, dereferencing the returned pointer on the host is
 * not safe because the underlying memory lives on the device.
 *
 * @tparam T The element type to allocate on the GPU.
 * @param numElements The number of elements of type \c T to allocate.
 * @return A \c std::unique_ptr to the GPU-allocated array.
 */
template <typename T>
std::unique_ptr<T[], GpuArrayDeleter<T>>
make_gpu_unique_ptr_array(std::size_t numElements)
{
    T* ptr = nullptr;
    OPM_GPU_SAFE_CALL(cudaMalloc(&ptr, numElements * sizeof(T)));
    return std::unique_ptr<T[], GpuArrayDeleter<T>>(ptr);
}

/**
 * @brief Creates a unique pointer managing GPU unified (managed) memory for a single object.
 *
 * Allocates one \c T worth of unified memory with \c cudaMallocManaged, then constructs a
 * \c T in-place using the supplied arguments via placement-new. The returned
 * \c std::unique_ptr destroys the contained object and releases the unified memory with
 * \c cudaFree when it goes out of scope.
 *
 * Unified memory is accessible from both host and device, so the returned pointer can be
 * dereferenced directly on either side. This is useful when an object must be visible to
 * GPU kernels but also needs to be touched from host code (e.g. when a device-side member
 * holds a pointer into the same unified-memory region).
 *
 * @tparam T The type of the object to construct in unified memory.
 * @tparam Args Argument types forwarded to \c T's constructor.
 * @param args Arguments forwarded to \c T's constructor.
 * @return A \c std::unique_ptr owning the unified-memory allocation.
 */
template <typename T, class... Args>
std::unique_ptr<T, GpuManagedDeleter<T>>
make_gpu_managed_unique_ptr(Args&&... args)
{
    void* raw = nullptr;
    OPM_GPU_SAFE_CALL(cudaMallocManaged(&raw, sizeof(T)));
    T* ptr = nullptr;
    try {
        ptr = new (raw) T(std::forward<Args>(args)...);
    } catch (...) {
        OPM_GPU_WARN_IF_ERROR(cudaFree(raw));
        throw;
    }
    return std::unique_ptr<T, GpuManagedDeleter<T>>(ptr);
}

/**
 * @brief Copies a value from GPU-allocated memory to the host.
 *
 * @param value A pointer to the value on the GPU.
 *
 * @return The value copied from the GPU.
 *
 * @note This function is involves a sychronization point, and should be used with care.
 */
template <class T>
T
copyFromGPU(const T* value)
{
#ifndef NDEBUG
    OPM_ERROR_IF(!Opm::gpuistl::detail::isGPUPointer(value), "The pointer is not associated with GPU memory.");
#endif
    T result;
    OPM_GPU_SAFE_CALL(cudaMemcpy(&result, value, sizeof(T), cudaMemcpyDeviceToHost));
    return result;
}

/**
 * @brief Copies a value from GPU-allocated memory to the host.
 *
 * @param value A shared pointer to the value on the GPU.
 *
 * @return The value copied from the GPU.
 *
 * @note This function is involves a sychronization point, and should be used with care.
 */
template <class T>
T
copyFromGPU(const std::shared_ptr<T>& value)
{
    return copyFromGPU(value.get());
}

/**
 * @brief Copies a value from GPU-allocated memory to the host.
 *
 * @tparam Deleter The custom deleter type.
 * @param value A unique pointer to the value on the GPU (with a custom deleter).
 *
 * @return The value copied from the GPU.
 *
 * @note This function is involves a sychronization point, and should be used with care.
 */
template <class T, class Deleter>
T
copyFromGPU(const std::unique_ptr<T, Deleter>& value)
{
    return copyFromGPU(value.get());
}

/**
 * @brief Copies a value from the host to GPU-allocated memory.
 *
 * @param value The value to copy to the GPU.
 * @param ptr A pointer to the GPU-allocated memory.
 *
 * @note This function is involves a sychronization point, and should be used with care.
 */
template <class T>
void
copyToGPU(const T& value, T* ptr)
{
#ifndef NDEBUG
    OPM_ERROR_IF(!Opm::gpuistl::detail::isGPUPointer(ptr), "The pointer is not associated with GPU memory.");
#endif
    OPM_GPU_SAFE_CALL(cudaMemcpy(ptr, &value, sizeof(T), cudaMemcpyHostToDevice));
}

/**
 * @brief Copies a value from the host to GPU-allocated memory using a shared_ptr.
 *
 * @param value The value to copy to the GPU.
 * @param ptr A shared_ptr to the GPU-allocated memory.
 *
 * @note This function involves a synchronization point, and should be used with care.
 */
template <class T>
void
copyToGPU(const T& value, const std::shared_ptr<T>& ptr)
{
    copyToGPU(value, ptr.get());
}

/**
 * @brief Copies a value from the host to GPU-allocated memory using a unique_ptr.
 *
 * @tparam Deleter The custom deleter type.
 * @param value The value to copy to the GPU.
 * @param ptr A unique_ptr to the GPU-allocated memory (with a custom deleter).
 *
 * @note This function involves a synchronization point, and should be used with care.
 */
template <class T, class Deleter>
void
copyToGPU(const T& value, const std::unique_ptr<T, Deleter>& ptr)
{
    copyToGPU(value, ptr.get());
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

    OPM_HOST_DEVICE const T& operator*() const
    {
        return *ptr_;
    }

    OPM_HOST_DEVICE T& operator*()
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

/**
 * @brief Specialization of PointerView for void type
 * This is needed beause we cannot have a PointerView<void> specialization
 * due to dereferincing a void ptr
 */
template <>
class PointerView<void>
{
public:
    PointerView(const PointerView& other) = default;

    PointerView(const std::shared_ptr<void>& ptr)
        : ptr_(ptr.get())
    {
    }

    template <class Deleter>
    PointerView(const std::unique_ptr<void, Deleter>& ptr)
        : ptr_(ptr.get())
    {
    }

    PointerView(void* ptr)
        : ptr_(ptr)
    {
    }

    OPM_HOST_DEVICE void* get() const
    {
        return ptr_;
    }

    OPM_HOST_DEVICE void* operator->() const
    {
        return ptr_;
    }

private:
    void* ptr_;
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

/**
 * @brief A value stored with a pointer interface.
 * Can be used to wrap objects in GPU kernels that were otherwise stored as pointers
 */
template<class T>
class ValueAsPointer {
public:
    using element_type = T;

    OPM_HOST_DEVICE ValueAsPointer() = default;

    OPM_HOST_DEVICE explicit ValueAsPointer(const T& t) : value(t) {}

    OPM_HOST_DEVICE T* operator->() {
        return &value;
    }

    OPM_HOST_DEVICE T* get() {
        return &value;
    }

    OPM_HOST_DEVICE const T* operator->() const {
        return &value;
    }

    OPM_HOST_DEVICE const T* get() const {
        return &value;
    }

    OPM_HOST_DEVICE T& operator*() {
        return value;
    }

    OPM_HOST_DEVICE const T& operator*() const {
        return value;
    }
private:
    T value;
};
} // namespace Opm::gpuistl
#endif
