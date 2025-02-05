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

#ifndef OPM_SIMULATORS_LINALG_GPUISTL_GPU_RESOURCES_HPP
#define OPM_SIMULATORS_LINALG_GPUISTL_GPU_RESOURCES_HPP
#include <cuda.h>
#include <cuda_runtime.h>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>
#include <type_traits>

/**
 * @brief Creates an RAII-style GPU resource class.
 *
 * This macro automatically generates a class that:
 * - Defines a resource type and a corresponding instance.
 * - Initializes the resource in the constructor via the provided creation function.
 * - Cleans up the resource in the destructor via the provided destruction function.
 *
 * @tparam name   The class name to be defined.
 * @tparam type   The underlying resource type to manage.
 * @tparam create The function used to create or allocate the resource.
 * @tparam destroy The function used to destroy or deallocate the resource.
 * @param ...     Additional arguments forwarded to the creation function.
 *
 * Usage example:
 * @code
 * OPM_CREATE_GPU_RESOURCE(MyResource, ResourceType, createFunction, destroyFunction, args...);
 * @endcode
 *
 * This ensures correct setup and teardown of GPU resources without clutter in user code.
 */
#define OPM_CREATE_GPU_RESOURCE(name, type, create, destroy, ...)                                                      \
    class name                                                                                                         \
    {                                                                                                                  \
    public:                                                                                                            \
        using resource_type = type;                                                                                    \
                                                                                                                       \
        name()                                                                                                         \
        {                                                                                                              \
            OPM_GPU_SAFE_CALL(create(&resource_, ##__VA_ARGS__));                                                      \
        }                                                                                                              \
        name(const name&) = delete;                                                                                    \
        name& operator=(const name&) = delete;                                                                         \
                                                                                                                       \
        ~name()                                                                                                        \
        {                                                                                                              \
            OPM_GPU_WARN_IF_ERROR(destroy(resource_));                                                                 \
        }                                                                                                              \
                                                                                                                       \
        const resource_type& get() const                                                                               \
        {                                                                                                              \
            return resource_;                                                                                          \
        }                                                                                                              \
        resource_type& get()                                                                                           \
        {                                                                                                              \
            return resource_;                                                                                          \
        }                                                                                                              \
                                                                                                                       \
    private:                                                                                                           \
        resource_type resource_ = nullptr;                                                                             \
    }

/**
 * @brief A helper macro for creating a GPU resource wrapper without invoking creation logic.
 *
 * This macro, closely related to \c OPM_CREATE_GPU_RESOURCE, provides a class that manages a GPU
 * resource of the specified type. Upon destruction, it automatically destroys the underlying
 * resource to ensure proper cleanup of GPU memory or handles.
 *
 * Usage:
 * - Construct an instance of this wrapper to store a GPU resource.
 * - The resource is destroyed when the wrapper goes out of scope.
 *
 * Use \c OPM_CREATE_GPU_RESOURCE for variants that require resource creation logic.
 *
 * @tparam name Name of the wrapper class to manage the GPU resource.
 * @tparam type The type of the GPU resource to manage.
 * @tparam destroy The function or macro used to destroy the resource.
 */
#define OPM_CREATE_GPU_RESOURCE_NO_CREATE(name, type, destroy)                                                         \
    class name                                                                                                         \
    {                                                                                                                  \
    public:                                                                                                            \
        using resource_type = type;                                                                                    \
                                                                                                                       \
        name() = default;                                                                                              \
                                                                                                                       \
        ~name()                                                                                                        \
        {                                                                                                              \
            OPM_GPU_WARN_IF_ERROR(destroy(resource_));                                                                 \
        }                                                                                                              \
        name(const name&) = delete;                                                                                    \
        name& operator=(const name&) = delete;                                                                         \
                                                                                                                       \
        const resource_type& get() const                                                                               \
        {                                                                                                              \
            return resource_;                                                                                          \
        }                                                                                                              \
        resource_type& get()                                                                                           \
        {                                                                                                              \
            return resource_;                                                                                          \
        }                                                                                                              \
                                                                                                                       \
    private:                                                                                                           \
        resource_type resource_;                                                                                       \
    }

namespace Opm::gpuistl
{
/**
 * @brief Manages a CUDA stream resource.
 *
 * This resource encapsulates a cudaStream_t handle and provides
 * automatic creation and destruction of the CUDA stream.
 * Use this resource to schedule and synchronize GPU kernels
 * or other asynchronous operations.
 */
OPM_CREATE_GPU_RESOURCE(GPUStream, cudaStream_t, cudaStreamCreate, cudaStreamDestroy);

/**
 * @brief Manages a CUDA event resource.
 *
 * This resource encapsulates a cudaEvent_t handle and provides
 * automatic creation and destruction of the CUDA event.
 * Use this resource to measure elapsed time or synchronize
 * GPU executions between different streams.
 */
OPM_CREATE_GPU_RESOURCE(GPUEvent, cudaEvent_t, cudaEventCreate, cudaEventDestroy);

/**
 * @brief Manages a CUDA graph resource.
 *
 * This resource encapsulates a cudaGraph_t handle and provides
 * automatic creation and destruction of a CUDA graph.
 * It represents a series of operations captured for
 * efficient replay, execution, or modification.
 */
OPM_CREATE_GPU_RESOURCE(
    GPUGraph, cudaGraph_t, cudaGraphCreate, cudaGraphDestroy, 0); // Per documentation, needs 0 as last argument

/**
 * @brief Manages a CUDA graph execution resource.
 *
 * This resource encapsulates a cudaGraphExec_t handle and provides
 * automatic destruction of the CUDA graph execution object.
 * It represents the compiled and optimized version of a CUDA graph
 * ready for efficient execution.
 */
OPM_CREATE_GPU_RESOURCE_NO_CREATE(GPUGraphExec, cudaGraphExec_t, cudaGraphExecDestroy);
} // namespace Opm::gpuistl
#endif
