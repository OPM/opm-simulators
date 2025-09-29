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

#ifndef OPM_HYPRE_GPU_TRANSFERS_HPP
#define OPM_HYPRE_GPU_TRANSFERS_HPP

#include <opm/simulators/linalg/gpuistl/hypreinterface/HypreDataStructures.hpp>
#include <opm/simulators/linalg/gpuistl/hypreinterface/HypreErrorHandling.hpp>

#if HAVE_CUDA
#if USE_HIP
#include <opm/simulators/linalg/gpuistl_hip/GpuSparseMatrixWrapper.hpp>
#include <opm/simulators/linalg/gpuistl_hip/GpuVector.hpp>
#else
#include <opm/simulators/linalg/gpuistl/GpuSparseMatrixWrapper.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>
#endif
#endif // HAVE_CUDA

#include <HYPRE.h>
#include <_hypre_utilities.h>

namespace Opm::gpuistl::HypreInterface
{

template <typename VectorType>
void setContinuousGpuVectorForHypre(const VectorType& v,
                                    std::vector<HYPRE_Real>& continuous_vector_values,
                                    const std::vector<int>& local_hypre_to_local_dune);

template <typename VectorType>
void setGpuVectorFromContinuousVector(VectorType& v,
                                      const std::vector<HYPRE_Real>& continuous_vector_values,
                                      const std::vector<int>& local_hypre_to_local_dune);
#if HYPRE_USING_CUDA || HYPRE_USING_HIP

/**
 * @brief Transfer GPU vector to Hypre vector
 */
template <typename VectorType>
void
transferGpuVectorToHypre(const VectorType& gpu_vec,
                         HYPRE_IJVector hypre_vec,
                         HypreHostDataArrays& host_arrays,
                         const HypreDeviceDataArrays& device_arrays,
                         const ParallelInfo& par_info,
                         bool use_gpu_backend)
{
    const int N = static_cast<int>(host_arrays.indices.size());
    using T = typename VectorType::field_type;

    if (use_gpu_backend) {
        // GPU backend with GPU input: use pre-allocated device arrays
        if (par_info.owner_first) {
            // Direct device-to-device transfer for owner-first ordering
            const T* device_ptr = gpu_vec.data();
            OPM_HYPRE_SAFE_CALL(
                HYPRE_IJVectorSetValues(hypre_vec, N, device_arrays.indices_device, const_cast<T*>(device_ptr)));
        } else {
            // Use continuous storage and device buffer for non-owner-first ordering
            setContinuousGpuVectorForHypre(
                gpu_vec, host_arrays.continuous_vector_values, par_info.local_hypre_to_local_dune);
            hypre_TMemcpy(device_arrays.vector_buffer_device,
                          host_arrays.continuous_vector_values.data(),
                          HYPRE_Real,
                          N,
                          HYPRE_MEMORY_DEVICE,
                          HYPRE_MEMORY_HOST);
            OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorSetValues(
                hypre_vec, N, device_arrays.indices_device, device_arrays.vector_buffer_device));
        }
    } else {
        // CPU backend with GPU input: copy via host memory
        if (par_info.owner_first) {
            // Get values to host and then set
            auto host_values = gpu_vec.asStdVector();
            OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorSetValues(hypre_vec,
                                                        N,
                                                        const_cast<HYPRE_BigInt*>(host_arrays.indices.data()),
                                                        reinterpret_cast<HYPRE_Real*>(host_values.data())));
        } else {
            // Use continuous storage for non-owner-first ordering
            setContinuousGpuVectorForHypre(
                gpu_vec, host_arrays.continuous_vector_values, par_info.local_hypre_to_local_dune);
            OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorSetValues(hypre_vec,
                                                        N,
                                                        const_cast<HYPRE_BigInt*>(host_arrays.indices.data()),
                                                        host_arrays.continuous_vector_values.data()));
        }
    }
    OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorAssemble(hypre_vec));
}

/**
 * @brief Transfer Hypre vector to GPU vector
 */
template <typename VectorType>
void
transferHypreToGpuVector(HYPRE_IJVector hypre_vec,
                         VectorType& gpu_vec,
                         HypreHostDataArrays& host_arrays,
                         const HypreDeviceDataArrays& device_arrays,
                         const ParallelInfo& par_info,
                         bool use_gpu_backend)
{
    const int N = static_cast<int>(host_arrays.indices.size());
    using T = typename VectorType::field_type;

    if (use_gpu_backend) {
        // GPU backend with GPU input: use pre-allocated device arrays
        if (par_info.owner_first) {
            // Direct device-to-device transfer for owner-first ordering
            T* device_ptr = gpu_vec.data();
            OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorGetValues(hypre_vec, N, device_arrays.indices_device, device_ptr));
        } else {
            // Use device buffer and then remap for non-owner-first ordering
            OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorGetValues(
                hypre_vec, N, device_arrays.indices_device, device_arrays.vector_buffer_device));
            hypre_TMemcpy(host_arrays.continuous_vector_values.data(),
                          device_arrays.vector_buffer_device,
                          HYPRE_Real,
                          N,
                          HYPRE_MEMORY_HOST,
                          HYPRE_MEMORY_DEVICE);
            setGpuVectorFromContinuousVector(
                gpu_vec, host_arrays.continuous_vector_values, par_info.local_hypre_to_local_dune);
        }
    } else {
        // CPU backend with GPU input: copy via host memory
        if (par_info.owner_first) {
            // Get values to host and then copy to GPU
            auto host_values = gpu_vec.asStdVector();
            OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorGetValues(hypre_vec,
                                                        N,
                                                        const_cast<HYPRE_BigInt*>(host_arrays.indices.data()),
                                                        reinterpret_cast<HYPRE_Real*>(host_values.data())));
            gpu_vec = VectorType(host_values);
        } else {
            // Use continuous storage for non-owner-first ordering
            OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorGetValues(hypre_vec,
                                                        N,
                                                        const_cast<HYPRE_BigInt*>(host_arrays.indices.data()),
                                                        host_arrays.continuous_vector_values.data()));
            setGpuVectorFromContinuousVector(
                gpu_vec, host_arrays.continuous_vector_values, par_info.local_hypre_to_local_dune);
        }
    }
}

/**
 * @brief Update Hypre matrix values from GPU matrix using pre-allocated helper arrays
 * Uses HYPRE_IJMatrixSetValues2 with pre-computed row_indexes, which allows us to use the original GPU matrix data
 * (with potential ghost values) directly.
 */
template <typename MatrixType>
void
updateMatrixFromGpuSparseMatrix(const MatrixType& gpu_matrix,
                                HYPRE_IJMatrix hypre_matrix,
                                const SparsityPattern& sparsity_pattern,
                                const HypreHostDataArrays& host_arrays,
                                const HypreDeviceDataArrays& device_arrays,
                                bool use_gpu_backend)
{
    const auto N = sparsity_pattern.rows.size();
    using T = typename MatrixType::field_type;
    const T* values = gpu_matrix.getNonZeroValues().data();

    if (use_gpu_backend) {
        // GPU backend with GPU input: use pre-allocated device arrays
        // Direct device-to-device transfer using smart row_indexes
        OPM_HYPRE_SAFE_CALL(HYPRE_IJMatrixSetValues2(hypre_matrix,
                                                     N,
                                                     device_arrays.ncols_device,
                                                     device_arrays.rows_device,
                                                     device_arrays.row_indexes_device,
                                                     device_arrays.cols_device,
                                                     values));
    } else {
        // CPU backend with GPU input: copy to host first
        auto host_values = gpu_matrix.getNonZeroValues().asStdVector();
        OPM_HYPRE_SAFE_CALL(HYPRE_IJMatrixSetValues2(hypre_matrix,
                                                     N,
                                                     const_cast<HYPRE_Int*>(sparsity_pattern.ncols.data()),
                                                     const_cast<HYPRE_BigInt*>(sparsity_pattern.rows.data()),
                                                     const_cast<HYPRE_Int*>(host_arrays.row_indexes.data()),
                                                     const_cast<HYPRE_BigInt*>(sparsity_pattern.cols.data()),
                                                     host_values.data()));
    }
    OPM_HYPRE_SAFE_CALL(HYPRE_IJMatrixAssemble(hypre_matrix));
}

template <typename VectorType>
void
setContinuousGpuVectorForHypre(const VectorType& v,
                               std::vector<HYPRE_Real>& continuous_vector_values,
                               const std::vector<int>& local_hypre_to_local_dune)
{
    // Get vector data to host first
    auto host_values = v.asStdVector();
    // Set values using the mapping
    for (size_t i = 0; i < local_hypre_to_local_dune.size(); ++i) {
        continuous_vector_values[i] = host_values[local_hypre_to_local_dune[i]];
    }
}

template <typename VectorType>
void
setGpuVectorFromContinuousVector(VectorType& v,
                                 const std::vector<HYPRE_Real>& continuous_vector_values,
                                 const std::vector<int>& local_hypre_to_local_dune)
{
    // Copy values to host and update values with mapping
    auto host_values = v.asStdVector();
    for (size_t i = 0; i < local_hypre_to_local_dune.size(); ++i) {
        host_values[local_hypre_to_local_dune[i]] = continuous_vector_values[i];
    }
    // Copy back to GPU
    v = VectorType(host_values);
}
#endif // HYPRE_USING_CUDA || HYPRE_USING_HIP

} // namespace Opm::gpuistl::HypreInterface

#endif // OPM_HYPRE_GPU_TRANSFERS_HPP
