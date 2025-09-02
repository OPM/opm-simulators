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

#ifndef OPM_HYPRE_CPU_TRANSFERS_HPP
#define OPM_HYPRE_CPU_TRANSFERS_HPP

#include <opm/simulators/linalg/gpuistl/hypreinterface/HypreDataStructures.hpp>
#include <opm/simulators/linalg/gpuistl/hypreinterface/HypreErrorHandling.hpp>

#include <HYPRE.h>
#include <_hypre_utilities.h>

namespace Opm::gpuistl::HypreInterface
{

/**
 * @brief Extract owned vector values in the order expected by HYPRE
 */
 template <typename VectorType>
 void
 setContinuousVectorForHypre(const VectorType& v,
                             std::vector<HYPRE_Real>& continuous_vector_values,
                             const std::vector<int>& local_hypre_to_local_dune)
 {
     // Set v values solution vectors
     for (size_t i = 0; i < local_hypre_to_local_dune.size(); ++i) {
         continuous_vector_values[i] = v[local_hypre_to_local_dune[i]][0];
     }
 }

 /**
  * @brief Distribute HYPRE vector values back to original vector positions
  */
 template <typename VectorType>
 void
 setDuneVectorFromContinuousVector(VectorType& v,
                                   const std::vector<HYPRE_Real>& continuous_vector_values,
                                   const std::vector<int>& local_hypre_to_local_dune)
 {
     // Place values back in original positions
     for (size_t i = 0; i < local_hypre_to_local_dune.size(); ++i) {
         v[local_hypre_to_local_dune[i]][0] = continuous_vector_values[i];
     }
 }

/**
 * @brief Transfer CPU vector to Hypre vector
 */
template <typename VectorType>
void
transferCpuVectorToHypre(const VectorType& cpu_vec,
                         HYPRE_IJVector hypre_vec,
                         HypreHostDataArrays& host_arrays,
                         const HypreDeviceDataArrays& device_arrays,
                         const ParallelInfo& par_info,
                         bool use_gpu_backend)
{
    const int N = static_cast<int>(host_arrays.indices.size());
    using T = typename VectorType::field_type;

    if (use_gpu_backend) {
#if HYPRE_USING_CUDA || HYPRE_USING_HIP
        // GPU backend with CPU input: use pre-allocated device arrays
        if (par_info.owner_first) {
            // Direct host-to-device transfer for owner-first ordering
            const T* values = &(cpu_vec[0][0]);
            hypre_TMemcpy(
                device_arrays.vector_buffer_device, values, HYPRE_Real, N, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
            OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorSetValues(
                hypre_vec, N, device_arrays.indices_device, device_arrays.vector_buffer_device));
        } else {
            // Use continuous storage and device buffer for non-owner-first ordering
            setContinuousVectorForHypre(
                cpu_vec, host_arrays.continuous_vector_values, par_info.local_hypre_to_local_dune);
            hypre_TMemcpy(device_arrays.vector_buffer_device,
                          host_arrays.continuous_vector_values.data(),
                          HYPRE_Real,
                          N,
                          HYPRE_MEMORY_DEVICE,
                          HYPRE_MEMORY_HOST);
            OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorSetValues(
                hypre_vec, N, device_arrays.indices_device, device_arrays.vector_buffer_device));
        }
#endif // HYPRE_USING_CUDA || HYPRE_USING_HIP
    } else {
        // CPU backend with CPU input: use host arrays directly
        if (par_info.owner_first) {
            // Direct transfer for owner-first ordering
            const T* values = &(cpu_vec[0][0]);
            OPM_HYPRE_SAFE_CALL(
                HYPRE_IJVectorSetValues(hypre_vec, N, const_cast<HYPRE_BigInt*>(host_arrays.indices.data()), values));
        } else {
            // Use continuous storage for non-owner-first ordering
            setContinuousVectorForHypre(
                cpu_vec, host_arrays.continuous_vector_values, par_info.local_hypre_to_local_dune);
            OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorSetValues(hypre_vec,
                                                        N,
                                                        const_cast<HYPRE_BigInt*>(host_arrays.indices.data()),
                                                        host_arrays.continuous_vector_values.data()));
        }
    }

    OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorAssemble(hypre_vec));
}

/**
 * @brief Transfer Hypre vector to CPU vector
 */
template <typename VectorType>
void
transferHypreToCpuVector(HYPRE_IJVector hypre_vec,
                         VectorType& cpu_vec,
                         HypreHostDataArrays& host_arrays,
                         const HypreDeviceDataArrays& device_arrays,
                         const ParallelInfo& par_info,
                         bool use_gpu_backend)
{
    const int N = static_cast<int>(host_arrays.indices.size());
    using T = typename VectorType::field_type;

    if (use_gpu_backend) {
#if HYPRE_USING_CUDA || HYPRE_USING_HIP
        // GPU backend with CPU input: use pre-allocated device arrays
        if (par_info.owner_first) {
            // Direct device-to-host transfer for owner-first ordering
            T* values = &(cpu_vec[0][0]);
            OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorGetValues(
                hypre_vec, N, device_arrays.indices_device, device_arrays.vector_buffer_device));
            hypre_TMemcpy(
                values, device_arrays.vector_buffer_device, HYPRE_Real, N, HYPRE_MEMORY_HOST, HYPRE_MEMORY_DEVICE);
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
            setDuneVectorFromContinuousVector(
                cpu_vec, host_arrays.continuous_vector_values, par_info.local_hypre_to_local_dune);
        }
#endif // HYPRE_USING_CUDA || HYPRE_USING_HIP
    } else {
        // CPU backend with CPU input: use host arrays directly
        if (par_info.owner_first) {
            // Direct transfer for owner-first ordering
            T* values = &(cpu_vec[0][0]);
            OPM_HYPRE_SAFE_CALL(
                HYPRE_IJVectorGetValues(hypre_vec, N, const_cast<HYPRE_BigInt*>(host_arrays.indices.data()), values));
        } else {
            // Use continuous storage for non-owner-first ordering
            OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorGetValues(hypre_vec,
                                                        N,
                                                        const_cast<HYPRE_BigInt*>(host_arrays.indices.data()),
                                                        host_arrays.continuous_vector_values.data()));
            setDuneVectorFromContinuousVector(
                cpu_vec, host_arrays.continuous_vector_values, par_info.local_hypre_to_local_dune);
        }
    }
}

/**
 * @brief Update Hypre matrix from CPU matrix
 * Uses HYPRE_IJMatrixSetValues2 with pre-computed row_indexes, which allows us to use the original CPU matrix data
 * (with potential ghost values) directly.
 */
template <typename MatrixType>
void
updateMatrixFromCpuMatrix(const MatrixType& cpu_matrix,
                          HYPRE_IJMatrix hypre_matrix,
                          const SparsityPattern& sparsity_pattern,
                          const HypreHostDataArrays& host_arrays,
                          const HypreDeviceDataArrays& device_arrays,
                          bool use_gpu_backend)
{
    const auto N = sparsity_pattern.rows.size();
    const auto nnz = cpu_matrix.nonzeroes(); // Total entries including ghost

    using T = typename MatrixType::field_type;
    const T* values = &(cpu_matrix[0][0][0][0]);

    if (use_gpu_backend) {
#if HYPRE_USING_CUDA || HYPRE_USING_HIP
        hypre_TMemcpy(
            device_arrays.matrix_buffer_device, values, HYPRE_Real, nnz, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
        OPM_HYPRE_SAFE_CALL(HYPRE_IJMatrixSetValues2(hypre_matrix,
                                                     N,
                                                     device_arrays.ncols_device,
                                                     device_arrays.rows_device,
                                                     device_arrays.row_indexes_device,
                                                     device_arrays.cols_device,
                                                     device_arrays.matrix_buffer_device));
#endif
    } else {
        OPM_HYPRE_SAFE_CALL(HYPRE_IJMatrixSetValues2(hypre_matrix,
                                                     N,
                                                     const_cast<HYPRE_Int*>(sparsity_pattern.ncols.data()),
                                                     const_cast<HYPRE_BigInt*>(sparsity_pattern.rows.data()),
                                                     const_cast<HYPRE_Int*>(host_arrays.row_indexes.data()),
                                                     const_cast<HYPRE_BigInt*>(sparsity_pattern.cols.data()),
                                                     values));
    }

    OPM_HYPRE_SAFE_CALL(HYPRE_IJMatrixAssemble(hypre_matrix));
}
} // namespace Opm::gpuistl::HypreInterface

#endif // OPM_HYPRE_CPU_TRANSFERS_HPP
