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

#include <opm/simulators/linalg/hypreinterface/HypreDataStructures.hpp>
#include <opm/simulators/linalg/hypreinterface/HypreErrorHandling.hpp>

#if HYPRE_USING_CUDA
#include <opm/simulators/linalg/gpuistl/hypreinterface/HypreCpuTransfers.hpp>
#elif HYPRE_USING_HIP
#include <opm/simulators/linalg/gpuistl_hip/hypreinterface/HypreCpuTransfers.hpp>
#endif

#include <HYPRE.h>
#include <_hypre_utilities.h>

namespace Opm::linalg::HypreInterface
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
                         linalg::HypreInterface::HostDataArrays& host_arrays,
                         [[maybe_unused]] const linalg::HypreInterface::DeviceDataArrays& device_arrays,
                         const linalg::HypreInterface::ParallelInfo& par_info,
                         bool use_gpu_backend)
{
    const int N = static_cast<int>(host_arrays.indices.size());
    using T = typename VectorType::field_type;

    if (use_gpu_backend) {
#if HYPRE_USING_CUDA || HYPRE_USING_HIP
        gpuistl::HypreInterface::transferCpuVectorToHypre(cpu_vec, hypre_vec, host_arrays,
                                                          device_arrays, par_info);
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
                         linalg::HypreInterface::HostDataArrays& host_arrays,
                         [[maybe_unused]] const linalg::HypreInterface::DeviceDataArrays& device_arrays,
                         const linalg::HypreInterface::ParallelInfo& par_info,
                         bool use_gpu_backend)
{
    const int N = static_cast<int>(host_arrays.indices.size());
    using T = typename VectorType::field_type;

    if (use_gpu_backend) {
#if HYPRE_USING_CUDA || HYPRE_USING_HIP
        gpuistl::HypreInterface::transferHypreToCpuVector(hypre_vec, cpu_vec, host_arrays,
                                                          device_arrays, par_info);
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
                          const linalg::HypreInterface::SparsityPattern& sparsity_pattern,
                          const linalg::HypreInterface::HostDataArrays& host_arrays,
                          [[maybe_unused]] const linalg::HypreInterface::DeviceDataArrays& device_arrays,
                          bool use_gpu_backend)
{
    const auto N = sparsity_pattern.rows.size();

    using T = typename MatrixType::field_type;
    const T* values = &(cpu_matrix[0][0][0][0]);

    if (use_gpu_backend) {
#if HYPRE_USING_CUDA || HYPRE_USING_HIP
        gpuistl::HypreInterface::updateMatrixFromCpuMatrix(cpu_matrix, hypre_matrix, sparsity_pattern, device_arrays);
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

} // namespace Opm::linalg::HypreInterface

#endif // OPM_HYPRE_CPU_TRANSFERS_HPP
