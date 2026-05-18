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

#ifndef OPM_HYPRE_INTERFACE_HPP
#define OPM_HYPRE_INTERFACE_HPP

#include <opm/simulators/linalg/hypreinterface/HypreCpuTransfers.hpp>
#include <opm/simulators/linalg/hypreinterface/HypreDataStructures.hpp>
#include <opm/simulators/linalg/hypreinterface/HypreErrorHandling.hpp>
#include <opm/simulators/linalg/hypreinterface/HypreSetup.hpp>
#include <opm/simulators/linalg/hypreinterface/HypreUtils.hpp>

// Fix conflict with Dune's matrixmarket.hh header - HYPRE defines MM_MAX_LINE_LENGTH as macro
// but Dune expects it as enum value
#ifdef MM_MAX_LINE_LENGTH
#undef MM_MAX_LINE_LENGTH
#endif

#if HYPRE_USING_CUDA
#include <opm/simulators/linalg/gpuistl/detail/gpu_type_detection.hpp>
#include <opm/simulators/linalg/gpuistl/hypreinterface/HypreGpuTransfers.hpp>
#elif HYPRE_USING_HIP
#include <opm/simulators/linalg/gpuistl_hip/detail/gpu_type_detection.hpp>
#include <opm/simulators/linalg/gpuistl_hip/hypreinterface/HypreGpuTransfers.hpp>
#endif

/**
 * @brief Unified interface for Hypre operations with both CPU and GPU data structures
 *
 * This namespace provides utilities for working with Hypre resources and transferring
 * data between CPU/GPU data structures and Hypre handles. It handles type detection
 * and automatically chooses the most efficient transfer method:
 * - For GPU types: Zero-copy device-to-device transfers when using GPU backend
 * - For CPU types: Host-to-device transfers when using GPU backend
 * - For CPU types with CPU backend: Direct host memory usage
 *
 * Supports four use cases:
 * 1. Input type is CPU and backend acceleration is CPU
 * 2. Input type is CPU and backend acceleration is GPU
 * 3. Input type is GPU and backend acceleration is GPU
 * 4. Input type is GPU and backend acceleration is CPU
 *
 * @note Error handling: All functions throw HypreError exceptions when Hypre operations fail.
 * @note This is a consolidated version that includes all functionality across multiple files.
 */
namespace Opm::linalg::HypreInterface
{

    /**
     * @brief Initialize the Hypre library and set memory/execution policy
     */
    void initialize(bool use_gpu_backend);

    /**
     * @brief Create Hypre solver (BoomerAMG)
     */
    HYPRE_Solver createAMGSolver();

    /**
     * @brief Set solver parameters from property tree
     */
    void setSolverParameters(HYPRE_Solver solver, const PropertyTree& prm, bool use_gpu_backend);

    /**
     * @brief Create Hypre matrix
     */
    template <typename CommType>
    HYPRE_IJMatrix createMatrix(HYPRE_Int N, HYPRE_Int dof_offset, const CommType& comm);

    /**
     * @brief Create Hypre vector
     */
    template <typename CommType>
    HYPRE_IJVector createVector(HYPRE_Int N, HYPRE_Int dof_offset, const CommType& comm);

    /**
     * @brief Destroy Hypre solver
     */
    void destroySolver(HYPRE_Solver solver);

    /**
     * @brief Destroy Hypre matrix
     */
    void destroyMatrix(HYPRE_IJMatrix matrix);

    /**
     * @brief Destroy Hypre vector
     */
    void destroyVector(HYPRE_IJVector vector);

    /**
     * @brief Setup parallel information for Hypre (automatically detects serial/parallel)
     */
    template <typename CommType, typename MatrixType>
    ParallelInfo
    setupHypreParallelInfo(const CommType& comm, const MatrixType& matrix);

    /**
     * @brief Setup sparsity pattern from matrix (automatically detects CPU/GPU type)
     */
    template <typename MatrixType>
    SparsityPattern
    setupSparsityPattern(const MatrixType& matrix,
                         const ParallelInfo& par_info,
                         bool owner_first);

    /**
     * @brief Compute row indexes for HYPRE_IJMatrixSetValues2
     */
    template <typename MatrixType>
    std::vector<HYPRE_Int> computeRowIndexes(const MatrixType& matrix,
                                             const std::vector<HYPRE_Int>& ncols,
                                             const std::vector<int>& local_dune_to_local_hypre,
                                             bool owner_first);

    /**
     * @brief Transfer vector to Hypre from any vector type (CPU or GPU)
     */
    template <typename VectorType>
    void transferVectorToHypre(const VectorType& vec,
                               HYPRE_IJVector hypre_vec,
                               HostDataArrays& host_arrays,
                               const DeviceDataArrays& device_arrays,
                               const ParallelInfo& par_info,
                               bool use_gpu_backend);

    /**
     * @brief Transfer vector from Hypre to any vector type (CPU or GPU)
     */
    template <typename VectorType>
    void transferVectorFromHypre(HYPRE_IJVector hypre_vec,
                                 VectorType& vec,
                                 HostDataArrays& host_arrays,
                                 const DeviceDataArrays& device_arrays,
                                 const ParallelInfo& par_info,
                                 bool use_gpu_backend);

    /**
     * @brief Update matrix values in Hypre
     */
    template <typename MatrixType>
    void updateMatrixValues(const MatrixType& matrix,
                            HYPRE_IJMatrix hypre_matrix,
                            const SparsityPattern& sparsity_pattern,
                            const HostDataArrays& host_arrays,
                            const DeviceDataArrays& device_arrays,
                            bool use_gpu_backend);

    template <typename VectorType>
    void transferVectorToHypre(const VectorType& vec,
                               HYPRE_IJVector hypre_vec,
                               HostDataArrays& host_arrays,
                               const DeviceDataArrays& device_arrays,
                               const ParallelInfo& par_info,
                               bool use_gpu_backend)
    {
#if HYPRE_USING_CUDA || HYPRE_USING_HIP
        if constexpr (gpuistl::is_gpu_type<VectorType>::value) {
            gpuistl::HypreInterface::transferGpuVectorToHypre(vec, hypre_vec, host_arrays,
                                                              device_arrays, par_info, use_gpu_backend);
        } else
#endif
        {
            transferCpuVectorToHypre(vec, hypre_vec, host_arrays,
                                     device_arrays, par_info, use_gpu_backend);
        }
    }

    template <typename VectorType>
    void transferVectorFromHypre(HYPRE_IJVector hypre_vec,
                                 VectorType& vec,
                                 HostDataArrays& host_arrays,
                                 const DeviceDataArrays& device_arrays,
                                 const ParallelInfo& par_info,
                                 bool use_gpu_backend)
    {
#if HYPRE_USING_CUDA || HYPRE_USING_HIP
        if constexpr (gpuistl::is_gpu_type<VectorType>::value) {
            gpuistl::HypreInterface::transferHypreToGpuVector(hypre_vec, vec, host_arrays,
                                                              device_arrays, par_info, use_gpu_backend);
        } else
#endif
        {
            transferHypreToCpuVector(hypre_vec, vec, host_arrays, device_arrays,
                                     par_info, use_gpu_backend);
        }
    }

    template <typename MatrixType>
    void updateMatrixValues(const MatrixType& matrix,
                            HYPRE_IJMatrix hypre_matrix,
                            const SparsityPattern& sparsity_pattern,
                            const HostDataArrays& host_arrays,
                            const DeviceDataArrays& device_arrays,
                            bool use_gpu_backend)
    {
#if HYPRE_USING_CUDA || HYPRE_USING_HIP
        if constexpr (gpuistl::is_gpu_type<MatrixType>::value) {
            gpuistl::HypreInterface::updateMatrixFromGpuSparseMatrix(matrix,
                                                                     hypre_matrix,
                                                                     sparsity_pattern,
                                                                     host_arrays,
                                                                     device_arrays,
                                                                     use_gpu_backend);
        } else
#endif
        {
            updateMatrixFromCpuMatrix(matrix, hypre_matrix, sparsity_pattern,
                                      host_arrays, device_arrays, use_gpu_backend);
        }
    }

} // namespace Opm::linalg::HypreInterface

#endif // OPM_HYPRE_INTERFACE_HPP
