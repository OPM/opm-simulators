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

#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_type_detection.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>

#include <dune/istl/paamg/pinfo.hh>

#if HAVE_CUDA
#if USE_HIP
#include <opm/simulators/linalg/gpuistl_hip/GpuVector.hpp>
#include <opm/simulators/linalg/gpuistl_hip/GpuSparseMatrix.hpp>
#include <opm/simulators/linalg/gpuistl_hip/detail/safe_conversion.hpp>
#else
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>
#include <opm/simulators/linalg/gpuistl/GpuSparseMatrix.hpp>
#include <opm/simulators/linalg/gpuistl/detail/safe_conversion.hpp>
#endif
#endif // HAVE_CUDA

#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>
#include <HYPRE_krylov.h>
#include <_hypre_utilities.h>

// Fix conflict with Dune's matrixmarket.hh header - HYPRE defines MM_MAX_LINE_LENGTH as macro
// but Dune expects it as enum value
#ifdef MM_MAX_LINE_LENGTH
#undef MM_MAX_LINE_LENGTH
#endif

#include <fmt/core.h>

#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace Opm::gpuistl {

/**
 * @brief Exception class for Hypre errors
 */
class HypreError : public std::runtime_error {
public:
    explicit HypreError(const std::string& msg) : std::runtime_error(msg) {}
};

/**
 * @brief Get a descriptive error message for a Hypre error code
 *
 * @param err The Hypre error code
 * @param expression The Hypre expression that caused the error
 * @param file The source file where the error occurred
 * @param function The function where the error occurred
 * @param line The line number where the error occurred
 * @return std::string The formatted error message
 */
inline std::string getHypreErrorMessage(HYPRE_Int err,
                                        const std::string& expression,
                                        const std::string& file,
                                        const std::string& function,
                                        int line) {
    return fmt::format("Hypre error in expression: {}\nError code: {}\nLocation: {}:{} in function {}",
                       expression,
                       err,
                       file,
                       line,
                       function);
}

/**
 * @brief Safe call wrapper for Hypre functions
 *
 * Checks the return code from Hypre functions and throws a HypreError if an error occurred.
 *
 * @param rc The Hypre return code to check
 * @param expression The expression being evaluated (for error reporting)
 * @param file The source file (typically __FILE__)
 * @param function The function name (typically __func__)
 * @param line The line number (typically __LINE__)
 * @throws HypreError if the return code indicates an error
 */
inline void hypreSafeCall(HYPRE_Int rc,
                          const std::string& expression,
                          const std::string& file,
                          const std::string& function,
                          int line) {
    if (rc != 0) {
        throw HypreError(getHypreErrorMessage(rc, expression, file, function, line));
    }
}

/**
 * @brief Macro to wrap Hypre function calls with error checking
 *
 * Example usage:
 * @code
 * OPM_HYPRE_SAFE_CALL(HYPRE_BoomerAMGCreate(&solver));
 * @endcode
 */
#define OPM_HYPRE_SAFE_CALL(expr) \
    ::Opm::gpuistl::hypreSafeCall((expr), #expr, __FILE__, __func__, __LINE__)

/**
 * @brief Short form macro for Hypre function calls (for backward compatibility)
 */
#ifndef HYPRE_SAFE_CALL
#define HYPRE_SAFE_CALL(expr) OPM_HYPRE_SAFE_CALL(expr)
#endif



/**
 * @brief Unified interface for Hypre operations with both CPU and GPU data structures
 *
 * This class provides utilities for working with Hypre resources and transferring
 * data between CPU/GPU data structures and Hypre handles. It handles type detection
 * and automatically chooses the most efficient transfer method:
 * - For GPU types: Zero-copy device-to-device transfers when using GPU backend
 * - For CPU types: Host-to-device transfers when using GPU backend
 * - For CPU types with CPU backend: Direct host memory usage
 *
 * The class supports three use cases:
 * 1. Input type is CPU and backend acceleration is CPU
 * 2. Input type is CPU and backend acceleration is GPU
 * 3. Input type is GPU and backend acceleration is GPU
 *
 * @note Error handling: All methods throw HypreError exceptions when Hypre operations fail.
 */
class HypreInterface {
public:
    //--------------------------------------------------------------------------
    // Hypre Library Lifecycle Management
    //--------------------------------------------------------------------------

    /**
     * @brief Initialize the Hypre library and set memory/execution policy
     *
     * @param use_gpu_backend Whether to use GPU backend acceleration
     * @throws HypreError if initialization fails
     */
    static void initialize(bool use_gpu_backend) {
        // Set memory location and execution policy
#if HYPRE_USING_CUDA || HYPRE_USING_HIP
        if (use_gpu_backend) {
            OPM_HYPRE_SAFE_CALL(HYPRE_SetMemoryLocation(HYPRE_MEMORY_DEVICE));
            OPM_HYPRE_SAFE_CALL(HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE));
            // use hypre's SpGEMM instead of vendor implementation
            OPM_HYPRE_SAFE_CALL(HYPRE_SetSpGemmUseVendor(false));
            // use cuRand for PMIS
            OPM_HYPRE_SAFE_CALL(HYPRE_SetUseGpuRand(1));
            OPM_HYPRE_SAFE_CALL(HYPRE_DeviceInitialize());
            OPM_HYPRE_SAFE_CALL(HYPRE_PrintDeviceInfo());
        }
        else
#endif
        {
            OPM_HYPRE_SAFE_CALL(HYPRE_SetMemoryLocation(HYPRE_MEMORY_HOST));
            OPM_HYPRE_SAFE_CALL(HYPRE_SetExecutionPolicy(HYPRE_EXEC_HOST));
        }
    }

    //--------------------------------------------------------------------------
    // Hypre Resource Creation and Management
    //--------------------------------------------------------------------------

    /**
     * @brief Create Hypre solver (BoomerAMG)
     *
     * @return HYPRE_Solver The created solver handle
     * @throws HypreError if solver creation fails
     */
    static HYPRE_Solver createSolver() {
        HYPRE_Solver solver;
        OPM_HYPRE_SAFE_CALL(HYPRE_BoomerAMGCreate(&solver));
        return solver;
    }

    static void setSolverParameters(HYPRE_Solver solver, const PropertyTree& prm, bool use_gpu_backend) {
        // Set parameters from property tree with defaults
        HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetPrintLevel(solver, prm.get<int>("print_level", 0)));
        HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetMaxIter(solver, prm.get<int>("max_iter", 1)));
        HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetStrongThreshold(solver, prm.get<double>("strong_threshold", 0.5)));
        HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetAggTruncFactor(solver, prm.get<double>("agg_trunc_factor", 0.3)));
        HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetInterpType(solver, prm.get<int>("interp_type", 6)));
        HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetMaxLevels(solver, prm.get<int>("max_levels", 15)));
        HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetTol(solver, prm.get<double>("tolerance", 0.0)));

        if (use_gpu_backend) {
            HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetRelaxType(solver, prm.get<int>("relax_type", 16)));
            HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetCoarsenType(solver, prm.get<int>("coarsen_type", 8)));
            HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetAggNumLevels(solver, prm.get<int>("agg_num_levels", 0)));
            HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetAggInterpType(solver, prm.get<int>("agg_interp_type", 6)));
            // Keep transpose to avoid SpMTV
            HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetKeepTranspose(solver, true));
        }
        else {
            HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetRelaxType(solver, prm.get<int>("relax_type", 13)));
            HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetCoarsenType(solver, prm.get<int>("coarsen_type", 10)));
            HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetAggNumLevels(solver, prm.get<int>("agg_num_levels", 1)));
            HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetAggInterpType(solver, prm.get<int>("agg_interp_type", 4)));
        }
    }

    /**
     * @brief Create Hypre matrix
     *
     * @param N Number of rows/columns
     * @param dof_offset Offset for the matrix
     * @param comm The communicator
     * @return HYPRE_IJMatrix The created matrix handle
     * @throws HypreError if matrix creation fails
     */
    template<typename Commun>
    static HYPRE_IJMatrix createMatrix(HYPRE_Int N, HYPRE_Int dof_offset, const Commun& comm) {
        HYPRE_IJMatrix matrix;
        MPI_Comm mpi_comm;
        if constexpr (std::is_same_v<Commun, Dune::Amg::SequentialInformation>) {
            mpi_comm = MPI_COMM_SELF;
        } else {
            mpi_comm = comm.communicator();
        }
        OPM_HYPRE_SAFE_CALL(HYPRE_IJMatrixCreate(mpi_comm, dof_offset, dof_offset + (N-1), dof_offset, dof_offset + (N-1), &matrix));
        OPM_HYPRE_SAFE_CALL(HYPRE_IJMatrixSetObjectType(matrix, HYPRE_PARCSR));
        OPM_HYPRE_SAFE_CALL(HYPRE_IJMatrixInitialize(matrix));
        return matrix;
    }

    /**
     * @brief Create Hypre vector
     *
     * @param N Size of the vector
     * @param dof_offset Offset for the vector
     * @param comm The communicator
     * @return HYPRE_IJVector The created vector handle
     * @throws HypreError if vector creation fails
     */
    template<typename Commun>
    static HYPRE_IJVector createVector(HYPRE_Int N, HYPRE_Int dof_offset, const Commun& comm) {
        HYPRE_IJVector vector;
        MPI_Comm mpi_comm;
        if constexpr (std::is_same_v<Commun, Dune::Amg::SequentialInformation>) {
            mpi_comm = MPI_COMM_SELF;
        } else {
            mpi_comm = comm.communicator();
        }
        OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorCreate(mpi_comm, dof_offset, dof_offset + (N-1), &vector));
        OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorSetObjectType(vector, HYPRE_PARCSR));
        OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorInitialize(vector));
        return vector;
    }

    /**
     * @brief Destroy Hypre solver
     *
     * @param solver The solver handle to destroy
     * @throws HypreError if solver destruction fails
     */
    static void destroySolver(HYPRE_Solver solver) {
        if (solver) {
            OPM_HYPRE_SAFE_CALL(HYPRE_BoomerAMGDestroy(solver));
        }
    }

    /**
     * @brief Destroy Hypre matrix
     *
     * @param matrix The matrix handle to destroy
     * @throws HypreError if matrix destruction fails
     */
    static void destroyMatrix(HYPRE_IJMatrix matrix) {
        if (matrix) {
            OPM_HYPRE_SAFE_CALL(HYPRE_IJMatrixDestroy(matrix));
        }
    }

    /**
     * @brief Destroy Hypre vector
     *
     * @param vector The vector handle to destroy
     * @throws HypreError if vector destruction fails
     */
    static void destroyVector(HYPRE_IJVector vector) {
        if (vector) {
            OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorDestroy(vector));
        }
    }

    //--------------------------------------------------------------------------
    // Vector Operations
    //--------------------------------------------------------------------------

    /**
     * @brief Transfer vector to Hypre from any vector type (CPU or GPU)
     *
     * Selects the transfer method based on the vector type and backend configuration.
     *
     * @param vec Source vector (typically BlockVector (CPU) or GpuVector (GPU))
     * @param hypre_vec Destination Hypre vector
     * @param use_gpu_backend Whether the Hypre backend is using GPU acceleration
     * @param indices Pre-allocated host array for owned DOF indices
     * @param indices_device Pre-allocated device array for owned DOF indices (for GPU backend)
     * @param vector_buffer_device Pre-allocated device buffer for vector values (for GPU backend)
     * @param continuous_vector_values Host buffer for continuous vector values
     * @param local_hypre_to_local_dune Mapping from HYPRE indices to Dune indices (used for non-owner-first)
     * @param owner_first Whether owned DOFs come first in the ordering
     * @throws HypreError if the transfer fails
     */
    template<typename VectorType>
    static void transferVectorToHypre(const VectorType& vec, HYPRE_IJVector hypre_vec, bool use_gpu_backend,
                                    const std::vector<HYPRE_BigInt>& indices,
                                    HYPRE_BigInt* indices_device,
                                    HYPRE_Real* vector_buffer_device,
                                    std::vector<HYPRE_Real>& continuous_vector_values,
                                    const std::vector<int>& local_hypre_to_local_dune,
                                    bool owner_first) {
#if HAVE_CUDA
        if constexpr (is_gpu_type<VectorType>::value) {
            // GPU input type
            transferGpuVectorToHypre(vec, hypre_vec, use_gpu_backend, indices, indices_device, vector_buffer_device, continuous_vector_values, local_hypre_to_local_dune, owner_first);
        } else
#endif
        {
            // CPU input type
            transferCpuVectorToHypre(vec, hypre_vec, use_gpu_backend, indices, indices_device, vector_buffer_device, continuous_vector_values, local_hypre_to_local_dune, owner_first);
        }
    }

    /**
     * @brief Transfer vector from Hypre to any vector type (CPU or GPU)
     *
     * Selects the transfer method based on the vector type and backend configuration.
     *
     * @param hypre_vec Source Hypre vector
     * @param vec Destination vector (typically BlockVector (CPU) or GpuVector (GPU))
     * @param use_gpu_backend Whether the Hypre backend is using GPU acceleration
     * @param indices Pre-allocated host array for owned DOF indices
     * @param indices_device Pre-allocated device array for owned DOF indices (for GPU backend)
     * @param vector_buffer_device Pre-allocated device buffer for vector values (for GPU backend)
     * @param continuous_vector_values Host buffer for continuous vector values
     * @param local_hypre_to_local_dune Mapping from HYPRE indices to Dune indices (used for non-owner-first)
     * @param owner_first Whether owned DOFs come first in the ordering
     * @throws HypreError if the transfer fails
     */
    template<typename VectorType>
    static void transferVectorFromHypre(HYPRE_IJVector hypre_vec, VectorType& vec, bool use_gpu_backend,
                                      const std::vector<HYPRE_BigInt>& indices,
                                      HYPRE_BigInt* indices_device,
                                      HYPRE_Real* vector_buffer_device,
                                      std::vector<HYPRE_Real>& continuous_vector_values,
                                      const std::vector<int>& local_hypre_to_local_dune,
                                      bool owner_first) {
#if HAVE_CUDA
        if constexpr (is_gpu_type<VectorType>::value) {
            // GPU input type
            transferHypreToGpuVector(hypre_vec, vec, use_gpu_backend, indices, indices_device, vector_buffer_device, continuous_vector_values, local_hypre_to_local_dune, owner_first);
        } else
#endif
        {
            // CPU input type
            transferHypreToCpuVector(hypre_vec, vec, use_gpu_backend, indices, indices_device, vector_buffer_device, continuous_vector_values, local_hypre_to_local_dune, owner_first);
        }
    }

    //--------------------------------------------------------------------------
    // Matrix Operations
    //--------------------------------------------------------------------------



    /**
     * @brief Initialize Hypre matrix from any matrix type (CPU or GPU)
     *
     * Selects the transfer method based on the matrix type and backend configuration.
     * Uses HYPRE_IJMatrixSetValues2 with pre-computed row_indexes for optimal performance.
     *
     * @param matrix Source matrix (typically BCRSMatrix (CPU) or GpuSparseMatrix (GPU))
     * @param hypre_matrix Destination Hypre matrix
     * @param use_gpu_backend Whether the Hypre backend is using GPU acceleration
     * @param ncols Pre-allocated host array for number of columns per row
     * @param rows Pre-allocated host array for row indices
     * @param cols Pre-allocated host array for column indices
     * @param row_indexes Pre-computed row indexes for HYPRE_IJMatrixSetValues2
     * @param ncols_device Pre-allocated device array for number of columns per row (for GPU backend)
     * @param rows_device Pre-allocated device array for row indices (for GPU backend)
     * @param cols_device Pre-allocated device array for column indices (for GPU backend)
     * @param row_indexes_device Pre-allocated device array for row indexes (for GPU backend)
     * @param matrix_buffer_device Pre-allocated device buffer for matrix values (for GPU backend)
     * @throws HypreError if initialization fails
     */
    template<typename MatrixType>
    static void initializeMatrix(const MatrixType& matrix, HYPRE_IJMatrix hypre_matrix, bool use_gpu_backend,
                                const std::vector<HYPRE_Int>& ncols,
                                const std::vector<HYPRE_BigInt>& rows,
                                const std::vector<HYPRE_BigInt>& cols,
                                const std::vector<HYPRE_Int>& row_indexes,
                                HYPRE_Int* ncols_device,
                                HYPRE_BigInt* rows_device,
                                HYPRE_BigInt* cols_device,
                                HYPRE_Int* row_indexes_device,
                                HYPRE_Real* matrix_buffer_device) {
#if HAVE_CUDA
        if constexpr (is_gpu_type<MatrixType>::value) {
            // GPU input type
            initializeMatrixFromGpuSparseMatrix(matrix, hypre_matrix, use_gpu_backend,
                                              ncols, rows, cols, row_indexes,
                                              ncols_device, rows_device, cols_device, row_indexes_device);
        } else
#endif
        {
            // CPU input type
            initializeMatrixFromCpuMatrix(matrix, hypre_matrix, use_gpu_backend,
                                        ncols, rows, cols, row_indexes,
                                        ncols_device, rows_device, cols_device, row_indexes_device,
                                        matrix_buffer_device);
        }
    }



    /**
     * @brief Update matrix values in Hypre
     *
     * Selects the update method based on the matrix type and backend configuration.
     * Updates the coefficients of a Hypre matrix, without changing its sparsity pattern.
     * Uses HYPRE_IJMatrixSetValues2 with pre-computed row_indexes for optimal performance.
     *
     * @param matrix Source matrix with updated values
     * @param hypre_matrix Hypre matrix to update
     * @param use_gpu_backend Whether the Hypre backend is using GPU acceleration
     * @param ncols Pre-allocated host array for number of columns per row
     * @param rows Pre-allocated host array for row indices
     * @param cols Pre-allocated host array for column indices
     * @param row_indexes Pre-computed row indexes for HYPRE_IJMatrixSetValues2
     * @param ncols_device Pre-allocated device array for number of columns per row (for GPU backend)
     * @param rows_device Pre-allocated device array for row indices (for GPU backend)
     * @param cols_device Pre-allocated device array for column indices (for GPU backend)
     * @param row_indexes_device Pre-allocated device array for row indexes (for GPU backend)
     * @param matrix_buffer_device Pre-allocated device buffer for matrix values (for GPU backend)
     * @throws HypreError if the update fails
     */
    template<typename MatrixType>
    static void updateMatrixValues(const MatrixType& matrix, HYPRE_IJMatrix hypre_matrix, bool use_gpu_backend,
                                 const std::vector<HYPRE_Int>& ncols,
                                 const std::vector<HYPRE_BigInt>& rows,
                                 const std::vector<HYPRE_BigInt>& cols,
                                 const std::vector<HYPRE_Int>& row_indexes,
                                 HYPRE_Int* ncols_device,
                                 HYPRE_BigInt* rows_device,
                                 HYPRE_BigInt* cols_device,
                                 HYPRE_Int* row_indexes_device,
                                 HYPRE_Real* matrix_buffer_device) {
#if HAVE_CUDA
        if constexpr (is_gpu_type<MatrixType>::value) {
            // GPU input type
            updateMatrixFromGpuSparseMatrix(matrix, hypre_matrix, use_gpu_backend,
                                          ncols, rows, cols, row_indexes,
                                          ncols_device, rows_device, cols_device, row_indexes_device);
        } else
#endif
        {
            // CPU input type
            updateMatrixFromCpuMatrix(matrix, hypre_matrix, use_gpu_backend,
                                    ncols, rows, cols, row_indexes,
                                    ncols_device, rows_device, cols_device, row_indexes_device,
                                    matrix_buffer_device);
        }
    }

    //--------------------------------------------------------------------------
    // Matrix Value Retrieval Operations
    //--------------------------------------------------------------------------

        /**
     * @brief Get matrix values from Hypre matrix
     *
     * Retrieves the values from a Hypre matrix for verification or debugging purposes.
     *
     * @param hypre_matrix The Hypre matrix to read values from
     * @param ncols Array containing number of columns per row
     * @param rows Array containing row indices
     * @param cols Array containing column indices
     * @param use_gpu_backend Whether the Hypre backend is using GPU acceleration
     * @return Vector containing the retrieved matrix values
     * @throws HypreError if the operation fails
     * @note For GPU backend, uses HYPRE_ParCSRMatrixGetRow since HYPRE_IJMatrixGetValues not supported on GPU
     */
    static std::vector<HYPRE_Real> getMatrixValues(HYPRE_IJMatrix hypre_matrix,
                                                  const std::vector<HYPRE_Int>& ncols,
                                                  const std::vector<HYPRE_BigInt>& rows,
                                                  const std::vector<HYPRE_BigInt>& cols,
                                                  bool use_gpu_backend = false) {
        const auto N = rows.size();
        std::vector<HYPRE_Real> values(cols.size());

        if (use_gpu_backend) {
            // HYPRE_IJMatrixGetValues not supported on GPU, use HYPRE_ParCSRMatrixGetRow instead
            void* object;
            OPM_HYPRE_SAFE_CALL(HYPRE_IJMatrixGetObject(hypre_matrix, &object));
            auto par_csr_matrix = static_cast<HYPRE_ParCSRMatrix>(object);
            HYPRE_Int max_row_size = *std::max_element(ncols.begin(), ncols.end());
            HYPRE_Complex* row_values_host = hypre_CTAlloc(HYPRE_Complex, max_row_size, HYPRE_MEMORY_HOST);

            size_t value_idx = 0;
            for (size_t row_idx = 0; row_idx < N; ++row_idx) {
                HYPRE_Int row_ncols;
                HYPRE_BigInt* row_cols_device;
                HYPRE_Complex* row_values_device;
                OPM_HYPRE_SAFE_CALL(HYPRE_ParCSRMatrixGetRow(par_csr_matrix, rows[row_idx], &row_ncols, &row_cols_device, &row_values_device));
                if (row_ncols > 0) {
                    hypre_TMemcpy(row_values_host, row_values_device, HYPRE_Complex, row_ncols, HYPRE_MEMORY_HOST, HYPRE_MEMORY_DEVICE);
                    for (HYPRE_Int j = 0; j < row_ncols; ++j, ++value_idx) {
                        values[value_idx] = static_cast<HYPRE_Real>(row_values_host[j]);
                    }
                }
                OPM_HYPRE_SAFE_CALL(HYPRE_ParCSRMatrixRestoreRow(par_csr_matrix, rows[row_idx], &row_ncols, &row_cols_device, &row_values_device));
            }
            hypre_TFree(row_values_host, HYPRE_MEMORY_HOST);
        } else {
            OPM_HYPRE_SAFE_CALL(HYPRE_IJMatrixGetValues(hypre_matrix, N,
                                                        const_cast<HYPRE_Int*>(ncols.data()),
                                                        const_cast<HYPRE_BigInt*>(rows.data()),
                                                        const_cast<HYPRE_BigInt*>(cols.data()),
                                                        values.data()));
        }
        return values;
    }

    /**
     * @brief Compute row indexes for HYPRE_IJMatrixSetValues2
     *
     * Chooses the appropriate row index computation method based on owner_first and matrix type.
     * For owner_first=true: Uses simple prefix sum (contiguous data)
     * For owner_first=false: Uses mapping-based approach (handles gaps in data)
     *
     * @param matrix The matrix to analyze
     * @param ncols Array containing number of columns per row
     * @param local_dune_to_local_hypre Mapping from Dune indices to Hypre indices
     * @param owner_first Whether owned DOFs come first in the ordering
     * @return Vector containing row indexes for HYPRE_IJMatrixSetValues2
     */
    template<typename MatrixType>
    static std::vector<HYPRE_Int> computeRowIndexes(const MatrixType& matrix,
                                                   const std::vector<HYPRE_Int>& ncols,
                                                   const std::vector<int>& local_dune_to_local_hypre,
                                                   bool owner_first) {
        if (owner_first) {
            // Simple contiguous case: prefix sum of ncols
            std::vector<HYPRE_Int> row_indexes(ncols.size());
            row_indexes[0] = 0;
            for (size_t i = 1; i < ncols.size(); ++i) {
                row_indexes[i] = row_indexes[i-1] + ncols[i-1];
            }
            return row_indexes;
        } else {
            // We need to compute the row indexes with mapping since we have gaps in the data
#if HAVE_CUDA
            if constexpr (is_gpu_type<MatrixType>::value) {
                return computeRowIndexesWithMappingGpu(matrix, local_dune_to_local_hypre);
            } else
#endif
            {
                return computeRowIndexesWithMappingCpu(matrix, local_dune_to_local_hypre);
            }
        }
    }

    /**
     * @brief Compute row indexes for CPU matrix with ownership mapping
     *
     * Computes row_indexes that point directly into the FULL (including ghost) matrix,
     * allowing us to use the original matrix without any data copying.
     *
     * @param matrix Matrix to analyze
     * @param local_dune_to_local_hypre Mapping from Dune indices to Hypre indices (-1 for non-owned)
     * @return Vector containing row indexes that map to full matrix data positions
     */
    template<typename MatrixType>
    static std::vector<HYPRE_Int> computeRowIndexesWithMappingCpu(const MatrixType& matrix,
                                                               const std::vector<int>& local_dune_to_local_hypre) {
        const int N = std::count_if(local_dune_to_local_hypre.begin(), local_dune_to_local_hypre.end(),
                                   [](int val) { return val >= 0; });
        std::vector<HYPRE_Int> row_indexes(N);
        int cumulative_full_nnz = 0;  // Counter for entries in the FULL (including ghost) matrix
        // Iterate through matrix rows in order
        for (auto row = matrix.begin(); row != matrix.end(); ++row) {
            const int dune_row_idx = row.index();
            const int hypre_row_idx = local_dune_to_local_hypre[dune_row_idx];

            if (hypre_row_idx >= 0) {
                // This is an owned row - record where its data starts in the FULL (including ghost) matrix
                row_indexes[hypre_row_idx] = cumulative_full_nnz;
            }
            // Always advance the cumulative count (including non-owned rows - this creates gaps)
            cumulative_full_nnz += row->size();
        }
        return row_indexes;
    }

#if HAVE_CUDA
    /**
     * @brief Compute row indexes for GPU matrix with ownership mapping
     *
     * Specialized version for GPU matrices that computes row_indexes pointing
     * directly into the FULL (including ghost) GPU matrix data, allowing zero-copy operation.
     */
    template<typename T>
    static std::vector<HYPRE_Int> computeRowIndexesWithMappingGpu(const GpuSparseMatrix<T>& gpu_matrix,
                                                                 const std::vector<int>& local_dune_to_local_hypre) {
        const int N = std::count_if(local_dune_to_local_hypre.begin(), local_dune_to_local_hypre.end(),
                                   [](int val) { return val >= 0; });
        std::vector<HYPRE_Int> row_indexes(N);

        // Get row pointers from GPU matrix
        auto host_row_ptrs = gpu_matrix.getRowIndices().asStdVector();

        // Map each owned Hypre row to its starting position in the FULL (including ghost) GPU matrix
        for (int dune_row_idx = 0; dune_row_idx < static_cast<int>(gpu_matrix.N()); ++dune_row_idx) {
            const int hypre_row_idx = local_dune_to_local_hypre[dune_row_idx];

            if (hypre_row_idx >= 0) {
                // This is an owned row - record where its data starts in the FULL (including ghost) GPU matrix
                row_indexes[hypre_row_idx] = host_row_ptrs[dune_row_idx];
            }
            // Non-owned rows create natural gaps in the indexing
        }
        return row_indexes;
    }
#endif

private:
#if HAVE_CUDA
    //--------------------------------------------------------------------------
    // GPU Vector Operations (GpuVector input)
    //--------------------------------------------------------------------------

    /**
     * @brief Transfer GpuVector to Hypre vector
     */
    template<typename T>
    static void transferGpuVectorToHypre(const GpuVector<T>& gpu_vec, HYPRE_IJVector hypre_vec, bool use_gpu_backend,
                                       const std::vector<HYPRE_BigInt>& indices, HYPRE_BigInt* indices_device,
                                       HYPRE_Real* vector_buffer_device,
                                       std::vector<HYPRE_Real>& continuous_vector_values,
                                       const std::vector<int>& local_hypre_to_local_dune,
                                       bool owner_first) {
        const int N = static_cast<int>(indices.size());

        if (use_gpu_backend) {
#if HYPRE_USING_CUDA || HYPRE_USING_HIP
            // GPU backend with GPU input: use pre-allocated device arrays
            if (owner_first) {
                // Direct device-to-device transfer for owner-first ordering
                const T* device_ptr = gpu_vec.data();
                OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorSetValues(hypre_vec, N, indices_device, const_cast<T*>(device_ptr)));
            } else {
                // Use continuous storage and device buffer for non-owner-first ordering
                setContinuousGpuVectorForHypre(gpu_vec, continuous_vector_values, local_hypre_to_local_dune);
                hypre_TMemcpy(vector_buffer_device, continuous_vector_values.data(), HYPRE_Real, N, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
                OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorSetValues(hypre_vec, N, indices_device, vector_buffer_device));
            }
#endif
        } else {
            // CPU backend with GPU input: copy via host memory
            if (owner_first) {
                // Get values to host and then set
                auto host_values = gpu_vec.asStdVector();
                OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorSetValues(hypre_vec, N, const_cast<HYPRE_BigInt*>(indices.data()), reinterpret_cast<HYPRE_Real*>(host_values.data())));
            } else {
                // Use continuous storage for non-owner-first ordering
                setContinuousGpuVectorForHypre(gpu_vec, continuous_vector_values, local_hypre_to_local_dune);
                OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorSetValues(hypre_vec, N, const_cast<HYPRE_BigInt*>(indices.data()), continuous_vector_values.data()));
            }
        }
        OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorAssemble(hypre_vec));
    }

    /**
     * @brief Transfer Hypre vector to GpuVector
     */
    template<typename T>
    static void transferHypreToGpuVector(HYPRE_IJVector hypre_vec, GpuVector<T>& gpu_vec, bool use_gpu_backend,
                                       const std::vector<HYPRE_BigInt>& indices, HYPRE_BigInt* indices_device,
                                       HYPRE_Real* vector_buffer_device,
                                       std::vector<HYPRE_Real>& continuous_vector_values,
                                       const std::vector<int>& local_hypre_to_local_dune,
                                       bool owner_first) {
        const int N = static_cast<int>(indices.size());

        if (use_gpu_backend) {
#if HYPRE_USING_CUDA || HYPRE_USING_HIP
            // GPU backend with GPU input: use pre-allocated device arrays
            if (owner_first) {
                // Direct device-to-device transfer for owner-first ordering
                T* device_ptr = gpu_vec.data();
                OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorGetValues(hypre_vec, N, indices_device, device_ptr));
            } else {
                // Use device buffer and then remap for non-owner-first ordering
                OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorGetValues(hypre_vec, N, indices_device, vector_buffer_device));
                hypre_TMemcpy(continuous_vector_values.data(), vector_buffer_device, HYPRE_Real, N, HYPRE_MEMORY_HOST, HYPRE_MEMORY_DEVICE);
                setGpuVectorFromContinuousVector(gpu_vec, continuous_vector_values, local_hypre_to_local_dune);
            }
#endif
        } else {
            // CPU backend with GPU input: copy via host memory
            if (owner_first) {
                // Get values to host and then copy to GPU
                auto host_values = gpu_vec.asStdVector();
                OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorGetValues(hypre_vec, N, const_cast<HYPRE_BigInt*>(indices.data()), reinterpret_cast<HYPRE_Real*>(host_values.data())));
                gpu_vec = GpuVector<T>(host_values);
            } else {
                // Use continuous storage for non-owner-first ordering
                OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorGetValues(hypre_vec, N, const_cast<HYPRE_BigInt*>(indices.data()), continuous_vector_values.data()));
                setGpuVectorFromContinuousVector(gpu_vec, continuous_vector_values, local_hypre_to_local_dune);
            }
        }
    }

#endif // HAVE_CUDA

    //--------------------------------------------------------------------------
    // GPU Vector Helper Functions
    //--------------------------------------------------------------------------

#if HAVE_CUDA
    template<typename T>
    static void setContinuousGpuVectorForHypre(const GpuVector<T>& v, std::vector<HYPRE_Real>& continuous_vector_values, const std::vector<int>& local_hypre_to_local_dune)
    {
        // Get vector data to host first
        auto host_values = v.asStdVector();
        // Set values using the mapping
        for (size_t i = 0; i < local_hypre_to_local_dune.size(); ++i) {
            continuous_vector_values[i] = host_values[local_hypre_to_local_dune[i]];
        }
    }

    template<typename T>
    static void setGpuVectorFromContinuousVector(GpuVector<T>& v, const std::vector<HYPRE_Real>& continuous_vector_values, const std::vector<int>& local_hypre_to_local_dune)
    {
        // Copy values to host and update values with mapping
        auto host_values = v.asStdVector();
        for (size_t i = 0; i < local_hypre_to_local_dune.size(); ++i) {
            host_values[local_hypre_to_local_dune[i]] = continuous_vector_values[i];
        }
        // Copy back to GPU
        v = GpuVector<T>(host_values);
    }
#endif

    //--------------------------------------------------------------------------
    // CPU Vector Operations (BlockVector input)
    //--------------------------------------------------------------------------

    /**
     * @brief Extract owned vector values in the order expected by HYPRE
     */
    template<typename VectorType>
    static void setContinuousVectorForHypre(const VectorType& v, std::vector<HYPRE_Real>& continuous_vector_values, const std::vector<int>& local_hypre_to_local_dune)
    {
        // Set v values solution vectors
        for (size_t i = 0; i < local_hypre_to_local_dune.size(); ++i) {
            continuous_vector_values[i] = v[local_hypre_to_local_dune[i]][0];
        }
    }

    /**
     * @brief Distribute HYPRE vector values back to original vector positions
     */
    template<typename VectorType>
    static void setDuneVectorFromContinuousVector(VectorType& v, const std::vector<HYPRE_Real>& continuous_vector_values, const std::vector<int>& local_hypre_to_local_dune)
    {
        // Place values back in original positions
        for (size_t i = 0; i < local_hypre_to_local_dune.size(); ++i) {
            v[local_hypre_to_local_dune[i]][0] = continuous_vector_values[i];
        }
    }

    /**
     * @brief Transfer CPU vector to Hypre vector
     */
    template<typename VectorType>
    static void transferCpuVectorToHypre(const VectorType& cpu_vec, HYPRE_IJVector hypre_vec, bool use_gpu_backend,
                                       const std::vector<HYPRE_BigInt>& indices,
                                       HYPRE_BigInt* indices_device,
                                       HYPRE_Real* vector_buffer_device,
                                       std::vector<HYPRE_Real>& continuous_vector_values,
                                       const std::vector<int>& local_hypre_to_local_dune,
                                       bool owner_first) {
        const int N = static_cast<int>(indices.size());
        using T = typename VectorType::field_type;

        if (use_gpu_backend) {
#if HYPRE_USING_CUDA || HYPRE_USING_HIP
            // GPU backend with CPU input: use pre-allocated device arrays
            if (owner_first) {
                // Direct host-to-device transfer for owner-first ordering
                const T* values = &(cpu_vec[0][0]);
                hypre_TMemcpy(vector_buffer_device, values, HYPRE_Real, N, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
                OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorSetValues(hypre_vec, N, indices_device, vector_buffer_device));
            } else {
                // Use continuous storage and device buffer for non-owner-first ordering
                setContinuousVectorForHypre(cpu_vec, continuous_vector_values, local_hypre_to_local_dune);
                hypre_TMemcpy(vector_buffer_device, continuous_vector_values.data(), HYPRE_Real, N, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
                OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorSetValues(hypre_vec, N, indices_device, vector_buffer_device));
            }
#endif
        } else {
            // CPU backend with CPU input: use host arrays directly
            if (owner_first) {
                // Direct transfer for owner-first ordering
                const T* values = &(cpu_vec[0][0]);
                OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorSetValues(hypre_vec, N, const_cast<HYPRE_BigInt*>(indices.data()), values));
            } else {
                // Use continuous storage for non-owner-first ordering
                setContinuousVectorForHypre(cpu_vec, continuous_vector_values, local_hypre_to_local_dune);
                OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorSetValues(hypre_vec, N, const_cast<HYPRE_BigInt*>(indices.data()), continuous_vector_values.data()));
            }
        }

        OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorAssemble(hypre_vec));
    }

    /**
     * @brief Transfer Hypre vector to CPU vector
     */
    template<typename VectorType>
    static void transferHypreToCpuVector(HYPRE_IJVector hypre_vec, VectorType& cpu_vec, bool use_gpu_backend,
                                       const std::vector<HYPRE_BigInt>& indices,
                                       HYPRE_BigInt* indices_device,
                                       HYPRE_Real* vector_buffer_device,
                                       std::vector<HYPRE_Real>& continuous_vector_values,
                                       const std::vector<int>& local_hypre_to_local_dune,
                                       bool owner_first) {
        const int N = static_cast<int>(indices.size());
        using T = typename VectorType::field_type;

        if (use_gpu_backend) {
#if HYPRE_USING_CUDA || HYPRE_USING_HIP
            // GPU backend with CPU input: use pre-allocated device arrays
            if (owner_first) {
                // Direct device-to-host transfer for owner-first ordering
                T* values = &(cpu_vec[0][0]);
                OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorGetValues(hypre_vec, N, indices_device, vector_buffer_device));
                hypre_TMemcpy(values, vector_buffer_device, HYPRE_Real, N, HYPRE_MEMORY_HOST, HYPRE_MEMORY_DEVICE);
            } else {
                // Use device buffer and then remap for non-owner-first ordering
                OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorGetValues(hypre_vec, N, indices_device, vector_buffer_device));
                hypre_TMemcpy(continuous_vector_values.data(), vector_buffer_device, HYPRE_Real, N, HYPRE_MEMORY_HOST, HYPRE_MEMORY_DEVICE);
                setDuneVectorFromContinuousVector(cpu_vec, continuous_vector_values, local_hypre_to_local_dune);
            }
#endif
        } else {
            // CPU backend with CPU input: use host arrays directly
            if (owner_first) {
                // Direct transfer for owner-first ordering
                T* values = &(cpu_vec[0][0]);
                OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorGetValues(hypre_vec, N, const_cast<HYPRE_BigInt*>(indices.data()), values));
            } else {
                // Use continuous storage for non-owner-first ordering
                OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorGetValues(hypre_vec, N, const_cast<HYPRE_BigInt*>(indices.data()), continuous_vector_values.data()));
                setDuneVectorFromContinuousVector(cpu_vec, continuous_vector_values, local_hypre_to_local_dune);
            }
        }
    }

#if HAVE_CUDA
    //--------------------------------------------------------------------------
    // GPU Matrix Operations (GpuSparseMatrix input)
    //--------------------------------------------------------------------------

    /**
     * @brief Initialize Hypre matrix from GpuSparseMatrix using pre-allocated helper arrays
     * Uses HYPRE_IJMatrixSetValues2 with pre-computed row_indexes, which allows us to use the original GPU matrix data (with potential ghost values) directly.
     */
    template<typename T>
    static void initializeMatrixFromGpuSparseMatrix(const GpuSparseMatrix<T>& gpu_matrix, HYPRE_IJMatrix hypre_matrix, bool use_gpu_backend,
                                                   const std::vector<HYPRE_Int>& ncols,
                                                   const std::vector<HYPRE_BigInt>& rows,
                                                   const std::vector<HYPRE_BigInt>& cols,
                                                   const std::vector<HYPRE_Int>& row_indexes,
                                                   HYPRE_Int* ncols_device,
                                                   HYPRE_BigInt* rows_device,
                                                   HYPRE_BigInt* cols_device,
                                                   HYPRE_Int* row_indexes_device) {
        const auto N = rows.size();
        const T* values = gpu_matrix.getNonZeroValues().data();

        if (use_gpu_backend) {
            // GPU backend with GPU input: use pre-allocated device arrays
#if HYPRE_USING_CUDA || HYPRE_USING_HIP
            // Direct device-to-device transfer using smart row_indexes
            OPM_HYPRE_SAFE_CALL(HYPRE_IJMatrixSetValues2(hypre_matrix, N, ncols_device, rows_device, row_indexes_device, cols_device, values));
#endif
        } else {
            // CPU backend with GPU input: copy to host first
            auto host_values = gpu_matrix.getNonZeroValues().asStdVector();
            OPM_HYPRE_SAFE_CALL(HYPRE_IJMatrixSetValues2(hypre_matrix, N, const_cast<HYPRE_Int*>(ncols.data()), const_cast<HYPRE_BigInt*>(rows.data()), const_cast<HYPRE_Int*>(row_indexes.data()), const_cast<HYPRE_BigInt*>(cols.data()), host_values.data()));
        }
        OPM_HYPRE_SAFE_CALL(HYPRE_IJMatrixAssemble(hypre_matrix));
    }

    /**
     * @brief Update Hypre matrix values from GpuSparseMatrix using pre-allocated helper arrays
     * Uses HYPRE_IJMatrixSetValues2 with pre-computed row_indexes, which allows us to use the original GPU matrix data (with potential ghost values) directly.
     */
    template<typename T>
    static void updateMatrixFromGpuSparseMatrix(const GpuSparseMatrix<T>& gpu_matrix, HYPRE_IJMatrix hypre_matrix, bool use_gpu_backend,
                                               const std::vector<HYPRE_Int>& ncols,
                                               const std::vector<HYPRE_BigInt>& rows,
                                               const std::vector<HYPRE_BigInt>& cols,
                                               const std::vector<HYPRE_Int>& row_indexes,
                                               HYPRE_Int* ncols_device,
                                               HYPRE_BigInt* rows_device,
                                               HYPRE_BigInt* cols_device,
                                               HYPRE_Int* row_indexes_device) {
        // For value updates, we need to re-upload the full matrix since Hypre doesn't have
        // a separate coefficient update API
        initializeMatrixFromGpuSparseMatrix(gpu_matrix, hypre_matrix, use_gpu_backend,
                                          ncols, rows, cols, row_indexes,
                                          ncols_device, rows_device, cols_device, row_indexes_device);
    }

#endif // HAVE_CUDA

    /**
     * @brief Initialize Hypre matrix from CPU matrix
     * Uses HYPRE_IJMatrixSetValues2 with pre-computed row_indexes, which allows us to use the original CPU matrix data (with potential ghost values) directly.
     */
    template<typename MatrixType>
    static void initializeMatrixFromCpuMatrix(const MatrixType& cpu_matrix, HYPRE_IJMatrix hypre_matrix, bool use_gpu_backend,
                                            const std::vector<HYPRE_Int>& ncols,
                                            const std::vector<HYPRE_BigInt>& rows,
                                            const std::vector<HYPRE_BigInt>& cols,
                                            const std::vector<HYPRE_Int>& row_indexes,
                                            HYPRE_Int* ncols_device,
                                            HYPRE_BigInt* rows_device,
                                            HYPRE_BigInt* cols_device,
                                            HYPRE_Int* row_indexes_device,
                                            HYPRE_Real* matrix_buffer_device) {
        const auto N = rows.size();
        const auto nnz = cpu_matrix.nonzeroes();  // Total entries including ghost

        using T = typename MatrixType::field_type;
        const T* values = &(cpu_matrix[0][0][0][0]);

        if (use_gpu_backend) {
#if HYPRE_USING_CUDA || HYPRE_USING_HIP
            hypre_TMemcpy(matrix_buffer_device, values, HYPRE_Real, nnz, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
            OPM_HYPRE_SAFE_CALL(HYPRE_IJMatrixSetValues2(hypre_matrix, N, ncols_device, rows_device, row_indexes_device, cols_device, matrix_buffer_device));
#endif
        } else {
            OPM_HYPRE_SAFE_CALL(HYPRE_IJMatrixSetValues2(hypre_matrix, N, const_cast<HYPRE_Int*>(ncols.data()), const_cast<HYPRE_BigInt*>(rows.data()), const_cast<HYPRE_Int*>(row_indexes.data()), const_cast<HYPRE_BigInt*>(cols.data()), values));
        }

        OPM_HYPRE_SAFE_CALL(HYPRE_IJMatrixAssemble(hypre_matrix));
    }

    /**
     * @brief Update Hypre matrix values from CPU matrix
     * Uses HYPRE_IJMatrixSetValues2 with pre-computed row_indexes for optimal performance.
     */
    template<typename MatrixType>
    static void updateMatrixFromCpuMatrix(const MatrixType& cpu_matrix, HYPRE_IJMatrix hypre_matrix, bool use_gpu_backend,
                                        const std::vector<HYPRE_Int>& ncols,
                                        const std::vector<HYPRE_BigInt>& rows,
                                        const std::vector<HYPRE_BigInt>& cols,
                                        const std::vector<HYPRE_Int>& row_indexes,
                                        HYPRE_Int* ncols_device,
                                        HYPRE_BigInt* rows_device,
                                        HYPRE_BigInt* cols_device,
                                        HYPRE_Int* row_indexes_device,
                                        HYPRE_Real* matrix_buffer_device) {
        // For value updates, we need to re-upload the full matrix since Hypre doesn't have
        // a separate coefficient update API like AMGX
        initializeMatrixFromCpuMatrix(cpu_matrix, hypre_matrix, use_gpu_backend,
                                    ncols, rows, cols, row_indexes,
                                    ncols_device, rows_device, cols_device, row_indexes_device,
                                    matrix_buffer_device);
    }
};

} // namespace Opm::gpuistl

#endif // OPM_HYPRE_INTERFACE_HPP