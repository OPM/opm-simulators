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

#ifndef OPM_AMGX_INTERFACE_HPP
#define OPM_AMGX_INTERFACE_HPP

#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/gpuistl/GpuSparseMatrix.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_type_detection.hpp>

#include <amgx_c.h>
#include <cuda_runtime.h>

#include <fmt/core.h>

#include <string>
#include <vector>

namespace Opm::gpuistl
{

/**
 * @brief Exception class for AMGX errors
 */
class AmgxError : public std::runtime_error
{
public:
    explicit AmgxError(const std::string& msg)
        : std::runtime_error(msg)
    {
    }
};

/**
 * @brief Get a descriptive error message for an AMGX error code
 *
 * @param err The AMGX error code
 * @param expression The AMGX expression that caused the error
 * @param file The source file where the error occurred
 * @param function The function where the error occurred
 * @param line The line number where the error occurred
 * @return std::string The formatted error message
 */
inline std::string
getAmgxErrorMessage(
    AMGX_RC err, const std::string& expression, const std::string& file, const std::string& function, int line)
{
    char amgx_err_msg[4096];
    AMGX_get_error_string(err, amgx_err_msg, sizeof(amgx_err_msg));

    return fmt::format(
        "AMGX error in expression: {}\nError code: {}\nError message: {}\nLocation: {}:{} in function {}",
        expression,
        static_cast<int>(err),
        amgx_err_msg,
        file,
        line,
        function);
}

/**
 * @brief Safe call wrapper for AMGX functions
 *
 * Checks the return code from AMGX functions and throws an AmgxError if an error occurred.
 *
 * @param rc The AMGX return code to check
 * @param expression The expression being evaluated (for error reporting)
 * @param file The source file (typically __FILE__)
 * @param function The function name (typically __func__)
 * @param line The line number (typically __LINE__)
 * @throws AmgxError if the return code indicates an error
 */
inline void
amgxSafeCall(AMGX_RC rc, const std::string& expression, const std::string& file, const std::string& function, int line)
{
    if (rc != AMGX_RC_OK) {
        throw AmgxError(getAmgxErrorMessage(rc, expression, file, function, line));
    }
}

/**
 * @brief Macro to wrap AMGX function calls with error checking
 *
 * Example usage:
 * @code
 * OPM_AMGX_SAFE_CALL(AMGX_initialize());
 * @endcode
 */
#define OPM_AMGX_SAFE_CALL(expr) ::Opm::gpuistl::amgxSafeCall((expr), #expr, __FILE__, __func__, __LINE__)


/**
 * @brief Unified interface for AMGX operations with both CPU and GPU data structures
 *
 * This class provides utilities for working with AMGX resources and transferring
 * data between CPU/GPU data structures and AMGX handles. It handles type detection
 * and automatically chooses the most efficient transfer method:
 * - For GPU types: Device-to-device transfers
 * - For CPU types: Host-to-device transfers with memory pinning
 *
 * @note Error handling: All methods throw AmgxError exceptions when AMGX operations fail.
 */
class AmgxInterface
{
public:
    /**
     * @brief Initialize the AMGX library
     *
     * This should be called once at the start of the program before using any AMGX functionality.
     * @throws AmgxError if initialization fails
     */
    static void initialize()
    {
        OPM_AMGX_SAFE_CALL(AMGX_initialize());
    }

    /**
     * @brief Finalize the AMGX library
     *
     * This should be called once at the end of the program to release AMGX resources.
     * @throws AmgxError if finalization fails
     */
    static void finalize()
    {
        OPM_AMGX_SAFE_CALL(AMGX_finalize());
    }

    /**
     * @brief Create an AMGX config handle from a configuration string
     *
     * @param config_string Configuration string for AMGX
     * @return AMGX_config_handle The created config handle
     * @throws AmgxError if config creation fails
     */
    static AMGX_config_handle createConfig(const std::string& config_string)
    {
        AMGX_config_handle config;
        OPM_AMGX_SAFE_CALL(AMGX_config_create(&config, config_string.c_str()));
        return config;
    }

    /**
     * @brief Create AMGX resources from a config
     *
     * @param config The AMGX config handle
     * @return AMGX_resources_handle The created resources handle
     * @throws AmgxError if resource creation fails
     */
    static AMGX_resources_handle createResources(AMGX_config_handle config)
    {
        AMGX_resources_handle resources;
        OPM_AMGX_SAFE_CALL(AMGX_resources_create_simple(&resources, config));
        return resources;
    }

    /**
     * @brief Create an AMGX solver
     *
     * @param resources AMGX resources handle
     * @param mode AMGX mode (precision configuration)
     * @param config AMGX config handle
     * @return AMGX_solver_handle The created solver handle
     * @throws AmgxError if solver creation fails
     */
    static AMGX_solver_handle createSolver(AMGX_resources_handle resources, AMGX_Mode mode, AMGX_config_handle config)
    {
        AMGX_solver_handle solver;
        OPM_AMGX_SAFE_CALL(AMGX_solver_create(&solver, resources, mode, config));
        return solver;
    }

    /**
     * @brief Create an AMGX matrix
     *
     * @param resources AMGX resources handle
     * @param mode AMGX mode (precision configuration)
     * @return AMGX_matrix_handle The created matrix handle
     * @throws AmgxError if matrix creation fails
     */
    static AMGX_matrix_handle createMatrix(AMGX_resources_handle resources, AMGX_Mode mode)
    {
        AMGX_matrix_handle matrix;
        OPM_AMGX_SAFE_CALL(AMGX_matrix_create(&matrix, resources, mode));
        return matrix;
    }

    /**
     * @brief Create an AMGX vector
     *
     * @param resources AMGX resources handle
     * @param mode AMGX mode (precision configuration)
     * @return AMGX_vector_handle The created vector handle
     * @throws AmgxError if vector creation fails
     */
    static AMGX_vector_handle createVector(AMGX_resources_handle resources, AMGX_Mode mode)
    {
        AMGX_vector_handle vector;
        OPM_AMGX_SAFE_CALL(AMGX_vector_create(&vector, resources, mode));
        return vector;
    }

    /**
     * @brief Destroy an AMGX config handle
     *
     * @param config The config handle to destroy
     * @throws AmgxError if config destruction fails
     */
    static void destroyConfig(AMGX_config_handle config)
    {
        if (config) {
            OPM_AMGX_SAFE_CALL(AMGX_config_destroy(config));
        }
    }

    /**
     * @brief Destroy an AMGX resources handle
     *
     * @param resources The resources handle to destroy
     * @throws AmgxError if resource destruction fails
     *
     * @note There is a known issue where destroying resources and then
     * reinitializing AMGX can cause crashes in some versions of the library.
     * If you encounter this issue, consider keeping the resources alive for
     * the entire program lifetime.
     */
    static void destroyResources(AMGX_resources_handle resources)
    {
        if (resources) {
            OPM_AMGX_SAFE_CALL(AMGX_resources_destroy(resources));
        }
    }

    /**
     * @brief Destroy an AMGX solver handle
     *
     * @param solver The solver handle to destroy
     * @throws AmgxError if solver destruction fails
     */
    static void destroySolver(AMGX_solver_handle solver)
    {
        if (solver) {
            OPM_AMGX_SAFE_CALL(AMGX_solver_destroy(solver));
        }
    }

    /**
     * @brief Destroy an AMGX matrix handle
     *
     * @param matrix The matrix handle to destroy
     * @throws AmgxError if matrix destruction fails
     */
    template <typename MatrixType>
    static void destroyMatrix(AMGX_matrix_handle amgx_matrix, const MatrixType& matrix)
    {
        if constexpr (!is_gpu_type<MatrixType>::value) {
            // Unpin memory for CPU matrices
            using T = typename MatrixType::field_type;
            const T* values = &(matrix[0][0][0][0]);
            OPM_AMGX_SAFE_CALL(AMGX_unpin_memory(const_cast<T*>(values)));
        }

        if (amgx_matrix) {
            OPM_AMGX_SAFE_CALL(AMGX_matrix_destroy(amgx_matrix));
        }
    }

    /**
     * @brief Destroy an AMGX vector handle
     *
     * @param vector The vector handle to destroy
     * @throws AmgxError if vector destruction fails
     */
    static void destroyVector(AMGX_vector_handle vector)
    {
        if (vector) {
            OPM_AMGX_SAFE_CALL(AMGX_vector_destroy(vector));
        }
    }

    /**
     * @brief Update an AMGX vector from a GpuVector (device-to-device transfer)
     *
     * Updates the AMGX vector with the contents of the GpuVector using direct device memory access.
     * The AMGX vector must already be created with appropriate resources and mode.
     *
     * @param gpu_vec The source GpuVector
     * @param amgx_vec The AMGX vector to update
     * @throws AmgxError if the transfer fails or sizes mismatch
     */
    template <typename T>
    static void updateAmgxFromGpuVector(const GpuVector<T>& gpu_vec, AMGX_vector_handle amgx_vec)
    {
        // Get vector size from AMGX to verify compatibility
        int n, block_dim;
        OPM_AMGX_SAFE_CALL(AMGX_vector_get_size(amgx_vec, &n, &block_dim));

        if (n > 0 && static_cast<size_t>(n * block_dim) != gpu_vec.dim()) {
            throw AmgxError(fmt::format("Vector size mismatch in updateAmgxFromGpuVector: "
                                        "AMGX vector size {} vs. GpuVector size {}",
                                        n * block_dim,
                                        gpu_vec.dim()));
        }

        // Get raw device pointer directly from GpuVector
        const T* device_ptr = gpu_vec.data();

        // Update AMGX vector directly with the device pointer
        OPM_AMGX_SAFE_CALL(AMGX_vector_upload(amgx_vec, n, block_dim, const_cast<T*>(device_ptr)));
    }

    /**
     * @brief Update a GpuVector from an AMGX vector (device-to-device transfer)
     *
     * Updates the GpuVector with the contents of the AMGX vector using direct device memory access.
     *
     * @param amgx_vec The source AMGX vector
     * @param gpu_vec The GpuVector to update
     * @throws AmgxError if the transfer fails or sizes mismatch
     */
    template <typename T>
    static void updateGpuVectorFromAmgx(AMGX_vector_handle amgx_vec, GpuVector<T>& gpu_vec)
    {
        // Get vector size from AMGX to verify compatibility
        int n, block_dim;
        OPM_AMGX_SAFE_CALL(AMGX_vector_get_size(amgx_vec, &n, &block_dim));

        if (static_cast<size_t>(n * block_dim) != gpu_vec.dim()) {
            throw AmgxError(fmt::format("Vector size mismatch in updateGpuVectorFromAmgx: "
                                        "AMGX vector size {} vs. GpuVector size {}",
                                        n * block_dim,
                                        gpu_vec.dim()));
        }

        // Get destination device pointer from GpuVector
        T* dst_device_ptr = gpu_vec.data();

        // Download data directly from AMGX vector to GpuVector's device memory
        OPM_AMGX_SAFE_CALL(AMGX_vector_download(amgx_vec, dst_device_ptr));
    }

    /**
     * @brief Transfer vector to AMGX from any vector type (CPU or GPU)
     *
     * Selects the transfer method based on the vector type.
     *
     * @param vec Source vector (typically BlockVector (CPU) or GpuVector (GPU))
     * @param amgx_vec Destination AMGX vector
     * @throws AmgxError if the transfer fails
     */
    template <typename VectorType>
    static void transferVectorToAmgx(const VectorType& vec, AMGX_vector_handle amgx_vec)
    {
        if constexpr (is_gpu_type<VectorType>::value) {
            // GPU implementation - device-to-device transfer
            updateAmgxFromGpuVector(vec, amgx_vec);
        } else {
            // CPU implementation - host-to-device transfer
            // For BlockVector, get dimensions
            const int N = vec.size();
            const int block_size = 1; // Assumes scalar vectors for CPU case

            // Upload vector to AMGX
            OPM_AMGX_SAFE_CALL(AMGX_vector_upload(amgx_vec, N, block_size, &vec[0][0]));
        }
    }

    /**
     * @brief Transfer vector from AMGX to any vector type (CPU or GPU)
     *
     * Selects the transfer method based on the vector type.
     *
     * @param amgx_vec Source AMGX vector
     * @param vec Destination vector (typically BlockVector (CPU) or GpuVector (GPU))
     * @throws AmgxError if the transfer fails
     */
    template <typename VectorType>
    static void transferVectorFromAmgx(AMGX_vector_handle amgx_vec, VectorType& vec)
    {
        if constexpr (is_gpu_type<VectorType>::value) {
            // GPU implementation - device-to-device transfer
            updateGpuVectorFromAmgx(amgx_vec, vec);
        } else {
            // CPU implementation - device-to-host transfer
            OPM_AMGX_SAFE_CALL(AMGX_vector_download(amgx_vec, &vec[0][0]));
        }
    }

    /**
     * @brief Update an AMGX matrix from a GpuSparseMatrix (device-to-device transfer)
     *
     * Uploads the entire matrix structure and values from GpuSparseMatrix to AMGX.
     *
     * @param gpuSparseMatrix The source GpuSparseMatrix
     * @param amgxMatrix The AMGX matrix to update
     * @throws AmgxError if the transfer fails
     */
    template <typename T>
    static void updateAmgxMatrixFromGpuSparseMatrix(const GpuSparseMatrix<T>& gpuSparseMatrix,
                                                    AMGX_matrix_handle amgxMatrix)
    {
        // Get matrix dimensions and sparsity information
        auto n = detail::to_int(gpuSparseMatrix.N());
        auto nnz = detail::to_int(gpuSparseMatrix.nonzeroes());
        auto block_size = detail::to_int(gpuSparseMatrix.blockSize());

        // Get device pointers directly from GpuSparseMatrix
        const T* values = gpuSparseMatrix.getNonZeroValues().data();
        const int* row_ptrs = gpuSparseMatrix.getRowIndices().data();
        const int* col_indices = gpuSparseMatrix.getColumnIndices().data();

        // Update AMGX matrix with the device pointers
        OPM_AMGX_SAFE_CALL(
            AMGX_matrix_upload_all(amgxMatrix, n, nnz, block_size, block_size, row_ptrs, col_indices, values, nullptr));
    }

    /**
     * @brief Update only the coefficient values of an AMGX matrix from a GpuSparseMatrix
     *
     * Updates just the coefficient values of an AMGX matrix without changing its sparsity pattern.
     * This is more efficient when the matrix structure remains the same.
     *
     * @param gpuSparseMatrix The source GpuSparseMatrix with updated values
     * @param amgxMatrix The AMGX matrix to update
     * @throws AmgxError if the update fails
     */
    template <typename T>
    static void updateAmgxMatrixCoefficientsFromGpuSparseMatrix(const GpuSparseMatrix<T>& gpuSparseMatrix,
                                                                AMGX_matrix_handle amgxMatrix)
    {
        // Get matrix dimensions and sparsity information
        auto n = detail::to_int(gpuSparseMatrix.N());
        auto nnz = detail::to_int(gpuSparseMatrix.nonzeroes());

        // Get device pointer to values directly from GpuSparseMatrix
        const T* values = gpuSparseMatrix.getNonZeroValues().data();

        // Update coefficients in AMGX matrix with the device pointer
        OPM_AMGX_SAFE_CALL(AMGX_matrix_replace_coefficients(amgxMatrix, n, nnz, values, nullptr));
    }

    /**
     * @brief Update a GpuSparseMatrix from an AMGX matrix (device-to-device transfer)
     *
     * Downloads the matrix data from AMGX and updates the GpuSparseMatrix with the values.
     * The sparsity pattern is assumed to be identical, so only values are updated.
     *
     * @note This function is intended primarily for testing purposes and should not be used in performance-critical code paths,
     *       as it may involve unnecessary device-to-device transfers and temporary allocations.
     *
     * @param amgxMatrix The source AMGX matrix
     * @param gpuSparseMatrix The GpuSparseMatrix to update
     * @throws AmgxError if the transfer fails
     */
    template <typename T>
    static void updateGpuSparseMatrixFromAmgxMatrix(AMGX_matrix_handle amgxMatrix, GpuSparseMatrix<T>& gpuSparseMatrix)
    {
        // Get matrix dimensions from AMGX
        int n, nnz, block_sizex, block_sizey;
        OPM_AMGX_SAFE_CALL(AMGX_matrix_get_size(amgxMatrix, &n, &block_sizex, &block_sizey));
        OPM_AMGX_SAFE_CALL(AMGX_matrix_get_nnz(amgxMatrix, &nnz));

        // Allocate temporary device memory for row pointers and column indices (required by AMGX API)
        int* temp_row_ptrs;
        int* temp_col_indices;
        OPM_GPU_SAFE_CALL(cudaMalloc(&temp_row_ptrs, (n + 1) * sizeof(int)));
        OPM_GPU_SAFE_CALL(cudaMalloc(&temp_col_indices, nnz * sizeof(int)));

        // Get device pointer to values in the GpuSparseMatrix
        T* gpu_values = gpuSparseMatrix.getNonZeroValues().data();

        // AMGX requires a valid pointer for diagonal data (even if unused)
        void* diag_data_ptr = nullptr;

        // Download matrix values directly from AMGX to GpuSparseMatrix
        OPM_AMGX_SAFE_CALL(AMGX_matrix_download_all(amgxMatrix,
                                                    temp_row_ptrs,
                                                    temp_col_indices,
                                                    gpu_values,
                                                    &diag_data_ptr));

        // Clean up temporary device memory
        OPM_GPU_SAFE_CALL(cudaFree(temp_row_ptrs));
        OPM_GPU_SAFE_CALL(cudaFree(temp_col_indices));
    }

    /**
     * @brief Initialize an AMGX matrix from any matrix type (CPU or GPU)
     *
     * Selects the transfer method based on the matrix type.
     *
     * @param matrix Source matrix (typically BCRSMatrix (CPU) or GpuSparseMatrix (GPU))
     * @param amgx_matrix Destination AMGX matrix
     * @throws AmgxError if initialization fails
     */
    template <typename MatrixType>
    static void initializeMatrix(const MatrixType& matrix, AMGX_matrix_handle amgx_matrix)
    {
        if constexpr (is_gpu_type<MatrixType>::value) {
            // GPU implementation - device-to-device transfer
            updateAmgxMatrixFromGpuSparseMatrix(matrix, amgx_matrix);
        } else {
            // CPU implementation - host-to-device transfer
            auto N = detail::to_int(matrix.N());
            auto nnz = detail::to_int(matrix.nonzeroes());
            const int block_size = 1; // Assumes scalar matrices for CPU case

            // Setup sparsity pattern
            std::vector<int> row_ptrs(N + 1);
            std::vector<int> col_indices(nnz);

            int pos = 0;
            row_ptrs[0] = 0;
            for (auto row = matrix.begin(); row != matrix.end(); ++row) {
                for (auto col = row->begin(); col != row->end(); ++col) {
                    col_indices[pos++] = col.index();
                }
                row_ptrs[row.index() + 1] = pos;
            }

            // Get values pointer
            using T = typename MatrixType::field_type;
            const T* values = &(matrix[0][0][0][0]);
            // Indexing explanation:
            // matrix[0]             - First row of the matrix
            //          [0]          - First block in that row
            //            [0]        - First row within the 1x1 block
            //              [0]      - First column within the 1x1 block

            // Pin memory for CPU-to-GPU transfers
            OPM_AMGX_SAFE_CALL(AMGX_pin_memory(const_cast<T*>(values), sizeof(T) * nnz * block_size * block_size));

            // Upload to AMGX
            OPM_AMGX_SAFE_CALL(AMGX_matrix_upload_all(amgx_matrix,
                                                      N,
                                                      nnz,
                                                      block_size,
                                                      block_size,
                                                      row_ptrs.data(),
                                                      col_indices.data(),
                                                      const_cast<T*>(values),
                                                      nullptr));
        }
    }

    /**
     * @brief Initialize an AMGX vector with zeros
     *
     * Creates a zero-initialized vector of the appropriate size.
     *
     * @param N Size of the vector
     * @param block_size Block size (typically 1 for scalar vectors)
     * @param amgx_vector The AMGX vector to initialize
     * @throws AmgxError if initialization fails
     */
    static void initializeVector(int N, int block_size, AMGX_vector_handle amgx_vector)
    {
        OPM_AMGX_SAFE_CALL(AMGX_vector_set_zero(amgx_vector, N, block_size));
    }

    /**
     * @brief Update matrix values in AMGX
     *
     * Selects the update method based on the matrix type.
     * Updates the coefficients of an AMGX matrix, without changing its sparsity pattern.
     *
     * @param matrix Source matrix with updated values
     * @param amgx_matrix AMGX matrix to update
     * @throws AmgxError if the update fails
     */
    template <typename MatrixType>
    static void updateMatrixValues(const MatrixType& matrix, AMGX_matrix_handle amgx_matrix)
    {
        if constexpr (is_gpu_type<MatrixType>::value) {
            // GPU implementation - device-to-device transfer
            updateAmgxMatrixCoefficientsFromGpuSparseMatrix(matrix, amgx_matrix);
        } else {
            // CPU implementation - host-to-device transfer
            using T = typename MatrixType::field_type;
            const T* values = &(matrix[0][0][0][0]);
            OPM_AMGX_SAFE_CALL(AMGX_matrix_replace_coefficients(
                amgx_matrix, matrix.N(), matrix.nonzeroes(), const_cast<T*>(values), nullptr));
        }
    }

    /**
     * @brief Determine the appropriate AMGX mode based on matrix and vector field types
     *
     * @tparam MatrixFieldType The field type of the matrix
     * @tparam VectorFieldType The field type of the vector
     * @return AMGX_Mode The appropriate AMGX mode
     * @throws std::runtime_error if the type combination is not supported
     */
    template <typename MatrixFieldType, typename VectorFieldType>
    static AMGX_Mode determineAmgxMode()
    {
        if constexpr (std::is_same_v<MatrixFieldType, double> && std::is_same_v<VectorFieldType, double>) {
            return AMGX_mode_dDDI;
        } else if constexpr (std::is_same_v<MatrixFieldType, float> && std::is_same_v<VectorFieldType, double>) {
            return AMGX_mode_dDFI;
        } else if constexpr (std::is_same_v<MatrixFieldType, float> && std::is_same_v<VectorFieldType, float>) {
            return AMGX_mode_dFFI;
        } else {
            OPM_THROW(std::runtime_error,
                      "Unsupported combination of matrix and vector types for AMGX: "
                          + std::string(typeid(MatrixFieldType).name()) + " and "
                          + std::string(typeid(VectorFieldType).name()));
        }
    }
};

} // namespace Opm::gpuistl

#endif // OPM_AMGX_INTERFACE_HPP
