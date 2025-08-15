/*
  Copyright 2024 SINTEF AS
  Copyright 2024 Equinor ASA

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

#ifndef OPM_HYPRE_PRECONDITIONER_HEADER_INCLUDED
#define OPM_HYPRE_PRECONDITIONER_HEADER_INCLUDED

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/TimingMacros.hpp>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/gpuistl/HypreInterface.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_type_detection.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>
#include <HYPRE_krylov.h>
#include <_hypre_utilities.h>
#include <numeric>
#include <vector>

namespace Hypre {

using HypreInterface = Opm::gpuistl::HypreInterface;

/**
 * @brief Wrapper for Hypre's BoomerAMG preconditioner.
 *
 * This class provides an interface to the BoomerAMG preconditioner from the Hypre library.
 * It is designed to work with matrices, update vectors, and defect vectors specified by the template parameters.
 * The HypreInterface class provides a unified interface to Hypre's functionality, allowing for easy 
 * switching between CPU and GPU input data types and backend acceleration.
 *
 * Supports three use cases:
 * 1. Input type is CPU and backend acceleration is CPU
 * 2. Input type is CPU and backend acceleration is GPU
 * 3. Input type is GPU and backend acceleration is GPU
 *
 * @tparam M The matrix type
 * @tparam X The vector type for the solution
 * @tparam Y The vector type for the right-hand side
 */
template<class M, class X, class Y>
class HyprePreconditioner : public Dune::PreconditionerWithUpdate<X,Y> {
public:
    //! \brief The matrix type the preconditioner is for
    using matrix_type = M;
    //! \brief The field type of the matrix
    using matrix_field_type = typename M::field_type;
    //! \brief The domain type of the preconditioner
    using domain_type = X;
    //! \brief The range type of the preconditioner
    using range_type = Y;
    //! \brief The field type of the vectors
    using vector_field_type = typename X::field_type;

    /**
     * @brief Constructor for the HyprePreconditioner class.
     *
     * Initializes the preconditioner with the given matrix and property tree.
     *
     * @param A The matrix for which the preconditioner is constructed.
     * @param prm The property tree containing configuration parameters.
     */
    HyprePreconditioner(const M& A, const Opm::PropertyTree prm)
        : A_(A)
    {
        OPM_TIMEBLOCK(prec_construct);

        int size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        if (size > 1) {
            OPM_THROW(std::runtime_error, "HyprePreconditioner is currently only implemented for sequential runs");
        }

        use_gpu_backend_ = prm.get<bool>("use_gpu", false);

        // Initialize Hypre library with backend configuration
        HypreInterface::initialize(use_gpu_backend_);

        // Create solver
        solver_ = HypreInterface::createSolver();
        HypreInterface::setSolverParameters(solver_, prm, use_gpu_backend_);

        // Cache matrix dimensions
        N_ = static_cast<HYPRE_Int>(A_.N());
        nnz_ = static_cast<HYPRE_Int>(A_.nonzeroes());

        // Setup helper arrays for both CPU and GPU matrices
        setupHelperArrays();

        // Create Hypre matrix and vectors
        A_hypre_ = HypreInterface::createMatrix(N_);
        x_hypre_ = HypreInterface::createVector(N_);
        b_hypre_ = HypreInterface::createVector(N_);

        // Initialize matrix structure and values using pre-allocated helper arrays
        HypreInterface::initializeMatrix(A_, A_hypre_, use_gpu_backend_,
                                       ncols_, rows_, cols_,
                                       ncols_device_, rows_device_, cols_device_,
                                       matrix_buffer_device_);

        // Perform initial update
        update();
    }

    /**
     * @brief Destructor for HyprePreconditioner
     *
     * Cleans up resources allocated by the preconditioner.
     */
    ~HyprePreconditioner()
    {
        // Clean up device arrays if allocated
        cleanupHelperArrays();

        HypreInterface::destroySolver(solver_);
        HypreInterface::destroyVector(x_hypre_);
        HypreInterface::destroyVector(b_hypre_);
        HypreInterface::destroyMatrix(A_hypre_);
    }

    /**
     * @brief Updates the preconditioner with the current matrix values.
     *
     * This method should be called whenever the matrix values change.
     */
    void update() override {
        OPM_TIMEBLOCK(prec_update);

        // Update matrix values using pre-allocated helper arrays
        HypreInterface::updateMatrixValues(A_, A_hypre_, use_gpu_backend_,
                                         ncols_, rows_, cols_,
                                         ncols_device_, rows_device_, cols_device_,
                                         matrix_buffer_device_);

        // Get the underlying ParCSR matrix for setup
        HYPRE_ParCSRMatrix parcsr_A;
        HYPRE_SAFE_CALL(HYPRE_IJMatrixGetObject(A_hypre_, reinterpret_cast<void**>(&parcsr_A)));

        // Get the underlying ParVector objects
        HYPRE_ParVector par_x, par_b;
        HYPRE_SAFE_CALL(HYPRE_IJVectorGetObject(x_hypre_, reinterpret_cast<void**>(&par_x)));
        HYPRE_SAFE_CALL(HYPRE_IJVectorGetObject(b_hypre_, reinterpret_cast<void**>(&par_b)));

        // Setup the solver
        HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetup(solver_, parcsr_A, par_b, par_x));
    }

    /**
     * @brief Pre-processing step before applying the preconditioner.
     *
     * This method is currently a no-op.
     *
     * @param v The update vector.
     * @param d The defect vector.
     */
    void pre(X& /*v*/, Y& /*d*/) override {
    }

    /**
     * @brief Applies the preconditioner to a vector.
     *
     * Performs one AMG V-cycle to solve the system.
     * Involves uploading vectors to Hypre, applying the preconditioner,
     * and transferring the result back to the vector.
     *
     * @param v The update vector.
     * @param d The defect vector.
     */
    void apply(X& v, const Y& d) override {
        OPM_TIMEBLOCK(prec_apply);

        // Transfer vectors to Hypre
        HypreInterface::transferVectorToHypre(v, x_hypre_, use_gpu_backend_, indices_, indices_device_, vector_buffer_device_);
        HypreInterface::transferVectorToHypre(d, b_hypre_, use_gpu_backend_, indices_, indices_device_, vector_buffer_device_);

        // Get the underlying ParCSR matrix and ParVector objects
        HYPRE_ParCSRMatrix parcsr_A;
        HYPRE_ParVector par_x, par_b;
        HYPRE_SAFE_CALL(HYPRE_IJMatrixGetObject(A_hypre_, reinterpret_cast<void**>(&parcsr_A)));
        HYPRE_SAFE_CALL(HYPRE_IJVectorGetObject(x_hypre_, reinterpret_cast<void**>(&par_x)));
        HYPRE_SAFE_CALL(HYPRE_IJVectorGetObject(b_hypre_, reinterpret_cast<void**>(&par_b)));

        // Apply the preconditioner (one AMG V-cycle)
        HYPRE_SAFE_CALL(HYPRE_BoomerAMGSolve(solver_, parcsr_A, par_b, par_x));

        // Transfer result back
        HypreInterface::transferVectorFromHypre(x_hypre_, v, use_gpu_backend_, indices_, indices_device_, vector_buffer_device_);
    }

    /**
     * @brief Post-processing step after applying the preconditioner.
     *
     * This method is currently a no-op.
     *
     * @param v The update vector.
     */
    void post(X& /*v*/) override {
    }

    /**
     * @brief Returns the solver category.
     *
     * @return The solver category, which is sequential.
     */
    Dune::SolverCategory::Category category() const override {
        return Dune::SolverCategory::sequential;
    }

    /**
     * @brief Checks if the preconditioner has a perfect update.
     *
     * @return True, indicating that the preconditioner can be perfectly updated.
     */
    bool hasPerfectUpdate() const override
    {
        // The Hypre preconditioner can depend on the values of the matrix so it does not have perfect update.
        // However, copying the matrix to Hypre requires to setup the solver again, so this is handled internally.
        // So for ISTLSolver, we can return true.
        return true;
    }

private:
    // Reference to the input matrix
    const M& A_;

    // Hypre handles
    HYPRE_Solver solver_ = nullptr;
    HYPRE_IJMatrix A_hypre_ = nullptr;
    HYPRE_IJVector x_hypre_ = nullptr;
    HYPRE_IJVector b_hypre_ = nullptr;

    // Backend configuration
    bool use_gpu_backend_ = false;

    // Matrix dimensions (cached for performance)
    HYPRE_Int N_;
    HYPRE_Int nnz_;

    // Host helper arrays for matrix operations
    std::vector<HYPRE_Int> ncols_;
    std::vector<HYPRE_BigInt> rows_;
    std::vector<HYPRE_BigInt> cols_;
    std::vector<HYPRE_BigInt> indices_;

    // Device helper arrays (only allocated when use_gpu_backend_ is true)
    HYPRE_Int* ncols_device_ = nullptr;
    HYPRE_BigInt* rows_device_ = nullptr;
    HYPRE_BigInt* cols_device_ = nullptr;
    HYPRE_BigInt* indices_device_ = nullptr;
    HYPRE_Real* vector_buffer_device_ = nullptr;
    HYPRE_Real* matrix_buffer_device_ = nullptr;

    /**
     * @brief Setup helper arrays for matrix and vector operations
     */
    void setupHelperArrays() {
        // Setup host arrays
        ncols_.resize(N_);
        rows_.resize(N_);
        cols_.resize(nnz_);
        indices_.resize(N_);

        // Setup sparsity pattern based on matrix type
#if HYPRE_USING_CUDA || HYPRE_USING_HIP
        if constexpr (Opm::gpuistl::is_gpu_type<M>::value) {
            setupSparsityPatternFromGpuMatrix();
        } else
#endif
        {
            setupSparsityPatternFromCpuMatrix();
        }

        // Setup indices array for vectors
        std::iota(indices_.begin(), indices_.end(), 0);

        // Allocate device arrays if using GPU backend
        if (use_gpu_backend_) {
#if HYPRE_USING_CUDA || HYPRE_USING_HIP
            ncols_device_ = hypre_CTAlloc(HYPRE_Int, N_, HYPRE_MEMORY_DEVICE);
            rows_device_ = hypre_CTAlloc(HYPRE_BigInt, N_, HYPRE_MEMORY_DEVICE);
            cols_device_ = hypre_CTAlloc(HYPRE_BigInt, nnz_, HYPRE_MEMORY_DEVICE);
            indices_device_ = hypre_CTAlloc(HYPRE_BigInt, N_, HYPRE_MEMORY_DEVICE);
            vector_buffer_device_ = hypre_CTAlloc(HYPRE_Real, N_, HYPRE_MEMORY_DEVICE);
            matrix_buffer_device_ = hypre_CTAlloc(HYPRE_Real, nnz_, HYPRE_MEMORY_DEVICE);
            // Copy data to device
            hypre_TMemcpy(ncols_device_, ncols_.data(), HYPRE_Int, N_, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
            hypre_TMemcpy(rows_device_, rows_.data(), HYPRE_BigInt, N_, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
            hypre_TMemcpy(cols_device_, cols_.data(), HYPRE_BigInt, nnz_, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
            hypre_TMemcpy(indices_device_, indices_.data(), HYPRE_BigInt, N_, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
#endif
        }
    }

    /**
     * @brief Setup sparsity pattern from CPU matrix (BCRSMatrix)
     */
    void setupSparsityPatternFromCpuMatrix() {
        int pos = 0;
        for (auto row = A_.begin(); row != A_.end(); ++row) {
            const int rowIdx = row.index();
            rows_[rowIdx] = rowIdx;
            ncols_[rowIdx] = row->size();

            for (auto col = row->begin(); col != row->end(); ++col) {
                cols_[pos++] = col.index();
            }
        }
    }

    #if HYPRE_USING_CUDA || HYPRE_USING_HIP
    /**
     * @brief Setup sparsity pattern from GPU matrix (GpuSparseMatrix)
     */
    void setupSparsityPatternFromGpuMatrix() {
        // Get row pointers from GPU matrix (one-time host copy during setup)
        auto host_row_ptrs = A_.getRowIndices().asStdVector();
        auto host_col_indices = A_.getColumnIndices().asStdVector();

        // Setup row information and ncols
        for (int i = 0; i < N_; ++i) {
            rows_[i] = i;
            ncols_[i] = host_row_ptrs[i + 1] - host_row_ptrs[i];
        }

        // Convert column indices to HYPRE_BigInt format
        for (int i = 0; i < nnz_; ++i) {
            cols_[i] = host_col_indices[i];
        }
    }
#endif

    /**
     * @brief Clean up device helper arrays
     */
    void cleanupHelperArrays() {
        if (use_gpu_backend_) {
#if HYPRE_USING_CUDA || HYPRE_USING_HIP
            if (ncols_device_) {
                hypre_TFree(ncols_device_, HYPRE_MEMORY_DEVICE);
                ncols_device_ = nullptr;
            }
            if (rows_device_) {
                hypre_TFree(rows_device_, HYPRE_MEMORY_DEVICE);
                rows_device_ = nullptr;
            }
            if (cols_device_) {
                hypre_TFree(cols_device_, HYPRE_MEMORY_DEVICE);
                cols_device_ = nullptr;
            }
            if (indices_device_) {
                hypre_TFree(indices_device_, HYPRE_MEMORY_DEVICE);
                indices_device_ = nullptr;
            }
            if (vector_buffer_device_) {
                hypre_TFree(vector_buffer_device_, HYPRE_MEMORY_DEVICE);
                vector_buffer_device_ = nullptr;
            }
            if (matrix_buffer_device_) {
                hypre_TFree(matrix_buffer_device_, HYPRE_MEMORY_DEVICE);
                matrix_buffer_device_ = nullptr;
            }
#endif
        }
    }
};

} // namespace Hypre

#endif
