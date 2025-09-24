/*
  Copyright 2024 SINTEF AS
  Copyright 2024-2025 Equinor ASA

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
#include <HYPRE_krylov.h>
#include <HYPRE_parcsr_ls.h>
#include <_hypre_utilities.h>

#include <numeric>
#include <vector>

namespace Hypre
{

namespace HypreInterface = Opm::gpuistl::HypreInterface;

/**
 * @brief Wrapper for Hypre's BoomerAMG preconditioner.
 *
 * This class provides an interface to the BoomerAMG preconditioner from the Hypre library.
 * It is designed to work with matrices, update vectors, and defect vectors specified by the template parameters.
 * The HypreInterface class provides a unified interface to Hypre's functionality, allowing for easy
 * switching between CPU and GPU input data types and backend acceleration.
 *
 * Supports four use cases:
 * 1. Input type is CPU and backend acceleration is CPU
 * 2. Input type is CPU and backend acceleration is GPU
 * 3. Input type is GPU and backend acceleration is GPU
 * 4. Input type is GPU and backend acceleration is CPU
 *
 * @tparam M The matrix type
 * @tparam X The vector type for the solution
 * @tparam Y The vector type for the right-hand side
 */
template <class M, class X, class Y, class Comm>
class HyprePreconditioner : public Dune::PreconditionerWithUpdate<X, Y>
{
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
     * @param comm Parallel communicator.
     */
    HyprePreconditioner(const M& A, const Opm::PropertyTree prm, const Comm& comm)
        : A_(A)
        , comm_(comm)
    {
        OPM_TIMEBLOCK(prec_construct);
        int size;
        int rank;
        MPI_Comm mpi_comm;
        if constexpr (std::is_same_v<Comm, Dune::Amg::SequentialInformation>) {
            mpi_comm = MPI_COMM_SELF;
        } else {
            mpi_comm = comm.communicator();
        }
        MPI_Comm_size(mpi_comm, &size);
        MPI_Comm_rank(mpi_comm, &rank);
        if (size > 1) {
            assert(size == comm.communicator().size());
            assert(rank == comm.communicator().rank());
        }
        // Set use_gpu_backend_ to user value if specified, otherwise match input type
        use_gpu_backend_ = prm.get<bool>("use_gpu", Opm::gpuistl::is_gpu_type<M>::value);

        // Initialize Hypre library with backend configuration
        HypreInterface::initialize(use_gpu_backend_);

        // Create solver
        solver_ = HypreInterface::createAMGSolver();
        HypreInterface::setSolverParameters(solver_, prm, use_gpu_backend_);

        // Setup parallel info and mappings
        par_info_ = HypreInterface::setupHypreParallelInfo(comm_, A_);

        // Setup sparsity pattern
        sparsity_pattern_ = HypreInterface::setupSparsityPattern(A_, par_info_, par_info_.owner_first);

        // Setup host arrays
        host_arrays_.row_indexes = HypreInterface::computeRowIndexes(
            A_, sparsity_pattern_.ncols, par_info_.local_dune_to_local_hypre, par_info_.owner_first);

        // Create indices for vector operations - simple sequential indices for owned DOFs
        host_arrays_.indices.resize(par_info_.N_owned);
        std::iota(host_arrays_.indices.begin(), host_arrays_.indices.end(), par_info_.dof_offset);

        // Setup continuous vector values buffer - only needed for non-owner-first
        if (!par_info_.owner_first) {
            host_arrays_.continuous_vector_values.resize(par_info_.N_owned);
        }

        // Allocate device arrays if using GPU backend
        if (use_gpu_backend_) {
#if HYPRE_USING_CUDA || HYPRE_USING_HIP
            device_arrays_.ncols_device = hypre_CTAlloc(HYPRE_Int, par_info_.N_owned, HYPRE_MEMORY_DEVICE);
            device_arrays_.rows_device = hypre_CTAlloc(HYPRE_BigInt, par_info_.N_owned, HYPRE_MEMORY_DEVICE);
            device_arrays_.cols_device = hypre_CTAlloc(HYPRE_BigInt, sparsity_pattern_.nnz, HYPRE_MEMORY_DEVICE);
            device_arrays_.row_indexes_device = hypre_CTAlloc(HYPRE_Int, par_info_.N_owned, HYPRE_MEMORY_DEVICE);
            device_arrays_.indices_device = hypre_CTAlloc(HYPRE_BigInt, par_info_.N_owned, HYPRE_MEMORY_DEVICE);
            device_arrays_.vector_buffer_device = hypre_CTAlloc(HYPRE_Real, par_info_.N_owned, HYPRE_MEMORY_DEVICE);
            if constexpr (!Opm::gpuistl::is_gpu_type<M>::value) {
                // For CPU input and GPU backend we need to allocate space for transfering the matrix values
                // Note that the buffer must be allocated with the number of nonzeroes in the matrix, not the
                // sparsity_pattern.nnz because we need to copy the entire matrix values from the host to the device.
                device_arrays_.matrix_buffer_device = hypre_CTAlloc(HYPRE_Real, A_.nonzeroes(), HYPRE_MEMORY_DEVICE);
            }
            // Copy data to device
            hypre_TMemcpy(device_arrays_.ncols_device, sparsity_pattern_.ncols.data(), HYPRE_Int, par_info_.N_owned, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
            hypre_TMemcpy(device_arrays_.rows_device, sparsity_pattern_.rows.data(), HYPRE_BigInt, par_info_.N_owned, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
            hypre_TMemcpy(device_arrays_.cols_device, sparsity_pattern_.cols.data(), HYPRE_BigInt, sparsity_pattern_.nnz, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
            hypre_TMemcpy(device_arrays_.row_indexes_device, host_arrays_.row_indexes.data(), HYPRE_Int, par_info_.N_owned, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
            hypre_TMemcpy(device_arrays_.indices_device, host_arrays_.indices.data(), HYPRE_BigInt, par_info_.N_owned, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
#endif
        }

        // Create Hypre matrix and vectors
        A_hypre_ = HypreInterface::createMatrix(par_info_.N_owned, par_info_.dof_offset, comm_);
        x_hypre_ = HypreInterface::createVector(par_info_.N_owned, par_info_.dof_offset, comm_);
        b_hypre_ = HypreInterface::createVector(par_info_.N_owned, par_info_.dof_offset, comm_);

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
        if (use_gpu_backend_) {
#if HYPRE_USING_CUDA || HYPRE_USING_HIP
            if (device_arrays_.ncols_device) {
                hypre_TFree(device_arrays_.ncols_device, HYPRE_MEMORY_DEVICE);
            }
            if (device_arrays_.rows_device) {
                hypre_TFree(device_arrays_.rows_device, HYPRE_MEMORY_DEVICE);
            }
            if (device_arrays_.cols_device) {
                hypre_TFree(device_arrays_.cols_device, HYPRE_MEMORY_DEVICE);
            }
            if (device_arrays_.row_indexes_device) {
                hypre_TFree(device_arrays_.row_indexes_device, HYPRE_MEMORY_DEVICE);
            }
            if (device_arrays_.indices_device) {
                hypre_TFree(device_arrays_.indices_device, HYPRE_MEMORY_DEVICE);
            }
            if (device_arrays_.vector_buffer_device) {
                hypre_TFree(device_arrays_.vector_buffer_device, HYPRE_MEMORY_DEVICE);
            }
            if (device_arrays_.matrix_buffer_device) {
                hypre_TFree(device_arrays_.matrix_buffer_device, HYPRE_MEMORY_DEVICE);
            }
#endif
        }

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
    void update() override
    {
        OPM_TIMEBLOCK(prec_update);

        // Update matrix values using pre-allocated helper arrays
        HypreInterface::updateMatrixValues(
            A_, A_hypre_, sparsity_pattern_, host_arrays_, device_arrays_, use_gpu_backend_);

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
    void pre(X& v, Y& /*d*/) override
    {
        comm_.copyOwnerToAll(v, v); // From dune: make dirichlet values consistent ??
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
    void apply(X& v, const Y& d) override
    {
        OPM_TIMEBLOCK(prec_apply);

        // Transfer vectors to Hypre
        HypreInterface::transferVectorToHypre(v, x_hypre_, host_arrays_, device_arrays_, par_info_, use_gpu_backend_);
        HypreInterface::transferVectorToHypre(d, b_hypre_, host_arrays_, device_arrays_, par_info_, use_gpu_backend_);

        // Get the underlying ParCSR matrix and ParVector objects
        HYPRE_ParCSRMatrix parcsr_A;
        HYPRE_ParVector par_x, par_b;
        HYPRE_SAFE_CALL(HYPRE_IJMatrixGetObject(A_hypre_, reinterpret_cast<void**>(&parcsr_A)));
        HYPRE_SAFE_CALL(HYPRE_IJVectorGetObject(x_hypre_, reinterpret_cast<void**>(&par_x)));
        HYPRE_SAFE_CALL(HYPRE_IJVectorGetObject(b_hypre_, reinterpret_cast<void**>(&par_b)));

        // Apply the preconditioner (one AMG V-cycle)
        HYPRE_SAFE_CALL(HYPRE_BoomerAMGSolve(solver_, parcsr_A, par_b, par_x));

        // Transfer result back
        HypreInterface::transferVectorFromHypre(x_hypre_, v, host_arrays_, device_arrays_, par_info_, use_gpu_backend_);
        // NB do we need to sync values to get correct values since a preconditioner
        // consistent result (a operator apply should give unique).
        comm_.copyOwnerToAll(v, v);
    }

    /**
     * @brief Post-processing step after applying the preconditioner.
     *
     * This method is currently a no-op.
     *
     * @param v The update vector.
     */
    void post(X& /*v*/) override
    {
    }

    /**
     * @brief Returns the solver category.
     *
     * @return The solver category, which is sequential.
     */
    Dune::SolverCategory::Category category() const override
    {
        return std::is_same_v<Comm, Dune::Amg::SequentialInformation> ? Dune::SolverCategory::sequential
                                                                      : Dune::SolverCategory::overlapping;
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
    const Comm& comm_; //!< The communication object for parallel operations.

    // Parallel information and sparsity pattern from HypreInterface
    HypreInterface::ParallelInfo par_info_;
    HypreInterface::SparsityPattern sparsity_pattern_;
    HypreInterface::HostArrays host_arrays_;
    HypreInterface::DeviceArrays device_arrays_;

    // Hypre handles
    HYPRE_Solver solver_ = nullptr;
    HYPRE_IJMatrix A_hypre_ = nullptr;
    HYPRE_IJVector x_hypre_ = nullptr;
    HYPRE_IJVector b_hypre_ = nullptr;

    // Backend configuration
    bool use_gpu_backend_ = false;
};

} // namespace Hypre

#endif
