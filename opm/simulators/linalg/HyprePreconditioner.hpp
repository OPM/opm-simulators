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
#include <dune/istl/repartition.hh>
#include <dune/istl/paamg/pinfo.hh>

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
template<class M, class X, class Y,class Comm>
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
    // template <typename Prm = std::enable_if_t<std::is_same_v<Comm, Dune::Amg::SequentialInformation>, ::Opm::PropertyTree>>
    // HyprePreconditioner (const M& A, const Prm& prm)
    //     : HyprePreconditioner(A, prm, Dune::Amg::SequentialInformation())
    // {
    //     //NB if this is used comm_ can never be used
    // }

    HyprePreconditioner(const M& A, const Opm::PropertyTree prm,const Comm& comm)
        : A_(A),comm_(comm)
    {
        OPM_TIMEBLOCK(prec_construct);
        int size;
        int rank;
        MPI_Comm mpi_comm;
        if constexpr (std::is_same_v<Comm,Dune::Amg::SequentialInformation>){
            mpi_comm = MPI_COMM_SELF;
        }else{
            mpi_comm = comm.communicator();
        }
        MPI_Comm_size(mpi_comm, &size);
        MPI_Comm_rank(mpi_comm, &rank);
        if (size > 1) {
            assert(size == comm.communicator().size());
            assert(rank == comm.communicator().rank());
        }

        use_gpu_backend_ = prm.get<bool>("use_gpu", false);

        // Initialize Hypre library with backend configuration
        HypreInterface::initialize(use_gpu_backend_);

        // Create solver
        solver_ = HypreInterface::createSolver();
        HypreInterface::setSolverParameters(solver_, prm, use_gpu_backend_);

        // Setup parallel info and mappings first
        setupHypreParallelInfo(comm_);
        // Setup helper arrays for matrix operations
        setupHelperArrays();

        // Create Hypre matrix and vectors
        A_hypre_ = HypreInterface::createMatrix(N_, dof_offset_, comm_);
        x_hypre_ = HypreInterface::createVector(N_, dof_offset_, comm_);
        b_hypre_ = HypreInterface::createVector(N_, dof_offset_, comm_);

        // Initialize matrix structure and values using pre-allocated helper arrays
        HypreInterface::initializeMatrix(A_, A_hypre_, use_gpu_backend_,
                                       ncols_, rows_, cols_,
                                       ncols_device_, rows_device_, cols_device_,
                                       matrix_buffer_device_,
                                       continuous_matrix_entries_,
                                       local_dune_to_local_hypre_,
                                       owner_first_);

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
                                         matrix_buffer_device_,
                                         continuous_matrix_entries_,
                                         local_dune_to_local_hypre_,
                                         owner_first_);

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
    void pre(X& v, Y& /*d*/) override {
      comm_.copyOwnerToAll(v,v); // From dune: make dirichlet values consistent ??
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
        HypreInterface::transferVectorToHypre(v, x_hypre_, use_gpu_backend_, indices_, indices_device_, vector_buffer_device_, continuous_vector_values_, local_hypre_to_local_dune_, owner_first_);
        HypreInterface::transferVectorToHypre(d, b_hypre_, use_gpu_backend_, indices_, indices_device_, vector_buffer_device_, continuous_vector_values_, local_hypre_to_local_dune_, owner_first_);

        // Get the underlying ParCSR matrix and ParVector objects
        HYPRE_ParCSRMatrix parcsr_A;
        HYPRE_ParVector par_x, par_b;
        HYPRE_SAFE_CALL(HYPRE_IJMatrixGetObject(A_hypre_, reinterpret_cast<void**>(&parcsr_A)));
        HYPRE_SAFE_CALL(HYPRE_IJVectorGetObject(x_hypre_, reinterpret_cast<void**>(&par_x)));
        HYPRE_SAFE_CALL(HYPRE_IJVectorGetObject(b_hypre_, reinterpret_cast<void**>(&par_b)));

        // Apply the preconditioner (one AMG V-cycle)
        HYPRE_SAFE_CALL(HYPRE_BoomerAMGSolve(solver_, parcsr_A, par_b, par_x));

        // Transfer result back
        HypreInterface::transferVectorFromHypre(x_hypre_, v, use_gpu_backend_, indices_, indices_device_, vector_buffer_device_, continuous_vector_values_, local_hypre_to_local_dune_, owner_first_);
        // NB do we need to sync values to get correct values since a preconditioner
        // consistent result (a operator apply should give unique).
        comm_.copyOwnerToAll(v,v);
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
        return std::is_same_v<Comm, Dune::Amg::SequentialInformation> ?
           Dune::SolverCategory::sequential : Dune::SolverCategory::overlapping;
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


    // mappings between dune and hypre use int to have negative values as not owned
    std::vector<int> local_dune_to_local_hypre_;
    std::vector<int> local_dune_to_global_hypre_;
    std::vector<int> local_hypre_to_local_dune_;// will only be needed if not owner first
    std::vector<HYPRE_Real> continuous_matrix_entries_;
    std::vector<HYPRE_Real> continuous_vector_values_;
    int dof_offset_;
    bool owner_first_;

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
    HYPRE_Real* values_device_ = nullptr; //!< Device array for matrix values.
    HYPRE_Real* x_values_device_ = nullptr; //!< Device array for solution vector values.
    HYPRE_Real* b_values_device_ = nullptr; //!< Device array for right-hand side vector values.
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

        // Create indices vector
        indices_.resize(N_);
        std::iota(indices_.begin(), indices_.end(), dof_offset_);

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
        // Setup arrays and fill column indices
        int pos = 0;
        for (auto row = A_.begin(); row != A_.end(); ++row) {
            const int rind = row.index();
            if (local_dune_to_local_hypre_[rind] < 0) {
                // This is a ghost dof, skip it
                continue;
            }

            const int local_rowIdx = local_dune_to_local_hypre_[rind];
            const int global_rowIdx = local_dune_to_global_hypre_[rind];

            rows_[local_rowIdx] = global_rowIdx;
            ncols_[local_rowIdx] = row->size();

            for (auto col = row->begin(); col != row->end(); ++col) {
                const int global_colIdx = local_dune_to_global_hypre_[col.index()];
                assert(global_colIdx >= 0);
                cols_[pos++] = global_colIdx;
            }
        }
    }

    #if HYPRE_USING_CUDA || HYPRE_USING_HIP
    /**
     * @brief Setup sparsity pattern from GPU matrix (GpuSparseMatrix)
     */
    void setupSparsityPatternFromGpuMatrix() {
        // Get row pointers and column indices from GPU matrix (one-time host copy during setup)
        auto host_row_ptrs = A_.getRowIndices().asStdVector();
        auto host_col_indices = A_.getColumnIndices().asStdVector();

        int pos = 0;
        for (int rind = 0; rind < static_cast<int>(A_.N()); ++rind) {
            if (local_dune_to_local_hypre_[rind] < 0) {
                // This is a ghost dof, skip it
                continue;
            }
            const int row_start = host_row_ptrs[rind];
            const int row_end = host_row_ptrs[rind + 1];
            const int num_cols = row_end - row_start;

            const int local_rowIdx = local_dune_to_local_hypre_[rind];
            const int global_rowIdx = local_dune_to_global_hypre_[rind];

            rows_[local_rowIdx] = global_rowIdx;
            ncols_[local_rowIdx] = num_cols;

            // Extract column indices for this row and map them to global Hypre indices
            for (int col_idx = row_start; col_idx < row_end; ++col_idx) {
                const int colIdx = host_col_indices[col_idx];
                const int global_colIdx = local_dune_to_global_hypre_[colIdx];
                assert(global_colIdx >= 0);
                cols_[pos++] = global_colIdx;
            }
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



    /**
     * @brief Set up helping structure for hypre in serial case
     */
    void setupHypreParallelInfo(const Dune::Amg::SequentialInformation& /*comm*/)
    {
        N_ = A_.N();
        nnz_ = A_.nonzeroes();

        local_dune_to_local_hypre_.resize(N_);
        local_dune_to_global_hypre_.resize(N_);
        local_hypre_to_local_dune_.resize(N_);
        for (int i = 0; i < N_; ++i) {
            local_dune_to_local_hypre_[i] = i;
            local_hypre_to_local_dune_[i] = i;
            local_dune_to_global_hypre_[i] = i;
        }
        dof_offset_ = 0;
        owner_first_ = true;
    }

    /**
     * @brief Set up helping structure for hypre in parallel case
     */

    template <class Commun>
    void setupHypreParallelInfo(const Commun& comm)
    {

        const auto& collective_comm = comm.communicator();
        // make global numbering and mapping for matrix

        size_t count = 0;
        local_dune_to_local_hypre_.resize(comm.indexSet().size(), -1);
        local_dune_to_global_hypre_.resize(comm.indexSet().size(), -1);

        // in case index set is not full dune we fix it
        if(!(A_.N() == comm.indexSet().size())){
            // in OPM this will likely not be trigged
            // ensure no holes in index sett
            const_cast<Commun&>(comm).buildGlobalLookup(A_.N());// need?
            Dune::Amg::MatrixGraph<M> graph(const_cast<M&>(A_));// do not know why not const ref is sufficient
            Dune::fillIndexSetHoles(graph, const_cast<Commun&>(comm));
            assert(A_.N() == comm.indexSet().size());
        }
        // NB iterations in index set is not in order of indexes
        // first find owners
        for (const auto& ind : comm.indexSet()) {
            int local_ind = ind.local().local();
            if (ind.local().attribute() == Dune::OwnerOverlapCopyAttributeSet::owner) {
                local_dune_to_local_hypre_[local_ind] = 1;//count;
                count += 1;

            } else {
               local_dune_to_local_hypre_[local_ind] = -1;
            }
        }
        N_ = count;

         // make order of variables
         // NB owner first will fail for cprw...
        bool owner_first = true;
        bool visited_copy = false;
        count = 0;
        for (size_t i=0; i < local_dune_to_local_hypre_.size(); ++i) {
            if (local_dune_to_local_hypre_[i] < 0) {
                visited_copy = true;
                assert(local_dune_to_local_hypre_[i] == -1);
            } else {
                local_dune_to_local_hypre_[i] = count;
                local_hypre_to_local_dune_.push_back(i);
                owner_first = owner_first && !visited_copy;
                count +=1;
            }
        }

         // owner first need other copying of data
        owner_first_ = owner_first;
        if (!owner_first){  // only need this storage if not owner first
            continuous_vector_values_.resize(N_);
        }

        // Gather the number of DOFs from each process
        std::vector<int> dof_counts_per_process(collective_comm.size());
        collective_comm.allgather(&N_, 1, dof_counts_per_process.data());

        // Calculate our process's offset (sum of DOFs in processes before us)
        dof_offset_ = std::accumulate(dof_counts_per_process.begin(),
                                      dof_counts_per_process.begin() + collective_comm.rank(), 0);

        // Create global DOF indices by adding offset to local indices
        for (size_t i = 0; i < local_dune_to_local_hypre_.size(); ++i) {
            if (local_dune_to_local_hypre_[i] >= 0) {
                local_dune_to_global_hypre_[i] = local_dune_to_local_hypre_[i] + dof_offset_;
            } else {
              local_dune_to_global_hypre_[i] = -1;
            }
         }

        if (collective_comm.rank() > 0) {
            assert(dof_offset_>0);
        }
        // communicate global numbering of dofs
        comm.copyOwnerToAll(local_dune_to_global_hypre_, local_dune_to_global_hypre_);

        // calculate matrix propeties which can be used to optimize memmory
        int nnz = 0;
        for (auto row = A_.begin(); row != A_.end(); ++row) {
            const int rowIdx = row.index();
            for (auto col = row->begin(); col != row->end(); ++col) {
                if (!(local_dune_to_local_hypre_[rowIdx] < 0) ){
                    nnz++;
                }
            }
        }
        nnz_ = nnz;
        if (!owner_first){
            continuous_matrix_entries_.resize(nnz_);
        }
    }
};

} // namespace Hypre

#endif
