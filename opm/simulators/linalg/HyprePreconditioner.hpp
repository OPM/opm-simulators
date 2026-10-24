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
#include <opm/simulators/linalg/hypreinterface/HypreInterface.hpp>
#include <opm/simulators/linalg/hypreinterface/HypreNumericRap.hpp>

#if HYPRE_USING_CUDA
#include <opm/simulators/linalg/gpuistl/detail/gpu_type_detection.hpp>
#elif HYPRE_USING_HIP
#include <opm/simulators/linalg/gpuistl_hip/detail/gpu_type_detection.hpp>
#endif

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <HYPRE.h>
#include <HYPRE_krylov.h>
#include <HYPRE_parcsr_ls.h>
#include <_hypre_parcsr_ls.h>
#include <_hypre_utilities.h>

#include <algorithm>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

namespace Opm::linalg
{

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
#if HYPRE_USING_CUDA || HYPRE_USING_HIP
        use_gpu_backend_ = prm.get<bool>("use_gpu", gpuistl::is_gpu_type<M>::value);
#endif

        // Initialize Hypre library with backend configuration
        HypreInterface::initialize(use_gpu_backend_);

        // Create solver
        solver_ = HypreInterface::createAMGSolver();
        HypreInterface::setSolverParameters(solver_, prm, use_gpu_backend_);
        // >1: keep the hierarchy for that many updates; 0: always re-setup cheaply and let the
        // linear solver decide when to recreate (no perfect update then)
        const int interval = prm.get<int>("reuse_interval", 1);
        external_reuse_ = interval == 0;
        reuse_interval_ = external_reuse_ ? std::numeric_limits<int>::max() : std::max(1, interval);
        // "frozen": keep the whole hierarchy; "galerkin": keep coarsening and P, recompute P^T A P
        const auto mode = prm.get<std::string>("reuse_mode", "frozen");
        if (mode != "frozen" && mode != "galerkin") {
            OPM_THROW(std::invalid_argument, "HyprePreconditioner: unknown reuse_mode " + mode);
        }
        reuse_galerkin_ = mode == "galerkin";
        numeric_rap_ = prm.get<bool>("numeric_rap", true);
        relax_type_ = prm.get<int>("relax_type", use_gpu_backend_ ? 16 : 13);

        // Setup parallel info and mappings
        par_info_ = HypreInterface::setupHypreParallelInfo(comm_, A_);

        // Setup sparsity pattern
        sparsity_pattern_ = HypreInterface::setupSparsityPattern(A_, par_info_, par_info_.owner_first);

        // Setup host arrays
        host_arrays_.row_indexes = HypreInterface::computeRowIndexes(A_,
                                                                     sparsity_pattern_.ncols,
                                                                     par_info_.local_dune_to_local_hypre,
                                                                     par_info_.owner_first);

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

        // Update matrix values: after the first update, a direct gather into the ParCSR arrays
        if (!direct_slots_.empty()) {
            OPM_TIMEBLOCK(prec_update_copy);
            copyValuesDirect_();
        } else {
            HypreInterface::updateMatrixValues(A_, A_hypre_, sparsity_pattern_,
                                               host_arrays_, device_arrays_, use_gpu_backend_);
            if (!use_gpu_backend_) {
                buildDirectSlots_();
            }
        }

        // Get the underlying ParCSR matrix for setup
        HYPRE_ParCSRMatrix parcsr_A;
        HYPRE_SAFE_CALL(HYPRE_IJMatrixGetObject(A_hypre_, reinterpret_cast<void**>(&parcsr_A)));

        // Get the underlying ParVector objects
        HYPRE_ParVector par_x, par_b;
        HYPRE_SAFE_CALL(HYPRE_IJVectorGetObject(x_hypre_, reinterpret_cast<void**>(&par_x)));
        HYPRE_SAFE_CALL(HYPRE_IJVectorGetObject(b_hypre_, reinterpret_cast<void**>(&par_b)));

        // The values were updated in place, so the fine level of a kept hierarchy already
        // sees the new matrix; only its l1 smoother norms are stale.
        if (is_setup_ && !use_gpu_backend_ && ++updates_since_setup_ < reuse_interval_) {
            OPM_TIMEBLOCK(prec_update_reuse);
            if (reuse_galerkin_) {
                resetupGalerkin_(parcsr_A);
            } else {
                refreshFineL1Norms_(parcsr_A);
            }
            return;
        }
        HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetup(solver_, parcsr_A, par_b, par_x));
        rap_maps_.clear();
        is_setup_ = true;
        updates_since_setup_ = 0;
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
        HypreInterface::transferVectorToHypre(v, x_hypre_, host_arrays_,
                                              device_arrays_, par_info_, use_gpu_backend_);
        HypreInterface::transferVectorToHypre(d, b_hypre_, host_arrays_,
                                              device_arrays_, par_info_, use_gpu_backend_);

        // Get the underlying ParCSR matrix and ParVector objects
        HYPRE_ParCSRMatrix parcsr_A;
        HYPRE_ParVector par_x, par_b;
        HYPRE_SAFE_CALL(HYPRE_IJMatrixGetObject(A_hypre_, reinterpret_cast<void**>(&parcsr_A)));
        HYPRE_SAFE_CALL(HYPRE_IJVectorGetObject(x_hypre_, reinterpret_cast<void**>(&par_x)));
        HYPRE_SAFE_CALL(HYPRE_IJVectorGetObject(b_hypre_, reinterpret_cast<void**>(&par_b)));

        // Apply the preconditioner (one AMG V-cycle)
        HYPRE_SAFE_CALL(HYPRE_BoomerAMGSolve(solver_, parcsr_A, par_b, par_x));

        // Transfer result back
        HypreInterface::transferVectorFromHypre(x_hypre_, v, host_arrays_,
                                                device_arrays_, par_info_, use_gpu_backend_);
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
        // Full setups are handled internally, unless the linear solver decides (reuse_interval 0).
        return !external_reuse_;
    }

private:
    // Where each Dune value (as passed to HYPRE_IJMatrixSetValues2) lives in the assembled
    // ParCSR matrix: diag slot h, or offd slot nnz(diag) + h.
    void buildDirectSlots_()
    {
        HYPRE_ParCSRMatrix parcsr_A;
        HYPRE_SAFE_CALL(HYPRE_IJMatrixGetObject(A_hypre_, reinterpret_cast<void**>(&parcsr_A)));
        auto* A = reinterpret_cast<hypre_ParCSRMatrix*>(parcsr_A);
        hypre_CSRMatrix* Ad = hypre_ParCSRMatrixDiag(A);
        hypre_CSRMatrix* Ao = hypre_ParCSRMatrixOffd(A);
        const HYPRE_Int nnzD = hypre_CSRMatrixNumNonzeros(Ad);
        const HYPRE_BigInt firstRow = hypre_ParCSRMatrixFirstRowIndex(A);
        const HYPRE_BigInt firstCol = hypre_ParCSRMatrixFirstColDiag(A);
        const HYPRE_Int nLocCols = hypre_CSRMatrixNumCols(Ad);
        const HYPRE_BigInt* cmap = hypre_ParCSRMatrixColMapOffd(A);
        const HYPRE_Int nOffCols = hypre_CSRMatrixNumCols(Ao);
        const auto& sp = sparsity_pattern_;
        std::vector<HYPRE_Int> src, dst;
        for (std::size_t i = 0; i < sp.rows.size(); ++i) {
            const HYPRE_Int r = static_cast<HYPRE_Int>(sp.rows[i] - firstRow);
            // values and columns share the offsets given to HYPRE_IJMatrixSetValues2
            const HYPRE_Int start = host_arrays_.row_indexes[i];
            for (HYPRE_Int t = 0; t < sp.ncols[i]; ++t) {
                const HYPRE_BigInt c = sp.cols[start + t];
                HYPRE_Int slot = -1;
                if (c >= firstCol && c < firstCol + nLocCols) {
                    const HYPRE_Int lc = static_cast<HYPRE_Int>(c - firstCol);
                    for (HYPRE_Int h = hypre_CSRMatrixI(Ad)[r]; h < hypre_CSRMatrixI(Ad)[r + 1]; ++h) {
                        if (hypre_CSRMatrixJ(Ad)[h] == lc) { slot = h; break; }
                    }
                } else if (nOffCols > 0) {
                    const HYPRE_BigInt* it = std::lower_bound(cmap, cmap + nOffCols, c);
                    if (it != cmap + nOffCols && *it == c) {
                        const HYPRE_Int oc = static_cast<HYPRE_Int>(it - cmap);
                        for (HYPRE_Int h = hypre_CSRMatrixI(Ao)[r]; h < hypre_CSRMatrixI(Ao)[r + 1]; ++h) {
                            if (hypre_CSRMatrixJ(Ao)[h] == oc) { slot = nnzD + h; break; }
                        }
                    }
                }
                if (slot < 0) {
                    return; // not found: keep using HYPRE_IJMatrixSetValues2
                }
                src.push_back(start + t);
                dst.push_back(slot);
            }
        }
        direct_src_ = std::move(src);
        direct_slots_ = std::move(dst);
        direct_nnz_diag_ = nnzD;
    }

    void copyValuesDirect_()
    {
        HYPRE_ParCSRMatrix parcsr_A;
        HYPRE_SAFE_CALL(HYPRE_IJMatrixGetObject(A_hypre_, reinterpret_cast<void**>(&parcsr_A)));
        auto* A = reinterpret_cast<hypre_ParCSRMatrix*>(parcsr_A);
        HYPRE_Real* d = hypre_CSRMatrixData(hypre_ParCSRMatrixDiag(A));
        HYPRE_Real* o = hypre_CSRMatrixData(hypre_ParCSRMatrixOffd(A));
        const auto* values = &(A_[0][0][0][0]);
        for (std::size_t t = 0; t < direct_slots_.size(); ++t) {
            const HYPRE_Int h = direct_slots_[t];
            (h < direct_nnz_diag_ ? d[h] : o[h - direct_nnz_diag_]) = values[direct_src_[t]];
        }
    }

    // Same l1 norms (option 4, no CF ordering) that BoomerAMGSetup builds for l1 smoothers.
    void refreshL1Norms_(HYPRE_Int level)
    {
        if (relax_type_ != 8 && relax_type_ != 13 && relax_type_ != 14) {
            return;
        }
        auto* amg = reinterpret_cast<hypre_ParAMGData*>(solver_);
        hypre_Vector** l1 = hypre_ParAMGDataL1Norms(amg);
        if (l1 == nullptr || l1[level] == nullptr) {
            return;
        }
        hypre_ParCSRMatrix* A = hypre_ParAMGDataAArray(amg)[level];
        HYPRE_Real* d = nullptr;
        HYPRE_SAFE_CALL(hypre_ParCSRComputeL1Norms(A, 4, nullptr, &d));
        hypre_SeqVectorDestroy(l1[level]);
        l1[level] = hypre_SeqVectorCreate(hypre_ParCSRMatrixNumRows(A));
        hypre_VectorData(l1[level]) = d;
        hypre_SeqVectorInitialize_v2(l1[level], hypre_ParCSRMatrixMemoryLocation(A));
    }

    void refreshFineL1Norms_(HYPRE_ParCSRMatrix /*parcsr_A*/)
    {
        refreshL1Norms_(0);
    }

    // Keep the coarsening and interpolation, rebuild the coarse operators from the current
    // matrix as BoomerAMGSetup does, then the coarsest elimination and the smoother norms.
    void resetupGalerkin_(HYPRE_ParCSRMatrix parcsr_A)
    {
        auto* amg = reinterpret_cast<hypre_ParAMGData*>(solver_);
        const HYPRE_Int L = hypre_ParAMGDataNumLevels(amg);
        hypre_ParCSRMatrix** A_array = hypre_ParAMGDataAArray(amg);
        hypre_ParCSRMatrix** P_array = hypre_ParAMGDataPArray(amg);
        A_array[0] = reinterpret_cast<hypre_ParCSRMatrix*>(parcsr_A);
        int num_procs = 1;
        MPI_Comm_size(hypre_ParCSRMatrixComm(A_array[0]), &num_procs);
        const HYPRE_Int keepT = hypre_ParAMGDataKeepTranspose(amg);
        if (numeric_rap_ && rap_maps_.empty()) {
            rap_maps_.resize(L - 1);
            for (HYPRE_Int l = 0; l < L - 1; ++l) {
                rap_maps_[l].build(A_array[l], P_array[l], A_array[l + 1]);
            }
        }
        for (HYPRE_Int l = 0; l < L - 1; ++l) {
            if (numeric_rap_ && rap_maps_[l].valid()) {
                rap_maps_[l].apply(A_array[l], P_array[l], A_array[l + 1]);
                continue;
            }
            hypre_ParCSRMatrix* AH = nullptr;
            if (hypre_ParAMGDataModularizedMatMat(amg)) {
                AH = hypre_ParCSRMatrixRAPKT(P_array[l], A_array[l], P_array[l], keepT, 1);
            } else {
                hypre_BoomerAMGBuildCoarseOperatorKT(P_array[l], A_array[l], P_array[l], keepT, &AH);
            }
            if (num_procs > 1 && hypre_ParCSRMatrixCommPkg(AH) == nullptr) {
                hypre_MatvecCommPkgCreate(AH);
            }
            hypre_ParCSRMatrixSetNumNonzeros(AH);
            hypre_ParCSRMatrixSetDNumNonzeros(AH);
            hypre_ParCSRMatrixDestroy(A_array[l + 1]);
            A_array[l + 1] = AH;
        }
        HYPRE_Int* grt = hypre_ParAMGDataGridRelaxType(amg);
        if (L > 1 && (grt[3] == 9 || grt[3] == 19 || grt[3] == 98 || grt[3] == 99
                      || grt[3] == 198 || grt[3] == 199)) {
            // free what hypre_GaussElimSetup allocated for the previous coarsest matrix
            const HYPRE_MemoryLocation gml = hypre_ParAMGDataGEMemoryLocation(amg);
            hypre_TFree(hypre_ParAMGDataAPiv(amg), gml);
            hypre_TFree(hypre_ParAMGDataAMat(amg), gml);
            hypre_TFree(hypre_ParAMGDataAWork(amg), gml);
            hypre_TFree(hypre_ParAMGDataBVec(amg), gml);
            hypre_TFree(hypre_ParAMGDataUVec(amg), gml);
            hypre_TFree(hypre_ParAMGDataCommInfo(amg), HYPRE_MEMORY_HOST);
            if (hypre_ParAMGDataNewComm(amg) != hypre_MPI_COMM_NULL) {
                MPI_Comm nc = hypre_ParAMGDataNewComm(amg);
                hypre_MPI_Comm_free(&nc);
                hypre_ParAMGDataNewComm(amg) = hypre_MPI_COMM_NULL;
            }
            hypre_GaussElimSetup(amg, L - 1, grt[3]);
        }
        for (HYPRE_Int l = 0; l < L - 1; ++l) {
            refreshL1Norms_(l);
        }
    }

    // Reference to the input matrix
    const M& A_;
    const Comm& comm_; //!< The communication object for parallel operations.

    // Parallel information and sparsity pattern from HypreInterface
    HypreInterface::ParallelInfo par_info_;
    HypreInterface::SparsityPattern sparsity_pattern_;
    HypreInterface::HostDataArrays host_arrays_;
    HypreInterface::DeviceDataArrays device_arrays_;

    // Hypre handles
    HYPRE_Solver solver_ = nullptr;
    HYPRE_IJMatrix A_hypre_ = nullptr;
    HYPRE_IJVector x_hypre_ = nullptr;
    HYPRE_IJVector b_hypre_ = nullptr;

    // Backend configuration
    bool use_gpu_backend_ = false;

    int reuse_interval_ = 1;
    bool external_reuse_ = false;
    bool reuse_galerkin_ = false;
    bool numeric_rap_ = true;
    std::vector<HYPRE_Int> direct_src_, direct_slots_;
    HYPRE_Int direct_nnz_diag_ = 0;
    std::vector<HypreInterface::NumericRap> rap_maps_;
    int relax_type_ = 13;
    int updates_since_setup_ = 0;
    bool is_setup_ = false;
};

} // namespace Hypre

#endif
