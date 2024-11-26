/*
  Copyright 2024 SINTEF AS
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

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>
#include <HYPRE_krylov.h>
#include <_hypre_utilities.h>

#include <memory>
#include <vector>
#include <numeric>

namespace Hypre {

/// Wrapper for Hypre's BoomerAMG preconditioner
template<class M, class X, class Y>
class HyprePreconditioner : public Dune::PreconditionerWithUpdate<X,Y> {
public:
    //! \brief The matrix type the preconditioner is for
    using matrix_type = M;
    //! \brief The domain type of the preconditioner
    using domain_type = X;
    //! \brief The range type of the preconditioner
    using range_type = Y;
    //! \brief The field type of the preconditioner
    using field_type = typename X::field_type;

    // Constructor
    HyprePreconditioner (const M& A, const Opm::PropertyTree& prm)
        : A_(A), prm_(prm)
    {
        OPM_TIMEBLOCK(prec_construct);

        use_gpu_ = prm_.get<bool>("use_gpu", false);

        // Set memory location and execution policy
        if (use_gpu_) {
            HYPRE_SetMemoryLocation(HYPRE_MEMORY_DEVICE);
            HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE);
            // use hypre's SpGEMM instead of vendor implementation
            HYPRE_SetSpGemmUseVendor(false);
            // use cuRand for PMIS
            HYPRE_SetUseGpuRand(1);
            HYPRE_DeviceInitialize();
            HYPRE_PrintDeviceInfo();
        }
        else {
            HYPRE_SetMemoryLocation(HYPRE_MEMORY_HOST);
            HYPRE_SetExecutionPolicy(HYPRE_EXEC_HOST);
        }

        // Create the solver (BoomerAMG)
        HYPRE_BoomerAMGCreate(&solver_);

        // Set parameters from property tree with defaults
        HYPRE_BoomerAMGSetPrintLevel(solver_, prm_.get<int>("print_level", 0));
        HYPRE_BoomerAMGSetMaxIter(solver_, prm_.get<int>("max_iter", 1));
        HYPRE_BoomerAMGSetStrongThreshold(solver_, prm_.get<double>("strong_threshold", 0.5));
        HYPRE_BoomerAMGSetAggTruncFactor(solver_, prm_.get<double>("agg_trunc_factor", 0.3));
        HYPRE_BoomerAMGSetInterpType(solver_, prm_.get<int>("interp_type", 6));
        HYPRE_BoomerAMGSetMaxLevels(solver_, prm_.get<int>("max_levels", 15));
        HYPRE_BoomerAMGSetTol(solver_, prm_.get<double>("tolerance", 0.0));

        if (use_gpu_) {
            HYPRE_BoomerAMGSetRelaxType(solver_, 16);
            HYPRE_BoomerAMGSetCoarsenType(solver_, 8);
            HYPRE_BoomerAMGSetAggNumLevels(solver_, 0);
            HYPRE_BoomerAMGSetAggInterpType(solver_, 6);
            // Keep transpose to avoid SpMTV
            HYPRE_BoomerAMGSetKeepTranspose(solver_, true);
        }
        else {
            HYPRE_BoomerAMGSetRelaxType(solver_, prm_.get<int>("relax_type", 13));
            HYPRE_BoomerAMGSetCoarsenType(solver_, prm_.get<int>("coarsen_type", 10));
            HYPRE_BoomerAMGSetAggNumLevels(solver_, prm_.get<int>("agg_num_levels", 1));
            HYPRE_BoomerAMGSetAggInterpType(solver_, prm_.get<int>("agg_interp_type", 4));
        }

        // Create Hypre vectors
        N_ = A_.N();
        nnz_ = A_.nonzeroes();
        HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, N_-1, &x_hypre_);
        HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, N_-1, &b_hypre_);
        HYPRE_IJVectorSetObjectType(x_hypre_, HYPRE_PARCSR);
        HYPRE_IJVectorSetObjectType(b_hypre_, HYPRE_PARCSR);
        HYPRE_IJVectorInitialize(x_hypre_);
        HYPRE_IJVectorInitialize(b_hypre_);
        // Create indices vector
        indices_.resize(N_);
        std::iota(indices_.begin(), indices_.end(), 0);
        if (use_gpu_) {
            indices_device_ = hypre_CTAlloc(HYPRE_BigInt, N_, HYPRE_MEMORY_DEVICE);
            hypre_TMemcpy(indices_device_, indices_.data(), HYPRE_BigInt, N_, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
            // Allocate device vectors
            x_values_device_ = hypre_CTAlloc(HYPRE_Real, N_, HYPRE_MEMORY_DEVICE);
            b_values_device_ = hypre_CTAlloc(HYPRE_Real, N_, HYPRE_MEMORY_DEVICE);
        }

        // Create Hypre matrix
        HYPRE_IJMatrixCreate(MPI_COMM_WORLD, 0, N_-1, 0, N_-1, &A_hypre_);
        HYPRE_IJMatrixSetObjectType(A_hypre_, HYPRE_PARCSR);
        HYPRE_IJMatrixInitialize(A_hypre_);

        setupSparsityPattern();
        update();
    }

    // Destructor
    ~HyprePreconditioner() {
        if (solver_) {
            HYPRE_BoomerAMGDestroy(solver_);
        }
        if (A_hypre_) {
            HYPRE_IJMatrixDestroy(A_hypre_);
        }
        if (x_hypre_) {
            HYPRE_IJVectorDestroy(x_hypre_);
        }
        if (b_hypre_) {
            HYPRE_IJVectorDestroy(b_hypre_);
        }
        if (values_device_) {
            hypre_TFree(values_device_, HYPRE_MEMORY_DEVICE);
        }
        if (x_values_device_) {
            hypre_TFree(x_values_device_, HYPRE_MEMORY_DEVICE);
        }
        if (b_values_device_) {
            hypre_TFree(b_values_device_, HYPRE_MEMORY_DEVICE);
        }
        if (indices_device_) {
            hypre_TFree(indices_device_, HYPRE_MEMORY_DEVICE);
        }
    }

    void update() override {
        OPM_TIMEBLOCK(prec_update);
        copyMatrixToHypre();
        HYPRE_BoomerAMGSetup(solver_, parcsr_A_, par_b_, par_x_);
    }

    void pre(X& x, Y& b) override {
        DUNE_UNUSED_PARAMETER(x);
        DUNE_UNUSED_PARAMETER(b);
    }

    void apply(X& v, const Y& d) override {
        OPM_TIMEBLOCK(prec_apply);

        // Copy vectors to Hypre format
        copyVectorsToHypre(v, d);

        // Apply the preconditioner (one AMG V-cycle)
        HYPRE_BoomerAMGSolve(solver_, parcsr_A_, par_b_, par_x_);

        // Copy result back
        copyVectorFromHypre(v);
    }

    void post(X& x) override {
        DUNE_UNUSED_PARAMETER(x);
    }

    Dune::SolverCategory::Category category() const override {
        return Dune::SolverCategory::sequential;
    }

    bool hasPerfectUpdate() const override
    {
        // The Hypre preconditioner can depend on the values of the matrix, so it must be recreated
        return false;
    }

private:
    void setupSparsityPattern() {
        // Allocate arrays required by Hypre
        ncols_.resize(N_);
        rows_.resize(N_);
        cols_.resize(nnz_);

        // Setup arrays and fill column indices
        int pos = 0;
        for (auto row = A_.begin(); row != A_.end(); ++row) {
            const int rowIdx = row.index();
            rows_[rowIdx] = rowIdx;
            ncols_[rowIdx] = row->size();

            for (auto col = row->begin(); col != row->end(); ++col) {
                cols_[pos++] = col.index();
            }
        }
        if (use_gpu_) {
            // Allocate device arrays
            ncols_device_ = hypre_CTAlloc(HYPRE_Int, N_, HYPRE_MEMORY_DEVICE);
            rows_device_ = hypre_CTAlloc(HYPRE_BigInt, N_, HYPRE_MEMORY_DEVICE);
            cols_device_ = hypre_CTAlloc(HYPRE_BigInt, nnz_, HYPRE_MEMORY_DEVICE);
            values_device_ = hypre_CTAlloc(HYPRE_Real, nnz_, HYPRE_MEMORY_DEVICE);

            // Copy to device
            hypre_TMemcpy(ncols_device_, ncols_.data(), HYPRE_Int, N_, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
            hypre_TMemcpy(rows_device_, rows_.data(), HYPRE_BigInt, N_, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
            hypre_TMemcpy(cols_device_, cols_.data(), HYPRE_BigInt, nnz_, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
        }
    }

    void copyMatrixToHypre() {
        OPM_TIMEBLOCK(prec_copy_matrix);
        // Get pointer to matrix values array
        const HYPRE_Real* values = &(A_[0][0][0][0]);
        // Indexing explanation:
        // A_[row]           - First row of the matrix
        //     [0]           - First block in that row
        //        [0]        - First row within the 1x1 block
        //           [0]     - First column within the 1x1 block

        if (use_gpu_) {
            hypre_TMemcpy(values_device_, values, HYPRE_Real, nnz_, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
            HYPRE_IJMatrixSetValues(A_hypre_, N_, ncols_device_, rows_device_, cols_device_, values_device_);
        }
        else {
            HYPRE_IJMatrixSetValues(A_hypre_, N_, ncols_.data(), rows_.data(), cols_.data(), values);
        }

        HYPRE_IJMatrixAssemble(A_hypre_);
        HYPRE_IJMatrixGetObject(A_hypre_, (void**)&parcsr_A_);
    }

    void copyVectorsToHypre(const X& v, const Y& d) {
        OPM_TIMEBLOCK(prec_copy_vectors_to_hypre);
        const HYPRE_Real* x_vals = &(v[0][0]);
        const HYPRE_Real* b_vals = &(d[0][0]);

        if (use_gpu_) {
            hypre_TMemcpy(x_values_device_, x_vals, HYPRE_Real, N_, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);
            hypre_TMemcpy(b_values_device_, b_vals, HYPRE_Real, N_, HYPRE_MEMORY_DEVICE, HYPRE_MEMORY_HOST);

            HYPRE_IJVectorSetValues(x_hypre_, N_, indices_device_, x_values_device_);
            HYPRE_IJVectorSetValues(b_hypre_, N_, indices_device_, b_values_device_);
        }
        else {
            HYPRE_IJVectorSetValues(x_hypre_, N_, indices_.data(), x_vals);
            HYPRE_IJVectorSetValues(b_hypre_, N_, indices_.data(), b_vals);
        }

        HYPRE_IJVectorAssemble(x_hypre_);
        HYPRE_IJVectorAssemble(b_hypre_);
        HYPRE_IJVectorGetObject(x_hypre_, (void**)&par_x_);
        HYPRE_IJVectorGetObject(b_hypre_, (void**)&par_b_);
    }

    void copyVectorFromHypre(X& v) {
        OPM_TIMEBLOCK(prec_copy_vector_from_hypre);
        HYPRE_Real* values = &(v[0][0]);
        if (use_gpu_) {
            HYPRE_IJVectorGetValues(x_hypre_, N_, indices_device_, x_values_device_);
            hypre_TMemcpy(values, x_values_device_, HYPRE_Real, N_, HYPRE_MEMORY_HOST, HYPRE_MEMORY_DEVICE);
        }
        else {
            HYPRE_IJVectorGetValues(x_hypre_, N_, indices_.data(), values);
        }
    }

    const M& A_;
    const Opm::PropertyTree& prm_;
    bool use_gpu_ = false;

    HYPRE_Solver solver_ = nullptr;
    HYPRE_IJMatrix A_hypre_ = nullptr;
    HYPRE_ParCSRMatrix parcsr_A_ = nullptr;
    HYPRE_IJVector x_hypre_ = nullptr;
    HYPRE_IJVector b_hypre_ = nullptr;
    HYPRE_ParVector par_x_ = nullptr;
    HYPRE_ParVector par_b_ = nullptr;

    // Store sparsity pattern
    std::vector<HYPRE_Int> ncols_;
    std::vector<HYPRE_BigInt> rows_;
    std::vector<HYPRE_BigInt> cols_;
    HYPRE_Int* ncols_device_ = nullptr;
    HYPRE_BigInt* rows_device_ = nullptr;
    HYPRE_BigInt* cols_device_ = nullptr;
    HYPRE_Real* values_device_ = nullptr;
    // Store indices vector
    std::vector<int> indices_;
    HYPRE_BigInt* indices_device_ = nullptr;
    HYPRE_Int N_ = -1;
    HYPRE_Int nnz_ = -1;

    HYPRE_Real* x_values_device_ = nullptr;
    HYPRE_Real* b_values_device_ = nullptr;
};

} // namespace Hypre

#endif
