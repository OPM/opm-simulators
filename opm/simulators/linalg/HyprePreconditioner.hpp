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

        // Initialize Hypre
        HYPRE_Init();

        // Create the solver (BoomerAMG)
        HYPRE_BoomerAMGCreate(&solver_);

        // Set some default parameters
        HYPRE_BoomerAMGSetPrintLevel(solver_, 0);         // Reduce output
        HYPRE_BoomerAMGSetMaxIter(solver_, 1);            // Only one V-cycle, as used as preconditioner
        HYPRE_BoomerAMGSetCoarsenType(solver_, 10);       // HMIS coarsening
        HYPRE_BoomerAMGSetStrongThreshold(solver_, 0.5);  // Strength threshold for 3D
        HYPRE_BoomerAMGSetAggNumLevels(solver_, 1);       // Aggressive coarsening on first level
        HYPRE_BoomerAMGSetAggTruncFactor(solver_, 0.3);   // Remove weak connections
        HYPRE_BoomerAMGSetInterpType(solver_, 6);         // ext+i interpolation
        HYPRE_BoomerAMGSetMaxLevels(solver_, 15);         // Maximum number of levels
        HYPRE_BoomerAMGSetTol(solver_, 0);              // Convergence tolerance, 0 as used as preconditioner with one V-cycle


        // Create Hypre vectors
        const int N = A_.N();
        HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, N-1, &x_hypre_);
        HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, N-1, &b_hypre_);
        HYPRE_IJVectorSetObjectType(x_hypre_, HYPRE_PARCSR);
        HYPRE_IJVectorSetObjectType(b_hypre_, HYPRE_PARCSR);
        HYPRE_IJVectorInitialize(x_hypre_);
        HYPRE_IJVectorInitialize(b_hypre_);
        // Create indices vector
        indices_.resize(A_.N());
        std::iota(indices_.begin(), indices_.end(), 0);

        // Create Hypre matrix
        HYPRE_IJMatrixCreate(MPI_COMM_WORLD, 0, N-1, 0, N-1, &A_hypre_);
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
        if (parcsr_A_) {
            HYPRE_IJMatrixDestroy(A_hypre_);
        }
        HYPRE_Finalize();
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
        const int N = A_.N();
        const int nnz = A_.nonzeroes();

        // Allocate arrays required by Hypre
        ncols_.resize(N);
        rows_.resize(N);
        cols_.resize(nnz);

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
    }

    void copyMatrixToHypre() {
        // Get pointer to matrix values array
        const double* vals = &(A_[0][0][0][0]);  // Indexing explanation:
                           // A_[row]              - First row of the matrix
                           //     [0]              - First block in that row
                           //        [0]           - First row within the 1x1 block
                           //           [0]        - First column within the 1x1 block

        // Set all values at once using stored sparsity pattern
        HYPRE_IJMatrixSetValues(A_hypre_, A_.N(), ncols_.data(), rows_.data(), cols_.data(), vals);

        HYPRE_IJMatrixAssemble(A_hypre_);
        HYPRE_IJMatrixGetObject(A_hypre_, (void**)&parcsr_A_);
    }

    void copyVectorsToHypre(const X& v, const Y& d) {
        const double* x_vals = &(v[0][0]);
        const double* b_vals = &(d[0][0]);

        HYPRE_IJVectorSetValues(x_hypre_, A_.N(), indices_.data(), x_vals);
        HYPRE_IJVectorSetValues(b_hypre_, A_.N(), indices_.data(), b_vals);

        HYPRE_IJVectorGetObject(x_hypre_, (void**)&par_x_);
        HYPRE_IJVectorGetObject(b_hypre_, (void**)&par_b_);
    }

    void copyVectorFromHypre(X& v) {
        double* vals = &(v[0][0]);
        HYPRE_IJVectorGetValues(x_hypre_, A_.N(), indices_.data(), vals);
    }

    const M& A_;
    const Opm::PropertyTree& prm_;

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

    // Store indices vector
    std::vector<int> indices_;
};

} // namespace Dune

#endif
