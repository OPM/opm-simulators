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

namespace Dune {

/// Wrapper for Hypre's BoomerAMG preconditioner
template<class M, class X, class Y>
class HyprePreconditioner : public PreconditionerWithUpdate<X,Y> {
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
    HyprePreconditioner (const M& A)
        : A_(A)
    {
        OPM_TIMEBLOCK(prec_construct);

        // Initialize Hypre
        HYPRE_Init();

        // Create the solver (BoomerAMG)
        HYPRE_BoomerAMGCreate(&solver_);

        // Set some default parameters
        HYPRE_BoomerAMGSetPrintLevel(solver_, 0);  // Reduce output
        HYPRE_BoomerAMGSetOldDefault(solver_);     // Falgout coarsening with modified classical interpolation
        HYPRE_BoomerAMGSetRelaxType(solver_, 3);   // G-S/Jacobi hybrid relaxation
        HYPRE_BoomerAMGSetRelaxOrder(solver_, 1);  // Uses C/F relaxation
        HYPRE_BoomerAMGSetNumSweeps(solver_, 1);   // Sweeps on each level
        HYPRE_BoomerAMGSetMaxLevels(solver_, 20);  // Maximum number of levels
        HYPRE_BoomerAMGSetTol(solver_, 1e-7);      // Convergence tolerance


        // Create the vectors
        const int N = A_.N();
        HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, N-1, &x_hypre_);
        HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, N-1, &b_hypre_);
        HYPRE_IJVectorSetObjectType(x_hypre_, HYPRE_PARCSR);
        HYPRE_IJVectorSetObjectType(b_hypre_, HYPRE_PARCSR);
        HYPRE_IJVectorInitialize(x_hypre_);
        HYPRE_IJVectorInitialize(b_hypre_);

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

    SolverCategory::Category category() const override {
        return SolverCategory::sequential;
    }

    bool hasPerfectUpdate() const override
    {
        // The Hypre preconditioner can depend on the values of the matrix, so it must be recreated
        return false;
    }

private:
    void copyMatrixToHypre() {
        const int N = A_.N();
        HYPRE_IJMatrixCreate(MPI_COMM_WORLD, 0, N-1, 0, N-1, &A_hypre_);
        HYPRE_IJMatrixSetObjectType(A_hypre_, HYPRE_PARCSR);
        HYPRE_IJMatrixInitialize(A_hypre_);

        // Copy values from Dune matrix to Hypre matrix
        for (auto row = A_.begin(); row != A_.end(); ++row) {
            int i = row.index();
            std::vector<int> cols;
            std::vector<double> vals;

            for (auto col = row->begin(); col != row->end(); ++col) {
                cols.push_back(col.index());
                vals.push_back(*col);
            }

            int ncols = cols.size();
            HYPRE_IJMatrixSetValues(A_hypre_, 1, &ncols, &i, cols.data(), vals.data());
        }

        HYPRE_IJMatrixAssemble(A_hypre_);
        HYPRE_IJMatrixGetObject(A_hypre_, (void**)&parcsr_A_);

        // Setup the solver with the new matrix
        HYPRE_BoomerAMGSetup(solver_, parcsr_A_, par_b_, par_x_);
    }

    void copyVectorsToHypre(const X& v, const Y& d) {
        const int N = A_.N();

        // Copy values
        std::vector<int> indices(N);
        std::vector<double> x_vals(N);
        std::vector<double> b_vals(N);

        for (int i = 0; i < N; ++i) {
            indices[i] = i;
            x_vals[i] = v[i];
            b_vals[i] = d[i];
        }

        HYPRE_IJVectorSetValues(x_hypre_, N, indices.data(), x_vals.data());
        HYPRE_IJVectorSetValues(b_hypre_, N, indices.data(), b_vals.data());

        HYPRE_IJVectorGetObject(x_hypre_, (void**)&par_x_);
        HYPRE_IJVectorGetObject(b_hypre_, (void**)&par_b_);
    }

    void copyVectorFromHypre(X& v) {
        const int N = A_.N();
        std::vector<int> indices(N);
        std::vector<double> vals(N);

        for (int i = 0; i < N; ++i) {
            indices[i] = i;
        }

        HYPRE_IJVectorGetValues(x_hypre_, N, indices.data(), vals.data());

        for (int i = 0; i < N; ++i) {
            v[i] = vals[i];
        }
    }

    const M& A_;
    HYPRE_Solver solver_ = nullptr;
    HYPRE_IJMatrix A_hypre_ = nullptr;
    HYPRE_ParCSRMatrix parcsr_A_ = nullptr;
    HYPRE_IJVector x_hypre_ = nullptr;
    HYPRE_IJVector b_hypre_ = nullptr;
    HYPRE_ParVector par_x_ = nullptr;
    HYPRE_ParVector par_b_ = nullptr;
};

} // namespace Dune

#endif
