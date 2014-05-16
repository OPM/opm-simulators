/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.

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

#include <config.h>

#include <opm/autodiff/NewtonIterationBlackoilCPR.hpp>
#include <opm/autodiff/CPRPreconditioner.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/linalg/LinearSolverFactory.hpp>

#include "disable_warning_pragmas.h"

#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/io.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/schwarz.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/kamg.hh>
#include <dune/istl/paamg/pinfo.hh>

#include "reenable_warning_pragmas.h"


namespace Opm
{


    typedef AutoDiffBlock<double> ADB;
    typedef ADB::V V;
    typedef ADB::M M;

    typedef Dune::FieldVector<double, 1   > VectorBlockType;
    typedef Dune::FieldMatrix<double, 1, 1> MatrixBlockType;
    typedef Dune::BCRSMatrix <MatrixBlockType>        Mat;
    typedef Dune::BlockVector<VectorBlockType>        Vector;


    namespace {
        /// Eliminate a variable via Schur complement.
        /// \param[in]  eqs  set of equations with Jacobians
        /// \param[in]  n    index of equation/variable to eliminate.
        /// \return          new set of equations, one smaller than eqs.
        /// Note: this method requires the eliminated variable to have the same size
        /// as the equation in the corresponding position (that also will be eliminated).
        /// It also required the jacobian block n of equation n to be diagonal.
        std::vector<ADB> eliminateVariable(const std::vector<ADB>& eqs, const int n);

        /// Recover that value of a variable previously eliminated.
        /// \param[in]  equation          previously eliminated equation.
        /// \param[in]  partial_solution  solution to the remainder system after elimination.
        /// \param[in]  n                 index of equation/variable that was eliminated.
        /// \return                       solution to complete system.
        V recoverVariable(const ADB& equation, const V& partial_solution, const int n);

        /// Determine diagonality of a sparse matrix.
        /// If there are off-diagonal elements in the sparse
        /// structure, this function returns true if they are all
        /// equal to zero.
        /// \param[in]  matrix  the matrix under consideration
        /// \return             true if matrix is diagonal
        bool isDiagonal(const M& matrix);

        /// Create a dune-istl matrix from an Eigen matrix.
        /// \param[in]  matrix       input Eigen::SparseMatrix
        /// \return                  output Dune::BCRSMatrix
        Mat makeIstlMatrix(const Eigen::SparseMatrix<double, Eigen::RowMajor>& matrix);

    } // anonymous namespace





    /// Construct a system solver.
    /// \param[in] linsolver   linear solver to use
    NewtonIterationBlackoilCPR::NewtonIterationBlackoilCPR(const parameter::ParameterGroup& param)
    {
        parameter::ParameterGroup cpr_elliptic = param.getDefault("cpr_elliptic", param);
        linsolver_elliptic_.reset(new LinearSolverFactory(cpr_elliptic));
        parameter::ParameterGroup cpr_full = param.getDefault("cpr_full", param);
        linsolver_full_.reset(new LinearSolverFactory(cpr_full));
    }





    /// Solve the linear system Ax = b, with A being the
    /// combined derivative matrix of the residual and b
    /// being the residual itself.
    /// \param[in] residual   residual object containing A and b.
    /// \return               the solution x
    NewtonIterationBlackoilCPR::SolutionVector
    NewtonIterationBlackoilCPR::computeNewtonIncrement(const LinearisedBlackoilResidual& residual) const
    {
        // Build the vector of equations.
        const int np = residual.material_balance_eq.size();
        std::vector<ADB> eqs;
        eqs.reserve(np + 2);
        for (int phase = 0; phase < np; ++phase) {
            eqs.push_back(residual.material_balance_eq[phase]);
        }
        eqs.push_back(residual.well_flux_eq);
        eqs.push_back(residual.well_eq);

        // Eliminate the well-related unknowns, and corresponding equations.
        std::vector<ADB> elim_eqs;
        elim_eqs.reserve(2);
        elim_eqs.push_back(eqs[np]);
        eqs = eliminateVariable(eqs, np); // Eliminate well flux unknowns.
        elim_eqs.push_back(eqs[np]);
        eqs = eliminateVariable(eqs, np); // Eliminate well bhp unknowns.
        assert(int(eqs.size()) == np);

        // Scale material balance equations.
        const double matbalscale[3] = { 1.1169, 1.0031, 0.0031 }; // HACK hardcoded instead of computed.
        for (int phase = 0; phase < np; ++phase) {
            eqs[phase] = eqs[phase] * matbalscale[phase];
        }

        // Add material balance equations to form pressure equation
        // in place of first balance equation.
        for (int phase = 1; phase < np; ++phase) {
            eqs[0] += eqs[phase];
        }

        // Scale pressure equation.
        const double pscale = 200*unit::barsa;
        eqs[0] = eqs[0] * pscale;

        // Combine in single block.
        ADB total_residual = eqs[0];
        for (int phase = 1; phase < np; ++phase) {
            total_residual = vertcat(total_residual, eqs[phase]);
        }
        total_residual = collapseJacs(total_residual);

        // Solve reduced system.
        const Eigen::SparseMatrix<double, Eigen::RowMajor> matr = total_residual.derivative()[0];
        SolutionVector dx(SolutionVector::Zero(total_residual.size()));
        /*
        Opm::LinearSolverInterface::LinearSolverReport rep
            = linsolver_full_->solve(matr.rows(), matr.nonZeros(),
                                     matr.outerIndexPtr(), matr.innerIndexPtr(), matr.valuePtr(),
                                     total_residual.value().data(), dx.data());
        if (!rep.converged) {
            OPM_THROW(std::runtime_error,
                      "FullyImplicitBlackoilSolver::solveJacobianSystem(): "
                      "Linear solver convergence failure.");
        }
        */


        // Create ISTL matrix.
        Mat A = makeIstlMatrix(matr);

        // Create ISTL matrix for elliptic part.
        const int nc = residual.material_balance_eq[0].size();
        Mat Ae = makeIstlMatrix(matr.topLeftCorner(nc, nc));

        // Construct operator, scalar product and vectors needed.
        typedef Dune::MatrixAdapter<Mat,Vector,Vector> Operator;
        Operator opA(A);
        Dune::SeqScalarProduct<Vector> sp;
        // Right hand side.
        Vector b(opA.getmat().N());
        std::copy_n(total_residual.value().data(), b.size(), b.begin());
        // System solution
        Vector x(opA.getmat().M());
        x = 0.0;

        // Construct preconditioner.
        // typedef Dune::SeqILU0<Mat,Vector,Vector> Preconditioner;
        typedef Opm::CPRPreconditioner<Mat,Vector,Vector> Preconditioner;
        const double relax = 1.0;
        Preconditioner precond(A, Ae, relax);

        // Construct linear solver.
        const double tolerance = 1e-3;
        const int maxit = 5000;
        const int verbosity = 1;
        // Dune::BiCGSTABSolver<Vector> linsolve(opA, sp, precond, tolerance, maxit, verbosity);
        const int restart = 40;
        Dune::RestartedGMResSolver<Vector> linsolve(opA, sp, precond, tolerance, restart, maxit, verbosity);

        // Solve system.
        Dune::InverseOperatorResult result;
        linsolve.apply(x, b, result);

        // Output results.
        LinearSolverInterface::LinearSolverReport res;
        res.converged = result.converged;
        res.iterations = result.iterations;
        res.residual_reduction = result.reduction;

        std::copy(x.begin(), x.end(), dx.data());






        // Compute full solution using the eliminated equations.
        // Recovery in inverse order of elimination.
        dx = recoverVariable(elim_eqs[1], dx, np);
        dx = recoverVariable(elim_eqs[0], dx, np);
        return dx;
    }




    namespace
    {


        std::vector<ADB> eliminateVariable(const std::vector<ADB>& eqs, const int n)
        {
            // Check that the variable index to eliminate is within bounds.
            const int num_eq = eqs.size();
            const int num_vars = eqs[0].derivative().size();
            if (num_eq != num_vars) {
                OPM_THROW(std::logic_error, "eliminateVariable() requires the same number of variables and equations.");
            }
            if (n >= num_eq) {
                OPM_THROW(std::logic_error, "Trying to eliminate variable from too small set of equations.");
            }

            // Schur complement of (A B ; C D) wrt. D is A - B*inv(D)*C.
            // This is applied to all 2x2 block submatrices.
            // We require that D is diagonal.
            const M& D = eqs[n].derivative()[n];
            if (!isDiagonal(D)) {
                // std::cout << "++++++++++++++++++++++++++++++++++++++++++++\n"
                //           << D
                //           << "++++++++++++++++++++++++++++++++++++++++++++\n" << std::endl;
                OPM_THROW(std::logic_error, "Cannot do Schur complement with respect to non-diagonal block.");
            }
            V diag = D.diagonal();
            Eigen::DiagonalMatrix<double, Eigen::Dynamic> invD = (1.0 / diag).matrix().asDiagonal();
            std::vector<V> vals(num_eq);              // Number n will remain empty.
            std::vector<std::vector<M>> jacs(num_eq); // Number n will remain empty.
            for (int eq = 0; eq < num_eq; ++eq) {
                if (eq == n) {
                    continue;
                }
                jacs[eq].reserve(num_eq - 1);
                const M& B = eqs[eq].derivative()[n];
                for (int var = 0; var < num_eq; ++var) {
                    if (var == n) {
                        continue;
                    }
                    // Create new jacobians.
                    M schur_jac = eqs[eq].derivative()[var] - B * (invD * eqs[n].derivative()[var]);
                    jacs[eq].push_back(schur_jac);
                }
                // Update right hand side.
                vals[eq] = eqs[eq].value().matrix() - B * (invD * eqs[n].value().matrix());
            }

            // Create return value.
            std::vector<ADB> retval;
            retval.reserve(num_eq - 1);
            for (int eq = 0; eq < num_eq; ++eq) {
                if (eq == n) {
                    continue;
                }
                retval.push_back(ADB::function(vals[eq], jacs[eq]));
            }
            return retval;
        }





        V recoverVariable(const ADB& equation, const V& partial_solution, const int n)
        {
            // The equation to solve for the unknown y (to be recovered) is
            //    Cx + Dy = b
            //    y = inv(D) (b - Cx)
            // where D is the eliminated block, C is the jacobian of
            // the eliminated equation with respect to the
            // non-eliminated unknowms, b is the right-hand side of
            // the eliminated equation, and x is the partial solution
            // of the non-eliminated unknowns.
            // We require that D is diagonal.

            // Find inv(D).
            const M& D = equation.derivative()[n];
            if (!isDiagonal(D)) {
                OPM_THROW(std::logic_error, "Cannot do Schur complement with respect to non-diagonal block.");
            }
            V diag = D.diagonal();
            Eigen::DiagonalMatrix<double, Eigen::Dynamic> invD = (1.0 / diag).matrix().asDiagonal();

            // Build C.
            std::vector<M> C_jacs = equation.derivative();
            C_jacs.erase(C_jacs.begin() + n);
            ADB eq_coll = collapseJacs(ADB::function(equation.value(), C_jacs));
            const M& C = eq_coll.derivative()[0];

            // Compute value of eliminated variable.
            V elim_var = invD * (equation.value().matrix() - C * partial_solution.matrix());

            // Find the relevant sizes to use when reconstructing the full solution.
            const int nelim = equation.size();
            const int npart = partial_solution.size();
            assert(C.cols() == npart);
            const int full_size = nelim + npart;
            int start = 0;
            for (int i = 0; i < n; ++i) {
                start += equation.derivative()[i].cols();
            }
            assert(start < full_size);

            // Reconstruct complete solution vector.
            V sol(full_size);
            std::copy_n(partial_solution.data(), start, sol.data());
            std::copy_n(elim_var.data(), nelim, sol.data() + start);
            std::copy_n(partial_solution.data() + start, npart - start, sol.data() + start + nelim);
            return sol;
        }





        bool isDiagonal(const M& matr)
        {
            M matrix = matr;
            matrix.makeCompressed();
            for (int k = 0; k < matrix.outerSize(); ++k) {
                for (M::InnerIterator it(matrix, k); it; ++it) {
                    if (it.col() != it.row()) {
                        // Off-diagonal element.
                        if (it.value() != 0.0) {
                            // Nonzero off-diagonal element.
                            // std::cout << "off-diag: " << it.row() << ' ' << it.col() << std::endl;
                            return false;
                        }
                    }
                }
            }
            return true;
        }




        Mat makeIstlMatrix(const Eigen::SparseMatrix<double, Eigen::RowMajor>& matrix)
        {
            // Create ISTL matrix.
            const int size = matrix.rows();
            const int nonzeros = matrix.nonZeros();
            const int* ia = matrix.outerIndexPtr();
            const int* ja = matrix.innerIndexPtr();
            const double* sa = matrix.valuePtr();
            Mat A(size, size, nonzeros, Mat::row_wise);
            for (Mat::CreateIterator row = A.createbegin(); row != A.createend(); ++row) {
                const int ri = row.index();
                for (int i = ia[ri]; i < ia[ri + 1]; ++i) {
                    row.insert(ja[i]);
                }
            }
            for (int ri = 0; ri < size; ++ri) {
                for (int i = ia[ri]; i < ia[ri + 1]; ++i) {
                    A[ri][ja[i]] = sa[i];
                }
            }
            return A;
        }



    } // anonymous namespace


} // namespace Opm

