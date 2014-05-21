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

        /// Form an elliptic system of equations.
        /// \param[in]       num_phases  the number of fluid phases
        /// \param[in]       eqs         the equations
        /// \param[out]      A           the resulting full system matrix
        /// \param[out]      b           the right hand side
        /// This function will deal with the first num_phases
        /// equations in eqs, and return a matrix A for the full
        /// system that has a elliptic upper left corner, if possible.
        void formEllipticSystem(const int num_phases,
                                const std::vector<ADB>& eqs,
                                Eigen::SparseMatrix<double, Eigen::RowMajor>& A,
                                V& b);

        /// Create a dune-istl matrix from an Eigen matrix.
        /// \param[in]  matrix       input Eigen::SparseMatrix
        /// \return                  output Dune::BCRSMatrix
        Mat makeIstlMatrix(const Eigen::SparseMatrix<double, Eigen::RowMajor>& matrix);

    } // anonymous namespace





    /// Construct a system solver.
    /// \param[in] linsolver   linear solver to use
    NewtonIterationBlackoilCPR::NewtonIterationBlackoilCPR(const parameter::ParameterGroup& /*param*/)
    {
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

        // Add material balance equations (or other manipulations) to
        // form pressure equation in top left of full system.
        Eigen::SparseMatrix<double, Eigen::RowMajor> A;
        V b;
        formEllipticSystem(np, eqs, A, b);

        // Scale pressure equation.
        const double pscale = 200*unit::barsa;
        const int nc = residual.material_balance_eq[0].size();
        A.topRows(nc) *= pscale;
        b.topRows(nc) *= pscale;

        // Solve reduced system.
        SolutionVector dx(SolutionVector::Zero(b.size()));

        // Create ISTL matrix.
        Mat istlA = makeIstlMatrix(A);

        // Create ISTL matrix for elliptic part.
        Mat istlAe = makeIstlMatrix(A.topLeftCorner(nc, nc));

        // Construct operator, scalar product and vectors needed.
        typedef Dune::MatrixAdapter<Mat,Vector,Vector> Operator;
        Operator opA(istlA);
        Dune::SeqScalarProduct<Vector> sp;
        // Right hand side.
        Vector istlb(opA.getmat().N());
        std::copy_n(b.data(), istlb.size(), istlb.begin());
        // System solution
        Vector x(opA.getmat().M());
        x = 0.0;

        // Construct preconditioner.
        // typedef Dune::SeqILU0<Mat,Vector,Vector> Preconditioner;
        typedef Opm::CPRPreconditioner<Mat,Vector,Vector> Preconditioner;
        const double relax = 1.0;
        Preconditioner precond(istlA, istlAe, relax);

        // Construct linear solver.
        const double tolerance = 1e-3;
        const int maxit = 5000;
        const int verbosity = 1;
        const int restart = 40;
        Dune::RestartedGMResSolver<Vector> linsolve(opA, sp, precond, tolerance, restart, maxit, verbosity);

        // Solve system.
        Dune::InverseOperatorResult result;
        linsolve.apply(x, istlb, result);

        // Check for failure of linear solver.
        if (!result.converged) {
            OPM_THROW(std::runtime_error, "Convergence failure for linear solver.");
        }

        // Copy solver output to dx.
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




        /// Form an elliptic system of equations.
        /// \param[in]       num_phases  the number of fluid phases
        /// \param[in]       eqs         the equations
        /// \param[out]      A           the resulting full system matrix
        /// \param[out]      b           the right hand side
        /// This function will deal with the first num_phases
        /// equations in eqs, and return a matrix A for the full
        /// system that has a elliptic upper left corner, if possible.
        void formEllipticSystem(const int num_phases,
                                const std::vector<ADB>& eqs_in,
                                Eigen::SparseMatrix<double, Eigen::RowMajor>& A,
                                V& b)
        {
            if (num_phases != 3) {
                OPM_THROW(std::logic_error, "formEllipticSystem() requires 3 phases.");
            }

            // A concession to MRST, to obtain more similar behaviour:
            // swap the first two equations, so that oil is first, then water.
            auto eqs = eqs_in;
            std::swap(eqs[0], eqs[1]);

            // Characterize the material balance equations.
            const int n = eqs[0].size();
            const double ratio_limit = 0.01;
            typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> Block;
            // The l1 block indicates if the equation for a given cell and phase is
            // sufficiently strong on the diagonal.
            Block l1 = Block::Zero(n, num_phases);
            for (int phase = 0; phase < num_phases; ++phase) {
                const M& J = eqs[phase].derivative()[0];
                V dj = J.diagonal().cwiseAbs();
                V sod = V::Zero(n);
                for (int elem = 0; elem < n; ++elem) {
                    sod(elem) = J.col(elem).cwiseAbs().sum() - dj(elem);
                }
                l1.col(phase) = (dj/sod > ratio_limit).cast<double>();
            }

            // By default, replace first equation with sum of all phase equations.
            // Build helper vectors.
            V l21 = V::Zero(n);
            V l22 = V::Ones(n);
            V l31 = V::Zero(n);
            V l33 = V::Ones(n);

            // If the first phase diagonal is not strong enough, we need further treatment.
            // Then the first equation will be the sum of the remaining equations,
            // and we swap the first equation into one of their slots.
            for (int elem = 0; elem < n; ++elem) {
                if (!l1(elem, 0)) {
                    const double l12x = l1(elem, 1);
                    const double l13x = l1(elem, 2);
                    const bool allzero = (l12x + l13x == 0);
                    if (allzero) {
                        l1(elem, 0) = 1;
                    } else {
                        if (l12x >= l13x) {
                            l21(elem) = 1;
                            l22(elem) = 0;
                        } else {
                            l31(elem) = 1;
                            l33(elem) = 0;
                        }
                    }
                }
            }

            // Construct the sparse matrix L that does the swaps and sums.
            Span i1(n, 1, 0);
            Span i2(n, 1, n);
            Span i3(n, 1, 2*n);
            std::vector< Eigen::Triplet<double> > t;
            t.reserve(7*n);
            for (int ii = 0; ii < n; ++ii) {
                t.emplace_back(i1[ii], i1[ii], l1(ii));
                t.emplace_back(i1[ii], i2[ii], l1(ii+n));
                t.emplace_back(i1[ii], i3[ii], l1(ii+2*n));
                t.emplace_back(i2[ii], i1[ii], l21(ii));
                t.emplace_back(i2[ii], i2[ii], l22(ii));
                t.emplace_back(i3[ii], i1[ii], l31(ii));
                t.emplace_back(i3[ii], i3[ii], l33(ii));
            }
            M L(3*n, 3*n);
            L.setFromTriplets(t.begin(), t.end());

            // Combine in single block.
            ADB total_residual = eqs[0];
            for (int phase = 1; phase < num_phases; ++phase) {
                total_residual = vertcat(total_residual, eqs[phase]);
            }
            total_residual = collapseJacs(total_residual);

            // Create output as product of L with equations.
            A = L * total_residual.derivative()[0];
            b = L * total_residual.value().matrix();
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

