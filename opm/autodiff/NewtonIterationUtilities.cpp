/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
  Copyright 2014 IRIS AS

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

#include <opm/autodiff/NewtonIterationUtilities.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/core/linalg/ParallelIstlInformation.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#if HAVE_UMFPACK
#include <Eigen/UmfPackSupport>
#else
#include <Eigen/SparseLU>
#endif
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

namespace Opm
{


    typedef AutoDiffBlock<double> ADB;
    typedef ADB::V V;
    typedef ADB::M M;
    typedef Eigen::SparseMatrix<double> S;


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
        // This is applied to all 2x2 block submatrices
        // The right hand side is modified accordingly. bi = bi - B * inv(D)* bn;
        // We do not explicitly compute inv(D) instead Du = C is solved

        // Extract the submatrix
        const std::vector<M>& Jn = eqs[n].derivative();

        // Use sparse LU to solve the block submatrices i.e compute inv(D)
        typedef Eigen::SparseMatrix<double> Sp;
        Sp Jnn;
        Jn[n].toSparse(Jnn);
#if HAVE_UMFPACK
        const Eigen::UmfPackLU<Sp> solver(Jnn);
#else
        const Eigen::SparseLU<Sp> solver(Jnn);
#endif
        Sp id(Jn[n].rows(), Jn[n].cols());
        id.setIdentity();
        const Sp Di = solver.solve(id);

        // compute inv(D)*bn for the update of the right hand side
        // Note: Eigen version > 3.2 requires a non-const reference to solve.
        ADB::V eqs_n_v = eqs[n].value();
        const Eigen::VectorXd& Dibn = solver.solve(eqs_n_v.matrix());

        std::vector<V> vals(num_eq);              // Number n will remain empty.
        std::vector<std::vector<M>> jacs(num_eq); // Number n will remain empty.
        for (int eq = 0; eq < num_eq; ++eq) {
            jacs[eq].reserve(num_eq - 1);
            const std::vector<M>& Je = eqs[eq].derivative();
            const M& B = Je[n];
            // Update right hand side.
            vals[eq] = eqs[eq].value().matrix() - B * Dibn;
        }
        for (int var = 0; var < num_eq; ++var) {
            if (var == n) {
                continue;
            }
            // solve Du = C
            // const M u = Di * Jn[var]; // solver.solve(Jn[var]);
            M u;
            fastSparseProduct(Di, Jn[var], u); // solver.solve(Jn[var]);
            for (int eq = 0; eq < num_eq; ++eq) {
                if (eq == n) {
                    continue;
                }
                const std::vector<M>& Je = eqs[eq].derivative();
                const M& B = Je[n];

                // Create new jacobians.
                // Add A
                jacs[eq].push_back(Je[var]);
                M& J = jacs[eq].back();
                // Subtract Bu (B*inv(D)*C)
                M Bu;
                fastSparseProduct(B, u, Bu);
                J = J + (Bu * -1.0);
            }
        }

        // Create return value.
        std::vector<ADB> retval;
        retval.reserve(num_eq - 1);
        for (int eq = 0; eq < num_eq; ++eq) {
            if (eq == n) {
                continue;
            }
            retval.push_back(ADB::function(std::move(vals[eq]), std::move(jacs[eq])));
        }
        return retval;
    }





    V recoverVariable(const ADB& equation, const V& partial_solution, const int n)
    {
        // The equation to solve for the unknown y (to be recovered) is
        //    Cx + Dy = b
        //    Dy = (b - Cx)
        // where D is the eliminated block, C is the jacobian of
        // the eliminated equation with respect to the
        // non-eliminated unknowms, b is the right-hand side of
        // the eliminated equation, and x is the partial solution
        // of the non-eliminated unknowns.

        const M& D1 = equation.derivative()[n];
        // Build C.
        std::vector<M> C_jacs = equation.derivative();
        C_jacs.erase(C_jacs.begin() + n);
        V equation_value = equation.value();
        ADB eq_coll = collapseJacs(ADB::function(std::move(equation_value), std::move(C_jacs)));
        const M& C = eq_coll.derivative()[0];

        // Use sparse LU to solve the block submatrices
        typedef Eigen::SparseMatrix<double> Sp;
        Sp D;
        D1.toSparse(D);
#if HAVE_UMFPACK
        const Eigen::UmfPackLU<Sp> solver(D);
#else
        const Eigen::SparseLU<Sp> solver(D);
#endif

        // Compute value of eliminated variable.
        const Eigen::VectorXd b = (equation.value().matrix() - C * partial_solution.matrix());
        const Eigen::VectorXd elim_var = solver.solve(b);

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
        eqs[0].swap(eqs[1]);

        // Characterize the material balance equations.
        const int n = eqs[0].size();
        const double ratio_limit = 0.01;
        typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> Block;
        // The l1 block indicates if the equation for a given cell and phase is
        // sufficiently strong on the diagonal.
        Block l1 = Block::Zero(n, num_phases);
        {
            S J;
            for (int phase = 0; phase < num_phases; ++phase) {
                eqs[phase].derivative()[0].toSparse(J);
                V dj = J.diagonal().cwiseAbs();
                V sod = V::Zero(n);
                for (int elem = 0; elem < n; ++elem) {
                    sod(elem) = J.col(elem).cwiseAbs().sum() - dj(elem);
                }
                l1.col(phase) = (dj/sod > ratio_limit).cast<double>();
            }
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
        S L(3*n, 3*n);
        L.setFromTriplets(t.begin(), t.end());

        // Combine in single block.
        ADB total_residual = vertcatCollapseJacs(eqs);

        S derivative;
        total_residual.derivative()[0].toSparse(derivative);

        // Create output as product of L with equations.
        A = L * derivative;
        b = L * total_residual.value().matrix();
    }




    /// Return true if this is a serial run, or rank zero on an MPI run.
    bool isIORank(const boost::any& parallel_info)
    {
#if HAVE_MPI
        if (parallel_info.type() == typeid(ParallelISTLInformation)) {
            const ParallelISTLInformation& info =
                boost::any_cast<const ParallelISTLInformation&>(parallel_info);
            return info.communicator().rank() == 0;
        } else {
            return true;
        }
#else
        static_cast<void>(parallel_info); // Suppress unused argument warning.
        return true;
#endif
    }


} // namespace Opm

