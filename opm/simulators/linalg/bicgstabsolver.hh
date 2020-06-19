// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::Linear::BiCGStabSolver
 */
#ifndef EWOMS_BICG_STAB_SOLVER_HH
#define EWOMS_BICG_STAB_SOLVER_HH

#include "convergencecriterion.hh"
#include "residreductioncriterion.hh"
#include "linearsolverreport.hh"

#include <opm/models/utils/timer.hh>
#include <opm/models/utils/timerguard.hh>

#include <opm/material/common/Exceptions.hpp>

#include <memory>

namespace Opm {
namespace Linear {
/*!
 * \brief Implements a preconditioned stabilized BiCG linear solver.
 *
 * This solves a linear system of equations Ax = b, where the matrix A is sparse and may
 * be unsymmetric.
 *
 * See https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method, (article
 * date: December 19, 2016)
 */
template <class LinearOperator, class Vector, class Preconditioner>
class BiCGStabSolver
{
    using ConvergenceCriterion = Opm::Linear::ConvergenceCriterion<Vector>;
    using Scalar = typename LinearOperator::field_type;

public:
    BiCGStabSolver(Preconditioner& preconditioner,
                   ConvergenceCriterion& convergenceCriterion,
                   Dune::ScalarProduct<Vector>& scalarProduct)
        : preconditioner_(preconditioner)
        , convergenceCriterion_(convergenceCriterion)
        , scalarProduct_(scalarProduct)
    {
        A_ = nullptr;
        b_ = nullptr;

        maxIterations_ = 1000;
    }

    /*!
     * \brief Set the maximum number of iterations before we give up without achieving
     *        convergence.
     */
    void setMaxIterations(unsigned value)
    { maxIterations_ = value; }

    /*!
     * \brief Return the maximum number of iterations before we give up without achieving
     *        convergence.
     */
    unsigned maxIterations() const
    { return maxIterations_; }

    /*!
     * \brief Set the verbosity level of the linear solver
     *
     * The levels correspont to those used by the dune-istl solvers:
     *
     * - 0: no output
     * - 1: summary output at the end of the solution proceedure (if no exception was
     *      thrown)
     * - 2: detailed output after each iteration
     */
    void setVerbosity(unsigned value)
    { verbosity_ = value; }

    /*!
     * \brief Return the verbosity level of the linear solver.
     */
    unsigned verbosity() const
    { return verbosity_; }

    /*!
     * \brief Set the matrix "A" of the linear system.
     */
    void setLinearOperator(const LinearOperator* A)
    { A_ = A; }

    /*!
     * \brief Set the right hand side "b" of the linear system.
     */
    void setRhs(const Vector* b)
    { b_ = b; }

    /*!
     * \brief Run the stabilized BiCG solver and store the result into the "x" vector.
     */
    bool apply(Vector& x)
    {
        // epsilon used for detecting breakdowns
        const Scalar breakdownEps = std::numeric_limits<Scalar>::min() * Scalar(1e10);

        // start the stop watch for the solution proceedure, but make sure that it is
        // turned off regardless of how we leave the stadium. (i.e., that the timer gets
        // stopped in case exceptions are thrown as well as if the method returns
        // regularly.)
        report_.reset();
        Opm::TimerGuard reportTimerGuard(report_.timer());
        report_.timer().start();

        // preconditioned stabilized biconjugate gradient method
        //
        // See https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method,
        // (article date: December 19, 2016)

        // set the initial solution to the zero vector
        x = 0.0;

        // prepare the preconditioner. to allow some optimizations, we assume that the
        // preconditioner does not change the initial solution x if the initial solution
        // is a zero vector.
        Vector r = *b_;
        preconditioner_.pre(x, r);

#ifndef NDEBUG
        // ensure that the preconditioner does not change the initial solution. since
        // this is a debugging check, we don't care if it does not work properly in
        // parallel. (because this goes wrong, it should be considered to be a bug in the
        // code anyway.)
        for (unsigned i = 0; i < x.size(); ++i) {
            const auto& u = x[i];
            if (u*u != 0.0)
                throw std::logic_error("The preconditioner is assumed not to modify the initial solution!");
        }
#endif // NDEBUG

        convergenceCriterion_.setInitial(x, r);
        if (convergenceCriterion_.converged()) {
            report_.setConverged(true);
            return report_.converged();
        }

        if (verbosity_ > 0) {
            std::cout << "-------- BiCGStabSolver --------" << std::endl;
            convergenceCriterion_.printInitial();
        }

        // r0 = b - Ax (i.e., r -= A*x_0 = b, because x_0 == 0)
        //A_->applyscaleadd(/*alpha=*/-1.0, x, r);

        // r0hat = r0
        const Vector& r0hat = *b_;

        // rho0 = alpha = omega0 = 1
        Scalar rho = 1.0;
        Scalar alpha = 1.0;
        Scalar omega = 1.0;

        // v_0 = p_0 = 0;
        Vector v(r);
        v = 0.0;
        Vector p(v);

        // create all the temporary vectors which we need. Be aware that some of them
        // actually point to the same object because they are not needed at the same time!
        Vector y(x);
        Vector& h(x);
        Vector& s(r);
        Vector z(x);
        Vector& t(y);
        unsigned n = x.size();

        for (; report_.iterations() < maxIterations_; report_.increment()) {
            // rho_i = (r0hat,r_(i-1))
            Scalar rho_i = scalarProduct_.dot(r0hat, r);

            // beta = (rho_i/rho_(i-1))*(alpha/omega_(i-1))
            if (std::abs(rho) <= breakdownEps || std::abs(omega) <= breakdownEps)
                throw Opm::NumericalIssue("Breakdown of the BiCGStab solver (division by zero)");
            Scalar beta = (rho_i/rho)*(alpha/omega);

            // make rho correspond to the current iteration (i.e., forget rho_(i-1))
            rho = rho_i;

            // this loop conflates the following operations:
            //
            // p_i = r_(i-1) + beta*(p_(i-1) - omega_(i-1)*v_(i-1))
            // y = p
            for (unsigned i = 0; i < n; ++i) {
                // p_i = r_(i-1) + beta*(p_(i-1) - omega_(i-1)*v_(i-1))
                auto tmp = v[i];
                tmp *= omega;
                tmp -= p[i];
                tmp *= -beta;
                p[i] = r[i];
                p[i] += tmp;

                // y = p; not required because the precontioner overwrites y anyway...
                // y[i] = p[i];
            }

            // y = K^-1 * p_i
            preconditioner_.apply(y, p);

            // v_i = A*y
            A_->apply(y, v);

            // alpha = rho_i/(r0hat,v_i)
            Scalar denom = scalarProduct_.dot(r0hat, v);
            if (std::abs(denom) <= breakdownEps)
                throw Opm::NumericalIssue("Breakdown of the BiCGStab solver (division by zero)");
            alpha = rho_i/denom;
            if (std::abs(alpha) <= breakdownEps)
                throw Opm::NumericalIssue("Breakdown of the BiCGStab solver (stagnation detected)");

            // h = x_(i-1) + alpha*y
            // s = r_(i-1) - alpha*v_i
            for (unsigned i = 0; i < n; ++i) {
                auto tmp = y[i];
                tmp *= alpha;
                tmp += x[i];
                h[i] = tmp;

                //s[i] = r[i]; // not necessary because r and s are the same object
                tmp = v[i];
                tmp *= alpha;
                s[i] -= tmp;
            }

            // do convergence check and print terminal output
            convergenceCriterion_.update(/*curSol=*/h, /*delta=*/y, s);
            if (convergenceCriterion_.converged()) {
                if (verbosity_ > 0) {
                    convergenceCriterion_.print(report_.iterations() + 0.5);
                    std::cout << "-------- /BiCGStabSolver --------" << std::endl;
                }

                // x = h; // not necessary because x and h are the same object
                preconditioner_.post(x);
                report_.setConverged(true);
                return report_.converged();
            }
            else if (convergenceCriterion_.failed()) {
                if (verbosity_ > 0) {
                    convergenceCriterion_.print(report_.iterations() + 0.5);
                    std::cout << "-------- /BiCGStabSolver --------" << std::endl;
                }

                report_.setConverged(false);
                return report_.converged();
            }

            if (verbosity_ > 1)
                convergenceCriterion_.print(report_.iterations() + 0.5);

            // z = K^-1*s
            z = s;
            preconditioner_.apply(z, s);

            // t = Az
            t = z;
            A_->apply(z, t);

            // omega_i = (t*s)/(t*t)
            denom = scalarProduct_.dot(t, t);
            if (std::abs(denom) <= breakdownEps)
                throw Opm::NumericalIssue("Breakdown of the BiCGStab solver (division by zero)");
            omega = scalarProduct_.dot(t, s)/denom;
            if (std::abs(omega) <= breakdownEps)
                throw Opm::NumericalIssue("Breakdown of the BiCGStab solver (stagnation detected)");

            // x_i = h + omega_i*z
            // x = h; // not necessary because x and h are the same object
            x.axpy(/*a=*/omega, /*y=*/z);

            // do convergence check and print terminal output
            convergenceCriterion_.update(/*curSol=*/x, /*delta=*/z, r);
            if (convergenceCriterion_.converged()) {
                if (verbosity_ > 0) {
                    convergenceCriterion_.print(1.0 + report_.iterations());
                    std::cout << "-------- /BiCGStabSolver --------" << std::endl;
                }

                preconditioner_.post(x);
                report_.setConverged(true);
                return report_.converged();
            }
            else if (convergenceCriterion_.failed()) {
                if (verbosity_ > 0) {
                    convergenceCriterion_.print(1.0 + report_.iterations());
                    std::cout << "-------- /BiCGStabSolver --------" << std::endl;
                }

                report_.setConverged(false);
                return report_.converged();
            }

            if (verbosity_ > 1)
                convergenceCriterion_.print(1.0 + report_.iterations());

            // r_i = s - omega*t
            // r = s; // not necessary because r and s are the same object
            r.axpy(/*a=*/-omega, /*y=*/t);
        }

        report_.setConverged(false);
        return report_.converged();
    }

    void setConvergenceCriterion(ConvergenceCriterion& crit)
    {
        convergenceCriterion_ = &crit;
    }

    const Opm::Linear::SolverReport& report() const
    { return report_; }

private:
    const LinearOperator* A_;
    const Vector* b_;

    Preconditioner& preconditioner_;
    ConvergenceCriterion& convergenceCriterion_;
    Dune::ScalarProduct<Vector>& scalarProduct_;
    Opm::Linear::SolverReport report_;

    unsigned maxIterations_;
    unsigned verbosity_;
};

} // namespace Linear
} // namespace Opm

#endif
