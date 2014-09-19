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

#ifndef OPM_CPRPRECONDITIONER_HEADER_INCLUDED
#define OPM_CPRPRECONDITIONER_HEADER_INCLUDED


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

#include <opm/core/utility/ErrorMacros.hpp>

namespace Opm
{
    /*!
      \brief Sequential CPR preconditioner.

      This is a two-stage preconditioner, combining an elliptic-type
      partial solution with ILU0 for the whole system.

      \tparam M The matrix type to operate on
      \tparam X Type of the update
      \tparam Y Type of the defect
    */
    template<class M, class X, class Y>
    class CPRPreconditioner : public Dune::Preconditioner<X,Y> {
    public:
        //! \brief The matrix type the preconditioner is for.
        typedef typename Dune::remove_const<M>::type matrix_type;
        //! \brief The domain type of the preconditioner.
        typedef X domain_type;
        //! \brief The range type of the preconditioner.
        typedef Y range_type;
        //! \brief The field type of the preconditioner.
        typedef typename X::field_type field_type;

        // define the category
        enum {
            //! \brief The category the preconditioner is part of.
            category = Dune::SolverCategory::sequential
        };

        /*! \brief Constructor.

          Constructor gets all parameters to operate the prec.
          \param A  The matrix to operate on.
          \param Ae The top-left elliptic part of A.
          \param w  The ILU0 relaxation factor.
          \param amg  If true, use AMG preconditioner for elliptic part.
        */
        CPRPreconditioner (const M& A, const M& Ae, const field_type relax, const bool amg)
            : A_(A),
              ILU_(A), // copy A (will be overwritten by ILU decomp)
              Ae_(Ae),
              relax_(relax),
              use_amg_(amg)
        {
            Dune::bilu0_decomposition(ILU_);
        }

        /*!
          \brief Prepare the preconditioner.

          \copydoc Preconditioner::pre(X&,Y&)
        */
        virtual void pre (X& /*x*/, Y& /*b*/)
        {
        }

        /*!
          \brief Apply the preconditoner.

          \copydoc Preconditioner::apply(X&,const Y&)
        */
        virtual void apply (X& v, const Y& d)
        {
            // Extract part of d corresponding to elliptic part.
            Y de(Ae_.N());
            // Note: Assumes that the elliptic part comes first.
            std::copy_n(d.begin(), Ae_.N(), de.begin());

            // Solve elliptic part, extend solution to full.
            Y ve = solveElliptic(de);
            Y vfull(ILU_.N());
            vfull = 0.0;
            // Again assuming that the elliptic part comes first.
            std::copy(ve.begin(), ve.end(), vfull.begin());

            // Subtract elliptic residual from initial residual.
            // dmodified = d - A * vfull
            Y dmodified = d;
            A_.mmv(vfull, dmodified);

            // Apply ILU0.
            Y vilu(ILU_.N());
            Dune::bilu_backsolve(ILU_, vilu, dmodified);
            v = vfull;
            v += vilu;
            v *= relax_;
        }

        /*!
          \brief Clean up.

          \copydoc Preconditioner::post(X&)
        */
        virtual void post (X& /*x*/)
        {
        }

    private:
        Y solveElliptic(Y& de)
        {
            // std::cout << "solveElliptic()" << std::endl;
            // Construct operator, scalar product and vectors needed.
            typedef Dune::MatrixAdapter<M,X,X> Operator;
            Operator opAe(Ae_);
            Dune::SeqScalarProduct<X> sp;
            // Right hand side.
            // System solution
            X x(opAe.getmat().M());
            x = 0.0;

            Dune::InverseOperatorResult result;
            if (use_amg_) {
                typedef Dune::Amg::FirstDiagonal CouplingMetric;
                // typedef Dune::Amg::RowSum        CouplingMetric;
                // typedef Dune::Amg::AggregationCriterion<Dune::Amg::SymmetricMatrixDependency<M,CouplingMetric> > CriterionBase;
                typedef Dune::Amg::UnSymmetricCriterion<M,CouplingMetric> CriterionBase;
                typedef Dune::Amg::CoarsenCriterion<CriterionBase> Criterion;
                // typedef Dune::SeqILU0<M,X,X>        Smoother;
                typedef Dune::SeqSOR<M,X,X>        Smoother;

                typedef Dune::Amg::AMG<Operator,X,Smoother>   Precond;

                // Construct preconditioner.
                Criterion criterion;
                const double linsolver_prolongate_factor = 1.6;
                const int verbosity = 0;
                const int smooth_steps = 1;
                criterion.setDebugLevel(verbosity);
                criterion.setDefaultValuesAnisotropic(3, 2);
                criterion.setProlongationDampingFactor(linsolver_prolongate_factor);
                criterion.setNoPreSmoothSteps(smooth_steps);
                criterion.setNoPostSmoothSteps(smooth_steps);
                criterion.setGamma(1); // V-cycle; this is the default
                typename Precond::SmootherArgs smootherArgs;
                Precond precond(opAe, criterion, smootherArgs);

                // Construct linear solver.
                const double tolerance = 1e-4;
                const int maxit = 5000;
                Dune::BiCGSTABSolver<X> linsolve(opAe, sp, precond, tolerance, maxit, verbosity);
                // const int restart = 40;
                // Dune::RestartedGMResSolver<X> linsolve(opAe, sp, precond, tolerance, restart, maxit, verbosity);
                // Solve system.
                linsolve.apply(x, de, result);
            } else {
                // Construct preconditioner.
                typedef typename Dune::SeqILU0<M,X,X> Preconditioner;
                const double relax = 1.0;
                Preconditioner precond(Ae_, relax);

                // Construct linear solver.
                const double tolerance = 1e-4;
                const int maxit = 5000;
                const int verbosity = 0;
                Dune::BiCGSTABSolver<X> linsolve(opAe, sp, precond, tolerance, maxit, verbosity);
                // Solve system.
                linsolve.apply(x, de, result);
            }
            if (!result.converged) {
                OPM_THROW(std::runtime_error, "CPRPreconditioner failed to solve elliptic subsystem.");
            }
            return x;
        }

        //! \brief The matrix for the full linear problem.
        const matrix_type& A_;
        //! \brief The ILU0 decomposition of the matrix.
        matrix_type ILU_;
        //! \brief The elliptic part of the matrix.
        matrix_type Ae_;
        //! \brief The relaxation factor to use.
        field_type relax_;
        //! \brief If true, use an AMG preconditioner for the elliptic part.
        bool use_amg_;
    };

} // namespace Opm

#endif // OPM_CPRPRECONDITIONER_HEADER_INCLUDED
