/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
  Copyright 2014 IRIS AS.

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

#include <memory>

#include <opm/core/utility/platform_dependent/disable_warnings.h>

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

#include <opm/core/utility/platform_dependent/reenable_warnings.h>

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
    class CPRPreconditioner : public Dune::Preconditioner<X,Y>
    {
        // prohibit copying for now
        CPRPreconditioner( const CPRPreconditioner& );

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

        //! \brief Elliptic Operator
        typedef Dune::MatrixAdapter<M,X,X> Operator;

        //! \brief ilu-0 preconditioner for the elliptic system
        typedef Dune::SeqILU0<M,X,X> Preconditioner;

        //! \brief amg preconditioner for the elliptic system
        typedef Preconditioner Smoother;
        typedef Dune::Amg::AMG<Operator, X, Smoother> AMG;

        /*! \brief Constructor.

          Constructor gets all parameters to operate the prec.
          \param A       The matrix to operate on.
          \param Ae      The top-left elliptic part of A.
          \param relax   The ILU0 relaxation factor.
          \param useAMG  if true, AMG is used as a preconditioner for the elliptic sub-system, otherwise ilu-0 (default)
          \param useBiCG if true, BiCG solver is used (default), otherwise CG solver
        */
        CPRPreconditioner (const M& A, const M& Ae, const field_type relax,
                           const bool useAMG  = false,
                           const bool useBiCG = true )
            : A_(A),
              Ae_(Ae),
              de_( Ae_.N() ),
              ve_( Ae_.M() ),
              dmodified_( A_.N() ),
              opAe_( Ae_ ),
              precond_(), // ilu0 preconditioner for elliptic system
              amg_(),     // amg  preconditioner for elliptic system
              ILU_(A),    // copy A (will be overwritten by ILU decomp)
              vilu_( ILU_.N() ),
              relax_(relax),
              use_bicg_solver_( useBiCG )
        {
            // create appropriate preconditioner for elliptic system
            createPreconditioner( useAMG );

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
            // Note: Assumes that the elliptic part comes first.
            std::copy_n(d.begin(), de_.size(), de_.begin());

            // Solve elliptic part, extend solution to full.
            // reset result
            ve_ = 0;
            solveElliptic( ve_, de_ );

            //reset return value
            v = 0.0;
            // Again assuming that the elliptic part comes first.
            std::copy(ve_.begin(), ve_.end(), v.begin());

            // Subtract elliptic residual from initial residual.
            // dmodified = d - A * vfull
            dmodified_ = d;
            A_.mmv(v, dmodified_);

            // Apply ILU0.
            Dune::bilu_backsolve(ILU_, vilu_, dmodified_);
            v += vilu_;

            // don't apply relaxation if relax_ == 1 
            if( std::abs( relax_ - 1.0 ) < 1e-12 ) return;

            v *= relax_;
        }

        /*!
          \brief Clean up.

          \copydoc Preconditioner::post(X&)
        */
        virtual void post (X& /*x*/)
        {
        }

     protected:
        void solveElliptic(Y& x, Y& de)
        {
            // Linear solver parameters
            const double tolerance = 1e-4;
            const int maxit = 5000;
            const int verbosity = 0;

            // operator result containing iterations etc.
            Dune::InverseOperatorResult result;

            // sequential scalar product
            Dune::SeqScalarProduct<X> sp;
            if( amg_ )
            {
                // Solve system with AMG
                if( use_bicg_solver_ ) {
                    Dune::BiCGSTABSolver<X> linsolve(opAe_, sp, (*amg_), tolerance, maxit, verbosity);
                    linsolve.apply(x, de, result);
                }
                else {
                    Dune::CGSolver<X> linsolve(opAe_, sp, (*amg_), tolerance, maxit, verbosity);
                    linsolve.apply(x, de, result);
                }
            }
            else
            {
                assert( precond_ );
                // Solve system with ILU-0
                if( use_bicg_solver_ ) {
                    Dune::BiCGSTABSolver<X> linsolve(opAe_, sp, (*precond_), tolerance, maxit, verbosity);
                    linsolve.apply(x, de, result);
                }
                else {
                    Dune::CGSolver<X> linsolve(opAe_, sp, (*precond_), tolerance, maxit, verbosity);
                    linsolve.apply(x, de, result);
                }

            }

            if (!result.converged) {
                OPM_THROW(std::runtime_error, "CPRPreconditioner failed to solve elliptic subsystem.");
            }
        }

        //! \brief The matrix for the full linear problem.
        const matrix_type& A_;
        //! \brief The elliptic part of the matrix.
        const matrix_type& Ae_;

        //! \brief temporary variables for elliptic solve
        Y de_, ve_, dmodified_;

        //! \brief elliptic operator
        Operator opAe_;

        //! \brief ILU0 preconditioner for the elliptic system
        std::unique_ptr< Preconditioner > precond_;
        //! \brief AMG preconditioner with ILU0 smoother
        std::unique_ptr< AMG > amg_;

        //! \brief The ILU0 decomposition of the matrix.
        matrix_type ILU_;

        //! \brief temporary variables for ILU solve
        Y vilu_;

        //! \brief The relaxation factor to use.
        field_type relax_;

        //! \brief true if ISTL BiCGSTABSolver is used, otherwise ISTL CGSolver is used
        const bool use_bicg_solver_;

     protected:
        void createPreconditioner( const bool amg )
        {
            if( amg )
            {
              typedef Dune::Amg::CoarsenCriterion< Dune::Amg::SymmetricCriterion<M, Dune::Amg::FirstDiagonal> > Criterion;
              typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments SmootherArgs;

              SmootherArgs smootherArgs;

              smootherArgs.iterations = 1;
              smootherArgs.relaxationFactor = relax_;

              int coarsenTarget=1200;
              Criterion criterion(15,coarsenTarget);
              criterion.setDebugLevel( 0 ); // no debug information, 1 for printing hierarchy information
              criterion.setDefaultValuesIsotropic(2);
              criterion.setAlpha(.67);
              criterion.setBeta(1.0e-6);
              criterion.setMaxLevel(10);
              amg_ = std::unique_ptr< AMG > (new AMG(opAe_, criterion, smootherArgs));
            }
            else
              precond_ = std::unique_ptr< Preconditioner > (new Preconditioner( Ae_, relax_ ));
       }
    };

} // namespace Opm

#endif // OPM_CPRPRECONDITIONER_HEADER_INCLUDED
