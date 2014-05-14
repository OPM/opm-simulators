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
          \param ne The size of the elliptic top-left part
          \param w  The ILU0 relaxation factor.
        */
        CPRPreconditioner (const M& A, const int ne, const field_type w)
            : ILU(A) // copy A
        {
            _w =w;
            Dune::bilu0_decomposition(ILU);
            // Ae = A(pindx, pindx);
        }

        /*!
          \brief Prepare the preconditioner.

          \copydoc Preconditioner::pre(X&,Y&)
        */
        virtual void pre (X& x, Y& b)
        {
            DUNE_UNUSED_PARAMETER(x);
            DUNE_UNUSED_PARAMETER(b);
        }

        /*!
          \brief Apply the preconditoner.

          \copydoc Preconditioner::apply(X&,const Y&)
        */
        virtual void apply (X& v, const Y& d)
        {
            Dune::bilu_backsolve(ILU,v,d);
            v *= _w;
        }

        /*!
          \brief Clean up.

          \copydoc Preconditioner::post(X&)
        */
        virtual void post (X& x)
        {
            DUNE_UNUSED_PARAMETER(x);
        }

    private:
        //! \brief The relaxation factor to use.
        field_type _w;
        //! \brief The ILU0 decomposition of the matrix.
        matrix_type ILU;
        //! \brief The elliptic part of the matrix.
        matrix_type Ae;
    };

} // namespace Opm

#endif // OPM_CPRPRECONDITIONER_HEADER_INCLUDED
