/*
  Copyright SINTEF AS

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
#ifndef OPM_PRECONDITIONERADAPTER_HEADER_INCLUDED
#define OPM_PRECONDITIONERADAPTER_HEADER_INCLUDED
#include <cusparse.h>
#include <dune/common/simd.hh>
#include <dune/istl/preconditioner.hh>
#include <opm/simulators/linalg/cuistl/CuMatrixDescription.hpp>
#include <opm/simulators/linalg/cuistl/CuSparseHandle.hpp>
#include <opm/simulators/linalg/cuistl/CuSparseMatrix.hpp>
#include <opm/simulators/linalg/cuistl/CuSparseResource.hpp>
#include <opm/simulators/linalg/cuistl/cusparse_constants.hpp>
#include <opm/simulators/linalg/cuistl/cusparse_safe_call.hpp>


namespace Opm::cuistl
{
//!\brief Makes a CUDA preconditioner available to a CPU simulator.
template <class CudaPreconditionerType, class M, class X, class Y, int l = 1>
class PreconditionerAdapter : public Dune::Preconditioner<X, Y>
{
public:
    //! \brief The matrix type the preconditioner is for.
    typedef typename std::remove_const<M>::type matrix_type;
    //! \brief The domain type of the preconditioner.
    typedef X domain_type;
    //! \brief The range type of the preconditioner.
    typedef Y range_type;
    //! \brief The field type of the preconditioner.
    typedef typename X::field_type field_type;
    //! \brief scalar type underlying the field_type
    typedef Dune::SimdScalar<field_type> scalar_field_type;

    /*! \brief Constructor.

    Constructor gets all parameters to operate the prec.
       \param A The matrix to operate on.
       \param w The relaxation factor.
            */
    PreconditionerAdapter(std::shared_ptr<CudaPreconditionerType> preconditioner_)
        : underlyingPeconditioner(preconditioner_)
    {

    }

    /*!
       \brief Prepare the preconditioner.

    \copydoc Preconditioner::pre(X&,Y&)
        */
    virtual void pre(X& x, Y& b)
    {
        DUNE_UNUSED_PARAMETER(x);
        DUNE_UNUSED_PARAMETER(b);
    }

    /*!
       \brief Apply the preconditoner.

    \copydoc Preconditioner::apply(X&,const Y&)
        */
    virtual void apply(X& v, const Y& d)
    {
        if (!inputBuffer) {
            inputBuffer.reset(new CuVector<field_type>(v.dim() * v.N()));
            outputBuffer.reset(new CuVector<field_type>(v.dim() * v.N()));
        }
        inputBuffer->copyFrom(d);
        underlyingPeconditioner->apply(*inputBuffer, *outputBuffer);
        outputBuffer->copyTo(v);
    }

    /*!
       \brief Clean up.

    \copydoc Preconditioner::post(X&)
        */
    virtual void post(X& x)
    {
        DUNE_UNUSED_PARAMETER(x);
    }

    //! Category of the preconditioner (see SolverCategory::Category)
    virtual Dune::SolverCategory::Category category() const
    {
        return Dune::SolverCategory::sequential;
    }

private:
    //! \brief the underlying preconditioner to use
    std::shared_ptr<CudaPreconditionerType> underlyingPeconditioner;

    std::unique_ptr<CuVector<field_type>> inputBuffer;
    std::unique_ptr<CuVector<field_type>> outputBuffer;
};
} // end namespace Opm::cuistl

#endif
