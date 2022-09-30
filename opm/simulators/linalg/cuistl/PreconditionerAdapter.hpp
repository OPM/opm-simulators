/*
  Copyright SINTEF AS 2022

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
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/cuistl/CuSparseMatrix.hpp>
#include <opm/simulators/linalg/cuistl/PreconditionerHolder.hpp>
#include <opm/simulators/linalg/cuistl/impl/CuMatrixDescription.hpp>
#include <opm/simulators/linalg/cuistl/impl/CuSparseHandle.hpp>
#include <opm/simulators/linalg/cuistl/impl/CuSparseResource.hpp>
#include <opm/simulators/linalg/cuistl/impl/cusparse_constants.hpp>
#include <opm/simulators/linalg/cuistl/impl/cusparse_safe_call.hpp>
#include <opm/simulators/linalg/cuistl/impl/preconditioner_should_call_post_pre.hpp>


namespace Opm::cuistl
{
//!\brief Makes a CUDA preconditioner available to a CPU simulator.
template <class X, class Y, class CudaPreconditionerType>
class PreconditionerAdapter
    : public Dune::PreconditionerWithUpdate<X, Y>,
      public PreconditionerHolder<CuVector<typename X::field_type>, CuVector<typename Y::field_type>>
{
public:
    //! \brief The domain type of the preconditioner.
    using domain_type = X;
    //! \brief The range type of the preconditioner.
    using range_type = Y;
    //! \brief The field type of the preconditioner.
    using field_type = typename X::field_type;
    //! \brief scalar type underlying the field_type
    typedef Dune::SimdScalar<field_type> scalar_field_type;

    //! \brief Constructor.
    //!
    //! Constructor gets all parameters to operate the prec.
    //!    \param A The matrix to operate on.
    //!    \param w The relaxation factor.
    //!
    PreconditionerAdapter(std::shared_ptr<CudaPreconditionerType> preconditioner_)
        : m_underlyingPreconditioner(preconditioner_)
    {
    }


    //! \brief Prepare the preconditioner.
    //!
    //! \copydoc Preconditioner::pre(X&,Y&)
    virtual void pre([[maybe_unused]] X& x, [[maybe_unused]] Y& b) override
    {
        static_assert(!impl::shouldCallPreconditionerPre<CudaPreconditionerType>(),
                      "We currently do not support Preconditioner::pre().");
    }


    //! \brief Apply the preconditoner.
    //!
    //! \copydoc Preconditioner::apply(X&,const Y&)
    virtual void apply(X& v, const Y& d) override
    {
        if (!m_inputBuffer) {
            m_inputBuffer.reset(new CuVector<field_type>(v.dim()));
            m_outputBuffer.reset(new CuVector<field_type>(v.dim()));
        }
        m_inputBuffer->copyFromHost(d);
        m_underlyingPreconditioner->apply(*m_outputBuffer, *m_inputBuffer);
        m_outputBuffer->copyToHost(v);
    }


    //! \brief Clean up.
    //!
    //! \copydoc Preconditioner::post(X&)
    virtual void post([[maybe_unused]] X& x) override
    {
        static_assert(!impl::shouldCallPreconditionerPost<CudaPreconditionerType>(),
                      "We currently do not support Preconditioner::post().");
    }

    //! Category of the preconditioner (see SolverCategory::Category)
    virtual Dune::SolverCategory::Category category() const
    {
        return m_underlyingPreconditioner->category();
    }

    virtual void update() override
    {
        m_underlyingPreconditioner->update();
    }

    std::shared_ptr<Dune::PreconditionerWithUpdate<CuVector<field_type>, CuVector<field_type>>>
    getUnderlyingPreconditioner() override
    {
        return m_underlyingPreconditioner;
    }

    static constexpr bool shouldCallPre()
    {
        return impl::shouldCallPreconditionerPost<CudaPreconditionerType>();
    }
    static constexpr bool shouldCallPost()
    {
        return impl::shouldCallPreconditionerPre<CudaPreconditionerType>();
    }

private:
    //! \brief the underlying preconditioner to use
    std::shared_ptr<CudaPreconditionerType> m_underlyingPreconditioner;

    std::unique_ptr<CuVector<field_type>> m_inputBuffer;
    std::unique_ptr<CuVector<field_type>> m_outputBuffer;
};
} // end namespace Opm::cuistl

#endif
