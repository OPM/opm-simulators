/*
  Copyright 2022-2023 SINTEF AS

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
#ifndef OPM_PRECONDITIONERADAPTER_HPP
#define OPM_PRECONDITIONERADAPTER_HPP
#include <cusparse.h>
#include <dune/istl/preconditioner.hh>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>
#include <opm/simulators/linalg/gpuistl/PreconditionerHolder.hpp>
#include <opm/simulators/linalg/gpuistl/detail/preconditioner_should_call_post_pre.hpp>


namespace Opm::gpuistl
{
//!\brief Makes a CUDA preconditioner available to a CPU simulator.
//!
//! The use case for this adapter is to use a CUDA preconditioner during a linear
//! solver that works on the CPU. The motivation for this is benchmarking new preconditioners on the GPU.
//!
//! \tparam X the domain type (should be on the CPU). Typicall a Dune::BlockVector
//! \tparam Y the range type (should be on the CPU). Typicall a Dune::BlockVector
//! \tparam CudaPreconditionerType the preconditioner taking GpuVector<real_type> as arguments to apply
template <class X, class Y, class CudaPreconditionerType>
class PreconditionerAdapter
    : public Dune::PreconditionerWithUpdate<X, Y>,
      public PreconditionerHolder<GpuVector<typename X::field_type>, GpuVector<typename Y::field_type>>
{
public:
    //! \brief The domain type of the preconditioner.
    using domain_type = X;
    //! \brief The range type of the preconditioner.
    using range_type = Y;
    //! \brief The field type of the preconditioner.
    using field_type = typename X::field_type;

    //! \brief Constructor.
    //!
    //! Constructor gets all parameters to operate the prec.
    //!    \param preconditioner Preconditioner to use in adapter
    //!
    explicit PreconditionerAdapter(std::shared_ptr<CudaPreconditionerType> preconditioner)
        : m_underlyingPreconditioner(preconditioner)
    {
    }


    //! \brief Prepare the preconditioner.
    //!
    //! Currently not supported.
    virtual void pre([[maybe_unused]] X& x, [[maybe_unused]] Y& b) override
    {
        static_assert(!detail::shouldCallPreconditionerPre<CudaPreconditionerType>(),
                      "We currently do not support Preconditioner::pre().");
    }


    //! \brief Apply the preconditoner.
    //!
    //! \copydoc Preconditioner::apply(X&,const Y&)
    virtual void apply(X& v, const Y& d) override
    {
        if (!m_inputBuffer) {
            m_inputBuffer.reset(new GpuVector<field_type>(v.dim()));
            m_outputBuffer.reset(new GpuVector<field_type>(v.dim()));
        }
        m_inputBuffer->copyFromHost(d);
        m_underlyingPreconditioner->apply(*m_outputBuffer, *m_inputBuffer);
        m_outputBuffer->copyToHost(v);
    }


    //! \brief Clean up.
    //!
    //! Currently not supported.
    virtual void post([[maybe_unused]] X& x) override
    {
        static_assert(!detail::shouldCallPreconditionerPost<CudaPreconditionerType>(),
                      "We currently do not support Preconditioner::post().");
    }


    //! Category of the preconditioner (see SolverCategory::Category)
    Dune::SolverCategory::Category category() const override
    {
        return m_underlyingPreconditioner->category();
    }

    //! Calls update on the underlying CUDA preconditioner
    virtual void update() override
    {
        m_underlyingPreconditioner->update();
    }

    static constexpr bool shouldCallPre()
    {
        return detail::shouldCallPreconditionerPost<CudaPreconditionerType>();
    }
    static constexpr bool shouldCallPost()
    {
        return detail::shouldCallPreconditionerPre<CudaPreconditionerType>();
    }

    virtual std::shared_ptr<Dune::PreconditionerWithUpdate<GpuVector<field_type>, GpuVector<field_type>>>
    getUnderlyingPreconditioner() override
    {
        return m_underlyingPreconditioner;
    }

    virtual bool hasPerfectUpdate() const override {
        return m_underlyingPreconditioner->hasPerfectUpdate();
    }

private:
    //! \brief the underlying preconditioner to use
    std::shared_ptr<CudaPreconditionerType> m_underlyingPreconditioner;

    std::unique_ptr<GpuVector<field_type>> m_inputBuffer;
    std::unique_ptr<GpuVector<field_type>> m_outputBuffer;
};
} // end namespace Opm::gpuistl

#endif
