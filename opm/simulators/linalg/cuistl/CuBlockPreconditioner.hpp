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
#ifndef OPM_CUISTL_CUBLOCKPRECONDITIONER_HPP
#define OPM_CUISTL_CUBLOCKPRECONDITIONER_HPP

#include <dune/common/shared_ptr.hh>
#include <memory>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/cuistl/PreconditionerHolder.hpp>
#include <opm/simulators/linalg/cuistl/detail/preconditioner_should_call_post_pre.hpp>

namespace Opm::cuistl
{
//! @brief Is an adaptation of Dune::BlockPreconditioner that works within the CuISTL framework.
//!
//! @note We aim to intgrate this into OwningBlockPreconditioner (or a relative thereof).
template <class X, class Y, class C, class P = Dune::PreconditionerWithUpdate<X, Y>>
class CuBlockPreconditioner : public Dune::PreconditionerWithUpdate<X, Y>, public PreconditionerHolder<X, Y>
{
public:
    using domain_type = X;
    using range_type = Y;
    using field_type = typename X::field_type;
    using communication_type = C;


    //! @brief Constructor.
    //!
    //! constructor gets all parameters to operate the prec.
    //! @param p The sequential preconditioner.
    //! @param c The communication object for syncing overlap and copy
    //! data points. (E.~g. OwnerOverlapCopyCommunication )
    //!
    CuBlockPreconditioner(const std::shared_ptr<P>& p, const std::shared_ptr<const communication_type>& c)
        : m_preconditioner(p)
        , m_communication(c)
    {
    }

    CuBlockPreconditioner(const std::shared_ptr<P>& p, const communication_type& c)
        : m_preconditioner(p)
        , m_communication(Dune::stackobject_to_shared_ptr(c))
    {
    }

    //! \brief Prepare the preconditioner.
    //!
    //! \copydoc Preconditioner::pre(X&,Y&)
    virtual void pre(X& x, Y& b) override
    {
        // TODO: [perf] Do we always need to copy, or only when we need the pre step?
        m_communication->copyOwnerToAll(x, x); // make Dirichlet values consistent
        if constexpr (detail::shouldCallPreconditionerPre<P>()) {
            m_preconditioner->pre(x, b);
        }
    }

    //! \brief Apply the preconditioner
    //!
    //! \copydoc Preconditioner::apply(X&,const Y&)
    virtual void apply(X& v, const Y& d) override
    {
        m_preconditioner->apply(v, d);
        m_communication->copyOwnerToAll(v, v);
    }


    virtual void update() override
    {
        m_preconditioner->update();
    }

    virtual void post(X& x) override
    {
        if constexpr (detail::shouldCallPreconditionerPost<P>()) {
            m_preconditioner->post(x);
        }
    }

    //! Category of the preconditioner (see SolverCategory::Category)
    virtual Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::overlapping;
    }

    static constexpr bool shouldCallPre()
    {
        return detail::shouldCallPreconditionerPost<P>();
    }
    static constexpr bool shouldCallPost()
    {
        return detail::shouldCallPreconditionerPre<P>();
    }

    virtual std::shared_ptr<Dune::PreconditionerWithUpdate<X, Y>> getUnderlyingPreconditioner() override
    {
        return m_preconditioner;
    }


private:
    //! \brief a sequential preconditioner
    std::shared_ptr<P> m_preconditioner;

    //! \brief the communication object
    std::shared_ptr<const communication_type> m_communication;
};
} // namespace Opm::cuistl
#endif
