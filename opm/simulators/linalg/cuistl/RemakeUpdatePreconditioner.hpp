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
#ifndef OPM_REMAKEUPDATEPRECONDITIONER_HEADER_INCLUDED
#define OPM_REMAKEUPDATEPRECONDITIONER_HEADER_INCLUDED
#include <cusparse.h>
#include <dune/common/simd.hh>
#include <dune/istl/preconditioner.hh>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/cuistl/CuSparseMatrix.hpp>
#include <opm/simulators/linalg/cuistl/impl/CuMatrixDescription.hpp>
#include <opm/simulators/linalg/cuistl/impl/CuSparseHandle.hpp>
#include <opm/simulators/linalg/cuistl/impl/CuSparseResource.hpp>
#include <opm/simulators/linalg/cuistl/impl/cusparse_constants.hpp>
#include <opm/simulators/linalg/cuistl/impl/cusparse_safe_call.hpp>


namespace Opm::cuistl
{
//! @brief Enables updating on a preconditioner that does not have an update function (but needs one)
//!
//! This works by re-creating the underlying preconditioner at every update call.
//!
//! @note This is meant mostly for debugging and testing older Dune preconditioners
//!
//! @note This is currently not in use.
template <class M, class X, class Y, int l = 1>
class RemakeUpdatePreconditioner : public Dune::PreconditionerWithUpdate<X, Y>
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

    //!\brief Constructor.
    //!
    //! Constructor gets all parameters to operate the prec.
    //! \param creator A functor to create the preconditioner (the functor needs to hold a reference to a matrix if
    //! needed)
    RemakeUpdatePreconditioner(std::function<std::unique_ptr<Dune::Preconditioner<X, Y>>()> creator)
        : m_creator(creator)
        , m_underlyingPreconditioner(creator())
    {
    }

    //! \brief Prepare the preconditioner.
    virtual void pre(X& x, Y& b) override
    {
        m_underlyingPreconditioner->pre(x, b);
    }

    //! \brief Apply the preconditoner.
    virtual void apply(X& v, const Y& d) override
    {

        m_underlyingPreconditioner->apply(v, d);
    }

    //! \brief Clean up.
    virtual void post(X& x) override
    {
        m_underlyingPreconditioner->post(x);
    }

    //! Category of the preconditioner (see SolverCategory::Category)
    virtual Dune::SolverCategory::Category category() const
    {
        return m_underlyingPreconditioner->category();
    }

    virtual void update() override
    {
        m_underlyingPreconditioner = m_creator();
    }

private:
    std::function<std::unique_ptr<Dune::Preconditioner<X, Y>>()> m_creator;
    //! \brief the underlying preconditioner to use
    std::unique_ptr<Dune::Preconditioner<X, Y>> m_underlyingPreconditioner;
};
} // end namespace Opm::cuistl

#endif
