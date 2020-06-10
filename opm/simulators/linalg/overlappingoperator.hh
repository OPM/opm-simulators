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
 * \copydoc Opm::Linear::OverlappingOperator
 */
#ifndef EWOMS_OVERLAPPING_OPERATOR_HH
#define EWOMS_OVERLAPPING_OPERATOR_HH

#include <dune/istl/operators.hh>
#include <dune/common/version.hh>

namespace Opm {
namespace Linear {

/*!
 * \brief An overlap aware linear operator usable by ISTL.
 */
template <class OverlappingMatrix, class DomainVector, class RangeVector>
class OverlappingOperator
    : public Dune::AssembledLinearOperator<OverlappingMatrix, DomainVector, RangeVector>
{
    using Overlap = typename OverlappingMatrix::Overlap;

public:
    //! export types
    using domain_type = DomainVector;
    using field_type = typename domain_type::field_type;

    OverlappingOperator(const OverlappingMatrix& A) : A_(A)
    {}

    //! the kind of computations supported by the operator. Either overlapping or non-overlapping
    Dune::SolverCategory::Category category() const override
    { return Dune::SolverCategory::overlapping; }

    //! apply operator to x:  \f$ y = A(x) \f$
    virtual void apply(const DomainVector& x, RangeVector& y) const override
    {
        A_.mv(x, y);
        y.sync();
    }

    //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
    virtual void applyscaleadd(field_type alpha, const DomainVector& x,
                               RangeVector& y) const override
    {
        A_.usmv(alpha, x, y);
        y.sync();
    }

    //! returns the matrix
    virtual const OverlappingMatrix& getmat() const override
    { return A_; }

    const Overlap& overlap() const
    { return A_.overlap(); }

private:
    const OverlappingMatrix& A_;
};

} // namespace Linear
} // namespace Opm

#endif
