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
 * \copydoc Opm::Linear::OverlappingScalarProduct
 */
#ifndef EWOMS_OVERLAPPING_SCALAR_PRODUCT_HH
#define EWOMS_OVERLAPPING_SCALAR_PRODUCT_HH

#include <dune/common/version.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/istl/scalarproducts.hh>

namespace Opm {
namespace Linear {

/*!
 * \brief An overlap aware ISTL scalar product.
 */
template <class OverlappingBlockVector, class Overlap>
class OverlappingScalarProduct
    : public Dune::ScalarProduct<OverlappingBlockVector>
{
public:
    using field_type = typename OverlappingBlockVector::field_type;
    using CollectiveCommunication = typename Dune::CollectiveCommunication<typename Dune::MPIHelper::MPICommunicator>;

    using real_type = typename Dune::ScalarProduct<OverlappingBlockVector>::real_type;

    //! the kind of computations supported by the operator. Either overlapping or non-overlapping
    Dune::SolverCategory::Category category() const override
    { return Dune::SolverCategory::overlapping; }

    OverlappingScalarProduct(const Overlap& overlap)
        : overlap_(overlap), comm_( Dune::MPIHelper::getCollectiveCommunication() )
    {}

#if DUNE_VERSION_NEWER(DUNE_ISTL, 2,7)
    field_type dot(const OverlappingBlockVector& x,
                   const OverlappingBlockVector& y) const override
#else
    field_type dot(const OverlappingBlockVector& x,
                   const OverlappingBlockVector& y) override
#endif
    {
        field_type sum = 0;
        size_t numLocal = overlap_.numLocal();
        for (unsigned localIdx = 0; localIdx < numLocal; ++localIdx) {
            if (overlap_.iAmMasterOf(static_cast<int>(localIdx)))
                sum += x[localIdx] * y[localIdx];
        }

        // return the global sum
        return comm_.sum( sum );
    }

#if DUNE_VERSION_NEWER(DUNE_ISTL, 2,7)
    real_type norm(const OverlappingBlockVector& x) const override
#else
    real_type norm(const OverlappingBlockVector& x) override
#endif
    { return std::sqrt(dot(x, x)); }

private:
    const Overlap& overlap_;
    const CollectiveCommunication comm_;
};

} // namespace Linear
} // namespace Opm

#endif
