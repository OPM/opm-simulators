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
 * \copydoc PffGridVector
 */
#ifndef EWOMS_PFF_GRID_VECTOR_HH
#define EWOMS_PFF_GRID_VECTOR_HH

#include <opm/models/utils/prefetch.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/common/version.hh>

#include <vector>

namespace Opm {
/*!
 * \brief A random-access container which stores data attached to a grid's degrees of
 *        freedom in a prefetch friendly manner.
 *
 * This container often reduces the number of cache faults considerably, thus improving
 * performance. On the flipside data cannot be written to on an individual basis and it
 * requires significantly more memory than a plain array. PffVector stands for "PreFetch
 * Friendly Grid Vector".
 */
template <class GridView, class Stencil, class Data, class DofMapper>
class PffGridVector
{
    typedef typename GridView::template Codim<0>::Entity Element;

#if DUNE_VERSION_NEWER(DUNE_GRID, 2,6)
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView> ElementMapper;
#else
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, Dune::MCMGElementLayout> ElementMapper;
#endif

public:
    PffGridVector(const GridView& gridView, const DofMapper& dofMapper)
        : gridView_(gridView)
#if DUNE_VERSION_NEWER(DUNE_GRID, 2,6)
        , elementMapper_(gridView_, Dune::mcmgElementLayout())
#else
        , elementMapper_(gridView_)
#endif
        , dofMapper_(dofMapper)
    { }

    template <class DistFn>
    void update(const DistFn& distFn)
    {
        unsigned numElements = gridView_.size(/*codim=*/0);
        unsigned numLocalDofs = computeNumLocalDofs_();

        elemData_.resize(numElements);
        data_.resize(numLocalDofs);

        // update the pointers for the element data: for this, we need to loop over the
        // whole grid and update a stencil for each element
        Data *curElemDataPtr = &data_[0];
        Stencil stencil(gridView_, dofMapper_);
        auto elemIt = gridView_.template begin</*codim=*/0>();
        const auto& elemEndIt = gridView_.template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            // set the DOF data pointer for the current element
            const auto& elem = *elemIt;
            unsigned elemIdx = elementMapper_.index(elem);
            elemData_[elemIdx] = curElemDataPtr;

            stencil.update(elem);
            unsigned numDof = stencil.numDof();
            for (unsigned localDofIdx = 0; localDofIdx < numDof; ++ localDofIdx)
                distFn(curElemDataPtr[localDofIdx], stencil, localDofIdx);

            // update the element data pointer to make it point to the beginning of the
            // data for DOFs of the next element
            curElemDataPtr += numDof;
        }
    }

    void prefetch(const Element& elem) const
    {
        unsigned elemIdx = elementMapper_.index(elem);

        // we use 0 as the temporal locality, because it is reasonable to assume that an
        // entry will only be accessed once.
        Opm::prefetch</*temporalLocality=*/0>(elemData_[elemIdx]);
    }

    const Data& get(const Element& elem, unsigned localDofIdx) const
    {
        unsigned elemIdx = elementMapper_.index(elem);
        return elemData_[elemIdx][localDofIdx];
    }

private:
    unsigned computeNumLocalDofs_() const
    {
        unsigned result = 0;

        // loop over the whole grid and sum up the number of local DOFs of all Stencils
        Stencil stencil(gridView_, dofMapper_);
        auto elemIt = gridView_.template begin</*codim=*/0>();
        const auto& elemEndIt = gridView_.template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            stencil.update(*elemIt);
            result += stencil.numDof();
        }

        return result;
    }

    GridView gridView_;
    ElementMapper elementMapper_;
    const DofMapper& dofMapper_;
    std::vector<Data> data_;
    std::vector<Data*> elemData_;
};

} // namespace Opm

#endif
