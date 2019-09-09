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
 * \copydoc Opm::Linear::VertexBorderListFromGrid
 */
#ifndef EWOMS_VERTEX_BORDER_LIST_FROM_GRID_HH
#define EWOMS_VERTEX_BORDER_LIST_FROM_GRID_HH

#include "overlaptypes.hh"
#include "blacklist.hh"

#include <opm/material/common/Unused.hpp>

#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/operators.hh>
#include <dune/common/version.hh>

#include <algorithm>

namespace Opm {
namespace Linear {
/*!
 * \brief Uses communication on the grid to find the initial seed list
 *        of indices.
 *
 * \todo implement this class generically. For this, it must be
 *       possible to query the mapper whether it contains entities of
 *       a given codimension without the need to hand it an actual
 *       entity.
 */
template <class GridView, class VertexMapper>
class VertexBorderListFromGrid
    : public Dune::CommDataHandleIF<VertexBorderListFromGrid<GridView, VertexMapper>,
                                    int>
{
    static const int dimWorld = GridView::dimensionworld;

public:
    VertexBorderListFromGrid(const GridView& gridView, const VertexMapper& map)
        : gridView_(gridView), map_(map)
    {
        gridView.communicate(*this,
                             Dune::InteriorBorder_InteriorBorder_Interface,
                             Dune::ForwardCommunication);

        auto vIt = gridView.template begin<dimWorld>();
        const auto& vEndIt = gridView.template end<dimWorld >();
        for (; vIt != vEndIt; ++vIt) {
            if (vIt->partitionType() != Dune::InteriorEntity
                && vIt->partitionType() != Dune::BorderEntity)
            {
                Index vIdx = static_cast<Index>(map_.index(*vIt));
                blackList_.addIndex(vIdx);
            }
        }
    }

    // data handle methods
    bool contains(int dim, int codim) const
    { return dim == codim; }

    bool fixedsize(int dim OPM_UNUSED, int codim OPM_UNUSED) const
    { return true; }

    template <class EntityType>
    size_t size(const EntityType& e OPM_UNUSED) const
    { return 2; }

    template <class MessageBufferImp, class EntityType>
    void gather(MessageBufferImp& buff, const EntityType& e) const
    {
        buff.write(static_cast<int>(gridView_.comm().rank()));
        buff.write(static_cast<int>(map_.index(e)));
    }

    template <class MessageBufferImp, class EntityType>
    void scatter(MessageBufferImp& buff, const EntityType& e, size_t n OPM_UNUSED)
    {
        BorderIndex bIdx;

        bIdx.localIdx = static_cast<Index>(map_.index(e));
        {
            int tmp;
            buff.read(tmp);
            bIdx.peerRank = static_cast<ProcessRank>(tmp);
        }
        {
            int tmp;
            buff.read(tmp);
            bIdx.peerIdx = static_cast<Index>(tmp);
        }
        bIdx.borderDistance = 0;

        borderList_.push_back(bIdx);
    }

    // Access to the border list.
    const BorderList& borderList() const
    { return borderList_; }

    // Access to the black-list indices.
    const BlackList& blackList() const
    { return blackList_; }

private:
    const GridView gridView_;
    const VertexMapper& map_;
    BorderList borderList_;
    BlackList blackList_;
};

} // namespace Linear
} // namespace Opm

#endif
