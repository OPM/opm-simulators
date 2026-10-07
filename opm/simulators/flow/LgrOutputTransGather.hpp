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
#ifndef OPM_LGR_OUTPUT_TRANS_GATHER_HPP
#define OPM_LGR_OUTPUT_TRANS_GATHER_HPP

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/partitionset.hh>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/common/CommunicationUtils.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <utility>
#include <vector>

namespace Opm {

/// Key-sorted index over gathered LGR connection records.
///
/// The I/O rank holds every connection of the global grid and looks each one up
/// once per output pass. The records are sorted on their key only and searched
/// with std::lower_bound. The key is the pair of global cell ids of a
/// connection (see lgrTransKey); values are transmissibilities.
template <std::size_t N>
class LgrTransIndex
{
public:
    using Key = std::array<int, N>;

    LgrTransIndex() = default;

    /// Build the index from the gathered (key, value) records. Keys are unique by
    /// construction (every connection is recorded by exactly one rank).
    explicit LgrTransIndex(std::vector<std::pair<Key, double>> records)
        : records_(std::move(records))
    {
        std::sort(records_.begin(), records_.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });
        assert(std::adjacent_find(records_.begin(), records_.end(),
                                  [](const auto& a, const auto& b) { return a.first == b.first; })
               == records_.end() && "duplicate LGR transmissibility key");
    }

    /// \return pointer to the value for \p key, or nullptr if absent.
    const double* find(const Key& key) const
    {
        const auto it = std::lower_bound(records_.begin(), records_.end(), key,
                                         [](const auto& record, const Key& k) { return record.first < k; });
        return (it != records_.end() && it->first == key) ? &it->second : nullptr;
    }

private:
    std::vector<std::pair<Key, double>> records_;
};

/// Connection transmissibilities gathered for parallel LGR INIT output, keyed by
/// lgrTransKey. Held on the I/O rank; empty on all other ranks.
using GatheredLgrOutputTrans = LgrTransIndex<2>;

/// Key of the connection between two cells of a CpGrid: their global ids,
/// smaller first. Every rank and the global grid on the I/O rank agree on
/// these ids. A global id is the cell's level index plus the entity counts of
/// the levels before it, so it fits in int as those counts do.
template <class Element>
std::array<int,2> lgrTransKey(const Dune::CpGrid& grid,
                              const Element& cell1,
                              const Element& cell2)
{
    const int id1 = static_cast<int>(grid.globalIdSet().id(cell1));
    const int id2 = static_cast<int>(grid.globalIdSet().id(cell2));
    return {std::min(id1, id2), std::max(id1, id2)};
}

/// Gather the simulator's own (distributed) transmissibilities for parallel LGR INIT output.
///
/// Each rank walks its interior leaf cells and records every connection it owns, keyed by
/// the global ids of its two cells (see lgrTransKey), then the records are gathered on the
/// I/O rank (rank 0), whose output walk over the equil grid looks values up by the same ids.
///
/// Every connection is recorded exactly once, by the rank that owns the cell with the smaller
/// global id; that comparison is rank-independent, so at a rank boundary only one of the two
/// owner ranks records the connection -- which it can, because its partner cell is present in
/// its overlap layer (this requires at least one overlap layer). The ids of a level come
/// after those of the levels before it, so a connection between two levels is recorded from
/// the smaller-level side.
///
/// This reuses the values the simulation itself computed in parallel instead of recomputing a
/// whole-grid transmissibility on the I/O rank.
///
/// \return the complete records indexed on rank 0; empty index on all other ranks.
template <class GridView, class TransFn>
GatheredLgrOutputTrans
gatherLgrOutputTrans(const Dune::CpGrid& grid,
                     const GridView& gridView,
                     TransFn&& transFn)
{
    // Build the final (key, value) records directly -- no separate flat key/value
    // buffers, no zip pass afterwards.
    std::vector<std::pair<std::array<int,2>, double>> records;

    const Dune::MultipleCodimMultipleGeomTypeMapper<GridView>
        elemMapper(gridView, Dune::mcmgElementLayout());

    for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
        // The inside cell is the same for every intersection of this element.
        const auto idxIn = elemMapper.index(elem);
        const int idIn = static_cast<int>(grid.globalIdSet().id(elem));

        for (const auto& is : intersections(gridView, elem)) {
            if (!is.neighbor()) {
                continue;
            }

            const auto outside = is.outside();
            const auto key = lgrTransKey(grid, elem, outside);

            if (key[0] != idIn) {
                continue; // recorded by the owner of the cell with the smaller id
            }

            records.emplace_back(key, transFn(idxIn, elemMapper.index(outside)));
        }
    }

    // Gather the records directly to the I/O rank. gatherv is generic; Dune's
    // MPITraits<std::pair<std::array<int,N>,double>> composes a byte-blob array with
    // MPI_DOUBLE, so the record moves as a single MPI datatype. gatherv is collective
    // and returns an empty vector on the non-root ranks. MPI's int counts/displacements
    // bound the total record count across ALL ranks at ~2^31 records -- a shared
    // MPI-wide ceiling, not a per-rank one.
    return GatheredLgrOutputTrans(gatherv(records, grid.comm(), 0).first);
}

} // namespace Opm

#endif // OPM_LGR_OUTPUT_TRANS_GATHER_HPP
