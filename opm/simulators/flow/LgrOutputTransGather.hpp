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

#include <opm/grid/common/CommunicationUtils.hpp>
#include <opm/grid/cpgrid/LevelCartesianIndexMapper.hpp>

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
/// with std::lower_bound. The key is the level-Cartesian index tuple (see
/// GatheredLgrOutputTrans); values are transmissibilities.
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

/// Connection transmissibilities gathered for parallel LGR INIT output.
///
/// Held on the I/O rank as key-indexed tables (see LgrTransIndex); empty on all
/// other ranks.
///
/// sameLevel:  key (level, min level-Cartesian index, max level-Cartesian index)
///             -- every connection between two cells of the same level grid
///                (TRANX/TRANY/TRANZ and same-level NNCs).
/// crossLevel: key (smaller level, its level-Cartesian index,
///                  larger level, its level-Cartesian index)
///             -- every connection between cells of different level grids
///                (global<->LGR and LGR<->LGR NNCs; TRANGL/TRANLL).
struct GatheredLgrOutputTrans
{
    LgrTransIndex<3> sameLevel;
    LgrTransIndex<4> crossLevel;
};

/// Gather the simulator's own (distributed) transmissibilities for parallel LGR INIT output.
///
/// Each rank walks its interior leaf cells and records every connection it owns, keyed by
/// level-Cartesian indices (see GatheredLgrOutputTrans), then the records are gathered on the
/// I/O rank (rank 0). The level-Cartesian key is geometrically canonical -- defined by the LGR
/// specification, identical on the distributed grid and the undistributed (equil) copy -- so
/// the I/O rank's output walk over the equil grid can look values up directly.
///
/// Every connection is recorded exactly once. A same-level connection is recorded by the rank
/// that owns the cell with the smaller level-Cartesian index; that comparison is
/// rank-independent, so at a rank boundary only one of the two owner ranks records the
/// connection -- which it can, because its partner cell is present in its overlap layer (this
/// requires at least one overlap layer; the caller guards --num-overlap). A cross-level
/// connection is recorded from the smaller-level side only.
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
    std::vector<std::pair<std::array<int,3>, double>> same;   // level, minCart, maxCart
    std::vector<std::pair<std::array<int,4>, double>> cross;  // smallLevel, smallCart, largeLevel, largeCart

    const LevelCartesianIndexMapper<Dune::CpGrid> levelCartMapp(grid);
    const Dune::MultipleCodimMultipleGeomTypeMapper<GridView>
        elemMapper(gridView, Dune::mcmgElementLayout());

    for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
        // The inside cell is the same for every intersection of this element.
        const int levelIn = elem.level();
        const int cartIn = levelCartMapp.cartesianIndex(elem.getLevelElem().index(), levelIn);
        const auto idxIn = elemMapper.index(elem);

        for (const auto& is : intersections(gridView, elem)) {
            if (!is.neighbor()) {
                continue;
            }

            const auto outside = is.outside();
            const int levelOut = outside.level();

            if (levelIn != levelOut) {
                if (levelIn > levelOut) {
                    continue; // recorded exactly once, from the smaller-level side
                }

                cross.emplace_back(
                    std::array<int,4>{levelIn, cartIn, levelOut,
                                      levelCartMapp.cartesianIndex(outside.getLevelElem().index(), levelOut)},
                    transFn(idxIn, elemMapper.index(outside)));
                continue;
            }

            const int cartOut = levelCartMapp.cartesianIndex(
                outside.getLevelElem().index(), levelIn);

            if (cartIn > cartOut) {
                // Record each connection once, in canonical direction. The comparison
                // is rank-independent, so at a rank boundary exactly one of the two
                // owner ranks records the connection.
                continue;
            }

            same.emplace_back(std::array<int,3>{levelIn, cartIn, cartOut},
                              transFn(idxIn, elemMapper.index(outside)));
        }
    }

    // Gather each record vector directly to the I/O rank -- one collective per kind.
    // gatherv is generic; Dune's
    // MPITraits<std::pair<std::array<int,N>,double>> composes a byte-blob array with
    // MPI_DOUBLE, so the record moves as a single MPI datatype. gatherv is collective
    // and returns empty vectors on the non-root ranks. MPI's int counts/displacements
    // bound the total record count across ALL ranks at ~2^31 records -- a shared
    // MPI-wide ceiling, not a per-rank one.
    const auto& comm = grid.comm();
    auto allSame  = gatherv(same,  comm, 0).first;
    auto allCross = gatherv(cross, comm, 0).first;

    // Index the gathered records on the I/O rank; empty on the others.
    GatheredLgrOutputTrans gathered;
    gathered.sameLevel  = LgrTransIndex<3>(std::move(allSame));
    gathered.crossLevel = LgrTransIndex<4>(std::move(allCross));
    return gathered;
}

} // namespace Opm

#endif // OPM_LGR_OUTPUT_TRANS_GATHER_HPP
