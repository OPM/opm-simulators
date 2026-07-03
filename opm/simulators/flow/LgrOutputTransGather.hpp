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

/// Connection transmissibilities gathered for parallel LGR INIT output.
///
/// Records are sorted by key (established by gatherLgrOutputTrans below) for
/// binary-search lookup. The I/O rank holds every connection of the global
/// grid, so a node-based container is deliberately avoided.
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
    std::vector<std::pair<std::array<int,3>, double>> sameLevel;
    std::vector<std::pair<std::array<int,4>, double>> crossLevel;
};

namespace detail {

/// Zip flattened (key, value) records into key-sorted records ready for
/// binary-search lookup. keys holds N ints per record, values one double.
template <std::size_t N>
std::vector<std::pair<std::array<int,N>, double>>
sortLgrTransRecords(const std::vector<int>& keys,
                    const std::vector<double>& values)
{
    std::vector<std::pair<std::array<int,N>, double>> records;
    records.reserve(values.size());
    for (std::size_t i = 0; i < values.size(); ++i) {
        std::array<int,N> key;
        for (std::size_t j = 0; j < N; ++j) {
            key[j] = keys[N*i + j];
        }
        records.emplace_back(key, values[i]);
    }
    std::sort(records.begin(), records.end());
    // Keys are unique by construction: every connection is recorded by exactly
    // one rank (see gatherLgrOutputTrans below).
    assert(std::adjacent_find(records.begin(), records.end(),
                              [](const auto& a, const auto& b)
                              { return a.first == b.first; }) == records.end());
    return records;
}

} // namespace detail

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
/// \return the complete, key-sorted records on rank 0; empty records on all other ranks.
template <class GridView, class TransFn>
GatheredLgrOutputTrans
gatherLgrOutputTrans(const Dune::CpGrid& grid,
                     const GridView& gridView,
                     TransFn&& transFn)
{
    std::vector<int> sameKeys;      // flattened: level, minCart, maxCart per record
    std::vector<double> sameValues;
    std::vector<int> crossKeys;     // flattened: smallLevel, smallCart, largeLevel, largeCart
    std::vector<double> crossValues;

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

                crossKeys.push_back(levelIn);
                crossKeys.push_back(cartIn);
                crossKeys.push_back(levelOut);
                crossKeys.push_back(levelCartMapp.cartesianIndex(
                    outside.getLevelElem().index(), levelOut));
                crossValues.push_back(transFn(idxIn, elemMapper.index(outside)));
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

            sameKeys.push_back(levelIn);
            sameKeys.push_back(cartIn);
            sameKeys.push_back(cartOut);
            sameValues.push_back(transFn(idxIn, elemMapper.index(outside)));
        }
    }

    // Gather every rank's records on the I/O rank; gatherv is collective and
    // returns empty vectors on the other ranks. MPI's int counts/displacements
    // bound the total record count across ALL ranks at ~2^31 ints -- a shared
    // MPI-wide ceiling (~500M same-level records), not a per-rank one.
    const auto& comm = grid.comm();
    const auto allSameKeys    = gatherv(sameKeys, comm, 0).first;
    const auto allSameValues  = gatherv(sameValues, comm, 0).first;
    const auto allCrossKeys   = gatherv(crossKeys, comm, 0).first;
    const auto allCrossValues = gatherv(crossValues, comm, 0).first;

    GatheredLgrOutputTrans gathered;
    gathered.sameLevel  = detail::sortLgrTransRecords<3>(allSameKeys, allSameValues);
    gathered.crossLevel = detail::sortLgrTransRecords<4>(allCrossKeys, allCrossValues);
    return gathered;
}

} // namespace Opm

#endif // OPM_LGR_OUTPUT_TRANS_GATHER_HPP
