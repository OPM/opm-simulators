/*
  Copyright 2013, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014, 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2014, 2015 Statoil ASA.
  Copyright 2015 NTNU
  Copyright 2015, 2016, 2017 IRIS AS

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

#ifndef AQUIFER_GRID_UTILS_HEADER_INCLUDED
#define AQUIFER_GRID_UTILS_HEADER_INCLUDED

#include <opm/grid/CpGrid.hpp>

#include <opm/input/eclipse/EclipseState/Aquifer/Aquancon.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>

#include <dune/grid/common/rangegenerators.hh>

#include <algorithm>
#include <map>
#include <set>
#include <type_traits>
#include <unordered_map>
#include <vector>

namespace Opm {

template<class Grid>
struct IsNumericalAquiferCell {
    explicit IsNumericalAquiferCell(const Grid&)
    {}

    template<class T>
    bool operator()(const T&) const { return false; }
};

template<>
struct IsNumericalAquiferCell<Dune::CpGrid> {
    explicit IsNumericalAquiferCell(const Dune::CpGrid& grid)
        : grid_(grid)
    {}

    template<class T>
    bool operator()(const T& elem) const
    {
        const auto& aquiferCells = grid_.sortedNumAquiferCells();
        if (aquiferCells.empty())
        {
          return false;
        }
        auto candidate = std::lower_bound(aquiferCells.begin(),
                                          aquiferCells.end(), elem.index());
        return candidate != aquiferCells.end() && *candidate == elem.index();
    }

private:
    const Dune::CpGrid& grid_;
};

/// Aquifer connections as the leaf grid has them.
struct LeafAquiferConnections
{
    std::vector<Aquancon::AquancCell> connections;
    std::vector<int> leafCell;           //!< -1: found by the connection's global index
    std::vector<std::size_t> origin;     //!< index of the deck connection
};

/// A connection on a refined host, which the deck can only name by the host, is shared
/// among the host's children with a boundary face on that side, by face area.
template<class Simulator>
LeafAquiferConnections
aquiferConnectionsOnLeaf(const std::vector<Aquancon::AquancCell>& connections,
                         const Simulator& simulator)
{
    auto out = LeafAquiferConnections{};
    const auto& vanguard = simulator.vanguard();
    const auto& gridView = vanguard.gridView();
    using Element = typename std::remove_cvref_t<decltype(gridView)>::template Codim<0>::Entity;

    auto children = std::unordered_map<std::size_t, std::vector<std::pair<int, Element>>>{};
    if constexpr (std::is_same_v<std::remove_cvref_t<decltype(vanguard.grid())>, Dune::CpGrid>) {
        if (vanguard.grid().maxLevel() > 0) {
            const auto& mapper = simulator.model().dofMapper();
            for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
                if (elem.level() > 0) {
                    const auto idx = static_cast<int>(mapper.index(elem));
                    children[vanguard.cartesianIndex(idx)].emplace_back(idx, elem);
                }
            }
        }
    }

    for (std::size_t c = 0; c < connections.size(); ++c) {
        const auto& conn = connections[c];
        const auto host = children.find(conn.global_index);
        auto area = std::map<int, double>{};
        auto total = 0.0;
        if ((host != children.end()) && (vanguard.compressedIndex(conn.global_index) < 0)) {
            for (const auto& [idx, elem] : host->second) {
                for (const auto& is : intersections(gridView, elem)) {
                    if (is.boundary() &&
                        (FaceDir::FromIntersectionIndex(is.indexInInside()) == conn.face_dir)) {
                        area[idx] += is.geometry().volume();
                        total += is.geometry().volume();
                    }
                }
            }
        }
        if (total > 0.0) {
            for (const auto& [idx, a] : area) {
                auto child = conn;
                child.influx_coeff *= a / total;
                child.effective_facearea *= a / total;
                out.connections.push_back(child);
                out.leafCell.push_back(idx);
                out.origin.push_back(c);
            }
        }
        else {
            out.connections.push_back(conn);
            out.leafCell.push_back(-1);
            out.origin.push_back(c);
        }
    }
    return out;
}

} // namespace Opm

#endif
