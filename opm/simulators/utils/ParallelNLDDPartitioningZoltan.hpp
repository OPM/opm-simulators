/*
  Copyright 2023 Equinor ASA

  This file is part of the Open Porous Media Project (OPM).

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

#ifndef OPM_UTILITY_PARALLEL_NLDD_PARTITIONING_ZOLTAN_HPP
#define OPM_UTILITY_PARALLEL_NLDD_PARTITIONING_ZOLTAN_HPP

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <cstddef>
#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace Opm {

    /// Partition connectivity graph into non-overlapping domains using the
    /// Zoltan graph partitioning software package.  Primarily intended for
    /// use in NLDD-based non-linear solves.
    class ParallelNLDDPartitioningZoltan
    {
    public:
        /// Callback type for mapping a vertex/cell ID to a globally unique
        /// ID.
        using GlobalCellID = std::function<int(int)>;

        /// Collection of parameters to control the partitioning procedure.
        using ZoltanParamMap = std::map<std::string, std::string>;

        /// Constructor.
        ///
        /// \param[in] comm MPI communication object.  Needed by Zoltan.
        ///
        /// \param[in] numElements Number of potential vertices in
        ///   connectivity graph.  Typically the total number of cells on
        ///   the current rank, i.e., both owned cells and overlap cells.
        ///
        /// \param[in] globalCell Callback for mapping (local) vertex IDs to
        ///   globally unique vertex IDs.
        explicit ParallelNLDDPartitioningZoltan(const Parallel::Communication comm,
                                                const std::size_t             numElements,
                                                const GlobalCellID&           globalCell)
            : comm_        { comm }
            , numElements_ { numElements }
            , globalCell_  { globalCell }
        {}

        /// Insert directed graph edge between two vertices.
        ///
        /// \param[in] c1 Source vertex.
        ///
        /// \param]in] c2 Sink/destination vertex.
        void registerConnection(std::size_t c1, std::size_t c2)
        {
            this->conns_.emplace_back(c1, c2);
        }

        /// Force collection of cells to be in same result domain.
        ///
        /// Mostly as a means to ensuring wells do not intersect multiple
        /// domains/blocks.
        ///
        /// \param[in] cells Cell collection.  Typically those cells which
        ///   are intersected by a single well.
        void forceSameDomain(std::vector<int>&& cells)
        {
            this->sameDomain_.emplace_back(std::move(cells));
        }

        /// Partition connectivity graph using Zoltan graph partitioning
        /// package.
        ///
        /// Honours any prescribed requirement that certain cells be placed
        /// in a single domain/block.
        ///
        /// \param[in] params Parameters for Zoltan.  Override default
        ///   settings.
        ///
        /// \return Partition vector of size \p numElements.  Reachable
        ///   vertices/cells are partitioned into N blocks numbered 0..N-1.
        ///   Unreachable vertices get a block ID of -1.
        std::vector<int>
        partitionElements(const ZoltanParamMap& params) const;

    private:
        /// Connection/graph edge.
        using Connection = std::pair<std::size_t, std::size_t>;

        /// MPI communication object.  Needed by Zoltan.
        Parallel::Communication comm_{};

        /// Maximum number of vertices in connectivity graph.  Needed to
        /// size partition vector return value.
        std::size_t numElements_{};

        /// Callback for mapping (local) vertex/cell IDs to globally unique
        /// vertex/cell IDs.
        GlobalCellID globalCell_{};

        /// Connectivity graph edges.
        std::vector<Connection> conns_{};

        /// Collections of vertices/cells which must be coalesced to the
        /// same domain/block.  All vertices within a single collection will
        /// be placed on the same domain, but cells from different
        /// collections may be placed on different domains.
        std::vector<std::vector<int>> sameDomain_{};
    };

} // namespace Opm

#endif  // OPM_UTILITY_PARALLEL_NLDD_PARTITIONING_ZOLTAN_HPP
