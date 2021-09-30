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
#ifndef EWOMS_COLLECT_TO_IO_RANK_HH
#define EWOMS_COLLECT_TO_IO_RANK_HH

#include <opm/output/data/Aquifer.hpp>
#include <opm/output/data/Cells.hpp>
#include <opm/output/data/Solution.hpp>
#include <opm/output/data/Wells.hpp>
#include <opm/output/data/Groups.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/WellTestState.hpp>

#include <opm/grid/common/p2pcommunicator.hh>

#include <map>
#include <utility>
#include <vector>

namespace Dune {
template<class Grid> class CartesianIndexMapper;
}

namespace Opm {

template <class Grid, class EquilGrid, class GridView>
class CollectDataToIORank
{
public:
    using CollectiveCommunication = typename Grid::CollectiveCommunication;
    using P2PCommunicatorType = Dune::Point2PointCommunicator<Dune::SimpleMessageBuffer>;
    using IndexMapType = std::vector<int>;
    using IndexMapStorageType = std::vector<IndexMapType>;

    static const int dimension = Grid::dimension;

    enum { ioRank = 0 };

    static const bool needsReordering =
        !std::is_same<Grid, EquilGrid>::value;

    CollectDataToIORank(const Grid& grid,
                        const EquilGrid* equilGrid,
                        const GridView& gridView,
                        const Dune::CartesianIndexMapper<Grid>& cartMapper,
                        const Dune::CartesianIndexMapper<EquilGrid>* equilCartMapper);

    // gather solution to rank 0 for EclipseWriter
    void collect(const data::Solution& localCellData,
                 const std::map<std::pair<std::string, int>, double>& localBlockData,
                 const std::map<std::size_t, double>& localWBPData,
                 const data::Wells& localWellData,
                 const data::GroupAndNetworkValues& localGroupAndNetworkData,
                 const data::Aquifers& localAquiferData,
                 const WellTestState& localWellTestState);

    const std::map<std::size_t, double>& globalWBPData() const
    { return this->globalWBPData_; }

    const std::map<std::pair<std::string, int>, double>& globalBlockData() const
    { return globalBlockData_; }

    const data::Solution& globalCellData() const
    { return globalCellData_; }

    const data::Wells& globalWellData() const
    { return globalWellData_; }

    const WellTestState& globalWellTestState() const
    { return this->globalWellTestState_; }

    const data::GroupAndNetworkValues& globalGroupAndNetworkData() const
    { return globalGroupAndNetworkData_; }

    const data::Aquifers& globalAquiferData() const
    { return globalAquiferData_; }

    bool isIORank() const
    { return toIORankComm_.rank() == ioRank; }

    bool isParallel() const
    { return toIORankComm_.size() > 1; }

    int localIdxToGlobalIdx(unsigned localIdx) const;

    size_t numCells () const
    { return globalCartesianIndex_.size(); }

    const std::vector<int>& globalRanks() const
    { return globalRanks_; }

    bool isCartIdxOnThisRank(int cartIdx) const;

protected:
    P2PCommunicatorType toIORankComm_;
    IndexMapType globalCartesianIndex_;
    IndexMapType localIndexMap_;
    IndexMapStorageType indexMaps_;
    std::vector<int> globalRanks_;
    data::Solution globalCellData_;
    std::map<std::pair<std::string, int>, double> globalBlockData_;
    std::map<std::size_t, double> globalWBPData_;
    data::Wells globalWellData_;
    data::GroupAndNetworkValues globalGroupAndNetworkData_;
    data::Aquifers globalAquiferData_;
    WellTestState globalWellTestState_;
    std::vector<int> localIdxToGlobalIdx_;
    /// \brief sorted list of cartesian indices present-
    ///
    /// non-empty only when running in parallel
    std::vector<int> sortedCartesianIdx_;
};

} // end namespace Opm

#endif
