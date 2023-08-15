/*
  Copyright 2018 Andreas Thune

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

#ifndef OPM_FINDOVERLAPROWSANDCOLUMNS_HEADER_INCLUDED
#define OPM_FINDOVERLAPROWSANDCOLUMNS_HEADER_INCLUDED

#include <opm/grid/common/WellConnections.hpp>
#include <opm/grid/common/CartesianIndexMapper.hpp>

#include <cstddef>
#include <utility>
#include <vector>

namespace Dune {
template<class Grid> class CartesianIndexMapper;
}

namespace Opm
{
namespace detail
{
    /// \brief Find cell IDs for wells contained in local grid.
    ///
    /// Cell IDs of wells stored in a graph, so it can be used to create
    /// an adjacency pattern. Only relevant when the UseWellContribusion option is set to true
    /// \tparam The type of the DUNE grid.
    /// \tparam Well vector type
    /// \param grid The grid where we look for overlap cells.
    /// \param wells List of wells contained in grid.
    /// \param useWellConn Boolean that is true when UseWellContribusion is true
    /// \param wellGraph Cell IDs of well cells stored in a graph.
    template<class Grid, class CartMapper, class W>
    void setWellConnections(const Grid& grid, const CartMapper& cartMapper, const W& wells, bool useWellConn, std::vector<std::set<int>>& wellGraph, int numJacobiBlocks)
    {
        if ( grid.comm().size() > 1 || numJacobiBlocks > 1)
        {
            const int numCells = cartMapper.compressedSize(); // grid.numCells()
            wellGraph.resize(numCells);

            if (useWellConn) {
                const auto& cpgdim = cartMapper.cartesianDimensions();

                std::vector<int> cart(cpgdim[0]*cpgdim[1]*cpgdim[2], -1);

                for( int i=0; i < numCells; ++i )
                    cart[ cartMapper.cartesianIndex( i ) ] = i;

                Dune::cpgrid::WellConnections well_indices;
                well_indices.init(wells, cpgdim, cart);

                for (auto& well : well_indices)
                {
                    for (auto perf = well.begin(); perf != well.end(); ++perf)
                    {
                        auto perf2 = perf;
                        for (++perf2; perf2 != well.end(); ++perf2)
                        {
                            wellGraph[*perf].insert(*perf2);
                            wellGraph[*perf2].insert(*perf);
                        }
                    }
                }
            }
        }
    }

    /// \brief Find the rows corresponding to overlap cells
    ///
    /// Loop over grid and store cell ids of rows
    /// corresponding to overlap cells.
    /// \tparam The type of the DUNE grid.
    /// \param grid The grid where we look for overlap cells.
    /// \param overlapRows List where overlap rows are stored.
    /// \param interiorRows List where overlap rows are stored.
    template<class Grid, class Mapper>
    void findOverlapAndInterior(const Grid& grid, const Mapper& mapper, std::vector<int>& overlapRows,
                                std::vector<int>& interiorRows)
    {
        //only relevant in parallel case.
        if ( grid.comm().size() > 1)
        {
            //Numbering of cells
            const auto& gridView = grid.leafGridView();
            //loop over cells in mesh
            for (const auto& elem : elements(gridView))
            {
                int lcell = mapper.index(elem);

                if (elem.partitionType() != Dune::InteriorEntity)
                {
                    //add row to list
                    overlapRows.push_back(lcell);
                } else {
                    interiorRows.push_back(lcell);
                }
            }
        }
    }

    /// \brief If ownerFirst=true, returns the number of interior cells in grid, else just numCells().
    ///
    /// If cells in grid is ordered so that interior/owner cells come before overlap/copy cells, the method
    /// returns the number of interior cells numInterior. In the linear solver only the first numInterior rows of
    /// the matrix are needed.
    template <class Grid>
    std::size_t numMatrixRowsToUseInSolver(const Grid& grid, bool ownerFirst)
    {
        std::size_t numInterior = 0;
        if (!ownerFirst || grid.comm().size()==1)
            return grid.leafGridView().size(0);
        const auto& gridView = grid.leafGridView();

        // loop over cells in mesh
        const auto& range = elements(gridView, Dune::Partitions::interior);
        numInterior = std::distance(range.begin(), range.end());

        return numInterior;
    }
} // namespace detail
} // namespace Opm

#endif // OPM_FINDOVERLAPROWSANDCOLUMNS_HEADER_INCLUDED
