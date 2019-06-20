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

#include <vector>
#include <utility>

namespace Opm
{
namespace detail
{

    /// \brief Find the rows corresponding to overlap cells
    ///
    /// Loop over grid and store cell ids of row-column pairs
    /// corresponding to overlap cells.
    /// \tparam The type of the DUNE grid.
    /// \param grid The grid where we look for overlap cells.
    /// \param overlapRowAndColumns List where overlap rows and columns are stored.
    template <class Grid>
    void findOverlapRowsAndColumns(const Grid& grid,
                                   std::vector<std::pair<int, std::vector<int>>>& overlapRowAndColumns)
    {
        // only relevant in parallel case.
        if (grid.comm().size() > 1) {
            // Numbering of cells
            auto lid = grid.localIdSet();

            const auto& gridView = grid.leafGridView();
            auto elemIt = gridView.template begin<0>();
            const auto& elemEndIt = gridView.template end<0>();

            // loop over cells in mesh
            for (; elemIt != elemEndIt; ++elemIt) {
                const auto& elem = *elemIt;

                // If cell has partition type not equal to interior save row
                if (elem.partitionType() != Dune::InteriorEntity) {
                    // local id of overlap cell
                    int lcell = lid.id(elem);

                    std::vector<int> columns;
                    // loop over faces of cell
                    auto isend = gridView.iend(elem);
                    for (auto is = gridView.ibegin(elem); is != isend; ++is) {
                        // check if face has neighbor
                        if (is->neighbor()) {
                            // get index of neighbor cell
                            int ncell = lid.id(is->outside());
                            columns.push_back(ncell);
                        }
                    }
                    // add row to list
                    overlapRowAndColumns.push_back(std::pair<int, std::vector<int>>(lcell, columns));
                }
            }
        }
    }

} // namespace detail
} // namespace Opm

#endif // OPM_FINDOVERLAPROWSANDCOLUMNS_HEADER_INCLUDED
