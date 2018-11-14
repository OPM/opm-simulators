/*
  Copyright 2015 Andreas Lauser

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
#ifndef OPM_CREATE_GLOBAL_CELL_ARRAY_HPP
#define OPM_CREATE_GLOBAL_CELL_ARRAY_HPP

#include <opm/grid/GridHelpers.hpp>

namespace Opm
{
/*!
 * \brief Create a mapping from a global cell index of a grid to the logically
 *        Cartesian index of the ECL deck.
 */
template <class Grid>
void createGlobalCellArray(const Grid &grid, std::vector<int>& dest)
{
    int numCells = Opm::UgGridHelpers::numCells(grid);
    dest.resize(numCells);
    const auto& globalCell = Opm::UgGridHelpers::globalCell(grid);
    std::vector<int> compressedToCartesianIdx(numCells);
    for (int cellIdx = 0; cellIdx < numCells; ++cellIdx) {
        if (globalCell) {
            dest[cellIdx] = globalCell[cellIdx];
        }
        else {
            dest[cellIdx] = cellIdx;
        }
    }

}
}

#endif
