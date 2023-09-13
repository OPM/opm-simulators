/*
  Copyright 2013, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014, 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2014, 2015 Statoil ASA.
  Copyright 2015 NTNU
  Copyright 2015 IRIS AS

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

#ifndef OPM_COUNTGLOBALCELLS_HEADER_INCLUDED
#define OPM_COUNTGLOBALCELLS_HEADER_INCLUDED

#include <opm/core/props/BlackoilPhases.hpp>

#include <dune/grid/common/gridview.hh>

#include <algorithm>
#include <vector>

namespace Opm {
namespace detail {

        /// \brief Get the number of local interior cells in a grid view.
        ///
        /// \tparam GridView Type of DUNE grid view.
        ///
        /// \param[in] gridView Grid view for which to count the number of
        ///    interior cells.
        ///
        /// \return The number of interior cell in the partition of the grid
        ///    stored on this process.
        template <class GridView>
        std::size_t countLocalInteriorCellsGridView(const GridView& gridView)
        {
            if (gridView.comm().size() == 1) {
                return gridView.size(0);
            }

            return std::distance(gridView.template begin<0, Dune::Interior_Partition>(),
                                 gridView.template end<0, Dune::Interior_Partition>());
        }

        /// \brief Get the number of local interior cells in a grid.
        ///
        /// \tparam Grid Type of DUNE grid.
        ///
        /// \param[in] grid Grid for which to count interior cells.
        ///
        /// \return The number of interior cells in the partition of the
        ///    grid stored on this process.
        template<class Grid>
        std::size_t countLocalInteriorCells(const Grid& grid)
        {
            return countLocalInteriorCellsGridView(grid.leafGridView());
        }

        /// \brief Get the number of cells of a global grid.
        ///
        /// In a parallel run this is the number of cells that a grid would
        /// have if the whole grid was stored on one process only.
        ///
        /// \tparam Grid Type of DUNE grid.
        ///
        /// \param[in] grid Grid for which to count global cells.
        ///
        /// \return The global number of cells.
        template<class Grid>
        std::size_t countGlobalCells(const Grid& grid)
        {
            if ( grid.comm().size() == 1)
            {
                return grid.size(0);
            }
            const std::size_t count = countLocalInteriorCellsGridView(grid.leafGridView());
            return grid.comm().sum(count);
        }

    } // namespace detail
} // namespace Opm

#endif // OPM_BLACKOILDETAILS_HEADER_INCLUDED
