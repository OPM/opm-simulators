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

#include <any>
#include <vector>

namespace Opm {
namespace detail {


    std::vector<int> buildAllCells(const int nc);



    template <class PU>
    std::vector<bool>
    activePhases(const PU& pu)
    {
        const int maxnp = BlackoilPhases::MaxNumPhases;
        std::vector<bool> active(maxnp, false);

        for (int p = 0; p < pu.MaxNumPhases; ++p) {
            active[ p ] = pu.phase_used[ p ] != 0;
        }

        return active;
    }



    template <class PU>
    std::vector<int>
    active2Canonical(const PU& pu)
    {
        const int maxnp = BlackoilPhases::MaxNumPhases;
        std::vector<int> act2can(maxnp, -1);

        for (int phase = 0; phase < maxnp; ++phase) {
            if (pu.phase_used[ phase ]) {
                act2can[ pu.phase_pos[ phase ] ] = phase;
            }
        }

        return act2can;
    }



    double getGravity(const double* g, const int dim);

        /// \brief Get the number of local interior cells in a grid.
        /// \tparam The type of the DUNE grid.
        /// \param grid The grid which cells we count
        /// \return The number of interior cell in the partition of the
        /// grid stored on this process.
        template<class Grid>
        std::size_t countLocalInteriorCells(const Grid& grid)
        {
            if ( grid.comm().size() == 1)
            {
                return grid.size(0);
            }
            std::size_t count = 0;
            const auto& gridView = grid.leafGridView();
            for(auto cell = gridView.template begin<0, Dune::Interior_Partition>(),
                    endCell = gridView.template end<0, Dune::Interior_Partition>();
                cell != endCell; ++cell)
            {
                    ++count;
            }
            return count;
        }

        /// \brief Get the number of cells of a global grid.
        ///
        /// In a parallel run this is the number of cells that a grid would
        /// have if the whole grid was stored on one process only.
        /// \tparam The type of the DUNE grid.
        /// \param grid The grid which cells we count
        /// \return The global number of cells.
        template<class Grid>
        std::size_t countGlobalCells(const Grid& grid)
        {
            if ( grid.comm().size() == 1)
            {
                return grid.size(0);
            }
            std::size_t count = countLocalInteriorCells(grid);
            return grid.comm().sum(count);
        }

    } // namespace detail
} // namespace Opm

#endif // OPM_BLACKOILDETAILS_HEADER_INCLUDED
