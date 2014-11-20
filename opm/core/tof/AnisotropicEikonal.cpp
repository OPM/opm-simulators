/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.

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

#include <opm/core/tof/AnisotropicEikonal.hpp>
#include <opm/core/grid/GridUtilities.hpp>
#include <opm/core/grid.h>

namespace Opm
{

    /// Construct solver.
    /// \param[in] grid      A 2d grid.
    AnisotropicEikonal2d::AnisotropicEikonal2d(const UnstructuredGrid& grid)
	: grid_(grid)
    {
	if (grid.dimensions != 2) {
	    OPM_THROW(std::logic_error, "Grid for AnisotropicEikonal2d must be 2d.");
	}
	cell_neighbours_ = vertexNeighbours(grid);
	orderCounterClockwise(grid, cell_neighbours_);
    }

    /// Solve the eikonal equation.
    /// \param[in]  metric            Array of metric tensors, M, for each cell.
    /// \param[in]  startcells        Array of cells where u = 0 at the centroid.
    /// \param[out] solution          Array of solution to the eikonal equation.
    void AnisotropicEikonal2d::solve(const double* metric,
				     const std::vector<int>& startcells,
				     std::vector<double>& solution) const
    {
	// The algorithm used is described in J.A. Sethian and A. Vladimirsky,
	// "Ordered Upwind Methods for Static Hamilton-Jacobi Equations".

	// 1. 
    }

} // namespace Opm
