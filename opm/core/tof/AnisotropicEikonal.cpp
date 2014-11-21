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
#include <set>

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
	considered_.reserve(100);
    }

    /// Solve the eikonal equation.
    /// \param[in]  metric            Array of metric tensors, M, for each cell.
    /// \param[in]  startcells        Array of cells where u = 0 at the centroid.
    /// \param[out] solution          Array of solution to the eikonal equation.
    void AnisotropicEikonal2d::solve(const double* metric,
				     const std::vector<int>& startcells,
				     std::vector<double>& solution)
    {
	// The algorithm used is described in J.A. Sethian and A. Vladimirsky,
	// "Ordered Upwind Methods for Static Hamilton-Jacobi Equations".
	// Notation in comments is as used in that paper: U is the solution,
	// and q is the boundary condition. One difference is that we talk about
	// grid cells instead of mesh points.
	//
	// Algorithm summary:
	// 1. Put all cells in Far. U_i = \inf.
	// 2. Move the startcells to Accepted. U_i = q(x_i)
	// 3. Move cells adjacent to startcells to Considered, evaluate
	//    U_i = min_{(x_j,x_k) \in NF(x_i)} G_{j,k}
	// 4. Find the Considered cell with the smallest value: r.
	// 5. Move cell r to Accepted. Update AcceptedFront.
	// 6. Move cells adjacent to r from Far to Considered.
	// 7. Recompute the value for all Considered cells within
	//    distance h * F_2/F1 from x_r. Use min of previous and new.
	// 8. If Considered is not empty, go to step 4.

	// 1. Put all cells in Far. U_i = \inf.
	const int num_cells = grid_.number_of_cells;
	const double inf = 1e100;
	solution.clear();
	solution.resize(num_cells, inf);
	considered_.clear();
	is_considered_.clear();
	is_considered_.resize(num_cells, false);

	// 2. Move the startcells to Accepted. U_i = q(x_i)
	std::vector<char> accepted(num_cells, false);
	const int num_startcells = startcells.size();
	for (int ii = 0; ii < num_startcells; ++ii) {
	    accepted[startcells[ii]] = true;
	    solution[startcells[ii]] = 0.0;
	}
	std::set<int> accepted_front(startcells.begin(), startcells.end());

	// 3. Move cells adjacent to startcells to Considered, evaluate
	//    U_i = min_{(x_j,x_k) \in NF(x_i)} G_{j,k}
	for (int ii = 0; ii < num_startcells; ++ii) {
	    const int scell = startcells[ii];
	    const int num_nb = cell_neighbours_[scell].size();
	    for (int nb = 0; nb < num_nb; ++nb) {
		const int nb_cell = cell_neighbours_[scell][nb];
		if (!is_considered_[nb_cell]) {
		    const double value = computeValue(nb_cell);
		    pushConsidered(std::make_pair(value, nb_cell));
		    is_considered_[nb_cell] = true;
		}
	    }
	}

	// 4. Find the Considered cell with the smallest value: r.

	// 5. Move cell r to Accepted. Update AcceptedFront.

	// 6. Move cells adjacent to r from Far to Considered.

	// 7. Recompute the value for all Considered cells within
	//    distance h * F_2/F1 from x_r. Use min of previous and new.

	// 8. If Considered is not empty, go to step 4.

    }





    double AnisotropicEikonal2d::computeValue(const int cell) const
    {
	const auto& nbs = cell_neighbours_[cell];
	const int num_nbs = nbs.size();
	double val = 1e100;
	for (int ii = 0; ii < num_nbs; ++ii) {
	    const int n[2] = { nbs[ii], nbs[(ii+1) % num_nbs] };
	    // if ... accepted front
	}
	return val;
    }





    const AnisotropicEikonal2d::ValueAndCell& AnisotropicEikonal2d::topConsidered() const
    {
	return considered_.front();
    }





    void AnisotropicEikonal2d::pushConsidered(const ValueAndCell& vc)
    {
	considered_.push_back(vc);
	std::push_heap(considered_.begin(), considered_.end());
    }





    void AnisotropicEikonal2d::popConsidered()
    {
	std::pop_heap(considered_.begin(), considered_.end());
	considered_.pop_back();
    }





} // namespace Opm
