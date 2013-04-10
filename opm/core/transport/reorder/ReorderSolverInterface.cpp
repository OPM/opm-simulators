/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#include "config.h"
#include <opm/core/transport/reorder/ReorderSolverInterface.hpp>
#include <opm/core/transport/reorder/reordersequence.h>
#include <opm/core/grid.h>
#include <opm/core/utility/StopWatch.hpp>

#include <vector>
#include <cassert>


void Opm::ReorderSolverInterface::reorderAndTransport(const UnstructuredGrid& grid, const double* darcyflux)
{
    // Compute reordered sequence of single-cell problems
    sequence_.resize(grid.number_of_cells);
    components_.resize(grid.number_of_cells + 1);
    int ncomponents;
    time::StopWatch clock;
    clock.start();
    compute_sequence(&grid, darcyflux, &sequence_[0], &components_[0], &ncomponents);
    clock.stop();
    std::cout << "Topological sort took: " << clock.secsSinceStart() << " seconds." << std::endl;

    // Make vector's size match actual used data.
    components_.resize(ncomponents + 1);

    // Invoke appropriate solve method for each interdependent component.
    for (int comp = 0; comp < ncomponents; ++comp) {
#if 0
#ifdef MATLAB_MEX_FILE
	// \TODO replace this with general signal handling code, check if it costs performance.
        if (interrupt_signal) {
            mexPrintf("Reorder loop interrupted by user: %d of %d "
                      "cells finished.\n", i, grid.number_of_cells);
            break;
        }
#endif
#endif
	const int comp_size = components_[comp + 1] - components_[comp];
	if (comp_size == 1) {
	    solveSingleCell(sequence_[components_[comp]]);
	} else {
	    solveMultiCell(comp_size, &sequence_[components_[comp]]);
	}
    }
}


const std::vector<int>& Opm::ReorderSolverInterface::sequence() const
{
    return sequence_;
}


const std::vector<int>& Opm::ReorderSolverInterface::components() const
{
    return components_;
}
