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

#include "countGlobalCells.hpp"

#include <cassert>

namespace Opm {
namespace detail {

std::vector<int> buildAllCells(const int nc) {
    std::vector<int> all_cells(nc);
    std::iota(all_cells.begin(), all_cells.end(), 0);
    return all_cells;
}

double getGravity(const double* g, const int dim) {
    double grav = 0.0;
    if (g) {
        // Guard against gravity in anything but last dimension.
        for (int dd = 0; dd < dim - 1; ++dd) {
            assert(g[dd] == 0.0);
        }
        grav = g[dim - 1];
    }
    return grav;
}

} // namespace detail
} // namespace Opm
