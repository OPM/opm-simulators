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

#include <opm/core/transport/reorder/TransportModelTracerTof.hpp>
#include <opm/core/grid.h>
#include <opm/core/utility/ErrorMacros.hpp>
#include <algorithm>

namespace Opm
{


    /// Construct solver.
    /// \param[in] grid      A 2d or 3d grid.
    TransportModelTracerTof::TransportModelTracerTof(const UnstructuredGrid& grid)
        : grid_(grid)
    {
    }




    /// Solve for time-of-flight at next timestep.
    /// \param[in] darcyflux         Array of signed face fluxes.
    /// \param[in] porevolume        Array of pore volumes.
    /// \param[in, out] tof          Array of time-of-flight values.
    void TransportModelTracerTof::solveTof(const double* darcyflux,
                                           const double* porevolume,
                                           std::vector<double>& tof)
    {
        darcyflux_ = darcyflux;
        porevolume_ = porevolume;
        tof.resize(grid_.number_of_cells);
        std::fill(tof.begin(), tof.end(), 0.0);
        tof_ = &tof[0];
        reorderAndTransport(grid_, darcyflux);
    }




    void TransportModelTracerTof::solveSingleCell(const int cell)
    {
        // Compute flux terms.
        double upwind_term = 0.0;
        double downwind_flux = 0.0;
        for (int i = grid_.cell_facepos[cell]; i < grid_.cell_facepos[cell+1]; ++i) {
            int f = grid_.cell_faces[i];
            double flux;
            int other;
            // Compute cell flux
            if (cell == grid_.face_cells[2*f]) {
                flux  = darcyflux_[f];
                other = grid_.face_cells[2*f+1];
            } else {
                flux  =-darcyflux_[f];
                other = grid_.face_cells[2*f];
            }
            // Add flux to upwind_term or downwind_flux, if interior.
            if (other != -1) {
                if (flux < 0.0) {
                    upwind_term  += flux*tof_[other];
                } else {
                    downwind_flux += flux;
                }
            }
        }

        // Compute tof.
        tof_[cell] = (porevolume_[cell] - upwind_term)/downwind_flux;
    }




    void TransportModelTracerTof::solveMultiCell(const int num_cells, const int* cells)
    {
        std::cout << "Pretending to solve multi-cell dependent equation with " << num_cells << " cells." << std::endl;
        for (int i = 0; i < num_cells; ++i) {
            solveSingleCell(cells[i]);
        }
    }




} // namespace Opm
