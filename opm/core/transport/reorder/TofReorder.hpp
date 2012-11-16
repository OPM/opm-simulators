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

#ifndef OPM_TOFREORDER_HEADER_INCLUDED
#define OPM_TOFREORDER_HEADER_INCLUDED

#include <opm/core/transport/reorder/ReorderSolverInterface.hpp>
#include <vector>
#include <map>
#include <ostream>
struct UnstructuredGrid;

namespace Opm
{

    class IncompPropertiesInterface;

    /// Implements a first-order finite volume solver for
    /// (single-phase) time-of-flight using reordering.
    /// The equation solved is:
    ///     v \cdot \grad\tau = \phi
    /// where v is the fluid velocity, \tau is time-of-flight and
    /// \phi is the porosity. This is a boundary value problem, where
    /// \tau is specified to be zero on all inflow boundaries.
    class TofReorder : public ReorderSolverInterface
    {
    public:
        /// Construct solver.
        /// \param[in] grid      A 2d or 3d grid.
        TofReorder(const UnstructuredGrid& grid);

        /// Solve for time-of-flight.
        /// \param[in]  darcyflux         Array of signed face fluxes.
        /// \param[in]  porevolume        Array of pore volumes.
        /// \param[in]  source            Source term. Sign convention is:
        ///                                 (+) inflow flux,
        ///                                 (-) outflow flux.
        /// \param[out] tof               Array of time-of-flight values.
        void solveTof(const double* darcyflux,
                      const double* porevolume,
                      const double* source,
                      std::vector<double>& tof);

    private:
        virtual void solveSingleCell(const int cell);
        virtual void solveMultiCell(const int num_cells, const int* cells);

    private:
        const UnstructuredGrid& grid_;
        const double* darcyflux_;   // one flux per grid face
        const double* porevolume_;  // one volume per cell
        const double* source_;      // one volumetric source term per cell
        double* tof_;
    };

} // namespace Opm

#endif // OPM_TRANSPORTMODELTRACERTOF_HEADER_INCLUDED
