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


/// \page tutorial2 Flow Solver for a single phase 
/// \details The flow equations consist of the mass conservation equation
/// \f[\nabla\cdot u=q\f] and the Darcy law \f[u=-\frac{1}{\mu}K\nabla p.\f] Here,
/// \f$u\f$ denotes the velocity and \f$p\f$ the pressure. The permeability tensor is
/// given by \f$K\f$ and \f$\mu\f$ denotes the viscosity. 
/// 
/// We solve the flow equations for a carthesian grid and we set the source term
/// \f$q\f$ be zero except at the left-lower and right-upper corner, where it is equal 
/// with opposite sign (inflow equal to outflow).


#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <opm/core/grid.h>
#include <opm/core/GridManager.hpp>
#include <opm/core/utility/writeVtkData.hpp>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <opm/core/linalg/LinearSolverUmfpack.hpp>
#include <opm/core/pressure/IncompTpfa.hpp>
#include <opm/core/pressure/FlowBCManager.hpp>
#include <opm/core/utility/miscUtilities.hpp>

/// \page tutorial2
/// \section commentedcode Commented code:
/// \code
int main()
{
    /// \endcode
    /// \page tutorial2
    /// We construct a carthesian grid
    /// \code
    int dim = 3;
    int nx = 40;
    int ny = 40;
    int nz = 1;
    Opm::GridManager grid(nx, ny, nz);
    /// \endcode
    /// \page tutorial2
    /// \details We access the unstructured grid  through
    /// the pointer given by \c grid.c_grid(). For more  details on unstructured
    /// grid, see  grid.h.
    /// \code
    int num_cells = grid.c_grid()->number_of_cells;
    int num_faces = grid.c_grid()->number_of_faces;
    /// \endcode
    /// \page tutorial2
    /// We define the viscosity (unit: cP).
    /// \code
    double mu = 1.0;
    /// \endcode
    /// \page tutorial2
    /// We define the permeability (unit: mD).
    /// \code
    double k = 100.0;
    /// \endcode
    /// \page tutorial2
    /// \details
    /// We set up the permeability tensor and compute the mobility for each cell.
    /// The permeability tensor is flattened in a vector.
    /// \code
    std::vector<double> permeability(num_cells*dim*dim, 0.);
    std::vector<double> mob(num_cells);
    for (int cell = 0; cell < num_cells; ++cell) {
        permeability[9*cell + 0] = k;
        permeability[9*cell + 4] = k;
        permeability[9*cell + 8] = k;
        mob[cell] = 1/mu;
    }
    /// \endcode

    /// \page tutorial2
    /// We choose the UMFPACK linear solver for the pressure solver.
    /// \code
    Opm::LinearSolverUmfpack linsolver;
    /// \endcode
    /// \page tutorial2
    /// We set up the pressure solver
    /// The third argument which corresponds to gravity is set to
    /// zero (no gravity).
    /// \code
    Opm::IncompTpfa psolver(*grid.c_grid(), &permeability[0], 0, linsolver);
    /// \endcode
    /// \page tutorial2
    /// We define the source term.
    /// \code
    std::vector<double> src(num_cells, 0.0);
    src[0] = 100.;
    src[num_cells-1] = -100.;
    /// \endcode
    /// \page tutorial2
    /// \details We set up the boundary conditions. We do not modify them. 
    /// By default, we obtain no outflow boundary conditions.
    /// \code
    /// \code
    Opm::FlowBCManager bcs;
    /// \endcode

    /// \page tutorial2
    /// We declare the solution vectors, i.e., the pressure and face
    /// flux vectors we are going to compute.  The well solution
    /// vectors are needed for interface compatibility with the
    /// <CODE>solve()</CODE> method of class
    /// <CODE>Opm::IncompTPFA</CODE>.
    /// \code
    std::vector<double> pressure(num_cells);
    std::vector<double> faceflux(num_faces);
    std::vector<double> well_bhp;
    std::vector<double> well_flux;
    /// \endcode 
    /// \page tutorial2
    /// \details
    /// We declare the gravity term which is required by the pressure solver (see
    /// Opm::IncompTpfa.solve()). In the absence of gravity, an empty vector is required.
    /// \code
    std::vector<double> omega; 
    /// \endcode
    
    /// \page tutorial2
    /// \details
    /// We declare the wdp term which is required by the pressure solver (see
    /// Opm::IncompTpfa.solve()). In the absence of wells, an empty vector is required.
    /// \code
    std::vector<double> wdp; 
    /// \endcode

    /// \page tutorial2
    /// We call the pressure solver.
    /// \code
    psolver.solve(mob, omega, src, wdp, bcs.c_bcs(),
                  pressure, faceflux, well_bhp, well_flux);
    /// \endcode

    /// \page tutorial2
    /// We write the results in a file in VTK format.
    /// \code
    std::ofstream vtkfile("tutorial2.vtu");
    Opm::DataMap dm;
    dm["pressure"] = &pressure;
    std::vector<double> cell_velocity;
    Opm::estimateCellVelocity(*grid.c_grid(), faceflux, cell_velocity);
    dm["velocity"] = &cell_velocity;
    Opm::writeVtkData(*grid.c_grid(), dm, vtkfile);
}
/// \endcode
/// \page tutorial2
/// We read the the vtu output file in \a Paraview and obtain the following pressure
/// distribution. \image html tutorial2.png


/// \page tutorial2
/// \section sourcecode Complete source code.
/// \include tutorial2.cpp 
