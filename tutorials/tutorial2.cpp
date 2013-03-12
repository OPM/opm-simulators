/*!
\cond SKIP
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
\endcond
*/
/// \page tutorial2 Flow Solver for a single phase
/// \details The flow equations consist of the mass conservation equation
/// \f[\nabla\cdot {\bf u}=q\f] and the Darcy law \f[{\bf u} =- \frac{1}{\mu}K\nabla p.\f] Here,
/// \f${\bf u}\f$ denotes the velocity and \f$p\f$ the pressure. The permeability tensor is
/// given by \f$K\f$ and \f$\mu\f$ denotes the viscosity.
///
/// We solve the flow equations for a Cartesian grid and we set the source term
/// \f$q\f$ be zero except at the left-lower and right-upper corner, where it is equal
/// with opposite sign (inflow equal to outflow).
#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <opm/core/grid.h>
#include <opm/core/GridManager.hpp>
#include <opm/core/io/vtk/writeVtkData.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <opm/core/fluid/IncompPropertiesBasic.hpp>
#include <opm/core/linalg/LinearSolverUmfpack.hpp>
#include <opm/core/pressure/IncompTpfa.hpp>
#include <opm/core/pressure/FlowBCManager.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/simulator/TwophaseState.hpp>
#include <opm/core/simulator/WellState.hpp>

/// \page tutorial2
/// \section commentedcode2 Program walk-through.
/// 

int main()
{
    
    /// \page tutorial2
    /// We construct a Cartesian grid
    /// \snippet tutorial2.cpp cartesian grid
    /// \internal [cartesian grid]
    int dim = 3;
    int nx = 40;
    int ny = 40;
    int nz = 1;
    Opm::GridManager grid(nx, ny, nz);
    /// \internal [cartesian grid]
    /// \endinternal
    /// \page tutorial2
    /// \details We access the unstructured grid  through
    /// the pointer given by \c grid.c_grid(). For more  details on the
    /// UnstructuredGrid data structure, see  grid.h.
    /// \snippet tutorial2.cpp access grid
    /// \internal [access grid]
    int num_cells = grid.c_grid()->number_of_cells;
    int num_faces = grid.c_grid()->number_of_faces;
    /// \internal [access grid]
    /// endinternal


    /// \page tutorial2
    /// \details
    /// We define a fluid viscosity equal to 1 cP and density equal
    /// to 1000 kg/m^3.
    /// The <opm/core/utility/Units.hpp> header contains support
    /// for common units and prefixes, in the namespaces Opm::unit
    /// and Opm::prefix.
    /// \snippet tutorial2.cpp fluid
    /// \internal [fluid]
    using namespace Opm::unit;
    using namespace Opm::prefix;
    int num_phases = 1;
    std::vector<double> mu(num_phases, 1.0*centi*Poise);
    std::vector<double> rho(num_phases, 1000.0*kilogram/cubic(meter));
    /// \internal [fluid]
    /// \endinternal
    /// \page tutorial2
    /// \details
    /// We define a permeability equal to 100 mD.
    /// \snippet tutorial2.cpp perm
    /// \internal [perm]
    double k = 100.0*milli*darcy;
    /// \internal [perm]
    /// \endinternal

    /// \page tutorial2
    /// \details
    /// We set up a simple property object for a single-phase situation.
    /// \snippet tutorial2.cpp single-phase property
    /// \internal [single-phase property]
    Opm::IncompPropertiesBasic props(1, Opm::SaturationPropsBasic::Constant, rho,
                                     mu, 1.0, k, dim, num_cells);
    /// \internal [single-phase property]
    /// /endinternal

    /// \page tutorial2
    /// \details
    /// We take UMFPACK as the linear solver for the pressure solver
    /// (this library has therefore to be installed).
    /// \snippet tutorial2.cpp linsolver
    /// \internal [linsolver] 
    Opm::LinearSolverUmfpack linsolver;
    /// \internal [linsolver] 
    /// \endinternal

    /// \page tutorial2
    /// We define the source term.
    /// \snippet tutorial2.cpp source
    /// \internal [source]
    std::vector<double> src(num_cells, 0.0);
    src[0] = 100.;
    src[num_cells-1] = -100.;
    /// \internal [source]
    /// \endinternal
    
    /// \page tutorial2
    /// \details We set up the boundary conditions.
    /// By default, we obtain no-flow boundary conditions.
    /// \snippet tutorial2.cpp boundary
    /// \internal [boundary]
    Opm::FlowBCManager bcs;
    /// \internal [boundary]
    /// \endinternal

    /// We set up a pressure solver for the incompressible problem,
    /// using the two-point flux approximation discretization.  The
    /// null pointers correspond to arguments for gravity, wells and
    /// boundary conditions, which are all defaulted (to zero gravity,
    /// no wells, and no-flow boundaries).
    /// \snippet tutorial2.cpp tpfa
    /// \internal [tpfa]
    Opm::IncompTpfa psolver(*grid.c_grid(), props, linsolver, NULL, NULL, src, NULL);
    /// \internal [tpfa]
    /// \endinternal
    
    /// \page tutorial2
    /// We declare the state object, that will contain the pressure and face
    /// flux vectors we are going to compute.  The well state
    /// object is needed for interface compatibility with the
    /// <CODE>solve()</CODE> method of class
    /// <CODE>Opm::IncompTPFA</CODE>.
    /// \snippet tutorial2.cpp state
    /// \internal [state]
    Opm::TwophaseState state;
    state.pressure().resize(num_cells, 0.0);
    state.faceflux().resize(num_faces, 0.0);
    state.saturation().resize(num_cells, 1.0);
    Opm::WellState well_state;
    /// \internal [state]
    /// \endinternal

    /// \page tutorial2
    /// We call the pressure solver.
    /// The first (timestep) argument does not matter for this
    /// incompressible case.
    /// \snippet tutorial2.cpp pressure solver
    /// \internal [pressure solver]
    psolver.solve(1.0*day, state, well_state);
    /// \internal [pressure solver]
    /// \endinternal

    /// \page tutorial2
    /// We write the results to a file in VTK format.
    /// The data vectors added to the Opm::DataMap must
    /// contain cell data. They may be a scalar per cell
    /// (pressure) or a vector per cell (cell_velocity).
    /// \snippet tutorial2.cpp write output
    /// \internal [write output]
    std::ofstream vtkfile("tutorial2.vtu");
    Opm::DataMap dm;
    dm["pressure"] = &state.pressure();
    std::vector<double> cell_velocity;
    Opm::estimateCellVelocity(*grid.c_grid(), state.faceflux(), cell_velocity);
    dm["velocity"] = &cell_velocity;
    Opm::writeVtkData(*grid.c_grid(), dm, vtkfile);
    /// \internal [write output]
    /// \endinternal
}

/// \page tutorial2
/// We read the vtu output file in \a Paraview and obtain the following pressure
/// distribution. \image html tutorial2.png


/// \page tutorial2
/// \section completecode2 Complete source code:
/// \include tutorial2.cpp

/// \page tutorial2
/// \details
/// \section pythonscript2 python script to generate figures: 
/// \snippet generate_doc_figures.py tutorial2
