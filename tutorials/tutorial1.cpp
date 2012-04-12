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






#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

/// \page tutorial1 A simple carthesian grid
/// This tutorial explains how to construct a simple carthesian grid.\n\n
/// We construct a 2x2 two dimensional carthesian grid with 4 blocks of equal size.

#include <opm/core/grid.h>
#include <opm/core/GridManager.hpp>
#include <opm/core/utility/writeVtkData.hpp>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

// ----------------- Main program -----------------

/// \page tutorial1 
/// \code
int main()
{
    /// \endcode
    /// \page tutorial1
    /// By setting <code>nz = 1</code>, we make the grid two dimensional
    /// \code
    int nx = 3;
    int ny = 3;
    int nz = 1;
    /// \endcode
    /// The size of each block is 1x1x1. We use standard units (SI)
    /// \code
    double dx = 1.0;
    double dy = 1.0;
    double dz = 1.0;
    /// \endcode 
    /// \page tutorial1 
    /// One of the constructors of the class Opm::GridManager takes <code>nx,ny,nz,dx,dy,dz</code> 
    /// and construct the corresponding carthesian grid.
    /// \code
    Opm::GridManager grid(nx, ny, nz, dx, dy, dz);
    /// \endcode
    /// \page tutorial1
    /// We open a file to write down the output
    /// \code
    std::ofstream vtkfile("tutorial1.vtu");
    /// \endcode
    /// \page tutorial1
    /// The Opm::writeVtkData() function writes output data. Here, we just want to visualize the
    /// grid. We construct an empty Opm::DataMap object, which we send to Opm::writeVtkData() together with the grid
    /// \code
    Opm::DataMap dm;
    /// \endcode
    /// \page tutorial1
    /// The function Opm::writeVtkData() writes down the output.
    /// \code
    Opm::writeVtkData(*grid.c_grid(), dm, vtkfile);
}
/// \endcode
/// \page tutorial1
/// We read the the vtu output file in \a Paraview and obtain the following grid.
/// \image html tutorial1.png

/// \page tutorial1 
/// \section sourcecode Source code.
/// \include tutorial1.cpp 

