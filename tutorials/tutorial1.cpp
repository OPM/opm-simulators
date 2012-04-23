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

/// \page tutorial1 A simple cartesian grid
/// This tutorial explains how to construct a simple cartesian grid,
/// and we will take a look at some output facilities.

/// \page tutorial1
/// \section commentedsource1 Program walkthrough.
/// All headers from opm-core are found in the opm/core/ directory.
/// Some important headers are at the root, other headers are found
/// in subdirectories.
#include <opm/core/grid.h>
#include <opm/core/GridManager.hpp>
#include <opm/core/utility/writeVtkData.hpp>
#include <iostream>
#include <fstream>
#include <vector>

/**
\code
#include <opm/core/grid.h>
#include <opm/core/GridManager.hpp>
#include <opm/core/utility/writeVtkData.hpp>
#include <iostream>
#include <fstream>
#include <vector>
\endcode
*/

// ----------------- Main program -----------------

int main()
{
    /// \page tutorial1
    /// We set the number of blocks in each direction.
    /// \code
    int nx = 4;
    int ny = 3;
    int nz = 2;
    /// \endcode
    /// The size of each block is 1m x 1m x 1m. The default units are always the
    /// standard units (SI). But other units can easily be dealt with, see Opm::unit.
    /// \code
    double dx = 1.0;
    double dy = 1.0;
    double dz = 1.0;
    /// \endcode
    /// \page tutorial1
    /// In opm-core, grid information is accessed via the UnstructuredGrid data structure.
    /// This data structure has a pure C API, including helper functions to construct and
    /// destroy the data structure. In this tutorial however, we will use Opm::GridManager,
    /// which is a C++ class that wraps the UnstructuredGrid and takes care of
    /// object lifetime issues.
    /// One of the constructors of the class Opm::GridManager takes <code>nx, ny, nz, dx, dy, dz</code>
    /// and construct the corresponding cartesian grid.
    /// \code
    Opm::GridManager grid(nx, ny, nz, dx, dy, dz);
    /// \endcode
    /// \page tutorial1
    /// We open an output file stream for the output
    /// \code
    std::ofstream vtkfile("tutorial1.vtu");
    /// \endcode
    /// \page tutorial1
    /// The Opm::writeVtkData() function writes a grid together with
    /// data to a stream. Here, we just want to visualize the grid. We
    /// construct an empty Opm::DataMap object, which we send to
    /// Opm::writeVtkData() together with the grid
    /// \code
    Opm::DataMap dm;
    /// \endcode
    /// \page tutorial1
    /// Call Opm::writeVtkData() to write the output file.
    /// \code
    Opm::writeVtkData(*grid.c_grid(), dm, vtkfile);
}
/// \endcode
/// \page tutorial1
/// We read the vtu output file in \a Paraview and obtain the following grid.
/// \image html tutorial1.png

/// \page tutorial1
/// \section completecode1 Complete source code:
/// \include tutorial1.cpp

