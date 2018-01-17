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

#ifndef OPM_FLOWBCMANAGER_HEADER_INCLUDED
#define OPM_FLOWBCMANAGER_HEADER_INCLUDED

#include <opm/core/pressure/flow_bc.h>

struct UnstructuredGrid;

namespace Opm
{

    /// This class manages a FlowBoundaryConditions struct in the
    /// sense that it encapsulates creation and destruction of the
    /// data structure.
    /// The resulting struct is available through the c_bcs() method.
    class FlowBCManager
    {
    public:
	/// Default constructor sets up empty boundary conditions.
	/// By convention, this is equivalent to all-noflow conditions.
	FlowBCManager();

	/// Destructor.
	~FlowBCManager();

	/// Remove all appended BCs.
	/// By convention, BCs are now equivalent to all-noflow conditions.
	void clear();

	/// Append a single boundary condition.
	/// If the type is BC_NOFLOW the value argument is not used.
	/// If the type is BC_PRESSURE the value argument is a pressure value.
	/// If the type is BC_FLUX_TOTVOL the value argument is a total flux value (m^3/s).
	/// Note: unset boundary conditions are noflow by convention,
	/// so it is normally not necessary to explicitly append
	/// BC_NOFLOW conditions. However, it may make sense to do so
	/// if the bc will change during a simulation run.
	/// Note: if normal velocity bcs are desired, convert to
	/// fluxes by multiplying with face area.
	void append(const FlowBCType type,
		    const int face,
		    const double value);

	/// Defines the canonical sides for logical cartesian grids.
	enum Side { Xmin, Xmax, Ymin, Ymax, Zmin, Zmax };

	/// Add BC_PRESSURE boundary conditions to all faces on a given side.
	/// The grid must have a logical cartesian structure, and grid
	/// faces must be tagged (i.e. grid.cell_facetag must be
	/// non-null). Only the set of faces adjacent to cells with
	/// minimum/maximum I/J/K coordinate (depending on side) are
	/// considered.
	void pressureSide(const UnstructuredGrid& grid,
			  const Side side,
			  const double pressure);

	/// Add BC_FLUX_TOTVOL boundary conditions to all faces on a given side.
	/// The grid must have a logical cartesian structure, and grid
	/// faces must be tagged (i.e. grid.cell_facetag must be
	/// non-null). Only the set of faces adjacent to cells with
	/// minimum/maximum I/J/K coordinate (depending on side) are
	/// considered.
	/// The flux specified is taken to be the total flux through
	/// the side, each individual face receiving a part of the
	/// total flux in proportion to its area, so that all faces
	/// will have identical normal velocities.
	void fluxSide(const UnstructuredGrid& grid,
		      const Side side,
		      const double flux);

	/// Access the managed boundary conditions.
	/// The method is named similarly to c_str() in std::string,
	/// to make it clear that we are returning a C-compatible struct.
	const FlowBoundaryConditions* c_bcs() const;

    private:
	// Disable copying and assignment.
	FlowBCManager(const FlowBCManager& other);
	FlowBCManager& operator=(const FlowBCManager& other);
	// The managed struct.
	FlowBoundaryConditions* bc_;
    };

} // namespace Opm

#endif // OPM_FLOWBCMANAGER_HEADER_INCLUDED
