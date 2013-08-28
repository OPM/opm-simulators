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
#include <opm/core/pressure/FlowBCManager.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/grid.h>
#include <vector>

namespace Opm
{

    namespace
    {
	std::string sideString(FlowBCManager::Side s);
	void findSideFaces(const UnstructuredGrid& grid,
			   const FlowBCManager::Side side,
			   std::vector<int>& faces);
    } // anon namespace


    /// Default constructor sets up empty boundary conditions.
    /// By convention, this is equivalent to all-noflow conditions.
    FlowBCManager::FlowBCManager()
	: bc_(0)
    {
	bc_ = flow_conditions_construct(0);
	if (!bc_) {
	    OPM_THROW(std::runtime_error, "Failed to construct FlowBoundaryConditions struct.");
	}
    }


    /// Destructor.
    FlowBCManager::~FlowBCManager()
    {
	flow_conditions_destroy(bc_);
    }


    /// Remove all appended BCs.
    /// By convention, BCs are now equivalent to all-noflow conditions.
    void FlowBCManager::clear()
    {
	flow_conditions_clear(bc_);
    }


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
    void FlowBCManager::append(const FlowBCType type,
			       const int face,
			       const double value)
    {
	int ok = flow_conditions_append(type, face, value, bc_);
	if (!ok) {
	    OPM_THROW(std::runtime_error, "Failed to append boundary condition for face " << face);
	}
    }


    /// Add BC_PRESSURE boundary conditions to all faces on a given side.
    /// The grid must have a logical cartesian structure, and grid
    /// faces must be tagged (i.e. grid.cell_facetag must be
    /// non-null). Only the set of faces adjacent to cells with
    /// minimum/maximum I/J/K coordinate (depending on side) are
    /// considered.
    void FlowBCManager::pressureSide(const UnstructuredGrid& grid,
				     const Side side,
				     const double pressure)
    {
	std::vector<int> faces;
	findSideFaces(grid, side, faces);
	int ok = flow_conditions_append_multi(BC_PRESSURE, faces.size(), &faces[0], pressure, bc_);
	if (!ok) {
	    OPM_THROW(std::runtime_error, "Failed to append pressure boundary conditions for side " << sideString(side));
	}
    }


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
    void FlowBCManager::fluxSide(const UnstructuredGrid& grid,
				 const Side side,
				 const double flux)
    {
	// Find side faces.
	std::vector<int> faces;
	findSideFaces(grid, side, faces);

	// Compute total area of faces.
	double tot_area = 0.0;
	for (int fi = 0; fi < int(faces.size()); ++fi) {
	    tot_area += grid.face_areas[faces[fi]];
	}

	// Append flux conditions for all the faces individually.
	for (int fi = 0; fi < int(faces.size()); ++fi) {
	    const double face_flux = flux * grid.face_areas[faces[fi]] / tot_area;
	    int ok = flow_conditions_append(BC_FLUX_TOTVOL, faces[fi], face_flux, bc_);
	    if (!ok) {
		OPM_THROW(std::runtime_error, "Failed to append flux boundary conditions for face " << faces[fi] << " on side " << sideString(side));
	    }
	}
    }



    /// Access the managed boundary conditions.
    /// The method is named similarly to c_str() in std::string,
    /// to make it clear that we are returning a C-compatible struct.
    const FlowBoundaryConditions* FlowBCManager::c_bcs() const
    {
	return bc_;
    }




    // ------ Utility functions ------


    namespace
    {
	std::string sideString(FlowBCManager::Side s)
	{
	    switch (s) {
	    case FlowBCManager::Xmin: return "Xmin";
	    case FlowBCManager::Xmax: return "Xmax";
	    case FlowBCManager::Ymin: return "Ymin";
	    case FlowBCManager::Ymax: return "Ymax";
	    case FlowBCManager::Zmin: return "Zmin";
	    case FlowBCManager::Zmax: return "Zmax";
	    default: OPM_THROW(std::runtime_error, "Unknown side tag " << s);
	    }
	}



	void cartCoord(const int ndims,
                       const int log_cart_coord,
		       const int* dims,
		       int* ijk)
	{
            int ix = log_cart_coord;

            for (int dim = 0; dim < ndims; ++dim) {
                ijk[dim]  = ix % dims[dim];
                ix       /=      dims[dim];
            }

            ASSERT2 (ix == 0,
                     "Lexicographic index is not consistent "
                     "with grid dimensions.");
	}



	/// The grid must have a logical cartesian structure, and grid
	/// faces must be tagged (i.e. grid.cell_facetag must be
	/// non-null). Only the set of faces adjacent to cells with
	/// minimum/maximum I/J/K coordinate (depending on side) are
	/// considered.
	void findSideFaces(const UnstructuredGrid& grid,
			   const FlowBCManager::Side side,
			   std::vector<int>& faces)
	{
	    if (grid.cell_facetag == 0) {
		OPM_THROW(std::runtime_error, "Faces not tagged - cannot extract " << sideString(side) << " faces.");
	    }

            ASSERT2 (grid.dimensions <= 3,
                     "Grid must have three dimensions or less.");

            ASSERT2 (side < 2 * grid.dimensions,
                     "Boundary condition side not consistent with "
                     "number of physical grid dimensions.");

	    // Get all boundary faces with the correct tag and with
	    // min/max i/j/k (depending on side).
	    const int correct_ijk = (side % 2) ? grid.cartdims[side/2] - 1 : 0;
	    for (int c = 0; c < grid.number_of_cells; ++c) {
		int ijk[3] = { -1, -1, -1 };
                int gc = (grid.global_cell != 0) ? grid.global_cell[c] : c;
		cartCoord(grid.dimensions, gc, grid.cartdims, ijk);
		if (ijk[side/2] != correct_ijk) {
		    continue;
		}
		for (int hf = grid.cell_facepos[c]; hf < grid.cell_facepos[c + 1]; ++hf) {
		    if (grid.cell_facetag[hf] == side) {
			// Tag is correct.
			const int f = grid.cell_faces[hf];
			if (grid.face_cells[2*f] == -1 || grid.face_cells[2*f + 1] == -1) {
			    // Face is on boundary.
			    faces.push_back(f);
			} else {
			    OPM_THROW(std::runtime_error, "Face not on boundary, even with correct tag and boundary cell. This should not occur.");
			}
		    }
		}
	    }
	}

    } // anon namespace

} // namespace Opm
