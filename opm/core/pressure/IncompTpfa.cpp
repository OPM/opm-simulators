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

#include <opm/core/pressure/IncompTpfa.hpp>
#include <opm/core/pressure/tpfa/ifs_tpfa.h>
#include <opm/core/pressure/tpfa/trans_tpfa.h>
#include <opm/core/pressure/mimetic/mimetic.h>
#include <opm/core/pressure/flow_bc.h>
#include <opm/core/linalg/LinearSolverInterface.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

namespace Opm
{



    /// Construct solver.
    /// \param[in] g             A 2d or 3d grid.
    /// \param[in] permeability  Array of permeability tensors, the array
    ///                          should have size N*D^2, if D == g.dimensions
    ///                          and N == g.number_of_cells.
    /// \param[in] gravity       Gravity vector. If nonzero, the array should
    ///                          have D elements.
    IncompTpfa::IncompTpfa(const UnstructuredGrid& g,
			   const double* permeability,
			   const double* gravity,
                           const LinearSolverInterface& linsolver)
	: grid_(g),
          linsolver_(linsolver),
	  htrans_(g.cell_facepos[ g.number_of_cells ]),
	  trans_ (g.number_of_faces)
    {
	UnstructuredGrid* gg = const_cast<UnstructuredGrid*>(&grid_);
	tpfa_htrans_compute(gg, permeability, &htrans_[0]);
	if (gravity) {
	    gpress_.resize(g.cell_facepos[ g.number_of_cells ], 0.0);

	    mim_ip_compute_gpress(gg->number_of_cells, gg->dimensions, gravity,
				  gg->cell_facepos, gg->cell_faces,
				  gg->face_centroids, gg->cell_centroids,
				  &gpress_[0]);
	}
	// gpress_omegaweighted_ is sent to assembler always, and it dislikes
	// getting a zero pointer.
	gpress_omegaweighted_.resize(g.cell_facepos[ g.number_of_cells ], 0.0);
	h_ = ifs_tpfa_construct(gg, 0);
    }




    /// Destructor.
    IncompTpfa::~IncompTpfa()
    {
	ifs_tpfa_destroy(h_);
    }



    /// Assemble and solve pressure system.
    /// \param[in]  totmob     Must contain N total mobility values (one per cell).
    ///                        totmob = \sum_{p} kr_p/mu_p.
    /// \param[in]  omega      Must be empty if constructor gravity argument was null.
    ///                        Otherwise must contain N mobility-weighted density values (one per cell).
    ///                        omega = \frac{\sum_{p} mob_p rho_p}{\sum_p rho_p}.
    /// \param[in]  src        Must contain N source rates (one per cell).
    ///                        Positive values represent total inflow rates,
    ///                        negative values represent total outflow rates.
    /// \param[in]  bcs        If non-null, specifies boundary conditions.
    ///                        If null, noflow conditions are assumed.
    /// \param[out] pressure   Will contain N cell-pressure values.
    /// \param[out] faceflux   Will contain F signed face flux values.
    void IncompTpfa::solve(const std::vector<double>& totmob,
			   const std::vector<double>& omega,
			   const std::vector<double>& src,
			   const FlowBoundaryConditions* bcs,
			   std::vector<double>& pressure,
			   std::vector<double>& faceflux)
    {
	UnstructuredGrid* gg = const_cast<UnstructuredGrid*>(&grid_);
	tpfa_eff_trans_compute(gg, &totmob[0], &htrans_[0], &trans_[0]);

	if (!omega.empty()) {
	    if (gpress_.empty()) {
		THROW("Nozero omega argument given, but gravity was null in constructor.");
	    }
	    mim_ip_density_update(gg->number_of_cells, gg->cell_facepos,
				  &omega[0],
				  &gpress_[0], &gpress_omegaweighted_[0]);
	} else {
	    if (!gpress_.empty()) {
		THROW("Empty omega argument given, but gravity was non-null in constructor.");
	    }
	}

        const ifs_tpfa_forces F = { &src[0], bcs };

	ifs_tpfa_assemble(gg, &F, &trans_[0], &gpress_omegaweighted_[0], h_);

	linsolver_.solve(h_->A, h_->b, h_->x);

	pressure.resize(grid_.number_of_cells);
	faceflux.resize(grid_.number_of_faces);

	ifs_tpfa_press_flux(gg, &F, &trans_[0], h_,
			    &pressure[0],
			    &faceflux[0]);
    }




} // namespace Opm
