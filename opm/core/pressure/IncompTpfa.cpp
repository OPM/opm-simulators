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
#include <opm/core/linalg/sparse_sys.h>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/newwells.h>

namespace Opm
{



    /// Construct solver.
    /// \param[in] g             A 2d or 3d grid.
    /// \param[in] permeability  Array of permeability tensors, the array
    ///                          should have size N*D^2, if D == g.dimensions
    ///                          and N == g.number_of_cells.
    /// \param[in] gravity       Gravity vector. If nonzero, the array should
    ///                          have D elements.
    /// \param[in] wells         The wells argument. Will be used in solution, 
    ///                          is ignored if NULL
    IncompTpfa::IncompTpfa(const UnstructuredGrid& g,
			   const double* permeability,
			   const double* gravity,
                           const LinearSolverInterface& linsolver,
                           const struct Wells* wells)
	: grid_(g),
          linsolver_(linsolver),
	  htrans_(g.cell_facepos[ g.number_of_cells ]),
	  trans_ (g.number_of_faces),
          wells_(wells)
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
    /// \param[out] well_bhp   Will contain bhp values for each well passed
    ///                        in the constructor.
    /// \param[out] well_rate  Will contain rate values for each well passed
    ///                        in the constructor.
    void IncompTpfa::solve(const std::vector<double>& totmob,
			   const std::vector<double>& omega,
			   const std::vector<double>& src,
			   const FlowBoundaryConditions* bcs,
			   std::vector<double>& pressure,
			   std::vector<double>& faceflux,
                           std::vector<double>& well_bhp,
                           std::vector<double>& well_rate)
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

        ifs_tpfa_forces F = { NULL, NULL, wells_, NULL, NULL };
        if (! src.empty()) { F.src = &src[0]; }
        F.bc = bcs;

	ifs_tpfa_assemble(gg, &F, &trans_[0], &gpress_omegaweighted_[0], h_);

	linsolver_.solve(h_->A, h_->b, h_->x);

	pressure.resize(grid_.number_of_cells);
	faceflux.resize(grid_.number_of_faces);

        ifs_tpfa_solution soln = { NULL, NULL, NULL, NULL };
        soln.cell_press = &pressure[0];
        soln.face_flux  = &faceflux[0];

        if(wells_ != NULL) {
            well_bhp.resize(wells_->number_of_wells);
            well_rate.resize(wells_->number_of_wells);
            soln.well_flux = &well_rate[0];
            soln.well_press = &well_bhp[0];
        }
        
        ifs_tpfa_press_flux(gg, &F, &trans_[0], h_, &soln);
    }


    /// Assemble and solve pressure system with rock compressibility (assumed constant per cell).
    /// \param[in]  totmob     Must contain N total mobility values (one per cell).
    ///                        totmob = \sum_{p} kr_p/mu_p.
    /// \param[in]  omega      Must be empty if constructor gravity argument was null.
    ///                        Otherwise must contain N fractional-flow-weighted density
    ///                        values (one per cell).
    /// \param[in]  src        Must contain N source rates (one per cell).
    ///                        Positive values represent total inflow rates,
    ///                        negative values represent total outflow rates.
    /// \param[in]  bcs        If non-null, specifies boundary conditions.
    ///                        If null, noflow conditions are assumed.
    /// \param[in]  porevol    Must contain N pore volumes.
    /// \param[in]  rock_comp  Must contain N rock compressibilities.
    ///                        rock_comp = (d poro / d p)*(1/poro).
    /// \param[in]  dt         Timestep.
    /// \param[out] pressure   Will contain N cell-pressure values.
    /// \param[out] faceflux   Will contain F signed face flux values.
    /// \param[out] well_bhp   Will contain bhp values for each well passed
    ///                        in the constructor
    /// \param[out] well_rate  Will contain rate values for each well passed
    ///                        in the constructor
    void IncompTpfa::solve(const std::vector<double>& totmob,
                           const std::vector<double>& omega,
                           const std::vector<double>& src,
                           const FlowBoundaryConditions* bcs,
                           const std::vector<double>& porevol,
                           const std::vector<double>& rock_comp,
                           const double dt,
                           std::vector<double>& pressure,
                           std::vector<double>& faceflux,
                           std::vector<double>& well_bhp,
                           std::vector<double>& well_rate)
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

        ifs_tpfa_forces F = { NULL, NULL, wells_, NULL, NULL };
        if (! src.empty()) { F.src = &src[0]; }
        F.bc = bcs;

	ifs_tpfa_assemble(gg, &F, &trans_[0], &gpress_omegaweighted_[0], h_);

        // TODO: this is a hack, it would be better to handle this in a
        // (variant of) ifs_tpfa_assemble().
        if (!rock_comp.empty()) {
            // We must compensate for adjustment made in ifs_tpfa_assemble()
            // to make the system nonsingular.
            h_->A->sa[0] *= 0.5;

            // The extra term of the equation is
            //
            //     porevol*rock_comp*(p - p0)/dt.
            //
            // The p part goes on the diagonal, the p0 on the rhs.
            for (int c = 0; c < gg->number_of_cells; ++c) {
                // Find diagonal
                size_t j = csrmatrix_elm_index(c, c, h_->A);
                double d = porevol[c] * rock_comp[c] / dt;
                h_->A->sa[j] += d;
                h_->b[c]     += d * pressure[c];
            }
        }

	linsolver_.solve(h_->A, h_->b, h_->x);

	pressure.resize(grid_.number_of_cells);
	faceflux.resize(grid_.number_of_faces);

        ifs_tpfa_solution soln = { NULL, NULL, NULL, NULL };
        soln.cell_press = &pressure[0];
        soln.face_flux  = &faceflux[0];

        if(wells_ != NULL) {
            well_bhp.resize(wells_->number_of_wells);
            well_rate.resize(wells_->number_of_wells);
            soln.well_flux = &well_rate[0];
            soln.well_press = &well_bhp[0];
        }
        ifs_tpfa_press_flux(gg, &F, &trans_[0], h_, &soln);
    }



} // namespace Opm
