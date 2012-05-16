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


#include <opm/core/pressure/CompressibleTpfa.hpp>
#include <opm/core/pressure/tpfa/cfs_tpfa_residual.h>
#include <opm/core/pressure/tpfa/compr_quant_general.h>
#include <opm/core/pressure/tpfa/compr_source.h>
#include <opm/core/pressure/tpfa/trans_tpfa.h>
#include <opm/core/linalg/LinearSolverInterface.hpp>
#include <opm/core/linalg/sparse_sys.h>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/newwells.h>
#include <opm/core/BlackoilState.hpp>

#include <algorithm>

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
    CompressibleTpfa::CompressibleTpfa(const UnstructuredGrid& g,
                                       const double* permeability,
                                       const double* /* gravity */, // ???
                                       const LinearSolverInterface& linsolver,
                                       const struct Wells* wells,
                                       const int num_phases)
	: grid_(g),
          linsolver_(linsolver),
	  htrans_(g.cell_facepos[ g.number_of_cells ]),
	  trans_ (g.number_of_faces),
          wells_(wells)
    {
        if (wells_ && (wells_->number_of_phases != num_phases)) {
            THROW("Inconsistent number of phases specified: "
                  << wells_->number_of_phases << " != " << num_phases);
        }
        const int num_dofs = g.number_of_cells + (wells ? wells->number_of_wells : 0);
        pressure_increment_.resize(num_dofs);
	UnstructuredGrid* gg = const_cast<UnstructuredGrid*>(&grid_);
	tpfa_htrans_compute(gg, permeability, &htrans_[0]);
	tpfa_trans_compute(gg, &htrans_[0], &trans_[0]);
        cfs_tpfa_res_wells w;
        w.W = const_cast<struct Wells*>(wells_);
        w.data = NULL;
	h_ = cfs_tpfa_res_construct(gg, &w, num_phases);
    }




    /// Destructor.
    CompressibleTpfa::~CompressibleTpfa()
    {
	cfs_tpfa_res_destroy(h_);
    }




    /// Solve pressure equation, by Newton iterations.
    void CompressibleTpfa::solve(const double dt,
                                 BlackoilState& state)
    {
        // Set up dynamic data.
        computeDynamicData();

        // Assemble J and F.
        assemble(dt, state);

        bool residual_ok = false; // Replace with tolerance check.
        while (!residual_ok) {
            // Solve for increment in Newton method:
            //   incr = x_{n+1} - x_{n} = -J^{-1}F
            // (J is Jacobian matrix, F is residual)
            solveIncrement();

            // Update pressure vars with increment.

            // Set up dynamic data.
            computeDynamicData();

            // Assemble J and F.
            assemble(dt, state);

            // Check for convergence.
            // Include both tolerance check for residual
            // and solution change.
        }

        // Write to output parameters.
        // computeResults(...);
    }




    /// Solve pressure equation, by Newton iterations.
    void CompressibleTpfa::computeDynamicData()
    {
    }




    /// Compute the residual and Jacobian.
    void CompressibleTpfa::assemble(const double dt,
                                    const BlackoilState& state)
    {
        // Arguments or members?
        const double* z = &state.surfacevol()[0];
        const double* cell_press = &state.pressure()[0];
        const double* well_bhp = NULL;
        const double* porevol = NULL;

	UnstructuredGrid* gg = const_cast<UnstructuredGrid*>(&grid_);
        CompletionData completion_data;
        completion_data.gpot = 0; // TODO
        completion_data.A = 0; // TODO
        completion_data.phasemob = 0; // TODO
        cfs_tpfa_res_wells wells_tmp;
        wells_tmp.W = const_cast<Wells*>(wells_);
        wells_tmp.data = &completion_data;
        cfs_tpfa_res_forces forces;
        forces.wells = &wells_tmp;
        forces.src = NULL; // Check if it is legal to leave it as NULL.
        compr_quantities_gen cq;
        cq.Ac = 0; // TODO
        cq.dAc = 0; // TODO
        cq.Af = 0; // TODO
        cq.phasemobf = 0; // TODO
        cq.voldiscr = 0; // TODO
        // TODO: gravcapf_ must be set already.
        cfs_tpfa_res_assemble(gg, dt, &forces, z, &cq, &trans_[0],
                              &gravcapf_[0], &cell_press[0], &well_bhp[0],
                              &porevol[0], h_);
    }




    /// Computes pressure_increment_.
    void CompressibleTpfa::solveIncrement()
    {
        // Increment is equal to -J^{-1}F
	linsolver_.solve(h_->J, h_->F, &pressure_increment_[0]);        
        std::transform(pressure_increment_.begin(), pressure_increment_.end(),
                       pressure_increment_.begin(), std::negate<double>());
    }




    /// Compute the output.
    void CompressibleTpfa::computeResults(std::vector<double>& // pressure
                                          ,
                                          std::vector<double>& // faceflux
                                          ,
                                          std::vector<double>& // well_bhp
                                          ,
                                          std::vector<double>& // well_rate
                                          )
    {
    }

} // namespace Opm
