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
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/newwells.h>
#include <opm/core/BlackoilState.hpp>
#include <opm/core/WellState.hpp>

#include <algorithm>

namespace Opm
{


    /// Construct solver.
    /// \param[in] grid          A 2d or 3d grid.
    /// \param[in] props         Rock and fluid properties.
    /// \param[in] linsolver     Linear solver to use.
    /// \param[in] gravity       Gravity vector. If nonzero, the array should
    ///                          have D elements.
    /// \param[in] wells         The wells argument. Will be used in solution, 
    ///                          is ignored if NULL.
    ///                          Note: this class observes the well object, and
    ///                                makes the assumption that the well topology
    ///                                and completions does not change during the
    ///                                run. However, controls (only) are allowed
    ///                                to change.
    CompressibleTpfa::CompressibleTpfa(const UnstructuredGrid& grid,
                                       const BlackoilPropertiesInterface& props,
                                       const LinearSolverInterface& linsolver,
                                       const double* gravity,
                                       const struct Wells* wells)
	: grid_(grid),
          props_(props),
          linsolver_(linsolver),
          gravity_(gravity),
          wells_(wells),
	  htrans_(grid.cell_facepos[ grid.number_of_cells ]),
	  trans_ (grid.number_of_faces)
    {
        if (wells_ && (wells_->number_of_phases != props.numPhases())) {
            THROW("Inconsistent number of phases specified (wells vs. props): "
                  << wells_->number_of_phases << " != " << props.numPhases());
        }
        const int num_dofs = grid.number_of_cells + (wells ? wells->number_of_wells : 0);
        pressure_increment_.resize(num_dofs);
	UnstructuredGrid* gg = const_cast<UnstructuredGrid*>(&grid_);
	tpfa_htrans_compute(gg, props.permeability(), &htrans_[0]);
	tpfa_trans_compute(gg, &htrans_[0], &trans_[0]);
        computePorevolume(grid_, props.porosity(), porevol_);
        cfs_tpfa_res_wells w;
        w.W = const_cast<struct Wells*>(wells_);
        w.data = NULL;
	h_ = cfs_tpfa_res_construct(gg, &w, props.numPhases());
    }




    /// Destructor.
    CompressibleTpfa::~CompressibleTpfa()
    {
	cfs_tpfa_res_destroy(h_);
    }




    /// Solve pressure equation, by Newton iterations.
    void CompressibleTpfa::solve(const double dt,
                                 BlackoilState& state,
                                 WellState& well_state)
    {
        // Set up dynamic data.
        computePerSolveDynamicData(state);
        computePerIterationDynamicData();

        // Assemble J and F.
        assemble(dt, state, well_state);

        bool residual_ok = false; // Replace with tolerance check.
        while (!residual_ok) {
            // Solve for increment in Newton method:
            //   incr = x_{n+1} - x_{n} = -J^{-1}F
            // (J is Jacobian matrix, F is residual)
            solveIncrement();

            // Update pressure vars with increment.

            // Set up dynamic data.
            computePerIterationDynamicData();

            // Assemble J and F.
            assemble(dt, state, well_state);

            // Check for convergence.
            // Include both tolerance check for residual
            // and solution change.
        }

        // Write to output parameters.
        // computeResults(...);
    }





    /// Compute well potentials.
    void CompressibleTpfa::computeWellPotentials(const BlackoilState& state)
    {
        if (wells_ == NULL) return;

        const int nw = wells_->number_of_wells;
        const int np = props_.numPhases();
        const int nperf = wells_->well_connpos[nw];
        const int dim = grid_.dimensions;
        const double grav = gravity_ ? gravity_[dim - 1] : 0.0;
        if (grav == 0.0) {
            wellperf_gpot_.clear();
            wellperf_gpot_.resize(np*nperf, 0.0);
            return;
        }
        // Temporary storage for perforation A matrices and densities.
        std::vector<double> A(np*np, 0.0);
        std::vector<double> rho(np, 0.0);

        // Main loop, iterate over all perforations,
        // using the following formula (by phase):
        //    gpot(perf) = g*(perf_z - well_ref_z)*rho(perf)
        // where the phase densities rho(perf) are taken to be
        // the densities in the perforation cell.
        for (int w = 0; w < nw; ++w) {
            const double ref_depth = wells_->depth_ref[w];
            for (int j = wells_->well_connpos[w]; j < wells_->well_connpos[w + 1]; ++j) {
                const int cell = wells_->well_cells[j];
                const double cell_depth = grid_.cell_centroids[dim * cell + dim - 1];
                props_.matrix(1, &state.pressure()[cell], &state.surfacevol()[np*cell], &cell, &A[0], 0);
                props_.density(1, &A[0], &rho[0]);
                for (int phase = 0; phase < np; ++phase) {
                    wellperf_gpot_[np*j + phase] = rho[phase]*grav*(cell_depth - ref_depth);
                }
            }
        }
    }




    /// Compute per-solve dynamic properties.
    void CompressibleTpfa::computePerSolveDynamicData(const BlackoilState& state)
    {
        computeWellPotentials(state);
    }




    /// Compute per-iteration dynamic properties.
    void CompressibleTpfa::computePerIterationDynamicData()
    {
        // std::vector<double> face_gravcap_;
        // std::vector<double> wellperf_A_;
        // std::vector<double> wellperf_phasemob_;
        // std::vector<double> cell_A_;
        // std::vector<double> cell_dA_;
        // std::vector<double> face_A_;
        // std::vector<double> face_phasemob_;
        // std::vector<double> cell_voldisc_;
    }




    /// Compute the residual and Jacobian.
    void CompressibleTpfa::assemble(const double dt,
                                    const BlackoilState& state,
                                    const WellState& well_state)
    {
        const double* z = &state.surfacevol()[0];
        const double* cell_press = &state.pressure()[0];
        const double* well_bhp = well_state.bhp().empty() ? NULL : &well_state.bhp()[0];
	UnstructuredGrid* gg = const_cast<UnstructuredGrid*>(&grid_);

        CompletionData completion_data;
        completion_data.gpot = &wellperf_gpot_[0];
        completion_data.A = &wellperf_A_[0];
        completion_data.phasemob = &wellperf_phasemob_[0];
        cfs_tpfa_res_wells wells_tmp;
        wells_tmp.W = const_cast<Wells*>(wells_);
        wells_tmp.data = &completion_data;
        cfs_tpfa_res_forces forces;
        forces.wells = &wells_tmp;
        forces.src = NULL; // Check if it is legal to leave it as NULL.
        compr_quantities_gen cq;
        cq.Ac = &cell_A_[0];
        cq.dAc = &cell_dA_[0];
        cq.Af = &face_A_[0];
        cq.phasemobf = &face_phasemob_[0];
        cq.voldiscr = &cell_voldisc_[0];
        cfs_tpfa_res_assemble(gg, dt, &forces, z, &cq, &trans_[0],
                              &face_gravcap_[0], cell_press, well_bhp,
                              &porevol_[0], h_);
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
