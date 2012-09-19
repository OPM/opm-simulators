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
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/fluid/RockCompressibility.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <numeric>

namespace Opm
{


    /// Construct solver.
    /// \param[in] grid          A 2d or 3d grid.
    /// \param[in] props         Rock and fluid properties.
    /// \param[in] linsolver     Linear solver to use.
    /// \param[in] residual_tol  Solution accepted if inf-norm of residual is smaller.
    /// \param[in] change_tol    Solution accepted if inf-norm of change in pressure is smaller.
    /// \param[in] maxiter       Maximum acceptable number of iterations.
    /// \param[in] gravity       Gravity vector. If non-null, the array should
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
                                       const RockCompressibility* rock_comp_props,
                                       const LinearSolverInterface& linsolver,
                                       const double residual_tol,
                                       const double change_tol,
                                       const int maxiter,
                                       const double* gravity,
                                       const struct Wells* wells)
        : grid_(grid),
          props_(props),
          rock_comp_props_(rock_comp_props),
          linsolver_(linsolver),
          residual_tol_(residual_tol),
          change_tol_(change_tol),
          maxiter_(maxiter),
          gravity_(gravity),
          wells_(wells),
          htrans_(grid.cell_facepos[ grid.number_of_cells ]),
          trans_ (grid.number_of_faces),
          allcells_(grid.number_of_cells),
          singular_(false)
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
        // If we have rock compressibility, pore volumes are updated
        // in the compute*() methods, otherwise they are constant and
        // hence may be computed here.
        if (rock_comp_props_ == NULL || !rock_comp_props_->isActive()) {
            computePorevolume(grid_, props.porosity(), porevol_);
        }
        for (int c = 0; c < grid.number_of_cells; ++c) {
            allcells_[c] = c;
        }
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
        const int nc = grid_.number_of_cells;
        const int nw = wells_->number_of_wells;

        // Set up dynamic data.
        computePerSolveDynamicData(dt, state, well_state);
        computePerIterationDynamicData(dt, state, well_state);

        // Assemble J and F.
        assemble(dt, state, well_state);

        double inc_norm = 0.0;
        int iter = 0;
        double res_norm = residualNorm();
        std::cout << "\nIteration         Residual        Change in p\n"
                  << std::setw(9) << iter
                  << std::setw(18) << res_norm
                  << std::setw(18) << '*' << std::endl;
        while ((iter < maxiter_) && (res_norm > residual_tol_)) {
            // Solve for increment in Newton method:
            //   incr = x_{n+1} - x_{n} = -J^{-1}F
            // (J is Jacobian matrix, F is residual)
            solveIncrement();
            ++iter;

            // Update pressure vars with increment.
            for (int c = 0; c < nc; ++c) {
                state.pressure()[c] += pressure_increment_[c];
            }
            for (int w = 0; w < nw; ++w) {
                well_state.bhp()[w] += pressure_increment_[nc + w];
            }

            // Stop iterating if increment is small.
            inc_norm = incrementNorm();
            if (inc_norm <= change_tol_) {
                std::cout << std::setw(9) << iter
                          << std::setw(18) << '*'
                          << std::setw(18) << inc_norm << std::endl;
                break;
            }

            // Set up dynamic data.
            computePerIterationDynamicData(dt, state, well_state);

            // Assemble J and F.
            assemble(dt, state, well_state);

            // Update residual norm.
            res_norm = residualNorm();

            std::cout << std::setw(9) << iter
                      << std::setw(18) << res_norm
                      << std::setw(18) << inc_norm << std::endl;
        }

        if ((iter == maxiter_) && (res_norm > residual_tol_) && (inc_norm > change_tol_)) {
            THROW("CompressibleTpfa::solve() failed to converge in " << maxiter_ << " iterations.");
        }

        std::cout << "Solved pressure in " << iter << " iterations." << std::endl;

        // Compute fluxes and face pressures.
        computeResults(state, well_state);
    }





    /// @brief After solve(), was the resulting pressure singular.
    /// Returns true if the pressure is singular in the following
    /// sense: if everything is incompressible and there are no
    /// pressure conditions, the absolute values of the pressure
    /// solution are arbitrary. (But the differences in pressure
    /// are significant.)
    bool CompressibleTpfa::singularPressure() const
    {
        return singular_;
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
        wellperf_gpot_.clear();
        wellperf_gpot_.resize(np*nperf, 0.0);
        if (grav == 0.0) {
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
    void CompressibleTpfa::computePerSolveDynamicData(const double /*dt*/,
                                                      const BlackoilState& state,
                                                      const WellState& /*well_state*/)
    {
        computeWellPotentials(state);
        if (rock_comp_props_ && rock_comp_props_->isActive()) {
            computePorevolume(grid_, props_.porosity(), *rock_comp_props_, state.pressure(), initial_porevol_);
        }
    }




    /// Compute per-iteration dynamic properties.
    void CompressibleTpfa::computePerIterationDynamicData(const double dt,
                                                          const BlackoilState& state,
                                                          const WellState& well_state)
    {
        // These are the variables that get computed by this function:
        //
        // std::vector<double> cell_A_;
        // std::vector<double> cell_dA_;
        // std::vector<double> cell_viscosity_;
        // std::vector<double> cell_phasemob_;
        // std::vector<double> cell_voldisc_;
        // std::vector<double> face_A_;
        // std::vector<double> face_phasemob_;
        // std::vector<double> face_gravcap_;
        // std::vector<double> wellperf_A_;
        // std::vector<double> wellperf_phasemob_;
        // std::vector<double> porevol_;   // Only modified if rock_comp_props_ is non-null.
        // std::vector<double> rock_comp_; // Empty unless rock_comp_props_ is non-null.
        computeCellDynamicData(dt, state, well_state);
        computeFaceDynamicData(dt, state, well_state);
        computeWellDynamicData(dt, state, well_state);
    }





    /// Compute per-iteration dynamic properties for cells.
    void CompressibleTpfa::computeCellDynamicData(const double /*dt*/,
                                                  const BlackoilState& state,
                                                  const WellState& /*well_state*/)
    {
        // These are the variables that get computed by this function:
        //
        // std::vector<double> cell_A_;
        // std::vector<double> cell_dA_;
        // std::vector<double> cell_viscosity_;
        // std::vector<double> cell_phasemob_;
        // std::vector<double> cell_voldisc_;
        // std::vector<double> porevol_;   // Only modified if rock_comp_props_ is non-null.
        // std::vector<double> rock_comp_; // Empty unless rock_comp_props_ is non-null.
        const int nc = grid_.number_of_cells;
        const int np = props_.numPhases();
        const double* cell_p = &state.pressure()[0];
        const double* cell_z = &state.surfacevol()[0];
        const double* cell_s = &state.saturation()[0];
        cell_A_.resize(nc*np*np);
        cell_dA_.resize(nc*np*np);
        props_.matrix(nc, cell_p, cell_z, &allcells_[0], &cell_A_[0], &cell_dA_[0]);
        cell_viscosity_.resize(nc*np);
        props_.viscosity(nc, cell_p, cell_z, &allcells_[0], &cell_viscosity_[0], 0);
        cell_phasemob_.resize(nc*np);
        props_.relperm(nc, cell_s, &allcells_[0], &cell_phasemob_[0], 0);
        std::transform(cell_phasemob_.begin(), cell_phasemob_.end(),
                       cell_viscosity_.begin(),
                       cell_phasemob_.begin(),
                       std::divides<double>());
        // Volume discrepancy: we have that
        //     z = Au, voldiscr = sum(u) - 1,
        // but I am not sure it is actually needed.
        // Use zero for now.
        // TODO: Check this!
        cell_voldisc_.clear();
        cell_voldisc_.resize(nc, 0.0);

        if (rock_comp_props_ && rock_comp_props_->isActive()) {
            computePorevolume(grid_, props_.porosity(), *rock_comp_props_, state.pressure(), porevol_);
            rock_comp_.resize(nc);
            for (int cell = 0; cell < nc; ++cell) {
                rock_comp_[cell] = rock_comp_props_->rockComp(state.pressure()[cell]);
            }
        }
    }




    /// Compute per-iteration dynamic properties for faces.
    void CompressibleTpfa::computeFaceDynamicData(const double /*dt*/,
                                                  const BlackoilState& state,
                                                  const WellState& /*well_state*/)
    {
        // These are the variables that get computed by this function:
        //
        // std::vector<double> face_A_;
        // std::vector<double> face_phasemob_;
        // std::vector<double> face_gravcap_;
        const int np = props_.numPhases();
        const int nf = grid_.number_of_faces;
        const int dim = grid_.dimensions;
        const double grav = gravity_ ? gravity_[dim - 1] : 0.0;
        std::vector<double> gravcontrib[2];
        std::vector<double> pot[2];
        gravcontrib[0].resize(np);
        gravcontrib[1].resize(np);
        pot[0].resize(np);
        pot[1].resize(np);
        face_A_.resize(nf*np*np);
        face_phasemob_.resize(nf*np);
        face_gravcap_.resize(nf*np);
        for (int face = 0; face < nf; ++face) {
            // Obtain properties from both sides of the face.
            const double face_depth = grid_.face_centroids[face*dim + dim - 1];
            const int* c = &grid_.face_cells[2*face];

            // Get pressures and compute gravity contributions,
            // to decide upwind directions.
            double c_press[2];
            for (int j = 0; j < 2; ++j) {
                if (c[j] >= 0) {
                    // Pressure
                    c_press[j] = state.pressure()[c[j]];
                    // Gravity contribution, gravcontrib = rho*(face_z - cell_z) [per phase].
                    if (grav != 0.0) {
                        const double depth_diff = face_depth - grid_.cell_centroids[c[j]*dim + dim - 1];
                        props_.density(1, &cell_A_[np*np*c[j]], &gravcontrib[j][0]);
                        for (int p = 0; p < np; ++p) {
                            gravcontrib[j][p] *= depth_diff*grav;
                        }
                    } else {
                        std::fill(gravcontrib[j].begin(), gravcontrib[j].end(), 0.0);
                    }
                } else {
                    // Pressures
                    c_press[j] = state.facepressure()[face];
                    // Gravity contribution.
                    std::fill(gravcontrib[j].begin(), gravcontrib[j].end(), 0.0);
                }
            }

            // Gravity contribution:
            //    gravcapf = rho_1*g*(z_12 - z_1) - rho_2*g*(z_12 - z_2)
            // where _1 and _2 refers to two neigbour cells, z is the
            // z coordinate of the centroid, and z_12 is the face centroid.
            // Also compute the potentials.
            for (int phase = 0; phase < np; ++phase) {
                face_gravcap_[np*face + phase] = gravcontrib[0][phase] - gravcontrib[1][phase];
                pot[0][phase] = c_press[0] + face_gravcap_[np*face + phase];
                pot[1][phase] = c_press[1];
            }

            // Now we can easily find the upwind direction for every phase,
            // we can also tell which boundary faces are inflow bdys.

            // Get upwind mobilities by phase.
            // Get upwind A matrix rows by phase.
            // NOTE:
            // We should be careful to upwind the R factors,
            // the B factors are not that vital.
            //      z = Au = RB^{-1}u,
            // where (this example is for gas-oil)
            //      R = [1 RgL; RoV 1], B = [BL 0 ; 0 BV]
            // (RgL is gas in Liquid phase, RoV is oil in Vapour phase.)
            //      A = [1/BL RgL/BV; RoV/BL 1/BV]
            // This presents us with a dilemma, as V factors should be
            // upwinded according to V phase flow, same for L. What then
            // about the RgL/BV and RoV/BL numbers?
            // We give priority to R, and therefore upwind the rows of A
            // by phase (but remember, Fortran matrix ordering).
            // This prompts the question if we should split the matrix()
            // property method into formation volume and R-factor methods.
            for (int phase = 0; phase < np; ++phase) {
                int upwindc = -1;
                if (c[0] >=0 && c[1] >= 0) {
                    upwindc = (pot[0] < pot[1]) ? c[1] : c[0];
                } else {
                    upwindc = (c[0] >= 0) ? c[0] : c[1];
                }
                face_phasemob_[np*face + phase] = cell_phasemob_[np*upwindc + phase];
                for (int p2 = 0; p2 < np; ++p2) {
                    // Recall: column-major ordering.
                    face_A_[np*np*face + phase + np*p2]
                        = cell_A_[np*np*upwindc + phase + np*p2];
                }
            }
        }
    }




    /// Compute per-iteration dynamic properties for wells.
    void CompressibleTpfa::computeWellDynamicData(const double /*dt*/,
                                                  const BlackoilState& /*state*/,
                                                  const WellState& well_state)
    {
        // These are the variables that get computed by this function:
        //
        // std::vector<double> wellperf_A_;
        // std::vector<double> wellperf_phasemob_;
        const int np = props_.numPhases();
        const int nw = wells_->number_of_wells;
        const int nperf = wells_->well_connpos[nw];
        wellperf_A_.resize(nperf*np*np);
        wellperf_phasemob_.resize(nperf*np);
        // The A matrix is set equal to the perforation grid cells'
        // matrix for producers, computed from bhp and injection
        // component fractions from
        // The mobilities are set equal to the perforation grid cells'
        // mobilities for producers.
        std::vector<double> mu(np);
        for (int w = 0; w < nw; ++w) {
            bool producer = (wells_->type[w] == PRODUCER);
            const double* comp_frac = &wells_->comp_frac[np*w];
            for (int j = wells_->well_connpos[w]; j < wells_->well_connpos[w+1]; ++j) {
                const int c = wells_->well_cells[j];
                double* wpA = &wellperf_A_[np*np*j];
                double* wpM = &wellperf_phasemob_[np*j];
                if (producer) {
                    const double* cA = &cell_A_[np*np*c];
                    std::copy(cA, cA + np*np, wpA);
                    const double* cM = &cell_phasemob_[np*c];
                    std::copy(cM, cM + np, wpM);
                } else {
                    const double bhp = well_state.bhp()[w];
                    double perf_p = bhp;
                    for (int phase = 0; phase < np; ++phase) {
                        perf_p += wellperf_gpot_[np*j + phase]*comp_frac[phase];
                    }
                    // Hack warning: comp_frac is used as a component
                    // surface-volume variable in calls to matrix() and
                    // viscosity(), but as a saturation in the call to
                    // relperm(). This is probably ok as long as injectors
                    // only inject pure fluids.
                    props_.matrix(1, &perf_p, comp_frac, &c, wpA, NULL);
                    props_.viscosity(1, &perf_p, comp_frac, &c, &mu[0], NULL);
                    ASSERT(std::fabs(std::accumulate(comp_frac, comp_frac + np, 0.0) - 1.0) < 1e-6);
                    props_.relperm  (1, comp_frac, &c, wpM , NULL);
                    for (int phase = 0; phase < np; ++phase) {
                        wpM[phase] /= mu[phase];
                    }
                }
            }
        }
    }




    /// Compute the residual and Jacobian.
    void CompressibleTpfa::assemble(const double dt,
                                    const BlackoilState& state,
                                    const WellState& well_state)
    {
        const double* cell_press = &state.pressure()[0];
        const double* well_bhp = well_state.bhp().empty() ? NULL : &well_state.bhp()[0];
        const double* z = &state.surfacevol()[0];
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
        cq.nphases = props_.numPhases();
        cq.Ac = &cell_A_[0];
        cq.dAc = &cell_dA_[0];
        cq.Af = &face_A_[0];
        cq.phasemobf = &face_phasemob_[0];
        cq.voldiscr = &cell_voldisc_[0];
        int was_adjusted = 0;
        if (! (rock_comp_props_ && rock_comp_props_->isActive())) {
            was_adjusted =
                cfs_tpfa_res_assemble(gg, dt, &forces, z, &cq, &trans_[0],
                                      &face_gravcap_[0], cell_press, well_bhp,
                                      &porevol_[0], h_);
        } else {
            was_adjusted =
                cfs_tpfa_res_comprock_assemble(gg, dt, &forces, z, &cq, &trans_[0],
                                               &face_gravcap_[0], cell_press, well_bhp,
                                               &porevol_[0], &initial_porevol_[0],
                                               &rock_comp_[0], h_);
        }
        singular_ = (was_adjusted == 1);
    }




    /// Computes pressure_increment_.
    void CompressibleTpfa::solveIncrement()
    {
        // Increment is equal to -J^{-1}F
        linsolver_.solve(h_->J, h_->F, &pressure_increment_[0]);
        std::transform(pressure_increment_.begin(), pressure_increment_.end(),
                       pressure_increment_.begin(), std::negate<double>());
    }




    namespace {
        template <class FI>
        double infnorm(FI beg, FI end)
        {
            double norm = 0.0;
            for (; beg != end; ++beg) {
                norm = std::max(norm, std::fabs(*beg));
            }
            return norm;
        }
    } // anonymous namespace




    /// Computes the inf-norm of the residual.
    double CompressibleTpfa::residualNorm() const
    {
        const int ndof = pressure_increment_.size();
        return infnorm(h_->F, h_->F + ndof);
    }




    /// Computes the inf-norm of pressure_increment_.
    double CompressibleTpfa::incrementNorm() const
    {
        return infnorm(pressure_increment_.begin(), pressure_increment_.end());
    }




    /// Compute the output.
    void CompressibleTpfa::computeResults(BlackoilState& state,
                                          WellState& well_state) const
    {
        UnstructuredGrid* gg = const_cast<UnstructuredGrid*>(&grid_);
        CompletionData completion_data;
        completion_data.gpot = const_cast<double*>(&wellperf_gpot_[0]);
        completion_data.A = const_cast<double*>(&wellperf_A_[0]);
        completion_data.phasemob = const_cast<double*>(&wellperf_phasemob_[0]);
        cfs_tpfa_res_wells wells_tmp;
        wells_tmp.W = const_cast<Wells*>(wells_);
        wells_tmp.data = &completion_data;
        cfs_tpfa_res_forces forces;
        forces.wells = &wells_tmp;
        forces.src = NULL;

        cfs_tpfa_res_flux(gg,
                          &forces,
                          props_.numPhases(),
                          &trans_[0],
                          &cell_phasemob_[0],
                          &face_phasemob_[0],
                          &face_gravcap_[0],
                          &state.pressure()[0],
                          &well_state.bhp()[0],
                          &state.faceflux()[0],
                          &well_state.perfRates()[0]);
        cfs_tpfa_res_fpress(gg,
                            props_.numPhases(),
                            &htrans_[0],
                            &face_phasemob_[0],
                            &face_gravcap_[0],
                            h_,
                            &state.pressure()[0],
                            &state.faceflux()[0],
                            &state.facepressure()[0]);
    }

} // namespace Opm
