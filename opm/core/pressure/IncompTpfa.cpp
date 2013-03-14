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

#include <opm/core/props/IncompPropertiesInterface.hpp>
#include <opm/core/props/rock/RockCompressibility.hpp>
#include <opm/core/pressure/tpfa/ifs_tpfa.h>
#include <opm/core/pressure/tpfa/trans_tpfa.h>
#include <opm/core/pressure/mimetic/mimetic.h>
#include <opm/core/pressure/flow_bc.h>
#include <opm/core/linalg/LinearSolverInterface.hpp>
#include <opm/core/linalg/sparse_sys.h>
#include <opm/core/simulator/TwophaseState.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/newwells.h>
#include <iomanip>
#include <cmath>
#include <algorithm>

namespace Opm
{



    /// Construct solver for incompressible case.
    /// \param[in] grid             A 2d or 3d grid.
    /// \param[in] props            Rock and fluid properties.
    /// \param[in] linsolver        Linear solver to use.
    /// \param[in] gravity          Gravity vector. If non-null, the array should
    ///                             have D elements.
    /// \param[in] wells            The wells argument. Will be used in solution,
    ///                             is ignored if NULL.
    ///                             Note: this class observes the well object, and
    ///                                   makes the assumption that the well topology
    ///                                   and completions does not change during the
    ///                                   run. However, controls (only) are allowed
    ///                                   to change.
    /// \param[in] src              Source terms. May be empty().
    /// \param[in] bcs              Boundary conditions, treat as all noflow if null.
    IncompTpfa::IncompTpfa(const UnstructuredGrid& grid,
                           const IncompPropertiesInterface& props,
                           LinearSolverInterface& linsolver,
                           const double* gravity,
                           const Wells* wells,
                           const std::vector<double>& src,
                           const FlowBoundaryConditions* bcs)
        : grid_(grid),
          props_(props),
          rock_comp_props_(NULL),
          linsolver_(linsolver),
          residual_tol_(0.0),
          change_tol_(0.0),
          maxiter_(0),
          gravity_(gravity),
          wells_(wells),
          src_(src),
          bcs_(bcs),
          htrans_(grid.cell_facepos[ grid.number_of_cells ]),
          allcells_(grid.number_of_cells),
          trans_ (grid.number_of_faces)
    {
        computeStaticData();
    }




    /// Construct solver, possibly with rock compressibility.
    /// \param[in] grid             A 2d or 3d grid.
    /// \param[in] props            Rock and fluid properties.
    /// \param[in] rock_comp_props  Rock compressibility properties. May be null.
    /// \param[in] linsolver        Linear solver to use.
    /// \param[in] residual_tol     Solution accepted if inf-norm of residual is smaller.
    /// \param[in] change_tol       Solution accepted if inf-norm of change in pressure is smaller.
    /// \param[in] maxiter          Maximum acceptable number of iterations.
    /// \param[in] gravity          Gravity vector. If non-null, the array should
    ///                             have D elements.
    /// \param[in] wells            The wells argument. Will be used in solution,
    ///                             is ignored if NULL.
    ///                             Note: this class observes the well object, and
    ///                                   makes the assumption that the well topology
    ///                                   and completions does not change during the
    ///                                   run. However, controls (only) are allowed
    ///                                   to change.
    /// \param[in] src              Source terms. May be empty().
    /// \param[in] bcs              Boundary conditions, treat as all noflow if null.
    IncompTpfa::IncompTpfa(const UnstructuredGrid& grid,
                           const IncompPropertiesInterface& props,
                           const RockCompressibility* rock_comp_props,
                           LinearSolverInterface& linsolver,
                           const double residual_tol,
                           const double change_tol,
                           const int maxiter,
                           const double* gravity,
                           const Wells* wells,
                           const std::vector<double>& src,
                           const FlowBoundaryConditions* bcs)
        : grid_(grid),
          props_(props),
          rock_comp_props_(rock_comp_props),
          linsolver_(linsolver),
          residual_tol_(residual_tol),
          change_tol_(change_tol),
          maxiter_(maxiter),
          gravity_(gravity),
          wells_(wells),
          src_(src),
          bcs_(bcs),
          htrans_(grid.cell_facepos[ grid.number_of_cells ]),
          allcells_(grid.number_of_cells),
          trans_ (grid.number_of_faces)
    {
        computeStaticData();
    }





    /// Destructor.
    IncompTpfa::~IncompTpfa()
    {
        ifs_tpfa_destroy(h_);
    }




    /// Solve the pressure equation. If there is no pressure
    /// dependency introduced by rock compressibility effects,
    /// the equation is linear, and it is solved directly.
    /// Otherwise, the nonlinear equations ares solved by a
    /// Newton-Raphson scheme.
    /// May throw an exception if the number of iterations
    /// exceed maxiter (set in constructor).
    void IncompTpfa::solve(const double dt,
                           TwophaseState& state,
                           WellState& well_state)
    {
        if (rock_comp_props_ != 0 && rock_comp_props_->isActive()) {
            solveRockComp(dt, state, well_state);
        } else {
            solveIncomp(dt, state, well_state);
        }
    }



    // Solve with no rock compressibility (linear eqn).
    void IncompTpfa::solveIncomp(const double dt,
                                 TwophaseState& state,
                                 WellState& well_state)
    {
        // Set up properties.
        computePerSolveDynamicData(dt, state, well_state);

        // Assemble.
        UnstructuredGrid* gg = const_cast<UnstructuredGrid*>(&grid_);
        int ok = ifs_tpfa_assemble(gg, &forces_, &trans_[0], &gpress_omegaweighted_[0], h_);
        if (!ok) {
            THROW("Failed assembling pressure system.");
        }

        // Solve.
        linsolver_.solve(h_->A, h_->b, h_->x);

        // Obtain solution.
        ASSERT(int(state.pressure().size()) == grid_.number_of_cells);
        ASSERT(int(state.faceflux().size()) == grid_.number_of_faces);
        ifs_tpfa_solution soln = { NULL, NULL, NULL, NULL };
        soln.cell_press = &state.pressure()[0];
        soln.face_flux  = &state.faceflux()[0];
        if (wells_ != NULL) {
            ASSERT(int(well_state.bhp().size()) == wells_->number_of_wells);
            ASSERT(int(well_state.perfRates().size()) == wells_->well_connpos[ wells_->number_of_wells ]);
            soln.well_flux = &well_state.perfRates()[0];
            soln.well_press = &well_state.bhp()[0];
        }
        ifs_tpfa_press_flux(gg, &forces_, &trans_[0], h_, &soln);
    }






    // Solve with rock compressibility (nonlinear eqn).
    void IncompTpfa::solveRockComp(const double dt,
                                   TwophaseState& state,
                                   WellState& well_state)
    {
        // This function is identical to CompressibleTpfa::solve().
        // \TODO refactor?

        const int nc = grid_.number_of_cells;
        const int nw = (wells_) ? wells_->number_of_wells : 0;

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
                state.pressure()[c] += h_->x[c];
            }
            for (int w = 0; w < nw; ++w) {
                well_state.bhp()[w] += h_->x[nc + w];
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
            THROW("IncompTpfa::solve() failed to converge in " << maxiter_ << " iterations.");
        }

        std::cout << "Solved pressure in " << iter << " iterations." << std::endl;

        // Compute fluxes and face pressures.
        computeResults(state, well_state);
    }






    /// Compute data that never changes (after construction).
    void IncompTpfa::computeStaticData()
    {
        if (wells_ && (wells_->number_of_phases != props_.numPhases())) {
            THROW("Inconsistent number of phases specified (wells vs. props): "
                  << wells_->number_of_phases << " != " << props_.numPhases());
        }
        const int num_dofs = grid_.number_of_cells + (wells_ ? wells_->number_of_wells : 0);
        pressures_.resize(num_dofs);
        UnstructuredGrid* gg = const_cast<UnstructuredGrid*>(&grid_);
        tpfa_htrans_compute(gg, props_.permeability(), &htrans_[0]);
        if (gravity_) {
            gpress_.resize(gg->cell_facepos[ gg->number_of_cells ], 0.0);

            mim_ip_compute_gpress(gg->number_of_cells, gg->dimensions, gravity_,
                                  gg->cell_facepos, gg->cell_faces,
                                  gg->face_centroids, gg->cell_centroids,
                                  &gpress_[0]);
        }
        // gpress_omegaweighted_ is sent to assembler always, and it dislikes
        // getting a zero pointer.
        gpress_omegaweighted_.resize(gg->cell_facepos[ gg->number_of_cells ], 0.0);
        if (rock_comp_props_) {
            rock_comp_.resize(grid_.number_of_cells);
        }
        for (int c = 0; c < grid_.number_of_cells; ++c) {
            allcells_[c] = c;
        }
        h_ = ifs_tpfa_construct(gg, const_cast<struct Wells*>(wells_));
    }






    /// Compute per-solve dynamic properties.
    void IncompTpfa::computePerSolveDynamicData(const double /*dt*/,
                                                const TwophaseState& state,
                                                const WellState& /*well_state*/)
    {
        // Computed here:
        //
        // std::vector<double> wdp_;
        // std::vector<double> totmob_;
        // std::vector<double> omega_;
        // std::vector<double> trans_;
        // std::vector<double> gpress_omegaweighted_;
        // std::vector<double> initial_porevol_;
        // ifs_tpfa_forces forces_;

        // wdp_
        if (wells_) {
            Opm::computeWDP(*wells_, grid_, state.saturation(), props_.density(),
                            gravity_ ? gravity_[2] : 0.0, true, wdp_);
        }
        // totmob_, omega_, gpress_omegaweighted_
        if (gravity_) {
            computeTotalMobilityOmega(props_, allcells_, state.saturation(), totmob_, omega_);
            mim_ip_density_update(grid_.number_of_cells, grid_.cell_facepos,
                                  &omega_[0],
                                  &gpress_[0], &gpress_omegaweighted_[0]);
        } else {
            computeTotalMobility(props_, allcells_, state.saturation(), totmob_);
        }
        // trans_
        tpfa_eff_trans_compute(const_cast<UnstructuredGrid*>(&grid_), &totmob_[0], &htrans_[0], &trans_[0]);
        // initial_porevol_
        if (rock_comp_props_ && rock_comp_props_->isActive()) {
            computePorevolume(grid_, props_.porosity(), *rock_comp_props_, state.pressure(), initial_porevol_);
        }
        // forces_
        forces_.src = src_.empty() ? NULL : &src_[0];
        forces_.bc = bcs_;
        forces_.W = wells_;
        forces_.totmob = &totmob_[0];
        forces_.wdp = wdp_.empty() ? NULL : &wdp_[0];
    }






    /// Compute per-iteration dynamic properties.
    void IncompTpfa::computePerIterationDynamicData(const double /*dt*/,
                                                    const TwophaseState& state,
                                                    const WellState& well_state)
    {
        // These are the variables that get computed by this function:
        //
        // std::vector<double> porevol_
        // std::vector<double> rock_comp_
        // std::vector<double> pressures_

        computePorevolume(grid_, props_.porosity(), *rock_comp_props_, state.pressure(), porevol_);
        if (rock_comp_props_ && rock_comp_props_->isActive()) {
            for (int cell = 0; cell < grid_.number_of_cells; ++cell) {
                rock_comp_[cell] = rock_comp_props_->rockComp(state.pressure()[cell]);
            }
        }
        if (wells_) {
            std::copy(state.pressure().begin(), state.pressure().end(), pressures_.begin());
            std::copy(well_state.bhp().begin(), well_state.bhp().end(), pressures_.begin() + grid_.number_of_cells);
        }
    }





    /// Compute the residual in h_->b and Jacobian in h_->A.
    void IncompTpfa::assemble(const double dt,
                              const TwophaseState& state,
                              const WellState& /*well_state*/)
    {
        const double* pressures = wells_ ? &pressures_[0] : &state.pressure()[0];

        bool ok = ifs_tpfa_assemble_comprock_increment(const_cast<UnstructuredGrid*>(&grid_),
                                                       &forces_, &trans_[0], &gpress_omegaweighted_[0],
                                                       &porevol_[0], &rock_comp_[0], dt, pressures,
                                                       &initial_porevol_[0], h_);
        if (!ok) {
            THROW("Failed assembling pressure system.");
        }
    }






    /// Computes pressure increment, puts it in h_->x
    void IncompTpfa::solveIncrement()
    {
        // Increment is equal to -J^{-1}R.
        // The Jacobian is in h_->A, residual in h_->b.
        linsolver_.solve(h_->A, h_->b, h_->x);
        // It is not necessary to negate the increment,
        // apparently the system for the increment is generated,
        // not the Jacobian and residual as such.
        // std::transform(h_->x, h_->x + h_->A->m, h_->x, std::negate<double>());
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
    double IncompTpfa::residualNorm() const
    {
        return infnorm(h_->b, h_->b + h_->A->m);
    }




    /// Computes the inf-norm of pressure_increment_.
    double IncompTpfa::incrementNorm() const
    {
        return infnorm(h_->x, h_->x + h_->A->m);
    }




    /// Compute the output.
    void IncompTpfa::computeResults(TwophaseState& state,
                                    WellState& well_state) const
    {
        // Make sure h_ contains the direct-solution matrix
        // and right hand side (not jacobian and residual).
        // TODO: optimize by only adjusting b and diagonal of A.
        UnstructuredGrid* gg = const_cast<UnstructuredGrid*>(&grid_);
        ifs_tpfa_assemble(gg, &forces_, &trans_[0], &gpress_omegaweighted_[0], h_);


        // Make sure h_->x contains the direct solution vector.
        ASSERT(int(state.pressure().size()) == grid_.number_of_cells);
        ASSERT(int(state.faceflux().size()) == grid_.number_of_faces);
        std::copy(state.pressure().begin(), state.pressure().end(), h_->x);
        std::copy(well_state.bhp().begin(), well_state.bhp().end(), h_->x + grid_.number_of_cells);

        // Obtain solution.
        ifs_tpfa_solution soln = { NULL, NULL, NULL, NULL };
        soln.cell_press = &state.pressure()[0];
        soln.face_flux  = &state.faceflux()[0];
        if (wells_ != NULL) {
            ASSERT(int(well_state.bhp().size()) == wells_->number_of_wells);
            ASSERT(int(well_state.perfRates().size()) == wells_->well_connpos[ wells_->number_of_wells ]);
            soln.well_flux = &well_state.perfRates()[0];
            soln.well_press = &well_state.bhp()[0];
        }
        ifs_tpfa_press_flux(gg, &forces_, &trans_[0], h_, &soln); // TODO: Check what parts of h_ are used here.
    }


} // namespace Opm
