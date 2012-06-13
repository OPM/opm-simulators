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


#include <opm/polymer/IncompTpfaPolymer.hpp>

#include <opm/core/fluid/IncompPropertiesInterface.hpp>
#include <opm/core/fluid/RockCompressibility.hpp>
#include <opm/core/pressure/tpfa/ifs_tpfa.h>
#include <opm/core/pressure/tpfa/trans_tpfa.h>
#include <opm/core/pressure/mimetic/mimetic.h>
#include <opm/core/pressure/flow_bc.h>
#include <opm/core/linalg/LinearSolverInterface.hpp>
#include <opm/core/linalg/sparse_sys.h>
#include <opm/polymer/PolymerState.hpp>
#include <opm/polymer/PolymerUtilities.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/newwells.h>
#include <iomanip>
#include <cmath>
#include <algorithm>

namespace Opm
{



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
    IncompTpfaPolymer::IncompTpfaPolymer(const UnstructuredGrid& grid,
                                         const IncompPropertiesInterface& props,
                                         const RockCompressibility* rock_comp_props,
                                         const PolymerProperties& poly_props,
                                         LinearSolverInterface& linsolver,
                                         const double residual_tol,
                                         const double change_tol,
                                         const int maxiter,
                                         const double* gravity,
                                         const Wells* wells,
                                         const std::vector<double>& src,
                                         const FlowBoundaryConditions* bcs)
        : IncompTpfa(grid, props, rock_comp_props, linsolver,
                     residual_tol, change_tol, maxiter,
                     gravity, wells, src, bcs),
          poly_props_(poly_props),
          c_(0),
          cmax_(0)
    {
    }




    /// Solve the pressure equation. If there is no pressure
    /// dependency introduced by rock compressibility effects,
    /// the equation is linear, and it is solved directly.
    /// Otherwise, the nonlinear equations ares solved by a
    /// Newton-Raphson scheme.
    /// May throw an exception if the number of iterations
    /// exceed maxiter (set in constructor).
    void IncompTpfaPolymer::solve(const double dt,
                                  PolymerState& state,
                                  WellState& well_state)
    {
        c_ = &state.concentration();
        cmax_ = &state.maxconcentration();
        if (rock_comp_props_ != 0 && rock_comp_props_->isActive()) {
            solveRockComp(dt, state.twophaseState(), well_state);
        } else {
            solveIncomp(dt, state.twophaseState(), well_state);
        }
    }







    /// Compute per-solve dynamic properties.
    void IncompTpfaPolymer::computePerSolveDynamicData(const double /*dt*/,
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

        // The only difference from IncompTpfa::computePerSolveDynamicData() is that
        // we call the polymer-aware versions of the computeTotalMobility*() functions.

        // wdp_
        if (wells_) {
            Opm::computeWDP(*wells_, grid_, state.saturation(), props_.density(),
                            gravity_ ? gravity_[2] : 0.0, true, wdp_);
        }
        // totmob_, omega_, gpress_omegaweighted_
        if (gravity_) {
            computeTotalMobilityOmega(props_, poly_props_, allcells_, state.saturation(), *c_, *cmax_,
                                      totmob_, omega_);
            mim_ip_density_update(grid_.number_of_cells, grid_.cell_facepos,
                                  &omega_[0],
                                  &gpress_[0], &gpress_omegaweighted_[0]);
        } else {
            computeTotalMobility(props_, poly_props_, allcells_, state.saturation(), *c_, *cmax_, totmob_);
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







} // namespace Opm
