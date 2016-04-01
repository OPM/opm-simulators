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

#include <config.h>

#include <opm/polymer/PolymerBlackoilState.hpp>
#include <opm/polymer/CompressibleTpfaPolymer.hpp>
#include <opm/core/props/BlackoilPropertiesInterface.hpp>
#include <opm/core/props/rock/RockCompressibility.hpp>
#include <opm/core/pressure/tpfa/ifs_tpfa.h>
#include <opm/core/pressure/tpfa/trans_tpfa.h>
#include <opm/core/pressure/mimetic/mimetic.h>
#include <opm/core/pressure/flow_bc.h>
#include <opm/core/linalg/LinearSolverInterface.hpp>
#include <opm/core/linalg/sparse_sys.h>
#include <opm/polymer/polymerUtilities.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/wells.h>
#include <iomanip>
#include <cmath>
#include <algorithm>

namespace Opm
{



    /// Construct solver
    /// \param[in] grid             A 2d or 3d grid.
    /// \param[in] props            Rock and fluid properties.
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
    CompressibleTpfaPolymer::CompressibleTpfaPolymer(const UnstructuredGrid& grid,
                                                     const BlackoilPropertiesInterface& props,
                                                     const RockCompressibility* rock_comp_props,
                                                     const PolymerProperties& poly_props,
                                                     const LinearSolverInterface& linsolver,
                                                     const double residual_tol,
                                                     const double change_tol,
                                                     const int maxiter,
                                                     const double* gravity,
                                                     const Wells* wells)
        : CompressibleTpfa(grid, props, rock_comp_props, linsolver,
                           residual_tol, change_tol, maxiter,
                           gravity, wells),
          poly_props_(poly_props),
          c_(0),
          cmax_(0)
    {
    }




    /// Solve the pressure equation. The nonlinear equations ares
    /// solved by a Newton-Raphson scheme.  May throw an exception if
    /// the number of iterations exceed maxiter (set in constructor).
    void CompressibleTpfaPolymer::solve(const double dt,
                                  PolymerBlackoilState& state,
                                  WellState& well_state)
    {
        c_ = &state.getCellData( state.CONCENTRATION );
        cmax_ = &state.getCellData( state.CMAX );
        CompressibleTpfa::solve(dt, state, well_state);
    }

    /// Compute per-solve dynamic properties.
    void CompressibleTpfaPolymer::computePerSolveDynamicData(const double /* dt */,
                                                             const BlackoilState& state,
                                                             const WellState& /* well_state */)
    {
        // std::vector<double> cell_relperm__;
        // std::vector<double> cell_eff_relperm_;
        const int nc = grid_.number_of_cells;
        const int np = props_.numPhases();
        cell_relperm_.resize(nc*np);
        cell_eff_viscosity_.resize(nc*np);
        const double* cell_s = &state.saturation()[0];
        props_.relperm(nc, cell_s, &allcells_[0], &cell_relperm_[0], 0);
        computeWellPotentials(state);
        if (rock_comp_props_ && rock_comp_props_->isActive()) {
            computePorevolume(grid_, props_.porosity(), *rock_comp_props_, state.pressure(), initial_porevol_);
        }
    }


    /// Compute per-iteration dynamic properties for cells.
    void CompressibleTpfaPolymer::computeCellDynamicData(const double /*dt*/,
                                                  const BlackoilState& state,
                                                  const WellState& /*well_state*/)
    {
        // These are the variables that get computed by this function:
        //
        // std::vector<double> cell_A_;
        // std::vector<double> cell_dA_;
        // std::vector<double> cell_viscosity_;
        // std::vector<double> cell_eff_viscosity_;
        // std::vector<double> cell_phasemob_;
        // std::vector<double> cell_voldisc_;
        // std::vector<double> porevol_;   // Only modified if rock_comp_props_ is non-null.
        // std::vector<double> rock_comp_; // Empty unless rock_comp_props_ is non-null.
        const int nc = grid_.number_of_cells;
        const int np = props_.numPhases();
        const double* cell_p = &state.pressure()[0];
        const double* cell_T = &state.temperature()[0];
        const double* cell_z = &state.surfacevol()[0];
        cell_A_.resize(nc*np*np);
        cell_dA_.resize(nc*np*np);
        props_.matrix(nc, cell_p, cell_T, cell_z, &allcells_[0], &cell_A_[0], &cell_dA_[0]);
        cell_viscosity_.resize(nc*np);
        props_.viscosity(nc, cell_p, cell_T, cell_z, &allcells_[0], &cell_viscosity_[0], 0);
        cell_phasemob_.resize(nc*np);
        for (int cell = 0; cell < nc; ++cell) {
            poly_props_.effectiveVisc((*c_)[cell], &cell_viscosity_[np*cell + 0], cell_eff_viscosity_[np*cell + 0]);
            poly_props_.effectiveMobilities((*c_)[cell], (*cmax_)[cell], &cell_viscosity_[np*cell + 0], &cell_relperm_[np*cell + 0], &cell_phasemob_[np*cell + 0]);
        }

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






} // namespace Opm
