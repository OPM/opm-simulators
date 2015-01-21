/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.

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
#include <opm/core/pressure/IncompTpfaSinglePhase.hpp>

#include <opm/core/props/IncompPropertiesSinglePhase.hpp>
#include <opm/core/pressure/tpfa/ifs_tpfa.h>
#include <opm/core/pressure/tpfa/trans_tpfa.h>
// #include <opm/core/pressure/mimetic/mimetic.h>
// #include <opm/core/pressure/flow_bc.h>
#include <opm/core/linalg/LinearSolverInterface.hpp>
#include <opm/core/linalg/sparse_sys.h>
// #include <opm/core/simulator/TwophaseState.hpp>
// #include <opm/core/simulator/WellState.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
// #include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/wells.h>
// #include <iostream>
// #include <iomanip>
// #include <cmath>
// #include <algorithm>

namespace Opm
{



    /// Construct solver for incompressible case.
    /// \param[in] grid             A 2d or 3d grid.
    /// \param[in] props            Rock and fluid properties.
    /// \param[in] linsolver        Linear solver to use.
    /// \param[in] wells            The wells used as driving forces.
    IncompTpfaSinglePhase::IncompTpfaSinglePhase(const UnstructuredGrid& grid,
                                                 const IncompPropertiesSinglePhase& props,
                                                 const LinearSolverInterface& linsolver,
                                                 const Wells& wells)
        : grid_(grid),
          props_(props),
          linsolver_(linsolver),
          wells_(wells),
          htrans_(grid.cell_facepos[ grid.number_of_cells ]),
          trans_ (grid.number_of_faces),
          zeros_(grid.cell_facepos[ grid.number_of_cells ])
    {
        computeStaticData();
    }






    /// Destructor.
    IncompTpfaSinglePhase::~IncompTpfaSinglePhase()
    {
        ifs_tpfa_destroy(h_);
    }






    /// Solve the pressure equation.
    void IncompTpfaSinglePhase::solve(std::vector<double>& press,
                                      std::vector<double>& flux,
                                      std::vector<double>& bhp,
                                      std::vector<double>& wellrates)
    {
        // Set up properties.
        computePerSolveDynamicData();

        // Assemble.
        UnstructuredGrid* gg = const_cast<UnstructuredGrid*>(&grid_);
        int ok = ifs_tpfa_assemble(gg, &forces_, trans_.data(), zeros_.data(), h_);
        if (!ok) {
            OPM_THROW(std::runtime_error, "Failed assembling pressure system.");
        }

        // Solve.
        linsolver_.solve(h_->A, h_->b, h_->x);

        // Obtain solution.
        press.resize(grid_.number_of_cells);
        flux.resize(grid_.number_of_faces);
        wellrates.resize(wells_.well_connpos[ wells_.number_of_wells ]);
        bhp.resize(wells_.number_of_wells);
        ifs_tpfa_solution soln = { NULL, NULL, NULL, NULL };
        soln.cell_press = press.data();
        soln.face_flux  = flux.data();
        soln.well_press = bhp.data();
        soln.well_flux = wellrates.data();
        ifs_tpfa_press_flux(gg, &forces_, &trans_[0], h_, &soln);
    }






    /// Compute data that never changes (after construction).
    void IncompTpfaSinglePhase::computeStaticData()
    {
        UnstructuredGrid* gg = const_cast<UnstructuredGrid*>(&grid_);
        tpfa_htrans_compute(gg, props_.permeability(), &htrans_[0]);
        h_ = ifs_tpfa_construct(gg, const_cast<struct Wells*>(&wells_));
    }






    /// Compute per-solve dynamic properties.
    void IncompTpfaSinglePhase::computePerSolveDynamicData()
    {
        // Computed here:
        //
        // std::vector<double> totmob_;
        // std::vector<double> trans_;
        // ifs_tpfa_forces forces_;

        // totmob_
        totmob_.clear();
        totmob_.resize(grid_.number_of_cells, 1.0/(*props_.viscosity()));
        // trans_
        tpfa_eff_trans_compute(const_cast<UnstructuredGrid*>(&grid_), totmob_.data(), htrans_.data(), trans_.data());
        // forces_
        forces_.src = NULL;
        forces_.bc = NULL;
        forces_.W = &wells_;
        forces_.totmob = totmob_.data();
        forces_.wdp = zeros_.data();
    }


} // namespace Opm
