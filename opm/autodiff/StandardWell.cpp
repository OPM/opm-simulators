/*
  Copyright 2017 SINTEF ICT, Applied Mathematics.
  Copyright 2017 Statoil ASA.

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

#include <opm/autodiff/StandardWell.hpp>



namespace Opm
{
    StandardWell::
    StandardWell(const Well* well, const int time_step, const Wells* wells)
    : WellInterface(well, time_step, wells)
    , perf_densities_(number_of_perforations_)
    , perf_pressure_diffs_(number_of_perforations_)
    , well_variables_(numWellEq) // the number of the primary variables
    {
        dune_B_.setBuildMode( Mat::row_wise );
        dune_C_.setBuildMode( Mat::row_wise );
        inv_dune_D_.setBuildMode( Mat::row_wise );
    }





    const std::vector<double>&
    StandardWell::
    perfDensities() const
    {
        return perf_densities_;
    }





    std::vector<double>&
    StandardWell::
    perfDensities()
    {
        return perf_densities_;
    }





    const std::vector<double>&
    StandardWell::
    perfPressureDiffs() const
    {
        return perf_pressure_diffs_;
    }





    std::vector<double>&
    StandardWell::
    perfPressureDiffs()
    {
        return perf_pressure_diffs_;
    }





    void
    StandardWell::
    assembleWellEq(Simulator& ebos_simulator,
                   const double dt,
                   WellState& well_state,
                   bool only_wells)
    {
    }





    void StandardWell::
    setWellVariables(const WellState& well_state)
    {
        const int np = number_of_phases_;
        const int nw = well_state.bhp().size();
        // TODO: it should be the number of primary variables
        // TODO: this is from the old version of StandardWellsDense, it is a coincidence, 3 phases and 3 primary variables
        // TODO: it needs to be careful.
        // TODO: the following code has to be rewritten later for correctness purpose.
        for (int phase = 0; phase < np; ++phase) {
            well_variables_[phase] = 0.0;
            well_variables_[phase].setValue(well_state.wellSolutions()[index_of_well_ + nw * phase]);
            well_variables_[phase].setDerivative(numEq + phase, 1.0);
        }
    }





    StandardWell::EvalWell
    StandardWell::
    getBhp() const
    {
        const WellControls* wc = well_controls_;
        if (well_controls_get_current_type(wc) == BHP) {
            EvalWell bhp = 0.0;
            const double target_rate = well_controls_get_current_target(wc);
            bhp.setValue(target_rate);
            return bhp;
        } else if (well_controls_get_current_type(wc) == THP) {
            const int control = well_controls_get_current(wc);
            const double thp = well_controls_get_current_target(wc);
            const double alq = well_controls_iget_alq(wc, control);
            const int table_id = well_controls_iget_vfp(wc, control);
            EvalWell aqua = 0.0;
            EvalWell liquid = 0.0;
            EvalWell vapour = 0.0;
            EvalWell bhp = 0.0;
            double vfp_ref_depth = 0.0;

            const Opm::PhaseUsage& pu = phaseUsage();

            if (active()[ Water ]) {
                aqua = getQs(pu.phase_pos[ Water]);
            }
            if (active()[ Oil ]) {
                liquid = getQs(pu.phase_pos[ Oil ]);
            }
            if (active()[ Gas ]) {
                vapour = getQs(pu.phase_pos[ Gas ]);
            }
            if (wellType() == INJECTOR) {
                bhp = vfp_properties_->getInj()->bhp(table_id, aqua, liquid, vapour, thp);
                vfp_ref_depth = vfp_properties_->getInj()->getTable(table_id)->getDatumDepth();
            } else {
                bhp = vfp_properties_->getProd()->bhp(table_id, aqua, liquid, vapour, thp, alq);
                vfp_ref_depth = vfp_properties_->getProd()->getTable(table_id)->getDatumDepth();
            }

            // pick the density in the top layer
            const double rho = perf_densities_[0];
            // TODO: not sure whether it is always correct
            const double well_ref_depth = perf_depth_[0];
            // const double dp = wellhelpers::computeHydrostaticCorrection(wells(), wellIdx, vfp_ref_depth, rho, gravity_);
            const double dp = wellhelpers::computeHydrostaticCorrection(well_ref_depth, vfp_ref_depth, rho, gravity_);
            bhp -= dp;
            return bhp;
        }

        return well_variables_[XvarWell];
    }





    StandardWell::EvalWell
    StandardWell::
    getQs(const int phase) const
    {
        EvalWell qs = 0.0;

        return qs; // temporary

        /* const WellControls* wc = well_controls_;
        const int np = number_of_phases_;

        // the target from the well controls
        const double target = well_controls_get_current_target(wc);

        // TODO: the formulation for the injectors decides it only work with single phase
        // surface rate injection control. Improvement will be required.
        // Injectors
        if (wellType() == INJECTOR) {
            // TODO: we should not rely on comp_frac anymore, it should depend on the values in distr
            const double comp_frac = wells().comp_frac[np*wellIdx + phaseIdx];
            if (comp_frac == 0.0) {
                return qs;
            }

            if (well_controls_get_current_type(wc) == BHP || well_controls_get_current_type(wc) == THP) {
                return well_variables_[XvarWell];
            }

            // rate control
            // TODO: if it is reservoir volume rate, it should be wrong here
            qs.setValue(target);
            return qs;
        }


        // Producers
        if (well_controls_get_current_type(wc) == BHP || well_controls_get_current_type(wc) == THP ) {
            return well_variables_[XvarWell] * wellVolumeFractionScaled(phase);
        } */

    }






    StandardWell::EvalWell
    StandardWell::
    wellVolumeFractionScaled(const int phase) const
    {
        // TODO: we should be able to set the g for the well based on the control type
        // instead of using explicit code for g all the times
        const WellControls* wc = well_controls_;
        if (well_controls_get_current_type(wc) == RESERVOIR_RATE) {
            const double* distr = well_controls_get_current_distr(wc);
            if (distr[phase] > 0.) {
                return wellVolumeFraction(phase) / distr[phase];
            } else {
                // TODO: not sure why return EvalWell(0.) causing problem here
                // Probably due to the wrong Jacobians.
                return wellVolumeFraction(phase);
            }
        }
        std::vector<double> g = {1,1,0.01};
        return (wellVolumeFraction(phase) / g[phase]);
    }





    StandardWell::EvalWell
    StandardWell::
    wellVolumeFraction(const int phase) const
    {
        if (phase == Water) {
            return well_variables_[WFrac];
        }

        if (phase == Gas) {
            return well_variables_[GFrac];
        }

        // Oil fraction
        EvalWell well_fraction = 1.0;
        if (active()[Water]) {
            well_fraction -= well_variables_[WFrac];
        }

        if (active()[Gas]) {
            well_fraction -= well_variables_[GFrac];
        }
        return well_fraction;
    }





    StandardWell::EvalWell
    StandardWell::
    wellSurfaceVolumeFraction(const int phase) const
    {
        EvalWell sum_volume_fraction_scaled = 0.;
        const int np = number_of_phases_;
        for (int p = 0; p < np; ++p) {
            sum_volume_fraction_scaled += wellVolumeFractionScaled(p);
        }

        assert(sum_volume_fraction_scaled.value() != 0.);

        return wellVolumeFractionScaled(phase) / sum_volume_fraction_scaled;
     }

}
