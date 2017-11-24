/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 Statoil ASA.

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


#ifndef OPM_WELLHELPERS_HEADER_INCLUDED
#define OPM_WELLHELPERS_HEADER_INCLUDED

#include <opm/core/wells.h>
// #include <opm/autodiff/AutoDiffHelpers.hpp>

#include <vector>

namespace Opm {




    namespace wellhelpers
    {


        inline
        double rateToCompare(const std::vector<double>& well_phase_flow_rate,
                             const int well,
                             const int num_phases,
                             const double* distr)
        {
            double rate = 0.0;
            for (int phase = 0; phase < num_phases; ++phase) {
                // Important: well_phase_flow_rate is ordered with all phase rates for first
                // well first, then all phase rates for second well etc.
                rate += well_phase_flow_rate[well*num_phases + phase] * distr[phase];
            }
            return rate;
        }

        inline
        bool constraintBroken(const std::vector<double>& bhp,
                              const std::vector<double>& thp,
                              const std::vector<double>& well_phase_flow_rate,
                              const int well,
                              const int num_phases,
                              const WellType& well_type,
                              const WellControls* wc,
                              const int ctrl_index)
        {
            const WellControlType ctrl_type = well_controls_iget_type(wc, ctrl_index);
            const double target = well_controls_iget_target(wc, ctrl_index);
            const double* distr = well_controls_iget_distr(wc, ctrl_index);

            bool broken = false;

            switch (well_type) {
            case INJECTOR:
            {
                switch (ctrl_type) {
                case BHP:
                    broken = bhp[well] > target;
                    break;

                case THP:
                    broken = thp[well] > target;
                    break;

                case RESERVOIR_RATE: // Intentional fall-through
                case SURFACE_RATE:
                    broken = rateToCompare(well_phase_flow_rate,
                                           well, num_phases, distr) > target;
                    break;
                }
            }
            break;

            case PRODUCER:
            {
                switch (ctrl_type) {
                case BHP:
                    broken = bhp[well] < target;
                    break;

                case THP:
                    broken = thp[well] < target;
                    break;

                case RESERVOIR_RATE: // Intentional fall-through
                case SURFACE_RATE:
                    // Note that the rates compared below are negative,
                    // so breaking the constraints means: too high flow rate
                    // (as for injection).
                    broken = rateToCompare(well_phase_flow_rate,
                                           well, num_phases, distr) < target;
                    break;
                }
            }
            break;

            default:
                OPM_THROW(std::logic_error, "Can only handle INJECTOR and PRODUCER wells.");
            }

            return broken;
        }


        // ---------      Types      ---------

        /**
         * Simple hydrostatic correction for VFP table
         * @param wells - wells struct
         * @param w Well number
         * @param vfp_table VFP table
         * @param well_perforation_densities Densities at well perforations
         * @param gravity Gravitational constant (e.g., 9.81...)
         */
        inline
        double computeHydrostaticCorrection(const Wells& wells, const int w, double vfp_ref_depth,
                                            const double& rho, const double gravity) {
            if ( wells.well_connpos[w] == wells.well_connpos[w+1] )
            {
                // This is a well with no perforations.
                // If this is the last well we would subscript over the
                // bounds below.
                // we assume well_perforation_densities to be 0
                return 0;
            }
            const double well_ref_depth = wells.depth_ref[w];
            const double dh = vfp_ref_depth - well_ref_depth;
            const double dp = rho*gravity*dh;

            return dp;
        }

        inline
        double computeHydrostaticCorrection(const double well_ref_depth, const double vfp_ref_depth,
                                            const double rho, const double gravity) {
            const double dh = vfp_ref_depth - well_ref_depth;
            const double dp = rho * gravity * dh;

            return dp;
        }

        template <class Vector>
        inline
        Vector computeHydrostaticCorrection(const Wells& wells, const Vector vfp_ref_depth,
                const Vector& well_perforation_densities, const double gravity) {
            const int nw = wells.number_of_wells;
            Vector retval = Vector::Zero(nw);

#if HAVE_OPENMP
#pragma omp parallel for schedule(static)
#endif // HAVE_OPENMP
            for (int i=0; i<nw; ++i) {
                const int perf = wells.well_connpos[i];
                retval[i] = computeHydrostaticCorrection(wells, i, vfp_ref_depth[i], well_perforation_densities[perf], gravity);
            }

            return retval;
        }

        inline
        std::vector<double> computeHydrostaticCorrection(const Wells& wells, const std::vector<double>& vfp_ref_depth,
                                                         const std::vector<double>& well_perforation_densities, const double gravity) {
            const int nw = wells.number_of_wells;
            std::vector<double> retval(nw,0.0);

#if HAVE_OPENMP
#pragma omp parallel for schedule(static)
#endif // HAVE_OPENMP
            for (int i=0; i<nw; ++i) {
                const int perf = wells.well_connpos[i];
                retval[i] = computeHydrostaticCorrection(wells, i, vfp_ref_depth[i], well_perforation_densities[perf], gravity);
            }

            return retval;
        }

    } // namespace wellhelpers

}

#endif
