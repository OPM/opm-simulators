/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.

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

#include <opm/autodiff/WellDensitySegmented.hpp>
#include <opm/core/wells.h>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <numeric>
#include <cmath>


std::vector<double> Opm::WellDensitySegmented::computeConnectionPressureDelta(const Wells& wells,
                                                                              const WellState& wstate,
                                                                              const PhaseUsage& phase_usage,
                                                                              const std::vector<double>& b_perf,
                                                                              const std::vector<double>& rsmax_perf,
                                                                              const std::vector<double>& rvmax_perf,
                                                                              const std::vector<double>& z_perf,
                                                                              const std::vector<double>& surf_dens,
                                                                              const double gravity)
{
    // Verify that we have consistent input.
    const int np = wells.number_of_phases;
    const int nw = wells.number_of_wells;
    const int nperf = wells.well_connpos[nw];
    if (wells.number_of_phases != phase_usage.num_phases) {
        OPM_THROW(std::logic_error, "Inconsistent input: wells vs. phase_usage.");
    }
    if (surf_dens.size() != size_t(wells.number_of_phases)) {
        OPM_THROW(std::logic_error, "Inconsistent input: surf_dens vs. phase_usage.");
    }
    if (nperf*np != int(wstate.perfRates().size())) {
        OPM_THROW(std::logic_error, "Inconsistent input: wells vs. wstate.");
    }
    if (nperf*np != int(b_perf.size())) {
        OPM_THROW(std::logic_error, "Inconsistent input: wells vs. b_perf.");
    }
    if (nperf != int(z_perf.size())) {
        OPM_THROW(std::logic_error, "Inconsistent input: wells vs. z_perf.");
    }
    if ((!rsmax_perf.empty()) || (!rvmax_perf.empty())) {
        // Need both oil and gas phases.
        if (!phase_usage.phase_used[BlackoilPhases::Liquid]) {
            OPM_THROW(std::logic_error, "Oil phase inactive, but non-empty rsmax_perf or rvmax_perf.");
        }
        if (!phase_usage.phase_used[BlackoilPhases::Vapour]) {
            OPM_THROW(std::logic_error, "Gas phase inactive, but non-empty rsmax_perf or rvmax_perf.");
        }
    }

    // Algorithm:

    // We'll assume the perforations are given in order from top to
    // bottom for each well.  By top and bottom we do not necessarily
    // mean in a geometric sense (depth), but in a topological sense:
    // the 'top' perforation is nearest to the surface topologically.
    // Our goal is to compute a pressure delta for each perforation.

    // 1. Compute the flow (in surface volume units for each
    //    component) exiting up the wellbore from each perforation,
    //    taking into account flow from lower in the well, and
    //    in/out-flow at each perforation.
    std::vector<double> q_out_perf(nperf*np);
    for (int w = 0; w < nw; ++w) {
        // Iterate over well perforations from bottom to top.
        for (int perf = wells.well_connpos[w+1] - 1; perf >= wells.well_connpos[w]; --perf) {
            for (int phase = 0; phase < np; ++phase) {
                if (perf == wells.well_connpos[w+1] - 1) {
                    // This is the bottom perforation. No flow from below.
                    q_out_perf[perf*np + phase] = 0.0;
                } else {
                    // Set equal to flow from below.
                    q_out_perf[perf*np + phase] = q_out_perf[(perf+1)*np + phase];
                }
                // Subtract outflow through perforation.
                q_out_perf[perf*np + phase] -= wstate.perfRates()[perf*np + phase];
            }
        }
    }

    // 2. Compute the component mix at each perforation as the
    //    absolute values of the surface rates divided by their sum.
    //    Then compute volume ratios (formation factors) for each perforation.
    //    Finally compute densities for the segments associated with each perforation.
    const int gaspos = phase_usage.phase_pos[BlackoilPhases::Vapour];
    const int oilpos = phase_usage.phase_pos[BlackoilPhases::Liquid];
    std::vector<double> mix(np);
    std::vector<double> x(np);
    std::vector<double> dens(nperf);
    for (int w = 0; w < nw; ++w) {
        for (int perf = wells.well_connpos[w]; perf < wells.well_connpos[w+1]; ++perf) {
            // Find component mix.
            const double tot_surf_rate = std::accumulate(q_out_perf.begin() + np*perf,
                                                         q_out_perf.begin() + np*(perf+1), 0.0);
            if (tot_surf_rate != 0.0) {
                for (int phase = 0; phase < np; ++phase) {
                    mix[phase] = std::fabs(q_out_perf[perf*np + phase]/tot_surf_rate);
                }
            } else {
                // No flow => use well specified fractions for mix.
                std::copy(wells.comp_frac + w*np, wells.comp_frac + (w+1)*np, mix.begin());
            }
            // Compute volume ratio.
            x = mix;
            double rs = 0.0;
            double rv = 0.0;
            if (!rsmax_perf.empty() && mix[oilpos] > 0.0) {
                rs = std::min(mix[gaspos]/mix[oilpos], rsmax_perf[perf]);
            }
            if (!rvmax_perf.empty() && mix[gaspos] > 0.0) {
                rv = std::min(mix[oilpos]/mix[gaspos], rvmax_perf[perf]);
            }
            if (rs != 0.0) {
                // Subtract gas in oil from gas mixture
                x[gaspos] = (mix[gaspos] - mix[oilpos]*rs)/(1.0 - rs*rv);
            }
            if (rv != 0.0) {
                // Subtract oil in gas from oil mixture
                x[oilpos] = (mix[oilpos] - mix[gaspos]*rv)/(1.0 - rs*rv);;
            }
            double volrat = 0.0;
            for (int phase = 0; phase < np; ++phase) {
                volrat += x[phase] / b_perf[perf*np + phase];
            }
            // Compute segment density.
            dens[perf] = std::inner_product(surf_dens.begin(), surf_dens.end(), mix.begin(), 0.0) / volrat;
        }
    }

    // 3. Compute pressure differences between perforations.
    //    dp_perf will contain the pressure difference between a
    //    perforation and the one above it, except for the first
    //    perforation for each well, for which it will be the
    //    difference to the reference (bhp) depth.
    std::vector<double> dp_perf(nperf);
    for (int w = 0; w < nw; ++w) {
        for (int perf = wells.well_connpos[w]; perf < wells.well_connpos[w+1]; ++perf) {
            const double z_above = perf == wells.well_connpos[w] ? wells.depth_ref[w] : z_perf[perf - 1];
            const double dz = z_perf[perf] - z_above;
            dp_perf[perf] = dz * dens[perf] * gravity;
        }
    }

    // 4. Compute pressure differences to the reference point (bhp) by
    //    accumulating the already computed adjacent pressure
    //    differences, storing the result in dp_perf.
    //    This accumulation must be done per well.
    for (int w = 0; w < nw; ++w) {
        const auto beg = dp_perf.begin() + wells.well_connpos[w];
        const auto end = dp_perf.begin() + wells.well_connpos[w + 1];
        std::partial_sum(beg, end, beg);
    }

    return dp_perf;
}
