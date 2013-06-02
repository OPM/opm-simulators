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

#ifndef OPM_WELLSTATE_HEADER_INCLUDED
#define OPM_WELLSTATE_HEADER_INCLUDED

#include <opm/core/wells.h>
#include <vector>

namespace Opm
{

    /// The state of a set of wells.
    class WellState
    {
    public:
        /// Allocate and initialize if wells is non-null.
        /// Also tries to give useful initial values to the bhp() and
        /// wellRates() fields, depending on controls.  The
        /// perfRates() field is filled with zero, and perfPress()
        /// with -1e100.
        template <class State>
        void init(const Wells* wells, const State& state)
        {
            if (wells) {
                const int nw = wells->number_of_wells;
                const int np = wells->number_of_phases;
                bhp_.resize(nw);
                wellrates_.resize(nw * np, 0.0);
                for (int w = 0; w < nw; ++w) {
                    const WellControls* ctrl = wells->ctrls[w];
                    // Initialize bhp to be target pressure if
                    // bhp-controlled well, otherwise set to a little
                    // above or below (depending on if the well is an
                    // injector or producer) pressure in first perforation
                    // cell.
                    if ((ctrl->current < 0) || // SHUT
                        (ctrl->type[ctrl->current] != BHP)) {
                        const int first_cell = wells->well_cells[wells->well_connpos[w]];
                        const double safety_factor = (wells->type[w] == INJECTOR) ? 1.01 : 0.99;
                        bhp_[w] = safety_factor*state.pressure()[first_cell];
                    } else {
                        bhp_[w] = ctrl->target[ctrl->current];
                    }
                    // Initialize well rates to match controls if type is SURFACE_RATE
                    if ((ctrl->current >= 0) && // open well
                        (ctrl->type[ctrl->current] != SURFACE_RATE)) {
                        const double rate_target = ctrl->target[ctrl->current];
                        for (int p = 0; p < np; ++p) {
                            const double phase_distr = ctrl->distr[np * ctrl->current + p]
                                * wells->comp_frac[np * w + p];
                            wellrates_[np*w + p] = rate_target * phase_distr;
                        }
                    }
                }
                // The perforation rates and perforation pressures are
                // not expected to be consistent with bhp_ and wellrates_
                // after init().
                perfrates_.resize(wells->well_connpos[nw], 0.0);
                perfpress_.resize(wells->well_connpos[nw], -1e100);
            }
        }

        /// One bhp pressure per well.
        std::vector<double>& bhp() { return bhp_; }
        const std::vector<double>& bhp() const { return bhp_; }

        /// One rate per well and phase.
        std::vector<double>& wellRates() { return wellrates_; }
        const std::vector<double>& wellRates() const { return wellrates_; }

        /// One rate per well connection.
        std::vector<double>& perfRates() { return perfrates_; }
        const std::vector<double>& perfRates() const { return perfrates_; }

        /// One pressure per well connection.
        std::vector<double>& perfPress() { return perfpress_; }
        const std::vector<double>& perfPress() const { return perfpress_; }

    private:
        std::vector<double> bhp_;
        std::vector<double> wellrates_;
        std::vector<double> perfrates_;
        std::vector<double> perfpress_;
    };

} // namespace Opm

#endif // OPM_WELLSTATE_HEADER_INCLUDED
