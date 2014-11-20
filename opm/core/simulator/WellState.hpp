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
#include <opm/core/well_controls.h>
#include <vector>
#include <cassert>

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
                temperature_.resize(nw, 273.15 + 20); // standard temperature for now
                wellrates_.resize(nw * np, 0.0);
                for (int w = 0; w < nw; ++w) {
                    assert((wells->type[w] == INJECTOR) || (wells->type[w] == PRODUCER));
                    const WellControls* ctrl = wells->ctrls[w];
                    if (well_controls_well_is_shut(ctrl)) {
                        // Shut well:
                        // 1. Assign zero well rates.
                        for (int p = 0; p < np; ++p) {
                            wellrates_[np*w + p] = 0.0;
                        }
                        // 2. Assign bhp equal to bhp control, if
                        //    applicable, otherwise assign equal to
                        //    first perforation cell pressure.
                        if (well_controls_get_current_type(ctrl) == BHP) {
                            bhp_[w] = well_controls_get_current_target( ctrl );
                        } else {
                            const int first_cell = wells->well_cells[wells->well_connpos[w]];
                            bhp_[w] = state.pressure()[first_cell];
                        }
                    } else {
                        // Open well:
                        // 1. Initialize well rates to match controls
                        //    if type is SURFACE_RATE.  Otherwise, we
                        //    cannot set the correct value here, so we
                        //    assign a small rate with the correct
                        //    sign so that any logic depending on that
                        //    sign will work as expected.
                        if (well_controls_get_current_type(ctrl) == SURFACE_RATE) {
                            const double rate_target = well_controls_get_current_target(ctrl);
                            const double * distr = well_controls_get_current_distr( ctrl );
                            for (int p = 0; p < np; ++p) {
                                wellrates_[np*w + p] = rate_target * distr[p];
                            }
                        } else {
                            const double small_rate = 1e-14;
                            const double sign = (wells->type[w] == INJECTOR) ? 1.0 : -1.0;
                            for (int p = 0; p < np; ++p) {
                                wellrates_[np*w + p] = small_rate * sign;
                            }
                        }
                        // 2. Initialize bhp to be target pressure if
                        //    bhp-controlled well, otherwise set to a
                        //    little above or below (depending on if
                        //    the well is an injector or producer)
                        //    pressure in first perforation cell.
                        if (well_controls_get_current_type(ctrl) == BHP) {
                            bhp_[w] = well_controls_get_current_target( ctrl );
                        } else {
                            const int first_cell = wells->well_cells[wells->well_connpos[w]];
                            const double safety_factor = (wells->type[w] == INJECTOR) ? 1.01 : 0.99;
                            bhp_[w] = safety_factor*state.pressure()[first_cell];
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

        /// One temperature per well.
        std::vector<double>& temperature() { return temperature_; }
        const std::vector<double>& temperature() const { return temperature_; }

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
        std::vector<double> temperature_;
        std::vector<double> wellrates_;
        std::vector<double> perfrates_;
        std::vector<double> perfpress_;
    };

} // namespace Opm

#endif // OPM_WELLSTATE_HEADER_INCLUDED
