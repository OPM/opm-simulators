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

#ifndef OPM_WELLSTATEFULLYIMPLICITBLACKOIL_HEADER_INCLUDED
#define OPM_WELLSTATEFULLYIMPLICITBLACKOIL_HEADER_INCLUDED


#include <opm/core/wells.h>
#include <opm/core/well_controls.h>
#include <opm/core/simulator/WellState.hpp>
#include <vector>
#include <cassert>

namespace Opm
{

    /// The state of a set of wells, tailored for use by the fully
    /// implicit blackoil simulator.
    class WellStateFullyImplicitBlackoil
    {
    public:
        /// Allocate and initialize if wells is non-null.  Also tries
        /// to give useful initial values to the bhp(), wellRates()
        /// and perfPhaseRates() fields, depending on controls
        template <class State>
        void init(const Wells* wells, const State& state)
        {
            if (wells == 0) {
                return;
            }

            // We use the WellState::init() function to do bhp and well rates init.
            // The alternative would be to copy that function wholesale.
            basic_well_state_.init(wells, state);

            // Initialize perfphaserates_, which must be done here.
            const int nw = wells->number_of_wells;
            const int np = wells->number_of_phases;
            const int nperf = wells->well_connpos[nw];
            perfphaserates_.resize(nperf * np, 0.0);
            for (int w = 0; w < nw; ++w) {
                assert((wells->type[w] == INJECTOR) || (wells->type[w] == PRODUCER));
                const WellControls* ctrl = wells->ctrls[w];
                if (well_controls_well_is_shut(ctrl)) {
                    // Shut well: perfphaserates_ are all zero.
                } else {
                    // Open well: Initialize perfphaserates_ to well
                    // rates divided by the number of perforations.
                    const int num_perf_this_well = wells->well_connpos[w + 1] - wells->well_connpos[w];
                    for (int perf = wells->well_connpos[w]; perf < wells->well_connpos[w + 1]; ++perf) {
                        for (int p = 0; p < np; ++p) {
                            perfphaserates_[np*perf + p] = wellRates()[np*w + p] / double(num_perf_this_well);
                        }
                    }
                }
            }
        }

        /// One bhp pressure per well.
        std::vector<double>& bhp() { return basic_well_state_.bhp(); }
        const std::vector<double>& bhp() const { return basic_well_state_.bhp(); }

        /// One rate per well and phase.
        std::vector<double>& wellRates() { return basic_well_state_.wellRates(); }
        const std::vector<double>& wellRates() const { return basic_well_state_.wellRates(); }

        /// One rate per phase and well connection.
        std::vector<double>& perfPhaseRates() { return perfphaserates_; }
        const std::vector<double>& perfPhaseRates() const { return perfphaserates_; }

        /// For interfacing with functions that take a WellState.
        const WellState& basicWellState() const
        {
            return basic_well_state_;
        }

    private:
        WellState basic_well_state_;
        std::vector<double> perfphaserates_;
    };

} // namespace Opm


#endif // OPM_WELLSTATEFULLYIMPLICITBLACKOIL_HEADER_INCLUDED
