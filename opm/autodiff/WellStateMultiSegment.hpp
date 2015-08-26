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

#ifndef OPM_WELLSTATEMULTISEGMENT_HEADER_INCLUDED
#define OPM_WELLSTATEMULTISEGMENT_HEADER_INCLUDED


#include <opm/core/wells.h>
#include <opm/core/well_controls.h>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/autodiff/WellMultiSegment.hpp>
#include <vector>
#include <cassert>
#include <string>
#include <utility>
#include <map>
#include <algorithm>
#include <array>

namespace Opm
{

    /// The state of a set of multi-sgemnet wells
    class WellStateMultiSegment
        : public WellState
    {
    public:

        // typedef std::array< int, 3 >  mapentry_t;
        // typedef std::map< std::string, mapentry_t > WellMapType;
        // this map needs to change a little bit?
        /* struct mapentry {
            int well_number;
            int number_of_segments;
            std::vector<int> number_of_performations;
        } */
        // MAYNOT NEED THIS

        /// Allocate and initialize if wells is non-null.  Also tries
        /// to give useful initial values to the bhp(), wellRates()
        /// and perfPhaseRates() fields, depending on controls
        template <class State, class PrevState>
        void init(const std::vector<WellMultiSegment>& wells, const State& state, const PrevState& prevState)
        {
            const int nw = wells.size();
            if (nw == 0) {
                return;
            }

            const int np = wells[0].numberOfPhases(); // number of the phases

            int nperf = 0; // the number of the perforations
            int nseg = 0; // the nubmer of the segments

            for (int iw = 0; iw < nw; ++iw) {
                nperf += wells[i].numberOfPerforations();
                nseg += wells[i].numberOfSegment();
            }

            bhp_.resize(nw);
            thp_.resize(nw);
            wellrates_.resize(nw * np, 0.0);

            current_controls_.resize(nw);
            for(int iw = 0; iw < nw; ++iw) {
                current_controls_[iw] = well_controls_get_current(wells[iw].wellControls());
            }

            for (int iw = 0; iw < nw; ++iw) {
                assert((wells[i].wellType() == INJECTOR) || (wells[i].wellType() == PRODUCER));
                const WellControls* ctrl = wells[iw]->wellControls();
            }


            // Map is used to map the value from the previous state to the current state as the initial values
            // TODO: handle this later.
            // Trying to figure out the work flow first.
        }

    private:
        // pressure for the segment nodes
        std::vector<double> seg_pressure_;
        // phase rates for the segments
        std::vector<double> seg_phaserates_;
        // phase rates for the completions
        std::vector<double> perf_phaserates_;
        // fractions for each segments (W, O, G)
        std::vector<double> seg_phasefrac_;
        // total flow rates for each segments, G_T
        std::vector<double> seg_totalrate_;
        std::vector<int> current_controls_;
        WellMapType wellMap_;
    };

} // namespace Opm


#endif // OPM_WELLSTATEMULTISEGMENT_HEADER_INCLUDE
