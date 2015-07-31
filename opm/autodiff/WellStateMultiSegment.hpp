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

        /// Allocate and initialize if wells is non-null.  Also tries
        /// to give useful initial values to the bhp(), wellRates()
        /// and perfPhaseRates() fields, depending on controls
        template <class State, class PrevState>
        void init(const Wells* wells, const State& state, const PrevState& prevState)
        {
        }

    private:
        // pressure for the segment nodes
        // pressure for the top segment nodes are the bhp
        std::vector<double> segpressure_;
        // phase rates for the segments
        std::vector<double> segphaserates_;
        // phase rates for the completions
        std::vector<double> perfphaserates_;
        std::vector<int> current_controls_;
        // WellMapType wellMap_;
    };

} // namespace Opm


#endif // OPM_WELLSTATEMULTISEGMENT_HEADER_INCLUDE
