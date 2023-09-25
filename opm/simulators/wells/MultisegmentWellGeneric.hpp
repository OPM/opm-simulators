/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
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


#ifndef OPM_MULTISEGMENTWELL_GENERIC_HEADER_INCLUDED
#define OPM_MULTISEGMENTWELL_GENERIC_HEADER_INCLUDED

#include <functional>
#include <optional>
#include <vector>
#include <array>

namespace Opm
{

class DeferredLogger;
class SummaryState;
class WellInterfaceGeneric;
enum class WellSegmentCompPressureDrop;
class WellSegments;
class WellState;

template <typename Scalar>
class MultisegmentWellGeneric
{
public:
    // get the WellSegments from the well_ecl_
    const WellSegments& segmentSet() const;

    // segment number is an ID of the segment, it is specified in the deck
    // get the loation of the segment with a segment number in the segmentSet
    int segmentNumberToIndex(const int segment_number) const;

    /// number of segments for this well
    int numberOfSegments() const;

protected:
    MultisegmentWellGeneric(WellInterfaceGeneric& baseif);

    // scale the segment rates and pressure based on well rates and bhp
    void scaleSegmentRatesWithWellRates(const std::vector<std::vector<int>>& segment_inlets,
                                        const std::vector<std::vector<int>>& segment_perforations,
                                        WellState& well_state) const;
    void scaleSegmentPressuresWithBhp(WellState& well_state) const;

    // components of the pressure drop to be included
    WellSegmentCompPressureDrop compPressureDrop() const;

    /// Detect oscillation or stagnation based on the residual measure history
    void detectOscillations(const std::vector<double>& measure_history,
                            bool& oscillate,
                            bool& stagnate) const;

    bool accelerationalPressureLossConsidered() const;
    bool frictionalPressureLossConsidered() const;

    double getSegmentDp(const int seg,
                        const double density,
                        const std::vector<double>& seg_dp) const;

    const WellInterfaceGeneric& baseif_;
};

}

#endif // OPM_MULTISEGMENTWELL_GENERIC_HEADER_INCLUDED
