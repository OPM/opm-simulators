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

#include <opm/input/eclipse/Schedule/MSW/WellSegments.hpp>

#include <functional>
#include <optional>
#include <vector>
#include <array>

namespace Opm
{

class DeferredLogger;
class SummaryState;
class WellInterfaceGeneric;
class WellState;

template <typename Scalar>
class MultisegmentWellGeneric
{
public:
    //! \brief Returns the inlet segments for each segment.
    const std::vector<std::vector<int>>& segmentInlets() const;

    //! \brief Returns the perforation index for each segment.
    const std::vector<std::vector<int>>& segmentPerforations() const;

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
    void scaleSegmentRatesWithWellRates(WellState& well_state) const;
    void scaleSegmentPressuresWithBhp(WellState& well_state) const;

    // components of the pressure drop to be included
    WellSegments::CompPressureDrop compPressureDrop() const;

    /// Detect oscillation or stagnation based on the residual measure history
    void detectOscillations(const std::vector<double>& measure_history,
                            const int it,
                            bool& oscillate,
                            bool& stagnate) const;

    bool accelerationalPressureLossConsidered() const;
    bool frictionalPressureLossConsidered() const;

    const WellInterfaceGeneric& baseif_;

    // TODO: trying to use the information from the Well opm-parser as much
    // as possible, it will possibly be re-implemented later for efficiency reason.

    // the completions that is related to each segment
    // the completions's ids are their index in the vector well_index_, well_cell_
    // This is also assuming the order of the completions in Well is the same with
    // the order of the completions in wells.
    // it is for convinience reason. we can just calcuate the inforation for segment once then using it for all the perofrations
    // belonging to this segment
    std::vector<std::vector<int>> segment_perforations_;

    // the inlet segments for each segment. It is for convenience and efficiency reason
    std::vector<std::vector<int>> segment_inlets_;

    std::vector<double> segment_depth_diffs_;

    // depth difference between the segment and the perforation
    // or in another way, the depth difference between the perforation and
    // the segment the perforation belongs to
    std::vector<double> perforation_segment_depth_diffs_;
};

}

#endif // OPM_MULTISEGMENTWELL_GENERIC_HEADER_INCLUDED
