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

#include <vector>

namespace Opm
{

class DeferredLogger;
class SummaryState;
template<typename Scalar, typename IndexTraits> class WellInterfaceGeneric;
enum class WellSegmentCompPressureDrop;
class WellSegments;
template<typename Scalar, typename IndexTraits> class WellState;

template<typename Scalar, typename IndexTraits>
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
    explicit MultisegmentWellGeneric(WellInterfaceGeneric<Scalar, IndexTraits>& baseif);

    // scale the segment rates and pressure based on well rates and bhp
    void scaleSegmentRatesWithWellRates(const std::vector<std::vector<int>>& segment_inlets,
                                        const std::vector<std::vector<int>>& segment_perforations,
                                        WellState<Scalar, IndexTraits>& well_state) const;
    void scaleSegmentPressuresWithBhp(WellState<Scalar, IndexTraits>& well_state) const;

    // components of the pressure drop to be included
    WellSegmentCompPressureDrop compPressureDrop() const;

    /// Detect oscillation or stagnation based on the residual measure history
    bool update_relaxation_factor(const std::vector<Scalar>& measure_history, Scalar& relaxation_factor, bool& regularize, DeferredLogger& deferred_logger) const;
    bool repeatedStagnation(const std::vector<Scalar>& measure_history, bool& regularize, DeferredLogger& deferred_logger) const;

    bool accelerationalPressureLossConsidered() const;
    bool frictionalPressureLossConsidered() const;

    Scalar getSegmentDp(const int seg,
                        const Scalar density,
                        const std::vector<Scalar>& seg_dp) const;

    const WellInterfaceGeneric<Scalar, IndexTraits>& baseif_;
};

}

#endif // OPM_MULTISEGMENTWELL_GENERIC_HEADER_INCLUDED
