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

#include <config.h>
#include <opm/simulators/wells/MultisegmentWellGeneric.hpp>

#include <opm/common/utility/numeric/RootFinders.hpp>

#include <opm/input/eclipse/Schedule/VFPInjTable.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/VFPHelpers.hpp>
#include <opm/simulators/wells/VFPProperties.hpp>
#include <opm/simulators/wells/WellBhpThpCalculator.hpp>
#include <opm/simulators/wells/WellHelpers.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <cassert>
#include <cmath>
#include <stdexcept>

namespace Opm
{

template<typename Scalar>
MultisegmentWellGeneric<Scalar>::
MultisegmentWellGeneric(WellInterfaceGeneric& baseif)
    : baseif_(baseif)
    , segment_perforations_(numberOfSegments())
    , segment_inlets_(numberOfSegments())
    , segment_depth_diffs_(numberOfSegments(), 0.0)
    , perforation_segment_depth_diffs_(baseif_.numPerfs(), 0.0)
{
    // since we decide to use the WellSegments from the well parser. we can reuse a lot from it.
    // for other facilities needed but not available from parser, we need to process them here

    // initialize the segment_perforations_ and update perforation_segment_depth_diffs_
    const WellConnections& completion_set = baseif_.wellEcl().getConnections();
    // index of the perforation within wells struct
    // there might be some perforations not active, which causes the number of the perforations in
    // well_ecl_ and wells struct different
    // the current implementation is a temporary solution for now, it should be corrected from the parser
    // side
    int i_perf_wells = 0;
    baseif.perfDepth().resize(baseif_.numPerfs(), 0.);
    for (size_t perf = 0; perf < completion_set.size(); ++perf) {
        const Connection& connection = completion_set.get(perf);
        if (connection.state() == Connection::State::OPEN) {
            const int segment_index = segmentNumberToIndex(connection.segment());
            segment_perforations_[segment_index].push_back(i_perf_wells);
            baseif.perfDepth()[i_perf_wells] = connection.depth();
            const double segment_depth = segmentSet()[segment_index].depth();
            perforation_segment_depth_diffs_[i_perf_wells] = baseif.perfDepth()[i_perf_wells] - segment_depth;
            i_perf_wells++;
        }
    }

    // initialize the segment_inlets_
    for (int seg = 0; seg < numberOfSegments(); ++seg) {
        const Segment& segment = segmentSet()[seg];
        const int segment_number = segment.segmentNumber();
        const int outlet_segment_number = segment.outletSegment();
        if (outlet_segment_number > 0) {
            const int segment_index = segmentNumberToIndex(segment_number);
            const int outlet_segment_index = segmentNumberToIndex(outlet_segment_number);
            segment_inlets_[outlet_segment_index].push_back(segment_index);
        }
    }

    // calculating the depth difference between the segment and its oulet_segments
    // for the top segment, we will make its zero unless we find other purpose to use this value
    for (int seg = 1; seg < numberOfSegments(); ++seg) {
        const double segment_depth = segmentSet()[seg].depth();
        const int outlet_segment_number = segmentSet()[seg].outletSegment();
        const Segment& outlet_segment = segmentSet()[segmentNumberToIndex(outlet_segment_number)];
        const double outlet_depth = outlet_segment.depth();
        segment_depth_diffs_[seg] = segment_depth - outlet_depth;
    }
}

template<typename Scalar>
void
MultisegmentWellGeneric<Scalar>::
scaleSegmentRatesWithWellRates(WellState& well_state) const
{
    auto& ws = well_state.well(baseif_.indexOfWell());
    auto& segments = ws.segments;
    auto& segment_rates = segments.rates;
    for (int phase = 0; phase < baseif_.numPhases(); ++phase) {
        const double unscaled_top_seg_rate = segment_rates[phase];
        const double well_phase_rate = ws.surface_rates[phase];
        if (std::abs(unscaled_top_seg_rate) > 1e-12) {
            for (int seg = 0; seg < numberOfSegments(); ++seg) {
                segment_rates[baseif_.numPhases() * seg + phase] *= well_phase_rate / unscaled_top_seg_rate;
            }
        } else {
            // Due to various reasons, the well/top segment rate can be zero for this phase.
            // We can not scale this rate directly. The following approach is used to initialize the segment rates.
            double sumTw = 0;
            for (int perf = 0; perf < baseif_.numPerfs(); ++perf) {
                sumTw += baseif_.wellIndex()[perf];
            }

            // only handling this specific phase
            constexpr double num_single_phase = 1;
            std::vector<double> perforation_rates(num_single_phase * baseif_.numPerfs(), 0.0);
            const double perf_phaserate_scaled = ws.surface_rates[phase] / sumTw;
            for (int perf = 0; perf < baseif_.numPerfs(); ++perf) {
                perforation_rates[perf] = baseif_.wellIndex()[perf] * perf_phaserate_scaled;
            }

            std::vector<double> rates;
            WellState::calculateSegmentRates(segment_inlets_, segment_perforations_, perforation_rates, num_single_phase, 0, rates);
            for (int seg = 0; seg < numberOfSegments(); ++seg) {
                segment_rates[baseif_.numPhases() * seg + phase] = rates[seg];
            }
        }
    }
}

template <typename Scalar>
void
MultisegmentWellGeneric<Scalar>::
scaleSegmentPressuresWithBhp(WellState& well_state) const
{
    auto& ws = well_state.well(baseif_.indexOfWell());
    auto& segments = ws.segments;
    segments.scale_pressure(ws.bhp);
}

template<typename Scalar>
const WellSegments&
MultisegmentWellGeneric<Scalar>::
segmentSet() const
{
    return baseif_.wellEcl().getSegments();
}

template <typename Scalar>
int
MultisegmentWellGeneric<Scalar>::
numberOfSegments() const
{
    return segmentSet().size();
}

template <typename Scalar>
WellSegments::CompPressureDrop
MultisegmentWellGeneric<Scalar>::
compPressureDrop() const
{
    return segmentSet().compPressureDrop();
}

template<typename Scalar>
const std::vector<std::vector<int>>&
MultisegmentWellGeneric<Scalar>::
segmentInlets() const
{
    return segment_inlets_;
}

template<typename Scalar>
const std::vector<std::vector<int>>&
MultisegmentWellGeneric<Scalar>::
segmentPerforations() const
{
    return segment_perforations_;
}

template<typename Scalar>
int
MultisegmentWellGeneric<Scalar>::
segmentNumberToIndex(const int segment_number) const
{
    return segmentSet().segmentNumberToIndex(segment_number);
}

template<typename Scalar>
void
MultisegmentWellGeneric<Scalar>::
detectOscillations(const std::vector<double>& measure_history,
                   const int it,
                   bool& oscillate,
                   bool& stagnate) const
{
    if ( it < 2 ) {
        oscillate = false;
        stagnate = false;
        return;
    }

    stagnate = true;
    const double F0 = measure_history[it];
    const double F1 = measure_history[it - 1];
    const double F2 = measure_history[it - 2];
    const double d1 = std::abs((F0 - F2) / F0);
    const double d2 = std::abs((F0 - F1) / F0);

    const double oscillaton_rel_tol = 0.2;
    oscillate = (d1 < oscillaton_rel_tol) && (oscillaton_rel_tol < d2);

    const double stagnation_rel_tol = 1.e-2;
    stagnate = std::abs((F1 - F2) / F2) <= stagnation_rel_tol;
}

template<typename Scalar>
bool
MultisegmentWellGeneric<Scalar>::
frictionalPressureLossConsidered() const
{
    // HF- and HFA needs to consider frictional pressure loss
    return (segmentSet().compPressureDrop() != WellSegments::CompPressureDrop::H__);
}

template<typename Scalar>
bool
MultisegmentWellGeneric<Scalar>::
accelerationalPressureLossConsidered() const
{
    return (segmentSet().compPressureDrop() == WellSegments::CompPressureDrop::HFA);
}

template class MultisegmentWellGeneric<double>;

} // namespace Opm
