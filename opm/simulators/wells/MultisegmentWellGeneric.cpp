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
#include <opm/input/eclipse/Schedule/MSW/WellSegments.hpp>

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

#include <fmt/format.h>

namespace Opm
{

template<typename Scalar>
MultisegmentWellGeneric<Scalar>::
MultisegmentWellGeneric(WellInterfaceGeneric& baseif)
    : baseif_(baseif)
{
}

template<typename Scalar>
void
MultisegmentWellGeneric<Scalar>::
scaleSegmentRatesWithWellRates(const std::vector<std::vector<int>>& segment_inlets,
                               const std::vector<std::vector<int>>& segment_perforations,
                               WellState& well_state) const
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
            WellState::calculateSegmentRates(segment_inlets,
                                             segment_perforations,
                                             perforation_rates,
                                             num_single_phase, 0, rates);
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
int
MultisegmentWellGeneric<Scalar>::
segmentNumberToIndex(const int segment_number) const
{
    return segmentSet().segmentNumberToIndex(segment_number);
}

template<typename Scalar>
void
MultisegmentWellGeneric<Scalar>::
detectOscillations(const std::vector<double>& measure_history, bool& oscillate, bool& stagnate) const
{
    const auto it = measure_history.size() - 1;
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


template<class Scalar>
double
MultisegmentWellGeneric<Scalar>::getSegmentDp(const int seg,
                                              const double density,
                                              const std::vector<double>& seg_dp) const
{
    const double segment_depth = this->segmentSet()[seg].depth();
    const int outlet_segment_index = this->segmentNumberToIndex(this->segmentSet()[seg].outletSegment());
    const double segment_depth_outlet = seg == 0 ? baseif_.refDepth() : this->segmentSet()[outlet_segment_index].depth();
    double dp = wellhelpers::computeHydrostaticCorrection(segment_depth_outlet, segment_depth,
                                                          density, baseif_.gravity());
    // we add the hydrostatic correction from the outlet segment
    // in order to get the correction all the way to the bhp ref depth.
    if (seg > 0) {
        dp += seg_dp[outlet_segment_index];
    }

    return dp;
}

template class MultisegmentWellGeneric<double>;

} // namespace Opm
