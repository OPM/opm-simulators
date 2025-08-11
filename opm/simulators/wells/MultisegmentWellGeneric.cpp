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

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/models/blackoil/blackoilvariableandequationindices.hh>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/VFPHelpers.hpp>
#include <opm/simulators/wells/VFPProperties.hpp>
#include <opm/simulators/wells/WellBhpThpCalculator.hpp>
#include <opm/simulators/wells/WellHelpers.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <cmath>

#include <fmt/format.h>

namespace Opm {

template<typename FluidSystem, typename Indices>
MultisegmentWellGeneric<FluidSystem, Indices>::
MultisegmentWellGeneric(WellInterfaceGeneric<FluidSystem, Indices>& baseif)
    : baseif_(baseif)
{
}

template<typename FluidSystem, typename Indices>
void
MultisegmentWellGeneric<FluidSystem, Indices>::
scaleSegmentRatesWithWellRates(const std::vector<std::vector<int>>& segment_inlets,
                               const std::vector<std::vector<int>>& segment_perforations,
                               WellState<FluidSystem, Indices>& well_state) const
{
    auto& ws = well_state.well(baseif_.indexOfWell());
    auto& segments = ws.segments;
    auto& segment_rates = segments.rates;
    Scalar sumTw = 0;
    bool calculateSumTw = std::any_of(segment_rates.begin(), segment_rates.end(), [](const auto& unscaled_top_seg_rate) {
        return std::abs(unscaled_top_seg_rate) <= 1e-12;
    });
    if (calculateSumTw) {
        // Due to various reasons, the well/top segment rate can be zero for this phase.
        // We can not scale this rate directly. The following approach is used to initialize the segment rates.
        for (int perf = 0; perf < baseif_.numLocalPerfs(); ++perf) {
            sumTw += baseif_.wellIndex()[perf];
        }
        // We need to communicate here to scale the perf_phaserate_scaled with the contribution of all segments
        sumTw = ws.parallel_info.get().communication().sum(sumTw);
    }
    for (int phase = 0; phase < baseif_.numPhases(); ++phase) {
        const Scalar unscaled_top_seg_rate = segment_rates[phase];
        const Scalar well_phase_rate = ws.surface_rates[phase];
        if (std::abs(unscaled_top_seg_rate) > 1e-12) {
            for (int seg = 0; seg < numberOfSegments(); ++seg) {
                segment_rates[baseif_.numPhases() * seg + phase] *= well_phase_rate / unscaled_top_seg_rate;
            }
        } else { // In this case, we calculated sumTw
            // only handling this specific phase
            constexpr Scalar num_single_phase = 1;
            std::vector<Scalar> perforation_rates(num_single_phase * baseif_.numLocalPerfs(), 0.0);
            const Scalar perf_phaserate_scaled = ws.surface_rates[phase] / sumTw;
            for (int perf = 0; perf < baseif_.numLocalPerfs(); ++perf) {
                perforation_rates[perf] = baseif_.wellIndex()[perf] * perf_phaserate_scaled;
            }

            std::vector<Scalar> rates;
            WellState<FluidSystem, Indices>::calculateSegmentRates(ws.parallel_info,
                                                     segment_inlets,
                                                     segment_perforations,
                                                     perforation_rates,
                                                     num_single_phase, 0, rates);
            for (int seg = 0; seg < numberOfSegments(); ++seg) {
                segment_rates[baseif_.numPhases() * seg + phase] = rates[seg];
            }
        }
    }
}

template <typename FluidSystem, typename Indices>
void
MultisegmentWellGeneric<FluidSystem, Indices>::
scaleSegmentPressuresWithBhp(WellState<FluidSystem, Indices>& well_state) const
{
    auto& ws = well_state.well(baseif_.indexOfWell());
    auto& segments = ws.segments;
    segments.scale_pressure(ws.bhp);
}

template<typename FluidSystem, typename Indices>
const WellSegments&
MultisegmentWellGeneric<FluidSystem, Indices>::
segmentSet() const
{
    return baseif_.wellEcl().getSegments();
}

template <typename FluidSystem, typename Indices>
int
MultisegmentWellGeneric<FluidSystem, Indices>::
numberOfSegments() const
{
    return segmentSet().size();
}

template <typename FluidSystem, typename Indices>
WellSegments::CompPressureDrop
MultisegmentWellGeneric<FluidSystem, Indices>::
compPressureDrop() const
{
    return segmentSet().compPressureDrop();
}

template<typename FluidSystem, typename Indices>
int
MultisegmentWellGeneric<FluidSystem, Indices>::
segmentNumberToIndex(const int segment_number) const
{
    return segmentSet().segmentNumberToIndex(segment_number);
}

template<typename FluidSystem, typename Indices>
bool
MultisegmentWellGeneric<FluidSystem, Indices>::
update_relaxation_factor(const std::vector<Scalar>& measure_history, Scalar& relaxation_factor, bool& regularize, DeferredLogger& deferred_logger) const
{
    const auto it = measure_history.size() - 1;
    if ( it < 2 ) {
        return false;
    }
    const Scalar F0 = measure_history[it];
    const Scalar F1 = measure_history[it - 1];
    const Scalar F2 = measure_history[it - 2];
    const Scalar d1 = std::abs((F0 - F2) / F0);
    const Scalar d2 = std::abs((F0 - F1) / F0);

    const Scalar oscillaton_rel_tol = 0.2;
    bool oscillate = (d1 < oscillaton_rel_tol) && (oscillaton_rel_tol < d2);
    if (!oscillate)
        return false;
   
    const Scalar min_relaxation_factor = 0.6;
    std::ostringstream sstr;
    if (relaxation_factor == min_relaxation_factor) {
        sstr << " well " << baseif_.name() << " observes severe oscillation. Terminates after " << it << "iterations.\n";

    }
    const Scalar reduction_mutliplier = 0.9;
    relaxation_factor = std::max(relaxation_factor * reduction_mutliplier, min_relaxation_factor);

    if (relaxation_factor == min_relaxation_factor) {
        return true;
    }
    sstr << " well " << baseif_.name() << " observes oscillation in inner iteration " << it << "\n";
    sstr << " relaxation_factor is " << relaxation_factor << "\n";
    regularize = true;
    deferred_logger.debug(sstr.str());
    return false;
}

template<typename FluidSystem, typename Indices>
bool
MultisegmentWellGeneric<FluidSystem, Indices>::
repeatedStagnation(const std::vector<Scalar>& measure_history, bool& regularize, DeferredLogger& deferred_logger) const
{
    const auto it = measure_history.size() - 1;
    if ( it < 2 ) {
        return false;
    }
    const Scalar F0 = measure_history[it];
    const Scalar F1 = measure_history[it - 1];
    const Scalar F2 = measure_history[it - 2];
    const Scalar stagnation_rel_tol = 1.e-2;
    bool stagnate = std::abs((F0 - F1) / F1) <= stagnation_rel_tol && std::abs((F1 - F2) / F2) <= stagnation_rel_tol;
    if (!stagnate)
        return false;

    std::ostringstream sstr;
    sstr << " well " << baseif_.name() << " observes stagnation in inner iteration " << it << "\n";
    // we frist try to regularize and only stop iterating if it still stagnates
    if (regularize ) {
        sstr << "Inner iterations are terminated.\n";
        return true;
    }
    regularize = true;
    sstr << "Try to regularize the equation.\n";
    deferred_logger.debug(sstr.str());
    return false;
}


template<typename FluidSystem, typename Indices>
bool
MultisegmentWellGeneric<FluidSystem, Indices>::
frictionalPressureLossConsidered() const
{
    // HF- and HFA needs to consider frictional pressure loss
    return (segmentSet().compPressureDrop() != WellSegments::CompPressureDrop::H__);
}

template<typename FluidSystem, typename Indices>
bool
MultisegmentWellGeneric<FluidSystem, Indices>::
accelerationalPressureLossConsidered() const
{
    return (segmentSet().compPressureDrop() == WellSegments::CompPressureDrop::HFA);
}


template<typename FluidSystem, typename Indices>
typename MultisegmentWellGeneric<FluidSystem, Indices>::Scalar
MultisegmentWellGeneric<FluidSystem, Indices>::getSegmentDp(const int seg,
                                              const Scalar density,
                                              const std::vector<Scalar>& seg_dp) const
{
    const Scalar segment_depth = this->segmentSet()[seg].depth();
    const int outlet_segment_index = this->segmentNumberToIndex(this->segmentSet()[seg].outletSegment());
    const Scalar segment_depth_outlet = seg == 0 ? baseif_.refDepth() : this->segmentSet()[outlet_segment_index].depth();
    Scalar dp = wellhelpers::computeHydrostaticCorrection(segment_depth_outlet, segment_depth,
                                                          density, baseif_.gravity());
    // we add the hydrostatic correction from the outlet segment
    // in order to get the correction all the way to the bhp ref depth.
    if (seg > 0) {
        dp += seg_dp[outlet_segment_index];
    }

    return dp;
}

    template<class Scalar>
    using FS = BlackOilFluidSystem<Scalar, BlackOilDefaultFluidSystemIndices>;

#define INSTANTIATE(T,...) \
    template class MultisegmentWellGeneric<FS<T>, __VA_ARGS__>;

#define INSTANTIATE_TYPE(T)                                                  \
    INSTANTIATE(T,BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>) \
    INSTANTIATE(T,BlackOilOnePhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>) \
    INSTANTIATE(T,BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,5u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,0u,0u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,2u,0u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,false,0u,2u,0u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,2u,0u,false,false,0u,2u,0u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,2u,0u>)  \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,0u,0u>)  \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,0u,0u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,true,0u,0u,0u>)  \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<1u,0u,0u,0u,false,false,0u,0u,0u>) \
    INSTANTIATE(T,BlackOilVariableAndEquationIndices<0u,0u,0u,0u,false,false,0u,0u>)            \
    INSTANTIATE(T,BlackOilVariableAndEquationIndices<0u,0u,0u,0u,false,false,1u,0u>)            \
    INSTANTIATE(T,BlackOilVariableAndEquationIndices<0u,0u,0u,0u,true,false,0u,0u>)             \
    INSTANTIATE(T,BlackOilVariableAndEquationIndices<0u,0u,0u,0u,false,true,0u,0u>)             \
    INSTANTIATE(T,BlackOilVariableAndEquationIndices<0u,0u,0u,0u,false,true,2u,0u>)             \
    INSTANTIATE(T,BlackOilVariableAndEquationIndices<1u,0u,0u,0u,false,false,0u,0u>)            \
    INSTANTIATE(T,BlackOilVariableAndEquationIndices<0u,1u,0u,0u,false,false,0u,0u>)            \
    INSTANTIATE(T,BlackOilVariableAndEquationIndices<0u,0u,1u,0u,false,false,0u,0u>)            \
    INSTANTIATE(T,BlackOilVariableAndEquationIndices<0u,0u,0u,1u,false,false,0u,0u>)            \
    INSTANTIATE(T,BlackOilVariableAndEquationIndices<0u,0u,0u,1u,false,true,0u,0u>)             \
    INSTANTIATE(T,BlackOilVariableAndEquationIndices<1u,0u,0u,0u,true,false,0u,0u>)

    INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
    INSTANTIATE_TYPE(float)
#endif

} // namespace Opm
