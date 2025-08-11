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

#include <opm/simulators/wells/MultisegmentWellEval.hpp>

#include <opm/input/eclipse/Schedule/MSW/WellSegments.hpp>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/models/blackoil/blackoilvariableandequationindices.hh>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/simulators/timestepping/ConvergenceReport.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/MultisegmentWellAssemble.hpp>
#include <opm/simulators/wells/WellAssemble.hpp>
#include <opm/simulators/wells/WellConvergence.hpp>
#include <opm/simulators/wells/WellInterfaceIndices.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <numeric>
#include <utility>
#include <vector>

namespace Opm
{

template<typename FluidSystem, typename Indices>
MultisegmentWellEval<FluidSystem,Indices>::
MultisegmentWellEval(WellInterfaceIndices<FluidSystem,Indices>& baseif, const ParallelWellInfo<Scalar>& pw_info)
    :  MultisegmentWellGeneric<FluidSystem, Indices>(baseif)
    , pw_info_(pw_info)
    , baseif_(baseif)
    , linSys_(*this, pw_info)
    , primary_variables_(baseif)
    , segments_(this->numberOfSegments(), pw_info, baseif)
    , cell_perforation_depth_diffs_(baseif_.numLocalPerfs(), 0.0)
    , cell_perforation_pressure_diffs_(baseif_.numLocalPerfs(), 0.0)
{
}

template<typename FluidSystem, typename Indices>
void
MultisegmentWellEval<FluidSystem,Indices>::
initMatrixAndVectors()
{
    linSys_.init(baseif_.numLocalPerfs(),
                 baseif_.cells(), segments_.inlets(),
                 segments_.perforations());
    primary_variables_.resize(this->numberOfSegments());
}

template<typename FluidSystem, typename Indices>
ConvergenceReport
MultisegmentWellEval<FluidSystem,Indices>::
getWellConvergence(const WellState<FluidSystem, Indices>& well_state,
                   const std::vector<Scalar>& B_avg,
                   DeferredLogger& deferred_logger,
                   const Scalar max_residual_allowed,
                   const Scalar tolerance_wells,
                   const Scalar relaxed_inner_tolerance_flow_ms_well,
                   const Scalar tolerance_pressure_ms_wells,
                   const Scalar relaxed_inner_tolerance_pressure_ms_well,
                   const bool relax_tolerance, 
                   const bool well_is_stopped) const
{
    assert(int(B_avg.size()) == baseif_.numComponents());

    // checking if any residual is NaN or too large. The two large one is only handled for the well flux
    std::vector<std::vector<Scalar>> abs_residual(this->numberOfSegments(),
                                                  std::vector<Scalar>(numWellEq, 0.0));
    for (int seg = 0; seg < this->numberOfSegments(); ++seg) {
        for (int eq_idx = 0; eq_idx < numWellEq; ++eq_idx) {
            abs_residual[seg][eq_idx] = std::abs(linSys_.residual()[seg][eq_idx]);
        }
    }

    std::vector<Scalar> maximum_residual(numWellEq, 0.0);

    ConvergenceReport report;
    // TODO: the following is a little complicated, maybe can be simplified in some way?
    for (int eq_idx = 0; eq_idx < numWellEq; ++eq_idx) {
        for (int seg = 0; seg < this->numberOfSegments(); ++seg) {
            if (eq_idx < baseif_.numComponents()) { // phase or component mass equations
                const Scalar flux_residual = B_avg[eq_idx] * abs_residual[seg][eq_idx];
                if (flux_residual > maximum_residual[eq_idx]) {
                    maximum_residual[eq_idx] = flux_residual;
                }
            } else { // pressure or control equation
                // for the top segment (seg == 0), it is control equation, will be checked later separately
                if (seg > 0) {
                    // Pressure equation
                    const Scalar pressure_residual = abs_residual[seg][eq_idx];
                    if (pressure_residual > maximum_residual[eq_idx]) {
                        maximum_residual[eq_idx] = pressure_residual;
                    }
                }
            }
        }
    }

    using CR = ConvergenceReport;
    for (int eq_idx = 0; eq_idx < numWellEq; ++eq_idx) {
        if (eq_idx < baseif_.numComponents()) { // phase or component mass equations
            const Scalar flux_residual = maximum_residual[eq_idx];
            // TODO: the report can not handle the segment number yet.

            if (std::isnan(flux_residual)) {
                report.setWellFailed({CR::WellFailure::Type::MassBalance, CR::Severity::NotANumber, eq_idx, baseif_.name()});
            } else if (flux_residual > max_residual_allowed) {
                report.setWellFailed({CR::WellFailure::Type::MassBalance, CR::Severity::TooLarge, eq_idx, baseif_.name()});
            } else if (!relax_tolerance && flux_residual > tolerance_wells) {
                report.setWellFailed({CR::WellFailure::Type::MassBalance, CR::Severity::Normal, eq_idx, baseif_.name()});
            } else if (flux_residual > relaxed_inner_tolerance_flow_ms_well) {
                report.setWellFailed({CR::WellFailure::Type::MassBalance, CR::Severity::Normal, eq_idx, baseif_.name()});
            }
        } else { // pressure equation
            const Scalar pressure_residual = maximum_residual[eq_idx];
            const int dummy_component = -1;
            if (std::isnan(pressure_residual)) {
                report.setWellFailed({CR::WellFailure::Type::Pressure, CR::Severity::NotANumber, dummy_component, baseif_.name()});
            } else if (std::isinf(pressure_residual)) {
                report.setWellFailed({CR::WellFailure::Type::Pressure, CR::Severity::TooLarge, dummy_component, baseif_.name()});
            } else if (!relax_tolerance && pressure_residual > tolerance_pressure_ms_wells) {
                report.setWellFailed({CR::WellFailure::Type::Pressure, CR::Severity::Normal, dummy_component, baseif_.name()});
            } else if (pressure_residual > relaxed_inner_tolerance_pressure_ms_well) {
                report.setWellFailed({CR::WellFailure::Type::Pressure, CR::Severity::Normal, dummy_component, baseif_.name()});
            }
        }
    }

    WellConvergence(baseif_).
        checkConvergenceControlEq(well_state,
                                  {tolerance_pressure_ms_wells,
                                   tolerance_pressure_ms_wells,
                                   tolerance_wells,
                                   tolerance_wells,
                                   max_residual_allowed},
                                  std::abs(linSys_.residual()[0][SPres]),
                                  well_is_stopped,  
                                  report,
                                  deferred_logger);

    // for stopped well, we do not enforce the following checking to avoid dealing with sign of near-zero values
    // for BHP or THP controlled wells, we need to make sure the flow direction is correct
    if (!well_is_stopped && baseif_.isPressureControlled(well_state)) {
        // checking the flow direction
        const Scalar sign = baseif_.isProducer() ? -1. : 1.;
        const auto weight_total_flux = this->primary_variables_.getWQTotal() * sign;
        constexpr int dummy_phase = -1;
        if (weight_total_flux < 0.) {
            report.setWellFailed(
                    {CR::WellFailure::Type::WrongFlowDirection, CR::Severity::Normal, dummy_phase, baseif_.name()});
        }
    }

    return report;
}

template<typename FluidSystem, typename Indices>
typename MultisegmentWellEval<FluidSystem,Indices>::EvalWell
MultisegmentWellEval<FluidSystem,Indices>::
extendEval(const Eval& in) const
{
    EvalWell out = 0.0;
    out.setValue(in.value());
    for(int eq_idx = 0; eq_idx < Indices::numEq;++eq_idx) {
        out.setDerivative(eq_idx, in.derivative(eq_idx));
    }
    return out;
}

template<typename FluidSystem, typename Indices>
void
MultisegmentWellEval<FluidSystem,Indices>::
assembleAccelerationPressureLoss(const int seg,
                                 WellState<FluidSystem, Indices>& well_state)
{
    // Computes and assembles p-drop due to acceleration
    assert(seg != 0); // top segment can not enter here
    auto& segments = well_state.well(baseif_.indexOfWell()).segments;
    const auto& segment_set = this->segmentSet();

    // Add segment head
    const Scalar seg_area = segment_set[seg].crossArea();
    const EvalWell signed_velocity_head = segments_.accelerationPressureLossContribution(seg, seg_area);
    segments.pressure_drop_accel[seg] = signed_velocity_head.value();
    
    const int seg_upwind = segments_.upwinding_segment(seg);
    // acceleration term is *subtracted* from pressure equation
    MultisegmentWellAssemble(baseif_).
        assembleAccelerationTerm(seg, seg, seg_upwind, signed_velocity_head, linSys_);
    if (seg != seg_upwind) {// special treatment for reverse flow
        // extra derivatives are *added* to Jacobian (hence minus)
        const EvalWell extra_derivatives = -segments_.accelerationPressureLossContribution(seg, seg_area, /*extra_derivatives*/ true);
        MultisegmentWellAssemble(baseif_).
            assemblePressureEqExtraDerivatives(seg, seg_upwind, extra_derivatives, linSys_);
    }

    // Subtract inlet head(s), opposite signs from above
    for (const int inlet : segments_.inlets(seg)) {
        // area used in formula is max of areas
        const Scalar inlet_area = std::max(seg_area,
                                           static_cast<Scalar>(segment_set[inlet].crossArea()));
        const EvalWell signed_velocity_head_inlet = segments_.accelerationPressureLossContribution(inlet, inlet_area);
        segments.pressure_drop_accel[seg] -= signed_velocity_head_inlet.value();

        const int inlet_upwind = segments_.upwinding_segment(inlet); 
        MultisegmentWellAssemble(baseif_).
            assembleAccelerationTerm(seg, inlet, inlet_upwind, -signed_velocity_head_inlet, linSys_);
        if (inlet != inlet_upwind) {// special treatment for reverse flow
            // extra derivatives are *added* to Jacobian (hence minus minus)
            const EvalWell extra_derivatives_inlet = segments_.accelerationPressureLossContribution(inlet, inlet_area, /*extra_derivatives*/ true);
            // in this case inlet_upwind = seg
            MultisegmentWellAssemble(baseif_).
                assemblePressureEqExtraDerivatives(seg, inlet_upwind, extra_derivatives_inlet, linSys_);
        }
    }
}

template<typename FluidSystem, typename Indices>
void
MultisegmentWellEval<FluidSystem,Indices>::
assembleDefaultPressureEq(const int seg,
                          WellState<FluidSystem, Indices>& well_state,
                          const bool use_average_density)
{
    assert(seg != 0); // not top segment
    const int seg_upwind = segments_.upwinding_segment(seg);
    const bool reverseFlow = seg != seg_upwind; // special treatment for reverse flow

    // for top segment, the well control equation will be used.
    EvalWell pressure_equation = primary_variables_.getSegmentPressure(seg);
    EvalWell extra_derivatives;

    // we need to handle the pressure difference between the two segments
    // hydrostatic pressure loss is assembled seperately at the end
    // TODO: we might be able to add member variables to store these values, then we update well state
    // after converged

    auto& ws = well_state.well(baseif_.indexOfWell());
    auto& segments = ws.segments;

    if (this->frictionalPressureLossConsidered()) {
        const auto friction_pressure_drop = segments_.getFrictionPressureLoss(seg);
        if (reverseFlow){
            // call function once again to obtain/assemble remaining derivatives
            extra_derivatives = -segments_.getFrictionPressureLoss(seg, /*extra_reverse_flow_derivatives*/ true);
            MultisegmentWellAssemble(baseif_).
                assemblePressureEqExtraDerivatives(seg, seg_upwind, extra_derivatives, linSys_);
        }
        pressure_equation -= friction_pressure_drop;
        segments.pressure_drop_friction[seg] = friction_pressure_drop.value();
    }

    // contribution from the outlet segment
    const int outlet_segment_index = this->segmentNumberToIndex(this->segmentSet()[seg].outletSegment());
    const EvalWell outlet_pressure = primary_variables_.getSegmentPressure(outlet_segment_index);

    MultisegmentWellAssemble(baseif_).
        assemblePressureEq(seg, seg_upwind, outlet_segment_index,
                           pressure_equation, outlet_pressure, linSys_);

    assembleAccelerationAndHydroPressureLosses(seg, well_state, use_average_density);
}

template<typename FluidSystem, typename Indices>
void
MultisegmentWellEval<FluidSystem,Indices>::
assembleICDPressureEq(const int seg,
                      const UnitSystem& unit_system,
                      WellState<FluidSystem, Indices>& well_state,
                      const SummaryState& summary_state,
                      const bool use_average_density,
                      DeferredLogger& deferred_logger)
{
    // TODO: upwinding needs to be taken care of
    // top segment can not be a spiral ICD device
    assert(seg != 0);

    if (const auto& segment = this->segmentSet()[seg];
       (segment.segmentType() == Segment::SegmentType::VALVE) &&
       (segment.valve().status() == Opm::ICDStatus::SHUT) ) { // we use a zero rate equation to handle SHUT valve
        MultisegmentWellAssemble(baseif_).
            assembleTrivialEq(seg, this->primary_variables_.eval(seg)[WQTotal].value(), linSys_);

        auto& ws = well_state.well(baseif_.indexOfWell());
        ws.segments.pressure_drop_friction[seg] = 0.;
        return;
    }

    // the pressure equation is something like
    // p_seg - deltaP - p_outlet = 0.
    // the major part is how to calculate the deltaP
    const int seg_upwind = segments_.upwinding_segment(seg);
    const bool reverseFlow = seg != seg_upwind; // special treatment for reverse flow

    EvalWell pressure_equation = primary_variables_.getSegmentPressure(seg);

    EvalWell icd_pressure_drop;
    EvalWell extra_derivatives;
    switch(this->segmentSet()[seg].segmentType()) {
        case Segment::SegmentType::SICD :
            icd_pressure_drop = segments_.pressureDropSpiralICD(seg);
            if (reverseFlow){
                extra_derivatives = -segments_.pressureDropSpiralICD(seg, /*extra_reverse_flow_derivatives*/ true);
            }
            break;
        case Segment::SegmentType::AICD :
            icd_pressure_drop = segments_.pressureDropAutoICD(seg, unit_system);
            if (reverseFlow){
                extra_derivatives = -segments_.pressureDropAutoICD(seg, unit_system, /*extra_reverse_flow_derivatives*/ true);
            }
            break;
        case Segment::SegmentType::VALVE :
            icd_pressure_drop = segments_.pressureDropValve(seg, summary_state);
            if (reverseFlow){
                extra_derivatives = -segments_.pressureDropValve(seg, summary_state, /*extra_reverse_flow_derivatives*/ true);
            }
            break;
        default: {
            OPM_DEFLOG_THROW(std::runtime_error,
                             fmt::format("Segment {} for well {} is not of ICD type",
                                         this->segmentSet()[seg].segmentNumber(),
                                         baseif_.name()),
                             deferred_logger);
        }
    }
    if (reverseFlow){
        MultisegmentWellAssemble(baseif_).
            assemblePressureEqExtraDerivatives(seg, seg_upwind, extra_derivatives, linSys_);
    }

    pressure_equation = pressure_equation - icd_pressure_drop;
    auto& ws = well_state.well(baseif_.indexOfWell());
    ws.segments.pressure_drop_friction[seg] = icd_pressure_drop.value();

    // contribution from the outlet segment
    const int outlet_segment_index = this->segmentNumberToIndex(this->segmentSet()[seg].outletSegment());
    const EvalWell outlet_pressure = primary_variables_.getSegmentPressure(outlet_segment_index);

    MultisegmentWellAssemble(baseif_).
        assemblePressureEq(seg, seg_upwind, outlet_segment_index,
                           pressure_equation, outlet_pressure,
                           linSys_);

    assembleAccelerationAndHydroPressureLosses(seg, well_state, use_average_density);
}

template<typename FluidSystem, typename Indices>
void
MultisegmentWellEval<FluidSystem,Indices>::
assembleAccelerationAndHydroPressureLosses(const int seg,
                                           WellState<FluidSystem, Indices>& well_state,
                                           const bool use_average_density)
{
    if (this->accelerationalPressureLossConsidered()) {
        assembleAccelerationPressureLoss(seg, well_state);
    }

    // Since density derivatives are organized differently than what is required for assemblePressureEq,
    // this part needs to be assembled separately. Optionally use average density variant.
    const auto hydro_pressure_drop_seg = segments_.getHydroPressureLoss(seg, seg);
    auto& ws = well_state.well(baseif_.indexOfWell());
    auto& segments = ws.segments;
    if (!use_average_density){
        MultisegmentWellAssemble(baseif_).
            assembleHydroPressureLoss(seg, seg, hydro_pressure_drop_seg, linSys_);
        segments.pressure_drop_hydrostatic[seg] = hydro_pressure_drop_seg.value();
    } else {
        const int seg_outlet = this->segmentNumberToIndex(this->segmentSet()[seg].outletSegment());
        const auto hydro_pressure_drop_outlet = segments_.getHydroPressureLoss(seg, seg_outlet);
        MultisegmentWellAssemble(baseif_).
            assembleHydroPressureLoss(seg, seg, 0.5*hydro_pressure_drop_seg, linSys_);
        MultisegmentWellAssemble(baseif_).
            assembleHydroPressureLoss(seg, seg_outlet, 0.5*hydro_pressure_drop_outlet, linSys_);
        segments.pressure_drop_hydrostatic[seg] = 0.5*hydro_pressure_drop_seg.value() + 0.5*hydro_pressure_drop_outlet.value();
    }
}

template<typename FluidSystem, typename Indices>
void
MultisegmentWellEval<FluidSystem,Indices>::
assemblePressureEq(const int seg,
                   const UnitSystem& unit_system,
                   WellState<FluidSystem, Indices>& well_state,
                   const SummaryState& summary_state,
                   const bool use_average_density,
                   DeferredLogger& deferred_logger)
{
    switch(this->segmentSet()[seg].segmentType()) {
        case Segment::SegmentType::SICD :
        case Segment::SegmentType::AICD :
        case Segment::SegmentType::VALVE : {
            assembleICDPressureEq(seg, unit_system, well_state, summary_state, use_average_density, deferred_logger);
            break;
        }
        default :
            assembleDefaultPressureEq(seg, well_state, use_average_density);
    }
}

template<typename FluidSystem, typename Indices>
std::pair<bool, std::vector<typename FluidSystem::Scalar> >
MultisegmentWellEval<FluidSystem,Indices>::
getFiniteWellResiduals(const std::vector<Scalar>& B_avg,
                       DeferredLogger& deferred_logger) const
{
    assert(int(B_avg.size() ) == baseif_.numComponents());
    std::vector<Scalar> residuals(numWellEq + 1, 0.0);

    for (int seg = 0; seg < this->numberOfSegments(); ++seg) {
        for (int eq_idx = 0; eq_idx < numWellEq; ++eq_idx) {
            Scalar residual = 0.;
            if (eq_idx < baseif_.numComponents()) {
                residual = std::abs(linSys_.residual()[seg][eq_idx]) * B_avg[eq_idx];
            } else {
                if (seg > 0) {
                    residual = std::abs(linSys_.residual()[seg][eq_idx]);
                }
            }
            if (std::isnan(residual) || std::isinf(residual)) {
                deferred_logger.debug(fmt::format("nan or inf value for residual for well {} segment {} eq_idx {}",
                                                  baseif_.name(), seg, eq_idx));
                return {false, residuals};
            }

            if (residual > residuals[eq_idx]) {
                residuals[eq_idx] = residual;
            }
        }
    }

    // handling the control equation residual
    {
        const Scalar control_residual = std::abs(linSys_.residual()[0][numWellEq - 1]);
        if (std::isnan(control_residual) || std::isinf(control_residual)) {
           deferred_logger.debug(fmt::format("nan or inf value for control residual for well {}", baseif_.name()));
           return {false, residuals};
        }
        residuals[numWellEq] = control_residual;
    }

    return {true, residuals};
}

template<typename FluidSystem, typename Indices>
typename MultisegmentWellEval<FluidSystem,Indices>::Scalar
MultisegmentWellEval<FluidSystem,Indices>::
getControlTolerance(const WellState<FluidSystem, Indices>& well_state,
                    const Scalar tolerance_wells,
                    const Scalar tolerance_pressure_ms_wells,
                    DeferredLogger& deferred_logger) const
{
    Scalar control_tolerance = 0.;

    const int well_index = baseif_.indexOfWell();
    const auto& ws = well_state.well(well_index);
    if (baseif_.isInjector() )
    {
        auto current = ws.injection_cmode;
        switch(current) {
        case Well::InjectorCMode::THP:
            control_tolerance = tolerance_pressure_ms_wells;
            break;
        case Well::InjectorCMode::BHP:
            control_tolerance = tolerance_wells;
            break;
        case Well::InjectorCMode::RATE:
        case Well::InjectorCMode::RESV:
            control_tolerance = tolerance_wells;
            break;
        case Well::InjectorCMode::GRUP:
            control_tolerance = tolerance_wells;
            break;
        default:
            OPM_DEFLOG_THROW(std::runtime_error,
                             fmt::format("Unknown well control control types for well {}", baseif_.name()),
                             deferred_logger);
        }
    }

    if (baseif_.isProducer() )
    {
        auto current = ws.production_cmode;
        switch(current) {
        case Well::ProducerCMode::THP:
            control_tolerance = tolerance_pressure_ms_wells; // 0.1 bar
            break;
        case Well::ProducerCMode::BHP:
            control_tolerance = tolerance_wells; // 0.01 bar
            break;
        case Well::ProducerCMode::ORAT:
        case Well::ProducerCMode::WRAT:
        case Well::ProducerCMode::GRAT:
        case Well::ProducerCMode::LRAT:
        case Well::ProducerCMode::RESV:
        case Well::ProducerCMode::CRAT:
            control_tolerance = tolerance_wells; // smaller tolerance for rate control
            break;
        case Well::ProducerCMode::GRUP:
            control_tolerance = tolerance_wells; // smaller tolerance for rate control
            break;
        default:
            OPM_DEFLOG_THROW(std::runtime_error,
                             fmt::format("Unknown well control control types for well {}", baseif_.name()),
                             deferred_logger);
        }
    }

    return control_tolerance;
}

template<typename FluidSystem, typename Indices>
typename MultisegmentWellEval<FluidSystem,Indices>::Scalar
MultisegmentWellEval<FluidSystem,Indices>::
getResidualMeasureValue(const WellState<FluidSystem, Indices>& well_state,
                        const std::vector<Scalar>& residuals,
                        const Scalar tolerance_wells,
                        const Scalar tolerance_pressure_ms_wells,
                        DeferredLogger& deferred_logger) const
{
    assert(int(residuals.size()) == numWellEq + 1);

    const Scalar rate_tolerance = tolerance_wells;
    [[maybe_unused]] int count = 0;
    Scalar sum = 0;
    for (int eq_idx = 0; eq_idx < numWellEq - 1; ++eq_idx) {
        if (residuals[eq_idx] > rate_tolerance) {
            sum += residuals[eq_idx] / rate_tolerance;
            ++count;
        }
    }

    const Scalar pressure_tolerance = tolerance_pressure_ms_wells;
    if (residuals[SPres] > pressure_tolerance) {
        sum += residuals[SPres] / pressure_tolerance;
        ++count;
    }

    const Scalar control_tolerance = getControlTolerance(well_state,
                                                         tolerance_wells,
                                                         tolerance_pressure_ms_wells,
                                                         deferred_logger);
    if (residuals[SPres + 1] > control_tolerance) {
        sum += residuals[SPres + 1] / control_tolerance;
        ++count;
    }

    return sum;
}

#include <opm/simulators/utils/InstantiationIndicesMacros.hpp>

INSTANTIATE_TYPE_INDICES(MultisegmentWellEval, double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE_INDICES(MultisegmentWellEval, float)
#endif

} // namespace Opm
