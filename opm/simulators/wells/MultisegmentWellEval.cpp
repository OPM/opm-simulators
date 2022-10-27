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

#include <dune/istl/umfpack.hh>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/models/blackoil/blackoilindices.hh>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/simulators/linalg/bda/WellContributions.hpp>
#include <opm/simulators/timestepping/ConvergenceReport.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/MSWellHelpers.hpp>
#include <opm/simulators/wells/WellInterfaceIndices.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <cassert>

namespace Opm
{

template<typename FluidSystem, typename Indices, typename Scalar>
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
MultisegmentWellEval(WellInterfaceIndices<FluidSystem,Indices,Scalar>& baseif)
    : MultisegmentWellGeneric<Scalar>(baseif)
    , baseif_(baseif)
    , upwinding_segments_(this->numberOfSegments(), 0)
    , segment_densities_(this->numberOfSegments(), 0.0)
    , segment_mass_rates_(this->numberOfSegments(), 0.0)
    , segment_viscosities_(this->numberOfSegments(), 0.0)
    , segment_phase_densities_(this->numberOfSegments(), std::vector<EvalWell>(baseif_.numComponents(), 0.0)) // number of phase here?
    , segment_phase_fractions_(this->numberOfSegments(), std::vector<EvalWell>(baseif_.numComponents(), 0.0)) // number of phase here?
    , segment_phase_viscosities_(this->numberOfSegments(), std::vector<EvalWell>(baseif_.numComponents(), 0.0)) // number of phase here?
    , cell_perforation_depth_diffs_(baseif_.numPerfs(), 0.0)
    , cell_perforation_pressure_diffs_(baseif_.numPerfs(), 0.0)
{
}

template<typename FluidSystem, typename Indices, typename Scalar>
void
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
initMatrixAndVectors(const int num_cells) const
{
    duneB_.setBuildMode(OffDiagMatWell::row_wise);
    duneC_.setBuildMode(OffDiagMatWell::row_wise);
    duneD_.setBuildMode(DiagMatWell::row_wise);

    // set the size and patterns for all the matrices and vectors
    // [A C^T    [x    =  [ res
    //  B D] x_well]      res_well]

    // calculatiing the NNZ for duneD_
    // NNZ = number_of_segments + 2 * (number_of_inlets / number_of_outlets)
    {
        int nnz_d = this->numberOfSegments();
        for (const std::vector<int>& inlets : this->segment_inlets_) {
            nnz_d += 2 * inlets.size();
        }
        duneD_.setSize(this->numberOfSegments(), this->numberOfSegments(), nnz_d);
    }
    duneB_.setSize(this->numberOfSegments(), num_cells, baseif_.numPerfs());
    duneC_.setSize(this->numberOfSegments(), num_cells, baseif_.numPerfs());

    // we need to add the off diagonal ones
    for (auto row = duneD_.createbegin(), end = duneD_.createend(); row != end; ++row) {
        // the number of the row corrspnds to the segment now
        const int seg = row.index();
        // adding the item related to outlet relation
        const Segment& segment = this->segmentSet()[seg];
        const int outlet_segment_number = segment.outletSegment();
        if (outlet_segment_number > 0) { // if there is a outlet_segment
            const int outlet_segment_index = this->segmentNumberToIndex(outlet_segment_number);
            row.insert(outlet_segment_index);
        }

        // Add nonzeros for diagonal
        row.insert(seg);

        // insert the item related to its inlets
        for (const int& inlet : this->segment_inlets_[seg]) {
            row.insert(inlet);
        }
    }

    // make the C matrix
    for (auto row = duneC_.createbegin(), end = duneC_.createend(); row != end; ++row) {
        // the number of the row corresponds to the segment number now.
        for (const int& perf : this->segment_perforations_[row.index()]) {
            const int cell_idx = baseif_.cells()[perf];
            row.insert(cell_idx);
        }
    }

    // make the B^T matrix
    for (auto row = duneB_.createbegin(), end = duneB_.createend(); row != end; ++row) {
        // the number of the row corresponds to the segment number now.
        for (const int& perf : this->segment_perforations_[row.index()]) {
            const int cell_idx = baseif_.cells()[perf];
            row.insert(cell_idx);
        }
    }

    resWell_.resize(this->numberOfSegments());

    primary_variables_.resize(this->numberOfSegments());
    primary_variables_evaluation_.resize(this->numberOfSegments());
}

template<typename FluidSystem, typename Indices, typename Scalar>
void
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
initPrimaryVariablesEvaluation() const
{
    for (int seg = 0; seg < this->numberOfSegments(); ++seg) {
        for (int eq_idx = 0; eq_idx < numWellEq; ++eq_idx) {
            primary_variables_evaluation_[seg][eq_idx] = 0.0;
            primary_variables_evaluation_[seg][eq_idx].setValue(primary_variables_[seg][eq_idx]);
            primary_variables_evaluation_[seg][eq_idx].setDerivative(eq_idx + Indices::numEq, 1.0);
        }
    }
}

template<typename FluidSystem, typename Indices, typename Scalar>
void
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
checkConvergenceControlEq(const WellState& well_state,
                          ConvergenceReport& report,
                          const double tolerance_pressure_ms_wells,
                          const double tolerance_wells,
                          const double max_residual_allowed,
                          DeferredLogger& deferred_logger) const
{
    double control_tolerance = 0.;
    using CR = ConvergenceReport;
    CR::WellFailure::Type ctrltype = CR::WellFailure::Type::Invalid;

    const int well_index = baseif_.indexOfWell();
    const auto& ws = well_state.well(well_index);
    if (baseif_.isInjector() )
    {
        auto current = ws.injection_cmode;
        switch(current) {
        case Well::InjectorCMode::THP:
            ctrltype = CR::WellFailure::Type::ControlTHP;
            control_tolerance = tolerance_pressure_ms_wells;
            break;
        case Well::InjectorCMode::BHP:
            ctrltype = CR::WellFailure::Type::ControlBHP;
            control_tolerance = tolerance_pressure_ms_wells;
            break;
        case Well::InjectorCMode::RATE:
        case Well::InjectorCMode::RESV:
            ctrltype = CR::WellFailure::Type::ControlRate;
            control_tolerance = tolerance_wells;
            break;
        case Well::InjectorCMode::GRUP:
            ctrltype = CR::WellFailure::Type::ControlRate;
            control_tolerance = tolerance_wells;
            break;
        default:
            OPM_DEFLOG_THROW(std::runtime_error, "Unknown well control control types for well " << baseif_.name(), deferred_logger);
        }
    }

    if (baseif_.isProducer() )
    {
        auto current = ws.production_cmode;
        switch(current) {
        case Well::ProducerCMode::THP:
            ctrltype = CR::WellFailure::Type::ControlTHP;
            control_tolerance = tolerance_pressure_ms_wells;
            break;
        case Well::ProducerCMode::BHP:
            ctrltype = CR::WellFailure::Type::ControlBHP;
            control_tolerance = tolerance_pressure_ms_wells;
            break;
        case Well::ProducerCMode::ORAT:
        case Well::ProducerCMode::WRAT:
        case Well::ProducerCMode::GRAT:
        case Well::ProducerCMode::LRAT:
        case Well::ProducerCMode::RESV:
        case Well::ProducerCMode::CRAT:
            ctrltype = CR::WellFailure::Type::ControlRate;
            control_tolerance = tolerance_wells;
            break;
        case Well::ProducerCMode::GRUP:
            ctrltype = CR::WellFailure::Type::ControlRate;
            control_tolerance = tolerance_wells;
            break;
        default:
            OPM_DEFLOG_THROW(std::runtime_error, "Unknown well control control types for well " << baseif_.name(), deferred_logger);
        }
    }

    const double well_control_residual = std::abs(resWell_[0][SPres])/this->bhp_control_scaling_;
    const int dummy_component = -1;
    if (std::isnan(well_control_residual)) {
        report.setWellFailed({ctrltype, CR::Severity::NotANumber, dummy_component, baseif_.name()});
    } else if (well_control_residual > max_residual_allowed * 10.) {
        report.setWellFailed({ctrltype, CR::Severity::TooLarge, dummy_component, baseif_.name()});
    } else if ( well_control_residual > control_tolerance) {
        report.setWellFailed({ctrltype, CR::Severity::Normal, dummy_component, baseif_.name()});
    }
}

template<typename FluidSystem, typename Indices, typename Scalar>
ConvergenceReport
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
getWellConvergence(const WellState& well_state,
                   const std::vector<double>& B_avg,
                   DeferredLogger& deferred_logger,
                   const double max_residual_allowed,
                   const double tolerance_wells,
                   const double relaxed_inner_tolerance_flow_ms_well,
                   const double tolerance_pressure_ms_wells,
                   const double relaxed_inner_tolerance_pressure_ms_well,
                   const bool relax_tolerance) const
{
    assert(int(B_avg.size()) == baseif_.numComponents());

    // checking if any residual is NaN or too large. The two large one is only handled for the well flux
    std::vector<std::vector<double>> abs_residual(this->numberOfSegments(), std::vector<double>(numWellEq, 0.0));
    for (int seg = 0; seg < this->numberOfSegments(); ++seg) {
        for (int eq_idx = 0; eq_idx < numWellEq; ++eq_idx) {
            abs_residual[seg][eq_idx] = std::abs(resWell_[seg][eq_idx]);
        }
    }

    std::vector<double> maximum_residual(numWellEq, 0.0);

    ConvergenceReport report;
    // TODO: the following is a little complicated, maybe can be simplified in some way?
    for (int eq_idx = 0; eq_idx < numWellEq; ++eq_idx) {
        for (int seg = 0; seg < this->numberOfSegments(); ++seg) {
            if (eq_idx < baseif_.numComponents()) { // phase or component mass equations
                const double flux_residual = B_avg[eq_idx] * abs_residual[seg][eq_idx];
                if (flux_residual > maximum_residual[eq_idx]) {
                    maximum_residual[eq_idx] = flux_residual;
                }
            } else { // pressure or control equation
                // for the top segment (seg == 0), it is control equation, will be checked later separately
                if (seg > 0) {
                    // Pressure equation
                    const double pressure_residual = abs_residual[seg][eq_idx];
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
            const double flux_residual = maximum_residual[eq_idx];
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
            const double pressure_residual = maximum_residual[eq_idx];
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

    checkConvergenceControlEq(well_state,
                              report,
                              tolerance_pressure_ms_wells,
                              tolerance_wells,
                              max_residual_allowed,
                              deferred_logger);

    return report;
}

template<typename FluidSystem, typename Indices, typename Scalar>
void
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
processFractions(const int seg) const
{
    static constexpr int Water = BlackoilPhases::Aqua;
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Gas = BlackoilPhases::Vapour;

    const PhaseUsage& pu = baseif_.phaseUsage();

    std::vector<double> fractions(baseif_.numPhases(), 0.0);

    assert( FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) );
    const int oil_pos = pu.phase_pos[Oil];
    fractions[oil_pos] = 1.0;

    if ( FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) ) {
        const int water_pos = pu.phase_pos[Water];
        fractions[water_pos] = primary_variables_[seg][WFrac];
        fractions[oil_pos] -= fractions[water_pos];
    }

    if ( FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) ) {
        const int gas_pos = pu.phase_pos[Gas];
        fractions[gas_pos] = primary_variables_[seg][GFrac];
        fractions[oil_pos] -= fractions[gas_pos];
    }

    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
        const int water_pos = pu.phase_pos[Water];
        if (fractions[water_pos] < 0.0) {
            if ( FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) ) {
                fractions[pu.phase_pos[Gas]] /= (1.0 - fractions[water_pos]);
            }
            fractions[oil_pos] /= (1.0 - fractions[water_pos]);
            fractions[water_pos] = 0.0;
        }
    }

    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
        const int gas_pos = pu.phase_pos[Gas];
        if (fractions[gas_pos] < 0.0) {
            if ( FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) ) {
                fractions[pu.phase_pos[Water]] /= (1.0 - fractions[gas_pos]);
            }
            fractions[oil_pos] /= (1.0 - fractions[gas_pos]);
            fractions[gas_pos] = 0.0;
        }
    }

    if (fractions[oil_pos] < 0.0) {
        if ( FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) ) {
            fractions[pu.phase_pos[Water]] /= (1.0 - fractions[oil_pos]);
        }
        if ( FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) ) {
            fractions[pu.phase_pos[Gas]] /= (1.0 - fractions[oil_pos]);
        }
        fractions[oil_pos] = 0.0;
    }

    if ( FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) ) {
        primary_variables_[seg][WFrac] = fractions[pu.phase_pos[Water]];
    }

    if ( FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) ) {
        primary_variables_[seg][GFrac] = fractions[pu.phase_pos[Gas]];
    }
}

template<typename FluidSystem, typename Indices, typename Scalar>
void
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
updatePrimaryVariablesNewton(const BVectorWell& dwells,
                             const double relaxation_factor,
                             const double dFLimit,
                             const double max_pressure_change) const
{
    const std::vector<std::array<double, numWellEq> > old_primary_variables = primary_variables_;

    for (int seg = 0; seg < this->numberOfSegments(); ++seg) {
        if (has_wfrac_variable) {
            const int sign = dwells[seg][WFrac] > 0. ? 1 : -1;
            const double dx_limited = sign * std::min(std::abs(dwells[seg][WFrac]) * relaxation_factor, dFLimit);
            primary_variables_[seg][WFrac] = old_primary_variables[seg][WFrac] - dx_limited;
        }

        if (has_gfrac_variable) {
            const int sign = dwells[seg][GFrac] > 0. ? 1 : -1;
            const double dx_limited = sign * std::min(std::abs(dwells[seg][GFrac]) * relaxation_factor, dFLimit);
            primary_variables_[seg][GFrac] = old_primary_variables[seg][GFrac] - dx_limited;
        }

        // handling the overshooting or undershooting of the fractions
        processFractions(seg);

        // update the segment pressure
        {
            const int sign = dwells[seg][SPres] > 0.? 1 : -1;
            const double dx_limited = sign * std::min(std::abs(dwells[seg][SPres]) * relaxation_factor, max_pressure_change/this->bhp_scaling_);
            primary_variables_[seg][SPres] = std::max( old_primary_variables[seg][SPres] - dx_limited, 1e5/this->bhp_scaling_);
        }

        // update the total rate // TODO: should we have a limitation of the total rate change?
        {
            primary_variables_[seg][WQTotal] = old_primary_variables[seg][WQTotal] - relaxation_factor * dwells[seg][WQTotal];

            // make sure that no injector produce and no producer inject
            if (seg == 0) {
                if (baseif_.isInjector()) {
                    primary_variables_[seg][WQTotal] = std::max( primary_variables_[seg][WQTotal], 0.0);
                } else {
                    primary_variables_[seg][WQTotal] = std::min( primary_variables_[seg][WQTotal], 0.0);
                }
            }
        }
    }
}

template<typename FluidSystem, typename Indices, typename Scalar>
void
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
updatePrimaryVariables(const WellState& well_state) const
{
    static constexpr int Water = BlackoilPhases::Aqua;
    static constexpr int Gas = BlackoilPhases::Vapour;

    // TODO: to test using rate conversion coefficients to see if it will be better than
    // this default one
    if (!baseif_.isOperableAndSolvable() && !baseif_.wellIsStopped()) return;

    const Well& well = baseif_.wellEcl();

    // the index of the top segment in the WellState
    const auto& ws = well_state.well(baseif_.indexOfWell());
    const auto& segments = ws.segments;
    const auto& segment_rates = segments.rates;
    const auto& segment_pressure = segments.pressure;
    const PhaseUsage& pu = baseif_.phaseUsage();

    for (int seg = 0; seg < this->numberOfSegments(); ++seg) {
        // calculate the total rate for each segment
        double total_seg_rate = 0.0;
        // the segment pressure
        primary_variables_[seg][SPres] = segment_pressure[seg]/this->bhp_scaling_;
        // TODO: under what kind of circustances, the following will be wrong?
        // the definition of g makes the gas phase is always the last phase
        for (int p = 0; p < baseif_.numPhases(); p++) {
            total_seg_rate += baseif_.scalingFactor(p) * segment_rates[baseif_.numPhases() * seg + p];
        }

        if (seg == 0) {
            if (baseif_.isInjector()) {
                total_seg_rate = std::max(total_seg_rate, 0.);
            } else {
                total_seg_rate = std::min(total_seg_rate, 0.);
            }
        }
        primary_variables_[seg][WQTotal] = total_seg_rate/this->rate_scaling_;
        if (std::abs(total_seg_rate) > 0.) {
            if (has_wfrac_variable) {
                const int water_pos = pu.phase_pos[Water];
                primary_variables_[seg][WFrac] = baseif_.scalingFactor(water_pos) * segment_rates[baseif_.numPhases() * seg + water_pos] / total_seg_rate;
            }
            if (has_gfrac_variable) {
                const int gas_pos = pu.phase_pos[Gas];
                primary_variables_[seg][GFrac] = baseif_.scalingFactor(gas_pos) * segment_rates[baseif_.numPhases() * seg + gas_pos] / total_seg_rate;
            }
        } else { // total_seg_rate == 0
            if (baseif_.isInjector()) {
                // only single phase injection handled
                auto phase = well.getInjectionProperties().injectorType;

                if (has_wfrac_variable) {
                    if (phase == InjectorType::WATER) {
                        primary_variables_[seg][WFrac] = 1.0;
                    } else {
                        primary_variables_[seg][WFrac] = 0.0;
                    }
                }

                if (has_gfrac_variable) {
                    if (phase == InjectorType::GAS) {
                        primary_variables_[seg][GFrac] = 1.0;
                    } else {
                        primary_variables_[seg][GFrac] = 0.0;
                    }
                }

            } else if (baseif_.isProducer()) { // producers
                if (has_wfrac_variable) {
                    primary_variables_[seg][WFrac] = 1.0 / baseif_.numPhases();
                }

                if (has_gfrac_variable) {
                    primary_variables_[seg][GFrac] = 1.0 / baseif_.numPhases();
                }
            }
        }
    }
}

template<typename FluidSystem, typename Indices, typename Scalar>
void
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
recoverSolutionWell(const BVector& x, BVectorWell& xw) const
{
    if (!baseif_.isOperableAndSolvable() && !baseif_.wellIsStopped()) return;

    BVectorWell resWell = resWell_;
    // resWell = resWell - B * x
    duneB_.mmv(x, resWell);
    // xw = D^-1 * resWell
    xw = mswellhelpers::applyUMFPack(duneD_, duneDSolver_, resWell);
}

template<typename FluidSystem, typename Indices, typename Scalar>
typename MultisegmentWellEval<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
volumeFraction(const int seg,
               const unsigned compIdx) const
{
    if (has_wfrac_variable && compIdx == Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx)) {
        return primary_variables_evaluation_[seg][WFrac];
    }

    if (has_gfrac_variable && compIdx == Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx)) {
        return primary_variables_evaluation_[seg][GFrac];
    }

    // Oil fraction
    EvalWell oil_fraction = 1.0;
    if (has_wfrac_variable) {
        oil_fraction -= primary_variables_evaluation_[seg][WFrac];
    }

    if (has_gfrac_variable) {
        oil_fraction -= primary_variables_evaluation_[seg][GFrac];
    }
    /* if (has_solvent) {
        oil_fraction -= primary_variables_evaluation_[seg][SFrac];
    } */
    return oil_fraction;
}

template<typename FluidSystem, typename Indices, typename Scalar>
typename MultisegmentWellEval<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
volumeFractionScaled(const int seg,
                     const int comp_idx) const
{
    // For reservoir rate control, the distr in well control is used for the
    // rate conversion coefficients. For the injection well, only the distr of the injection
    // phase is not zero.
    const double scale = baseif_.scalingFactor(baseif_.ebosCompIdxToFlowCompIdx(comp_idx));
    if (scale > 0.) {
        return volumeFraction(seg, comp_idx) / scale;
    }

    return volumeFraction(seg, comp_idx);
}

template<typename FluidSystem, typename Indices, typename Scalar>
typename MultisegmentWellEval<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
surfaceVolumeFraction(const int seg,
                      const int comp_idx) const
{
    EvalWell sum_volume_fraction_scaled = 0.;
    for (int idx = 0; idx < baseif_.numComponents(); ++idx) {
        sum_volume_fraction_scaled += volumeFractionScaled(seg, idx);
    }

    assert(sum_volume_fraction_scaled.value() != 0.);

    return volumeFractionScaled(seg, comp_idx) / sum_volume_fraction_scaled;
}

template<typename FluidSystem, typename Indices, typename Scalar>
typename MultisegmentWellEval<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
getSegmentRateUpwinding(const int seg,
                        const size_t comp_idx) const
{
    const int seg_upwind = upwinding_segments_[seg];
    // the result will contain the derivative with respect to WQTotal in segment seg,
    // and the derivatives with respect to WFrac GFrac in segment seg_upwind.
    // the derivative with respect to SPres should be zero.
    if (seg == 0 && baseif_.isInjector()) {
        const Well& well = baseif_.wellEcl();
        auto phase = well.getInjectionProperties().injectorType;

        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)
                && Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx) == comp_idx
                && phase == InjectorType::WATER)
            return primary_variables_evaluation_[seg][WQTotal]*this->rate_scaling_ / baseif_.scalingFactor(baseif_.ebosCompIdxToFlowCompIdx(comp_idx));


        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)
                && Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx) == comp_idx
                && phase == InjectorType::OIL)
            return primary_variables_evaluation_[seg][WQTotal]*this->rate_scaling_ / baseif_.scalingFactor(baseif_.ebosCompIdxToFlowCompIdx(comp_idx));

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)
                && Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx) == comp_idx
                && phase == InjectorType::GAS)
            return primary_variables_evaluation_[seg][WQTotal]*this->rate_scaling_ / baseif_.scalingFactor(baseif_.ebosCompIdxToFlowCompIdx(comp_idx));

        return 0.0;
    }

    const EvalWell segment_rate = primary_variables_evaluation_[seg][WQTotal] * volumeFractionScaled(seg_upwind, comp_idx);

    assert(segment_rate.derivative(SPres + Indices::numEq) == 0.);

    return segment_rate*this->rate_scaling_;
}

template<typename FluidSystem, typename Indices, typename Scalar>
typename MultisegmentWellEval<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
extendEval(const Eval& in) const
{
    EvalWell out = 0.0;
    out.setValue(in.value());
    for(int eq_idx = 0; eq_idx < Indices::numEq;++eq_idx) {
        out.setDerivative(eq_idx, in.derivative(eq_idx));
    }
    return out;
}

template<typename FluidSystem, typename Indices, typename Scalar>
void
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
computeSegmentFluidProperties(const EvalWell& temperature,
                              const EvalWell& saltConcentration,
                              int pvt_region_index,
                              DeferredLogger& deferred_logger)
{
    std::vector<double> surf_dens(baseif_.numComponents());
    // Surface density.
    for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
        if (!FluidSystem::phaseIsActive(phaseIdx)) {
            continue;
        }

        const unsigned compIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
        surf_dens[compIdx] = FluidSystem::referenceDensity( phaseIdx, pvt_region_index);
    }

    for (int seg = 0; seg < this->numberOfSegments(); ++seg) {
        // the compostion of the components inside wellbore under surface condition
        std::vector<EvalWell> mix_s(baseif_.numComponents(), 0.0);
        for (int comp_idx = 0; comp_idx < baseif_.numComponents(); ++comp_idx) {
            mix_s[comp_idx] = surfaceVolumeFraction(seg, comp_idx);
        }

        std::vector<EvalWell> b(baseif_.numComponents(), 0.0);
        std::vector<EvalWell> visc(baseif_.numComponents(), 0.0);
        std::vector<EvalWell>& phase_densities = segment_phase_densities_[seg];

        const EvalWell seg_pressure = getSegmentPressure(seg);
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
            b[waterCompIdx] =
                FluidSystem::waterPvt().inverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure, saltConcentration);
            visc[waterCompIdx] =
                FluidSystem::waterPvt().viscosity(pvt_region_index, temperature, seg_pressure, saltConcentration);
            // TODO: double check here
            // TODO: should not we use phaseIndex here?
            phase_densities[waterCompIdx] = b[waterCompIdx] * surf_dens[waterCompIdx];
        }

        EvalWell rv(0.0);
        EvalWell rvw(0.0);
        // gas phase
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
                const EvalWell rvmax = FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvt_region_index, temperature, seg_pressure);
                if (mix_s[oilCompIdx] > 0.0) {
                    if (mix_s[gasCompIdx] > 0.0) {
                        rv = mix_s[oilCompIdx] / mix_s[gasCompIdx];
                    }

                    if (rv > rvmax) {
                        rv = rvmax;
                    }
                    b[gasCompIdx] =
                        FluidSystem::gasPvt().inverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure, rv, rvw);
                    visc[gasCompIdx] =
                        FluidSystem::gasPvt().viscosity(pvt_region_index, temperature, seg_pressure, rv, rvw);
                    phase_densities[gasCompIdx] = b[gasCompIdx] * surf_dens[gasCompIdx]
                                                + rv * b[gasCompIdx] * surf_dens[oilCompIdx];
                } else { // no oil exists
                    b[gasCompIdx] =
                        FluidSystem::gasPvt().saturatedInverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure);
                    visc[gasCompIdx] =
                        FluidSystem::gasPvt().saturatedViscosity(pvt_region_index, temperature, seg_pressure);
                    phase_densities[gasCompIdx] = b[gasCompIdx] * surf_dens[gasCompIdx];
                }
            } else { // no Liquid phase
                // it is the same with zero mix_s[Oil]
                b[gasCompIdx] =
                    FluidSystem::gasPvt().saturatedInverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure);
                visc[gasCompIdx] =
                    FluidSystem::gasPvt().saturatedViscosity(pvt_region_index, temperature, seg_pressure);
            }
        }

        EvalWell rs(0.0);
        // oil phase
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                const EvalWell rsmax = FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvt_region_index, temperature, seg_pressure);
                if (mix_s[gasCompIdx] > 0.0) {
                    if (mix_s[oilCompIdx] > 0.0) {
                        rs = mix_s[gasCompIdx] / mix_s[oilCompIdx];
                    }

                    if (rs > rsmax) {
                        rs = rsmax;
                    }
                    b[oilCompIdx] =
                        FluidSystem::oilPvt().inverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure, rs);
                    visc[oilCompIdx] =
                        FluidSystem::oilPvt().viscosity(pvt_region_index, temperature, seg_pressure, rs);
                    phase_densities[oilCompIdx] = b[oilCompIdx] * surf_dens[oilCompIdx]
                                                + rs * b[oilCompIdx] * surf_dens[gasCompIdx];
                } else { // no oil exists
                    b[oilCompIdx] =
                        FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure);
                    visc[oilCompIdx] =
                        FluidSystem::oilPvt().saturatedViscosity(pvt_region_index, temperature, seg_pressure);
                    phase_densities[oilCompIdx] = b[oilCompIdx] * surf_dens[oilCompIdx];
                }
            } else { // no Liquid phase
                // it is the same with zero mix_s[Oil]
                b[oilCompIdx] =
                    FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure);
                visc[oilCompIdx] =
                    FluidSystem::oilPvt().saturatedViscosity(pvt_region_index, temperature, seg_pressure);
            }
        }

        segment_phase_viscosities_[seg] = visc;

        std::vector<EvalWell> mix(mix_s);
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
            const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);

            const EvalWell d = 1.0 - rs * rv;
            if (d <= 0.0) {
                std::ostringstream sstr;
                sstr << "Problematic d value " << d << " obtained for well " << baseif_.name()
                     << " during segment density calculations with rs " << rs
                     << ", rv " << rv << " and pressure " << seg_pressure
                     << " obtaining d " << d
                     << " Continue as if no dissolution (rs = 0) and vaporization (rv = 0) "
                     << " for this connection.";
                deferred_logger.debug(sstr.str());
            } else {
                if (rs > 0.0) {
                    mix[gasCompIdx] = (mix_s[gasCompIdx] - mix_s[oilCompIdx] * rs) / d;
                }
                if (rv > 0.0) {
                    mix[oilCompIdx] = (mix_s[oilCompIdx] - mix_s[gasCompIdx] * rv) / d;
                }
            }
        }

        EvalWell volrat(0.0);
        for (int comp_idx = 0; comp_idx < baseif_.numComponents(); ++comp_idx) {
            volrat += mix[comp_idx] / b[comp_idx];
        }

        this->segment_viscosities_[seg] = 0.;
        // calculate the average viscosity
        for (int comp_idx = 0; comp_idx < baseif_.numComponents(); ++comp_idx) {
            const EvalWell fraction =  mix[comp_idx] / b[comp_idx] / volrat;
            // TODO: a little more work needs to be done to handle the negative fractions here
            this->segment_phase_fractions_[seg][comp_idx] = fraction; // >= 0.0 ? fraction : 0.0;
            this->segment_viscosities_[seg] += visc[comp_idx] * this->segment_phase_fractions_[seg][comp_idx];
        }

        EvalWell density(0.0);
        for (int comp_idx = 0; comp_idx < baseif_.numComponents(); ++comp_idx) {
            density += surf_dens[comp_idx] * mix_s[comp_idx];
        }
        this->segment_densities_[seg] = density / volrat;

        // calculate the mass rates
        segment_mass_rates_[seg] = 0.;
        for (int comp_idx = 0; comp_idx < baseif_.numComponents(); ++comp_idx) {
            const EvalWell rate = getSegmentRateUpwinding(seg, comp_idx);
            this->segment_mass_rates_[seg] += rate * surf_dens[comp_idx];
        }
    }
}

template<typename FluidSystem, typename Indices, typename Scalar>
typename MultisegmentWellEval<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
getSegmentPressure(const int seg) const
{
    return primary_variables_evaluation_[seg][SPres]*this->bhp_scaling_;
}

template<typename FluidSystem, typename Indices, typename Scalar>
typename MultisegmentWellEval<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
getBhp() const
{
    return getSegmentPressure(0);
}

template<typename FluidSystem, typename Indices, typename Scalar>
typename MultisegmentWellEval<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
getSegmentRate(const int seg,
               const int comp_idx) const
{
    return primary_variables_evaluation_[seg][WQTotal] * volumeFractionScaled(seg, comp_idx)*this->rate_scaling_;
}

template<typename FluidSystem, typename Indices, typename Scalar>
typename MultisegmentWellEval<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
getQs(const int comp_idx) const
{
    return getSegmentRate(0, comp_idx);
}

template<typename FluidSystem, typename Indices, typename Scalar>
typename MultisegmentWellEval<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
getSegmentWQTotal(const int seg) const
{
    return primary_variables_evaluation_[seg][WQTotal];
}

template<typename FluidSystem, typename Indices, typename Scalar>
typename MultisegmentWellEval<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
getWQTotal() const
{
    return getSegmentWQTotal(0);
}

template<typename FluidSystem, typename Indices, typename Scalar>
typename MultisegmentWellEval<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
getHydroPressureLoss(const int seg) const
{
    return segment_densities_[seg] * baseif_.gravity() * this->segment_depth_diffs_[seg];
}

template<typename FluidSystem, typename Indices, typename Scalar>
typename MultisegmentWellEval<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
getFrictionPressureLoss(const int seg) const
{
    const EvalWell mass_rate = segment_mass_rates_[seg];
    const int seg_upwind = upwinding_segments_[seg];
    EvalWell density = segment_densities_[seg_upwind];
    EvalWell visc = segment_viscosities_[seg_upwind];
    // WARNING
    // We disregard the derivatives from the upwind density to make sure derivatives
    // wrt. to different segments dont get mixed.
    if (seg != seg_upwind) {
        density.clearDerivatives();
        visc.clearDerivatives();
    }
    const int outlet_segment_index = this->segmentNumberToIndex(this->segmentSet()[seg].outletSegment());
    const double length = this->segmentSet()[seg].totalLength() - this->segmentSet()[outlet_segment_index].totalLength();
    assert(length > 0.);
    const double roughness = this->segmentSet()[seg].roughness();
    const double area = this->segmentSet()[seg].crossArea();
    const double diameter = this->segmentSet()[seg].internalDiameter();

    const double sign = mass_rate < 0. ? 1.0 : - 1.0;

    return sign * mswellhelpers::frictionPressureLoss(length, diameter, area, roughness, density, mass_rate, visc);
}

template<typename FluidSystem, typename Indices, typename Scalar>
typename MultisegmentWellEval<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
pressureDropSpiralICD(const int seg) const
{
    const SICD& sicd = this->segmentSet()[seg].spiralICD();

    const int seg_upwind = upwinding_segments_[seg];
    const std::vector<EvalWell>& phase_fractions = segment_phase_fractions_[seg_upwind];
    const std::vector<EvalWell>& phase_viscosities = segment_phase_viscosities_[seg_upwind];

    EvalWell water_fraction = 0.;
    EvalWell water_viscosity = 0.;
    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
        const int water_pos = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
        water_fraction = phase_fractions[water_pos];
        water_viscosity = phase_viscosities[water_pos];
    }

    EvalWell oil_fraction = 0.;
    EvalWell oil_viscosity = 0.;
    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
        const int oil_pos = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
        oil_fraction = phase_fractions[oil_pos];
        oil_viscosity = phase_viscosities[oil_pos];
    }

    EvalWell gas_fraction = 0.;
    EvalWell gas_viscosity = 0.;
    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
        const int gas_pos = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
        gas_fraction = phase_fractions[gas_pos];
        gas_viscosity = phase_viscosities[gas_pos];
    }

    EvalWell density = segment_densities_[seg_upwind];
    // WARNING
    // We disregard the derivatives from the upwind density to make sure derivatives
    // wrt. to different segments dont get mixed.
    if (seg != seg_upwind) {
        water_fraction.clearDerivatives();
        water_viscosity.clearDerivatives();
        oil_fraction.clearDerivatives();
        oil_viscosity.clearDerivatives();
        gas_fraction.clearDerivatives();
        gas_viscosity.clearDerivatives();
        density.clearDerivatives();
    }

    const EvalWell liquid_emulsion_viscosity = mswellhelpers::emulsionViscosity(water_fraction, water_viscosity,
                                                                                oil_fraction, oil_viscosity, sicd);
    const EvalWell mixture_viscosity = (water_fraction + oil_fraction) * liquid_emulsion_viscosity + gas_fraction * gas_viscosity;

    const EvalWell reservoir_rate = segment_mass_rates_[seg] / density;

    const EvalWell reservoir_rate_icd = reservoir_rate * sicd.scalingFactor();

    const double viscosity_cali = sicd.viscosityCalibration();

    using MathTool = MathToolbox<EvalWell>;

    const double density_cali = sicd.densityCalibration();
    const EvalWell temp_value1 = MathTool::pow(density / density_cali, 0.75);
    const EvalWell temp_value2 = MathTool::pow(mixture_viscosity / viscosity_cali, 0.25);

    // formulation before 2016, base_strength is used
    // const double base_strength = sicd.strength() / density_cali;
    // formulation since 2016, strength is used instead
    const double strength = sicd.strength();

    const double sign = reservoir_rate_icd <= 0. ? 1.0 : -1.0;

    return sign * temp_value1 * temp_value2 * strength * reservoir_rate_icd * reservoir_rate_icd;
}

template<typename FluidSystem, typename Indices, typename Scalar>
typename MultisegmentWellEval<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
pressureDropAutoICD(const int seg,
                    const UnitSystem& unit_system) const
{
    const AutoICD& aicd = this->segmentSet()[seg].autoICD();

    const int seg_upwind = this->upwinding_segments_[seg];
    const std::vector<EvalWell>& phase_fractions = this->segment_phase_fractions_[seg_upwind];
    const std::vector<EvalWell>& phase_viscosities = this->segment_phase_viscosities_[seg_upwind];
    const std::vector<EvalWell>& phase_densities = this->segment_phase_densities_[seg_upwind];

    EvalWell water_fraction = 0.;
    EvalWell water_viscosity = 0.;
    EvalWell water_density = 0.;
    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
        const int water_pos = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
        water_fraction = phase_fractions[water_pos];
        water_viscosity = phase_viscosities[water_pos];
        water_density = phase_densities[water_pos];
    }

    EvalWell oil_fraction = 0.;
    EvalWell oil_viscosity = 0.;
    EvalWell oil_density = 0.;
    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
        const int oil_pos = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
        oil_fraction = phase_fractions[oil_pos];
        oil_viscosity = phase_viscosities[oil_pos];
        oil_density = phase_densities[oil_pos];
    }

    EvalWell gas_fraction = 0.;
    EvalWell gas_viscosity = 0.;
    EvalWell gas_density = 0.;
    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
        const int gas_pos = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
        gas_fraction = phase_fractions[gas_pos];
        gas_viscosity = phase_viscosities[gas_pos];
        gas_density = phase_densities[gas_pos];
    }

    EvalWell density = segment_densities_[seg_upwind];
    // WARNING
    // We disregard the derivatives from the upwind density to make sure derivatives
    // wrt. to different segments dont get mixed.
    if (seg != seg_upwind) {
        water_fraction.clearDerivatives();
        water_viscosity.clearDerivatives();
        water_density.clearDerivatives();
        oil_fraction.clearDerivatives();
        oil_viscosity.clearDerivatives();
        oil_density.clearDerivatives();
        gas_fraction.clearDerivatives();
        gas_viscosity.clearDerivatives();
        gas_density.clearDerivatives();
        density.clearDerivatives();
    }

    using MathTool = MathToolbox<EvalWell>;
    const EvalWell mixture_viscosity = MathTool::pow(water_fraction, aicd.waterViscExponent()) * water_viscosity
                                     + MathTool::pow(oil_fraction, aicd.oilViscExponent()) * oil_viscosity
                                     + MathTool::pow(gas_fraction, aicd.gasViscExponent()) * gas_viscosity;

    const EvalWell mixture_density = MathTool::pow(water_fraction, aicd.waterDensityExponent()) * water_density
                                   + MathTool::pow(oil_fraction, aicd.oilDensityExponent()) * oil_density
                                   + MathTool::pow(gas_fraction, aicd.gasDensityExponent()) * gas_density;

    const double rho_reference = aicd.densityCalibration();
    const double visc_reference = aicd.viscosityCalibration();
    const auto volume_rate_icd = this->segment_mass_rates_[seg] * aicd.scalingFactor() / mixture_density;
    const double sign = volume_rate_icd <= 0. ? 1.0 : -1.0;
    // convert 1 unit volume rate
    using M  = UnitSystem::measure;
    const double unit_volume_rate = unit_system.to_si(M::geometric_volume_rate, 1.);

    // TODO: we did not consider the maximum allowed rate here
    const auto result = sign / rho_reference * mixture_density * mixture_density
                      * MathTool::pow(visc_reference/mixture_viscosity, aicd.viscExponent())
                      * aicd.strength() * MathTool::pow( -sign * volume_rate_icd, aicd.flowRateExponent())
                      * std::pow(unit_volume_rate, (2. - aicd.flowRateExponent())) ;
    return result;
}

template<typename FluidSystem, typename Indices, typename Scalar>
typename MultisegmentWellEval<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
pressureDropValve(const int seg) const
{
    const Valve& valve = this->segmentSet()[seg].valve();

    const EvalWell& mass_rate = segment_mass_rates_[seg];
    const int seg_upwind = upwinding_segments_[seg];
    EvalWell visc = segment_viscosities_[seg_upwind];
    EvalWell density = segment_densities_[seg_upwind];
    // WARNING
    // We disregard the derivatives from the upwind density to make sure derivatives
    // wrt. to different segments dont get mixed.
    if (seg != seg_upwind) {
        visc.clearDerivatives();
        density.clearDerivatives();
    }

    const double additional_length = valve.pipeAdditionalLength();
    const double roughness = valve.pipeRoughness();
    const double diameter = valve.pipeDiameter();
    const double area = valve.pipeCrossArea();

    const EvalWell friction_pressure_loss =
        mswellhelpers::frictionPressureLoss(additional_length, diameter, area, roughness, density, mass_rate, visc);

    const double area_con = valve.conCrossArea();
    const double cv = valve.conFlowCoefficient();

    const EvalWell constriction_pressure_loss =
        mswellhelpers::valveContrictionPressureLoss(mass_rate, density, area_con, cv);

    const double sign = mass_rate <= 0. ? 1.0 : -1.0;
    return sign * (friction_pressure_loss + constriction_pressure_loss);
}

template<typename FluidSystem, typename Indices, typename Scalar>
typename MultisegmentWellEval<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
getSegmentSurfaceVolume(const EvalWell& temperature,
                        const EvalWell& saltConcentration,
                        const int pvt_region_index,
                        const int seg_idx) const
{
    const EvalWell seg_pressure = getSegmentPressure(seg_idx);

    std::vector<EvalWell> mix_s(baseif_.numComponents(), 0.0);
    for (int comp_idx = 0; comp_idx < baseif_.numComponents(); ++comp_idx) {
        mix_s[comp_idx] = surfaceVolumeFraction(seg_idx, comp_idx);
    }

    std::vector<EvalWell> b(baseif_.numComponents(), 0.);
    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
        const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
        b[waterCompIdx] =
            FluidSystem::waterPvt().inverseFormationVolumeFactor(pvt_region_index,
                                                                 temperature,
                                                                 seg_pressure,
                                                                 saltConcentration);
    }

    EvalWell rv(0.0);
    EvalWell rvw(0.0);
    // gas phase
    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
        const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
            EvalWell rvmax = FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvt_region_index,
                                                                                  temperature,
                                                                                  seg_pressure);
            if (rvmax < 0.0) { // negative rvmax can happen if the seg_pressure is outside the range of the table
                rvmax = 0.0;
            }
            if (mix_s[oilCompIdx] > 0.0) {
                if (mix_s[gasCompIdx] > 0.0) {
                    rv = mix_s[oilCompIdx] / mix_s[gasCompIdx];
                }

                if (rv > rvmax) {
                    rv = rvmax;
                }
                b[gasCompIdx] =
                    FluidSystem::gasPvt().inverseFormationVolumeFactor(pvt_region_index,
                                                                       temperature,
                                                                       seg_pressure,
                                                                       rv,
                                                                       rvw);
            } else { // no oil exists
                b[gasCompIdx] =
                    FluidSystem::gasPvt().saturatedInverseFormationVolumeFactor(pvt_region_index,
                                                                                temperature,
                                                                                seg_pressure);
            }
        } else { // no Liquid phase
            // it is the same with zero mix_s[Oil]
            b[gasCompIdx] =
                FluidSystem::gasPvt().saturatedInverseFormationVolumeFactor(pvt_region_index,
                                                                            temperature,
                                                                            seg_pressure);
        }
    }

    EvalWell rs(0.0);
    // oil phase
    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
        const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
            EvalWell rsmax = FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvt_region_index,
                                                                                 temperature,
                                                                                 seg_pressure);
            if (rsmax < 0.0) { // negative rsmax can happen if the seg_pressure is outside the range of the table
                rsmax = 0.0;
            }
            if (mix_s[gasCompIdx] > 0.0) {
                if (mix_s[oilCompIdx] > 0.0) {
                    rs = mix_s[gasCompIdx] / mix_s[oilCompIdx];
                }
                // std::cout << " rs " << rs.value() << " rsmax " << rsmax.value() << std::endl;

                if (rs > rsmax) {
                    rs = rsmax;
                }
                b[oilCompIdx] =
                    FluidSystem::oilPvt().inverseFormationVolumeFactor(pvt_region_index,
                                                                       temperature,
                                                                       seg_pressure,
                                                                       rs);
            } else { // no oil exists
                b[oilCompIdx] =
                    FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(pvt_region_index,
                                                                                temperature,
                                                                                seg_pressure);
            }
        } else { // no gas phase
            // it is the same with zero mix_s[Gas]
            b[oilCompIdx] =
                FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(pvt_region_index,
                                                                            temperature,
                                                                            seg_pressure);
        }
    }

    std::vector<EvalWell> mix(mix_s);
    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
        const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
        const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);

        const EvalWell d = 1.0 - rs * rv;
        if (d <= 0.0 || d > 1.0) {
            std::ostringstream sstr;
            sstr << "Problematic d value " << d << " obtained for well " << baseif_.name()
                 << " during conversion to surface volume with rs " << rs
                 << ", rv " << rv << " and pressure " << seg_pressure
                 << " obtaining d " << d
                 << " Continue as if no dissolution (rs = 0) and vaporization (rv = 0) "
                 << " for this connection.";
            OpmLog::debug(sstr.str());
        } else {
            if (rs > 0.0) {
                mix[gasCompIdx] = (mix_s[gasCompIdx] - mix_s[oilCompIdx] * rs) / d;
            }
            if (rv > 0.0) {
                mix[oilCompIdx] = (mix_s[oilCompIdx] - mix_s[gasCompIdx] * rv) / d;
            }
        }
    }

    EvalWell vol_ratio(0.0);
    for (int comp_idx = 0; comp_idx < baseif_.numComponents(); ++comp_idx) {
        vol_ratio += mix[comp_idx] / b[comp_idx];
    }

    // We increase the segment volume with a factor 10 to stabilize the system.
    const double volume = this->segmentSet()[seg_idx].volume();

    return volume / vol_ratio;
}

template<typename FluidSystem, typename Indices, typename Scalar>
void
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
assembleControlEq(const WellState& well_state,
                  const GroupState& group_state,
                  const Schedule& schedule,
                  const SummaryState& summaryState,
                  const Well::InjectionControls& inj_controls,
                  const Well::ProductionControls& prod_controls,
                  const double rho,
                  DeferredLogger& deferred_logger)
{
    static constexpr int Gas = BlackoilPhases::Vapour;
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Water = BlackoilPhases::Aqua;

    EvalWell control_eq(0.0);

    const auto& well = baseif_.wellEcl();

    auto getRates = [&]() {
        std::vector<EvalWell> rates(3, 0.0);
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            rates[Water] = getQs(Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx));
        }
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            rates[Oil] = getQs(Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx));
        }
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            rates[Gas] = getQs(Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx));
        }
        return rates;
    };

    if (baseif_.wellIsStopped()) {
        control_eq = getWQTotal();
    } else if (baseif_.isInjector() ) {
        // Find scaling factor to get injection rate,
        const InjectorType injectorType = inj_controls.injector_type;
        double scaling = 1.0;
        const auto& pu = baseif_.phaseUsage();
        switch (injectorType) {
        case InjectorType::WATER:
        {
            scaling = baseif_.scalingFactor(pu.phase_pos[BlackoilPhases::Aqua]);
            break;
        }
        case InjectorType::OIL:
        {
            scaling = baseif_.scalingFactor(pu.phase_pos[BlackoilPhases::Liquid]);
            break;
        }
        case InjectorType::GAS:
        {
            scaling = baseif_.scalingFactor(pu.phase_pos[BlackoilPhases::Vapour]);
            break;
        }
        default:
            throw("Expected WATER, OIL or GAS as type for injectors " + well.name());
        }
        const EvalWell injection_rate = getWQTotal() / scaling;
        // Setup function for evaluation of BHP from THP (used only if needed).
        auto bhp_from_thp = [&]() {
            const auto rates = getRates();
            return baseif_.calculateBhpFromThp(well_state, rates, well, summaryState, rho, deferred_logger);
        };
        // Call generic implementation.
        baseif_.assembleControlEqInj(well_state,
                                     group_state,
                                     schedule,
                                     summaryState,
                                     inj_controls,
                                     getBhp(),
                                     injection_rate,
                                     bhp_from_thp,
                                     control_eq,
                                     deferred_logger);
    } else {
        // Find rates.
        const auto rates = getRates();
        // Setup function for evaluation of BHP from THP (used only if needed).
        auto bhp_from_thp = [&]() {
            return baseif_.calculateBhpFromThp(well_state, rates, well, summaryState, rho, deferred_logger);
        };
        // Call generic implementation.
        baseif_.assembleControlEqProd(well_state,
                                      group_state,
                                      schedule,
                                      summaryState,
                                      prod_controls,
                                      getBhp(),
                                      rates,
                                      bhp_from_thp,
                                      control_eq,
                                      deferred_logger);
    }

    // using control_eq to update the matrix and residuals
    resWell_[0][SPres] = control_eq.value();
    for (int pv_idx = 0; pv_idx < numWellEq; ++pv_idx) {
        duneD_[0][0][SPres][pv_idx] = control_eq.derivative(pv_idx + Indices::numEq);
    }
}

template<typename FluidSystem, typename Indices, typename Scalar>
void
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
updateThp(WellState& well_state,
          const double rho,
          DeferredLogger& deferred_logger) const
{
    static constexpr int Gas = BlackoilPhases::Vapour;
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Water = BlackoilPhases::Aqua;
    auto& ws = well_state.well(baseif_.indexOfWell());

    // When there is no vaild VFP table provided, we set the thp to be zero.
    if (!baseif_.isVFPActive(deferred_logger) || baseif_.wellIsStopped()) {
        ws.thp = 0;
        return;
    }
    // For THP controlled wells, we know the thp value
    bool thp_controlled = baseif_.isInjector() ? ws.injection_cmode == Well::InjectorCMode::THP:
                                              ws.production_cmode == Well::ProducerCMode::THP;
    if (thp_controlled) {
        return;
    }

    // the well is under other control types, we calculate the thp based on bhp and rates
    std::vector<double> rates(3, 0.0);

    const PhaseUsage& pu = baseif_.phaseUsage();
    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
        rates[ Water ] = ws.surface_rates[pu.phase_pos[ Water ] ];
    }
    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
        rates[ Oil ] = ws.surface_rates[pu.phase_pos[ Oil ] ];
    }
    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
        rates[ Gas ] = ws.surface_rates[pu.phase_pos[ Gas ] ];
    }

    ws.thp = this->calculateThpFromBhp(rates, ws.bhp, rho, deferred_logger);
}

template<typename FluidSystem, typename Indices, typename Scalar>
void
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
handleAccelerationPressureLoss(const int seg,
                               WellState& well_state) const
{
    const double area = this->segmentSet()[seg].crossArea();
    const EvalWell mass_rate = segment_mass_rates_[seg];
    const int seg_upwind = upwinding_segments_[seg];
    EvalWell density = segment_densities_[seg_upwind];
    // WARNING
    // We disregard the derivatives from the upwind density to make sure derivatives
    // wrt. to different segments dont get mixed.
    if (seg != seg_upwind) {
        density.clearDerivatives();
    }

    EvalWell accelerationPressureLoss = mswellhelpers::velocityHead(area, mass_rate, density);
    // handling the velocity head of intlet segments
    for (const int inlet : this->segment_inlets_[seg]) {
        const int seg_upwind_inlet = upwinding_segments_[inlet];
        const double inlet_area = this->segmentSet()[inlet].crossArea();
        EvalWell inlet_density = this->segment_densities_[seg_upwind_inlet];
        // WARNING
        // We disregard the derivatives from the upwind density to make sure derivatives
        // wrt. to different segments dont get mixed.
        if (inlet != seg_upwind_inlet) {
            inlet_density.clearDerivatives();
        }
        const EvalWell inlet_mass_rate = segment_mass_rates_[inlet];
        accelerationPressureLoss -= mswellhelpers::velocityHead(std::max(inlet_area, area), inlet_mass_rate, inlet_density);
    }

    // We change the sign of the accelerationPressureLoss for injectors.
    // Is this correct? Testing indicates that this is what the reference simulator does
    const double sign = mass_rate < 0. ? 1.0 : - 1.0;
    accelerationPressureLoss *= sign;

    auto& segments = well_state.well(baseif_.indexOfWell()).segments;
    segments.pressure_drop_accel[seg] = accelerationPressureLoss.value();

    resWell_[seg][SPres] -= accelerationPressureLoss.value();
    duneD_[seg][seg][SPres][SPres] -= accelerationPressureLoss.derivative(SPres + Indices::numEq);
    duneD_[seg][seg][SPres][WQTotal] -= accelerationPressureLoss.derivative(WQTotal + Indices::numEq);
    if (has_wfrac_variable) {
        duneD_[seg][seg_upwind][SPres][WFrac] -= accelerationPressureLoss.derivative(WFrac + Indices::numEq);
    }
    if (has_gfrac_variable) {
        duneD_[seg][seg_upwind][SPres][GFrac] -= accelerationPressureLoss.derivative(GFrac + Indices::numEq);
    }
}

template<typename FluidSystem, typename Indices, typename Scalar>
void
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
assembleDefaultPressureEq(const int seg,
                          WellState& well_state) const
{
    assert(seg != 0); // not top segment

    // for top segment, the well control equation will be used.
    EvalWell pressure_equation = getSegmentPressure(seg);

    // we need to handle the pressure difference between the two segments
    // we only consider the hydrostatic pressure loss first
    // TODO: we might be able to add member variables to store these values, then we update well state
    // after converged
    const auto hydro_pressure_drop = getHydroPressureLoss(seg);
    auto& ws = well_state.well(baseif_.indexOfWell());
    auto& segments = ws.segments;
    segments.pressure_drop_hydrostatic[seg] = hydro_pressure_drop.value();
    pressure_equation -= hydro_pressure_drop;

    if (this->frictionalPressureLossConsidered()) {
        const auto friction_pressure_drop = getFrictionPressureLoss(seg);
        pressure_equation -= friction_pressure_drop;
        segments.pressure_drop_friction[seg] = friction_pressure_drop.value();
    }

    resWell_[seg][SPres] = pressure_equation.value();
    const int seg_upwind = upwinding_segments_[seg];
    duneD_[seg][seg][SPres][SPres] += pressure_equation.derivative(SPres + Indices::numEq);
    duneD_[seg][seg][SPres][WQTotal] += pressure_equation.derivative(WQTotal + Indices::numEq);
    if (has_wfrac_variable) {
        duneD_[seg][seg_upwind][SPres][WFrac] += pressure_equation.derivative(WFrac + Indices::numEq);
    }
    if (has_gfrac_variable) {
        duneD_[seg][seg_upwind][SPres][GFrac] += pressure_equation.derivative(GFrac + Indices::numEq);
    }

    // contribution from the outlet segment
    const int outlet_segment_index = this->segmentNumberToIndex(this->segmentSet()[seg].outletSegment());
    const EvalWell outlet_pressure = getSegmentPressure(outlet_segment_index);

    resWell_[seg][SPres] -= outlet_pressure.value();
    for (int pv_idx = 0; pv_idx < numWellEq; ++pv_idx) {
        duneD_[seg][outlet_segment_index][SPres][pv_idx] = -outlet_pressure.derivative(pv_idx + Indices::numEq);
    }

    if (this->accelerationalPressureLossConsidered()) {
        handleAccelerationPressureLoss(seg, well_state);
    }
}

template<typename FluidSystem, typename Indices, typename Scalar>
void
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
updateWellStateFromPrimaryVariables(WellState& well_state,
                                    const double rho,
                                    DeferredLogger& deferred_logger) const
{
    static constexpr int Gas = BlackoilPhases::Vapour;
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Water = BlackoilPhases::Aqua;

    const PhaseUsage& pu = baseif_.phaseUsage();
    assert( FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) );
    const int oil_pos = pu.phase_pos[Oil];

    auto& ws = well_state.well(baseif_.indexOfWell());
    auto& segments = ws.segments;
    auto& segment_rates = segments.rates;
    auto& segment_pressure = segments.pressure;
    for (int seg = 0; seg < this->numberOfSegments(); ++seg) {
        std::vector<double> fractions(baseif_.numPhases(), 0.0);
        fractions[oil_pos] = 1.0;

        if ( FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) ) {
            const int water_pos = pu.phase_pos[Water];
            fractions[water_pos] = primary_variables_[seg][WFrac];
            fractions[oil_pos] -= fractions[water_pos];
        }

        if ( FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) ) {
            const int gas_pos = pu.phase_pos[Gas];
            fractions[gas_pos] = primary_variables_[seg][GFrac];
            fractions[oil_pos] -= fractions[gas_pos];
        }

        // convert the fractions to be Q_p / G_total to calculate the phase rates
        for (int p = 0; p < baseif_.numPhases(); ++p) {
            const double scale = baseif_.scalingFactor(p);
            // for injection wells, there should only one non-zero scaling factor
            if (scale > 0.) {
                fractions[p] /= scale;
            } else {
                // this should only happens to injection wells
                fractions[p] = 0.;
            }
        }

        // calculate the phase rates based on the primary variables
        const double g_total = primary_variables_[seg][WQTotal];
        for (int p = 0; p < baseif_.numPhases(); ++p) {
            const double phase_rate = g_total * fractions[p];
            segment_rates[seg*baseif_.numPhases() + p] = phase_rate;
            if (seg == 0) { // top segment
                ws.surface_rates[p] = phase_rate;
            }
        }

        // update the segment pressure
        segment_pressure[seg] = primary_variables_[seg][SPres]*this->bhp_scaling_;
        if (seg == 0) { // top segment
            ws.bhp = segment_pressure[seg]*this->bhp_scaling_;
        }
    }
    updateThp(well_state, rho, deferred_logger);
}

template<typename FluidSystem, typename Indices, typename Scalar>
void
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
assembleICDPressureEq(const int seg,
                      const UnitSystem& unit_system,
                      WellState& well_state,
                      DeferredLogger& deferred_logger) const
{
    // TODO: upwinding needs to be taken care of
    // top segment can not be a spiral ICD device
    assert(seg != 0);

    if (const auto& segment = this->segmentSet()[seg];
       (segment.segmentType() == Segment::SegmentType::VALVE) &&
       (segment.valve().status() == Opm::ICDStatus::SHUT) ) { // we use a zero rate equation to handle SHUT valve
        resWell_[seg][SPres] = this->primary_variables_evaluation_[seg][WQTotal].value();
        duneD_[seg][seg][SPres][WQTotal] = 1.;

        auto& ws = well_state.well(baseif_.indexOfWell());
        ws.segments.pressure_drop_friction[seg] = 0.;
        return;
    }

    // the pressure equation is something like
    // p_seg - deltaP - p_outlet = 0.
    // the major part is how to calculate the deltaP

    EvalWell pressure_equation = getSegmentPressure(seg);

    EvalWell icd_pressure_drop;
    switch(this->segmentSet()[seg].segmentType()) {
        case Segment::SegmentType::SICD :
            icd_pressure_drop = pressureDropSpiralICD(seg);
            break;
        case Segment::SegmentType::AICD :
            icd_pressure_drop = pressureDropAutoICD(seg, unit_system);
            break;
        case Segment::SegmentType::VALVE :
            icd_pressure_drop = pressureDropValve(seg);
            break;
        default: {
            OPM_DEFLOG_THROW(std::runtime_error, "Segment " + std::to_string(this->segmentSet()[seg].segmentNumber())
                             + " for well " + baseif_.name() + " is not of ICD type", deferred_logger);
        }
    }
    pressure_equation = pressure_equation - icd_pressure_drop;
    auto& ws = well_state.well(baseif_.indexOfWell());
    ws.segments.pressure_drop_friction[seg] = icd_pressure_drop.value();

    const int seg_upwind = upwinding_segments_[seg];
    resWell_[seg][SPres] = pressure_equation.value();
    duneD_[seg][seg][SPres][SPres] += pressure_equation.derivative(SPres + Indices::numEq);
    duneD_[seg][seg][SPres][WQTotal] += pressure_equation.derivative(WQTotal + Indices::numEq);
    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
        duneD_[seg][seg_upwind][SPres][WFrac] += pressure_equation.derivative(WFrac + Indices::numEq);
    }
    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
        duneD_[seg][seg_upwind][SPres][GFrac] += pressure_equation.derivative(GFrac + Indices::numEq);
    }

    // contribution from the outlet segment
    const int outlet_segment_index = this->segmentNumberToIndex(this->segmentSet()[seg].outletSegment());
    const EvalWell outlet_pressure = getSegmentPressure(outlet_segment_index);

    resWell_[seg][SPres] -= outlet_pressure.value();
    for (int pv_idx = 0; pv_idx < numWellEq; ++pv_idx) {
        duneD_[seg][outlet_segment_index][SPres][pv_idx] = -outlet_pressure.derivative(pv_idx + Indices::numEq);
    }
}

template<typename FluidSystem, typename Indices, typename Scalar>
void
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
assemblePressureEq(const int seg,
                   const UnitSystem& unit_system,
                   WellState& well_state,
                   DeferredLogger& deferred_logger) const
{
    switch(this->segmentSet()[seg].segmentType()) {
        case Segment::SegmentType::SICD :
        case Segment::SegmentType::AICD :
        case Segment::SegmentType::VALVE : {
            assembleICDPressureEq(seg, unit_system, well_state,deferred_logger);
            break;
        }
        default :
            assembleDefaultPressureEq(seg, well_state);
    }
}

template<typename FluidSystem, typename Indices, typename Scalar>
std::pair<bool, std::vector<Scalar> >
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
getFiniteWellResiduals(const std::vector<Scalar>& B_avg,
                       DeferredLogger& deferred_logger) const
{
    assert(int(B_avg.size() ) == baseif_.numComponents());
    std::vector<Scalar> residuals(numWellEq + 1, 0.0);

    for (int seg = 0; seg < this->numberOfSegments(); ++seg) {
        for (int eq_idx = 0; eq_idx < numWellEq; ++eq_idx) {
            double residual = 0.;
            if (eq_idx < baseif_.numComponents()) {
                residual = std::abs(resWell_[seg][eq_idx]) * B_avg[eq_idx];
            } else {
                if (seg > 0) {
                    residual = std::abs(resWell_[seg][eq_idx]);
                }
            }
            if (std::isnan(residual) || std::isinf(residual)) {
                deferred_logger.debug("nan or inf value for residal get for well " + baseif_.name()
                                      + " segment " + std::to_string(seg) + " eq_idx " + std::to_string(eq_idx));
                return {false, residuals};
            }

            if (residual > residuals[eq_idx]) {
                residuals[eq_idx] = residual;
            }
        }
    }

    // handling the control equation residual
    {
        const double control_residual = std::abs(resWell_[0][numWellEq - 1]);
        if (std::isnan(control_residual) || std::isinf(control_residual)) {
           deferred_logger.debug("nan or inf value for control residal get for well " + baseif_.name());
           return {false, residuals};
        }
        residuals[numWellEq] = control_residual;
    }

    return {true, residuals};
}

template<typename FluidSystem, typename Indices, typename Scalar>
double
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
getControlTolerance(const WellState& well_state,
                    const double tolerance_wells,
                    const double tolerance_pressure_ms_wells,
                    DeferredLogger& deferred_logger) const
{
    double control_tolerance = 0.;

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
            OPM_DEFLOG_THROW(std::runtime_error, "Unknown well control control types for well " << baseif_.name(), deferred_logger);
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
            OPM_DEFLOG_THROW(std::runtime_error, "Unknown well control control types for well " << baseif_.name(), deferred_logger);
        }
    }

    return control_tolerance;
}

template<typename FluidSystem, typename Indices, typename Scalar>
double
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
getResidualMeasureValue(const WellState& well_state,
                        const std::vector<double>& residuals,
                        const double tolerance_wells,
                        const double tolerance_pressure_ms_wells,
                        DeferredLogger& deferred_logger) const
{
    assert(int(residuals.size()) == numWellEq + 1);

    const double rate_tolerance = tolerance_wells;
    int count = 0;
    double sum = 0;
    for (int eq_idx = 0; eq_idx < numWellEq - 1; ++eq_idx) {
        if (residuals[eq_idx] > rate_tolerance) {
            sum += residuals[eq_idx] / rate_tolerance;
            ++count;
        }
    }

    const double pressure_tolerance = tolerance_pressure_ms_wells;
    if (residuals[SPres] > pressure_tolerance) {
        sum += residuals[SPres] / pressure_tolerance;
        ++count;
    }

    const double control_tolerance = getControlTolerance(well_state,
                                                         tolerance_wells,
                                                         tolerance_pressure_ms_wells,
                                                         deferred_logger);
    if (residuals[SPres + 1] > control_tolerance) {
        sum += residuals[SPres + 1] / control_tolerance;
        ++count;
    }

    // if (count == 0), it should be converged.
    assert(count != 0);

    return sum;
}

template<typename FluidSystem, typename Indices, typename Scalar>
void
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
updateUpwindingSegments()
{
    for (int seg = 0; seg < this->numberOfSegments(); ++seg) {
        // special treatment is needed for segment 0
        if (seg == 0) {
            // we are not supposed to have injecting producers and producing injectors
            assert( ! (baseif_.isProducer() && primary_variables_evaluation_[seg][WQTotal] > 0.) );
            assert( ! (baseif_.isInjector() && primary_variables_evaluation_[seg][WQTotal] < 0.) );
            upwinding_segments_[seg] = seg;
            continue;
        }

        // for other normal segments
        if (primary_variables_evaluation_[seg][WQTotal] <= 0.) {
            upwinding_segments_[seg] = seg;
        } else {
            const int outlet_segment_index = this->segmentNumberToIndex(this->segmentSet()[seg].outletSegment());
            upwinding_segments_[seg] = outlet_segment_index;
        }
    }
}

template<typename FluidSystem, typename Indices, typename Scalar>
void
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
addWellContribution(WellContributions& wellContribs) const
{
    unsigned int Mb = duneB_.N();       // number of blockrows in duneB_, duneC_ and duneD_
    unsigned int BnumBlocks = duneB_.nonzeroes();
    unsigned int DnumBlocks = duneD_.nonzeroes();

    // duneC
    std::vector<unsigned int> Ccols;
    std::vector<double> Cvals;
    Ccols.reserve(BnumBlocks);
    Cvals.reserve(BnumBlocks * Indices::numEq * numWellEq);
    for (auto rowC = duneC_.begin(); rowC != duneC_.end(); ++rowC) {
        for (auto colC = rowC->begin(), endC = rowC->end(); colC != endC; ++colC) {
            Ccols.emplace_back(colC.index());
            for (int i = 0; i < numWellEq; ++i) {
                for (int j = 0; j < Indices::numEq; ++j) {
                    Cvals.emplace_back((*colC)[i][j]);
                }
            }
        }
    }

    // duneD
    Dune::UMFPack<DiagMatWell> umfpackMatrix(duneD_, 0);
    double *Dvals = umfpackMatrix.getInternalMatrix().getValues();
    auto *Dcols = umfpackMatrix.getInternalMatrix().getColStart();
    auto *Drows = umfpackMatrix.getInternalMatrix().getRowIndex();

    // duneB
    std::vector<unsigned int> Bcols;
    std::vector<unsigned int> Brows;
    std::vector<double> Bvals;
    Bcols.reserve(BnumBlocks);
    Brows.reserve(Mb+1);
    Bvals.reserve(BnumBlocks * Indices::numEq * numWellEq);
    Brows.emplace_back(0);
    unsigned int sumBlocks = 0;
    for (auto rowB = duneB_.begin(); rowB != duneB_.end(); ++rowB) {
        int sizeRow = 0;
        for (auto colB = rowB->begin(), endB = rowB->end(); colB != endB; ++colB) {
            Bcols.emplace_back(colB.index());
            for (int i = 0; i < numWellEq; ++i) {
                for (int j = 0; j < Indices::numEq; ++j) {
                    Bvals.emplace_back((*colB)[i][j]);
                }
            }
            sizeRow++;
        }
        sumBlocks += sizeRow;
        Brows.emplace_back(sumBlocks);
    }

    wellContribs.addMultisegmentWellContribution(Indices::numEq,
                                                 numWellEq,
                                                 Mb,
                                                 Bvals,
                                                 Bcols,
                                                 Brows,
                                                 DnumBlocks,
                                                 Dvals,
                                                 Dcols,
                                                 Drows,
                                                 Cvals);
}

#define INSTANCE(A,...) \
template class MultisegmentWellEval<BlackOilFluidSystem<double,A>,__VA_ARGS__,double>;

// One phase
INSTANCE(BlackOilDefaultIndexTraits,BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(BlackOilDefaultIndexTraits,BlackOilOnePhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(BlackOilDefaultIndexTraits,BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,5u>)

// Two phase
INSTANCE(BlackOilDefaultIndexTraits,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,0u,0u>)
INSTANCE(BlackOilDefaultIndexTraits,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(BlackOilDefaultIndexTraits,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,2u,0u>)
INSTANCE(BlackOilDefaultIndexTraits,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,2u,0u>)
INSTANCE(BlackOilDefaultIndexTraits,BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,false,0u,2u,0u>)
INSTANCE(BlackOilDefaultIndexTraits,BlackOilTwoPhaseIndices<0u,0u,2u,0u,false,false,0u,2u,0u>)
INSTANCE(BlackOilDefaultIndexTraits,BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(BlackOilDefaultIndexTraits,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,0u,0u>)

// Blackoil
INSTANCE(BlackOilDefaultIndexTraits,BlackOilIndices<0u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilDefaultIndexTraits,BlackOilIndices<0u,0u,0u,0u,true,false,0u,0u>)
INSTANCE(BlackOilDefaultIndexTraits,BlackOilIndices<0u,0u,0u,0u,false,true,0u,0u>)
INSTANCE(BlackOilDefaultIndexTraits,BlackOilIndices<0u,0u,0u,0u,false,true,2u,0u>)
INSTANCE(BlackOilDefaultIndexTraits,BlackOilIndices<1u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilDefaultIndexTraits,BlackOilIndices<0u,1u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilDefaultIndexTraits,BlackOilIndices<0u,0u,1u,0u,false,false,0u,0u>)
INSTANCE(BlackOilDefaultIndexTraits,BlackOilIndices<0u,0u,0u,1u,false,false,0u,0u>)
INSTANCE(BlackOilDefaultIndexTraits,BlackOilIndices<0u,0u,0u,0u,false,false,1u,0u>)
INSTANCE(BlackOilDefaultIndexTraits,BlackOilIndices<0u,0u,0u,1u,false,true,0u,0u>)

} // namespace Opm
