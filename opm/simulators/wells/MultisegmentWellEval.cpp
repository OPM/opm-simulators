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

#include <opm/simulators/timestepping/ConvergenceReport.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/MSWellHelpers.hpp>
#include <opm/simulators/wells/MultisegmentWellAssemble.hpp>
#include <opm/simulators/wells/RateConverter.hpp>
#include <opm/simulators/wells/WellAssemble.hpp>
#include <opm/simulators/wells/WellBhpThpCalculator.hpp>
#include <opm/simulators/wells/WellConvergence.hpp>
#include <opm/simulators/wells/WellInterfaceIndices.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <numeric>
#include <utility>
#include <vector>

namespace Opm
{

template<typename FluidSystem, typename Indices, typename Scalar>
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
MultisegmentWellEval(WellInterfaceIndices<FluidSystem,Indices,Scalar>& baseif)
    : MultisegmentWellGeneric<Scalar>(baseif)
    , baseif_(baseif)
    , linSys_(*this)
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
initMatrixAndVectors(const int num_cells)
{
    linSys_.init(num_cells, baseif_.numPerfs(), baseif_.cells());
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
            abs_residual[seg][eq_idx] = std::abs(linSys_.residual()[seg][eq_idx]);
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

    WellConvergence(baseif_).
        checkConvergenceControlEq(well_state,
                                  {tolerance_pressure_ms_wells,
                                   tolerance_pressure_ms_wells,
                                   tolerance_wells,
                                   tolerance_wells,
                                   max_residual_allowed},
                                  std::abs(linSys_.residual()[0][SPres]),
                                  report,
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
            const double dx_limited = sign * std::min(std::abs(dwells[seg][SPres]) * relaxation_factor, max_pressure_change);
            primary_variables_[seg][SPres] = std::max( old_primary_variables[seg][SPres] - dx_limited, 1e5);
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
    // maybe a throw for parallel running?
    assert(int(segments.size()) == this->numberOfSegments());
    const auto& segment_rates = segments.rates;
    const auto& segment_pressure = segments.pressure;
    const PhaseUsage& pu = baseif_.phaseUsage();

    for (int seg = 0; seg < this->numberOfSegments(); ++seg) {
        // calculate the total rate for each segment
        double total_seg_rate = 0.0;
        // the segment pressure
        primary_variables_[seg][SPres] = segment_pressure[seg];
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
        primary_variables_[seg][WQTotal] = total_seg_rate;
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
            return primary_variables_evaluation_[seg][WQTotal] / baseif_.scalingFactor(baseif_.ebosCompIdxToFlowCompIdx(comp_idx));


        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)
                && Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx) == comp_idx
                && phase == InjectorType::OIL)
            return primary_variables_evaluation_[seg][WQTotal] / baseif_.scalingFactor(baseif_.ebosCompIdxToFlowCompIdx(comp_idx));

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)
                && Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx) == comp_idx
                && phase == InjectorType::GAS)
            return primary_variables_evaluation_[seg][WQTotal] / baseif_.scalingFactor(baseif_.ebosCompIdxToFlowCompIdx(comp_idx));

        return 0.0;
    }

    const EvalWell segment_rate = primary_variables_evaluation_[seg][WQTotal] * volumeFractionScaled(seg_upwind, comp_idx);

    assert(segment_rate.derivative(SPres + Indices::numEq) == 0.);

    return segment_rate;
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
    return primary_variables_evaluation_[seg][SPres];
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
    return primary_variables_evaluation_[seg][WQTotal] * volumeFractionScaled(seg, comp_idx);
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


    const EvalWell liquid_fraction = water_fraction + oil_fraction;

    // viscosity contribution from the liquid
    const EvalWell liquid_viscosity_fraction = liquid_fraction < 1.e-30 ? oil_fraction * oil_viscosity + water_fraction * water_viscosity :
            liquid_fraction * mswellhelpers::emulsionViscosity(water_fraction, water_viscosity, oil_fraction, oil_viscosity, sicd);

    const EvalWell mixture_viscosity = liquid_viscosity_fraction + gas_fraction * gas_viscosity;

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
handleAccelerationPressureLoss(const int seg,
                               WellState& well_state)
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

    MultisegmentWellAssemble<FluidSystem,Indices,Scalar>(baseif_).
        assemblePressureLoss(seg, seg_upwind, accelerationPressureLoss, linSys_);
}

template<typename FluidSystem, typename Indices, typename Scalar>
void
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
assembleDefaultPressureEq(const int seg,
                          WellState& well_state)
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

    // contribution from the outlet segment
    const int outlet_segment_index = this->segmentNumberToIndex(this->segmentSet()[seg].outletSegment());
    const EvalWell outlet_pressure = getSegmentPressure(outlet_segment_index);

    const int seg_upwind = upwinding_segments_[seg];
    MultisegmentWellAssemble<FluidSystem,Indices,Scalar>(baseif_).
        assemblePressureEq(seg, seg_upwind, outlet_segment_index,
                           pressure_equation, outlet_pressure, linSys_);

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

    const auto pvtReg = std::max(this->baseif_.wellEcl().pvt_table_number() - 1, 0);

    const PhaseUsage& pu = baseif_.phaseUsage();
    assert( FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) );
    const int oil_pos = pu.phase_pos[Oil];

    auto& ws = well_state.well(baseif_.indexOfWell());
    auto& segments = ws.segments;
    auto& segment_rates = segments.rates;
    auto& disgas = segments.dissolved_gas_rate;
    auto& vapoil = segments.vaporized_oil_rate;
    auto& segment_pressure = segments.pressure;
    for (int seg = 0; seg < this->numberOfSegments(); ++seg) {
        std::vector<double> fractions(baseif_.numPhases(), 0.0);
        fractions[oil_pos] = 1.0;

        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            const int water_pos = pu.phase_pos[Water];
            fractions[water_pos] = primary_variables_[seg][WFrac];
            fractions[oil_pos] -= fractions[water_pos];
        }

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
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
        segment_pressure[seg] = primary_variables_[seg][SPres];

        if (seg == 0) { // top segment
            ws.bhp = segment_pressure[seg];
        }

        // Calculate other per-phase dynamic quantities.

        const auto temperature = 0.0; // Ignore thermal effects
        const auto saltConc = 0.0;    // Ignore salt precipitation
        const auto Rvw = 0.0;         // Ignore vaporised water.

        auto rsMax = 0.0;
        auto rvMax = 0.0;
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            // Both oil and gas active.
            rsMax = FluidSystem::oilPvt()
                .saturatedGasDissolutionFactor(pvtReg, temperature, segment_pressure[seg]);

            rvMax = FluidSystem::gasPvt()
                .saturatedOilVaporizationFactor(pvtReg, temperature, segment_pressure[seg]);
        }

        // 1) Infer phase splitting for oil/gas.
        const auto& [Rs, Rv] = this->baseif_.rateConverter().inferDissolvedVaporisedRatio
            (rsMax, rvMax, segment_rates.begin() + (seg + 0)*this->baseif_.numPhases());

        if (! FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            vapoil[seg] = disgas[seg] = 0.0;
        }
        else {
            const auto* qs = &segment_rates[seg*this->baseif_.numPhases() + 0];
            const auto denom = 1.0 - (Rs * Rv);
            const auto io = pu.phase_pos[Oil];
            const auto ig = pu.phase_pos[Gas];
            disgas[seg] = Rs * (qs[io] - Rv*qs[ig]) / denom;
            vapoil[seg] = Rv * (qs[ig] - Rs*qs[io]) / denom;
        }

        // 2) Local condition volume flow rates
        {
            // Use std::span<> in C++20 and beyond.
            const auto  rate_start = (seg + 0) * this->baseif_.numPhases();
            const auto* surf_rates = segment_rates.data()             + rate_start;
            auto*       resv_rates = segments.phase_resv_rates.data() + rate_start;

            this->baseif_.rateConverter().calcReservoirVoidageRates
                (pvtReg, segment_pressure[seg],
                 std::max(0.0, Rs),
                 std::max(0.0, Rv),
                 temperature, saltConc, surf_rates, resv_rates);
        }

        // 3) Local condition holdup fractions.
        const auto tot_resv =
            std::accumulate(segments.phase_resv_rates.begin() + (seg + 0)*this->baseif_.numPhases(),
                            segments.phase_resv_rates.begin() + (seg + 1)*this->baseif_.numPhases(),
                            0.0);

        std::transform(segments.phase_resv_rates.begin() + (seg + 0)*this->baseif_.numPhases(),
                       segments.phase_resv_rates.begin() + (seg + 1)*this->baseif_.numPhases(),
                       segments.phase_holdup.begin()     + (seg + 0)*this->baseif_.numPhases(),
                       [tot_resv](const auto qr) { return std::clamp(qr / tot_resv, 0.0, 1.0); });

        // 4) Local condition flow velocities for segments other than top segment.
        if (seg > 0) {
            // Possibly poor approximation
            //    Velocity = Flow rate / cross-sectional area.
            // Additionally ignores drift flux.
            const auto area = this->baseif_.wellEcl().getSegments()
                .getFromSegmentNumber(segments.segment_number()[seg]).crossArea();
            const auto velocity = (area > 0.0) ? tot_resv / area : 0.0;

            std::transform(segments.phase_holdup.begin()   + (seg + 0)*this->baseif_.numPhases(),
                           segments.phase_holdup.begin()   + (seg + 1)*this->baseif_.numPhases(),
                           segments.phase_velocity.begin() + (seg + 0)*this->baseif_.numPhases(),
                           [velocity](const auto hf) { return (hf > 0.0) ? velocity : 0.0; });
        }

        // 5) Local condition phase viscosities.
        segments.phase_viscosity[seg*this->baseif_.numPhases() + pu.phase_pos[Oil]] =
            FluidSystem::oilPvt().viscosity(pvtReg, temperature, segment_pressure[seg], Rs);

        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            segments.phase_viscosity[seg*this->baseif_.numPhases() + pu.phase_pos[Water]] =
                FluidSystem::waterPvt().viscosity(pvtReg, temperature, segment_pressure[seg], saltConc);
        }

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            segments.phase_viscosity[seg*this->baseif_.numPhases() + pu.phase_pos[Gas]] =
                FluidSystem::gasPvt().viscosity(pvtReg, temperature, segment_pressure[seg], Rv, Rvw);
        }
    }

    // Segment flow velocity in top segment.
    {
        const auto np = this->baseif_.numPhases();
        auto segVel = [&segments, np](const auto segmentNumber)
        {
            auto v = 0.0;
            const auto* vel = segments.phase_velocity.data() + segmentNumber*np;
            for (auto p = 0*np; p < np; ++p) {
                if (std::abs(vel[p]) > std::abs(v)) {
                    v = vel[p];
                }
            }

            return v;
        };

        const auto seg = 0;
        auto maxVel = 0.0;
        for (const auto& inlet : this->segmentSet()[seg].inletSegments()) {
            const auto v = segVel(this->segmentNumberToIndex(inlet));
            if (std::abs(v) > std::abs(maxVel)) {
                maxVel = v;
            }
        }

        std::transform(segments.phase_holdup.begin()   + (seg + 0)*this->baseif_.numPhases(),
                       segments.phase_holdup.begin()   + (seg + 1)*this->baseif_.numPhases(),
                       segments.phase_velocity.begin() + (seg + 0)*this->baseif_.numPhases(),
                       [maxVel](const auto hf) { return (hf > 0.0) ? maxVel : 0.0; });
    }

    WellBhpThpCalculator(this->baseif_)
        .updateThp(rho, [this]() { return this->baseif_.wellEcl().alq_value(); },
                   {FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx),
                    FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx),
                    FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)},
                   well_state, deferred_logger);
}

template<typename FluidSystem, typename Indices, typename Scalar>
void
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
assembleICDPressureEq(const int seg,
                      const UnitSystem& unit_system,
                      WellState& well_state,
                      DeferredLogger& deferred_logger)
{
    // TODO: upwinding needs to be taken care of
    // top segment can not be a spiral ICD device
    assert(seg != 0);

    if (const auto& segment = this->segmentSet()[seg];
       (segment.segmentType() == Segment::SegmentType::VALVE) &&
       (segment.valve().status() == Opm::ICDStatus::SHUT) ) { // we use a zero rate equation to handle SHUT valve
        MultisegmentWellAssemble<FluidSystem,Indices,Scalar>(baseif_).
            assembleTrivialEq(seg, this->primary_variables_evaluation_[seg][WQTotal].value(), linSys_);

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

    // contribution from the outlet segment
    const int outlet_segment_index = this->segmentNumberToIndex(this->segmentSet()[seg].outletSegment());
    const EvalWell outlet_pressure = getSegmentPressure(outlet_segment_index);

    const int seg_upwind = upwinding_segments_[seg];
    MultisegmentWellAssemble<FluidSystem,Indices,Scalar>(baseif_).
        assemblePressureEq(seg, seg_upwind, outlet_segment_index,
                           pressure_equation, outlet_pressure,
                           linSys_,
                           FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx),
                           FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx));
}

template<typename FluidSystem, typename Indices, typename Scalar>
void
MultisegmentWellEval<FluidSystem,Indices,Scalar>::
assemblePressureEq(const int seg,
                   const UnitSystem& unit_system,
                   WellState& well_state,
                   DeferredLogger& deferred_logger)
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
                residual = std::abs(linSys_.residual()[seg][eq_idx]) * B_avg[eq_idx];
            } else {
                if (seg > 0) {
                    residual = std::abs(linSys_.residual()[seg][eq_idx]);
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
        const double control_residual = std::abs(linSys_.residual()[0][numWellEq - 1]);
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

#define INSTANCE(...) \
template class MultisegmentWellEval<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,__VA_ARGS__,double>;

// One phase
INSTANCE(BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(BlackOilOnePhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,5u>)

// Two phase
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,0u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,false,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,2u,0u,false,false,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,0u,0u>)

// Blackoil
INSTANCE(BlackOilIndices<0u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,0u,true,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,0u,false,true,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,0u,false,true,2u,0u>)
INSTANCE(BlackOilIndices<1u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,1u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,1u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,1u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,0u,false,false,1u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,1u,false,true,0u,0u>)

} // namespace Opm
