/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2016 - 2017 IRIS AS.

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
#include <opm/simulators/wells/MultisegmentWellPrimaryVariables.hpp>

#include <opm/input/eclipse/Schedule/MSW/WellSegments.hpp>
#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/material/fluidsystems/BlackOilDefaultIndexTraits.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/models/blackoil/blackoilindices.hh>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/simulators/wells/MultisegmentWellGeneric.hpp>
#include <opm/simulators/wells/RateConverter.hpp>
#include <opm/simulators/wells/WellBhpThpCalculator.hpp>
#include <opm/simulators/wells/WellInterfaceIndices.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <algorithm>

namespace Opm {

template<class FluidSystem, class Indices, class Scalar>
void MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::
resize(const int numSegments)
{
    value_.resize(numSegments);
    evaluation_.resize(numSegments);
}

template<class FluidSystem, class Indices, class Scalar>
void MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::
init()
{
    for (std::size_t seg = 0; seg < value_.size(); ++seg) {
        for (int eq_idx = 0; eq_idx < numWellEq; ++eq_idx) {
            evaluation_[seg][eq_idx] = 0.0;
            evaluation_[seg][eq_idx].setValue(value_[seg][eq_idx]);
            evaluation_[seg][eq_idx].setDerivative(eq_idx + Indices::numEq, 1.0);
        }
    }
}

template<class FluidSystem, class Indices, class Scalar>
void MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::
update(const WellState& well_state, const bool stop_or_zero_rate_target)
{
    static constexpr int Water = BlackoilPhases::Aqua;
    static constexpr int Gas = BlackoilPhases::Vapour;

    // TODO: to test using rate conversion coefficients to see if it will be better than
    // this default one
    if (!well_.isOperableAndSolvable() && !well_.wellIsStopped())
        return;

    const Well& well = well_.wellEcl();

    // the index of the top segment in the WellState
    const auto& ws = well_state.well(well_.indexOfWell());
    const auto& segments = ws.segments;
    // maybe a throw for parallel running?
    assert(segments.size() == value_.size());
    const auto& segment_rates = segments.rates;
    const auto& segment_pressure = segments.pressure;
    const PhaseUsage& pu = well_.phaseUsage();

    for (std::size_t seg = 0; seg < value_.size(); ++seg) {
        // calculate the total rate for each segment
        double total_seg_rate = 0.0;
        // the segment pressure
        value_[seg][SPres] = segment_pressure[seg];
        // TODO: under what kind of circustances, the following will be wrong?
        // the definition of g makes the gas phase is always the last phase
        for (int p = 0; p < well_.numPhases(); p++) {
            total_seg_rate += well_.scalingFactor(p) * segment_rates[well_.numPhases() * seg + p];
        }

        if (seg == 0) {
            if (well_.isInjector()) {
                total_seg_rate = std::max(total_seg_rate, 0.);
            } else {
                total_seg_rate = std::min(total_seg_rate, 0.);
            }
        }
        value_[seg][WQTotal] = total_seg_rate;
        if (stop_or_zero_rate_target && seg == 0) {
            value_[seg][WQTotal] = 0;
        }
        if (std::abs(total_seg_rate) > 0.) {
            if (has_wfrac_variable) {
                const int water_pos = pu.phase_pos[Water];
                value_[seg][WFrac] = well_.scalingFactor(water_pos) * segment_rates[well_.numPhases() * seg + water_pos] / total_seg_rate;
            }
            if (has_gfrac_variable) {
                const int gas_pos = pu.phase_pos[Gas];
                value_[seg][GFrac] = well_.scalingFactor(gas_pos) * segment_rates[well_.numPhases() * seg + gas_pos] / total_seg_rate;
            }
        } else { // total_seg_rate == 0
            if (well_.isInjector()) {
                // only single phase injection handled
                auto phase = well.getInjectionProperties().injectorType;

                if (has_wfrac_variable) {
                    if (phase == InjectorType::WATER) {
                        value_[seg][WFrac] = 1.0;
                    } else {
                        value_[seg][WFrac] = 0.0;
                    }
                }

                if (has_gfrac_variable) {
                    if (phase == InjectorType::GAS) {
                        value_[seg][GFrac] = 1.0;
                    } else {
                        value_[seg][GFrac] = 0.0;
                    }
                }

            } else if (well_.isProducer()) { // producers
                if (has_wfrac_variable) {
                    value_[seg][WFrac] = 1.0 / well_.numPhases();
                }

                if (has_gfrac_variable) {
                    value_[seg][GFrac] = 1.0 / well_.numPhases();
                }
            }
        }
    }
}

template<class FluidSystem, class Indices, class Scalar>
void MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::
updateNewton(const BVectorWell& dwells,
             const double relaxation_factor,
             const double dFLimit,
             const bool stop_or_zero_rate_target,
             const double max_pressure_change)
{
    const std::vector<std::array<double, numWellEq>> old_primary_variables = value_;

    for (std::size_t seg = 0; seg < value_.size(); ++seg) {
        if (has_wfrac_variable) {
            const int sign = dwells[seg][WFrac] > 0. ? 1 : -1;
            const double dx_limited = sign * std::min(std::abs(dwells[seg][WFrac]) * relaxation_factor, dFLimit);
            value_[seg][WFrac] = old_primary_variables[seg][WFrac] - dx_limited;
        }

        if (has_gfrac_variable) {
            const int sign = dwells[seg][GFrac] > 0. ? 1 : -1;
            const double dx_limited = sign * std::min(std::abs(dwells[seg][GFrac]) * relaxation_factor, dFLimit);
            value_[seg][GFrac] = old_primary_variables[seg][GFrac] - dx_limited;
        }

        // handling the overshooting or undershooting of the fractions
        this->processFractions(seg);

        // update the segment pressure
        {
            const int sign = dwells[seg][SPres] > 0.? 1 : -1;
            const double dx_limited = sign * std::min(std::abs(dwells[seg][SPres]) * relaxation_factor, max_pressure_change);
            // some cases might have defaulted bhp constraint of 1 bar, we use a slightly smaller value as the bhp lower limit for Newton update
            // so that bhp constaint can be an active control when needed.
            // TODO: we might need to distinguish SPres of the top segment from the SPres of the other segments, since only the top segment
            // is involved in the well constraints/control checking
            constexpr double spres_lower_limit = 1. * unit::barsa - 1. * unit::Pascal;
            value_[seg][SPres] = std::max(old_primary_variables[seg][SPres] - dx_limited, spres_lower_limit);
        }

        // update the total rate // TODO: should we have a limitation of the total rate change?
        {
            value_[seg][WQTotal] = old_primary_variables[seg][WQTotal] - relaxation_factor * dwells[seg][WQTotal];

            // make sure that no injector produce and no producer inject
            if (seg == 0) {
                if (well_.isInjector()) {
                    value_[seg][WQTotal] = std::max(value_[seg][WQTotal], 0.0);
                } else {
                    value_[seg][WQTotal] = std::min(value_[seg][WQTotal], 0.0);
                }
            }
        }
    }

    if (stop_or_zero_rate_target) {
        value_[0][WQTotal] = 0.;
    }
}

template<class FluidSystem, class Indices, class Scalar>
void MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::
copyToWellState(const MultisegmentWellGeneric<Scalar>& mswell,
                const double rho,
                const bool stop_or_zero_rate_target,
                WellState& well_state,
                const SummaryState& summary_state,
                DeferredLogger& deferred_logger) const
{
    static constexpr int Gas = BlackoilPhases::Vapour;
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Water = BlackoilPhases::Aqua;

    const auto pvtReg = std::max(well_.wellEcl().pvt_table_number() - 1, 0);

    const PhaseUsage& pu = well_.phaseUsage();
    assert( FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) );
    const int oil_pos = pu.phase_pos[Oil];

    auto& ws = well_state.well(well_.indexOfWell());
    auto& segments = ws.segments;
    auto& segment_rates = segments.rates;
    auto& disgas = segments.dissolved_gas_rate;
    auto& vapoil = segments.vaporized_oil_rate;
    auto& segment_pressure = segments.pressure;
    for (std::size_t seg = 0; seg < value_.size(); ++seg) {
        std::vector<double> fractions(well_.numPhases(), 0.0);
        fractions[oil_pos] = 1.0;

        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            const int water_pos = pu.phase_pos[Water];
            fractions[water_pos] = value_[seg][WFrac];
            fractions[oil_pos] -= fractions[water_pos];
        }

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            const int gas_pos = pu.phase_pos[Gas];
            fractions[gas_pos] = value_[seg][GFrac];
            fractions[oil_pos] -= fractions[gas_pos];
        }

        // convert the fractions to be Q_p / G_total to calculate the phase rates
        for (int p = 0; p < well_.numPhases(); ++p) {
            const double scale = well_.scalingFactor(p);
            // for injection wells, there should only one non-zero scaling factor
            if (scale > 0.) {
                fractions[p] /= scale;
            } else {
                // this should only happens to injection wells
                fractions[p] = 0.;
            }
        }

        // calculate the phase rates based on the primary variables
        const double g_total = value_[seg][WQTotal];
        for (int p = 0; p < well_.numPhases(); ++p) {
            const double phase_rate = g_total * fractions[p];
            segment_rates[seg * well_.numPhases() + p] = phase_rate;
            if (seg == 0) { // top segment
                ws.surface_rates[p] = phase_rate;
            }
        }

        // update the segment pressure
        segment_pressure[seg] = value_[seg][SPres];

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
        const auto& [Rs, Rv] = well_.rateConverter().inferDissolvedVaporisedRatio
            (rsMax, rvMax, segment_rates.begin() + seg * well_.numPhases());

        if (! FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            vapoil[seg] = disgas[seg] = 0.0;
        }
        else {
            const auto* qs = &segment_rates[seg * well_.numPhases()];
            const auto denom = 1.0 - (Rs * Rv);
            const auto io = pu.phase_pos[Oil];
            const auto ig = pu.phase_pos[Gas];
            disgas[seg] = Rs * (qs[io] - Rv*qs[ig]) / denom;
            vapoil[seg] = Rv * (qs[ig] - Rs*qs[io]) / denom;
        }

        // 2) Local condition volume flow rates
        {
            // Use std::span<> in C++20 and beyond.
            const auto  rate_start = seg  * well_.numPhases();
            const auto* surf_rates = segment_rates.data()             + rate_start;
            auto*       resv_rates = segments.phase_resv_rates.data() + rate_start;

            well_.rateConverter().calcReservoirVoidageRates
                (pvtReg, segment_pressure[seg],
                 std::max(0.0, Rs),
                 std::max(0.0, Rv),
                 0.0, // Rsw
                 0.0, // Rvw
                 temperature, saltConc, surf_rates, resv_rates);
        }

        // 3) Local condition holdup fractions.
        const auto tot_resv =
            std::accumulate(segments.phase_resv_rates.begin() + (seg + 0) * well_.numPhases(),
                            segments.phase_resv_rates.begin() + (seg + 1) * well_.numPhases(),
                            0.0);

        std::transform(segments.phase_resv_rates.begin() + (seg + 0) * well_.numPhases(),
                       segments.phase_resv_rates.begin() + (seg + 1) * well_.numPhases(),
                       segments.phase_holdup.begin()     + (seg + 0) * well_.numPhases(),
                       [tot_resv](const auto qr) { return std::clamp(qr / tot_resv, 0.0, 1.0); });

        // 4) Local condition flow velocities for segments other than top segment.
        if (seg > 0) {
            // Possibly poor approximation
            //    Velocity = Flow rate / cross-sectional area.
            // Additionally ignores drift flux.
            const auto area = well_.wellEcl().getSegments()
                .getFromSegmentNumber(segments.segment_number()[seg]).crossArea();
            const auto velocity = (area > 0.0) ? tot_resv / area : 0.0;

            std::transform(segments.phase_holdup.begin()   + (seg + 0) * well_.numPhases(),
                           segments.phase_holdup.begin()   + (seg + 1) * well_.numPhases(),
                           segments.phase_velocity.begin() + (seg + 0) * well_.numPhases(),
                           [velocity](const auto hf) { return (hf > 0.0) ? velocity : 0.0; });
        }

        // 5) Local condition phase viscosities.
        segments.phase_viscosity[seg * well_.numPhases() + pu.phase_pos[Oil]] =
            FluidSystem::oilPvt().viscosity(pvtReg, temperature, segment_pressure[seg], Rs);

        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            segments.phase_viscosity[seg * well_.numPhases() + pu.phase_pos[Water]] =
                FluidSystem::waterPvt().viscosity(pvtReg, temperature, segment_pressure[seg],  0.0 /*Rsw*/, saltConc);
        }

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            segments.phase_viscosity[seg * well_.numPhases() + pu.phase_pos[Gas]] =
                FluidSystem::gasPvt().viscosity(pvtReg, temperature, segment_pressure[seg], Rv, Rvw);
        }
    }

    // Segment flow velocity in top segment.
    {
        const auto np = well_.numPhases();
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
        for (const auto& inlet : mswell.segmentSet()[seg].inletSegments()) {
            const auto v = segVel(mswell.segmentNumberToIndex(inlet));
            if (std::abs(v) > std::abs(maxVel)) {
                maxVel = v;
            }
        }

        std::transform(segments.phase_holdup.begin()   + (seg + 0) * well_.numPhases(),
                       segments.phase_holdup.begin()   + (seg + 1) * well_.numPhases(),
                       segments.phase_velocity.begin() + (seg + 0) * well_.numPhases(),
                       [maxVel](const auto hf) { return (hf > 0.0) ? maxVel : 0.0; });
    }

    // Note: for the ALQ value, in the StandardWell, WellInterfaceGeneric::getALQ(well_state) is used.
    // We might want to unify the way regarding AQL value.
    WellBhpThpCalculator(well_)
        .updateThp(rho, stop_or_zero_rate_target, [this, summary_state]() { return well_.wellEcl().alq_value(summary_state); },
                   {FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx),
                    FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx),
                    FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)},
                   well_state, summary_state, deferred_logger);
}

template<class FluidSystem, class Indices, class Scalar>
void MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::
processFractions(const int seg)
{
    static constexpr int Water = BlackoilPhases::Aqua;
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Gas = BlackoilPhases::Vapour;

    const PhaseUsage& pu = well_.phaseUsage();

    std::vector<double> fractions(well_.numPhases(), 0.0);

    assert(FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx));
    const int oil_pos = pu.phase_pos[Oil];
    fractions[oil_pos] = 1.0;

    if ( FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) ) {
        const int water_pos = pu.phase_pos[Water];
        fractions[water_pos] = value_[seg][WFrac];
        fractions[oil_pos] -= fractions[water_pos];
    }

    if ( FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) ) {
        const int gas_pos = pu.phase_pos[Gas];
        fractions[gas_pos] = value_[seg][GFrac];
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
        value_[seg][WFrac] = fractions[pu.phase_pos[Water]];
    }

    if ( FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) ) {
        value_[seg][GFrac] = fractions[pu.phase_pos[Gas]];
    }
}

template<typename FluidSystem, typename Indices, typename Scalar>
typename MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::
volumeFraction(const int seg,
               const unsigned compIdx) const
{
    if (has_wfrac_variable && compIdx == Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx)) {
        return evaluation_[seg][WFrac];
    }

    if (has_gfrac_variable && compIdx == Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx)) {
        return evaluation_[seg][GFrac];
    }

    // Oil fraction
    EvalWell oil_fraction = 1.0;
    if (has_wfrac_variable) {
        oil_fraction -= evaluation_[seg][WFrac];
    }

    if (has_gfrac_variable) {
        oil_fraction -= evaluation_[seg][GFrac];
    }
    /* if (has_solvent) {
        oil_fraction -= evaluation_[seg][SFrac];
    } */

    // oil_fraction may turn out negative due to round-off, in that case
    // set to zero (but keep derivatives)
    if (oil_fraction.value() < 0.0) {
        oil_fraction.setValue(0.0);
    }
    return oil_fraction;
}

template<class FluidSystem, class Indices, class Scalar>
typename MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::
volumeFractionScaled(const int seg,
                     const int comp_idx) const
{
    // For reservoir rate control, the distr in well control is used for the
    // rate conversion coefficients. For the injection well, only the distr of the injection
    // phase is not zero.
    const double scale = well_.scalingFactor(well_.ebosCompIdxToFlowCompIdx(comp_idx));
    if (scale > 0.) {
        return this->volumeFraction(seg, comp_idx) / scale;
    }

    return this->volumeFraction(seg, comp_idx);
}

template<class FluidSystem, class Indices, class Scalar>
typename MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::
surfaceVolumeFraction(const int seg,
                      const int comp_idx) const
{
    EvalWell sum_volume_fraction_scaled = 0.;
    for (int idx = 0; idx < well_.numComponents(); ++idx) {
        sum_volume_fraction_scaled += this->volumeFractionScaled(seg, idx);
    }

    assert(sum_volume_fraction_scaled.value() != 0.);

    return this->volumeFractionScaled(seg, comp_idx) / sum_volume_fraction_scaled;
}

template<class FluidSystem, class Indices, class Scalar>
typename MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::
getSegmentRateUpwinding(const int seg,
                        const int seg_upwind,
                        const std::size_t comp_idx) const
{
    // the result will contain the derivative with respect to WQTotal in segment seg,
    // and the derivatives with respect to WFrac GFrac in segment seg_upwind.
    // the derivative with respect to SPres should be zero.
    if (seg == 0 && well_.isInjector()) {
        const Well& well = well_.wellEcl();
        auto phase = well.getInjectionProperties().injectorType;

        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)
                && Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx) == comp_idx
                && phase == InjectorType::WATER)
            return evaluation_[seg][WQTotal] / well_.scalingFactor(well_.ebosCompIdxToFlowCompIdx(comp_idx));

        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)
                && Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx) == comp_idx
                && phase == InjectorType::OIL)
            return evaluation_[seg][WQTotal] / well_.scalingFactor(well_.ebosCompIdxToFlowCompIdx(comp_idx));

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)
                && Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx) == comp_idx
                && phase == InjectorType::GAS)
            return evaluation_[seg][WQTotal] / well_.scalingFactor(well_.ebosCompIdxToFlowCompIdx(comp_idx));

        return 0.0;
    }

    const EvalWell segment_rate = evaluation_[seg][WQTotal] *
                                  this->volumeFractionScaled(seg_upwind, comp_idx);

    assert(segment_rate.derivative(SPres + Indices::numEq) == 0.);

    return segment_rate;
}

template<class FluidSystem, class Indices, class Scalar>
typename MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::
getSegmentPressure(const int seg) const
{
    return evaluation_[seg][SPres];
}

template<class FluidSystem, class Indices, class Scalar>
typename MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::
getBhp() const
{
    return this->getSegmentPressure(0);
}

template<class FluidSystem, class Indices, class Scalar>
typename MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::
getSegmentRate(const int seg,
               const int comp_idx) const
{
    return evaluation_[seg][WQTotal] * this->volumeFractionScaled(seg, comp_idx);
}

template<class FluidSystem, class Indices, class Scalar>
typename MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::
getQs(const int comp_idx) const
{
    return this->getSegmentRate(0, comp_idx);
}

template<class FluidSystem, class Indices, class Scalar>
typename MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>::
getWQTotal() const
{
    return evaluation_[0][WQTotal];
}

#define INSTANCE(...) \
template class MultisegmentWellPrimaryVariables<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,__VA_ARGS__,double>;

// One phase
INSTANCE(BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(BlackOilOnePhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,5u>)

// Two phase
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,0u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,false,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,0u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,2u,0u,false,false,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,0u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,true,0u,0u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<1u,0u,0u,0u,false,false,0u,0u,0u>)
// Blackoil
INSTANCE(BlackOilIndices<0u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<1u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,1u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,1u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,1u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,0u,true,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,0u,false,true,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,1u,false,true,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,0u,false,false,1u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,0u,false,true,2u,0u>)

INSTANCE(BlackOilIndices<1u,0u,0u,0u,true,false,0u,0u>)

}
