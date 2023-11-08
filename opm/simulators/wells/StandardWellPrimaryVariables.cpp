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
#include <opm/simulators/wells/StandardWellPrimaryVariables.hpp>

#include <opm/common/Exceptions.hpp>
#include <opm/input/eclipse/Units/Units.hpp>

#include <dune/common/dynvector.hh>
#include <dune/istl/bvector.hh>

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/models/blackoil/blackoilindices.hh>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>

#include <opm/simulators/wells/WellInterfaceIndices.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <limits>
#include <string>

namespace {

/**
    * @brief Relaxation factor considering only one fraction value.
    * @param old_value The previous value before the change.
    * @param dx The solved newton update for the value.
    * @param well_name The name of the well that calls this function
    * @param value_name The name/meaning of the old_value.
    * @param deferred_logger The deferred logger used for logging results and messages.
    * @return The calculated relaxation factor.
 */
template<class Scalar>
Scalar relaxationFactorFraction(const Scalar old_value,
                                const Scalar dx,
                                const std::string& well_name,
                                const std::string& value_name,
                                Opm::DeferredLogger& deferred_logger)
{
    constexpr double epislon = std::numeric_limits<Scalar>::epsilon();
    if (old_value < -epislon ||  old_value > 1.0 + epislon) {
        const std::string msg = fmt::format(" illegal fraction value {} {} is found for well {}", value_name, old_value, well_name);
        OPM_DEFLOG_PROBLEM(Opm::NumericalProblem, msg, deferred_logger);
    }
    const Scalar& safe_old_value = std::clamp(old_value, 0.0, 1.0);

    Scalar relaxation_factor = 1.;

    // updated values without relaxation factor
    const Scalar possible_updated_value = safe_old_value - dx;

    // 0.95 is an experimental value remains to be optimized
    if (possible_updated_value < 0.0) {
        relaxation_factor = std::abs(safe_old_value / dx) * 0.95;
    } else if (possible_updated_value > 1.0) {
        relaxation_factor = std::abs((1. - safe_old_value) / dx) * 0.95;
    }
    // if possible_updated_value is between 0. and 1.0, then relaxation_factor
    // remains to be one

    assert(relaxation_factor >= 0. && relaxation_factor <= 1.);

    return relaxation_factor;
}

}

namespace Opm {

template<class FluidSystem, class Indices, class Scalar>
void StandardWellPrimaryVariables<FluidSystem,Indices,Scalar>::
init()
{
    for (int eqIdx = 0; eqIdx < numWellEq_; ++eqIdx) {
        evaluation_[eqIdx] =
            EvalWell::createVariable(numWellEq_ + Indices::numEq,
                                     value_[eqIdx],
                                     Indices::numEq + eqIdx);

    }
}

template<class FluidSystem, class Indices, class Scalar>
void StandardWellPrimaryVariables<FluidSystem,Indices,Scalar>::
resize(const int numWellEq)
{
    value_.resize(numWellEq, 0.0);
    evaluation_.resize(numWellEq, EvalWell{numWellEq + Indices::numEq, 0.0});
    numWellEq_ = numWellEq;
}

template<class FluidSystem, class Indices, class Scalar>
void StandardWellPrimaryVariables<FluidSystem,Indices,Scalar>::
update(const WellState& well_state,
       const bool stop_or_zero_rate_target,
       DeferredLogger& deferred_logger)
{
    static constexpr int Water = BlackoilPhases::Aqua;
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Gas = BlackoilPhases::Vapour;

    const int well_index = well_.indexOfWell();
    const int np = well_.numPhases();
    const auto& pu = well_.phaseUsage();
    const auto& ws = well_state.well(well_index);
    // the weighted total well rate
    double total_well_rate = 0.0;
    for (int p = 0; p < np; ++p) {
        total_well_rate += well_.scalingFactor(p) * ws.surface_rates[p];
    }

    // Not: for the moment, the first primary variable for the injectors is not G_total. The injection rate
    // under surface condition is used here
    if (well_.isInjector()) {
        switch (well_.wellEcl().injectorType()) {
        case InjectorType::WATER:
            value_[WQTotal] = ws.surface_rates[pu.phase_pos[Water]];
            break;
        case InjectorType::GAS:
            value_[WQTotal] = ws.surface_rates[pu.phase_pos[Gas]];
            break;
        case InjectorType::OIL:
            value_[WQTotal] = ws.surface_rates[pu.phase_pos[Oil]];
            break;
        case InjectorType::MULTI:
            // Not supported.
            deferred_logger.warning("MULTI_PHASE_INJECTOR_NOT_SUPPORTED",
                                    "Multi phase injectors are not supported, requested for well " + well_.name());
            break;
        }
    } else {
            value_[WQTotal] = total_well_rate;
            if (stop_or_zero_rate_target) {
                value_[WQTotal] = 0.;
            }
    }

    if (std::abs(total_well_rate) > 0.) {
        if constexpr (has_wfrac_variable) {
            value_[WFrac] = well_.scalingFactor(pu.phase_pos[Water]) * ws.surface_rates[pu.phase_pos[Water]] / total_well_rate;
        }
        if constexpr (has_gfrac_variable) {
            value_[GFrac] = well_.scalingFactor(pu.phase_pos[Gas]) *
                            (ws.surface_rates[pu.phase_pos[Gas]] -
                             (Indices::enableSolvent ? ws.sum_solvent_rates() : 0.0) ) / total_well_rate ;
        }
        if constexpr (Indices::enableSolvent) {
            value_[SFrac] = well_.scalingFactor(Indices::contiSolventEqIdx) * ws.sum_solvent_rates() / total_well_rate ;
        }

    } else { // total_well_rate == 0
        if (well_.isInjector()) {
            // only single phase injection handled
            if constexpr (has_wfrac_variable) {
                if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                    auto phase = well_.wellEcl().getInjectionProperties().injectorType;
                    if (phase == InjectorType::WATER) {
                        value_[WFrac] = 1.0;
                    } else {
                        value_[WFrac] = 0.0;
                    }
                }
            }
            if constexpr (has_gfrac_variable) {
                if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                    auto phase = well_.wellEcl().getInjectionProperties().injectorType;
                    if (phase == InjectorType::GAS) {
                        value_[GFrac] = (1.0 - well_.rsRvInj());
                        if constexpr (Indices::enableSolvent) {
                            value_[GFrac] = 1.0 - well_.rsRvInj() - well_.wsolvent();
                            value_[SFrac] = well_.wsolvent();
                        }
                    } else {
                        value_[GFrac] = 0.0;
                    }
                }
            }

            // TODO: it is possible to leave injector as a oil well,
            // when F_w and F_g both equals to zero, not sure under what kind of circumstance
            // this will happen.
        } else if (well_.isProducer()) { // producers
            // TODO: the following are not addressed for the solvent case yet
            if constexpr (has_wfrac_variable) {
                value_[WFrac] = 1.0 / np;
            }

            if constexpr (has_gfrac_variable) {
                value_[GFrac] = 1.0 / np;
            }
        } else {
            OPM_DEFLOG_THROW(std::logic_error, "Expected PRODUCER or INJECTOR type of well", deferred_logger);
        }
    }

    // BHP
    value_[Bhp] = ws.bhp;
}

template<class FluidSystem, class Indices, class Scalar>
void StandardWellPrimaryVariables<FluidSystem,Indices,Scalar>::
updatePolyMW(const WellState& well_state)
{
    if (well_.isInjector()) {
        const auto& ws = well_state.well(well_.indexOfWell());
        const auto& perf_data = ws.perf_data;
        const auto& water_velocity = perf_data.water_velocity;
        const auto& skin_pressure = perf_data.skin_pressure;
        for (int perf = 0; perf < well_.numPerfs(); ++perf) {
            value_[Bhp + 1 + perf] = water_velocity[perf];
            value_[Bhp + 1 + well_.numPerfs() + perf] = skin_pressure[perf];
        }
    }
}

template<class FluidSystem, class Indices, class Scalar>
void StandardWellPrimaryVariables<FluidSystem,Indices,Scalar>::
updateNewton(const BVectorWell& dwells,
             const bool stop_or_zero_rate_target,
             const double dFLimit,
             const double dBHPLimit,
             DeferredLogger& deferred_logger)
{
    // for injectors, very typical one of the fractions will be one, and it is easy to get zero value
    // fractions. not sure what is the best way to handle it yet, so we just use 1.0 here
    // The relaxationFactorFractionProducer code does not take into account solvent
    // so we use 1.0 for cases with solvent.
    [[maybe_unused]] const double relaxation_factor_fractions =
            (well_.isProducer() && !Indices::enableSolvent) ? this->relaxationFactorFractionsProducer(dwells, deferred_logger)
                                                            : 1.0;

    // update the second and third well variable (The flux fractions)
    if constexpr (has_wfrac_variable) {
        const int sign2 = dwells[0][WFrac] > 0 ? 1: -1;
        const double dx2_limited = sign2 * std::min(std::abs(dwells[0][WFrac] * relaxation_factor_fractions), dFLimit);
        value_[WFrac] = value_[WFrac] - dx2_limited;
    }

    if constexpr (has_gfrac_variable) {
        const int sign3 = dwells[0][GFrac] > 0 ? 1: -1;
        const double dx3_limited = sign3 * std::min(std::abs(dwells[0][GFrac] * relaxation_factor_fractions), dFLimit);
        value_[GFrac] = value_[GFrac] - dx3_limited;
    }

    if constexpr (Indices::enableSolvent) {
        const int sign4 = dwells[0][SFrac] > 0 ? 1: -1;
        const double dx4_limited = sign4 * std::min(std::abs(dwells[0][SFrac]) * relaxation_factor_fractions, dFLimit);
        value_[SFrac] = value_[SFrac] - dx4_limited;
    }

    this->processFractions();

    // updating the total rates Q_t
    value_[WQTotal] -= dwells[0][WQTotal];

    // here, we make sure it is zero for wells with zero rate target(including stopped wells)
    if (stop_or_zero_rate_target) {
        value_[WQTotal] = 0.;
    } else {
        // make sure that no injector produce and no producer inject
        if (well_.isInjector()) {
            value_[WQTotal] = std::max(value_[WQTotal], 0.0);
        } else {
            value_[WQTotal] = std::min(value_[WQTotal], 0.0);
        }
    }

    // updating the bottom hole pressure
    const int sign1 = dwells[0][Bhp] > 0 ? 1: -1;
    const double dx1_limited = sign1 * std::min(std::abs(dwells[0][Bhp]),
                                                std::abs(value_[Bhp]) * dBHPLimit);
    // some cases might have defaulted bhp constraint of 1 bar, we use a slightly smaller value as the bhp lower limit for Newton update
    // so that bhp constaint can be an active control when needed.
    constexpr double bhp_lower_limit = 1. * unit::barsa - 1. * unit::Pascal;
    value_[Bhp] = std::max(value_[Bhp] - dx1_limited, bhp_lower_limit);
}

template<class FluidSystem, class Indices, class Scalar>
void StandardWellPrimaryVariables<FluidSystem,Indices,Scalar>::
updateNewtonPolyMW(const BVectorWell& dwells)
{
    if (well_.isInjector()) {
        for (int perf = 0; perf < well_.numPerfs(); ++perf) {
            const int wat_vel_index = Bhp + 1 + perf;
            const int pskin_index = Bhp + 1 + well_.numPerfs() + perf;

            const double relaxation_factor = 0.9;
            const double dx_wat_vel = dwells[0][wat_vel_index];
            value_[wat_vel_index] -= relaxation_factor * dx_wat_vel;

            const double dx_pskin = dwells[0][pskin_index];
            value_[pskin_index] -= relaxation_factor * dx_pskin;
        }
    }
}

template<class FluidSystem, class Indices, class Scalar>
void StandardWellPrimaryVariables<FluidSystem,Indices,Scalar>::
copyToWellState(WellState& well_state,
                DeferredLogger& deferred_logger) const
{
    static constexpr int Water = BlackoilPhases::Aqua;
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Gas = BlackoilPhases::Vapour;

    const PhaseUsage& pu = well_.phaseUsage();
    std::vector<double> F(well_.numPhases(), 0.0);
    [[maybe_unused]] double F_solvent = 0.0;
    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
        const int oil_pos = pu.phase_pos[Oil];
        F[oil_pos] = 1.0;

        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            const int water_pos = pu.phase_pos[Water];
            F[water_pos] = value_[WFrac];
            F[oil_pos] -= F[water_pos];
        }

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            const int gas_pos = pu.phase_pos[Gas];
            F[gas_pos] = value_[GFrac];
            F[oil_pos] -= F[gas_pos];
        }

        if constexpr (Indices::enableSolvent) {
            F_solvent = value_[SFrac];
            F[oil_pos] -= F_solvent;
        }
    }
    else if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
        const int water_pos = pu.phase_pos[Water];
        F[water_pos] = 1.0;

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            const int gas_pos = pu.phase_pos[Gas];
            F[gas_pos] = value_[GFrac];
            F[water_pos] -= F[gas_pos];
        }

        if constexpr (Indices::enableSolvent) {
            F_solvent = value_[SFrac];
            F[water_pos] -= F_solvent;
        }
    }
    else if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
        const int gas_pos = pu.phase_pos[Gas];
        F[gas_pos] = 1.0;
    }

    // convert the fractions to be Q_p / G_total to calculate the phase rates
    for (int p = 0; p < well_.numPhases(); ++p) {
        const double scal = well_.scalingFactor(p);
        // for injection wells, there should only one non-zero scaling factor
        if (scal > 0) {
            F[p] /= scal ;
        } else {
            // this should only happens to injection wells
            F[p] = 0.;
        }
    }

    // F_solvent is added to F_gas. This means that well_rate[Gas] also contains solvent.
    // More testing is needed to make sure this is correct for well groups and THP.
    if constexpr (Indices::enableSolvent) {
        F_solvent /= well_.scalingFactor(Indices::contiSolventEqIdx);
        F[pu.phase_pos[Gas]] += F_solvent;
    }

    auto& ws = well_state.well(well_.indexOfWell());
    ws.bhp = value_[Bhp];

    // calculate the phase rates based on the primary variables
    // for producers, this is not a problem, while not sure for injectors here
    if (well_.isProducer()) {
        const double g_total = value_[WQTotal];
        for (int p = 0; p < well_.numPhases(); ++p) {
            ws.surface_rates[p] = g_total * F[p];
        }
    } else { // injectors
        for (int p = 0; p < well_.numPhases(); ++p) {
            ws.surface_rates[p] = 0.0;
        }
        switch (well_.wellEcl().injectorType()) {
        case InjectorType::WATER:
            ws.surface_rates[pu.phase_pos[Water]] = value_[WQTotal];
            break;
        case InjectorType::GAS:
            ws.surface_rates[pu.phase_pos[Gas]] = value_[WQTotal];
            break;
        case InjectorType::OIL:
            ws.surface_rates[pu.phase_pos[Oil]] = value_[WQTotal];
            break;
        case InjectorType::MULTI:
            // Not supported.
            deferred_logger.warning("MULTI_PHASE_INJECTOR_NOT_SUPPORTED",
                                    "Multi phase injectors are not supported, requested for well " + well_.name());
            break;
        }
    }
}

template<class FluidSystem, class Indices, class Scalar>
void StandardWellPrimaryVariables<FluidSystem,Indices,Scalar>::
copyToWellStatePolyMW(WellState& well_state) const
{
    if (well_.isInjector()) {
        auto& ws = well_state.well(well_.indexOfWell());
        auto& perf_data = ws.perf_data;
        auto& perf_water_velocity = perf_data.water_velocity;
        auto& perf_skin_pressure = perf_data.skin_pressure;
        for (int perf = 0; perf < well_.numPerfs(); ++perf) {
            perf_water_velocity[perf] = value_[Bhp + 1 + perf];
            perf_skin_pressure[perf] = value_[Bhp + 1 + well_.numPerfs() + perf];
        }
    }
}

template<class FluidSystem, class Indices, class Scalar>
typename StandardWellPrimaryVariables<FluidSystem,Indices,Scalar>::EvalWell
StandardWellPrimaryVariables<FluidSystem,Indices,Scalar>::
volumeFraction(const unsigned compIdx) const
{
    if (FluidSystem::numActivePhases() == 1) {
        return EvalWell(numWellEq_ + Indices::numEq, 1.0);
    }

    if (has_gfrac_variable && compIdx == Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx)) {
        return evaluation_[GFrac];
    }

    if (Indices::enableSolvent && compIdx == (unsigned)Indices::contiSolventEqIdx) {
        return evaluation_[SFrac];
    }

    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
        if (has_wfrac_variable && compIdx == Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx)) {
            return evaluation_[WFrac];
        }
    }
    // Compute the Oil fraction if oil is present
    // or the WATER fraction for a gas-water case
    EvalWell well_fraction(numWellEq_ + Indices::numEq, 1.0);
    if (has_wfrac_variable && FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
        well_fraction -= evaluation_[WFrac];
    }
    if (has_gfrac_variable) {
        well_fraction -= evaluation_[GFrac];
    }
    if constexpr (Indices::enableSolvent) {
        well_fraction -= evaluation_[SFrac];
    }
    return well_fraction;
}

template<class FluidSystem, class Indices, class Scalar>
typename StandardWellPrimaryVariables<FluidSystem,Indices,Scalar>::EvalWell
StandardWellPrimaryVariables<FluidSystem,Indices,Scalar>::
volumeFractionScaled(const int compIdx) const
{
    const int legacyCompIdx = well_.ebosCompIdxToFlowCompIdx(compIdx);
    const double scal = well_.scalingFactor(legacyCompIdx);
    if (scal > 0)
        return this->volumeFraction(compIdx) / scal;

    // the scaling factor may be zero for RESV controlled wells.
    return this->volumeFraction(compIdx);
}

template<class FluidSystem, class Indices, class Scalar>
typename StandardWellPrimaryVariables<FluidSystem,Indices,Scalar>::EvalWell
StandardWellPrimaryVariables<FluidSystem,Indices,Scalar>::
surfaceVolumeFraction(const int compIdx) const
{
    EvalWell sum_volume_fraction_scaled(numWellEq_ + Indices::numEq, 0.);
    for (int idx = 0; idx < well_.numComponents(); ++idx) {
        sum_volume_fraction_scaled += this->volumeFractionScaled(idx);
    }

    assert(sum_volume_fraction_scaled.value() != 0.);

    return this->volumeFractionScaled(compIdx) / sum_volume_fraction_scaled;
 }

template<class FluidSystem, class Indices, class Scalar>
typename StandardWellPrimaryVariables<FluidSystem,Indices,Scalar>::EvalWell
StandardWellPrimaryVariables<FluidSystem,Indices,Scalar>::
getQs(const int comp_idx) const
{
    // Note: currently, the WQTotal definition is still depends on Injector/Producer.
    assert(comp_idx < well_.numComponents());

    if (well_.isInjector()) { // only single phase injection
        double inj_frac = 0.0;
        switch (well_.wellEcl().injectorType()) {
        case InjectorType::WATER:
            if (comp_idx == int(Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx))) {
                inj_frac = 1.0;
            }
            break;
        case InjectorType::GAS:
            if (Indices::enableSolvent && comp_idx == Indices::contiSolventEqIdx) { // solvent
                inj_frac = well_.wsolvent();
            } else if (comp_idx == int(Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx))) {
                inj_frac = 1.0 - well_.rsRvInj();
                if constexpr (Indices::enableSolvent) {
                    inj_frac -= well_.wsolvent();
                }
            } else if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && comp_idx == int(Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx))) {
                inj_frac = well_.rsRvInj();
            }
            break;
        case InjectorType::OIL:
            if (comp_idx == int(Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx))) {
                inj_frac = 1.0 - well_.rsRvInj();
            } else if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && comp_idx == int(Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx))) {
                inj_frac = well_.rsRvInj();
            }
            break;
        case InjectorType::MULTI:
            // Not supported.
            // deferred_logger.warning("MULTI_PHASE_INJECTOR_NOT_SUPPORTED",
            //                         "Multi phase injectors are not supported, requested for well " + name());
            break;
        }
        return inj_frac * evaluation_[WQTotal];
    } else { // producers
        return evaluation_[WQTotal] * this->volumeFractionScaled(comp_idx);
    }
}

template<class FluidSystem, class Indices, class Scalar>
void StandardWellPrimaryVariables<FluidSystem,Indices,Scalar>::
processFractions()
{
    static constexpr int Water = BlackoilPhases::Aqua;
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Gas = BlackoilPhases::Vapour;

    const auto pu = well_.phaseUsage();
    std::vector<double> F(well_.numPhases(), 0.0);

    [[maybe_unused]] double F_solvent = 0.0;

    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
        F[pu.phase_pos[Oil]] = 1.0;

        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            F[pu.phase_pos[Water]] = value_[WFrac];
            F[pu.phase_pos[Oil]] -= F[pu.phase_pos[Water]];
        }

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            F[pu.phase_pos[Gas]] = value_[GFrac];
            F[pu.phase_pos[Oil]] -= F[pu.phase_pos[Gas]];
        }
        if constexpr (Indices::enableSolvent) {
            F_solvent = value_[SFrac];
            F[pu.phase_pos[Oil]] -= F_solvent;
        }
    }
    else if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
        F[pu.phase_pos[Water]] = 1.0;

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            F[pu.phase_pos[Gas]] = value_[GFrac];
            F[pu.phase_pos[Water]] -= F[pu.phase_pos[Gas]];
        }
        if constexpr (Indices::enableSolvent) {
            F_solvent = value_[SFrac];
            F[pu.phase_pos[Water]] -= F_solvent;
        }
    }
    else if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
        F[pu.phase_pos[Gas]] = 1.0;
    }

    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
        if (F[pu.phase_pos[Water]] < 0.0) {
            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                    F[pu.phase_pos[Gas]] /= (1.0 - F[pu.phase_pos[Water]]);
            }
            if constexpr (Indices::enableSolvent) {
                F_solvent /= (1.0 - F[pu.phase_pos[Water]]);
            }
            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                F[pu.phase_pos[Oil]] /= (1.0 - F[pu.phase_pos[Water]]);
            }
            F[pu.phase_pos[Water]] = 0.0;
        }
    }

    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
        if (F[pu.phase_pos[Gas]] < 0.0) {
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                F[pu.phase_pos[Water]] /= (1.0 - F[pu.phase_pos[Gas]]);
            }
            if constexpr (Indices::enableSolvent) {
                F_solvent /= (1.0 - F[pu.phase_pos[Gas]]);
            }
            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                F[pu.phase_pos[Oil]] /= (1.0 - F[pu.phase_pos[Gas]]);
            }
            F[pu.phase_pos[Gas]] = 0.0;
        }
    }

    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
        if (F[pu.phase_pos[Oil]] < 0.0) {
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                F[pu.phase_pos[Water]] /= (1.0 - F[pu.phase_pos[Oil]]);
            }
            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                F[pu.phase_pos[Gas]] /= (1.0 - F[pu.phase_pos[Oil]]);
            }
            if constexpr (Indices::enableSolvent) {
                F_solvent /= (1.0 - F[pu.phase_pos[Oil]]);
            }
            F[pu.phase_pos[Oil]] = 0.0;
        }
    }

    if constexpr (has_wfrac_variable) {
        value_[WFrac] = F[pu.phase_pos[Water]];
    }

    if constexpr (has_gfrac_variable) {
        value_[GFrac] = F[pu.phase_pos[Gas]];
    }
    if constexpr (Indices::enableSolvent) {
        value_[SFrac] = F_solvent;
    }
}

template<class FluidSystem, class Indices, class Scalar>
double StandardWellPrimaryVariables<FluidSystem,Indices,Scalar>::
relaxationFactorFractionsProducer(const BVectorWell& dwells, DeferredLogger& deferred_logger) const
{
    // TODO: not considering solvent yet
    // 0.95 is a experimental value, which remains to be optimized
    double relaxation_factor = 1.0;

    if (FluidSystem::numActivePhases() > 1) {
        if constexpr (has_wfrac_variable) {
            const double relaxation_factor_w = relaxationFactorFraction(value_[WFrac],
                                                                        dwells[0][WFrac],
                                                                        this->well_.name(),
                                                                        "WFrac",
                                                                        deferred_logger);
            relaxation_factor = std::min(relaxation_factor, relaxation_factor_w);
        }

        if constexpr (has_gfrac_variable) {
            const double relaxation_factor_g = relaxationFactorFraction(value_[GFrac],
                                                                        dwells[0][GFrac],
                                                                        this->well_.name(),
                                                                        "GFrac",
                                                                        deferred_logger);
            relaxation_factor = std::min(relaxation_factor, relaxation_factor_g);
        }


        if constexpr (has_wfrac_variable && has_gfrac_variable) {
            // We need to make sure the even with the relaxation_factor, the sum of F_w and F_g is below one, so there will
            // not be negative oil fraction later
            const double original_sum = value_[WFrac] + value_[GFrac];
            const double relaxed_update = (dwells[0][WFrac] + dwells[0][GFrac]) * relaxation_factor;
            const double possible_updated_sum = original_sum - relaxed_update;
            // We only relax if fraction is above 1.
            // The newton solver should handle the rest
            constexpr double epsilon = 0.001;
            if (possible_updated_sum > 1.0 + epsilon) {
                // since the orignal sum <= 1.0 the epsilon asserts that
                // the relaxed_update is non trivial.
                assert(relaxed_update != 0.);

                const double further_relaxation_factor = std::abs((1. - original_sum) / relaxed_update) * 0.95;
                relaxation_factor *= further_relaxation_factor;
            }
        }
        if (relaxation_factor < 0.0 || relaxation_factor > 1.0) {
            const std::string msg = fmt::format(" illegal relaxation factor {} is obtained for well {}", relaxation_factor, this->well_.name());
            OPM_DEFLOG_PROBLEM(NumericalProblem, msg, deferred_logger);
        }
    }
    return relaxation_factor;
}

template<class FluidSystem, class Indices, class Scalar>
void StandardWellPrimaryVariables<FluidSystem,Indices,Scalar>::
checkFinite(DeferredLogger& deferred_logger) const
{
    for (const Scalar v : value_) {
        if (!isfinite(v))
            OPM_DEFLOG_PROBLEM(NumericalProblem,
                               "Infinite primary variable after update from wellState, well: " + well_.name(),
                               deferred_logger);
    }
}

#define INSTANCE(...) \
template class StandardWellPrimaryVariables<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,__VA_ARGS__,double>;

// One phase
INSTANCE(BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(BlackOilOnePhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,5u>)

// Two phase
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,0u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,false,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,true,0u,2u,0u>)
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
INSTANCE(BlackOilIndices<0u,0u,0u,1u,false,false,1u,0u>)

INSTANCE(BlackOilIndices<1u,0u,0u,0u,true,false,0u,0u>)

}
