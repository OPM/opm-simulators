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
#include <opm/simulators/wells/StandardWellEval.hpp>

#include <opm/material/densead/DynamicEvaluation.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/models/blackoil/blackoilindices.hh>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/simulators/timestepping/ConvergenceReport.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>
#include <opm/simulators/wells/WellInterfaceIndices.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <cassert>
#include <cmath>

#if HAVE_CUDA || HAVE_OPENCL
#include <opm/simulators/linalg/bda/WellContributions.hpp>
#endif


namespace Opm
{

template<class FluidSystem, class Indices, class Scalar>
StandardWellEval<FluidSystem,Indices,Scalar>::
StandardWellEval(const WellInterfaceIndices<FluidSystem,Indices,Scalar>& baseif)
    : StandardWellGeneric<Scalar>(Bhp, baseif)
    , baseif_(baseif)
    , F0_(numWellConservationEq)
{
}

template<class FluidSystem, class Indices, class Scalar>
void StandardWellEval<FluidSystem,Indices,Scalar>::
initPrimaryVariablesEvaluation() const
{
    for (int eqIdx = 0; eqIdx < numWellEq_; ++eqIdx) {
        primary_variables_evaluation_[eqIdx] =
            EvalWell::createVariable(numWellEq_ + Indices::numEq, primary_variables_[eqIdx], Indices::numEq + eqIdx);
    }
}

template<class FluidSystem, class Indices, class Scalar>
typename StandardWellEval<FluidSystem,Indices,Scalar>::EvalWell
StandardWellEval<FluidSystem,Indices,Scalar>::
extendEval(const Eval& in) const
{
    EvalWell out(numWellEq_ + Indices::numEq, in.value());
    for(int eqIdx = 0; eqIdx < Indices::numEq;++eqIdx) {
        out.setDerivative(eqIdx, in.derivative(eqIdx));
    }
    return out;
}

template<class FluidSystem, class Indices, class Scalar>
double
StandardWellEval<FluidSystem,Indices,Scalar>::
relaxationFactorFractionsProducer(const std::vector<double>& primary_variables,
                                  const BVectorWell& dwells)
{
    // TODO: not considering solvent yet
    // 0.95 is a experimental value, which remains to be optimized
    double relaxation_factor = 1.0;

    if (FluidSystem::numActivePhases() > 1) {
        if constexpr (has_wfrac_variable) {
            const double relaxation_factor_w = StandardWellGeneric<Scalar>::
                                               relaxationFactorFraction(primary_variables[WFrac], dwells[0][WFrac]);
            relaxation_factor = std::min(relaxation_factor, relaxation_factor_w);
        }

        if constexpr (has_gfrac_variable) {
            const double relaxation_factor_g = StandardWellGeneric<Scalar>::
                                               relaxationFactorFraction(primary_variables[GFrac], dwells[0][GFrac]);
            relaxation_factor = std::min(relaxation_factor, relaxation_factor_g);
        }


        if constexpr (has_wfrac_variable && has_gfrac_variable) {
            // We need to make sure the even with the relaxation_factor, the sum of F_w and F_g is below one, so there will
            // not be negative oil fraction later
            const double original_sum = primary_variables[WFrac] + primary_variables[GFrac];
            const double relaxed_update = (dwells[0][WFrac] + dwells[0][GFrac]) * relaxation_factor;
            const double possible_updated_sum = original_sum - relaxed_update;
            // We only relax if fraction is above 1.
            // The newton solver should handle the rest
            const double epsilon = 0.001;
            if (possible_updated_sum > 1.0 + epsilon) {
                // since the orignal sum <= 1.0 the epsilon asserts that
                // the relaxed_update is non trivial.
                assert(relaxed_update != 0.);

                const double further_relaxation_factor = std::abs((1. - original_sum) / relaxed_update) * 0.95;
                relaxation_factor *= further_relaxation_factor;
            }
        }
        assert(relaxation_factor >= 0.0 && relaxation_factor <= 1.0);
    }
    return relaxation_factor;
}

template<class FluidSystem, class Indices, class Scalar>
typename StandardWellEval<FluidSystem,Indices,Scalar>::EvalWell
StandardWellEval<FluidSystem,Indices,Scalar>::
wellVolumeFraction(const unsigned compIdx) const
{
    if (FluidSystem::numActivePhases() == 1) {
        return EvalWell(numWellEq_ + Indices::numEq, 1.0);
    }

    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
        if (has_wfrac_variable && compIdx == Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx)) {
            return primary_variables_evaluation_[WFrac];
        }

        if (has_gfrac_variable && compIdx == Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx)) {
            return primary_variables_evaluation_[GFrac];
        }

        if (Indices::enableSolvent && compIdx == (unsigned)Indices::contiSolventEqIdx) {
            return primary_variables_evaluation_[SFrac];
        }
    }
    else if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
        if (has_gfrac_variable && compIdx == Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx)) {
            return primary_variables_evaluation_[GFrac];
        }
    }

    // Oil or WATER fraction
    EvalWell well_fraction(numWellEq_ + Indices::numEq, 1.0);
    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            well_fraction -= primary_variables_evaluation_[WFrac];
        }

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            well_fraction -= primary_variables_evaluation_[GFrac];
        }

        if (Indices::enableSolvent) {
            well_fraction -= primary_variables_evaluation_[SFrac];
        }
    }
    else if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) && (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx))) {

            well_fraction -= primary_variables_evaluation_[GFrac];
    }

    return well_fraction;
}

template<class FluidSystem, class Indices, class Scalar>
typename StandardWellEval<FluidSystem,Indices,Scalar>::EvalWell
StandardWellEval<FluidSystem,Indices,Scalar>::
getQs(const int comp_idx) const
{
    // Note: currently, the WQTotal definition is still depends on Injector/Producer.
    assert(comp_idx < baseif_.numComponents());

    if (baseif_.isInjector()) { // only single phase injection
        double inj_frac = 0.0;
        switch (baseif_.wellEcl().injectorType()) {
        case InjectorType::WATER:
            if (comp_idx == int(Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx))) {
                inj_frac = 1.0;
            }
            break;
        case InjectorType::GAS:
            if (Indices::enableSolvent && comp_idx == Indices::contiSolventEqIdx) { // solvent
                inj_frac = baseif_.wsolvent();
            } else if (comp_idx == int(Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx))) {
                inj_frac = 1.0 - baseif_.rsRvInj();
                if (Indices::enableSolvent) {
                    inj_frac -= baseif_.wsolvent();
                }
            } else if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && comp_idx == int(Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx))) {
                inj_frac = baseif_.rsRvInj();
            }
            break;
        case InjectorType::OIL:
            if (comp_idx == int(Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx))) {
                inj_frac = 1.0 - baseif_.rsRvInj();
            } else if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && comp_idx == int(Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx))) {
                inj_frac = baseif_.rsRvInj();
            }
            break;
        case InjectorType::MULTI:
            // Not supported.
            // deferred_logger.warning("MULTI_PHASE_INJECTOR_NOT_SUPPORTED",
            //                         "Multi phase injectors are not supported, requested for well " + name());
            break;
        }
        return inj_frac * primary_variables_evaluation_[WQTotal]*rate_scaling_;
    } else { // producers
        return primary_variables_evaluation_[WQTotal] * wellVolumeFractionScaled(comp_idx)*rate_scaling_;
    }
}

template<class FluidSystem, class Indices, class Scalar>
typename StandardWellEval<FluidSystem,Indices,Scalar>::EvalWell
StandardWellEval<FluidSystem,Indices,Scalar>::
wellVolumeFractionScaled(const int compIdx) const
{
    const int legacyCompIdx = baseif_.ebosCompIdxToFlowCompIdx(compIdx);
    const double scal = baseif_.scalingFactor(legacyCompIdx);
    if (scal > 0)
        return  wellVolumeFraction(compIdx) / scal;

    // the scaling factor may be zero for RESV controlled wells.
    return wellVolumeFraction(compIdx);
}

template<class FluidSystem, class Indices, class Scalar>
typename StandardWellEval<FluidSystem,Indices,Scalar>::EvalWell
StandardWellEval<FluidSystem,Indices,Scalar>::
wellSurfaceVolumeFraction(const int compIdx) const
{
    EvalWell sum_volume_fraction_scaled(numWellEq_ + Indices::numEq, 0.);
    for (int idx = 0; idx < baseif_.numComponents(); ++idx) {
        sum_volume_fraction_scaled += wellVolumeFractionScaled(idx);
    }

    assert(sum_volume_fraction_scaled.value() != 0.);

    return wellVolumeFractionScaled(compIdx) / sum_volume_fraction_scaled;
 }

template<class FluidSystem, class Indices, class Scalar>
void
StandardWellEval<FluidSystem,Indices,Scalar>::
updatePrimaryVariables(const WellState& well_state, DeferredLogger& deferred_logger) const
{
    static constexpr int Gas = WellInterfaceIndices<FluidSystem,Indices,Scalar>::Gas;
    static constexpr int Oil = WellInterfaceIndices<FluidSystem,Indices,Scalar>::Oil;
    static constexpr int Water = WellInterfaceIndices<FluidSystem,Indices,Scalar>::Water;

    if (!baseif_.isOperableAndSolvable() && !baseif_.wellIsStopped()) return;

    const int well_index = baseif_.indexOfWell();
    const int np = baseif_.numPhases();
    const auto& pu = baseif_.phaseUsage();
    const auto& ws = well_state.well(well_index);
    // the weighted total well rate
    double total_well_rate = 0.0;
    for (int p = 0; p < np; ++p) {
        total_well_rate += baseif_.scalingFactor(p) * ws.surface_rates[p];
    }

    // Not: for the moment, the first primary variable for the injectors is not G_total. The injection rate
    // under surface condition is used here
    if (baseif_.isInjector()) {
        switch (baseif_.wellEcl().injectorType()) {
        case InjectorType::WATER:
            primary_variables_[WQTotal] = ws.surface_rates[pu.phase_pos[Water]]/rate_scaling_;
            break;
        case InjectorType::GAS:
            primary_variables_[WQTotal] = ws.surface_rates[pu.phase_pos[Gas]]/rate_scaling_;
            break;
        case InjectorType::OIL:
            primary_variables_[WQTotal] = ws.surface_rates[pu.phase_pos[Oil]]/rate_scaling_;
            break;
        case InjectorType::MULTI:
            // Not supported.
            deferred_logger.warning("MULTI_PHASE_INJECTOR_NOT_SUPPORTED",
                                    "Multi phase injectors are not supported, requested for well " + baseif_.name());
            break;
        }
    } else {
            primary_variables_[WQTotal] = total_well_rate/rate_scaling_;
    }

    if (std::abs(total_well_rate) > 0.) {
        if constexpr (has_wfrac_variable) {
            primary_variables_[WFrac] = baseif_.scalingFactor(pu.phase_pos[Water]) * ws.surface_rates[pu.phase_pos[Water]] / total_well_rate;
        }
        if constexpr (has_gfrac_variable) {
            primary_variables_[GFrac] = baseif_.scalingFactor(pu.phase_pos[Gas]) * (ws.surface_rates[pu.phase_pos[Gas]]
                                                                                    - (Indices::enableSolvent ? ws.sum_solvent_rates() : 0.0) ) / total_well_rate ;
        }
        if constexpr (Indices::enableSolvent) {
            primary_variables_[SFrac] = baseif_.scalingFactor(pu.phase_pos[Gas]) * ws.sum_solvent_rates() / total_well_rate ;
        }
    } else { // total_well_rate == 0
        if (baseif_.isInjector()) {
            // only single phase injection handled
            if constexpr (has_wfrac_variable) {
                if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                    auto phase = baseif_.wellEcl().getInjectionProperties().injectorType;
                    if (phase == InjectorType::WATER) {
                        primary_variables_[WFrac] = 1.0;
                    } else {
                        primary_variables_[WFrac] = 0.0;
                    }
                }
            }
            if constexpr (has_gfrac_variable) {
                if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                    auto phase = baseif_.wellEcl().getInjectionProperties().injectorType;
                    if (phase == InjectorType::GAS) {
                        primary_variables_[GFrac] = (1.0 - baseif_.rsRvInj());
                        if constexpr (Indices::enableSolvent) {
                            primary_variables_[GFrac] = 1.0 - baseif_.rsRvInj() - baseif_.wsolvent();
                            primary_variables_[SFrac] = baseif_.wsolvent();
                        }
                    } else {
                        primary_variables_[GFrac] = 0.0;
                    }
                }
            }

            // TODO: it is possible to leave injector as a oil well,
            // when F_w and F_g both equals to zero, not sure under what kind of circumstance
            // this will happen.
        } else if (baseif_.isProducer()) { // producers
            // TODO: the following are not addressed for the solvent case yet
            if constexpr (has_wfrac_variable) {
                primary_variables_[WFrac] = 1.0 / np;
            }

            if constexpr (has_gfrac_variable) {
                primary_variables_[GFrac] = 1.0 / np;
            }
        } else {
            OPM_DEFLOG_THROW(std::logic_error, "Expected PRODUCER or INJECTOR type of well", deferred_logger);
        }
    }


    // BHP
    primary_variables_[Bhp] = ws.bhp/bhp_scaling_;
}

template<class FluidSystem, class Indices, class Scalar>
void
StandardWellEval<FluidSystem,Indices,Scalar>::
assembleControlEq(const WellState& well_state,
                  const GroupState& group_state,
                  const Schedule& schedule,
                  const SummaryState& summaryState,
                  DeferredLogger& deferred_logger)
{
    static constexpr int Gas = WellInterfaceIndices<FluidSystem,Indices,Scalar>::Gas;
    static constexpr int Oil = WellInterfaceIndices<FluidSystem,Indices,Scalar>::Oil;
    static constexpr int Water = WellInterfaceIndices<FluidSystem,Indices,Scalar>::Water;
    EvalWell control_eq(numWellEq_ + Indices::numEq, 0.0);

    const auto& well = baseif_.wellEcl();

    auto getRates = [&]() {
        std::vector<EvalWell> rates(3, EvalWell(numWellEq_ + Indices::numEq, 0.0));
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
    } else if (baseif_.isInjector()) {
        // Find injection rate.
        const EvalWell injection_rate = getWQTotal();
        // Setup function for evaluation of BHP from THP (used only if needed).
        auto bhp_from_thp = [&]() {
            const auto rates = getRates();
            return baseif_.calculateBhpFromThp(well_state, rates, well, summaryState, this->getRho(), deferred_logger);
        };
        // Call generic implementation.
        const auto& inj_controls = well.injectionControls(summaryState);
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
             return baseif_.calculateBhpFromThp(well_state, rates, well, summaryState, this->getRho(), deferred_logger);
        };
        // Call generic implementation.
        const auto& prod_controls = well.productionControls(summaryState);
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
    // TODO: we should use a different index system for the well equations
    this->resWell_[0][Bhp] = control_eq.value();
    for (int pv_idx = 0; pv_idx < numWellEq_; ++pv_idx) {
        this->duneD_[0][0][Bhp][pv_idx] = control_eq.derivative(pv_idx + Indices::numEq);
    }
}


template<class FluidSystem, class Indices, class Scalar>
void
StandardWellEval<FluidSystem,Indices,Scalar>::
updatePrimaryVariablesPolyMW(const BVectorWell& dwells) const
{
    if (baseif_.isInjector()) {
        for (int perf = 0; perf < baseif_.numPerfs(); ++perf) {
            const int wat_vel_index = Bhp + 1 + perf;
            const int pskin_index = Bhp + 1 + baseif_.numPerfs() + perf;

            const double relaxation_factor = 0.9;
            const double dx_wat_vel = dwells[0][wat_vel_index];
            primary_variables_[wat_vel_index] -= relaxation_factor * dx_wat_vel;

            const double dx_pskin = dwells[0][pskin_index];
            primary_variables_[pskin_index] -= relaxation_factor * dx_pskin;
        }
    }
}

template<class FluidSystem, class Indices, class Scalar>
void
StandardWellEval<FluidSystem,Indices,Scalar>::
processFractions() const
{
    static constexpr int Gas = WellInterfaceIndices<FluidSystem,Indices,Scalar>::Gas;
    static constexpr int Oil = WellInterfaceIndices<FluidSystem,Indices,Scalar>::Oil;
    static constexpr int Water = WellInterfaceIndices<FluidSystem,Indices,Scalar>::Water;
    const auto pu = baseif_.phaseUsage();
    std::vector<double> F(baseif_.numPhases(), 0.0);

    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
        F[pu.phase_pos[Oil]] = 1.0;

        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            F[pu.phase_pos[Water]] = primary_variables_[WFrac];
            F[pu.phase_pos[Oil]] -= F[pu.phase_pos[Water]];
        }

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            F[pu.phase_pos[Gas]] = primary_variables_[GFrac];
            F[pu.phase_pos[Oil]] -= F[pu.phase_pos[Gas]];
        }
    }
    else if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
        F[pu.phase_pos[Water]] = 1.0;

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            F[pu.phase_pos[Gas]] = primary_variables_[GFrac];
            F[pu.phase_pos[Water]] -= F[pu.phase_pos[Gas]];
        }
    }
    else if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
        F[pu.phase_pos[Gas]] = 1.0;
    }

    [[maybe_unused]] double F_solvent;
    if constexpr (Indices::enableSolvent) {
        F_solvent = primary_variables_[SFrac];
        F[pu.phase_pos[Oil]] -= F_solvent;
    }

    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
        if (F[Water] < 0.0) {
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
        primary_variables_[WFrac] = F[pu.phase_pos[Water]];
    }

    if constexpr (has_gfrac_variable) {
        primary_variables_[GFrac] = F[pu.phase_pos[Gas]];
    }
    if constexpr (Indices::enableSolvent) {
        primary_variables_[SFrac] = F_solvent;
    }
}


template<class FluidSystem, class Indices, class Scalar>
void
StandardWellEval<FluidSystem,Indices,Scalar>::
updateThp(WellState& well_state,
          DeferredLogger& deferred_logger) const
{
    static constexpr int Gas = WellInterfaceIndices<FluidSystem,Indices,Scalar>::Gas;
    static constexpr int Oil = WellInterfaceIndices<FluidSystem,Indices,Scalar>::Oil;
    static constexpr int Water = WellInterfaceIndices<FluidSystem,Indices,Scalar>::Water;
    auto& ws = well_state.well(baseif_.indexOfWell());

    // When there is no vaild VFP table provided, we set the thp to be zero.
    if (!baseif_.isVFPActive(deferred_logger) || baseif_.wellIsStopped()) {
        ws.thp = 0;
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

    ws.thp = this->calculateThpFromBhp(well_state,
                                         rates,
                                         ws.bhp,
                                         deferred_logger);
}

template<class FluidSystem, class Indices, class Scalar>
void
StandardWellEval<FluidSystem,Indices,Scalar>::
updateWellStateFromPrimaryVariables(WellState& well_state,
                                    DeferredLogger& deferred_logger) const
{
    static constexpr int Gas = WellInterfaceIndices<FluidSystem,Indices,Scalar>::Gas;
    static constexpr int Oil = WellInterfaceIndices<FluidSystem,Indices,Scalar>::Oil;
    static constexpr int Water = WellInterfaceIndices<FluidSystem,Indices,Scalar>::Water;

    const PhaseUsage& pu = baseif_.phaseUsage();
    std::vector<double> F(baseif_.numPhases(), 0.0);
    [[maybe_unused]] double F_solvent = 0.0;
    if ( FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) ) {
        const int oil_pos = pu.phase_pos[Oil];
        F[oil_pos] = 1.0;

        if ( FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) ) {
            const int water_pos = pu.phase_pos[Water];
            F[water_pos] = primary_variables_[WFrac];
            F[oil_pos] -= F[water_pos];
        }

        if ( FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) ) {
            const int gas_pos = pu.phase_pos[Gas];
            F[gas_pos] = primary_variables_[GFrac];
            F[oil_pos] -= F[gas_pos];
        }

        if constexpr (Indices::enableSolvent) {
            F_solvent = primary_variables_[SFrac];
            F[oil_pos] -= F_solvent;
        }
    }
    else if ( FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) ) {
        const int water_pos = pu.phase_pos[Water];
        F[water_pos] = 1.0;

        if ( FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) ) {
            const int gas_pos = pu.phase_pos[Gas];
            F[gas_pos] = primary_variables_[GFrac];
            F[water_pos] -= F[gas_pos];
        }
    }
    else if ( FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) ) {
        const int gas_pos = pu.phase_pos[Gas];
        F[gas_pos] = 1.0;
    }


    // convert the fractions to be Q_p / G_total to calculate the phase rates
    for (int p = 0; p < baseif_.numPhases(); ++p) {
        const double scal = baseif_.scalingFactor(p);
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
    if constexpr (Indices::enableSolvent){
        F_solvent /= baseif_.scalingFactor(Indices::contiSolventEqIdx);
        F[pu.phase_pos[Gas]] += F_solvent;
    }

    auto& ws = well_state.well(baseif_.indexOfWell());
    ws.bhp = primary_variables_[Bhp]*bhp_scaling_;

    // calculate the phase rates based on the primary variables
    // for producers, this is not a problem, while not sure for injectors here
    if (baseif_.isProducer()) {
        const double g_total = primary_variables_[WQTotal];
        for (int p = 0; p < baseif_.numPhases(); ++p) {
            ws.surface_rates[p] = g_total * F[p];
        }
    } else { // injectors
        for (int p = 0; p < baseif_.numPhases(); ++p) {
            ws.surface_rates[p] = 0.0;
        }
        switch (baseif_.wellEcl().injectorType()) {
        case InjectorType::WATER:
            ws.surface_rates[pu.phase_pos[Water]] = primary_variables_[WQTotal];
            break;
        case InjectorType::GAS:
            ws.surface_rates[pu.phase_pos[Gas]] = primary_variables_[WQTotal];
            break;
        case InjectorType::OIL:
            ws.surface_rates[pu.phase_pos[Oil]] = primary_variables_[WQTotal];
            break;
        case InjectorType::MULTI:
            // Not supported.
            deferred_logger.warning("MULTI_PHASE_INJECTOR_NOT_SUPPORTED",
                                    "Multi phase injectors are not supported, requested for well " + baseif_.name());
            break;
        }
    }

    updateThp(well_state, deferred_logger);
}

template<class FluidSystem, class Indices, class Scalar>
void
StandardWellEval<FluidSystem,Indices,Scalar>::
updateWellStateFromPrimaryVariablesPolyMW(WellState& well_state) const
{
    if (baseif_.isInjector()) {
        auto& ws = well_state.well(baseif_.indexOfWell());
        auto& perf_data = ws.perf_data;
        auto& perf_water_velocity = perf_data.water_velocity;
        auto& perf_skin_pressure = perf_data.skin_pressure;
        for (int perf = 0; perf < baseif_.numPerfs(); ++perf) {
            perf_water_velocity[perf] = primary_variables_[Bhp + 1 + perf];
            perf_skin_pressure[perf] = primary_variables_[Bhp + 1 + baseif_.numPerfs() + perf];
        }
    }
}

template<class FluidSystem, class Indices, class Scalar>
void
StandardWellEval<FluidSystem,Indices,Scalar>::
computeAccumWell()
{
    for (int eq_idx = 0; eq_idx < numWellConservationEq; ++eq_idx) {
        F0_[eq_idx] = wellSurfaceVolumeFraction(eq_idx).value();
    }
}

template<class FluidSystem, class Indices, class Scalar>
void
StandardWellEval<FluidSystem,Indices,Scalar>::
updatePrimaryVariablesNewton(const BVectorWell& dwells,
                             [[maybe_unused]] const double dFLimit,
                             const double dBHPLimit) const
{
    const std::vector<double> old_primary_variables = primary_variables_;

    // for injectors, very typical one of the fractions will be one, and it is easy to get zero value
    // fractions. not sure what is the best way to handle it yet, so we just use 1.0 here
    [[maybe_unused]] const double relaxation_factor_fractions =
        (baseif_.isProducer()) ? relaxationFactorFractionsProducer(old_primary_variables, dwells)
                               : 1.0;

    // update the second and third well variable (The flux fractions)

    if constexpr (has_wfrac_variable) {
        const int sign2 = dwells[0][WFrac] > 0 ? 1: -1;
        const double dx2_limited = sign2 * std::min(std::abs(dwells[0][WFrac] * relaxation_factor_fractions), dFLimit);
        // primary_variables_[WFrac] = old_primary_variables[WFrac] - dx2_limited;
        primary_variables_[WFrac] = old_primary_variables[WFrac] - dx2_limited;
    }

    if constexpr (has_gfrac_variable) {
        const int sign3 = dwells[0][GFrac] > 0 ? 1: -1;
        const double dx3_limited = sign3 * std::min(std::abs(dwells[0][GFrac] * relaxation_factor_fractions), dFLimit);
        primary_variables_[GFrac] = old_primary_variables[GFrac] - dx3_limited;
    }

    if constexpr (Indices::enableSolvent) {
        const int sign4 = dwells[0][SFrac] > 0 ? 1: -1;
        const double dx4_limited = sign4 * std::min(std::abs(dwells[0][SFrac]) * relaxation_factor_fractions, dFLimit);
        primary_variables_[SFrac] = old_primary_variables[SFrac] - dx4_limited;
    }

    processFractions();

    // updating the total rates Q_t
    const double relaxation_factor_rate = this->relaxationFactorRate(old_primary_variables, dwells);
    primary_variables_[WQTotal] = old_primary_variables[WQTotal] - dwells[0][WQTotal] * relaxation_factor_rate;

    // updating the bottom hole pressure
    {
        const int sign1 = dwells[0][Bhp] > 0 ? 1: -1;
        const double dx1_limited = sign1 * std::min(std::abs(dwells[0][Bhp]), std::abs(old_primary_variables[Bhp]) * dBHPLimit/bhp_scaling_);
        // 1e5 to make sure bhp will not be below 1bar
        primary_variables_[Bhp] = std::max(old_primary_variables[Bhp] - dx1_limited, 1e5/bhp_scaling_);
    }
}

template<class FluidSystem, class Indices, class Scalar>
ConvergenceReport
StandardWellEval<FluidSystem,Indices,Scalar>::
getWellConvergence(const WellState& well_state,
                   const std::vector<double>& B_avg,
                   const double maxResidualAllowed,
                   const double tol_wells,
                   const double relaxed_tolerance_flow,
                   const bool relax_tolerance,
                   std::vector<double>& res,
                   DeferredLogger& deferred_logger) const
{
    res.resize(numWellEq_);
    for (int eq_idx = 0; eq_idx < numWellEq_; ++eq_idx) {
        // magnitude of the residual matters
        res[eq_idx] = std::abs(this->resWell_[0][eq_idx]);
    }

    std::vector<double> well_flux_residual(baseif_.numComponents());

    // Finish computation
    for (int compIdx = 0; compIdx < baseif_.numComponents(); ++compIdx )
    {
        well_flux_residual[compIdx] = B_avg[compIdx] * res[compIdx];
    }

    ConvergenceReport report;
    using CR = ConvergenceReport;
    CR::WellFailure::Type type = CR::WellFailure::Type::MassBalance;
    // checking if any NaN or too large residuals found
    for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
        if (!FluidSystem::phaseIsActive(phaseIdx)) {
            continue;
        }

        const unsigned canonicalCompIdx = FluidSystem::solventComponentIndex(phaseIdx);
        const int compIdx = Indices::canonicalToActiveComponentIndex(canonicalCompIdx);

        if (std::isnan(well_flux_residual[compIdx])) {
            report.setWellFailed({type, CR::Severity::NotANumber, compIdx, baseif_.name()});
        } else if (well_flux_residual[compIdx] > maxResidualAllowed) {
            report.setWellFailed({type, CR::Severity::TooLarge, compIdx, baseif_.name()});
        } else if (!relax_tolerance && well_flux_residual[compIdx] > tol_wells) {
            report.setWellFailed({type, CR::Severity::Normal, compIdx, baseif_.name()});
        } else if (well_flux_residual[compIdx] > relaxed_tolerance_flow) {
            report.setWellFailed({type, CR::Severity::Normal, compIdx, baseif_.name()});
        }
    }

    this->checkConvergenceControlEq(well_state, report, deferred_logger, maxResidualAllowed);

    return report;
}

template<class FluidSystem, class Indices, class Scalar>
void
StandardWellEval<FluidSystem,Indices,Scalar>::
computeConnectionDensities(const std::vector<double>& perfComponentRates,
                           const std::vector<double>& b_perf,
                           const std::vector<double>& rsmax_perf,
                           const std::vector<double>& rvmax_perf,
                           const std::vector<double>& rvwmax_perf,
                           const std::vector<double>& surf_dens_perf,
                           DeferredLogger& deferred_logger)
{
    // Verify that we have consistent input.
    const int nperf = baseif_.numPerfs();
    const int num_comp = baseif_.numComponents();

    // 1. Compute the flow (in surface volume units for each
    //    component) exiting up the wellbore from each perforation,
    //    taking into account flow from lower in the well, and
    //    in/out-flow at each perforation.
    std::vector<double> q_out_perf((nperf)*num_comp, 0.0);

    // Step 1 depends on the order of the perforations. Hence we need to
    // do the modifications globally.
    // Create and get the global perforation information and do this sequentially
    // on each process

    const auto& factory = baseif_.parallelWellInfo().getGlobalPerfContainerFactory();
    auto global_q_out_perf = factory.createGlobal(q_out_perf, num_comp);
    auto global_perf_comp_rates = factory.createGlobal(perfComponentRates, num_comp);

    // TODO: investigate whether we should use the following techniques to calcuate the composition of flows in the wellbore
    // Iterate over well perforations from bottom to top.
    for (int perf = factory.numGlobalPerfs() - 1; perf >= 0; --perf) {
        for (int component = 0; component < num_comp; ++component) {
            auto index = perf * num_comp + component;
            if (perf == factory.numGlobalPerfs() - 1) {
                // This is the bottom perforation. No flow from below.
                global_q_out_perf[index] = 0.0;
            } else {
                // Set equal to flow from below.
                global_q_out_perf[index] = global_q_out_perf[index + num_comp];
            }
            // Subtract outflow through perforation.
            global_q_out_perf[index] -= global_perf_comp_rates[index];
        }
    }

    // Copy the data back to local view
    factory.copyGlobalToLocal(global_q_out_perf, q_out_perf, num_comp);

    // 2. Compute the component mix at each perforation as the
    //    absolute values of the surface rates divided by their sum.
    //    Then compute volume ratios (formation factors) for each perforation.
    //    Finally compute densities for the segments associated with each perforation.
    std::vector<double> mix(num_comp,0.0);
    std::vector<double> x(num_comp);
    std::vector<double> surf_dens(num_comp);

    for (int perf = 0; perf < nperf; ++perf) {
        // Find component mix.
        const double tot_surf_rate = std::accumulate(q_out_perf.begin() + num_comp*perf,
                                                     q_out_perf.begin() + num_comp*(perf+1), 0.0);
        if (tot_surf_rate != 0.0) {
            for (int component = 0; component < num_comp; ++component) {
                mix[component] = std::fabs(q_out_perf[perf*num_comp + component]/tot_surf_rate);
            }
        } else if (num_comp == 1) {
            mix[num_comp-1] = 1.0;
        } else {
            std::fill(mix.begin(), mix.end(), 0.0);
            // No flow => use well specified fractions for mix.
            if (baseif_.isInjector()) {
                switch (baseif_.wellEcl().injectorType()) {
                case InjectorType::WATER:
                    mix[FluidSystem::waterCompIdx] = 1.0;
                    break;
                case InjectorType::GAS:
                    mix[FluidSystem::gasCompIdx] = 1.0;
                    break;
                case InjectorType::OIL:
                    mix[FluidSystem::oilCompIdx] = 1.0;
                    break;
                case InjectorType::MULTI:
                    // Not supported.
                    // deferred_logger.warning("MULTI_PHASE_INJECTOR_NOT_SUPPORTED",
                    //                         "Multi phase injectors are not supported, requested for well " + name());
                    break;
                }
            } else {
                assert(baseif_.isProducer());
                // For the frist perforation without flow we use the preferred phase to decide the mix initialization.
                if (perf == 0) { //
                    switch (baseif_.wellEcl().getPreferredPhase()) {
                    case Phase::OIL:
                        mix[FluidSystem::oilCompIdx] = 1.0;
                        break;
                    case Phase::GAS:
                        mix[FluidSystem::gasCompIdx] = 1.0;
                        break;
                    case Phase::WATER:
                        mix[FluidSystem::waterCompIdx] = 1.0;
                        break;
                    default:
                        // No others supported.
                        break;
                    }
                // For the rest of the perforation without flow we use mix from the above perforation.
                } else {
                    mix = x;
                }

            }
        }
        // Compute volume ratio.
        x = mix;

        // Subtract dissolved gas from oil phase and vapporized oil from gas phase and vaporized water from gas phase
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            const unsigned gaspos = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
            const unsigned oilpos = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
            double rs = 0.0;
            double rv = 0.0;
            if (!rsmax_perf.empty() && mix[oilpos] > 1e-12) {
                rs = std::min(mix[gaspos]/mix[oilpos], rsmax_perf[perf]);
            }
            if (!rvmax_perf.empty() && mix[gaspos] > 1e-12) {
                rv = std::min(mix[oilpos]/mix[gaspos], rvmax_perf[perf]);
            }
            const double d = 1.0 - rs*rv;
            if (d <= 0.0) {
                std::ostringstream sstr;
                sstr << "Problematic d value " << d << " obtained for well " << baseif_.name()
                     << " during ccomputeConnectionDensities with rs " << rs
                     << ", rv " << rv
                     << " obtaining d " << d
                     << " Continue as if no dissolution (rs = 0) and vaporization (rv = 0) "
                     << " for this connection.";
                deferred_logger.debug(sstr.str());
            } else {
                if (rs > 0.0) {
                    // Subtract gas in oil from gas mixture
                    x[gaspos] = (mix[gaspos] - mix[oilpos]*rs)/d;
                }
                if (rv > 0.0) {
                    // Subtract oil in gas from oil mixture
                    x[oilpos] = (mix[oilpos] - mix[gaspos]*rv)/d;
                }
            }
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                //matrix system: (mix[oilpos] = q_os, x[oilpos] = bo*q_or, etc...)
                //┌             ┐   ┌                ┐  ┌           ┐
                //│mix[oilpos]  │   | 1     Rv     0 |  |x[oilpos]  |
                //│mix[gaspos]  │ = │ Rs    1      0 │  │x[gaspos]  │
                //│mix[waterpos]│   │ 0     Rvw    1 │  │x[waterpos │
                //└             ┘   └                ┘  └           ┘
                const unsigned waterpos = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
                double rvw = 0.0;
                if (!rvwmax_perf.empty() && mix[gaspos] > 1e-12) {
                    rvw = std::min(mix[waterpos]/mix[gaspos], rvwmax_perf[perf]);
                }
                if (rvw > 0.0) {
                    // Subtract water in gas from water mixture
                    x[waterpos] = mix[waterpos] - x[gaspos] * rvw;
                }
            }
        } else if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            //no oil
            const unsigned gaspos = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
            const unsigned waterpos = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
            double rvw = 0.0;
            if (!rvwmax_perf.empty() && mix[gaspos] > 1e-12) {
                rvw = std::min(mix[waterpos]/mix[gaspos], rvwmax_perf[perf]);
            }
            if (rvw > 0.0) {
               // Subtract water in gas from water mixture
               x[waterpos] = mix[waterpos] - mix[gaspos] * rvw;
            }
        }

        double volrat = 0.0;
        for (int component = 0; component < num_comp; ++component) {
            volrat += x[component] / b_perf[perf*num_comp+ component];
        }
        for (int component = 0; component < num_comp; ++component) {
            surf_dens[component] = surf_dens_perf[perf*num_comp+ component];
        }

        // Compute segment density.
        this->perf_densities_[perf] = std::inner_product(surf_dens.begin(), surf_dens.end(), mix.begin(), 0.0) / volrat;
    }
}


template<class FluidSystem, class Indices, class Scalar>
void
StandardWellEval<FluidSystem,Indices,Scalar>::
init(std::vector<double>& perf_depth,
     const std::vector<double>& depth_arg,
     const int num_cells,
     const bool has_polymermw)
{
    perf_depth.resize(baseif_.numPerfs(), 0.);
    for (int perf = 0; perf < baseif_.numPerfs(); ++perf) {
        const int cell_idx = baseif_.cells()[perf];
        perf_depth[perf] = depth_arg[cell_idx];
    }

    // counting/updating primary variable numbers
    if (has_polymermw) {
        if (baseif_.isInjector()) {
            // adding a primary variable for water perforation rate per connection
            numWellEq_ += baseif_.numPerfs();
            // adding a primary variable for skin pressure per connection
            numWellEq_ += baseif_.numPerfs();
        }
    }

    // with the updated numWellEq_, we can initialize the primary variables and matrices now
    primary_variables_.resize(numWellEq_, 0.0);
    primary_variables_evaluation_.resize(numWellEq_, EvalWell{numWellEq_ + Indices::numEq, 0.0});

    // setup sparsity pattern for the matrices
    //[A C^T    [x    =  [ res
    // B D] x_well]      res_well]
    // set the size of the matrices
    this->duneD_.setSize(1, 1, 1);
    this->duneB_.setSize(1, num_cells, baseif_.numPerfs());
    this->duneC_.setSize(1, num_cells, baseif_.numPerfs());

    for (auto row=this->duneD_.createbegin(), end = this->duneD_.createend(); row!=end; ++row) {
        // Add nonzeros for diagonal
        row.insert(row.index());
    }
    // the block size is run-time determined now
    this->duneD_[0][0].resize(numWellEq_, numWellEq_);

    for (auto row = this->duneB_.createbegin(), end = this->duneB_.createend(); row!=end; ++row) {
        for (int perf = 0 ; perf < baseif_.numPerfs(); ++perf) {
            const int cell_idx = baseif_.cells()[perf];
            row.insert(cell_idx);
        }
    }

    for (int perf = 0 ; perf < baseif_.numPerfs(); ++perf) {
        const int cell_idx = baseif_.cells()[perf];
         // the block size is run-time determined now
         this->duneB_[0][cell_idx].resize(numWellEq_, Indices::numEq);
    }

    // make the C^T matrix
    for (auto row = this->duneC_.createbegin(), end = this->duneC_.createend(); row!=end; ++row) {
        for (int perf = 0; perf < baseif_.numPerfs(); ++perf) {
            const int cell_idx = baseif_.cells()[perf];
            row.insert(cell_idx);
        }
    }

    for (int perf = 0; perf < baseif_.numPerfs(); ++perf) {
        const int cell_idx = baseif_.cells()[perf];
        this->duneC_[0][cell_idx].resize(numWellEq_, Indices::numEq);
    }

    this->resWell_.resize(1);
    // the block size of resWell_ is also run-time determined now
    this->resWell_[0].resize(numWellEq_);

    // resize temporary class variables
    this->Bx_.resize( this->duneB_.N() );
    for (unsigned i = 0; i < this->duneB_.N(); ++i) {
        this->Bx_[i].resize(numWellEq_);
    }

    this->invDrw_.resize( this->duneD_.N() );
    for (unsigned i = 0; i < this->duneD_.N(); ++i) {
        this->invDrw_[i].resize(numWellEq_);
    }
}

#if HAVE_CUDA || HAVE_OPENCL
template<class FluidSystem, class Indices, class Scalar>
void
StandardWellEval<FluidSystem,Indices,Scalar>::
addWellContribution(WellContributions& wellContribs) const
{
    std::vector<int> colIndices;
    std::vector<double> nnzValues;
    colIndices.reserve(this->duneB_.nonzeroes());
    nnzValues.reserve(this->duneB_.nonzeroes()*numStaticWellEq * Indices::numEq);

    // duneC
    for ( auto colC = this->duneC_[0].begin(), endC = this->duneC_[0].end(); colC != endC; ++colC )
    {
        colIndices.emplace_back(colC.index());
        for (int i = 0; i < numStaticWellEq; ++i) {
            for (int j = 0; j < Indices::numEq; ++j) {
                nnzValues.emplace_back((*colC)[i][j]);
            }
        }
    }
    wellContribs.addMatrix(WellContributions::MatrixType::C, colIndices.data(), nnzValues.data(), this->duneC_.nonzeroes());

    // invDuneD
    colIndices.clear();
    nnzValues.clear();
    colIndices.emplace_back(0);
    for (int i = 0; i < numStaticWellEq; ++i)
    {
        for (int j = 0; j < numStaticWellEq; ++j) {
            nnzValues.emplace_back(this->invDuneD_[0][0][i][j]);
        }
    }
    wellContribs.addMatrix(WellContributions::MatrixType::D, colIndices.data(), nnzValues.data(), 1);

    // duneB
    colIndices.clear();
    nnzValues.clear();
    for ( auto colB = this->duneB_[0].begin(), endB = this->duneB_[0].end(); colB != endB; ++colB )
    {
        colIndices.emplace_back(colB.index());
        for (int i = 0; i < numStaticWellEq; ++i) {
            for (int j = 0; j < Indices::numEq; ++j) {
                nnzValues.emplace_back((*colB)[i][j]);
            }
        }
    }
    wellContribs.addMatrix(WellContributions::MatrixType::B, colIndices.data(), nnzValues.data(), this->duneB_.nonzeroes());
}
#endif

#define INSTANCE(A,...) \
template class StandardWellEval<BlackOilFluidSystem<double,A>,__VA_ARGS__,double>;

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
INSTANCE(BlackOilDefaultIndexTraits,BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,true,0u,2u,0u>)
INSTANCE(BlackOilDefaultIndexTraits,BlackOilTwoPhaseIndices<0u,0u,2u,0u,false,false,0u,2u,0u>)
INSTANCE(BlackOilDefaultIndexTraits,BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(BlackOilDefaultIndexTraits,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,0u,0u>)

// Blackoil
INSTANCE(BlackOilDefaultIndexTraits,BlackOilIndices<0u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilDefaultIndexTraits,BlackOilIndices<0u,0u,0u,0u,true,false,0u,0u>)
INSTANCE(BlackOilDefaultIndexTraits,BlackOilIndices<0u,0u,0u,0u,false,true,0u,0u>)
INSTANCE(BlackOilDefaultIndexTraits,BlackOilIndices<1u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilDefaultIndexTraits,BlackOilIndices<0u,1u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilDefaultIndexTraits,BlackOilIndices<0u,0u,1u,0u,false,false,0u,0u>)
INSTANCE(BlackOilDefaultIndexTraits,BlackOilIndices<0u,0u,0u,1u,false,false,0u,0u>)
INSTANCE(BlackOilDefaultIndexTraits,BlackOilIndices<0u,0u,0u,1u,false,false,1u,0u>)
INSTANCE(BlackOilDefaultIndexTraits,BlackOilIndices<0u,0u,0u,1u,false,true,0u,0u>)

}
