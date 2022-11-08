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
#include <opm/simulators/wells/WellBhpThpCalculator.hpp>
#include <opm/simulators/wells/WellConvergence.hpp>
#include <opm/simulators/wells/WellInterfaceIndices.hpp>
#include <opm/simulators/wells/WellState.hpp>
#include <opm/simulators/linalg/bda/WellContributions.hpp>

#include <cassert>
#include <cmath>



namespace Opm
{

template<class FluidSystem, class Indices, class Scalar>
StandardWellEval<FluidSystem,Indices,Scalar>::
StandardWellEval(const WellInterfaceIndices<FluidSystem,Indices,Scalar>& baseif)
    : StandardWellGeneric<Scalar>(baseif)
    , baseif_(baseif)
    , primary_variables_(baseif_)
    , F0_(numWellConservationEq)
    , linSys_(baseif_.parallelWellInfo())
{
}

template<class FluidSystem, class Indices, class Scalar>
typename StandardWellEval<FluidSystem,Indices,Scalar>::EvalWell
StandardWellEval<FluidSystem,Indices,Scalar>::
extendEval(const Eval& in) const
{
    EvalWell out(primary_variables_.numWellEq() + Indices::numEq, in.value());
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
        return inj_frac * primary_variables_.evaluation_[WQTotal];
    } else { // producers
        return primary_variables_.evaluation_[WQTotal] *
               primary_variables_.volumeFractionScaled(comp_idx);
    }
}

template<class FluidSystem, class Indices, class Scalar>
void
StandardWellEval<FluidSystem,Indices,Scalar>::
updatePrimaryVariables(const WellState& well_state, DeferredLogger& deferred_logger) const
{
    if (!baseif_.isOperableAndSolvable() && !baseif_.wellIsStopped())
        return;

    this->primary_variables_.update(well_state, deferred_logger);
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
            F[pu.phase_pos[Water]] = primary_variables_.value_[WFrac];
            F[pu.phase_pos[Oil]] -= F[pu.phase_pos[Water]];
        }

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            F[pu.phase_pos[Gas]] = primary_variables_.value_[GFrac];
            F[pu.phase_pos[Oil]] -= F[pu.phase_pos[Gas]];
        }
    }
    else if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
        F[pu.phase_pos[Water]] = 1.0;

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            F[pu.phase_pos[Gas]] = primary_variables_.value_[GFrac];
            F[pu.phase_pos[Water]] -= F[pu.phase_pos[Gas]];
        }
    }
    else if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
        F[pu.phase_pos[Gas]] = 1.0;
    }

    [[maybe_unused]] double F_solvent;
    if constexpr (Indices::enableSolvent) {
        F_solvent = primary_variables_.value_[SFrac];
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
        primary_variables_.value_[WFrac] = F[pu.phase_pos[Water]];
    }

    if constexpr (has_gfrac_variable) {
        primary_variables_.value_[GFrac] = F[pu.phase_pos[Gas]];
    }
    if constexpr (Indices::enableSolvent) {
        primary_variables_.value_[SFrac] = F_solvent;
    }
}

template<class FluidSystem, class Indices, class Scalar>
void
StandardWellEval<FluidSystem,Indices,Scalar>::
updateWellStateFromPrimaryVariables(WellState& well_state,
                                    DeferredLogger& deferred_logger) const
{
    this->primary_variables_.copyToWellState(well_state, deferred_logger);

    WellBhpThpCalculator(baseif_).
            updateThp(this->getRho(),
                      [this,&well_state]() { return this->baseif_.getALQ(well_state); },
                      {FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx),
                       FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx),
                       FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)},
                      well_state, deferred_logger);
}

template<class FluidSystem, class Indices, class Scalar>
void
StandardWellEval<FluidSystem,Indices,Scalar>::
computeAccumWell()
{
    for (size_t eq_idx = 0; eq_idx < F0_.size(); ++eq_idx) {
        F0_[eq_idx] = this->primary_variables_.surfaceVolumeFraction(eq_idx).value();
    }
}

template<class FluidSystem, class Indices, class Scalar>
void
StandardWellEval<FluidSystem,Indices,Scalar>::
updatePrimaryVariablesNewton(const BVectorWell& dwells,
                             [[maybe_unused]] const double dFLimit,
                             const double dBHPLimit) const
{
    const std::vector<double> old_primary_variables = primary_variables_.value_;

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
        primary_variables_.value_[WFrac] = old_primary_variables[WFrac] - dx2_limited;
    }

    if constexpr (has_gfrac_variable) {
        const int sign3 = dwells[0][GFrac] > 0 ? 1: -1;
        const double dx3_limited = sign3 * std::min(std::abs(dwells[0][GFrac] * relaxation_factor_fractions), dFLimit);
        primary_variables_.value_[GFrac] = old_primary_variables[GFrac] - dx3_limited;
    }

    if constexpr (Indices::enableSolvent) {
        const int sign4 = dwells[0][SFrac] > 0 ? 1: -1;
        const double dx4_limited = sign4 * std::min(std::abs(dwells[0][SFrac]) * relaxation_factor_fractions, dFLimit);
        primary_variables_.value_[SFrac] = old_primary_variables[SFrac] - dx4_limited;
    }

    processFractions();

    // updating the total rates Q_t
    const double relaxation_factor_rate = this->relaxationFactorRate(old_primary_variables, dwells[0][WQTotal]);
    primary_variables_.value_[WQTotal] = old_primary_variables[WQTotal] - dwells[0][WQTotal] * relaxation_factor_rate;

    // updating the bottom hole pressure
    {
        const int sign1 = dwells[0][Bhp] > 0 ? 1: -1;
        const double dx1_limited = sign1 * std::min(std::abs(dwells[0][Bhp]), std::abs(old_primary_variables[Bhp]) * dBHPLimit);
        // 1e5 to make sure bhp will not be below 1bar
        primary_variables_.value_[Bhp] = std::max(old_primary_variables[Bhp] - dx1_limited, 1e5);
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
    res.resize(this->primary_variables_.numWellEq());
    for (int eq_idx = 0; eq_idx < this->primary_variables_.numWellEq(); ++eq_idx) {
        // magnitude of the residual matters
        res[eq_idx] = std::abs(this->linSys_.residual()[0][eq_idx]);
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

    WellConvergence(baseif_).
        checkConvergenceControlEq(well_state,
                                  {1.e3, 1.e4, 1.e-4, 1.e-6, maxResidualAllowed},
                                  std::abs(this->linSys_.residual()[0][Bhp]),
                                  report,
                                  deferred_logger);

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
    int numWellEq = primary_variables_.numWellEq();
    if (has_polymermw) {
        if (baseif_.isInjector()) {
            // adding a primary variable for water perforation rate per connection
            numWellEq += baseif_.numPerfs();
            // adding a primary variable for skin pressure per connection
            numWellEq += baseif_.numPerfs();
        }
    }

    // with the updated numWellEq, we can initialize the primary variables and matrices now
    primary_variables_.resize(numWellEq);

    // setup sparsity pattern for the matrices
    this->linSys_.init(num_cells, numWellEq, baseif_.numPerfs(), baseif_.cells());
}

#define INSTANCE(...) \
template class StandardWellEval<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,__VA_ARGS__,double>;

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
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,true,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,2u,0u,false,false,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,0u,0u>)

// Blackoil
INSTANCE(BlackOilIndices<0u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,0u,true,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,0u,false,true,0u,0u>)
INSTANCE(BlackOilIndices<1u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,1u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,1u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,1u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,1u,false,false,1u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,1u,false,true,0u,0u>)

}
