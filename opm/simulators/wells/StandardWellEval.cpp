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
#include <cstddef>

namespace Opm
{

template<class FluidSystem, class Indices, class Scalar>
StandardWellEval<FluidSystem,Indices,Scalar>::
StandardWellEval(const WellInterfaceIndices<FluidSystem,Indices,Scalar>& baseif)
    : baseif_(baseif)
    , primary_variables_(baseif_)
    , F0_(numWellConservationEq)
    , linSys_(baseif_.parallelWellInfo())
    , connections_(baseif)
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
void
StandardWellEval<FluidSystem,Indices,Scalar>::
updateWellStateFromPrimaryVariables(const bool stop_or_zero_rate_target,
                                    WellState& well_state,
                                    const SummaryState& summary_state,
                                    DeferredLogger& deferred_logger) const
{
    this->primary_variables_.copyToWellState(well_state, deferred_logger);

    WellBhpThpCalculator(baseif_).
            updateThp(connections_.rho(),
                      stop_or_zero_rate_target,
                      [this,&well_state]() { return this->baseif_.getALQ(well_state); },
                      {FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx),
                       FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx),
                       FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)},
                      well_state, summary_state, deferred_logger);
}

template<class FluidSystem, class Indices, class Scalar>
void
StandardWellEval<FluidSystem,Indices,Scalar>::
computeAccumWell()
{
    for (std::size_t eq_idx = 0; eq_idx < F0_.size(); ++eq_idx) {
        F0_[eq_idx] = this->primary_variables_.surfaceVolumeFraction(eq_idx).value();
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
                   const bool well_is_stopped, 
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
                                  well_is_stopped, 
                                  report,
                                  deferred_logger);

    // for stopped well, we do not enforce the following checking to avoid dealing with sign of near-zero values
    // for BHP or THP controlled wells, we need to make sure the flow direction is correct
    if (!well_is_stopped && baseif_.isPressureControlled(well_state)) {
        // checking the flow direction
        const double sign = baseif_.isProducer() ? -1. : 1.;
        const auto weight_total_flux = this->primary_variables_.value(PrimaryVariables::WQTotal) * sign;
        constexpr int dummy_phase = -1;
        if (weight_total_flux < 0.) {
            report.setWellFailed(
                    {CR::WellFailure::Type::WrongFlowDirection, CR::Severity::Normal, dummy_phase, baseif_.name()});
        }
    }

    return report;
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
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,0u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,true,0u,0u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<1u,0u,0u,0u,false,false,0u,0u,0u>)
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

INSTANCE(BlackOilIndices<1u,0u,0u,0u,true,false,0u,0u>)

}
