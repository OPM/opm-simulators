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

#include <opm/material/densead/EvaluationFormat.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/models/blackoil/blackoilvariableandequationindices.hh>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/simulators/timestepping/ConvergenceReport.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>
#include <opm/simulators/wells/WellBhpThpCalculator.hpp>
#include <opm/simulators/wells/WellConvergence.hpp>
#include <opm/simulators/wells/WellInterfaceIndices.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <cmath>
#include <cstddef>

#include <fmt/format.h>

namespace Opm {

template<class FluidSystem, class Indices>
StandardWellEval<FluidSystem,Indices>::
StandardWellEval(const WellInterfaceIndices<FluidSystem,Indices>& baseif)
    : baseif_(baseif)
    , primary_variables_(baseif_)
    , F0_(numWellConservationEq)
    , linSys_(baseif_.parallelWellInfo())
    , connections_(baseif)
{
}

template<class FluidSystem, class Indices>
typename StandardWellEval<FluidSystem,Indices>::EvalWell
StandardWellEval<FluidSystem,Indices>::
extendEval(const Eval& in) const
{
    EvalWell out(in.value());
    // total number of equations/derivatives (well + reservoir)
    const int totalNumEq = primary_variables_.numWellEq() + Indices::numEq;
    for(int eqIdx = 0; eqIdx < Indices::numEq; ++eqIdx) {
        out.setDerivative(eqIdx, in.derivative(eqIdx), totalNumEq);
    }
    return out;
}

template<class FluidSystem, class Indices>
void
StandardWellEval<FluidSystem,Indices>::
computeAccumWell()
{
    for (std::size_t eq_idx = 0; eq_idx < F0_.size(); ++eq_idx) {
        F0_[eq_idx] = this->primary_variables_.surfaceVolumeFraction(eq_idx).value();
    }
}

template<class FluidSystem, class Indices>
void
StandardWellEval<FluidSystem,Indices>::
addBCDMatrix(std::vector<BMatrix>& b_matrices,
             std::vector<CMatrix>& c_matrices,
             std::vector<DMatrix>& d_matrices,
             std::vector<std::vector<int>>& wcells,
             std::vector<WVector>& residual) const
{
    const auto& srcB = linSys_.getB();
    const auto& srcC = linSys_.getC();
    const auto& srcD = linSys_.getD();
    const size_t numPerfs = srcB.M();

    // Copy residual (DynamicVector -> FieldVector)
    WVector res(1);
    for (int i = 0; i < primary_variables_.numWellEq(); ++i) {
        res[0][i] = linSys_.residual()[0][i];
    }
    residual.push_back(res);

    // Copy D matrix (1x1 block, DynamicMatrix -> FieldMatrix)
    DMatrix duneD;
    duneD.setSize(1, 1, 1);
    for (auto row = duneD.createbegin(); row != duneD.createend(); ++row) {
        // Add nonzeros for diagonal
        row.insert(row.index());
    }
    {
        auto& dest = duneD[0][0];
        const auto& src = srcD[0][0];
        for (size_t i = 0; i < dest.N(); ++i) {
            for (size_t j = 0; j < dest.M(); ++j) {
                dest[i][j] = src[i][j];
            }
        }
    }

    // Copy B matrix (DynamicMatrix -> FieldMatrix)
    BMatrix duneB;
    assert(srcB.N() == 1);
    duneB.setSize(srcB.N(), srcB.M(), srcB.M());
    for (auto row = duneB.createbegin(); row != duneB.createend(); ++row) {
        for (size_t perf = 0; perf < numPerfs; ++perf) {
            row.insert(perf);
        }
    }
    for (size_t i = 0; i < srcB.M(); ++i) {
        const auto& src = srcB[0][i];
        auto& dest = duneB[0][i];
        for (size_t j = 0; j < src.N(); ++j) {
            for (size_t k = 0; k < src.M(); ++k) {
                dest[j][k] = src[j][k];
            }
        }
    }

    // Build C as transpose of internal C^T layout
    CMatrix duneC;
    duneC.setSize(srcC.M(), srcC.N(), srcC.M());
    for (auto row = duneC.createbegin(); row != duneC.createend(); ++row) {
        row.insert(0);  // Only one well dofs
    }
    for (size_t perf = 0; perf < numPerfs; ++perf) {
        auto& dest = duneC[perf][0];
        const auto& src = srcC[0][perf];
        for (size_t i = 0; i < src.N(); ++i) {
            for (size_t j = 0; j < src.M(); ++j) {
                dest[j][i] = src[i][j];
            }
        }
    }

    b_matrices.push_back(duneB);
    c_matrices.push_back(duneC);
    d_matrices.push_back(duneD);
    wcells.push_back(linSys_.cells());
}

template<class FluidSystem, class Indices>
ConvergenceReport
StandardWellEval<FluidSystem,Indices>::
getWellConvergence(const WellState<Scalar, IndexTraits>& well_state,
                   const std::vector<Scalar>& B_avg,
                   const Scalar maxResidualAllowed,
                   const Scalar tol_wells,
                   const Scalar relaxed_tolerance_flow,
                   const bool relax_tolerance,
                   const bool well_is_stopped,
                   std::vector<Scalar>& res,
                   DeferredLogger& deferred_logger) const
{
    res.resize(this->primary_variables_.numWellEq());
    for (int eq_idx = 0; eq_idx < this->primary_variables_.numWellEq(); ++eq_idx) {
        // magnitude of the residual matters
        res[eq_idx] = std::abs(this->linSys_.residual()[0][eq_idx]);
    }

    std::vector<Scalar> well_flux_residual(baseif_.numConservationQuantities());

    // Finish computation
    for (int compIdx = 0; compIdx < baseif_.numConservationQuantities(); ++compIdx) {
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
        const int compIdx = FluidSystem::canonicalToActiveCompIdx(canonicalCompIdx);

        if (std::isnan(well_flux_residual[compIdx])) {
            report.setWellFailed({type, CR::Severity::NotANumber, compIdx, baseif_.name()});
            report.setWellConvergenceMetric(type, CR::Severity::NotANumber, compIdx, well_flux_residual[compIdx], baseif_.name());
        } else if (well_flux_residual[compIdx] > maxResidualAllowed) {
            report.setWellFailed({type, CR::Severity::TooLarge, compIdx, baseif_.name()});
            report.setWellConvergenceMetric(type, CR::Severity::TooLarge, compIdx, well_flux_residual[compIdx], baseif_.name());
        } else if (!relax_tolerance && well_flux_residual[compIdx] > tol_wells) {
            report.setWellFailed({type, CR::Severity::Normal, compIdx, baseif_.name()});
            report.setWellConvergenceMetric(type, CR::Severity::Normal, compIdx, well_flux_residual[compIdx], baseif_.name());
        } else if (well_flux_residual[compIdx] > relaxed_tolerance_flow) {
            report.setWellFailed({type, CR::Severity::Normal, compIdx, baseif_.name()});
            report.setWellConvergenceMetric(type, CR::Severity::Normal, compIdx, well_flux_residual[compIdx], baseif_.name());
        } else {
            report.setWellConvergenceMetric(CR::WellFailure::Type::Invalid, CR::Severity::None, compIdx, well_flux_residual[compIdx], baseif_.name());
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
        const Scalar sign = baseif_.isProducer() ? -1. : 1.;
        const auto weight_total_flux = this->primary_variables_.value(PrimaryVariables::WQTotal) * sign;
        constexpr int dummy_phase = -1;
        if (weight_total_flux < 0.) {
            report.setWellFailed(
                    {CR::WellFailure::Type::WrongFlowDirection, CR::Severity::Normal, dummy_phase, baseif_.name()});
        }
    }

    return report;
}

template<class FluidSystem, class Indices>
void
StandardWellEval<FluidSystem,Indices>::
init(std::vector<Scalar>& perf_depth,
     const std::vector<Scalar>& depth_arg,
     const bool has_polymermw)
{
    perf_depth.resize(baseif_.numLocalPerfs(), 0.);
    for (int perf = 0; perf < baseif_.numLocalPerfs(); ++perf) {
        const int cell_idx = baseif_.cells()[perf];
        perf_depth[perf] = depth_arg[cell_idx];
    }

    // counting/updating primary variable numbers
    int numWellEq = primary_variables_.numWellEq();
    if (has_polymermw) {
        if (baseif_.isInjector()) {
            // adding a primary variable for water perforation rate per connection
            numWellEq += baseif_.numLocalPerfs();
            // adding a primary variable for skin pressure per connection
            numWellEq += baseif_.numLocalPerfs();
        }
    }

    // with the updated numWellEq, we can initialize the primary variables and matrices now
    primary_variables_.resize(numWellEq);

    // setup sparsity pattern for the matrices
    this->linSys_.init(numWellEq, baseif_.numLocalPerfs(), baseif_.cells());
}

#include <opm/simulators/utils/InstantiationIndicesMacros.hpp>

INSTANTIATE_TYPE_INDICES(StandardWellEval, double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE_INDICES(StandardWellEval, float)
#endif

}
