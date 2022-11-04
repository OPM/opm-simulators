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
#include <opm/simulators/wells/StandardWellAssemble.hpp>

#include <opm/core/props/BlackoilPhases.hpp>

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/models/blackoil/blackoilindices.hh>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/simulators/linalg/matrixblock.hh>

#include <opm/simulators/wells/WellAssemble.hpp>
#include <opm/simulators/wells/WellBhpThpCalculator.hpp>
#include <opm/simulators/wells/WellInterfaceFluidSystem.hpp>

namespace Opm {

template<class FluidSystem, class Indices, class Scalar>
template<class EvalWell>
void
StandardWellAssemble<FluidSystem,Indices,Scalar>::
assembleControlEq(const WellState& well_state,
                  const GroupState& group_state,
                  const Schedule& schedule,
                  const SummaryState& summaryState,
                  const int numWellEq,
                  const EvalWell& wqTotal,
                  const EvalWell& bhp,
                  const std::function<EvalWell(int)>& getQs,
                  const double rho,
                  StandardWellEquations<Indices,Scalar>& A,
                  DeferredLogger& deferred_logger) const
{
    static constexpr int Water = BlackoilPhases::Aqua;
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Gas = BlackoilPhases::Vapour;
    EvalWell control_eq(numWellEq + Indices::numEq, 0.0);

    const auto& well = well_.wellEcl();

    auto getRates = [&]() {
        std::vector<EvalWell> rates(3, EvalWell(numWellEq + Indices::numEq, 0.0));
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

    if (well_.wellIsStopped()) {
        control_eq = wqTotal;
    } else if (well_.isInjector()) {
        // Find injection rate.
        const EvalWell injection_rate = wqTotal;
        // Setup function for evaluation of BHP from THP (used only if needed).
        std::function<EvalWell()> bhp_from_thp = [&]() {
            const auto rates = getRates();
            return WellBhpThpCalculator(well_).calculateBhpFromThp(well_state,
                                                                   rates,
                                                                   well,
                                                                   summaryState,
                                                                   rho,
                                                                   deferred_logger);
        };

        // Call generic implementation.
        const auto& inj_controls = well.injectionControls(summaryState);
        WellAssemble(well_).
             assembleControlEqInj(well_state,
                                  group_state,
                                  schedule,
                                  summaryState,
                                  inj_controls,
                                  bhp,
                                  injection_rate,
                                  bhp_from_thp,
                                  control_eq,
                                  deferred_logger);
    } else {
        // Find rates.
        const auto rates = getRates();
        // Setup function for evaluation of BHP from THP (used only if needed).
        std::function<EvalWell()> bhp_from_thp = [&]() {
             return WellBhpThpCalculator(well_).calculateBhpFromThp(well_state,
                                                                    rates,
                                                                    well,
                                                                    summaryState,
                                                                    rho,
                                                                    deferred_logger);
        };
        // Call generic implementation.
        const auto& prod_controls = well.productionControls(summaryState);
        WellAssemble(well_).
            assembleControlEqProd(well_state,
                                  group_state,
                                  schedule,
                                  summaryState,
                                  prod_controls,
                                  bhp,
                                  rates,
                                  bhp_from_thp,
                                  control_eq,
                                  deferred_logger);
    }

    // using control_eq to update the matrix and residuals
    // TODO: we should use a different index system for the well equations
    A.resWell_[0][A.Bhp] = control_eq.value();
    for (int pv_idx = 0; pv_idx < numWellEq; ++pv_idx) {
        A.duneD_[0][0][A.Bhp][pv_idx] = control_eq.derivative(pv_idx + Indices::numEq);
    }
}

template<class FluidSystem, class Indices, class Scalar>
template<class EvalWell>
void StandardWellAssemble<FluidSystem,Indices,Scalar>::
assembleInjectivityEq(const EvalWell& eq_pskin,
                      const EvalWell& eq_wat_vel,
                      const int pskin_index,
                      const int wat_vel_index,
                      const int cell_idx,
                      const int numWellEq,
                      StandardWellEquations<Indices,Scalar>& eqns) const
{
    eqns.resWell_[0][pskin_index] = eq_pskin.value();
    eqns.resWell_[0][wat_vel_index] = eq_wat_vel.value();
    for (int pvIdx = 0; pvIdx < numWellEq; ++pvIdx) {
        eqns.duneD_[0][0][wat_vel_index][pvIdx] = eq_wat_vel.derivative(pvIdx+Indices::numEq);
        eqns.duneD_[0][0][pskin_index][pvIdx] = eq_pskin.derivative(pvIdx+Indices::numEq);
    }

    // the water velocity is impacted by the reservoir primary varaibles. It needs to enter matrix B
    for (int pvIdx = 0; pvIdx < Indices::numEq; ++pvIdx) {
        eqns.duneB_[0][cell_idx][wat_vel_index][pvIdx] = eq_wat_vel.derivative(pvIdx);
    }
}

template<class FluidSystem, class Indices, class Scalar>
template<class EvalWell>
void StandardWellAssemble<FluidSystem,Indices,Scalar>::
assemblePerforationEq(const EvalWell& cq_s_effective,
                      const int componentIdx,
                      const int cell_idx,
                      const int numWellEq,
                      StandardWellEquations<Indices,Scalar>& eqns) const
{
    // subtract sum of phase fluxes in the well equations.
    eqns.resWell_[0][componentIdx] += cq_s_effective.value();

    // assemble the jacobians
    for (int pvIdx = 0; pvIdx < numWellEq; ++pvIdx) {
        // also need to consider the efficiency factor when manipulating the jacobians.
        eqns.duneC_[0][cell_idx][pvIdx][componentIdx] -= cq_s_effective.derivative(pvIdx+Indices::numEq); // intput in transformed matrix
        eqns.duneD_[0][0][componentIdx][pvIdx] += cq_s_effective.derivative(pvIdx+Indices::numEq);
    }

    for (int pvIdx = 0; pvIdx < Indices::numEq; ++pvIdx) {
        eqns.duneB_[0][cell_idx][componentIdx][pvIdx] += cq_s_effective.derivative(pvIdx);
    }
}

template<class FluidSystem, class Indices, class Scalar>
template<class EvalWell>
void StandardWellAssemble<FluidSystem,Indices,Scalar>::
assembleZFracEq(const EvalWell& cq_s_zfrac_effective,
                const int cell_idx,
                const int numWellEq,
                StandardWellEquations<Indices,Scalar>& eqns) const
{
    for (int pvIdx = 0; pvIdx < numWellEq; ++pvIdx) {
        eqns.duneC_[0][cell_idx][pvIdx][Indices::contiZfracEqIdx] -= cq_s_zfrac_effective.derivative(pvIdx+Indices::numEq);
    }
}

template<class FluidSystem, class Indices, class Scalar>
template<class EvalWell>
void StandardWellAssemble<FluidSystem,Indices,Scalar>::
assembleSourceEq(const EvalWell& resWell_loc,
                 const int componentIdx,
                 const int numWellEq,
                 StandardWellEquations<Indices,Scalar>& eqns) const
{
    for (int pvIdx = 0; pvIdx < numWellEq; ++pvIdx) {
        eqns.duneD_[0][0][componentIdx][pvIdx] += resWell_loc.derivative(pvIdx+Indices::numEq);
    }
    eqns.resWell_[0][componentIdx] += resWell_loc.value();
}

template<class FluidSystem, class Indices, class Scalar>
template<class PressureMatrix, class BVector>
void StandardWellAssemble<FluidSystem,Indices,Scalar>::
addWellPressureEqs(const BVector& weights,
                   const int pressureVarIndex,
                   const bool use_well_weights,
                   const WellState& well_state,
                   const int numWellEq,
                   const StandardWellEquations<Indices,Scalar>& eqns,
                   PressureMatrix& jacobian) const

{
    // This adds pressure quation for cpr
    // For use_well_weights=true
    //    weights lamda = inv(D)'v  v = 0 v(bhpInd) = 1
    //    the well equations are summed i lambda' B(:,pressureVarINd) -> B  lambda'*D(:,bhpInd) -> D
    // For use_well_weights = false
    //    weights lambda = \sum_i w /n where ths sum is over weights of all perforation cells
    //    in the case of pressure controlled trivial equations are used and bhp  C=B=0
    //    then the flow part of the well equations are summed lambda'*B(1:n,pressureVarInd) -> B lambda'*D(1:n,bhpInd) -> D
    // For bouth
    //    C -> w'C(:,bhpInd) where w is weights of the perforation cell

    // Add the well contributions in cpr
    // use_well_weights is a quasiimpes formulation which is not implemented in multisegment
    int bhp_var_index = StandardWellEquations<Indices,Scalar>::Bhp;
    int nperf = 0;
    auto cell_weights = weights[0];// not need for not(use_well_weights)
    cell_weights = 0.0;
    assert(eqns.duneC_.M() == weights.size());
    const int welldof_ind = eqns.duneC_.M() + well_.indexOfWell();
    // do not assume anything about pressure controlled with use_well_weights (work fine with the assumtion also)
    if (!well_.isPressureControlled(well_state) || use_well_weights) {
        // make coupling for reservoir to well
        for (auto colC = eqns.duneC_[0].begin(),
                  endC = eqns.duneC_[0].end(); colC != endC; ++colC) {
            const auto row_ind = colC.index();
            const auto& bw = weights[row_ind];
            double matel = 0;
            assert((*colC).M() == bw.size());
            for (size_t i = 0; i < bw.size(); ++i) {
                matel += (*colC)[bhp_var_index][i] * bw[i];
            }

            jacobian[row_ind][welldof_ind] = matel;
            cell_weights += bw;
            nperf += 1;
        }
    }
    cell_weights /= nperf;

    using BVectorWell = typename StandardWellEquations<Indices,Scalar>::BVectorWell;

    BVectorWell  bweights(1);
    size_t blockSz = numWellEq;
    bweights[0].resize(blockSz);
    bweights[0] = 0.0;
    double diagElem = 0;
    {
        if (use_well_weights) {
            // calculate weighs and set diagonal element
            //NB! use this options without treating pressure controlled separated
            //NB! calculate quasiimpes well weights NB do not work well with trueimpes reservoir weights
            double abs_max = 0;
            BVectorWell rhs(1);
            rhs[0].resize(blockSz);
            rhs[0][bhp_var_index] = 1.0;
            auto inv_diag_block = eqns.invDuneD_[0][0];
            auto inv_diag_block_transpose = Opm::wellhelpers::transposeDenseDynMatrix(inv_diag_block);
            for (size_t i = 0; i < blockSz; ++i) {
                bweights[0][i] = 0;
                for (size_t j = 0; j < blockSz; ++j) {
                    bweights[0][i] += inv_diag_block_transpose[i][j]*rhs[0][j];
                }
                abs_max = std::max(abs_max, std::fabs(bweights[0][i]));
            }
            assert( abs_max > 0.0 );
            for (size_t i = 0; i < blockSz; ++i) {
                bweights[0][i] /= abs_max;
            }
            diagElem = 1.0/abs_max;
        } else {
            // set diagonal element
            if (well_.isPressureControlled(well_state)) {
                bweights[0][blockSz-1] = 1.0;
                diagElem = 1.0;// better scaling could have used the calculation below if weights were calculated
            } else {
                for (size_t i = 0; i < cell_weights.size(); ++i) {
                    bweights[0][i] = cell_weights[i];
                }
                bweights[0][blockSz-1] = 0.0;
                diagElem = 0.0;
                const auto& locmat = eqns.duneD_[0][0];
                for (size_t i = 0; i < cell_weights.size(); ++i) {
                    diagElem += locmat[i][bhp_var_index]*cell_weights[i];
                }
            }
        }
    }
    //
    jacobian[welldof_ind][welldof_ind] = diagElem;
    // set the matrix elements for well reservoir coupling
    if (!well_.isPressureControlled(well_state) || use_well_weights) {
        for (auto colB = eqns.duneB_[0].begin(),
                  endB = eqns.duneB_[0].end(); colB != endB; ++colB) {
            const auto col_index = colB.index();
            const auto& bw = bweights[0];
            double matel = 0;
            for (size_t i = 0; i < bw.size(); ++i) {
                 matel += (*colB)[i][pressureVarIndex] * bw[i];
            }
            jacobian[welldof_ind][col_index] = matel;
        }
    }
}

#define INSTANCE(Dim,Block,...) \
template class StandardWellAssemble<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,__VA_ARGS__,double>; \
template void \
StandardWellAssemble<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,__VA_ARGS__,double>:: \
assembleControlEq(const WellState&, \
                  const GroupState&, \
                  const Schedule&, \
                  const SummaryState&, \
                  const int, \
                  const DenseAd::Evaluation<double,-1,Dim>&, \
                  const DenseAd::Evaluation<double,-1,Dim>&, \
                  const std::function<DenseAd::Evaluation<double,-1,Dim>(int)>&, \
                  const double, \
                  StandardWellEquations<__VA_ARGS__,double>&, \
                  DeferredLogger&) const; \
template void \
StandardWellAssemble<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,__VA_ARGS__,double>:: \
assembleInjectivityEq(const DenseAd::Evaluation<double,-1,Dim>&, \
                      const DenseAd::Evaluation<double,-1,Dim>&, \
                      const int, \
                      const int, \
                      const int, \
                      const int, \
                      StandardWellEquations<__VA_ARGS__,double>&) const; \
template void \
StandardWellAssemble<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,__VA_ARGS__,double>:: \
assemblePerforationEq(const DenseAd::Evaluation<double,-1,Dim>&, \
                      const int, \
                      const int, \
                      const int, \
                      StandardWellEquations<__VA_ARGS__,double>&) const; \
template void \
StandardWellAssemble<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,__VA_ARGS__,double>:: \
assembleZFracEq(const DenseAd::Evaluation<double,-1,Dim>&, \
                const int, \
                const int, \
                StandardWellEquations<__VA_ARGS__,double>&) const; \
template void \
StandardWellAssemble<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,__VA_ARGS__,double>:: \
assembleSourceEq(const DenseAd::Evaluation<double,-1,Dim>&, \
                      const int, \
                      const int, \
                      StandardWellEquations<__VA_ARGS__,double>&) const; \
template void \
StandardWellAssemble<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,__VA_ARGS__,double>:: \
addWellPressureEqs(const Dune::BlockVector<Dune::FieldVector<double,Block>>&, \
                   const int, \
                   const bool, \
                   const WellState&, \
                   const int, \
                   const StandardWellEquations<__VA_ARGS__,double>&, \
                   Dune::BCRSMatrix<MatrixBlock<double,1,1>>& jacobian) const;

// One phase
INSTANCE(4u, 1, BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(5u, 2, BlackOilOnePhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(9u, 6, BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,5u>)

// Two phase
INSTANCE(6u, 2, BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,0u,0u>)
INSTANCE(6u, 2, BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(6u, 2, BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,2u,0u>)
INSTANCE(7u, 3, BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,false,0u,2u,0u>)
INSTANCE(7u, 1, BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,true,0u,2u,0u>)
INSTANCE(7u, 3, BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(7u, 3, BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,0u,0u>)
INSTANCE(7u, 3, BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,2u,0u>)
INSTANCE(8u, 4, BlackOilTwoPhaseIndices<0u,0u,2u,0u,false,false,0u,2u,0u>)

// Blackoil
INSTANCE(8u, 3, BlackOilIndices<0u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(9u, 4, BlackOilIndices<0u,0u,0u,0u,true,false,0u,0u>)
INSTANCE(9u, 4, BlackOilIndices<0u,0u,0u,0u,false,true,0u,0u>)
INSTANCE(9u, 4, BlackOilIndices<0u,1u,0u,0u,false,false,0u,0u>)
INSTANCE(9u, 4, BlackOilIndices<0u,0u,1u,0u,false,false,0u,0u>)
INSTANCE(9u, 4, BlackOilIndices<0u,0u,0u,1u,false,false,0u,0u>)
INSTANCE(10u, 4, BlackOilIndices<1u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(10u, 5, BlackOilIndices<0u,0u,0u,1u,false,true,0u,0u>)
INSTANCE(10u, 1, BlackOilIndices<0u,0u,0u,1u,false,false,1u,0u>)

}
