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

#include <opm/simulators/wells/StandardWellEquations.hpp>
#include <opm/simulators/wells/StandardWellPrimaryVariables.hpp>
#include <opm/simulators/wells/WellAssemble.hpp>
#include <opm/simulators/wells/WellBhpThpCalculator.hpp>
#include <opm/simulators/wells/WellInterfaceFluidSystem.hpp>
#include <opm/simulators/wells/WellState.hpp>

namespace Opm {

//! \brief Class administering assembler access to equation system.
template<class Scalar, int numEq>
class StandardWellEquationAccess {
public:
    //! \brief Constructor initializes reference to the equation system.
    StandardWellEquationAccess(StandardWellEquations<Scalar,numEq>& eqns)
        : eqns_(eqns)
    {}

    using BVectorWell = typename StandardWellEquations<Scalar,numEq>::BVectorWell;
    using DiagMatWell = typename StandardWellEquations<Scalar,numEq>::DiagMatWell;
    using OffDiatMatWell = typename StandardWellEquations<Scalar,numEq>::OffDiagMatWell;

    //! \brief Returns a reference to residual vector.
    BVectorWell& residual()
    {
        return eqns_.resWell_;
    }

    //! \brief Returns a reference to B matrix.
    OffDiatMatWell& B()
    {
        return eqns_.duneB_;
    }

    //! \brief Returns a reference to C matrix.
    OffDiatMatWell& C()
    {
        return eqns_.duneC_;
    }

    //! \brief Returns a reference to D matrix.
    DiagMatWell& D()
    {
        return eqns_.duneD_;
    }

private:
    StandardWellEquations<Scalar,numEq>& eqns_; //!< Reference to equation system
};

template<class FluidSystem, class Indices, class Scalar>
void
StandardWellAssemble<FluidSystem,Indices,Scalar>::
assembleControlEq(const WellState& well_state,
                  const GroupState& group_state,
                  const Schedule& schedule,
                  const SummaryState& summaryState,
                  const Well::InjectionControls& inj_controls,
                  const Well::ProductionControls& prod_controls,
                  const PrimaryVariables& primary_variables,
                  const double rho,
                  StandardWellEquations<Scalar,Indices::numEq>& eqns1,
                  DeferredLogger& deferred_logger) const
{
    static constexpr int Water = BlackoilPhases::Aqua;
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Gas = BlackoilPhases::Vapour;
    EvalWell control_eq(primary_variables.numWellEq() + Indices::numEq, 0.0);

    const auto& well = well_.wellEcl();

    auto getRates = [&]() {
        std::vector<EvalWell> rates(3, EvalWell(primary_variables.numWellEq() + Indices::numEq, 0.0));
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            rates[Water] = primary_variables.getQs(Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx));
        }
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            rates[Oil] = primary_variables.getQs(Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx));
        }
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            rates[Gas] = primary_variables.getQs(Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx));
        }
        if constexpr (Indices::enableSolvent) {
            rates[Gas] += primary_variables.getQs(Indices::contiSolventEqIdx);
        }
        return rates;
    };

    if (well_.stopppedOrZeroRateTarget(summaryState, well_state)) {
        control_eq = primary_variables.eval(PrimaryVariables::WQTotal);
    } else if (well_.isInjector()) {
        // Find injection rate.
        const EvalWell injection_rate = primary_variables.eval(PrimaryVariables::WQTotal);
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

        WellAssemble(well_).
            assembleControlEqInj(well_state,
                                 group_state,
                                 schedule,
                                 summaryState,
                                 inj_controls,
                                 primary_variables.eval(PrimaryVariables::Bhp),
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
        WellAssemble(well_).
            assembleControlEqProd(well_state,
                                  group_state,
                                  schedule,
                                  summaryState,
                                  prod_controls,
                                  primary_variables.eval(PrimaryVariables::Bhp),
                                  rates,
                                  bhp_from_thp,
                                  control_eq,
                                  deferred_logger);
    }

    // using control_eq to update the matrix and residuals
    // TODO: we should use a different index system for the well equations
    StandardWellEquationAccess eqns(eqns1);
    eqns.residual()[0][PrimaryVariables::Bhp] = control_eq.value();
    for (int pv_idx = 0; pv_idx < primary_variables.numWellEq(); ++pv_idx) {
        eqns.D()[0][0][PrimaryVariables::Bhp][pv_idx] = control_eq.derivative(pv_idx + Indices::numEq);
    }
}

template<class FluidSystem, class Indices, class Scalar>
void StandardWellAssemble<FluidSystem,Indices,Scalar>::
assembleInjectivityEq(const EvalWell& eq_pskin,
                      const EvalWell& eq_wat_vel,
                      const int pskin_index,
                      const int wat_vel_index,
                      const int cell_idx,
                      const int numWellEq,
                      StandardWellEquations<Scalar,Indices::numEq>& eqns1) const
{
    StandardWellEquationAccess eqns(eqns1);
    eqns.residual()[0][pskin_index] = eq_pskin.value();
    eqns.residual()[0][wat_vel_index] = eq_wat_vel.value();
    for (int pvIdx = 0; pvIdx < numWellEq; ++pvIdx) {
        eqns.D()[0][0][wat_vel_index][pvIdx] = eq_wat_vel.derivative(pvIdx+Indices::numEq);
        eqns.D()[0][0][pskin_index][pvIdx] = eq_pskin.derivative(pvIdx+Indices::numEq);
    }

         // the water velocity is impacted by the reservoir primary varaibles. It needs to enter matrix B
    for (int pvIdx = 0; pvIdx < Indices::numEq; ++pvIdx) {
        eqns.B()[0][cell_idx][wat_vel_index][pvIdx] = eq_wat_vel.derivative(pvIdx);
    }
}

template<class FluidSystem, class Indices, class Scalar>
void StandardWellAssemble<FluidSystem,Indices,Scalar>::
assemblePerforationEq(const EvalWell& cq_s_effective,
                      const int componentIdx,
                      const int cell_idx,
                      const int numWellEq,
                      StandardWellEquations<Scalar,Indices::numEq>& eqns1) const
{
    StandardWellEquationAccess eqns(eqns1);

    // subtract sum of phase fluxes in the well equations.
    eqns.residual()[0][componentIdx] += cq_s_effective.value();

       // assemble the jacobians
    for (int pvIdx = 0; pvIdx < numWellEq; ++pvIdx) {
        // also need to consider the efficiency factor when manipulating the jacobians.
        eqns.C()[0][cell_idx][pvIdx][componentIdx] -= cq_s_effective.derivative(pvIdx+Indices::numEq); // intput in transformed matrix
        eqns.D()[0][0][componentIdx][pvIdx] += cq_s_effective.derivative(pvIdx+Indices::numEq);
    }

    for (int pvIdx = 0; pvIdx < Indices::numEq; ++pvIdx) {
        eqns.B()[0][cell_idx][componentIdx][pvIdx] += cq_s_effective.derivative(pvIdx);
    }
}

template<class FluidSystem, class Indices, class Scalar>
void StandardWellAssemble<FluidSystem,Indices,Scalar>::
assembleSourceEq(const EvalWell& resWell_loc,
                 const int componentIdx,
                 const int numWellEq,
                 StandardWellEquations<Scalar,Indices::numEq>& eqns1) const
{
    StandardWellEquationAccess eqns(eqns1);
    for (int pvIdx = 0; pvIdx < numWellEq; ++pvIdx) {
        eqns.D()[0][0][componentIdx][pvIdx] += resWell_loc.derivative(pvIdx+Indices::numEq);
    }
    eqns.residual()[0][componentIdx] += resWell_loc.value();
}

template<class FluidSystem, class Indices, class Scalar>
void StandardWellAssemble<FluidSystem,Indices,Scalar>::
assembleZFracEq(const EvalWell& cq_s_zfrac_effective,
                const int cell_idx,
                const int numWellEq,
                StandardWellEquations<Scalar,Indices::numEq>& eqns1) const
{
    StandardWellEquationAccess eqns(eqns1);
    for (int pvIdx = 0; pvIdx < numWellEq; ++pvIdx) {
        eqns.C()[0][cell_idx][pvIdx][Indices::contiZfracEqIdx] -= cq_s_zfrac_effective.derivative(pvIdx+Indices::numEq);
    }
}

#define INSTANCE(Dim,...) \
template class StandardWellAssemble<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,__VA_ARGS__,double>;

// One phase
INSTANCE(4u, BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(5u, BlackOilOnePhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(9u, BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,5u>)

// Two phase
INSTANCE(6u, BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,0u,0u>)
INSTANCE(6u, BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(6u, BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,2u,0u>)
INSTANCE(7u, BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,false,0u,2u,0u>)
INSTANCE(7u, BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,true,0u,2u,0u>)
INSTANCE(7u, BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(7u, BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,0u,0u>)
INSTANCE(7u, BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,2u,0u>)
INSTANCE(8u, BlackOilTwoPhaseIndices<0u,0u,2u,0u,false,false,0u,2u,0u>)
INSTANCE(7u, BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,0u,0u>)
INSTANCE(8u, BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,true,0u,0u,0u>)
INSTANCE(7u, BlackOilTwoPhaseIndices<1u,0u,0u,0u,false,false,0u,0u,0u>)
// Blackoil
INSTANCE(8u, BlackOilIndices<0u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(9u, BlackOilIndices<0u,0u,0u,0u,true,false,0u,0u>)
INSTANCE(9u, BlackOilIndices<0u,0u,0u,0u,false,true,0u,0u>)
INSTANCE(9u, BlackOilIndices<0u,1u,0u,0u,false,false,0u,0u>)
INSTANCE(9u, BlackOilIndices<0u,0u,1u,0u,false,false,0u,0u>)
INSTANCE(9u, BlackOilIndices<0u,0u,0u,1u,false,false,0u,0u>)
INSTANCE(10u, BlackOilIndices<1u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(10u, BlackOilIndices<0u,0u,0u,1u,false,true,0u,0u>)
INSTANCE(10u, BlackOilIndices<0u,0u,0u,1u,false,false,1u,0u>)
INSTANCE(11u, BlackOilIndices<1u,0u,0u,0u,true,false,0u,0u>)

}
