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

#include <opm/simulators/wells/MultisegmentWellAssemble.hpp>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/models/blackoil/blackoilvariableandequationindices.hh>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/simulators/utils/BlackoilPhases.hpp>

#include <opm/simulators/wells/MultisegmentWellEquations.hpp>
#include <opm/simulators/wells/MultisegmentWellPrimaryVariables.hpp>
#include <opm/simulators/wells/WellAssemble.hpp>
#include <opm/simulators/wells/WellBhpThpCalculator.hpp>
#include <opm/simulators/wells/WellHelpers.hpp>
#include <opm/simulators/wells/WellInterfaceIndices.hpp>
#include <opm/simulators/wells/WellState.hpp>

namespace Opm {

//! \brief Class administering assembler access to equation system.
template<typename FluidSystem, typename Indices>
class MultisegmentWellEquationAccess {
public:
    using Scalar = typename FluidSystem::Scalar;
    //! \brief Constructor initializes reference to the equation system.
    explicit MultisegmentWellEquationAccess(MultisegmentWellEquations<FluidSystem, Indices>& eqns)
        : eqns_(eqns)
    {}

    using BVectorWell = typename MultisegmentWellEquations<FluidSystem, Indices>::BVectorWell;
    using DiagMatWell = typename MultisegmentWellEquations<FluidSystem, Indices>::DiagMatWell;
    using OffDiatMatWell = typename MultisegmentWellEquations<FluidSystem, Indices>::OffDiagMatWell;

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
    MultisegmentWellEquations<FluidSystem, Indices>& eqns_; //!< Reference to equation system
};

template<class FluidSystem, class Indices>
void MultisegmentWellAssemble<FluidSystem,Indices>::
assembleControlEq(const WellState<FluidSystem, Indices>& well_state,
                  const GroupState<Scalar>& group_state,
                  const Schedule& schedule,
                  const SummaryState& summaryState,
                  const Well::InjectionControls& inj_controls,
                  const Well::ProductionControls& prod_controls,
                  const Scalar rho,
                  const PrimaryVariables& primary_variables,
                  Equations& eqns1,
                  const bool stopped_or_zero_target,
                  DeferredLogger& deferred_logger) const
{
    /*
        This function assembles the control equation, similar as for StandardWells.
        It does *not* need communication.
    */
    static constexpr int Gas = BlackoilPhases::Vapour;
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Water = BlackoilPhases::Aqua;

    EvalWell control_eq(0.0);

    const auto& well = well_.wellEcl();

    auto getRates = [&]() {
        std::vector<EvalWell> rates(3, 0.0);
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            rates[Water] = primary_variables.getQs(Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx));
        }
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            rates[Oil] = primary_variables.getQs(Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx));
        }
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            rates[Gas] = primary_variables.getQs(Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx));
        }
        return rates;
    };

    if (stopped_or_zero_target) {
        control_eq = primary_variables.getWQTotal();
    } else if (well_.isInjector() ) {
        // Find scaling factor to get injection rate,
        const InjectorType injectorType = inj_controls.injector_type;
        Scalar scaling;
        const auto& pu = well_.phaseUsage();
        switch (injectorType) {
        case InjectorType::WATER:
        {
            scaling = well_.scalingFactor(pu.phase_pos[BlackoilPhases::Aqua]);
            break;
        }
        case InjectorType::OIL:
        {
            scaling = well_.scalingFactor(pu.phase_pos[BlackoilPhases::Liquid]);
            break;
        }
        case InjectorType::GAS:
        {
            scaling = well_.scalingFactor(pu.phase_pos[BlackoilPhases::Vapour]);
            break;
        }
        default:
            throw("Expected WATER, OIL or GAS as type for injectors " + well.name());
        }
        const EvalWell injection_rate = primary_variables.getWQTotal() / scaling;
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
        WellAssemble(well_).assembleControlEqInj(well_state,
                                                 group_state,
                                                 schedule,
                                                 summaryState,
                                                 inj_controls,
                                                 primary_variables.getBhp(),
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
        WellAssemble(well_).assembleControlEqProd(well_state,
                                                  group_state,
                                                  schedule,
                                                  summaryState,
                                                  prod_controls,
                                                  primary_variables.getBhp(),
                                                  rates,
                                                  bhp_from_thp,
                                                  control_eq,
                                                  deferred_logger);
    }

    MultisegmentWellEquationAccess<FluidSystem, Indices> eqns(eqns1);
    // using control_eq to update the matrix and residuals
    eqns.residual()[0][SPres] = control_eq.value();
    for (int pv_idx = 0; pv_idx < numWellEq; ++pv_idx) {
        eqns.D()[0][0][SPres][pv_idx] = control_eq.derivative(pv_idx + Indices::numEq);
    }
}

template<class FluidSystem, class Indices>
void MultisegmentWellAssemble<FluidSystem,Indices>::
assembleAccelerationTerm(const int seg_target,
                         const int seg,
                         const int seg_upwind,
                         const EvalWell& accelerationTerm,
                         Equations& eqns1) const
{   // seg_target:  segment for which we are assembling the acc term
    // seg:         segment for wich we have computed the term
    // seg_upwind:  upwind segment to seg
    // acceleration term shold be
    //  * velocity head for seg_target if seg = seg_target
    //  * negative velocity head for seg if seg != seg_target   
    
    /*
        This method is called in MultisegmentWellEval::assembleAccelerationPressureLoss.
        It does *not* need communication.
    */

    MultisegmentWellEquationAccess<FluidSystem, Indices> eqns(eqns1);
    eqns.residual()[seg_target][SPres] -= accelerationTerm.value();
    eqns.D()[seg_target][seg][SPres][SPres] -= accelerationTerm.derivative(SPres + Indices::numEq);
    eqns.D()[seg_target][seg][SPres][WQTotal] -= accelerationTerm.derivative(WQTotal + Indices::numEq);
    if constexpr (has_wfrac_variable) {
        eqns.D()[seg_target][seg_upwind][SPres][WFrac] -= accelerationTerm.derivative(WFrac + Indices::numEq);
    }
    if constexpr (has_gfrac_variable) {
        eqns.D()[seg_target][seg_upwind][SPres][GFrac] -= accelerationTerm.derivative(GFrac + Indices::numEq);
    }
}                     

template<class FluidSystem, class Indices>
void MultisegmentWellAssemble<FluidSystem,Indices>::
assembleHydroPressureLoss(const int seg,
                          const int seg_density,
                          const EvalWell& hydro_pressure_drop_seg,
                          Equations& eqns1) const
{
    /*
        This method is called in MultisegmentWellEval::assembleAccelerationAndHydroPressureLosses.
        It does *not* need communication.
    */

    MultisegmentWellEquationAccess<FluidSystem, Indices> eqns(eqns1);
    eqns.residual()[seg][SPres] -= hydro_pressure_drop_seg.value();
    for (int pv_idx = 0; pv_idx < numWellEq; ++pv_idx) {
        eqns.D()[seg][seg_density][SPres][pv_idx] -= hydro_pressure_drop_seg.derivative(pv_idx + Indices::numEq);
    }

}

template<class FluidSystem, class Indices>
void MultisegmentWellAssemble<FluidSystem,Indices>::
assemblePressureEqExtraDerivatives(const int seg,
                                   const int seg_upwind,
                                   const EvalWell& extra_derivatives,
                                   Equations& eqns1) const
{
    /*
        This method does *not* need communication.
    */
    MultisegmentWellEquationAccess<FluidSystem, Indices> eqns(eqns1);
    // disregard residual
    // Frac - derivatives are zero (they belong to upwind^2)
    eqns.D()[seg][seg_upwind][SPres][SPres] += extra_derivatives.derivative(SPres + Indices::numEq);
    eqns.D()[seg][seg_upwind][SPres][WQTotal] += extra_derivatives.derivative(WQTotal + Indices::numEq);
}


template<class FluidSystem, class Indices>
void MultisegmentWellAssemble<FluidSystem,Indices>::
assemblePressureEq(const int seg,
                   const int seg_upwind,
                   const int outlet_segment_index,
                   const EvalWell& pressure_equation,
                   const EvalWell& outlet_pressure,
                   Equations& eqns1) const
{
    /*
        This method does *not* need communication.
    */
    MultisegmentWellEquationAccess<FluidSystem, Indices> eqns(eqns1);
    eqns.residual()[seg][SPres] += pressure_equation.value();
    eqns.D()[seg][seg][SPres][SPres] += pressure_equation.derivative(SPres + Indices::numEq);
    eqns.D()[seg][seg][SPres][WQTotal] += pressure_equation.derivative(WQTotal + Indices::numEq);
    if constexpr (has_wfrac_variable) {
        eqns.D()[seg][seg_upwind][SPres][WFrac] += pressure_equation.derivative(WFrac + Indices::numEq);
    }
    if constexpr (has_gfrac_variable) {
        eqns.D()[seg][seg_upwind][SPres][GFrac] += pressure_equation.derivative(GFrac + Indices::numEq);
    }

    // contribution from the outlet segment
    eqns.residual()[seg][SPres] -= outlet_pressure.value();
    for (int pv_idx = 0; pv_idx < numWellEq; ++pv_idx) {
        eqns.D()[seg][outlet_segment_index][SPres][pv_idx] -= outlet_pressure.derivative(pv_idx + Indices::numEq);
    }
}

template<class FluidSystem, class Indices>
void MultisegmentWellAssemble<FluidSystem,Indices>::
assembleTrivialEq(const int seg,
                  const Scalar value,
                  Equations& eqns1) const
{
    /*
        This method is called from MultisegmentWellEval::assembleICDPressureEq,
        which is called from MultisegmentWellEval::assemblePressureEq.
        This is the counterpart to assembleControlEquation, where assembleControlEquation is responsible for the top segment
        and assembleICDPressureEq is responsible for the remaining segments.
        This method does *not* need communication.
    */
    MultisegmentWellEquationAccess<FluidSystem, Indices> eqns(eqns1);
    eqns.residual()[seg][SPres] = value;
    eqns.D()[seg][seg][SPres][WQTotal] = 1.;
}

template<class FluidSystem, class Indices>
void MultisegmentWellAssemble<FluidSystem,Indices>::
assembleAccumulationTerm(const int seg,
                         const int comp_idx,
                         const EvalWell& accumulation_term,
                         Equations& eqns1) const
{
    /*
        This method is called from MultisegmentWell::assembleWellEqWithoutIteration.
        It only assembles on the diagonal of D and it does *not* need communication.
    */
    MultisegmentWellEquationAccess<FluidSystem, Indices> eqns(eqns1);
    eqns.residual()[seg][comp_idx] += accumulation_term.value();
    for (int pv_idx = 0; pv_idx < numWellEq; ++pv_idx) {
      eqns.D()[seg][seg][comp_idx][pv_idx] += accumulation_term.derivative(pv_idx + Indices::numEq);
    }
}

template<class FluidSystem, class Indices>
void MultisegmentWellAssemble<FluidSystem,Indices>::
assembleOutflowTerm(const int seg,
                    const int seg_upwind,
                    const int comp_idx,
                    const EvalWell& segment_rate,
                    Equations& eqns1) const
{
    /*
        This method is called from MultisegmentWell::assembleWellEqWithoutIteration.
        It does *not* need communication.
    */
    MultisegmentWellEquationAccess<FluidSystem, Indices> eqns(eqns1);
    eqns.residual()[seg][comp_idx] -= segment_rate.value();
    eqns.D()[seg][seg][comp_idx][WQTotal] -= segment_rate.derivative(WQTotal + Indices::numEq);
    if constexpr (has_wfrac_variable) {
        eqns.D()[seg][seg_upwind][comp_idx][WFrac] -= segment_rate.derivative(WFrac + Indices::numEq);
    }
    if constexpr (has_gfrac_variable) {
        eqns.D()[seg][seg_upwind][comp_idx][GFrac] -= segment_rate.derivative(GFrac + Indices::numEq);
    }
    // pressure derivative should be zero
}

template<class FluidSystem, class Indices>
void MultisegmentWellAssemble<FluidSystem,Indices>::
assembleInflowTerm(const int seg,
                   const int inlet,
                   const int inlet_upwind,
                   const int comp_idx,
                   const EvalWell& inlet_rate,
                   Equations& eqns1) const
 {
    /*
        This method is called from MultisegmentWell::assembleWellEqWithoutIteration.
        It does *not* need communication.
    */
    MultisegmentWellEquationAccess<FluidSystem, Indices> eqns(eqns1);
    eqns.residual()[seg][comp_idx] += inlet_rate.value();
    eqns.D()[seg][inlet][comp_idx][WQTotal] += inlet_rate.derivative(WQTotal + Indices::numEq);
    if constexpr (has_wfrac_variable) {
        eqns.D()[seg][inlet_upwind][comp_idx][WFrac] += inlet_rate.derivative(WFrac + Indices::numEq);
    }
    if constexpr (has_gfrac_variable) {
        eqns.D()[seg][inlet_upwind][comp_idx][GFrac] += inlet_rate.derivative(GFrac + Indices::numEq);
    }
    // pressure derivative should be zero
}

template<class FluidSystem, class Indices>
void MultisegmentWellAssemble<FluidSystem,Indices>::
assemblePerforationEq(const int seg,
                      const int local_perf_index,
                      const int comp_idx,
                      const EvalWell& cq_s_effective,
                      Equations& eqns1) const
{
    /*
        This method is called from MultisegmentWell::assembleWellEqWithoutIteration.
        It *does* need communication, i.e. this method only assembles the parts of the matrix this process is responsible for
        and after calling this function, the diagonal of the matrix D and the residual need to be combined by calling
        the function MultisegmentWellEquations::sumDistributed.
    */
    MultisegmentWellEquationAccess<FluidSystem, Indices> eqns(eqns1);
    // subtract sum of phase fluxes in the well equations.
    eqns.residual()[seg][comp_idx] += cq_s_effective.value();

    // assemble the jacobians
    for (int pv_idx = 0; pv_idx < numWellEq; ++pv_idx) {
        // also need to consider the efficiency factor when manipulating the jacobians.
        eqns.C()[seg][local_perf_index][pv_idx][comp_idx] -= cq_s_effective.derivative(pv_idx + Indices::numEq); // input in transformed matrix

        // the index name for the D should be eq_idx / pv_idx
        eqns.D()[seg][seg][comp_idx][pv_idx] += cq_s_effective.derivative(pv_idx + Indices::numEq);
    }

    for (int pv_idx = 0; pv_idx < Indices::numEq; ++pv_idx) {
        // also need to consider the efficiency factor when manipulating the jacobians.
        eqns.B()[seg][local_perf_index][comp_idx][pv_idx] += cq_s_effective.derivative(pv_idx);
    }
}


#include <opm/simulators/utils/InstantiationIndicesMacros.hpp>

INSTANTIATE_TYPE_INDICES(MultisegmentWellAssemble, double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE_INDICES(MultisegmentWellAssemble, float)
#endif

}
