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


#ifndef OPM_MULTISEGMENTWELL_ASSEMBLE_HEADER_INCLUDED
#define OPM_MULTISEGMENTWELL_ASSEMBLE_HEADER_INCLUDED

#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/material/densead/Evaluation.hpp>

namespace Opm
{

class DeferredLogger;
class GroupState;
template<class Scalar, int numWellEq, int numEq> class MultisegmentWellEquations;
template<class FluidSystem, class Indices, class Scalar> class MultisegmentWellPrimaryVariables;
class Schedule;
class SummaryState;
template<class FluidSystem, class Indices, class Scalar> class WellInterfaceIndices;
class WellState;

//! \brief Class handling assemble of the equation system for MultisegmentWell.
template<class FluidSystem, class Indices, class Scalar>
class MultisegmentWellAssemble
{
    static constexpr bool has_water = (Indices::waterSwitchIdx >= 0);
    static constexpr bool has_gas = (Indices::compositionSwitchIdx >= 0);
    static constexpr bool has_oil = (Indices::numPhases - has_gas - has_water) > 0;

    // In the implementation, one should use has_wfrac_variable
    // rather than has_water to check if you should do something
    // with the variable at the WFrac location, similar for GFrac.
    static constexpr bool has_wfrac_variable = has_water && Indices::numPhases > 1;
    static constexpr bool has_gfrac_variable = has_gas && has_oil;

    static constexpr int WQTotal = 0;
    static constexpr int WFrac = has_wfrac_variable ? 1 : -1000;
    static constexpr int GFrac = has_gfrac_variable ? has_wfrac_variable + 1 : -1000;
    static constexpr int SPres = has_wfrac_variable + has_gfrac_variable + 1;

public:
    static constexpr int numWellEq = Indices::numPhases+1;
    using Equations = MultisegmentWellEquations<Scalar,numWellEq,Indices::numEq>;
    using PrimaryVariables = MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>;
    using EvalWell = DenseAd::Evaluation<Scalar, numWellEq+Indices::numEq>;

    //! \brief Constructor initializes reference to well.
    MultisegmentWellAssemble(const WellInterfaceIndices<FluidSystem,Indices,Scalar>& well)
        : well_(well)
    {}

    //! \brief Assemble control equation.
    void assembleControlEq(const WellState& well_state,
                           const GroupState& group_state,
                           const Schedule& schedule,
                           const SummaryState& summaryState,
                           const Well::InjectionControls& inj_controls,
                           const Well::ProductionControls& prod_controls,
                           const double rho,
                           const PrimaryVariables& primary_variables,
                           Equations& eqns,
                           DeferredLogger& deferred_logger) const;

    //! \brief Assemble piece of the acceleration term
    void assembleAccelerationTerm(const int seg_target,
                                  const int seg,
                                  const int seg_upwing,
                                  const EvalWell& accelerationTerm,
                                  Equations& eqns1) const;

    //! \brief Assemble hydraulic pressure term
    void assembleHydroPressureLoss(const int seg,
                                   const int seg_density,
                                   const EvalWell& hydro_pressure_drop_seg,
                                   Equations& eqns1) const;

    //! \brief Assemble additional derivatives due to reverse flow
    void assemblePressureEqExtraDerivatives(const int seg,
                                            const int seg_upwind,
                                            const EvalWell& extra_derivatives,
                                            Equations& eqns1) const;

    //! \brief Assemble pressure terms.
    void assemblePressureEq(const int seg,
                            const int seg_upwind,
                            const int outlet_segment_index,
                            const EvalWell& pressure_equation,
                            const EvalWell& outlet_pressure,
                            Equations& eqns,
                            bool wfrac = has_wfrac_variable,
                            bool gfrac = has_gfrac_variable) const;

    //! \brief Assembles a trivial equation.
    void assembleTrivialEq(const int seg,
                           const Scalar value,
                           Equations& eqns) const;

    //! \brief Assemble accumulation term.
    void assembleAccumulationTerm(const int seg,
                                  const int comp_idx,
                                  const EvalWell& accumulation_term,
                                  Equations& eqns1) const;

    //! \brief Assemble outflow term.
    void assembleOutflowTerm(const int seg,
                        const int seg_upwind,
                        const int comp_idx,
                        const EvalWell& segment_rate,
                        Equations& eqns1) const;

    //! \brief Assemble inflow term.
    void assembleInflowTerm(const int seg,
                            const int inlet,
                            const int inlet_upwind,
                            const int comp_idx,
                            const EvalWell& inlet_rate,
                            Equations& eqns) const;

    //! \brief Assemble equation for a perforation.
    void assemblePerforationEq(const int seg,
                               const int cell_idx,
                               const int comp_idx,
                               const EvalWell& cq_s_effective,
                               Equations& eqns) const;

private:
    const WellInterfaceIndices<FluidSystem,Indices,Scalar>& well_; //!< Reference to well
};

}

#endif // OPM_STANDARDWELL_ASSEMBLE_HEADER_INCLUDED
