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
template<class Scalar> class GroupState;
template<class Scalar, int numWellEq, int numEq> class MultisegmentWellEquations;
template<class FluidSystem, class Indices> class MultisegmentWellPrimaryVariables;
class Schedule;
class SummaryState;
template<class FluidSystem, class Indices> class WellInterfaceIndices;
template<class Scalar> class WellState;

//! \brief Class handling assemble of the equation system for MultisegmentWell.
template<class FluidSystem, class Indices>
class MultisegmentWellAssemble
{
    using PrimaryVariables = MultisegmentWellPrimaryVariables<FluidSystem,Indices>;
    static constexpr int WQTotal = PrimaryVariables::WQTotal;
    static constexpr bool has_wfrac_variable = PrimaryVariables::has_wfrac_variable;
    static constexpr bool has_gfrac_variable = PrimaryVariables::has_gfrac_variable;
    static constexpr int WFrac = PrimaryVariables::WFrac;
    static constexpr int GFrac = PrimaryVariables::GFrac;
    static constexpr int SPres = PrimaryVariables::SPres;

public:
    static constexpr int numWellEq = Indices::numPhases+1;
    using Scalar = typename FluidSystem::Scalar;
    using Equations = MultisegmentWellEquations<Scalar,numWellEq,Indices::numEq>;
    using EvalWell = DenseAd::Evaluation<Scalar, numWellEq+Indices::numEq>;

    //! \brief Constructor initializes reference to well.
    explicit MultisegmentWellAssemble(const WellInterfaceIndices<FluidSystem,Indices>& well)
        : well_(well)
    {}

    //! \brief Assemble control equation.
    void assembleControlEq(const WellState<Scalar>& well_state,
                           const GroupState<Scalar>& group_state,
                           const Schedule& schedule,
                           const SummaryState& summaryState,
                           const Well::InjectionControls& inj_controls,
                           const Well::ProductionControls& prod_controls,
                           const Scalar rho,
                           const PrimaryVariables& primary_variables,
                           Equations& eqns,
                           const bool stopped_or_zero_target,
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
                            Equations& eqns) const;

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
                               const int local_perf_index,
                               const int comp_idx,
                               const EvalWell& cq_s_effective,
                               Equations& eqns) const;

private:
    const WellInterfaceIndices<FluidSystem,Indices>& well_; //!< Reference to well
};

}

#endif // OPM_STANDARDWELL_ASSEMBLE_HEADER_INCLUDED
