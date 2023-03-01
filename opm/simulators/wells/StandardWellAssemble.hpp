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


#ifndef OPM_STANDARDWELL_ASSEMBLE_HEADER_INCLUDED
#define OPM_STANDARDWELL_ASSEMBLE_HEADER_INCLUDED

#include <opm/input/eclipse/Schedule/Well/Well.hpp>

namespace Opm
{

class DeferredLogger;
class GroupState;
class Schedule;
template<class Scalar, int numEq> class StandardWellEquations;
template<class FluidSystem, class Indices, class Scalar> class StandardWellPrimaryVariables;
class SummaryState;
template<class FluidSystem> class WellInterfaceFluidSystem;
class WellState;

//! \brief Class handling assemble of the equation system for StandardWell.
template<class FluidSystem, class Indices, class Scalar>
class StandardWellAssemble
{
public:
    using PrimaryVariables = StandardWellPrimaryVariables<FluidSystem,Indices,Scalar>;
    using EvalWell = typename PrimaryVariables::EvalWell;

    //! \brief Constructor initializes reference to well.
    StandardWellAssemble(const WellInterfaceFluidSystem<FluidSystem>& well)
        : well_(well)
    {}

    //! \brief Assemble control equation.
    void assembleControlEq(const WellState& well_state,
                           const GroupState& group_state,
                           const Schedule& schedule,
                           const SummaryState& summaryState,
                           const Well::InjectionControls& inj_controls,
                           const Well::ProductionControls& prod_controls,
                           const PrimaryVariables& primary_variables,
                           const double rho,
                           StandardWellEquations<Scalar,Indices::numEq>& eqns,
                           DeferredLogger& deferred_logger) const;

    //! \brief Assemble injectivity equation.
    void assembleInjectivityEq(const EvalWell& eq_pskin,
                               const EvalWell& eq_wat_vel,
                               const int pskin_index,
                               const int wat_vel_index,
                               const int cell_idx,
                               const int numWellEq,
                               StandardWellEquations<Scalar,Indices::numEq>& eqns) const;

    //! \brief Assemble equation for a perforation.
    void assemblePerforationEq(const EvalWell& cq_s_effective,
                               const int componentIdx,
                               const int cell_idx,
                               const int numWellEq,
                               StandardWellEquations<Scalar,Indices::numEq>& eqns) const;

    //! \brief Assemble equation for Z fraction.
    void assembleZFracEq(const EvalWell& cq_s_zfrac_effective,
                         const int cell_idx,
                         const int numWellEq,
                         StandardWellEquations<Scalar,Indices::numEq>& eqns) const;

    //! \brief Assemble a source term.
    void assembleSourceEq(const EvalWell& resWell_loc,
                          const int componentIdx,
                          const int numWellEq,
                          StandardWellEquations<Scalar,Indices::numEq>& eqns) const;

private:
    const WellInterfaceFluidSystem<FluidSystem>& well_; //!< Reference to well
};

}

#endif // OPM_STANDARDWELL_ASSEMBLE_HEADER_INCLUDED
