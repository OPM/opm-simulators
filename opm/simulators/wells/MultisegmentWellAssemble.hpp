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

#include <functional>

namespace Opm
{

class DeferredLogger;
class GroupState;
template<class Scalar, int numWellEq, int numEq> class MultisegmentWellEquations;
class Schedule;
class SummaryState;
template<class FluidSystem, class Indices, class Scalar> class WellInterfaceIndices;
class WellState;

//! \brief Class handling assemble of the equation system for MultisegmentWell.
template<class FluidSystem, class Indices, class Scalar>
class MultisegmentWellAssemble
{
public:
    static constexpr int numWellEq = Indices::numPhases+1;
    using Equations = MultisegmentWellEquations<Scalar,numWellEq,Indices::numEq>;
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
                           const EvalWell& wqTotal,
                           const EvalWell& bhp,
                           const int SPres,
                           const std::function<EvalWell(const int)>& getQs,
                           Equations& eqns,
                           DeferredLogger& deferred_logger) const;

private:
    const WellInterfaceIndices<FluidSystem,Indices,Scalar>& well_; //!< Reference to well
};

}

#endif // OPM_STANDARDWELL_ASSEMBLE_HEADER_INCLUDED
