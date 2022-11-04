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

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>

#include <opm/material/densead/DynamicEvaluation.hpp>

#include <opm/simulators/wells/StandardWellEquations.hpp>
#include <opm/simulators/wells/WellHelpers.hpp>

#include <functional>
#include <optional>
#include <vector>

namespace Opm
{

class DeferredLogger;
class GroupState;
class Schedule;
class SummaryState;
class WellContributions;
template<class T> class WellInterfaceFluidSystem;
class WellState;

template<class FluidSystem, class Indices, class Scalar>
class StandardWellAssemble {
public:
    StandardWellAssemble(const WellInterfaceFluidSystem<FluidSystem>& well)
        : well_(well)
    {}

    template<class EvalWell>
    void assembleControlEq(const WellState& well_state,
                           const GroupState& group_state,
                           const Schedule& schedule,
                           const SummaryState& summaryState,
                           const int numWellEq,
                           const EvalWell& wqTotal,
                           const EvalWell& bhp,
                           const std::function<EvalWell(int)>& getQs,
                           const double rho,
                           StandardWellEquations<Indices,Scalar>& A,
                           DeferredLogger& deferred_logger) const;

    template<class EvalWell>
    void assembleInjectivityEq(const EvalWell& eq_pskin,
                               const EvalWell& eq_wat_vel,
                               const int pskin_index,
                               const int wat_vel_index,
                               const int cell_idx,
                               const int numWellEq,
                               StandardWellEquations<Indices,Scalar>& eqns) const;

    template<class EvalWell>
    void assemblePerforationEq(const EvalWell& cq_s_effective,
                               const int componentIdx,
                               const int cell_idx,
                               const int numWellEq,
                               StandardWellEquations<Indices,Scalar>& eqns) const;

    template<class EvalWell>
    void assembleZFracEq(const EvalWell& cq_s_zfrac_effective,
                         const int cell_idx,
                         const int numWellEq,
                         StandardWellEquations<Indices,Scalar>& eqns) const;

    template<class EvalWell>
    void assembleSourceEq(const EvalWell& resWell_loc,
                          const int componentIdx,
                          const int numWellEq,
                          StandardWellEquations<Indices,Scalar>& eqns) const;

    template<class PressureMatrix, class BVector>
    void addWellPressureEqs(const BVector& x,
                            const int pressureVarIndex,
                            const bool use_well_weights,
                            const WellState& well_state,
                            const int numWellEq,
                            const StandardWellEquations<Indices,Scalar>& eqns,
                            PressureMatrix& mat) const;

private:
    const WellInterfaceFluidSystem<FluidSystem>& well_;
};

}

#endif // OPM_STANDARDWELL_ASSEMBLE_HEADER_INCLUDED
