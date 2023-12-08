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


#ifndef OPM_STANDARDWELL_EVAL_HEADER_INCLUDED
#define OPM_STANDARDWELL_EVAL_HEADER_INCLUDED

#include <opm/simulators/wells/StandardWellConnections.hpp>
#include <opm/simulators/wells/StandardWellEquations.hpp>
#include <opm/simulators/wells/StandardWellPrimaryVariables.hpp>

#include <opm/material/densead/Evaluation.hpp>

#include <optional>
#include <vector>

namespace Opm
{

class ConvergenceReport;
class DeferredLogger;
class GroupState;
class Schedule;
class SummaryState;
class WellContributions;
template<class FluidSystem, class Indices, class Scalar> class WellInterfaceIndices;
class WellState;

template<class FluidSystem, class Indices, class Scalar>
class StandardWellEval
{
protected:
    using PrimaryVariables = StandardWellPrimaryVariables<FluidSystem,Indices,Scalar>;
    using StdWellConnections = StandardWellConnections<FluidSystem,Indices,Scalar>;
    static constexpr int Bhp = PrimaryVariables::Bhp;
    static constexpr int WQTotal= PrimaryVariables::WQTotal;
    static constexpr int numWellConservationEq = PrimaryVariables::numWellConservationEq;

    static constexpr bool has_wfrac_variable = PrimaryVariables::has_wfrac_variable;
    static constexpr bool has_gfrac_variable = PrimaryVariables::has_gfrac_variable;
    static constexpr int WFrac = PrimaryVariables::WFrac;
    static constexpr int GFrac = PrimaryVariables::GFrac;
    static constexpr int SFrac = PrimaryVariables::SFrac;

public:
    using EvalWell = typename PrimaryVariables::EvalWell;
    using Eval = DenseAd::Evaluation<Scalar, Indices::numEq>;
    using BVectorWell = typename StandardWellEquations<Scalar,Indices::numEq>::BVectorWell;

    //! \brief Returns a const reference to equation system.
    const StandardWellEquations<Scalar,Indices::numEq>& linSys() const
    { return linSys_; }

protected:
    StandardWellEval(const WellInterfaceIndices<FluidSystem,Indices,Scalar>& baseif);

    const WellInterfaceIndices<FluidSystem,Indices,Scalar>& baseif_;

    EvalWell extendEval(const Eval& in) const;

    // computing the accumulation term for later use in well mass equations
    void computeAccumWell();

    ConvergenceReport getWellConvergence(const WellState& well_state,
                                         const std::vector<double>& B_avg,
                                         const double maxResidualAllowed,
                                         const double tol_wells,
                                         const double relaxed_tolerance_flow,
                                         const bool relax_tolerance,
                                         const bool well_is_stopped, 
                                         std::vector<double>& res,
                                         DeferredLogger& deferred_logger) const;

    void init(std::vector<double>& perf_depth,
              const std::vector<double>& depth_arg,
              const int num_cells,
              const bool has_polymermw);

    void updateWellStateFromPrimaryVariables(const bool stop_or_zero_rate_target,
                                             WellState& well_state,
                                             const SummaryState& summary_state,
                                             DeferredLogger& deferred_logger) const;

    PrimaryVariables primary_variables_; //!< Primary variables for well

    // the saturations in the well bore under surface conditions at the beginning of the time step
    std::vector<double> F0_;

    StandardWellEquations<Scalar,Indices::numEq> linSys_; //!< Linear equation system
    StdWellConnections connections_; //!< Connection level values
};

}

#endif // OPM_STANDARDWELL_EVAL_HEADER_INCLUDED
