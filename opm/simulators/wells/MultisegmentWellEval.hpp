/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.

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


#ifndef OPM_MULTISEGMENTWELL_EVAL_HEADER_INCLUDED
#define OPM_MULTISEGMENTWELL_EVAL_HEADER_INCLUDED

#include <opm/simulators/wells/MultisegmentWellEquations.hpp>
#include <opm/simulators/wells/MultisegmentWellGeneric.hpp>
#include <opm/simulators/wells/MultisegmentWellPrimaryVariables.hpp>
#include <opm/simulators/wells/MultisegmentWellSegments.hpp>

#include <opm/material/densead/Evaluation.hpp>

#include <array>
#include <memory>
#include <utility>
#include <vector>

namespace Opm
{

class ConvergenceReport;
class GroupState;
class Schedule;
class WellContributions;
template<class FluidSystem, class Indices, class Scalar> class WellInterfaceIndices;
class WellState;

template<typename FluidSystem, typename Indices, typename Scalar>
class MultisegmentWellEval : public MultisegmentWellGeneric<Scalar>
{
protected:
    using PrimaryVariables = MultisegmentWellPrimaryVariables<FluidSystem,Indices,Scalar>;
    static constexpr int numWellEq = PrimaryVariables::numWellEq;
    static constexpr int SPres = PrimaryVariables::SPres;
    static constexpr int WQTotal = PrimaryVariables::WQTotal;

    using Equations = MultisegmentWellEquations<Scalar,numWellEq,Indices::numEq>;
    using MSWSegments = MultisegmentWellSegments<FluidSystem,Indices,Scalar>;

    using BVector = typename Equations::BVector;
    using BVectorWell = typename Equations::BVectorWell;

    // TODO: for more efficient implementation, we should have EvalReservoir, EvalWell, and EvalRerservoirAndWell
    //                                                         EvalR (Eval), EvalW, EvalRW
    // TODO: for now, we only use one type to save some implementation efforts, while improve later.
    using EvalWell = typename PrimaryVariables::EvalWell;
    using Eval = DenseAd::Evaluation<Scalar, /*size=*/Indices::numEq>;

public:
    //! \brief Returns a const reference to equation system.
    const Equations& linSys() const
    { return linSys_; }

protected:
    MultisegmentWellEval(WellInterfaceIndices<FluidSystem,Indices,Scalar>& baseif);

    void initMatrixAndVectors(const int num_cells);

    void assembleDefaultPressureEq(const int seg,
                                   WellState& well_state,
                                   const bool use_average_density);

    // assemble pressure equation for ICD segments
    void assembleICDPressureEq(const int seg,
                               const UnitSystem& unit_system,
                               WellState& well_state,
                               const bool use_average_density,
                               DeferredLogger& deferred_logger);

    void assembleAccelerationAndHydroPressureLosses(const int seg,
                                                    WellState& well_state,
                                                    const bool use_average_density);


    void assemblePressureEq(const int seg,
                            const UnitSystem& unit_system,
                            WellState& well_state,
                            const bool use_average_density,
                            DeferredLogger& deferred_logger);

    /// check whether the well equations get converged for this well
    ConvergenceReport getWellConvergence(const WellState& well_state,
                                         const std::vector<double>& B_avg,
                                         DeferredLogger& deferred_logger,
                                         const double max_residual_allowed,
                                         const double tolerance_wells,
                                         const double relaxed_inner_tolerance_flow_ms_well,
                                         const double tolerance_pressure_ms_wells,
                                         const double relaxed_inner_tolerance_pressure_ms_well,
                                         const bool relax_tolerance, 
                                         const bool well_is_stopped) const;

    std::pair<bool, std::vector<Scalar> >
    getFiniteWellResiduals(const std::vector<Scalar>& B_avg,
                           DeferredLogger& deferred_logger) const;

    double getControlTolerance(const WellState& well_state,
                               const double tolerance_wells,
                               const double tolerance_pressure_ms_wells,
                               DeferredLogger& deferred_logger) const;

    double getResidualMeasureValue(const WellState& well_state,
                                   const std::vector<double>& residuals,
                                   const double tolerance_wells,
                                   const double tolerance_pressure_ms_wells,
                                   DeferredLogger& deferred_logger) const;

    void handleAccelerationPressureLoss(const int seg,
                                        WellState& well_state);

    EvalWell pressureDropAutoICD(const int seg,
                                 const UnitSystem& unit_system) const;

    // convert a Eval from reservoir to contain the derivative related to wells
    EvalWell extendEval(const Eval& in) const;

    const WellInterfaceIndices<FluidSystem,Indices,Scalar>& baseif_;

    Equations linSys_; //!< The equation system
    PrimaryVariables primary_variables_; //!< The primary variables
    MSWSegments segments_; //!< Segment properties

    // depth difference between perforations and the perforated grid cells
    std::vector<double> cell_perforation_depth_diffs_;
    // pressure correction due to the different depth of the perforation and
    // center depth of the grid block
    std::vector<double> cell_perforation_pressure_diffs_;
};

}

#endif // OPM_MULTISEGMENTWELL_GENERIC_HEADER_INCLUDED
