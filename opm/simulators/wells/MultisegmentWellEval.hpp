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

#include <opm/material/densead/Evaluation.hpp>

#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <array>
#include <memory>
#include <utility>
#include <vector>

namespace Dune {
template<class Matrix> class UMFPack;
}

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
    // TODO: for now, not considering the polymer, solvent and so on to simplify the development process.

    // TODO: we need to have order for the primary variables and also the order for the well equations.
    // sometimes, they are similar, while sometimes, they can have very different forms.

    // Table showing the primary variable indices, depending on what phases are present:
    //
    //         WOG     OG     WG     WO    W/O/G (single phase)
    // WQTotal   0      0      0      0                       0
    // WFrac     1  -1000      1      1                   -1000
    // GFrac     2      1  -1000  -1000                   -1000
    // Spres     3      2      2      2                       1

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

    //  the number of well equations  TODO: it should have a more general strategy for it
    static constexpr int numWellEq = Indices::numPhases + 1;

    using Equations = MultisegmentWellEquations<Scalar,numWellEq,Indices::numEq>;

    using BVector = typename Equations::BVector;
    using BVectorWell = typename Equations::BVectorWell;

    // TODO: for more efficient implementation, we should have EvalReservoir, EvalWell, and EvalRerservoirAndWell
    //                                                         EvalR (Eval), EvalW, EvalRW
    // TODO: for now, we only use one type to save some implementation efforts, while improve later.
    using EvalWell = DenseAd::Evaluation<double, /*size=*/Indices::numEq + numWellEq>;
    using Eval = DenseAd::Evaluation<Scalar, /*size=*/Indices::numEq>;

public:
    //! \brief Returns a const reference to equation system.
    const Equations& linSys() const
    { return linSys_; }

protected:
    MultisegmentWellEval(WellInterfaceIndices<FluidSystem,Indices,Scalar>& baseif);

    void initMatrixAndVectors(const int num_cells);
    void initPrimaryVariablesEvaluation() const;

    void assembleDefaultPressureEq(const int seg,
                                   WellState& well_state);


    // assemble pressure equation for ICD segments
    void assembleICDPressureEq(const int seg,
                               const UnitSystem& unit_system,
                               WellState& well_state,
                               DeferredLogger& deferred_logger);


    void assemblePressureEq(const int seg,
                            const UnitSystem& unit_system,
                            WellState& well_state,
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
                                         const bool relax_tolerance) const;

    // handling the overshooting and undershooting of the fractions
    void processFractions(const int seg) const;

    void updatePrimaryVariables(const WellState& well_state) const;

    void updateUpwindingSegments();

    // updating the well_state based on well solution dwells
    void updatePrimaryVariablesNewton(const BVectorWell& dwells,
                                      const double relaxation_factor,
                                      const double DFLimit,
                                      const double max_pressure_change) const;

    void computeSegmentFluidProperties(const EvalWell& temperature,
                                       const EvalWell& saltConcentration,
                                       int pvt_region_index,
                                       DeferredLogger& deferred_logger);

    EvalWell getBhp() const;
    EvalWell getFrictionPressureLoss(const int seg) const;
    EvalWell getHydroPressureLoss(const int seg) const;
    EvalWell getQs(const int comp_idx) const;
    EvalWell getSegmentWQTotal(const int seg) const;
    EvalWell getSegmentPressure(const int seg) const;
    EvalWell getSegmentRate(const int seg, const int comp_idx) const;
    EvalWell getSegmentRateUpwinding(const int seg,
                                     const size_t comp_idx) const;
    EvalWell getSegmentSurfaceVolume(const EvalWell& temperature,
                                     const EvalWell& saltConcentration,
                                     const int pvt_region_index,
                                     const int seg_idx) const;
    EvalWell getWQTotal() const;


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

    // pressure drop for Autonomous ICD segment (WSEGAICD)
    EvalWell pressureDropAutoICD(const int seg,
                                 const UnitSystem& unit_system) const;

    // pressure drop for Spiral ICD segment (WSEGSICD)
    EvalWell pressureDropSpiralICD(const int seg) const;

    // pressure drop for sub-critical valve (WSEGVALV)
    EvalWell pressureDropValve(const int seg) const;

    void updateWellStateFromPrimaryVariables(WellState& well_state,
                                             const double rho,
                                             DeferredLogger& deferred_logger) const;

    // fraction value of the primary variables
    // should we just use member variables to store them instead of calculating them again and again
    EvalWell volumeFraction(const int seg,
                            const unsigned compIdx) const;

    // F_p / g_p, the basic usage of this value is because Q_p = G_t * F_p / G_p
    EvalWell volumeFractionScaled(const int seg,
                                  const int comp_idx) const;

    // basically Q_p / \sigma_p Q_p
    EvalWell surfaceVolumeFraction(const int seg,
                                   const int comp_idx) const;

    // convert a Eval from reservoir to contain the derivative related to wells
    EvalWell extendEval(const Eval& in) const;

    const WellInterfaceIndices<FluidSystem,Indices,Scalar>& baseif_;

    Equations linSys_; //!< The equation system

    // the values for the primary varibles
    // based on different solutioin strategies, the wells can have different primary variables
    mutable std::vector<std::array<double, numWellEq> > primary_variables_;

    // the Evaluation for the well primary variables, which contain derivativles and are used in AD calculation
    mutable std::vector<std::array<EvalWell, numWellEq> > primary_variables_evaluation_;

    // the upwinding segment for each segment based on the flow direction
    std::vector<int> upwinding_segments_;

    // the densities of segment fluids
    // we should not have this member variable
    std::vector<EvalWell> segment_densities_;

    // the mass rate of the segments
    std::vector<EvalWell> segment_mass_rates_;

    // the viscosity of the segments
    std::vector<EvalWell> segment_viscosities_;

    std::vector<std::vector<EvalWell>> segment_phase_densities_;
    std::vector<std::vector<EvalWell>> segment_phase_fractions_;
    std::vector<std::vector<EvalWell>> segment_phase_viscosities_;

    // depth difference between perforations and the perforated grid cells
    std::vector<double> cell_perforation_depth_diffs_;
    // pressure correction due to the different depth of the perforation and
    // center depth of the grid block
    std::vector<double> cell_perforation_pressure_diffs_;
};

}

#endif // OPM_MULTISEGMENTWELL_GENERIC_HEADER_INCLUDED
