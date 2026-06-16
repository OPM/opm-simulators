/*
  Copyright 2026 SINTEF Digital, Mathematics and Cybernetics.

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

#ifndef OPM_GRAPH_MULTISEGMENTWELL_HEADER_INCLUDED
#define OPM_GRAPH_MULTISEGMENTWELL_HEADER_INCLUDED

#include <opm/simulators/wells/MultisegmentWell.hpp>

#include <flowexperimental/graphwell/GraphWellAssembler.hpp>
#include <flowexperimental/graphwell/GraphWellEquations.hpp>
#include <flowexperimental/graphwell/GraphWellFluidProperties.hpp>
#include <flowexperimental/graphwell/GraphWellPrimaryVariables.hpp>
#include <flowexperimental/graphwell/GraphWellTopology.hpp>

#include <array>
#include <utility>
#include <vector>

namespace Opm {

/// \brief Multisegment well built on the GraphWell entity-split formulation with a
/// native AD face-based assembler.
///
/// This derives from the production MultisegmentWell and reuses its control / THP /
/// group machinery (the well control equation is borrowed from the base assembly).
/// Everything else -- the mass-conservation equations on segments, the momentum
/// (pressure-drop) equations on connections, and the perforation inflow -- is built
/// from scratch by GraphWellAssembler using fixed-size AD with upwind weighting, into
/// a Dune::MultiTypeBlockMatrix. The whole well (inner control iterations and the
/// global reservoir coupling) is then driven by that system.
///
/// The GraphWell unknowns map onto the MultisegmentWell unknowns by the tree
/// bijection (segment pressure/fractions shared; one flux per connection with
/// Q_c = -WQTotal_{down(c)}), which is used to read the linearisation point from the
/// base primary variables, to borrow the control row, and to push Newton increments
/// back through the base well-state update.
template<typename TypeTag>
class GraphMultisegmentWell : public MultisegmentWell<TypeTag>
{
public:
    using Base = MultisegmentWell<TypeTag>;
    using MSWEval = typename Base::MSWEval;

    using typename Base::Simulator;
    using typename Base::FluidSystem;
    using typename Base::Indices;
    using typename Base::Scalar;
    using typename Base::Eval;
    using typename Base::ModelParameters;
    using typename Base::RateConverterType;
    using typename Base::BVector;
    using typename Base::WellStateType;
    using typename Base::GroupStateHelperType;
    using typename Base::PressureMatrix;

    static constexpr int NP = Indices::numPhases;
    static constexpr int NumResEq = Indices::numEq;

    using GraphTopo = GraphWellTopology<Scalar>;
    using GraphEqns = GraphWellEquations<Scalar, NP, NumResEq>;

    // Compact AD type for the control / IPR equations: derivative slots are the top
    // segment DOFs [0,NP) plus the surface-connection flux at slot NP. Padded to at
    // least 3 so the generic WellAssemble / WellBhpThpCalculator helpers (instantiated
    // from size 3) cover single-phase wells too.
    static constexpr int QSlot = NP;
    static constexpr int CEvalSize = (NP + 1 < 3) ? 3 : NP + 1;
    using CEval = DenseAd::Evaluation<Scalar, CEvalSize>;
    using GPV = GraphWellPrimaryVariables<FluidSystem, NP, NumResEq>;
    using GFP = GraphWellFluidProperties<FluidSystem, NP, NumResEq>;
    using GAsm = GraphWellAssembler<FluidSystem, NP, NumResEq>;
    using GEval = typename GPV::Eval;

    GraphMultisegmentWell(const Well& well,
                          const ParallelWellInfo<Scalar>& pw_info,
                          const int time_step,
                          const ModelParameters& param,
                          const RateConverterType& rate_converter,
                          const int pvtRegionIdx,
                          const int num_conservation_quantities,
                          const int num_phases,
                          const int index_of_well,
                          const std::vector<PerforationData<Scalar>>& perf_data);

    void init(const std::vector<Scalar>& depth_arg,
              const Scalar gravity_arg,
              const std::vector<Scalar>& B_avg,
              const bool changed_to_open_this_step) override;

    void assembleWellEqWithoutIteration(const Simulator& simulator,
                                        const GroupStateHelperType& groupStateHelper,
                                        const double dt,
                                        const Well::InjectionControls& inj_controls,
                                        const Well::ProductionControls& prod_controls,
                                        WellStateType& well_state,
                                        const bool solving_with_zero_rate) override;

    void solveEqAndUpdateWellState(const Simulator& simulator,
                                   const GroupStateHelperType& groupStateHelper,
                                   WellStateType& well_state) override;

    ConvergenceReport getWellConvergence(const GroupStateHelperType& groupStateHelper,
                                         const std::vector<Scalar>& B_avg,
                                         const bool relax_tolerance) const override;

    void apply(const BVector& x, BVector& Ax) const override;
    void apply(BVector& r) const override;

    //! CPRW well-pressure coupling is not supported by the GraphWell formulation; the
    //! well is handled through the Schur apply() path instead (run with a non-CPRW
    //! linear solver, or accept a CPR preconditioner without well-pressure coupling).
    void addWellPressureEquations(PressureMatrix& mat,
                                  const BVector& x,
                                  const int pressureVarIndex,
                                  const bool use_well_weights,
                                  const WellStateType& well_state) const override;

    void recoverWellSolutionAndUpdateWellState(const Simulator& simulator,
                                               const BVector& x,
                                               const GroupStateHelperType& groupStateHelper,
                                               WellStateType& well_state) override;

    //! The production updateIPRImplicit solves the base MultisegmentWell linear system,
    //! which the GraphWell never assembles. Use the explicit IPR fallback instead.
    void updateIPRImplicit(const Simulator& simulator,
                           const GroupStateHelperType& groupStateHelper,
                           WellStateType& well_state) override;

protected:
    bool iterateWellEqWithControl(const Simulator& simulator,
                                  const double dt,
                                  const Well::InjectionControls& inj_controls,
                                  const Well::ProductionControls& prod_controls,
                                  const GroupStateHelperType& groupStateHelper,
                                  WellStateType& well_state) override;

    bool iterateWellEqWithSwitching(const Simulator& simulator,
                                    const double dt,
                                    const Well::InjectionControls& inj_controls,
                                    const Well::ProductionControls& prod_controls,
                                    const GroupStateHelperType& groupStateHelper,
                                    WellStateType& well_state,
                                    const bool fixed_control,
                                    const bool fixed_status,
                                    const bool solving_with_zero_rate) override;

private:
    //! Per-equation finite-residual check on the GraphWell system (native replacement
    //! for MSWEval::getFiniteWellResiduals used by the inner Newton loops).
    std::pair<bool, std::vector<Scalar>>
    getFiniteWellResidualsGraph(const std::vector<Scalar>& B_avg,
                                DeferredLogger& deferred_logger) const;
    //! Copy the current linearisation point from the base primary variables.
    void setGraphStateFromBase();
    //! Gather per-perforation reservoir quantities (with reservoir-AD derivatives).
    void gatherPerforations(const Simulator& simulator,
                            const GroupStateHelperType& groupStateHelper,
                            std::vector<typename GAsm::PerforationInput>& perfs) const;
    //! Assemble the well control equation natively into the surface connection,
    //! reusing the generic (StandardWell-shared) WellAssemble / THP helpers.
    void assembleControl(const Simulator& simulator,
                         const GroupStateHelperType& groupStateHelper,
                         WellStateType& well_state,
                         const Well::InjectionControls& inj_controls,
                         const Well::ProductionControls& prod_controls,
                         bool solving_with_zero_rate);
    //! Set the reservoir source terms (connectionRates_) and well-state perforation
    //! rates from the native perforation rates (replaces the base bookkeeping).
    void updateWellStateRates(const GFP& fp,
                              const std::vector<typename GAsm::PerforationInput>& perfs,
                              bool allow_cf,
                              WellStateType& well_state);

    //! Translate a GraphWell increment back to a MultisegmentWell increment and apply it.
    void applyGraphSolution(const Simulator& simulator,
                            const typename GraphEqns::BVectorWell& gxw,
                            const GroupStateHelperType& groupStateHelper,
                            WellStateType& well_state,
                            Scalar relaxation_factor = 1.0);

    //! Top-segment bhp as a compact control AD value (derivative in slot 0).
    CEval topBhpCEval() const;
    //! Top surface phase rate getQs(comp) = -Q_surface * volumeFractionScaled(top,comp)
    //! as a compact control AD value (derivatives in the top-segment + flux slots).
    CEval surfaceRateCEval(int comp) const;

    //! Embed a reservoir Eval into the GraphWell Eval (reservoir derivative slots).
    GEval fromRes(const Eval& e) const
    {
        GEval g(e.value());
        for (int k = 0; k < NumResEq; ++k)
            g.setDerivative(k, e.derivative(k));
        return g;
    }

    GraphTopo topo_;
    GraphEqns graph_eqns_;
    GPV gpv_;
    GAsm assembler_;

    // segment -> connection whose 'down' end is that segment (bijection on a tree)
    std::vector<int> outlet_conn_;
    // MultisegmentWell primary-variable index -> GraphWell segment DOF (-1 for WQTotal)
    std::array<int, MSWEval::numWellEq> var_to_gdof_;
    int msw_wqtotal_{0};

    // phase map for the native control equation
    int ref_comp_{0};
    std::array<int, NP> gdof_comp_{};   // active component stored at segment DOF g (g>=1)
    std::array<Scalar, NP> scale_comp_{};
};

} // namespace Opm

#include <flowexperimental/graphwell/GraphMultisegmentWell_impl.hpp>

#endif // OPM_GRAPH_MULTISEGMENTWELL_HEADER_INCLUDED
