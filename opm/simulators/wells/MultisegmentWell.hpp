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


#ifndef OPM_MULTISEGMENTWELL_HEADER_INCLUDED
#define OPM_MULTISEGMENTWELL_HEADER_INCLUDED

#include <opm/models/common/multiphasebaseproperties.hh>

#include <opm/simulators/wells/WellInterface.hpp>
#include <opm/simulators/wells/MultisegmentWellEval.hpp>

namespace Opm {

    class DeferredLogger;

    template<typename TypeTag>
    class MultisegmentWell : public WellInterface<TypeTag>
                           , public MultisegmentWellEval<GetPropType<TypeTag, Properties::FluidSystem>,
                                                         GetPropType<TypeTag, Properties::Indices>>
    {
    public:
        using Base = WellInterface<TypeTag>;
        using MSWEval = MultisegmentWellEval<GetPropType<TypeTag, Properties::FluidSystem>,
                                             GetPropType<TypeTag, Properties::Indices>>;

        using typename Base::Simulator;
        using typename Base::IntensiveQuantities;
        using typename Base::FluidSystem;
        using typename Base::ModelParameters;
        using typename Base::MaterialLaw;
        using typename Base::Indices;
        using typename Base::RateConverterType;
        using typename Base::SparseMatrixAdapter;
        using typename Base::FluidState;

        using Base::has_solvent;
        using Base::has_polymer;
        using Base::Water;
        using Base::Oil;
        using Base::Gas;

        using typename Base::Scalar;

        /// the matrix and vector types for the reservoir
        using typename Base::BVector;
        using typename Base::Eval;

        using typename MSWEval::Equations;
        using typename MSWEval::EvalWell;
        using typename MSWEval::BVectorWell;
        using MSWEval::SPres;
        using typename Base::PressureMatrix;

        MultisegmentWell(const Well& well,
                         const ParallelWellInfo<Scalar>& pw_info,
                         const int time_step,
                         const ModelParameters& param,
                         const RateConverterType& rate_converter,
                         const int pvtRegionIdx,
                         const int num_components,
                         const int num_phases,
                         const int index_of_well,
                         const std::vector<PerforationData<Scalar>>& perf_data);

        void init(const PhaseUsage* phase_usage_arg,
                  const std::vector<Scalar>& depth_arg,
                  const Scalar gravity_arg,
                  const std::vector<Scalar>& B_avg,
                  const bool changed_to_open_this_step) override;

        void initPrimaryVariablesEvaluation() override;

        /// updating the well state based the current control mode
        void updateWellStateWithTarget(const Simulator& simulator,
                                       const GroupState<Scalar>& group_state,
                                       WellState<Scalar>& well_state,
                                       DeferredLogger& deferred_logger, 
                                       const bool initialize = true) override;

        /// check whether the well equations get converged for this well
        ConvergenceReport getWellConvergence(const Simulator& simulator,
                                             const WellState<Scalar>& well_state,
                                             const std::vector<Scalar>& B_avg,
                                             DeferredLogger& deferred_logger,
                                             const bool relax_tolerance) const override;

        /// Ax = Ax - C D^-1 B x
        void apply(const BVector& x, BVector& Ax) const override;
        /// r = r - C D^-1 Rw
        void apply(BVector& r) const override;

        /// using the solution x to recover the solution xw for wells and applying
        /// xw to update Well State
        void recoverWellSolutionAndUpdateWellState(const Simulator& simulator,
                                                   const BVector& x,
                                                   WellState<Scalar>& well_state,
                                                   DeferredLogger& deferred_logger) override;

        /// computing the well potentials for group control
        void computeWellPotentials(const Simulator& simulator,
                                   const WellState<Scalar>& well_state,
                                   std::vector<Scalar>& well_potentials,
                                   DeferredLogger& deferred_logger) override;

        void updatePrimaryVariables(const Simulator& simulator,
                                    const WellState<Scalar>& well_state,
                                    DeferredLogger& deferred_logger) override;

        void solveEqAndUpdateWellState(const Simulator& simulator,
                                       WellState<Scalar>& well_state,
                                       DeferredLogger& deferred_logger) override; // const?

        void calculateExplicitQuantities(const Simulator& simulator,
                                         const WellState<Scalar>& well_state,
                                         DeferredLogger& deferred_logger) override; // should be const?

        void updateIPRImplicit(const Simulator& simulator,
                               WellState<Scalar>& well_state,
                               DeferredLogger& deferred_logger) override;

        void updateProductivityIndex(const Simulator& simulator,
                                     const WellProdIndexCalculator<Scalar>& wellPICalc,
                                     WellState<Scalar>& well_state,
                                     DeferredLogger& deferred_logger) const override;

        Scalar connectionDensity(const int globalConnIdx,
                                 const int openConnIdx) const override;

        void addWellContributions(SparseMatrixAdapter& jacobian) const override;

        void addWellPressureEquations(PressureMatrix& mat,
                                      const BVector& x,
                                      const int pressureVarIndex,
                                      const bool use_well_weights,
                                      const WellState<Scalar>& well_state) const override;

        std::vector<Scalar>
        computeCurrentWellRates(const Simulator& simulator,
                                DeferredLogger& deferred_logger) const override;

        std::optional<Scalar>
        computeBhpAtThpLimitProdWithAlq(const Simulator& simulator,
                                        const SummaryState& summary_state,
                                        const Scalar alq_value,
                                        DeferredLogger& deferred_logger,
                                        bool iterate_if_no_solution) const override;

        std::vector<Scalar> getPrimaryVars() const override;

        int setPrimaryVars(typename std::vector<Scalar>::const_iterator it) override;

    protected:
        // regularize msw equation
        bool regularize_;

        // the intial amount of fluids in each segment under surface condition
        std::vector<std::vector<Scalar> > segment_fluid_initial_;

        mutable int debug_cost_counter_ = 0;

        // updating the well_state based on well solution dwells
        void updateWellState(const Simulator& simulator,
                             const BVectorWell& dwells,
                             WellState<Scalar>& well_state,
                             DeferredLogger& deferred_logger,
                             const Scalar relaxation_factor = 1.0);

        // computing the accumulation term for later use in well mass equations
        void computeInitialSegmentFluids(const Simulator& simulator);

        // compute the pressure difference between the perforation and cell center
        void computePerfCellPressDiffs(const Simulator& simulator);

        template<class Value>
        void computePerfRate(const IntensiveQuantities& int_quants,
                             const std::vector<Value>& mob_perfcells,
                             const std::vector<Scalar>& Tw,
                             const int seg,
                             const int perf,
                             const Value& segment_pressure,
                             const bool& allow_cf,
                             std::vector<Value>& cq_s,
                             Value& perf_press,
                             PerforationRates<Scalar>& perf_rates,
                             DeferredLogger& deferred_logger) const;

        template<class Value>
        void computePerfRate(const Value& pressure_cell,
                        const Value& rs,
                        const Value& rv,
                        const std::vector<Value>& b_perfcells,
                        const std::vector<Value>& mob_perfcells,
                        const std::vector<Scalar>& Tw,
                        const int perf,
                        const Value& segment_pressure,
                        const Value& segment_density,
                        const bool& allow_cf,
                        const std::vector<Value>& cmix_s,
                        std::vector<Value>& cq_s,
                        Value& perf_press,
                        PerforationRates<Scalar>& perf_rates,
                        DeferredLogger& deferred_logger) const;

        // compute the fluid properties, such as densities, viscosities, and so on, in the segments
        // They will be treated implicitly, so they need to be of Evaluation type
        void computeSegmentFluidProperties(const Simulator& simulator,
                                           DeferredLogger& deferred_logger);

        // get the mobility for specific perforation
        template<class Value>
        void getMobility(const Simulator& simulator,
                         const int perf,
                         std::vector<Value>& mob,
                         DeferredLogger& deferred_logger) const;

        void computeWellRatesAtBhpLimit(const Simulator& simulator,
                                        std::vector<Scalar>& well_flux,
                                        DeferredLogger& deferred_logger) const;

        void computeWellRatesWithBhp(const Simulator& simulator,
                                     const Scalar& bhp,
                                     std::vector<Scalar>& well_flux,
                                     DeferredLogger& deferred_logger) const override;

        void computeWellRatesWithBhpIterations(const Simulator& simulator,
                                               const Scalar& bhp,
                                               std::vector<Scalar>& well_flux,
                                               DeferredLogger& deferred_logger) const override;

        std::vector<Scalar>
        computeWellPotentialWithTHP(const WellState<Scalar>& well_state,
                                    const Simulator& simulator,
                                    DeferredLogger& deferred_logger) const;

        bool computeWellPotentialsImplicit(const Simulator& simulator,
                                           const WellState<Scalar>& well_state,
                                           std::vector<Scalar>& well_potentials,
                                           DeferredLogger& deferred_logger) const;

        Scalar getRefDensity() const override;

        bool iterateWellEqWithControl(const Simulator& simulator,
                                      const double dt,
                                      const Well::InjectionControls& inj_controls,
                                      const Well::ProductionControls& prod_controls,
                                      WellState<Scalar>& well_state,
                                      const GroupState<Scalar>& group_state,
                                      DeferredLogger& deferred_logger) override;

        bool iterateWellEqWithSwitching(const Simulator& simulator,
                                        const double dt,
                                        const Well::InjectionControls& inj_controls,
                                        const Well::ProductionControls& prod_controls,
                                        WellState<Scalar>& well_state,
                                        const GroupState<Scalar>& group_state,
                                        DeferredLogger& deferred_logger,
                                        const bool fixed_control = false,
                                        const bool fixed_status = false) override;

        void assembleWellEqWithoutIteration(const Simulator& simulator,
                                            const double dt,
                                            const Well::InjectionControls& inj_controls,
                                            const Well::ProductionControls& prod_controls,
                                            WellState<Scalar>& well_state,
                                            const GroupState<Scalar>& group_state,
                                            DeferredLogger& deferred_logger) override;

        void updateWaterThroughput(const double dt, WellState<Scalar>& well_state) const override;

        EvalWell getSegmentSurfaceVolume(const Simulator& simulator, const int seg_idx) const;

        // turn on crossflow to avoid singular well equations
        // when the well is banned from cross-flow and the BHP is not properly initialized,
        // we turn on crossflow to avoid singular well equations. It can result in wrong-signed
        // well rates, it can cause problem for THP calculation
        // TODO: looking for better alternative to avoid wrong-signed well rates
        bool openCrossFlowAvoidSingularity(const Simulator& simulator) const;

        // for a well, when all drawdown are in the wrong direction, then this well will not
        // be able to produce/inject .
        bool allDrawDownWrongDirection(const Simulator& simulator) const;

        std::optional<Scalar>
        computeBhpAtThpLimitProd(const WellState<Scalar>& well_state,
                                 const Simulator& ebos_simulator,
                                 const SummaryState& summary_state,
                                 DeferredLogger& deferred_logger) const;

        std::optional<Scalar>
        computeBhpAtThpLimitInj(const Simulator& ebos_simulator,
                                const SummaryState& summary_state,
                                DeferredLogger& deferred_logger) const;

        Scalar maxPerfPress(const Simulator& simulator) const;

        // check whether the well is operable under BHP limit with current reservoir condition
        void checkOperabilityUnderBHPLimit(const WellState<Scalar>& well_state,
                                           const Simulator& ebos_simulator,
                                           DeferredLogger& deferred_logger) override;

        // check whether the well is operable under THP limit with current reservoir condition
        void checkOperabilityUnderTHPLimit(const Simulator& ebos_simulator,
                                           const WellState<Scalar>& well_state,
                                           DeferredLogger& deferred_logger) override;

        // updating the inflow based on the current reservoir condition
        void updateIPR(const Simulator& ebos_simulator,
                       DeferredLogger& deferred_logger) const override;
    };

}

#include "MultisegmentWell_impl.hpp"

#endif // OPM_MULTISEGMENTWELL_HEADER_INCLUDED
