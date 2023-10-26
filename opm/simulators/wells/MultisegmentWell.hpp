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

#include <opm/simulators/wells/WellInterface.hpp>
#include <opm/simulators/wells/MultisegmentWellEval.hpp>

namespace Opm
{
    class DeferredLogger;

    template<typename TypeTag>
    class MultisegmentWell : public WellInterface<TypeTag>
                           , public MultisegmentWellEval<GetPropType<TypeTag, Properties::FluidSystem>,
                                                         GetPropType<TypeTag, Properties::Indices>,
                                                         GetPropType<TypeTag, Properties::Scalar>>
    {
    public:
        using Base = WellInterface<TypeTag>;
        using MSWEval = MultisegmentWellEval<GetPropType<TypeTag, Properties::FluidSystem>,
                                             GetPropType<TypeTag, Properties::Indices>,
                                             GetPropType<TypeTag, Properties::Scalar>>;

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
                         const ParallelWellInfo& pw_info,
                         const int time_step,
                         const ModelParameters& param,
                         const RateConverterType& rate_converter,
                         const int pvtRegionIdx,
                         const int num_components,
                         const int num_phases,
                         const int index_of_well,
                         const std::vector<PerforationData>& perf_data);

        virtual void init(const PhaseUsage* phase_usage_arg,
                          const std::vector<double>& depth_arg,
                          const double gravity_arg,
                          const int num_cells,
                          const std::vector< Scalar >& B_avg,
                          const bool changed_to_open_this_step) override;

        void initPrimaryVariablesEvaluation() override;

        /// updating the well state based the current control mode
        virtual void updateWellStateWithTarget(const Simulator& ebos_simulator,
                                               const GroupState& group_state,
                                               WellState& well_state,
                                               DeferredLogger& deferred_logger) const override;

        /// check whether the well equations get converged for this well
        virtual ConvergenceReport getWellConvergence(const SummaryState& summary_state,
                                                     const WellState& well_state,
                                                     const std::vector<double>& B_avg,
                                                     DeferredLogger& deferred_logger,
                                                     const bool relax_tolerance) const override;

        /// Ax = Ax - C D^-1 B x
        virtual void apply(const BVector& x, BVector& Ax) const override;
        /// r = r - C D^-1 Rw
        virtual void apply(BVector& r) const override;

        /// using the solution x to recover the solution xw for wells and applying
        /// xw to update Well State
        void recoverWellSolutionAndUpdateWellState(const SummaryState& summary_state,
                                                   const BVector& x,
                                                   WellState& well_state,
                                                   DeferredLogger& deferred_logger) override;

        /// computing the well potentials for group control
        virtual void computeWellPotentials(const Simulator& ebosSimulator,
                                           const WellState& well_state,
                                           std::vector<double>& well_potentials,
                                           DeferredLogger& deferred_logger) override;

        void updatePrimaryVariables(const SummaryState& summary_state,
                                    const WellState& well_state,
                                    DeferredLogger& deferred_logger) override;

        virtual void solveEqAndUpdateWellState(const SummaryState& summary_state,
                                               WellState& well_state,
                                               DeferredLogger& deferred_logger) override; // const?

        virtual void calculateExplicitQuantities(const Simulator& ebosSimulator,
                                                 const WellState& well_state,
                                                 DeferredLogger& deferred_logger) override; // should be const?

        void updateIPRImplicit(const Simulator& ebos_simulator, WellState& well_state, DeferredLogger& deferred_logger);                                                 

        virtual void updateProductivityIndex(const Simulator& ebosSimulator,
                                             const WellProdIndexCalculator& wellPICalc,
                                             WellState& well_state,
                                             DeferredLogger& deferred_logger) const override;

        double connectionDensity(const int globalConnIdx,
                                 const int openConnIdx) const override;

        void addWellContributions(SparseMatrixAdapter& jacobian) const override;

        void addWellPressureEquations(PressureMatrix& mat,
                                      const BVector& x,
                                      const int pressureVarIndex,
                                      const bool use_well_weights,
                                      const WellState& well_state) const override;

        virtual std::vector<double> computeCurrentWellRates(const Simulator& ebosSimulator,
                                                            DeferredLogger& deferred_logger) const override;

        std::optional<double>
        computeBhpAtThpLimitProdWithAlq(const Simulator& ebos_simulator,
                                        const SummaryState& summary_state,
                                        const double alq_value,
                                        DeferredLogger& deferred_logger) const override;

        std::vector<double> getPrimaryVars() const override;

        int setPrimaryVars(std::vector<double>::const_iterator it) override;

    protected:

        // regularize msw equation
        bool regularize_;

        // the intial amount of fluids in each segment under surface condition
        std::vector<std::vector<double> > segment_fluid_initial_;

        mutable int debug_cost_counter_ = 0;

        // updating the well_state based on well solution dwells
        void updateWellState(const SummaryState& summary_state,
                             const BVectorWell& dwells,
                             WellState& well_state,
                             DeferredLogger& deferred_logger,
                             const double relaxation_factor = 1.0);


        // computing the accumulation term for later use in well mass equations
        void computeInitialSegmentFluids(const Simulator& ebos_simulator);

        // compute the pressure difference between the perforation and cell center
        void computePerfCellPressDiffs(const Simulator& ebosSimulator);

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
                             PerforationRates& perf_rates,
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
                        PerforationRates& perf_rates,
                        DeferredLogger& deferred_logger) const;

        // compute the fluid properties, such as densities, viscosities, and so on, in the segments
        // They will be treated implicitly, so they need to be of Evaluation type
        void computeSegmentFluidProperties(const Simulator& ebosSimulator,
                                           DeferredLogger& deferred_logger);

        // get the mobility for specific perforation
        template<class Value>
        void getMobility(const Simulator& ebosSimulator,
                         const int perf,
                         std::vector<Value>& mob,
                         DeferredLogger& deferred_logger) const;

        void computeWellRatesAtBhpLimit(const Simulator& ebosSimulator,
                                        std::vector<double>& well_flux,
                                        DeferredLogger& deferred_logger) const;

        void computeWellRatesWithBhp(const Simulator& ebosSimulator,
                                     const double& bhp,
                                     std::vector<double>& well_flux,
                                     DeferredLogger& deferred_logger) const override;

        void computeWellRatesWithBhpIterations(const Simulator& ebosSimulator,
                                               const Scalar& bhp,
                                               std::vector<double>& well_flux,
                                               DeferredLogger& deferred_logger) const override;

        std::vector<double> computeWellPotentialWithTHP(
                                 const WellState& well_state,
                                 const Simulator& ebos_simulator,
                                 DeferredLogger& deferred_logger) const;

        bool computeWellPotentialsImplicit(const Simulator& ebos_simulator,
                                           std::vector<double>& well_potentials,
                                           DeferredLogger& deferred_logger) const;                                      

        virtual double getRefDensity() const override;

        virtual bool iterateWellEqWithControl(const Simulator& ebosSimulator,
                                              const double dt,
                                              const Well::InjectionControls& inj_controls,
                                              const Well::ProductionControls& prod_controls,
                                              WellState& well_state,
                                              const GroupState& group_state,
                                              DeferredLogger& deferred_logger) override;

        virtual bool iterateWellEqWithSwitching(const Simulator& ebosSimulator,
                                                const double dt,
                                                const Well::InjectionControls& inj_controls,
                                                const Well::ProductionControls& prod_controls,
                                                WellState& well_state,
                                                const GroupState& group_state,
                                                DeferredLogger& deferred_logger, 
                                                const bool allow_switch = true) override;

        virtual void assembleWellEqWithoutIteration(const Simulator& ebosSimulator,
                                                    const double dt,
                                                    const Well::InjectionControls& inj_controls,
                                                    const Well::ProductionControls& prod_controls,
                                                    WellState& well_state,
                                                    const GroupState& group_state,
                                                    DeferredLogger& deferred_logger) override;

        virtual void updateWaterThroughput(const double dt, WellState& well_state) const override;

        EvalWell getSegmentSurfaceVolume(const Simulator& ebos_simulator, const int seg_idx) const;

        // turn on crossflow to avoid singular well equations
        // when the well is banned from cross-flow and the BHP is not properly initialized,
        // we turn on crossflow to avoid singular well equations. It can result in wrong-signed
        // well rates, it can cause problem for THP calculation
        // TODO: looking for better alternative to avoid wrong-signed well rates
        bool openCrossFlowAvoidSingularity(const Simulator& ebos_simulator) const;

        // for a well, when all drawdown are in the wrong direction, then this well will not
        // be able to produce/inject .
        bool allDrawDownWrongDirection(const Simulator& ebos_simulator) const;



        std::optional<double> computeBhpAtThpLimitProd(
            const WellState& well_state,
            const Simulator& ebos_simulator,
            const SummaryState& summary_state,
            DeferredLogger& deferred_logger) const;

        std::optional<double> computeBhpAtThpLimitInj(const Simulator& ebos_simulator,
                                                      const SummaryState& summary_state,
                                                      DeferredLogger& deferred_logger) const;

        double maxPerfPress(const Simulator& ebos_simulator) const;

        // check whether the well is operable under BHP limit with current reservoir condition
        virtual void checkOperabilityUnderBHPLimit(const WellState& well_state, const Simulator& ebos_simulator, DeferredLogger& deferred_logger) override;

        // check whether the well is operable under THP limit with current reservoir condition
        virtual void checkOperabilityUnderTHPLimit(const Simulator& ebos_simulator, const WellState& well_state, DeferredLogger& deferred_logger) override;

        // updating the inflow based on the current reservoir condition
        virtual void updateIPR(const Simulator& ebos_simulator, DeferredLogger& deferred_logger) const override;
    };

}

#include "MultisegmentWell_impl.hpp"

#endif // OPM_MULTISEGMENTWELL_HEADER_INCLUDED
