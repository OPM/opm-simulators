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

#include <opm/parser/eclipse/EclipseState/Runspec.hpp>

namespace Opm
{

    template<typename TypeTag>
    class MultisegmentWell: public WellInterface<TypeTag>
    {
    public:
        typedef WellInterface<TypeTag> Base;

        using typename Base::WellState;
        using typename Base::Simulator;
        using typename Base::IntensiveQuantities;
        using typename Base::FluidSystem;
        using typename Base::ModelParameters;
        using typename Base::MaterialLaw;
        using typename Base::Indices;
        using typename Base::RateConverterType;
        using typename Base::SparseMatrixAdapter;
        using typename Base::FluidState;
        using typename Base::GasLiftSingleWell;
        using typename Base::GLiftProdWells;
        using typename Base::GLiftOptWells;
        using typename Base::GLiftWellStateMap;

        /// the number of reservior equations
        using Base::numEq;
        using Base::numPhases;

        using Base::has_solvent;
        using Base::has_polymer;
        using Base::Water;
        using Base::Oil;
        using Base::Gas;

        // TODO: for now, not considering the polymer, solvent and so on to simplify the development process.

        // TODO: we need to have order for the primary variables and also the order for the well equations.
        // sometimes, they are similar, while sometimes, they can have very different forms.

        static constexpr bool has_gas = (Indices::compositionSwitchIdx >= 0);
        static constexpr bool has_water = (Indices::waterSaturationIdx >= 0);

        static constexpr int GTotal = 0;
        static constexpr int WFrac = has_water ? 1: -1000;
        static constexpr int GFrac = has_gas ? has_water + 1 : -1000;
        static constexpr int SPres = has_gas + has_water + 1;

        //  the number of well equations  TODO: it should have a more general strategy for it
        static const int numWellEq = numPhases + 1;

        using typename Base::Scalar;

        /// the matrix and vector types for the reservoir
        using typename Base::BVector;
        using typename Base::Eval;

        // sparsity pattern for the matrices
        // [A C^T    [x       =  [ res
        //  B  D ]   x_well]      res_well]

        // the vector type for the res_well and x_well
        typedef Dune::FieldVector<Scalar, numWellEq> VectorBlockWellType;
        typedef Dune::BlockVector<VectorBlockWellType> BVectorWell;

        // the matrix type for the diagonal matrix D
        typedef Dune::FieldMatrix<Scalar, numWellEq, numWellEq > DiagMatrixBlockWellType;
        typedef Dune::BCRSMatrix <DiagMatrixBlockWellType> DiagMatWell;

        // the matrix type for the non-diagonal matrix B and C^T
        typedef Dune::FieldMatrix<Scalar, numWellEq, numEq>  OffDiagMatrixBlockWellType;
        typedef Dune::BCRSMatrix<OffDiagMatrixBlockWellType> OffDiagMatWell;

        // TODO: for more efficient implementation, we should have EvalReservoir, EvalWell, and EvalRerservoirAndWell
        //                                                         EvalR (Eval), EvalW, EvalRW
        // TODO: for now, we only use one type to save some implementation efforts, while improve later.
        typedef DenseAd::Evaluation<double, /*size=*/numEq + numWellEq> EvalWell;

        MultisegmentWell(const Well& well,
                         const ParallelWellInfo& pw_info,
                         const int time_step,
                         const ModelParameters& param,
                         const RateConverterType& rate_converter,
                         const int pvtRegionIdx,
                         const int num_components,
                         const int num_phases,
                         const int index_of_well,
                         const int first_perf_index,
                         const std::vector<PerforationData>& perf_data);

        virtual void init(const PhaseUsage* phase_usage_arg,
                          const std::vector<double>& depth_arg,
                          const double gravity_arg,
                          const int num_cells,
                          const std::vector< Scalar >& B_avg) override;


        virtual void initPrimaryVariablesEvaluation() const override;

        virtual void gasLiftOptimizationStage1 (
            WellState&,
            const Simulator&,
            DeferredLogger&,
            GLiftProdWells &,
            GLiftOptWells &,
            GLiftWellStateMap &
        ) const override {
            // Not implemented yet
        }

        virtual void assembleWellEq(const Simulator& ebosSimulator,
                                    const double dt,
                                    WellState& well_state,
                                    const GroupState& group_state,
                                    Opm::DeferredLogger& deferred_logger) override;

        /// updating the well state based the current control mode
        virtual void updateWellStateWithTarget(const Simulator& ebos_simulator,
                                               WellState& well_state,
                                               Opm::DeferredLogger& deferred_logger) const override;

        /// check whether the well equations get converged for this well
        virtual ConvergenceReport getWellConvergence(const WellState& well_state, const std::vector<double>& B_avg, Opm::DeferredLogger& deferred_logger, const bool relax_tolerance = false) const override;

        /// Ax = Ax - C D^-1 B x
        virtual void apply(const BVector& x, BVector& Ax) const override;
        /// r = r - C D^-1 Rw
        virtual void apply(BVector& r) const override;

#if HAVE_CUDA || HAVE_OPENCL
        /// add the contribution (C, D, B matrices) of this Well to the WellContributions object
        void addWellContribution(WellContributions& wellContribs) const;
#endif

        /// using the solution x to recover the solution xw for wells and applying
        /// xw to update Well State
        virtual void recoverWellSolutionAndUpdateWellState(const BVector& x,
                                                           WellState& well_state,
                                                           Opm::DeferredLogger& deferred_logger) const override;

        /// computing the well potentials for group control
        virtual void computeWellPotentials(const Simulator& ebosSimulator,
                                           const WellState& well_state,
                                           std::vector<double>& well_potentials,
                                           Opm::DeferredLogger& deferred_logger) override;

        virtual void updatePrimaryVariables(const WellState& well_state, Opm::DeferredLogger& deferred_logger) const override;

        virtual void solveEqAndUpdateWellState(WellState& well_state, Opm::DeferredLogger& deferred_logger) override; // const?

        virtual void calculateExplicitQuantities(const Simulator& ebosSimulator,
                                                 const WellState& well_state,
                                                 Opm::DeferredLogger& deferred_logger) override; // should be const?

        virtual void updateProductivityIndex(const Simulator& ebosSimulator,
                                             const WellProdIndexCalculator& wellPICalc,
                                             WellState& well_state,
                                             DeferredLogger& deferred_logger) const override;

        virtual void  addWellContributions(SparseMatrixAdapter& jacobian) const override;

        /// number of segments for this well
        /// int number_of_segments_;
        int numberOfSegments() const;

        int numberOfPerforations() const;

        virtual std::vector<double> computeCurrentWellRates(const Simulator& ebosSimulator,
                                                            DeferredLogger& deferred_logger) const override;

        void computeConnLevelProdInd(const FluidState& fs,
                                     const std::function<double(const double)>& connPICalc,
                                     const std::vector<EvalWell>& mobility,
                                     double* connPI) const;

        void computeConnLevelInjInd(const FluidState& fs,
                                    const Phase preferred_phase,
                                    const std::function<double(const double)>& connIICalc,
                                    const std::vector<EvalWell>& mobility,
                                    double* connII,
                                    DeferredLogger& deferred_logger) const;
    protected:
        int number_segments_;

        // components of the pressure drop to be included
        WellSegments::CompPressureDrop compPressureDrop() const;
        // multi-phase flow model
        WellSegments::MultiPhaseModel multiphaseModel() const;

        // get the WellSegments from the well_ecl_
        const WellSegments& segmentSet() const;

        // protected member variables from the Base class
        using Base::well_ecl_;
        using Base::vfp_properties_;
        using Base::ref_depth_;
        using Base::number_of_perforations_; // TODO: can use well_ecl_?
        using Base::current_step_;
        using Base::index_of_well_;
        using Base::number_of_phases_;

        // TODO: the current implementation really relies on the order of the
        // perforation does not change from the parser to Wells structure.
        using Base::well_cells_;
        using Base::param_;
        using Base::well_index_;
        using Base::first_perf_;
        using Base::saturation_table_number_;
        using Base::well_efficiency_factor_;
        using Base::gravity_;
        using Base::perf_depth_;
        using Base::num_components_;
        using Base::connectionRates_;
        using Base::ipr_a_;
        using Base::ipr_b_;
        using Base::changed_to_stopped_this_step_;

        // protected functions from the Base class
        using Base::phaseUsage;
        using Base::name;
        using Base::flowPhaseToEbosCompIdx;
        using Base::flowPhaseToEbosPhaseIdx;
        using Base::ebosCompIdxToFlowCompIdx;
        using Base::getAllowCrossFlow;
        using Base::scalingFactor;
        using Base::wellIsStopped;
        using Base::updateWellOperability;
        using Base::checkWellOperability;

        // TODO: trying to use the information from the Well opm-parser as much
        // as possible, it will possibly be re-implemented later for efficiency reason.

        // the completions that is related to each segment
        // the completions's ids are their index in the vector well_index_, well_cell_
        // This is also assuming the order of the completions in Well is the same with
        // the order of the completions in wells.
        // it is for convinience reason. we can just calcuate the inforation for segment once then using it for all the perofrations
        // belonging to this segment
        std::vector<std::vector<int> > segment_perforations_;

        // the inlet segments for each segment. It is for convinience and efficiency reason
        std::vector<std::vector<int> > segment_inlets_;

        // segment number is an ID of the segment, it is specified in the deck
        // get the loation of the segment with a segment number in the segmentSet
        int segmentNumberToIndex(const int segment_number) const;

        // TODO, the following should go to a class for computing purpose
        // two off-diagonal matrices
        mutable OffDiagMatWell duneB_;
        mutable OffDiagMatWell duneC_;
        // "diagonal" matrix for the well. It has offdiagonal entries for inlets and outlets.
        mutable DiagMatWell duneD_;
        /// \brief solver for diagonal matrix
        ///
        /// This is a shared_ptr as MultisegmentWell is copied in computeWellPotentials...
        mutable std::shared_ptr<Dune::UMFPack<DiagMatWell> > duneDSolver_;

        // residuals of the well equations
        mutable BVectorWell resWell_;

        // the values for the primary varibles
        // based on different solutioin strategies, the wells can have different primary variables
        mutable std::vector<std::array<double, numWellEq> > primary_variables_;

        // the Evaluation for the well primary variables, which contain derivativles and are used in AD calculation
        mutable std::vector<std::array<EvalWell, numWellEq> > primary_variables_evaluation_;

        // depth difference between perforations and the perforated grid cells
        std::vector<double> cell_perforation_depth_diffs_;
        // pressure correction due to the different depth of the perforation and
        // center depth of the grid block
        std::vector<double> cell_perforation_pressure_diffs_;

        // depth difference between the segment and the peforation
        // or in another way, the depth difference between the perforation and
        // the segment the perforation belongs to
        std::vector<double> perforation_segment_depth_diffs_;

        // the intial amount of fluids in each segment under surface condition
        std::vector<std::vector<double> > segment_fluid_initial_;

        // the densities of segment fluids
        // we should not have this member variable
        std::vector<EvalWell> segment_densities_;

        // the viscosity of the segments
        std::vector<EvalWell> segment_viscosities_;

        // the mass rate of the segments
        std::vector<EvalWell> segment_mass_rates_;

        std::vector<double> segment_depth_diffs_;

        // the upwinding segment for each segment based on the flow direction
        std::vector<int> upwinding_segments_;

        mutable int debug_cost_counter_ = 0;

        std::vector<std::vector<EvalWell>> segment_phase_fractions_;

        std::vector<std::vector<EvalWell>> segment_phase_viscosities_;

        std::vector<std::vector<EvalWell>> segment_phase_densities_;


        void initMatrixAndVectors(const int num_cells) const;

        EvalWell getBhp() const;
        EvalWell getQs(const int comp_idx) const;
        EvalWell getWQTotal() const;

        // xw = inv(D)*(rw - C*x)
        void recoverSolutionWell(const BVector& x, BVectorWell& xw) const;

        // updating the well_state based on well solution dwells
        void updateWellState(const BVectorWell& dwells,
                             WellState& well_state,
                             Opm::DeferredLogger& deferred_logger,
                             const double relaxation_factor=1.0) const;


        // scale the segment rates and pressure based on well rates and bhp
        void scaleSegmentRatesWithWellRates(WellState& well_state) const;
        void scaleSegmentPressuresWithBhp(WellState& well_state) const;

        // computing the accumulation term for later use in well mass equations
        void computeInitialSegmentFluids(const Simulator& ebos_simulator);

        // compute the pressure difference between the perforation and cell center
        void computePerfCellPressDiffs(const Simulator& ebosSimulator);

        // fraction value of the primary variables
        // should we just use member variables to store them instead of calculating them again and again
        EvalWell volumeFraction(const int seg, const unsigned comp_idx) const;

        // F_p / g_p, the basic usage of this value is because Q_p = G_t * F_p / G_p
        EvalWell volumeFractionScaled(const int seg, const int comp_idx) const;

        // basically Q_p / \sigma_p Q_p
        EvalWell surfaceVolumeFraction(const int seg, const int comp_idx) const;

        void computePerfRatePressure(const IntensiveQuantities& int_quants,
                                     const std::vector<EvalWell>& mob_perfcells,
                                     const double Tw,
                                     const int seg,
                                     const int perf,
                                     const EvalWell& segment_pressure,
                                     const bool& allow_cf,
                                     std::vector<EvalWell>& cq_s,
                                     EvalWell& perf_press,
                                     double& perf_dis_gas_rate,
                                     double& perf_vap_oil_rate,
                                     Opm::DeferredLogger& deferred_logger) const;

        // convert a Eval from reservoir to contain the derivative related to wells
        EvalWell extendEval(const Eval& in) const;


        template <class ValueType>
        ValueType calculateBhpFromThp(const std::vector<ValueType>& rates, const Well& well, const SummaryState& summaryState, Opm::DeferredLogger& deferred_logger) const;

        double calculateThpFromBhp(const std::vector<double>& rates, const double bhp, Opm::DeferredLogger& deferred_logger) const;
        void updateThp(WellState& well_state, Opm::DeferredLogger& deferred_logger) const;

        // compute the fluid properties, such as densities, viscosities, and so on, in the segments
        // They will be treated implicitly, so they need to be of Evaluation type
        void computeSegmentFluidProperties(const Simulator& ebosSimulator);

        EvalWell getSegmentPressure(const int seg) const;

        EvalWell getSegmentRate(const int seg, const int comp_idx) const;

        EvalWell getSegmentRateUpwinding(const int seg, const size_t comp_idx) const;

        EvalWell getSegmentGTotal(const int seg) const;

        // get the mobility for specific perforation
        void getMobility(const Simulator& ebosSimulator,
                         const int perf,
                         std::vector<EvalWell>& mob) const;

        void computeWellRatesAtBhpLimit(const Simulator& ebosSimulator,
                                        std::vector<double>& well_flux,
                                        Opm::DeferredLogger& deferred_logger) const;

        void computeWellRatesWithBhp(const Simulator& ebosSimulator,
                                     const Scalar bhp,
                                     std::vector<double>& well_flux,
                                     Opm::DeferredLogger& deferred_logger) const;

        std::vector<double>
        computeWellPotentialWithTHP(const Simulator& ebos_simulator,
                                    Opm::DeferredLogger& deferred_logger) const;

        void assembleControlEq(const WellState& well_state,
                               const GroupState& group_state,
                               const Opm::Schedule& schedule,
                               const SummaryState& summaryState,
                               const Well::InjectionControls& inj_controls,
                               const Well::ProductionControls& prod_controls,
                               Opm::DeferredLogger& deferred_logger);

        void assemblePressureEq(const int seg, const UnitSystem& unit_system,
                                WellState& well_state, DeferredLogger& deferred_logger) const;

        void assembleDefaultPressureEq(const int seg, WellState& well_state) const;

        // hytrostatic pressure loss
        EvalWell getHydroPressureLoss(const int seg) const;

        // frictinal pressure loss
        EvalWell getFrictionPressureLoss(const int seg) const;

        void handleAccelerationPressureLoss(const int seg, WellState& well_state) const;

        // handling the overshooting and undershooting of the fractions
        void processFractions(const int seg) const;

        void updateWellStateFromPrimaryVariables(WellState& well_state, Opm::DeferredLogger& deferred_logger) const;

        bool frictionalPressureLossConsidered() const;

        bool accelerationalPressureLossConsidered() const;

        virtual bool iterateWellEqWithControl(const Simulator& ebosSimulator,
                                              const double dt,
                                              const Well::InjectionControls& inj_controls,
                                              const Well::ProductionControls& prod_controls,
                                              WellState& well_state,
                                              const GroupState& group_state,
                                              Opm::DeferredLogger& deferred_logger) override;

        virtual void assembleWellEqWithoutIteration(const Simulator& ebosSimulator,
                                                    const double dt,
                                                    const Well::InjectionControls& inj_controls,
                                                    const Well::ProductionControls& prod_controls,
                                                    WellState& well_state,
                                                    const GroupState& group_state,
                                                    Opm::DeferredLogger& deferred_logger) override;

        virtual void updateWaterThroughput(const double dt, WellState& well_state) const override;

        EvalWell getSegmentSurfaceVolume(const Simulator& ebos_simulator, const int seg_idx) const;

        std::vector<Scalar> getWellResiduals(const std::vector<Scalar>& B_avg,
                                             DeferredLogger& deferred_logger) const;

        void detectOscillations(const std::vector<double>& measure_history,
                                const int it, bool& oscillate, bool& stagnate) const;

        double getResidualMeasureValue(const WellState& well_state,
                                       const std::vector<double>& residuals,
                                       DeferredLogger& deferred_logger) const;

        double getControlTolerance(const WellState& well_state, DeferredLogger& deferred_logger) const;

        void checkConvergenceControlEq(const WellState& well_state,
                                       ConvergenceReport& report,
                                       DeferredLogger& deferred_logger) const;

        void updateUpwindingSegments();

        // turn on crossflow to avoid singular well equations
        // when the well is banned from cross-flow and the BHP is not properly initialized,
        // we turn on crossflow to avoid singular well equations. It can result in wrong-signed
        // well rates, it can cause problem for THP calculation
        // TODO: looking for better alternative to avoid wrong-signed well rates
        bool openCrossFlowAvoidSingularity(const Simulator& ebos_simulator) const;

        // for a well, when all drawdown are in the wrong direction, then this well will not
        // be able to produce/inject .
        bool allDrawDownWrongDirection(const Simulator& ebos_simulator) const;


        std::optional<double> computeBhpAtThpLimitProd(const Simulator& ebos_simulator,
                                                       const SummaryState& summary_state,
                                                       DeferredLogger& deferred_logger) const;

        std::optional<double> computeBhpAtThpLimitInj(const Simulator& ebos_simulator,
                                                      const SummaryState& summary_state,
                                                      DeferredLogger& deferred_logger) const;

        double maxPerfPress(const Simulator& ebos_simulator) const;

        // pressure drop for Spiral ICD segment (WSEGSICD)
        EvalWell pressureDropSpiralICD(const int seg) const;

        // pressure drop for Autonomous ICD segment (WSEGAICD)
        EvalWell pressureDropAutoICD(const int seg, const UnitSystem& unit_system) const;

        // pressure drop for sub-critical valve (WSEGVALV)
        EvalWell pressureDropValve(const int seg) const;

        // assemble pressure equation for ICD segments
        void assembleICDPressureEq(const int seg, const UnitSystem& unit_system,
                                   WellState& well_state, DeferredLogger& deferred_logger) const;

        // check whether the well is operable under BHP limit with current reservoir condition
        virtual void checkOperabilityUnderBHPLimitProducer(const WellState& well_state, const Simulator& ebos_simulator, Opm::DeferredLogger& deferred_logger) override;

        // check whether the well is operable under THP limit with current reservoir condition
        virtual void checkOperabilityUnderTHPLimitProducer(const Simulator& ebos_simulator, const WellState& well_state, Opm::DeferredLogger& deferred_logger) override;

        // updating the inflow based on the current reservoir condition
        virtual void updateIPR(const Simulator& ebos_simulator, Opm::DeferredLogger& deferred_logger) const override;

    };

}

#include "MultisegmentWell_impl.hpp"

#endif // OPM_MULTISEGMENTWELL_HEADER_INCLUDED
