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

#include <opm/common/Exceptions.hpp>
#include <opm/simulators/wells/WellInterface.hpp>
#include <opm/simulators/wells/MultisegmentWellEval.hpp>

#include <limits>
#include <string_view>

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
        using typename Base::IndexTraits;
        using typename Base::ModelParameters;
        using typename Base::MaterialLaw;
        using typename Base::Indices;
        using typename Base::RateConverterType;
        using typename Base::SparseMatrixAdapter;
        using typename Base::FluidState;
        using typename Base::WellStateType;
        using typename Base::GroupStateHelperType;

        using Base::has_solvent;
        using Base::has_polymer;
        using Base::has_energy;
        using Base::has_brine;
        using Base::Water;
        using Base::Oil;
        using Base::Gas;

        using typename Base::Scalar;

        // True when the composition switch primary variable is active, i.e. both oil and
        // gas phases are present so that Rs/Rv are stored in the fluid state. Matches the
        // fluid state's enableDissolution flag (see WellInterface::BlackOilFluidStateType).
        static constexpr bool compositionSwitchEnabled =
            Indices::compositionSwitchIdx != std::numeric_limits<unsigned>::max();

        // True when the segment fluid state stores a temperature (any thermal mode, not just
        // the fully implicit one). Matches the fluid state's enableTemperature flag (see
        // WellInterface::BlackOilFluidStateType); used to decide whether createFluidState()
        // must set the temperature explicitly.
        static constexpr bool enable_temperature =
            Base::energyModuleType != EnergyModules::NoTemperature;

        // Scales the well-side energy equation onto the mass-balance residual
        // scale. Reuses the reservoir energy factor so both live on the same scale.
        static constexpr Scalar energy_scaling_factor_ =
            getPropValue<TypeTag, Properties::BlackOilEnergyScalingFactor>();

        /// the matrix and vector types for the reservoir
        using typename Base::BVector;
        using typename Base::Eval;

        using typename MSWEval::Equations;
        using typename MSWEval::EvalWell;
        using typename MSWEval::BVectorWell;
        using MSWEval::SPres;
        using typename Base::PressureMatrix;
        using FSInfo = std::tuple<Scalar, Scalar>;

        using BMatrix = typename Base::BMatrix;
        using CMatrix = typename Base::CMatrix;
        using DMatrix = typename Base::DMatrix;
        using WVector = typename Base::WVector;

        // a fluid state to calculate the properties inside the wellbore for each segment
        // it will be probably used for more things, but at the moment, it is for the enthalpy
        // calculation in the wellbore.
        template <typename ValueType>
        using SegmentFluidState = Base::template BlackOilFluidStateType<ValueType>;

        MultisegmentWell(const Well& well,
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

        /// updating the well state based the current control mode
        void updateWellStateWithTarget(const Simulator& simulator,
                                       const GroupStateHelperType& groupStateHelper,
                                       WellStateType& well_state) const override;

        /// updating the segment pressure and rates based the current bhp and well rates
        void scaleSegmentRatesAndPressure(WellStateType& well_state) const override;

        /// check whether the well equations get converged for this well
        ConvergenceReport getWellConvergence(const GroupStateHelperType& groupStateHelper,
                                             const std::vector<Scalar>& B_avg,
                                             const bool relax_tolerance) const override;

        /// Ax = Ax - C D^-1 B x
        void apply(const BVector& x, BVector& Ax) const override;
        /// r = r - C D^-1 Rw
        void apply(BVector& r) const override;

        /// using the solution x to recover the solution xw for wells and applying
        /// xw to update Well State
        void recoverWellSolutionAndUpdateWellState(const Simulator& simulator,
                                                   const BVector& x,
                                                   const GroupStateHelperType& groupStateHelper,
                                                   WellStateType& well_state) override;

        /// computing the well potentials for group control
        void computeWellPotentials(const Simulator& simulator,
                                   const WellStateType& well_state,
                                   const GroupStateHelperType& groupStateHelper,
                                   std::vector<Scalar>& well_potentials) override;

        void updatePrimaryVariables(const GroupStateHelperType& groupStateHelper) override;

        void solveEqAndUpdateWellState(const Simulator& simulator,
                                       const GroupStateHelperType& groupStateHelper,
                                       WellStateType& well_state) override; // const?

        void calculateExplicitQuantities(const Simulator& simulator,
                                         const GroupStateHelperType& groupStateHelper) override; // should be const?

        void updateIPRImplicit(const Simulator& simulator,
                               const GroupStateHelperType& groupStateHelper,
                               WellStateType& well_state) override;

        void updateProductivityIndex(const Simulator& simulator,
                                     const WellProdIndexCalculator<Scalar>& wellPICalc,
                                     WellStateType& well_state,
                                     DeferredLogger& deferred_logger) const override;

        Scalar connectionDensity(const int globalConnIdx,
                                 const int openConnIdx) const override;

        void addWellContributions(SparseMatrixAdapter& jacobian) const override;

        void addWellPressureEquations(PressureMatrix& mat,
                                      const BVector& x,
                                      const int pressureVarIndex,
                                      const bool use_well_weights,
                                      const WellStateType& well_state) const override;

        std::vector<Scalar>
        computeCurrentWellRates(const Simulator& simulator,
                                DeferredLogger& deferred_logger) const override;

        std::optional<Scalar>
        computeBhpAtThpLimitProdWithAlq(const Simulator& simulator,
                                        const GroupStateHelperType& groupStateHelper,
                                        const SummaryState& summary_state,
                                        const Scalar alq_value,
                                        bool iterate_if_no_solution) const override;

        std::vector<Scalar> getPrimaryVars() const override;

        int setPrimaryVars(typename std::vector<Scalar>::const_iterator it) override;
        void addBCDMatrix(std::vector<typename Base::BMatrix>& b_matrices,
                          std::vector<typename Base::CMatrix>& c_matrices,
                          std::vector<typename Base::DMatrix>& d_matrices,
                          Opm::SparseTable<int>& wcells) const override
        {
            // System_cpr preconditioner is only supported when well DOF dimensions
            // match between WellInterface and MultisegmentWellEval (standard 3-phase blackoil).
            if constexpr (Base::numWellDofs == MSWEval::numWellDofs) {
                MSWEval::addBCDMatrix(b_matrices, c_matrices, d_matrices, wcells);
            } else {
                OPM_THROW(std::runtime_error,
                          "system_cpr preconditioner with multisegment wells is only supported for standard "
                          "3-phase blackoil (Indices::numEq == 3). This model has different equation count.");
            }
        }

        void getScaledWellFractions(std::vector<Scalar>& scaled_fractions,
                                    DeferredLogger& deferred_logger) const override;

    protected:
        // regularize msw equation
        bool regularize_;

        // the initial amount of fluids in each segment under surface condition
        std::vector<std::vector<Scalar> > segment_fluid_initial_;
        // total energy inside the segments at the beginning of the time step
        std::vector<Scalar> segment_initial_energy_;

        // segment fluid state
        std::vector<SegmentFluidState<EvalWell>> segment_fluid_state_;

        // fluid state under the wellhead condition, it is used to calculate the enthalpy
        // under operation condition for energy injection
        // because BHP will be involved, we use EvalWell type here
        SegmentFluidState<EvalWell> wellhead_fluid_state_;

        mutable int debug_cost_counter_ = 0;

        // updating the well_state based on well solution dwells
        void updateWellState(const Simulator& simulator,
                             const BVectorWell& dwells,
                             const GroupStateHelperType& groupStateHelper,
                             WellStateType& well_state,
                             const Scalar relaxation_factor = 1.0);

        // computing the accumulation term for later use in well mass equations
        void computeInitialSegmentFluids(const FSInfo& info, DeferredLogger& deferred_logger);

        // compute the pressure difference between the perforation and cell center
        void computePerfCellPressDiffs(const Simulator& simulator);

        template<class Value>
        void computePerfRate(const IntensiveQuantities& int_quants,
                             const std::vector<Value>& mob_perfcells,
                             const std::vector<Value>& Tw,
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
                        const std::vector<Value>& Tw,
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

        // get the transmissibility multiplier for specific perforation
        template<class Value>
        void getTransMult(Value& trans_mult,
                          const Simulator& simulator,
                          const int cell_indx) const;

        // get the mobility for specific perforation
        template<class Value>
        void getMobility(const Simulator& simulator,
                         const int local_perf_index,
                         std::vector<Value>& mob,
                         DeferredLogger& deferred_logger) const;

        void computeWellRatesAtBhpLimit(const Simulator& simulator,
                                        const GroupStateHelperType& groupStateHelper,
                                        std::vector<Scalar>& well_flux) const;

        void computeWellRatesWithBhp(const Simulator& simulator,
                                     const Scalar& bhp,
                                     std::vector<Scalar>& well_flux,
                                     DeferredLogger& deferred_logger) const override;

        void computeWellRatesWithBhpIterations(const Simulator& simulator,
                                               const Scalar& bhp,
                                               const GroupStateHelperType& groupStateHelper,
                                               std::vector<Scalar>& well_flux) const override;

        std::vector<Scalar>
        computeWellPotentialWithTHP(const WellStateType& well_state,
                                    const Simulator& simulator,
                                    const GroupStateHelperType& groupStateHelper) const;

        bool computeWellPotentialsImplicit(const Simulator& simulator,
                                           const GroupStateHelperType& groupStateHelper,
                                           std::vector<Scalar>& well_potentials) const;

        Scalar getRefDensity() const override;

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

        void assembleWellEqWithoutIteration(const Simulator& simulator,
                                            const GroupStateHelperType& groupStateHelper,
                                            const double dt,
                                            const Well::InjectionControls& inj_controls,
                                            const Well::ProductionControls& prod_controls,
                                            WellStateType& well_state,
                                            const bool solving_with_zero_rate) override;

        void updateWaterThroughput(const double dt, WellStateType& well_state) const override;

        EvalWell getSegmentSurfaceVolume(const int seg_idx,
                                         const FSInfo& info,
                                         DeferredLogger& deferred_logger) const;

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
        computeBhpAtThpLimitProd(const WellStateType& well_state,
                                 const Simulator& ebos_simulator,
                                 const GroupStateHelperType& groupStateHelper,
                                 const SummaryState& summary_state) const;

        std::optional<Scalar>
        computeBhpAtThpLimitInj(const Simulator& ebos_simulator,
                                const GroupStateHelperType& groupStateHelper,
                                const SummaryState& summary_state) const;

        Scalar maxPerfPress(const Simulator& simulator) const override;

        // check whether the well is operable under BHP limit with current reservoir condition
        void checkOperabilityUnderBHPLimit(const WellStateType& well_state,
                                           const Simulator& ebos_simulator,
                                           DeferredLogger& deferred_logger) override;

        // check whether the well is operable under THP limit with current reservoir condition
        void checkOperabilityUnderTHPLimit(const Simulator& ebos_simulator,
                                           const WellStateType& well_state,
                                           const GroupStateHelperType& groupStateHelper) override;

        // updating the inflow based on the current reservoir condition
        void updateIPR(const Simulator& ebos_simulator,
                       DeferredLogger& deferred_logger) const override;

        FSInfo getFirstPerforationFluidStateInfo(const Simulator& simulator) const;

        // this function can potentially be shared between multisegment wells and standard wells
        // The optional @p volume_ratio output receives the segment volume ratio (reservoir
        // volume per unit surface volume, i.e. the sum of the unnormalized phase saturations),
        // so callers can cache and reuse it for the segment surface volume.
        template <typename ValueType = EvalWell>
        SegmentFluidState<ValueType>
        createFluidState(const std::vector<ValueType>& fluid_composition,
                         const ValueType& pressure,
                         const ValueType& temperature,
                         const ValueType& saltConcentration,
                         DeferredLogger& deferred_logger,
                         ValueType* volume_ratio = nullptr) const;

        SegmentFluidState<EvalWell>
        createSegmentFluidState(int seg, const FSInfo& info, DeferredLogger& deferred_logger,
                                EvalWell* volume_ratio = nullptr) const;

        void computeInitialSegmentEnergy();

        // assemble the energy equation contribution for a single perforation/connection
        void assemblePerforationEnergyEq(const IntensiveQuantities& int_quants,
                                         const std::vector<EvalWell>& cq_s,
                                         const int seg,
                                         const int local_perf_index,
                                         DeferredLogger& deferred_logger);

        void updateWellHeadCondition(const Simulator& simulator,
                                     const Scalar first_perf_temperature,
                                     const Scalar first_perf_salt_concentration,
                                     DeferredLogger& deferred_logger);

        void updateSegmentFluidState(const FSInfo& info, DeferredLogger& deferred_logger);

        template <typename ValueType = EvalWell>
        ValueType computeSegmentEnergy(int seg) const;

        // Convert per-component surface volumetric rates to a phase reservoir
        // volumetric rate, using @p fs (upwind) for the rs/rv coupling. If
        // (1 - rs*rv) <= 0, falls back to rate / invB (drops the cross-terms but
        // keeps @p fs's invB; Rs/Rv unchanged) and logs a @p context debug message.
        // @p fs may be a wellbore SegmentFluidState (EvalWell) or a reservoir-cell
        // fluid state (Eval, extended to EvalWell on the fly).
        // @return the reservoir volumetric rate of @p phaseIdx.
        template <typename FluidStateT>
        EvalWell surfaceToReservoirRate(unsigned phaseIdx,
                                        const FluidStateT& fs,
                                        const std::vector<EvalWell>& surface_rates,
                                        int seg,
                                        std::string_view context,
                                        DeferredLogger& deferred_logger) const;

        // Compute the energy flux carried by fluid flowing from segment
        // @p seg toward its outlet, using @p upwind_fs as the upwind
        // fluid state. @p context is used in diagnostic messages.
        // @return the energy flux as sum_phases(reservoir_rate * enthalpy * density).
        EvalWell computeSegmentEnergyRate(int seg,
                                          int upwind_seg,
                                          const SegmentFluidState<EvalWell>& upwind_fs,
                                          std::string_view context,
                                          DeferredLogger& deferred_logger) const;
    };

} // namespace Opm

#include "MultisegmentWell_impl.hpp"

#endif // OPM_MULTISEGMENTWELL_HEADER_INCLUDED
