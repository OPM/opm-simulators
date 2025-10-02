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


#ifndef OPM_STANDARDWELL_HEADER_INCLUDED
#define OPM_STANDARDWELL_HEADER_INCLUDED

#include <opm/simulators/timestepping/ConvergenceReport.hpp>
#include <opm/simulators/wells/RateConverter.hpp>
#include <opm/simulators/wells/RatioCalculator.hpp>
#include <opm/simulators/wells/VFPInjProperties.hpp>
#include <opm/simulators/wells/VFPProdProperties.hpp>
#include <opm/simulators/wells/WellInterface.hpp>
#include <opm/simulators/wells/WellProdIndexCalculator.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>

#include <opm/models/blackoil/blackoilpolymermodules.hh>
#include <opm/models/blackoil/blackoilsolventmodules.hh>
#include <opm/models/blackoil/blackoilextbomodules.hh>
#include <opm/models/blackoil/blackoilfoammodules.hh>
#include <opm/models/blackoil/blackoilbrinemodules.hh>
#include <opm/models/blackoil/blackoilbioeffectsmodules.hh>

#include <opm/material/densead/Evaluation.hpp>
#include <opm/input/eclipse/Schedule/ScheduleTypes.hpp>

#include <opm/simulators/wells/StandardWellEval.hpp>

#include <dune/common/dynvector.hh>
#include <dune/common/dynmatrix.hh>

#include <memory>
#include <optional>

namespace Opm
{

    template<typename TypeTag>
    class StandardWell : public WellInterface<TypeTag>
                       , public StandardWellEval<GetPropType<TypeTag, Properties::FluidSystem>,
                                                 GetPropType<TypeTag, Properties::Indices>>
    {

    public:
        using Base = WellInterface<TypeTag>;
        using StdWellEval = StandardWellEval<GetPropType<TypeTag, Properties::FluidSystem>,
                                             GetPropType<TypeTag, Properties::Indices>>;

        // TODO: some functions working with AD variables handles only with values (double) without
        // dealing with derivatives. It can be beneficial to make functions can work with either AD or scalar value.
        // And also, it can also be beneficial to make these functions hanle different types of AD variables.
        using typename Base::Simulator;
        using typename Base::IntensiveQuantities;
        using typename Base::FluidSystem;
        using typename Base::MaterialLaw;
        using typename Base::ModelParameters;
        using typename Base::Indices;
        using typename Base::RateConverterType;
        using typename Base::SparseMatrixAdapter;
        using typename Base::FluidState;
        using typename Base::RateVector;

        using Base::has_solvent;
        using Base::has_zFraction;
        using Base::has_polymer;
        using Base::has_polymermw;
        using Base::has_foam;
        using Base::has_brine;
        using Base::has_energy;
        using Base::has_bioeffects;
        using Base::has_micp;

        using PolymerModule =  BlackOilPolymerModule<TypeTag>;
        using FoamModule = BlackOilFoamModule<TypeTag>;
        using typename Base::PressureMatrix;

        // number of the conservation equations
        static constexpr int numWellConservationEq = Indices::numPhases + Indices::numSolvents;
        // number of the well control equations
        static constexpr int numWellControlEq = 1;
        // number of the well equations that will always be used
        // based on the solution strategy, there might be other well equations be introduced
        static constexpr int numStaticWellEq = numWellConservationEq + numWellControlEq;

        // the index for Bhp in primary variables and also the index of well control equation
        // they both will be the last one in their respective system.
        // TODO: we should have indices for the well equations and well primary variables separately
        static constexpr int Bhp = numStaticWellEq - numWellControlEq;

        using StdWellEval::WQTotal;

        using typename Base::Scalar;

        using Base::name;
        using Base::Water;
        using Base::Oil;
        using Base::Gas;

        using typename Base::BVector;

        using Eval = typename StdWellEval::Eval;
        using EvalWell = typename StdWellEval::EvalWell;
        using BVectorWell = typename StdWellEval::BVectorWell;

        using IndexTraits = typename FluidSystem::IndexTraitsType;
        using WellStateType = WellState<Scalar, IndexTraits>;

        StandardWell(const Well& well,
                     const ParallelWellInfo<Scalar>& pw_info,
                     const int time_step,
                     const ModelParameters& param,
                     const RateConverterType& rate_converter,
                     const int pvtRegionIdx,
                     const int num_conservation_quantities,
                     const int num_phases,
                     const int index_of_well,
                     const std::vector<PerforationData<Scalar>>& perf_data);

        virtual void init(const std::vector<Scalar>& depth_arg,
                          const Scalar gravity_arg,
                          const std::vector<Scalar>& B_avg,
                          const bool changed_to_open_this_step) override;

        /// check whether the well equations get converged for this well
        virtual ConvergenceReport getWellConvergence(const Simulator& simulator,
                                                     const WellStateType& well_state,
                                                     const std::vector<Scalar>& B_avg,
                                                     DeferredLogger& deferred_logger,
                                                     const bool relax_tolerance) const override;

        /// Ax = Ax - C D^-1 B x
        virtual void apply(const BVector& x, BVector& Ax) const override;
        /// r = r - C D^-1 Rw
        virtual void apply(BVector& r) const override;

        /// using the solution x to recover the solution xw for wells and applying
        /// xw to update Well State
        void recoverWellSolutionAndUpdateWellState(const Simulator& simulator,
                                                   const BVector& x,
                                                   WellStateType& well_state,
                                                   DeferredLogger& deferred_logger) override;

        /// computing the well potentials for group control
        void computeWellPotentials(const Simulator& simulator,
                                   const WellStateType& well_state,
                                   std::vector<Scalar>& well_potentials,
                                   DeferredLogger& deferred_logger) /* const */ override;

        void updatePrimaryVariables(const Simulator& simulator,
                                    const WellStateType& well_state,
                                    DeferredLogger& deferred_logger) override;

        void solveEqAndUpdateWellState(const Simulator& simulator,
                                       WellStateType& well_state,
                                       DeferredLogger& deferred_logger) override;

        void calculateExplicitQuantities(const Simulator& simulator,
                                         const WellStateType& well_state,
                                         DeferredLogger& deferred_logger) override; // should be const?

        void updateProductivityIndex(const Simulator& simulator,
                                     const WellProdIndexCalculator<Scalar>& wellPICalc,
                                     WellStateType& well_state,
                                     DeferredLogger& deferred_logger) const override;

        Scalar connectionDensity(const int globalConnIdx,
                                 const int openConnIdx) const override;

        void addWellContributions(SparseMatrixAdapter& mat) const override;

        void addWellPressureEquations(PressureMatrix& mat,
                                      const BVector& x,
                                      const int pressureVarIndex,
                                      const bool use_well_weights,
                                      const WellStateType& well_state) const override;

        // iterate well equations with the specified control until converged
        bool iterateWellEqWithControl(const Simulator& simulator,
                                      const double dt,
                                      const Well::InjectionControls& inj_controls,
                                      const Well::ProductionControls& prod_controls,
                                      WellStateType& well_state,
                                      const GroupState<Scalar>& group_state,
                                      DeferredLogger& deferred_logger) override;

        // iterate well equations including control switching
        bool iterateWellEqWithSwitching(const Simulator& simulator,
                                        const double dt,
                                        const Well::InjectionControls& inj_controls,
                                        const Well::ProductionControls& prod_controls,
                                        WellStateType& well_state,
                                        const GroupState<Scalar>& group_state,
                                        DeferredLogger& deferred_logger,
                                        const bool fixed_control = false,
                                        const bool fixed_status = false) override;

        /* returns BHP */
        Scalar computeWellRatesAndBhpWithThpAlqProd(const Simulator& ebos_simulator,
                                                    const SummaryState &summary_state,
                                                    DeferredLogger& deferred_logger,
                                                    std::vector<Scalar>& potentials,
                                                    Scalar alq) const;

        void computeWellRatesWithThpAlqProd(const Simulator& ebos_simulator,
                                            const SummaryState& summary_state,
                                            DeferredLogger& deferred_logger,
                                            std::vector<Scalar>& potentials,
                                            Scalar alq) const;

        std::optional<Scalar>
        computeBhpAtThpLimitProdWithAlq(const Simulator& ebos_simulator,
                                        const SummaryState& summary_state,
                                        const Scalar alq_value,
                                        DeferredLogger& deferred_logger,
                                        bool iterate_if_no_solution) const override;

        void updateIPRImplicit(const Simulator& simulator,
                               WellStateType& well_state,
                               DeferredLogger& deferred_logger) override;

        void computeWellRatesWithBhp(const Simulator& ebosSimulator,
                                     const Scalar& bhp,
                                     std::vector<Scalar>& well_flux,
                                     DeferredLogger& deferred_logger) const override;

        // NOTE: These cannot be protected since they are used by GasLiftRuntime
        using Base::vfp_properties_;

        std::vector<Scalar>
        computeCurrentWellRates(const Simulator& ebosSimulator,
                                DeferredLogger& deferred_logger) const override;

        std::vector<Scalar> getPrimaryVars() const override;

        int setPrimaryVars(typename std::vector<Scalar>::const_iterator it) override;

    protected:
        bool regularize_;

        // updating the well_state based on well solution dwells
        void updateWellState(const Simulator& simulator,
                             const BVectorWell& dwells,
                             WellStateType& well_state,
                             DeferredLogger& deferred_logger);

        using WellConnectionProps = typename StdWellEval::StdWellConnections::Properties;

        // Compute connection level PVT properties needed to calulate the
        // pressure difference between well connections.
        WellConnectionProps
        computePropertiesForWellConnectionPressures(const Simulator& simulator,
                                                    const WellStateType& well_state) const;

        void computeWellConnectionDensitesPressures(const Simulator& simulator,
                                                    const WellStateType& well_state,
                                                    const WellConnectionProps& props,
                                                    DeferredLogger& deferred_logger);

        void computeWellConnectionPressures(const Simulator& simulator,
                                            const WellStateType& well_state,
                                            DeferredLogger& deferred_logger);

        template<class Value>
        void computePerfRate(const IntensiveQuantities& intQuants,
                             const std::vector<Value>& mob,
                             const Value& bhp,
                             const std::vector<Scalar>& Tw,
                             const int perf,
                             const bool allow_cf,
                             std::vector<Value>& cq_s,
                             PerforationRates<Scalar>& perf_rates,
                             DeferredLogger& deferred_logger) const;

        template<class Value>
        void computePerfRate(const std::vector<Value>& mob,
                             const Value& pressure,
                             const Value& bhp,
                             const Value& rs,
                             const Value& rv,
                             const Value& rvw,
                             const Value& rsw,
                             std::vector<Value>& b_perfcells_dense,
                             const std::vector<Scalar>& Tw,
                             const int perf,
                             const bool allow_cf,
                             const Value& skin_pressure,
                             const std::vector<Value>& cmix_s,
                             std::vector<Value>& cq_s,
                             PerforationRates<Scalar>& perf_rates,
                             DeferredLogger& deferred_logger) const;

        void computeWellRatesWithBhpIterations(const Simulator& ebosSimulator,
                                               const Scalar& bhp,
                                               std::vector<Scalar>& well_flux,
                                               DeferredLogger& deferred_logger) const override;

        std::vector<Scalar>
        computeWellPotentialWithTHP(const Simulator& ebosSimulator,
                                    DeferredLogger& deferred_logger,
                                    const WellStateType& well_state) const;

        bool computeWellPotentialsImplicit(const Simulator& ebos_simulator,
                                           const WellStateType& well_state,
                                           std::vector<Scalar>& well_potentials,
                                           DeferredLogger& deferred_logger) const;

        // return the density at the perforation[0] of the rank owning this well,
        // value is cached to minimize the number of broadcasts
        Scalar getRefDensity() const override;

        // get the mobility for specific perforation
        template<class Value>
        void getMobility(const Simulator& simulator,
                         const int perf,
                         std::vector<Value>& mob,
                         DeferredLogger& deferred_logger) const;

        void updateWaterMobilityWithPolymer(const Simulator& simulator,
                                            const int perf,
                                            std::vector<EvalWell>& mob_water,
                                            DeferredLogger& deferred_logger) const;

        void updatePrimaryVariablesNewton(const BVectorWell& dwells,
                                          const bool stop_or_zero_rate_target,
                                          DeferredLogger& deferred_logger);

        void updateWellStateFromPrimaryVariables(WellStateType& well_state,
                                                 const SummaryState& summary_state,
                                                 DeferredLogger& deferred_logger) const;

        void assembleWellEqWithoutIteration(const Simulator& simulator,
                                            const double dt,
                                            const Well::InjectionControls& inj_controls,
                                            const Well::ProductionControls& prod_controls,
                                            WellStateType& well_state,
                                            const GroupState<Scalar>& group_state,
                                            DeferredLogger& deferred_logger) override;

        void assembleWellEqWithoutIterationImpl(const Simulator& simulator,
                                                const double dt,
                                                const Well::InjectionControls& inj_controls,
                                                const Well::ProductionControls& prod_controls,
                                                WellStateType& well_state,
                                                const GroupState<Scalar>& group_state,
                                                DeferredLogger& deferred_logger);

        void calculateSinglePerf(const Simulator& simulator,
                                 const int perf,
                                 WellStateType& well_state,
                                 std::vector<RateVector>& connectionRates,
                                 std::vector<EvalWell>& cq_s,
                                 EvalWell& water_flux_s,
                                 EvalWell& cq_s_zfrac_effective,
                                 DeferredLogger& deferred_logger) const;

        // check whether the well is operable under BHP limit with current reservoir condition
        void checkOperabilityUnderBHPLimit(const WellStateType& well_state,
                                           const Simulator& simulator,
                                           DeferredLogger& deferred_logger) override;

        // check whether the well is operable under THP limit with current reservoir condition
        void checkOperabilityUnderTHPLimit(const Simulator& simulator,
                                           const WellStateType& well_state,
                                           DeferredLogger& deferred_logger) override;

        // updating the inflow based on the current reservoir condition
        void updateIPR(const Simulator& simulator,
                       DeferredLogger& deferred_logger) const override;

        // for a well, when all drawdown are in the wrong direction, then this well will not
        // be able to produce/inject .
        bool allDrawDownWrongDirection(const Simulator& simulator) const;

        // turn on crossflow to avoid singular well equations
        // when the well is banned from cross-flow and the BHP is not properly initialized,
        // we turn on crossflow to avoid singular well equations. It can result in wrong-signed
        // well rates, it can cause problem for THP calculation
        // TODO: looking for better alternative to avoid wrong-signed well rates
        bool openCrossFlowAvoidSingularity(const Simulator& simulator) const;

        // calculate the skin pressure based on water velocity, throughput and polymer concentration.
        // throughput is used to describe the formation damage during water/polymer injection.
        // calculated skin pressure will be applied to the drawdown during perforation rate calculation
        // to handle the effect from formation damage.
        EvalWell pskin(const Scalar throughput,
                       const EvalWell& water_velocity,
                       const EvalWell& poly_inj_conc,
                       DeferredLogger& deferred_logger) const;

        // calculate the skin pressure based on water velocity, throughput during water injection.
        EvalWell pskinwater(const Scalar throughput,
                            const EvalWell& water_velocity,
                            DeferredLogger& deferred_logger) const;

        // calculate the injecting polymer molecular weight based on the througput and water velocity
        EvalWell wpolymermw(const Scalar throughput,
                            const EvalWell& water_velocity,
                            DeferredLogger& deferred_logger) const;

        // modify the water rate for polymer injectivity study
        void handleInjectivityRate(const Simulator& simulator,
                                   const int perf,
                                   std::vector<EvalWell>& cq_s) const;

        // handle the extra equations for polymer injectivity study
        void handleInjectivityEquations(const Simulator& simulator,
                                        const WellStateType& well_state,
                                        const int perf,
                                        const EvalWell& water_flux_s,
                                        DeferredLogger& deferred_logger);

        void updateWaterThroughput(const double dt,
                                   WellStateType& well_state) const override;

        // checking convergence of extra equations, if there are any
        void checkConvergenceExtraEqs(const std::vector<Scalar>& res,
                                      ConvergenceReport& report) const;

        // updating the connectionRates_ related polymer molecular weight
        void updateConnectionRatePolyMW(const EvalWell& cq_s_poly,
                                        const IntensiveQuantities& int_quants,
                                        const WellStateType& well_state,
                                        const int perf,
                                        std::vector<RateVector>& connectionRates,
                                        DeferredLogger& deferred_logger) const;

        std::optional<Scalar>
        computeBhpAtThpLimitProd(const WellStateType& well_state,
                                 const Simulator& simulator,
                                 const SummaryState& summary_state,
                                 DeferredLogger& deferred_logger) const;

        std::optional<Scalar>
        computeBhpAtThpLimitInj(const Simulator& simulator,
                                const SummaryState& summary_state,
                                DeferredLogger& deferred_logger) const;

    private:
        Eval connectionRateEnergy(const std::vector<EvalWell>& cq_s,
                                  const IntensiveQuantities& intQuants,
                                  DeferredLogger& deferred_logger) const;

        // density of the first perforation, might not be from this rank
        Scalar cachedRefDensity{0};
    };

}

#include "StandardWell_impl.hpp"

#endif // OPM_STANDARDWELL_HEADER_INCLUDED
