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
#include <opm/models/blackoil/blackoilmicpmodules.hh>

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
                                                 GetPropType<TypeTag, Properties::Indices>,
                                                 GetPropType<TypeTag, Properties::Scalar>>
    {

    public:
        using Base = WellInterface<TypeTag>;
        using StdWellEval = StandardWellEval<GetPropType<TypeTag, Properties::FluidSystem>,
                                             GetPropType<TypeTag, Properties::Indices>,
                                             GetPropType<TypeTag, Properties::Scalar>>;

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
        using Base::has_micp;

        using PolymerModule =  BlackOilPolymerModule<TypeTag>;
        using FoamModule = BlackOilFoamModule<TypeTag>;
        using BrineModule = BlackOilBrineModule<TypeTag>;
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

        StandardWell(const Well& well,
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
                                           DeferredLogger& deferred_logger) /* const */ override;

        void updatePrimaryVariables(const SummaryState& summary_state,
                                    const WellState& well_state,
                                    DeferredLogger& deferred_logger) override;

        virtual void solveEqAndUpdateWellState(const SummaryState& summary_state,
                                               WellState& well_state,
                                               DeferredLogger& deferred_logger) override;

        virtual void calculateExplicitQuantities(const Simulator& ebosSimulator,
                                                 const WellState& well_state,
                                                 DeferredLogger& deferred_logger) override; // should be const?

        virtual void updateProductivityIndex(const Simulator& ebosSimulator,
                                             const WellProdIndexCalculator& wellPICalc,
                                             WellState& well_state,
                                             DeferredLogger& deferred_logger) const override;

        virtual double connectionDensity(const int globalConnIdx,
                                         const int openConnIdx) const override;

        virtual void addWellContributions(SparseMatrixAdapter& mat) const override;

        virtual void addWellPressureEquations(PressureMatrix& mat,
                                              const BVector& x,
                                              const int pressureVarIndex,
                                              const bool use_well_weights,
                                              const WellState& well_state) const override;

        // iterate well equations with the specified control until converged
        bool iterateWellEqWithControl(const Simulator& ebosSimulator,
                                      const double dt,
                                      const Well::InjectionControls& inj_controls,
                                      const Well::ProductionControls& prod_controls,
                                      WellState& well_state,
                                      const GroupState& group_state,
                                      DeferredLogger& deferred_logger) override;

        // iterate well equations including control switching
        bool iterateWellEqWithSwitching(const Simulator& ebosSimulator,
                                        const double dt,
                                        const Well::InjectionControls& inj_controls,
                                        const Well::ProductionControls& prod_controls,
                                        WellState& well_state,
                                        const GroupState& group_state,
                                        DeferredLogger& deferred_logger, 
                                        const bool allow_switch = true) override;

        /// \brief Wether the Jacobian will also have well contributions in it.
        virtual bool jacobianContainsWellContributions() const override
        {
            return this->param_.matrix_add_well_contributions_;
        }

        /* returns BHP */
        double computeWellRatesAndBhpWithThpAlqProd(const Simulator &ebos_simulator,
                               const SummaryState &summary_state,
                               DeferredLogger &deferred_logger,
                               std::vector<double> &potentials,
                               double alq) const;

        void computeWellRatesWithThpAlqProd(
            const Simulator &ebos_simulator,
            const SummaryState &summary_state,
            DeferredLogger &deferred_logger,
            std::vector<double> &potentials,
            double alq) const;

        std::optional<double> computeBhpAtThpLimitProdWithAlq(
            const Simulator& ebos_simulator,
            const SummaryState& summary_state,
            const double alq_value,
            DeferredLogger& deferred_logger) const override;

        void updateIPRImplicit(const Simulator& ebosSimulator, WellState& well_state, DeferredLogger& deferred_logger);                

        virtual void computeWellRatesWithBhp(
            const Simulator& ebosSimulator,
            const double& bhp,
            std::vector<double>& well_flux,
            DeferredLogger& deferred_logger) const override;

        // NOTE: These cannot be protected since they are used by GasLiftRuntime
        using Base::phaseUsage;
        using Base::vfp_properties_;

        virtual std::vector<double> computeCurrentWellRates(const Simulator& ebosSimulator,
                                                            DeferredLogger& deferred_logger) const override;

        std::vector<double> getPrimaryVars() const override;

        int setPrimaryVars(std::vector<double>::const_iterator it) override;

    protected:
        bool regularize_;

        // updating the well_state based on well solution dwells
        void updateWellState(const SummaryState& summary_state,
                             const BVectorWell& dwells,
                             WellState& well_state,
                             DeferredLogger& deferred_logger);

        // calculate the properties for the well connections
        // to calulate the pressure difference between well connections.
        using WellConnectionProps = typename StdWellEval::StdWellConnections::Properties;
        void computePropertiesForWellConnectionPressures(const Simulator& ebosSimulator,
                                                         const WellState& well_state,
                                                         WellConnectionProps& props) const;

        void computeWellConnectionDensitesPressures(const Simulator& ebosSimulator,
                                                    const WellState& well_state,
                                                    const WellConnectionProps& props,
                                                    DeferredLogger& deferred_logger);

        void computeWellConnectionPressures(const Simulator& ebosSimulator,
                                            const WellState& well_state,
                                            DeferredLogger& deferred_logger);

        template<class Value>
        void computePerfRate(const IntensiveQuantities& intQuants,
                             const std::vector<Value>& mob,
                             const Value& bhp,
                             const std::vector<Scalar>& Tw,
                             const int perf,
                             const bool allow_cf,
                             std::vector<Value>& cq_s,
                             PerforationRates& perf_rates,
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
                             PerforationRates& perf_rates,
                             DeferredLogger& deferred_logger) const;

        void computeWellRatesWithBhpIterations(const Simulator& ebosSimulator,
                                              const double& bhp,
                                              std::vector<double>& well_flux,
                                              DeferredLogger& deferred_logger) const override;

        std::vector<double> computeWellPotentialWithTHP(
            const Simulator& ebosSimulator,
            DeferredLogger& deferred_logger,
            const WellState &well_state) const;

        bool computeWellPotentialsImplicit(const Simulator& ebos_simulator,
                                           std::vector<double>& well_potentials,
                                           DeferredLogger& deferred_logger) const;               

        virtual double getRefDensity() const override;

        // get the mobility for specific perforation
        template<class Value>
        void getMobility(const Simulator& ebosSimulator,
                         const int perf,
                         std::vector<Value>& mob,
                         DeferredLogger& deferred_logger) const;

        void updateWaterMobilityWithPolymer(const Simulator& ebos_simulator,
                                            const int perf,
                                            std::vector<EvalWell>& mob_water,
                                            DeferredLogger& deferred_logger) const;

        void updatePrimaryVariablesNewton(const BVectorWell& dwells,
                                          const bool stop_or_zero_rate_target,
                                          DeferredLogger& deferred_logger);

        void updateWellStateFromPrimaryVariables(const bool stop_or_zero_rate_target,
                                                 WellState& well_state,
                                                 const SummaryState& summary_state,
                                                 DeferredLogger& deferred_logger) const;

        virtual void assembleWellEqWithoutIteration(const Simulator& ebosSimulator,
                                                    const double dt,
                                                    const Well::InjectionControls& inj_controls,
                                                    const Well::ProductionControls& prod_controls,
                                                    WellState& well_state,
                                                    const GroupState& group_state,
                                                    DeferredLogger& deferred_logger) override;

        void assembleWellEqWithoutIterationImpl(const Simulator& ebosSimulator,
                                                const double dt,
                                                const Well::InjectionControls& inj_controls,
                                                const Well::ProductionControls& prod_controls,
                                                WellState& well_state,
                                                const GroupState& group_state,
                                                DeferredLogger& deferred_logger);

        void calculateSinglePerf(const Simulator& ebosSimulator,
                                 const int perf,
                                 WellState& well_state,
                                 std::vector<RateVector>& connectionRates,
                                 std::vector<EvalWell>& cq_s,
                                 EvalWell& water_flux_s,
                                 EvalWell& cq_s_zfrac_effective,
                                 DeferredLogger& deferred_logger) const;

        // check whether the well is operable under BHP limit with current reservoir condition
        virtual void checkOperabilityUnderBHPLimit(const WellState& well_state, const Simulator& ebos_simulator, DeferredLogger& deferred_logger) override;

        // check whether the well is operable under THP limit with current reservoir condition
        virtual void checkOperabilityUnderTHPLimit(const Simulator& ebos_simulator, const WellState& well_state, DeferredLogger& deferred_logger) override;

        // updating the inflow based on the current reservoir condition
        virtual void updateIPR(const Simulator& ebos_simulator, DeferredLogger& deferred_logger) const override;

        // for a well, when all drawdown are in the wrong direction, then this well will not
        // be able to produce/inject .
        bool allDrawDownWrongDirection(const Simulator& ebos_simulator) const;

        // whether the well can produce / inject based on the current well state (bhp)
        bool canProduceInjectWithCurrentBhp(const Simulator& ebos_simulator,
                                            const WellState& well_state,
                                            DeferredLogger& deferred_logger);

        // turn on crossflow to avoid singular well equations
        // when the well is banned from cross-flow and the BHP is not properly initialized,
        // we turn on crossflow to avoid singular well equations. It can result in wrong-signed
        // well rates, it can cause problem for THP calculation
        // TODO: looking for better alternative to avoid wrong-signed well rates
        bool openCrossFlowAvoidSingularity(const Simulator& ebos_simulator) const;

        // calculate the skin pressure based on water velocity, throughput and polymer concentration.
        // throughput is used to describe the formation damage during water/polymer injection.
        // calculated skin pressure will be applied to the drawdown during perforation rate calculation
        // to handle the effect from formation damage.
        EvalWell pskin(const double throuhgput,
                       const EvalWell& water_velocity,
                       const EvalWell& poly_inj_conc,
                       DeferredLogger& deferred_logger) const;

        // calculate the skin pressure based on water velocity, throughput during water injection.
        EvalWell pskinwater(const double throughput,
                            const EvalWell& water_velocity,
                            DeferredLogger& deferred_logger) const;

        // calculate the injecting polymer molecular weight based on the througput and water velocity
        EvalWell wpolymermw(const double throughput,
                            const EvalWell& water_velocity,
                            DeferredLogger& deferred_logger) const;

        // modify the water rate for polymer injectivity study
        void handleInjectivityRate(const Simulator& ebosSimulator,
                                   const int perf,
                                   std::vector<EvalWell>& cq_s) const;

        // handle the extra equations for polymer injectivity study
        void handleInjectivityEquations(const Simulator& ebosSimulator,
                                        const WellState& well_state,
                                        const int perf,
                                        const EvalWell& water_flux_s,
                                        DeferredLogger& deferred_logger);

        virtual void updateWaterThroughput(const double dt, WellState& well_state) const override;

        // checking convergence of extra equations, if there are any
        void checkConvergenceExtraEqs(const std::vector<double>& res,
                                      ConvergenceReport& report) const;

        // updating the connectionRates_ related polymer molecular weight
        void updateConnectionRatePolyMW(const EvalWell& cq_s_poly,
                                        const IntensiveQuantities& int_quants,
                                        const WellState& well_state,
                                        const int perf,
                                        std::vector<RateVector>& connectionRates,
                                        DeferredLogger& deferred_logger) const;


        std::optional<double> computeBhpAtThpLimitProd(const WellState& well_state,
                                                       const Simulator& ebos_simulator,
                                                       const SummaryState& summary_state,
                                                       DeferredLogger& deferred_logger) const;

        std::optional<double> computeBhpAtThpLimitInj(const Simulator& ebos_simulator,
                                                      const SummaryState& summary_state,
                                                      DeferredLogger& deferred_logger) const;

    private:
        Eval connectionRateEnergy(const double maxOilSaturation,
                                  const std::vector<EvalWell>& cq_s,
                                  const IntensiveQuantities& intQuants,
                                  DeferredLogger& deferred_logger) const;

        template<class Value>
        void gasOilPerfRateInj(const std::vector<Value>& cq_s,
                               PerforationRates& perf_rates,
                               const Value& rv,
                               const Value& rs,
                               const Value& pressure,
                               const Value& rvw,
                               DeferredLogger& deferred_logger) const;

        template<class Value>
        void gasOilPerfRateProd(std::vector<Value>& cq_s,
                                PerforationRates& perf_rates,
                                const Value& rv,
                                const Value& rs,
                                const Value& rvw) const;

        template<class Value>
        void gasWaterPerfRateProd(std::vector<Value>& cq_s,
                                  PerforationRates& perf_rates,
                                  const Value& rvw,
                                  const Value& rsw) const;

        template<class Value>
        void gasWaterPerfRateInj(const std::vector<Value>& cq_s,
                                 PerforationRates& perf_rates,
                                 const Value& rvw,
                                 const Value& rsw,
                                 const Value& pressure,
                                 DeferredLogger& deferred_logger) const;

        template<class Value>
        void disOilVapWatVolumeRatio(Value& volumeRatio,
                                     const Value& rvw,
                                     const Value& rsw,
                                     const Value& pressure,
                                     const std::vector<Value>& cmix_s,
                                     const std::vector<Value>& b_perfcells_dense,
                                     DeferredLogger& deferred_logger) const;

        template<class Value>
        void gasOilVolumeRatio(Value& volumeRatio,
                               const Value& rv,
                               const Value& rs,
                               const Value& pressure,
                               const std::vector<Value>& cmix_s,
                               const std::vector<Value>& b_perfcells_dense,
                               DeferredLogger& deferred_logger) const;
    };

}

#include "StandardWell_impl.hpp"

#endif // OPM_STANDARDWELL_HEADER_INCLUDED
