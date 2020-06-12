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

#if HAVE_CUDA || HAVE_OPENCL
#include <opm/simulators/linalg/bda/WellContributions.hpp>
#endif

#include <opm/simulators/wells/RateConverter.hpp>
#include <opm/simulators/wells/WellInterface.hpp>

#include <opm/models/blackoil/blackoilpolymermodules.hh>
#include <opm/models/blackoil/blackoilsolventmodules.hh>
#include <opm/models/blackoil/blackoilfoammodules.hh>
#include <opm/models/blackoil/blackoilbrinemodules.hh>

#include <opm/material/densead/DynamicEvaluation.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/ScheduleTypes.hpp>

#include <dune/common/dynvector.hh>
#include <dune/common/dynmatrix.hh>

#include <optional>

namespace Opm
{

    template<typename TypeTag>
    class StandardWell: public WellInterface<TypeTag>
    {

    public:
        typedef WellInterface<TypeTag> Base;

        // TODO: some functions working with AD variables handles only with values (double) without
        // dealing with derivatives. It can be beneficial to make functions can work with either AD or scalar value.
        // And also, it can also be beneficial to make these functions hanle different types of AD variables.
        using typename Base::Simulator;
        using typename Base::WellState;
        using typename Base::IntensiveQuantities;
        using typename Base::FluidSystem;
        using typename Base::MaterialLaw;
        using typename Base::ModelParameters;
        using typename Base::Indices;
        using typename Base::RateConverterType;
        using typename Base::SparseMatrixAdapter;
        using typename Base::FluidState;

        using Base::numEq;

        using Base::has_solvent;
        using Base::has_polymer;
        using Base::has_foam;
        using Base::has_brine;
        using Base::has_energy;

        using PolymerModule =  Opm::BlackOilPolymerModule<TypeTag>;
        using FoamModule = Opm::BlackOilFoamModule<TypeTag>;
        using BrineModule = Opm::BlackOilBrineModule<TypeTag>;

        // polymer concentration and temperature are already known by the well, so
        // polymer and energy conservation do not need to be considered explicitly
        static const int numPolymerEq = Indices::numPolymers;
        static const int numEnergyEq = Indices::numEnergy;
        static const int numFoamEq = Indices::numFoam;
        static const int numBrineEq = Indices::numBrine;

        // number of the conservation equations
        static const int numWellConservationEq = numEq - numPolymerEq - numEnergyEq - numFoamEq - numBrineEq;
        // number of the well control equations
        static const int numWellControlEq = 1;
        // number of the well equations that will always be used
        // based on the solution strategy, there might be other well equations be introduced
        static const int numStaticWellEq = numWellConservationEq + numWellControlEq;

        // the positions of the primary variables for StandardWell
        // the first one is the weighted total rate (WQ_t), the second and the third ones are F_w and F_g,
        // which represent the fraction of Water and Gas based on the weighted total rate, the last one is BHP.
        // correspondingly, we have four well equations for blackoil model, the first three are mass
        // converstation equations, and the last one is the well control equation.
        // primary variables related to other components, will be before the Bhp and after F_g.
        // well control equation is always the last well equation.
        // TODO: in the current implementation, we use the well rate as the first primary variables for injectors,
        // instead of G_t.
        static const bool gasoil = numEq == 2 && (Indices::compositionSwitchIdx >= 0);
        static const int WQTotal = 0;
        static const int WFrac = gasoil? -1000: 1;
        static const int GFrac = gasoil? 1: 2;
        static const int SFrac = !has_solvent ? -1000 : 3;
        // the index for Bhp in primary variables and also the index of well control equation
        // they both will be the last one in their respective system.
        // TODO: we should have indices for the well equations and well primary variables separately
        static const int Bhp = numStaticWellEq - numWellControlEq;

        using typename Base::Scalar;


        using Base::name;
        using Base::Water;
        using Base::Oil;
        using Base::Gas;

        using typename Base::BVector;
        using typename Base::Eval;

        // sparsity pattern for the matrices
        //[A C^T    [x       =  [ res
        // B  D ]   x_well]      res_well]

        // the vector type for the res_well and x_well
        typedef Dune::DynamicVector<Scalar> VectorBlockWellType;
        typedef Dune::BlockVector<VectorBlockWellType> BVectorWell;

        // the matrix type for the diagonal matrix D
        typedef Dune::DynamicMatrix<Scalar> DiagMatrixBlockWellType;
        typedef Dune::BCRSMatrix <DiagMatrixBlockWellType> DiagMatWell;

        // the matrix type for the non-diagonal matrix B and C^T
        typedef Dune::DynamicMatrix<Scalar> OffDiagMatrixBlockWellType;
        typedef Dune::BCRSMatrix<OffDiagMatrixBlockWellType> OffDiagMatWell;

        typedef DenseAd::DynamicEvaluation<Scalar, numStaticWellEq + numEq + 1> EvalWell;

        using Base::contiSolventEqIdx;
        using Base::contiPolymerEqIdx;
        using Base::contiFoamEqIdx;
        using Base::contiBrineEqIdx;
        static const int contiEnergyEqIdx = Indices::contiEnergyEqIdx;

        StandardWell(const Well& well, const int time_step,
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
                          const int num_cells) override;


        virtual void initPrimaryVariablesEvaluation() const override;

        virtual void assembleWellEq(const Simulator& ebosSimulator,
                                    const std::vector<Scalar>& B_avg,
                                    const double dt,
                                    WellState& well_state,
                                    Opm::DeferredLogger& deferred_logger) override;

        virtual void updateWellStateWithTarget(const Simulator& ebos_simulator,
                                               WellState& well_state,
                                               Opm::DeferredLogger& deferred_logger) const override;

        /// check whether the well equations get converged for this well
        virtual ConvergenceReport getWellConvergence(const WellState& well_state,
                                                     const std::vector<double>& B_avg,
                                                     Opm::DeferredLogger& deferred_logger,
                                                     const bool relax_tolerance = false) const override;

        /// Ax = Ax - C D^-1 B x
        virtual void apply(const BVector& x, BVector& Ax) const override;
        /// r = r - C D^-1 Rw
        virtual void apply(BVector& r) const override;

#if HAVE_CUDA || HAVE_OPENCL
        /// add the contribution (C, D^-1, B matrices) of this Well to the WellContributions object
        void addWellContribution(WellContributions& wellContribs) const;

        /// get the number of blocks of the C and B matrices, used to allocate memory in a WellContributions object
        void getNumBlocks(unsigned int& _nnzs) const;
#endif

        /// using the solution x to recover the solution xw for wells and applying
        /// xw to update Well State
        virtual void recoverWellSolutionAndUpdateWellState(const BVector& x,
                                                           WellState& well_state,
                                                           Opm::DeferredLogger& deferred_logger) const override;

        /// computing the well potentials for group control
        virtual void computeWellPotentials(const Simulator& ebosSimulator,
                                           const std::vector<Scalar>& B_avg,
                                           const WellState& well_state,
                                           std::vector<double>& well_potentials,
                                           Opm::DeferredLogger& deferred_logger) /* const */ override;

        virtual void updatePrimaryVariables(const WellState& well_state, Opm::DeferredLogger& deferred_logger) const override;

        virtual void solveEqAndUpdateWellState(WellState& well_state, Opm::DeferredLogger& deferred_logger) override;

        virtual void calculateExplicitQuantities(const Simulator& ebosSimulator,
                                                 const WellState& well_state,
                                                 Opm::DeferredLogger& deferred_logger) override; // should be const?

        virtual void  addWellContributions(SparseMatrixAdapter& mat) const override;

        // iterate well equations with the specified control until converged
        bool iterateWellEqWithControl(const Simulator& ebosSimulator,
                                      const std::vector<double>& B_avg,
                                      const double dt,
                                      const Well::InjectionControls& inj_controls,
                                      const Well::ProductionControls& prod_controls,
                                      WellState& well_state,
                                      Opm::DeferredLogger& deferred_logger) override;

        /// \brief Wether the Jacobian will also have well contributions in it.
        virtual bool jacobianContainsWellContributions() const override
        {
            return param_.matrix_add_well_contributions_;
        }

    protected:

        // protected functions from the Base class
        using Base::getAllowCrossFlow;
        using Base::phaseUsage;
        using Base::flowPhaseToEbosCompIdx;
        using Base::ebosCompIdxToFlowCompIdx;
        using Base::wsalt;
        using Base::wsolvent;
        using Base::wpolymer;
        using Base::wfoam;
        using Base::scalingFactor;
        using Base::scaleProductivityIndex;
        using Base::mostStrictBhpFromBhpLimits;

        // protected member variables from the Base class
        using Base::current_step_;
        using Base::well_ecl_;
        using Base::vfp_properties_;
        using Base::gravity_;
        using Base::param_;
        using Base::well_efficiency_factor_;
        using Base::first_perf_;
        using Base::ref_depth_;
        using Base::perf_depth_;
        using Base::well_cells_;
        using Base::number_of_perforations_;
        using Base::number_of_phases_;
        using Base::saturation_table_number_;
        using Base::well_index_;
        using Base::index_of_well_;
        using Base::num_components_;
        using Base::connectionRates_;

        using Base::perf_rep_radius_;
        using Base::perf_length_;
        using Base::bore_diameters_;

        using Base::wellIsStopped_;

        // total number of the well equations and primary variables
        // there might be extra equations be used, numWellEq will be updated during the initialization
        int numWellEq_ = numStaticWellEq;

        // densities of the fluid in each perforation
        std::vector<double> perf_densities_;
        // pressure drop between different perforations
        std::vector<double> perf_pressure_diffs_;

        // residuals of the well equations
        BVectorWell resWell_;

        // two off-diagonal matrices
        OffDiagMatWell duneB_;
        OffDiagMatWell duneC_;
        // diagonal matrix for the well
        DiagMatWell invDuneD_;

        // several vector used in the matrix calculation
        mutable BVectorWell Bx_;
        mutable BVectorWell invDrw_;

        // the values for the primary varibles
        // based on different solutioin strategies, the wells can have different primary variables
        mutable std::vector<double> primary_variables_;

        // the Evaluation for the well primary variables, which contain derivativles and are used in AD calculation
        mutable std::vector<EvalWell> primary_variables_evaluation_;

        // the saturations in the well bore under surface conditions at the beginning of the time step
        std::vector<double> F0_;

        // the vectors used to describe the inflow performance relationship (IPR)
        // Q = IPR_A - BHP * IPR_B
        // TODO: it minght need to go to WellInterface, let us implement it in StandardWell first
        // it is only updated and used for producers for now
        mutable std::vector<double> ipr_a_;
        mutable std::vector<double> ipr_b_;

        bool changed_to_stopped_this_step_ = false;

        const EvalWell& getBhp() const;

        EvalWell getQs(const int comp_idx) const;

        const EvalWell& getWQTotal() const;

        EvalWell wellVolumeFractionScaled(const int phase) const;

        EvalWell wellVolumeFraction(const unsigned compIdx) const;

        EvalWell wellSurfaceVolumeFraction(const int phase) const;

        EvalWell extendEval(const Eval& in) const;

        Eval getPerfCellPressure(const FluidState& fs) const;

        // xw = inv(D)*(rw - C*x)
        void recoverSolutionWell(const BVector& x, BVectorWell& xw) const;

        // updating the well_state based on well solution dwells
        void updateWellState(const BVectorWell& dwells,
                             WellState& well_state,
                             Opm::DeferredLogger& deferred_logger) const;

        // calculate the properties for the well connections
        // to calulate the pressure difference between well connections.
        void computePropertiesForWellConnectionPressures(const Simulator& ebosSimulator,
                                                         const WellState& well_state,
                                                         std::vector<double>& b_perf,
                                                         std::vector<double>& rsmax_perf,
                                                         std::vector<double>& rvmax_perf,
                                                         std::vector<double>& surf_dens_perf) const;

        // TODO: not total sure whether it is a good idea to put this function here
        // the major reason to put here is to avoid the usage of Wells struct
        void computeConnectionDensities(const std::vector<double>& perfComponentRates,
                                        const std::vector<double>& b_perf,
                                        const std::vector<double>& rsmax_perf,
                                        const std::vector<double>& rvmax_perf,
                                        const std::vector<double>& surf_dens_perf);

        void computeConnectionPressureDelta();

        void computeWellConnectionDensitesPressures(const WellState& well_state,
                                                    const std::vector<double>& b_perf,
                                                    const std::vector<double>& rsmax_perf,
                                                    const std::vector<double>& rvmax_perf,
                                                    const std::vector<double>& surf_dens_perf);

        // computing the accumulation term for later use in well mass equations
        void computeAccumWell();

        void computeWellConnectionPressures(const Simulator& ebosSimulator,
                                                    const WellState& well_state);

        struct PerfRates {
            // surface perforation rates
            std::vector<EvalWell> cq_s;
            // reservoir perforation rates
            std::vector<EvalWell> cq_r;
            // reservor perfration total rates
            EvalWell cq_r_t;

            double perf_dis_gas_rate;

            double perf_vap_oil_rate;

            PerfRates(int num_components, int num_welleq)
            : cq_s(num_components, {num_welleq + numEq, 0.})
            , cq_r(num_components, {num_welleq + numEq, 0.})
            , cq_r_t(num_welleq + numEq, 0.)
            , perf_dis_gas_rate(0.)
            , perf_vap_oil_rate(0.)
            {
            }
        };

        void computePerfRate(const IntensiveQuantities& intQuants,
                             const std::vector<EvalWell>& mob,
                             const EvalWell& bhp,
                             const double Tw,
                             const int perf,
                             const bool allow_cf,
                             PerfRates& perf_rates,
                             Opm::DeferredLogger& deferred_logger) const;

        void computeWellRatesWithBhp(const Simulator& ebosSimulator,
                                             const double& bhp,
                                             std::vector<double>& well_flux,
                                             Opm::DeferredLogger& deferred_logger) const;

        void computeWellRatesWithBhpPotential(const Simulator& ebosSimulator,
                                              const std::vector<Scalar>& B_avg,
                                              const double& bhp,
                                              std::vector<double>& well_flux,
                                              Opm::DeferredLogger& deferred_logger);

        std::vector<double> computeWellPotentialWithTHP(const Simulator& ebosSimulator,
                                                        Opm::DeferredLogger& deferred_logger) const;

        template <class ValueType>
        ValueType calculateBhpFromThp(const std::vector<ValueType>& rates, const Well& well, const SummaryState& summaryState, Opm::DeferredLogger& deferred_logger) const;


        double calculateThpFromBhp(const std::vector<double>& rates, const double bhp, Opm::DeferredLogger& deferred_logger) const;

        // get the mobility for specific perforation
        void getMobility(const Simulator& ebosSimulator,
                         const int perf,
                         std::vector<EvalWell>& mob,
                         Opm::DeferredLogger& deferred_logger) const;

        void updateWaterMobilityWithPolymer(const Simulator& ebos_simulator,
                                            const int perf,
                                            std::vector<EvalWell>& mob_water,
                                            Opm::DeferredLogger& deferred_logger) const;

        void updatePrimaryVariablesNewton(const BVectorWell& dwells,
                                          const WellState& well_state) const;

        // update extra primary vriables if there are any
        void updateExtraPrimaryVariables(const BVectorWell& dwells) const;


        void updateWellStateFromPrimaryVariables(WellState& well_state, Opm::DeferredLogger& deferred_logger) const;

        void updateThp(WellState& well_state, Opm::DeferredLogger& deferred_logger) const;

        void assembleControlEq(const WellState& well_state,
                               const Opm::Schedule& schedule,
                               const SummaryState& summaryState,
                               Opm::DeferredLogger& deferred_logger);

        // handle the non reasonable fractions due to numerical overshoot
        void processFractions() const;

        // updating the inflow based on the current reservoir condition
        void updateIPR(const Simulator& ebos_simulator, Opm::DeferredLogger& deferred_logger) const;

        // update the operability status of the well is operable under the current reservoir condition
        // mostly related to BHP limit and THP limit
        virtual void checkWellOperability(const Simulator& ebos_simulator,
                                          const WellState& well_state,
                                          Opm::DeferredLogger& deferred_logger) override;

        virtual void assembleWellEqWithoutIteration(const Simulator& ebosSimulator,
                                                    const double dt,
                                                    const Well::InjectionControls& inj_controls,
                                                    const Well::ProductionControls& prod_controls,
                                                    WellState& well_state,
                                                    Opm::DeferredLogger& deferred_logger) override;

        // check whether the well is operable under the current reservoir condition
        // mostly related to BHP limit and THP limit
        void updateWellOperability(const Simulator& ebos_simulator,
                                   const WellState& well_state,
                                   Opm::DeferredLogger& deferred_logger);

        // check whether the well is operable under BHP limit with current reservoir condition
        void checkOperabilityUnderBHPLimitProducer(const Simulator& ebos_simulator, Opm::DeferredLogger& deferred_logger);

        // check whether the well is operable under THP limit with current reservoir condition
        void checkOperabilityUnderTHPLimitProducer(const Simulator& ebos_simulator, Opm::DeferredLogger& deferred_logger);

        // for a well, when all drawdown are in the wrong direction, then this well will not
        // be able to produce/inject .
        bool allDrawDownWrongDirection(const Simulator& ebos_simulator) const;

        // whether the well can produce / inject based on the current well state (bhp)
        bool canProduceInjectWithCurrentBhp(const Simulator& ebos_simulator,
                                            const WellState& well_state,
                                            Opm::DeferredLogger& deferred_logger);

        // turn on crossflow to avoid singular well equations
        // when the well is banned from cross-flow and the BHP is not properly initialized,
        // we turn on crossflow to avoid singular well equations. It can result in wrong-signed
        // well rates, it can cause problem for THP calculation
        // TODO: looking for better alternative to avoid wrong-signed well rates
        bool openCrossFlowAvoidSingularity(const Simulator& ebos_simulator) const;

        // relaxation factor considering only one fraction value
        static double relaxationFactorFraction(const double old_value,
                                               const double dx);

        // calculate a relaxation factor to avoid overshoot of the fractions for producers
        // which might result in negative rates
        static double relaxationFactorFractionsProducer(const std::vector<double>& primary_variables,
                                                        const BVectorWell& dwells);

        // calculate a relaxation factor to avoid overshoot of total rates
        static double relaxationFactorRate(const std::vector<double>& primary_variables,
                                           const BVectorWell& dwells);

        virtual void wellTestingPhysical(const Simulator& simulator, const std::vector<double>& B_avg,
                                         const double simulation_time, const int report_step,
                                         WellState& well_state, WellTestState& welltest_state,
                                         Opm::DeferredLogger& deferred_logger) override;

        // calculate the skin pressure based on water velocity, throughput and polymer concentration.
        // throughput is used to describe the formation damage during water/polymer injection.
        // calculated skin pressure will be applied to the drawdown during perforation rate calculation
        // to handle the effect from formation damage.
        EvalWell pskin(const double throuhgput,
                       const EvalWell& water_velocity,
                       const EvalWell& poly_inj_conc,
                       Opm::DeferredLogger& deferred_logger) const;

        // calculate the skin pressure based on water velocity, throughput during water injection.
        EvalWell pskinwater(const double throughput,
                            const EvalWell& water_velocity,
                            Opm::DeferredLogger& deferred_logger) const;

        // calculate the injecting polymer molecular weight based on the througput and water velocity
        EvalWell wpolymermw(const double throughput,
                            const EvalWell& water_velocity,
                            Opm::DeferredLogger& deferred_logger) const;

        // handle the extra equations for polymer injectivity study
        void handleInjectivityRateAndEquations(const IntensiveQuantities& int_quants,
                                               const WellState& well_state,
                                               const int perf,
                                               std::vector<EvalWell>& cq_s,
                                               Opm::DeferredLogger& deferred_logger);

        virtual void updateWaterThroughput(const double dt, WellState& well_state) const override;

        // checking the convergence of the well control equations
        void checkConvergenceControlEq(const WellState& well_state,
                                       ConvergenceReport& report,
                                       DeferredLogger& deferred_logger) const;

        // checking convergence of extra equations, if there are any
        void checkConvergenceExtraEqs(const std::vector<double>& res,
                                      ConvergenceReport& report) const;

        // updating the connectionRates_ related polymer molecular weight
        void updateConnectionRatePolyMW(const EvalWell& cq_s_poly,
                                        const IntensiveQuantities& int_quants,
                                        const WellState& well_state,
                                        const int perf,
                                        DeferredLogger& deferred_logger);


        std::optional<double> computeBhpAtThpLimitProd(const Simulator& ebos_simulator,
                                                       const SummaryState& summary_state,
                                                       DeferredLogger& deferred_logger) const;

        std::optional<double> computeBhpAtThpLimitInj(const Simulator& ebos_simulator,
                                                      const SummaryState& summary_state,
                                                      DeferredLogger& deferred_logger) const;

        // assembling the control equation for sequential transport solution
        void assembleControlEqSeqTrans(const WellState& well_state);

    };

}

#include "StandardWell_impl.hpp"

#endif // OPM_STANDARDWELL_HEADER_INCLUDED
