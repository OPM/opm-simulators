/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2017 IRIS
  Copyright 2019 Norce

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


#ifndef OPM_WELLINTERFACE_HEADER_INCLUDED
#define OPM_WELLINTERFACE_HEADER_INCLUDED

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>

#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/WellTestState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group/GuideRate.hpp>

#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/simulators/timestepping/SimulatorReport.hpp>

#include <opm/simulators/wells/RateConverter.hpp>
#include <opm/simulators/wells/VFPProperties.hpp>
#include <opm/simulators/wells/WellHelpers.hpp>
#include <opm/simulators/wells/WellGroupHelpers.hpp>
#include <opm/simulators/wells/WellStateFullyImplicitBlackoil.hpp>
#include <opm/simulators/flow/BlackoilModelParametersEbos.hpp>

#include <opm/simulators/timestepping/ConvergenceReport.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>

#include<dune/common/fmatrix.hh>
#include<dune/istl/bcrsmatrix.hh>
#include<dune/istl/matrixmatrix.hh>

#include <opm/material/densead/Math.hpp>
#include <opm/material/densead/Evaluation.hpp>

#include <string>
#include <memory>
#include <vector>
#include <cassert>

namespace Opm
{


    template<typename TypeTag>
    class WellInterface
    {
    public:

        using WellState = WellStateFullyImplicitBlackoil;

        typedef BlackoilModelParametersEbos<TypeTag> ModelParameters;

        static const int Water = BlackoilPhases::Aqua;
        static const int Oil = BlackoilPhases::Liquid;
        static const int Gas = BlackoilPhases::Vapour;

        typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
        typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
        typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
        typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
        typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
        typedef typename GET_PROP_TYPE(TypeTag, SparseMatrixAdapter) SparseMatrixAdapter;
        typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;

        static const int numEq = Indices::numEq;
        typedef double Scalar;

        typedef Dune::FieldVector<Scalar, numEq    > VectorBlockType;
        typedef Dune::FieldMatrix<Scalar, numEq, numEq > MatrixBlockType;
        typedef Dune::BlockVector<VectorBlockType> BVector;
        typedef DenseAd::Evaluation<double, /*size=*/numEq> Eval;

        static const bool has_solvent = GET_PROP_VALUE(TypeTag, EnableSolvent);
        static const bool has_polymer = GET_PROP_VALUE(TypeTag, EnablePolymer);
        static const bool has_energy = GET_PROP_VALUE(TypeTag, EnableEnergy);
        static const bool has_temperature = GET_PROP_VALUE(TypeTag, EnableTemperature);
        // flag for polymer molecular weight related
        static const bool has_polymermw = GET_PROP_VALUE(TypeTag, EnablePolymerMW);
        static const bool has_foam = GET_PROP_VALUE(TypeTag, EnableFoam);
        static const bool has_brine = GET_PROP_VALUE(TypeTag, EnableBrine);
        static const int contiSolventEqIdx = Indices::contiSolventEqIdx;
        static const int contiPolymerEqIdx = Indices::contiPolymerEqIdx;
        // index for the polymer molecular weight continuity equation
        static const int contiPolymerMWEqIdx = Indices::contiPolymerMWEqIdx;
        static const int contiFoamEqIdx = Indices::contiFoamEqIdx;
        static const int contiBrineEqIdx = Indices::contiBrineEqIdx;

        // For the conversion between the surface volume rate and reservoir voidage rate
        using RateConverterType = RateConverter::
        SurfaceToReservoirVoidage<FluidSystem, std::vector<int> >;
        static const bool compositionSwitchEnabled = Indices::gasEnabled;
        using FluidState = Opm::BlackOilFluidState<Eval,
                                                   FluidSystem,
                                                   has_temperature,
                                                   has_energy,
                                                   compositionSwitchEnabled,
                                                   has_brine,
                                                   Indices::numPhases >;
        /// Constructor
        WellInterface(const Well& well, const int time_step,
                      const ModelParameters& param,
                      const RateConverterType& rate_converter,
                      const int pvtRegionIdx,
                      const int num_components,
                      const int num_phases,
                      const int index_of_well,
                      const int first_perf_index,
                      const std::vector<PerforationData>& perf_data);

        /// Virutal destructor
        virtual ~WellInterface() {}

        /// Well name.
        const std::string& name() const;

        /// True if the well is an injector.
        bool isInjector() const;

        /// True if the well is a producer.
        bool isProducer() const;

        /// Index of well in the wells struct and wellState
        int indexOfWell() const;

        /// Well cells.
        const std::vector<int>& cells() const {return well_cells_; }

        void setVFPProperties(const VFPProperties<VFPInjProperties,VFPProdProperties>* vfp_properties_arg);

        void setGuideRate(const GuideRate* guide_rate_arg);

        virtual void init(const PhaseUsage* phase_usage_arg,
                          const std::vector<double>& depth_arg,
                          const double gravity_arg,
                          const int num_cells);

        virtual void initPrimaryVariablesEvaluation() const = 0;

        virtual ConvergenceReport getWellConvergence(const WellState& well_state, const std::vector<double>& B_avg, Opm::DeferredLogger& deferred_logger, const bool relax_tolerance = false) const = 0;

        virtual void solveEqAndUpdateWellState(WellState& well_state, Opm::DeferredLogger& deferred_logger) = 0;

        virtual void assembleWellEq(const Simulator& ebosSimulator,
                                    const std::vector<Scalar>& B_avg,
                                    const double dt,
                                    WellState& well_state,
                                    Opm::DeferredLogger& deferred_logger
                                    ) = 0;

        void updateWellTestState(const WellState& well_state,
                                 const double& simulationTime,
                                 const bool& writeMessageToOPMLog,
                                 WellTestState& wellTestState,
                                 Opm::DeferredLogger& deferred_logger) const;

        void setWellEfficiencyFactor(const double efficiency_factor);

        void computeRepRadiusPerfLength(const Grid& grid, const std::vector<int>& cartesian_to_compressed, Opm::DeferredLogger& deferred_logger);

        /// using the solution x to recover the solution xw for wells and applying
        /// xw to update Well State
        virtual void recoverWellSolutionAndUpdateWellState(const BVector& x,
                                                           WellState& well_state,
                                                           Opm::DeferredLogger& deferred_logger) const = 0;

        /// Ax = Ax - C D^-1 B x
        virtual void apply(const BVector& x, BVector& Ax) const = 0;

        /// r = r - C D^-1 Rw
        virtual void apply(BVector& r) const = 0;

        // TODO: before we decide to put more information under mutable, this function is not const
        virtual void computeWellPotentials(const Simulator& ebosSimulator,
                                           const std::vector<Scalar>& B_avg,
                                           const WellState& well_state,
                                           std::vector<double>& well_potentials,
                                           Opm::DeferredLogger& deferred_logger) = 0;

        virtual void updateWellStateWithTarget(const Simulator& ebos_simulator,
                                               WellState& well_state,
                                               Opm::DeferredLogger& deferred_logger) const = 0;

        enum class IndividualOrGroup { Individual, Group, Both };
        bool updateWellControl(const Simulator& ebos_simulator,
                               const IndividualOrGroup iog,
                               WellState& well_state,
                               Opm::DeferredLogger& deferred_logger) /* const */;

        virtual void updatePrimaryVariables(const WellState& well_state, Opm::DeferredLogger& deferred_logger) const = 0;

        virtual void calculateExplicitQuantities(const Simulator& ebosSimulator,
                                                 const WellState& well_state,
                                                 Opm::DeferredLogger& deferred_logger) = 0; // should be const?

        /// \brief Wether the Jacobian will also have well contributions in it.
        virtual bool jacobianContainsWellContributions() const
        {
            return false;
        }

        // updating the voidage rates in well_state when requested
        void calculateReservoirRates(WellState& well_state) const;

        // Add well contributions to matrix
        virtual void addWellContributions(SparseMatrixAdapter&) const = 0;

        void addCellRates(RateVector& rates, int cellIdx) const;

        Scalar volumetricSurfaceRateForConnection(int cellIdx, int phaseIdx) const;


        template <class EvalWell>
        Eval restrictEval(const EvalWell& in) const
        {
            Eval out = 0.0;
            out.setValue(in.value());
            for(int eqIdx = 0; eqIdx < numEq;++eqIdx) {
                out.setDerivative(eqIdx, in.derivative(eqIdx));
            }
            return out;
        }

        void closeCompletions(WellTestState& wellTestState);

        const Well& wellEcl() const;

        // TODO: theoretically, it should be a const function
        // Simulator is not const is because that assembleWellEq is non-const Simulator
        void wellTesting(const Simulator& simulator, const std::vector<double>& B_avg,
                         const double simulation_time, const int report_step,
                         const WellTestConfig::Reason testing_reason,
                         /* const */ WellState& well_state, WellTestState& welltest_state,
                         Opm::DeferredLogger& deferred_logger);

        void updatePerforatedCell(std::vector<bool>& is_cell_perforated);

        virtual void checkWellOperability(const Simulator& ebos_simulator, const WellState& well_state, Opm::DeferredLogger& deferred_logger) = 0;

        // whether the well is operable
        bool isOperable() const;

        /// Returns true if the well has one or more THP limits/constraints.
        bool wellHasTHPConstraints(const SummaryState& summaryState) const;

        /// Returns true if the well is currently in prediction mode (i.e. not history mode).
        bool underPredictionMode() const;

        // update perforation water throughput based on solved water rate
        virtual void updateWaterThroughput(const double dt, WellState& well_state) const = 0;

        void stopWell() {
            wellIsStopped_ = true;
        }
        void openWell() {
            wellIsStopped_ = false;
        }

        bool wellIsStopped() const {
            return wellIsStopped_;
        }

        void setWsolvent(const double wsolvent);


    protected:

        // to indicate a invalid completion
        static const int INVALIDCOMPLETION = INT_MAX;

        Well well_ecl_;

        const int current_step_;

        // simulation parameters
        const ModelParameters& param_;

        // number of the perforations for this well
        int number_of_perforations_;

        // well index for each perforation
        std::vector<double> well_index_;

        // depth for each perforation
        std::vector<double> perf_depth_;

        // reference depth for the BHP
        double ref_depth_;

        double well_efficiency_factor_;

        // cell index for each well perforation
        std::vector<int> well_cells_;

        // saturation table nubmer for each well perforation
        std::vector<int> saturation_table_number_;

        // representative radius of the perforations, used in shear calculation
        std::vector<double> perf_rep_radius_;

        // length of the perforations, use in shear calculation
        std::vector<double> perf_length_;

        // well bore diameter
        std::vector<double> bore_diameters_;

        /*
         *  completions_ contains the mapping from completion id to connection indices
         *  {
         *      2 : [ConnectionIndex, ConnectionIndex],
         *      1 : [ConnectionIndex, ConnectionIndex, ConnectionIndex],
         *      5 : [ConnectionIndex],
         *      7 : [ConnectionIndex]
         *      ...
         *   }
         *   The integer IDs correspond to the COMPLETION id given by the COMPLUMP keyword.
         *   When there is no COMPLUMP keyword used, a default completion number will be assigned
         *   based on the order of the declaration of the connections.
         *   Since the connections not OPEN is not included in the Wells, so they will not be considered
         *   in this mapping relation.
         */
        std::map<int, std::vector<int>> completions_;

        const PhaseUsage* phase_usage_;

        bool getAllowCrossFlow() const;

        const VFPProperties<VFPInjProperties,VFPProdProperties>* vfp_properties_;

        const GuideRate* guide_rate_;

        double gravity_;

        // For the conversion between the surface volume rate and resrevoir voidage rate
        const RateConverterType& rateConverter_;

        // The pvt region of the well. We assume
        // We assume a well to not penetrate more than one pvt region.
        const int pvtRegionIdx_;

        const int num_components_;

        // number of phases
        int number_of_phases_;

        // the index of well in Wells struct
        int index_of_well_;

        // record the index of the first perforation
        // of states of individual well.
        int first_perf_;

        std::vector<int> originalConnectionIndex_;

        std::vector<RateVector> connectionRates_;

        // rates under reservoir condtion for each phase/component?
        // RateVector is based on Dune::FieldVector<Evaluation, numEq>
        std::vector<RateVector> conResRates_;

        std::vector<Eval> conTotalResRates_;

        bool wellIsStopped_;

        double wsolvent_;

        const PhaseUsage& phaseUsage() const;

        int flowPhaseToEbosCompIdx( const int phaseIdx ) const;

        int ebosCompIdxToFlowCompIdx( const unsigned compIdx ) const;

        double wsolvent() const;

        double wpolymer() const;

        double wfoam() const;

        double wsalt() const;

        bool checkRateEconLimits(const WellEconProductionLimits& econ_production_limits,
                                 const WellState& well_state,
                                 Opm::DeferredLogger& deferred_logger) const;

        double getTHPConstraint(const SummaryState& summaryState) const;

        // Component fractions for each phase for the well
        const std::vector<double>& compFrac() const;

        double mostStrictBhpFromBhpLimits(const SummaryState& summaryState) const;

        struct RatioLimitCheckReport;

        void checkMaxWaterCutLimit(const WellEconProductionLimits& econ_production_limits,
                              const WellState& well_state,
                              RatioLimitCheckReport& report) const;

        void checkMaxGORLimit(const WellEconProductionLimits& econ_production_limits,
                              const WellState& well_state,
                              RatioLimitCheckReport& report) const;

        void checkMaxWGRLimit(const WellEconProductionLimits& econ_production_limits,
                              const WellState& well_state,
                              RatioLimitCheckReport& report) const;

        void checkRatioEconLimits(const WellEconProductionLimits& econ_production_limits,
                                  const WellState& well_state,
                                  RatioLimitCheckReport& report,
                                  Opm::DeferredLogger& deferred_logger) const;


        template <typename RatioFunc>
        bool checkMaxRatioLimitWell(const WellState& well_state,
                                    const double max_ratio_limit,
                                    const RatioFunc& ratioFunc) const;

        template <typename RatioFunc>
        void checkMaxRatioLimitCompletions(const WellState& well_state,
                                           const double max_ratio_limit,
                                           const RatioFunc& ratioFunc,
                                           RatioLimitCheckReport& report) const;

        double scalingFactor(const int comp_idx) const;

        // whether a well is specified with a non-zero and valid VFP table number
        bool isVFPActive(Opm::DeferredLogger& deferred_logger) const;

        struct OperabilityStatus;

        OperabilityStatus operability_status_;

        void wellTestingEconomic(const Simulator& simulator, const std::vector<double>& B_avg,
                                 const double simulation_time, const WellState& well_state,
                                 WellTestState& welltest_state, Opm::DeferredLogger& deferred_logger);

        virtual void wellTestingPhysical(const Simulator& simulator, const std::vector<double>& B_avg,
                                 const double simulation_time, const int report_step,
                                         WellState& well_state, WellTestState& welltest_state, Opm::DeferredLogger& deferred_logger) = 0;


        virtual void assembleWellEqWithoutIteration(const Simulator& ebosSimulator,
                                                    const double dt,
                                                    const Well::InjectionControls& inj_controls,
                                                    const Well::ProductionControls& prod_controls,
                                                    WellState& well_state,
                                                    Opm::DeferredLogger& deferred_logger) = 0;

        // iterate well equations with the specified control until converged
        virtual bool iterateWellEqWithControl(const Simulator& ebosSimulator,
                                              const std::vector<double>& B_avg,
                                              const double dt,
                                              const Well::InjectionControls& inj_controls,
                                              const Well::ProductionControls& prod_controls,
                                              WellState& well_state,
                                              Opm::DeferredLogger& deferred_logger) = 0;

        bool iterateWellEquations(const Simulator& ebosSimulator,
                                  const std::vector<double>& B_avg,
                                  const double dt,
                                  WellState& well_state,
                                  Opm::DeferredLogger& deferred_logger);

        void updateWellTestStateEconomic(const WellState& well_state,
                                         const double simulation_time,
                                         const bool write_message_to_opmlog,
                                         WellTestState& well_test_state,
                                         Opm::DeferredLogger& deferred_logger) const;

        void updateWellTestStatePhysical(const WellState& well_state,
                                         const double simulation_time,
                                         const bool write_message_to_opmlog,
                                         WellTestState& well_test_state,
                                         Opm::DeferredLogger& deferred_logger) const;

        void solveWellForTesting(const Simulator& ebosSimulator, WellState& well_state,
                                 const std::vector<double>& B_avg,
                                 Opm::DeferredLogger& deferred_logger);

        void scaleProductivityIndex(const int perfIdx, double& productivity_index, const bool new_well, Opm::DeferredLogger& deferred_logger);

        void initCompletions();

        // count the number of times an output log message is created in the productivity
        // index calculations
        int well_productivity_index_logger_counter_;

        bool checkConstraints(WellState& well_state,
                              const Schedule& schedule,
                              const SummaryState& summaryState,
                              DeferredLogger& deferred_logger) const;

        bool checkIndividualConstraints(WellState& well_state,
                                        const SummaryState& summaryState) const;

        bool checkGroupConstraints(WellState& well_state,
                                   const Schedule& schedule,
                                   const SummaryState& summaryState,
                                   DeferredLogger& deferred_logger) const;

        std::pair<bool, double> checkGroupConstraintsProd(const Group& group,
                                       const WellState& well_state,
                                       const double efficiencyFactor,
                                       const Schedule& schedule,
                                       const SummaryState& summaryState,
                                       DeferredLogger& deferred_logger) const;

        std::pair<bool, double> checkGroupConstraintsInj(const Group& group,
                                      const WellState& well_state,
                                      const double efficiencyFactor,
                                      const Schedule& schedule,
                                      const SummaryState& summaryState,
                                      DeferredLogger& deferred_logger) const;

        template <class EvalWell>
        void getGroupInjectionControl(const Group& group,
                                      const WellState& well_state,
                                      const Opm::Schedule& schedule,
                                      const SummaryState& summaryState,
                                      const InjectorType& injectorType,
                                      const EvalWell& bhp,
                                      const EvalWell& injection_rate,
                                      EvalWell& control_eq,
                                      double efficiencyFactor);

        template <class EvalWell>
        void getGroupProductionControl(const Group& group,
                                       const WellState& well_state,
                                       const Opm::Schedule& schedule,
                                       const SummaryState& summaryState,
                                       const EvalWell& bhp,
                                       const std::vector<EvalWell>& rates,
                                       EvalWell& control_eq,
                                       double efficiencyFactor);

        template <class EvalWell, class BhpFromThpFunc>
        void assembleControlEqInj(const WellState& well_state,
                                  const Opm::Schedule& schedule,
                                  const SummaryState& summaryState,
                                  const Well::InjectionControls& controls,
                                  const EvalWell& bhp,
                                  const EvalWell& injection_rate,
                                  BhpFromThpFunc bhp_from_thp,
                                  EvalWell& control_eq,
                                  Opm::DeferredLogger& deferred_logger);

        template <class EvalWell, class BhpFromThpFunc>
        void assembleControlEqProd(const WellState& well_state,
                                   const Opm::Schedule& schedule,
                                   const SummaryState& summaryState,
                                   const Well::ProductionControls& controls,
                                   const EvalWell& bhp,
                                   const std::vector<EvalWell>& rates, // Always 3 canonical rates.
                                   BhpFromThpFunc bhp_from_thp,
                                   EvalWell& control_eq,
                                   Opm::DeferredLogger& deferred_logger);
    };





    // definition of the struct OperabilityStatus
    template<typename TypeTag>
    struct
    WellInterface<TypeTag>::
    OperabilityStatus {
        bool isOperable() const {
            if (!operable_under_only_bhp_limit) {
                return false;
            } else {
                return ( (isOperableUnderBHPLimit() || isOperableUnderTHPLimit()) );
            }
        }

        bool isOperableUnderBHPLimit() const {
            return operable_under_only_bhp_limit && obey_thp_limit_under_bhp_limit;
        }

        bool isOperableUnderTHPLimit() const {
            return can_obtain_bhp_with_thp_limit && obey_bhp_limit_with_thp_limit;
        }

        void reset() {
            operable_under_only_bhp_limit = true;
            obey_thp_limit_under_bhp_limit = true;
            can_obtain_bhp_with_thp_limit = true;
            obey_bhp_limit_with_thp_limit = true;
        }

        // whether the well can be operated under bhp limit
        // without considering other limits.
        // if it is false, then the well is not operable for sure.
        bool operable_under_only_bhp_limit = true;
        // if the well can be operated under bhp limit, will it obey(not violate)
        // the thp limit when operated under bhp limit
        bool obey_thp_limit_under_bhp_limit = true;
        // whether the well operate under the thp limit only
        bool can_obtain_bhp_with_thp_limit = true;
        // whether the well obey bhp limit when operated under thp limit
        bool obey_bhp_limit_with_thp_limit = true;

    };


    template<typename TypeTag>
    struct
    WellInterface<TypeTag>::
    RatioLimitCheckReport{
        bool ratio_limit_violated = false;
        int worst_offending_completion = INVALIDCOMPLETION;
        double violation_extent = 0.0;
    };

    const std::string modestring[4] = { "BHP", "THP", "RESERVOIR_RATE", "SURFACE_RATE" };

}

#include "WellInterface_impl.hpp"

#endif // OPM_WELLINTERFACE_HEADER_INCLUDED
