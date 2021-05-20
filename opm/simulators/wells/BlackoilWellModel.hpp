/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 - 2017 Statoil ASA.
  Copyright 2017 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2016 - 2018 IRIS AS

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


#ifndef OPM_BLACKOILWELLMODEL_HEADER_INCLUDED
#define OPM_BLACKOILWELLMODEL_HEADER_INCLUDED

#include <ebos/eclproblem.hh>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <cassert>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <stddef.h>

#include <opm/parser/eclipse/EclipseState/Runspec.hpp>

#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/WellTestState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group/GuideRate.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group/Group.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group/GConSale.hpp>

#include <opm/simulators/timestepping/SimulatorReport.hpp>
#include <opm/simulators/flow/countGlobalCells.hpp>
#include <opm/simulators/wells/GasLiftSingleWell.hpp>
#include <opm/simulators/wells/GasLiftStage2.hpp>
#include <opm/simulators/wells/GasLiftWellState.hpp>
#include <opm/simulators/wells/PerforationData.hpp>
#include <opm/simulators/wells/VFPInjProperties.hpp>
#include <opm/simulators/wells/VFPProdProperties.hpp>
#include <opm/simulators/wells/WellState.hpp>
#include <opm/simulators/wells/WGState.hpp>
#include <opm/simulators/wells/RateConverter.hpp>
#include <opm/simulators/wells/WellInterface.hpp>
#include <opm/simulators/wells/StandardWell.hpp>
#include <opm/simulators/wells/MultisegmentWell.hpp>
#include <opm/simulators/wells/WellGroupHelpers.hpp>
#include <opm/simulators/wells/WellProdIndexCalculator.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>
#include <opm/simulators/timestepping/gatherConvergenceReport.hpp>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixmatrix.hh>

#include <opm/material/densead/Math.hpp>

#include <opm/simulators/utils/DeferredLogger.hpp>

namespace Opm::Properties {

template<class TypeTag, class MyTypeTag>
struct EnableTerminalOutput {
    using type = UndefinedProperty;
};

} // namespace Opm::Properties

namespace Opm {

        /// Class for handling the blackoil well model.
        template<typename TypeTag>
        class BlackoilWellModel : public BaseAuxiliaryModule<TypeTag>
        {
        public:
            // ---------      Types      ---------
            typedef BlackoilModelParametersEbos<TypeTag> ModelParameters;

            using Grid = GetPropType<TypeTag, Properties::Grid>;
            using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
            using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
            using Indices = GetPropType<TypeTag, Properties::Indices>;
            using Simulator = GetPropType<TypeTag, Properties::Simulator>;
            using Scalar = GetPropType<TypeTag, Properties::Scalar>;
            using RateVector = GetPropType<TypeTag, Properties::RateVector>;
            using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
            using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;

            typedef typename BaseAuxiliaryModule<TypeTag>::NeighborSet NeighborSet;
            using GasLiftSingleWell = ::Opm::GasLiftSingleWell<TypeTag>;
            using GasLiftStage2 = ::Opm::GasLiftStage2<TypeTag>;
            using GLiftWellState = ::Opm::GasLiftWellState<TypeTag>;
            using GLiftWellStateMap =
                std::map<std::string,std::unique_ptr<GLiftWellState>>;
            using GLiftOptWells =
                std::map<std::string,std::unique_ptr<GasLiftSingleWell>>;
            using GLiftProdWells =
                std::map<std::string,const WellInterface<TypeTag> *>;

            static const int numEq = Indices::numEq;
            static const int solventSaturationIdx = Indices::solventSaturationIdx;
            static constexpr bool has_solvent_ = getPropValue<TypeTag, Properties::EnableSolvent>();
            static constexpr bool has_polymer_ = getPropValue<TypeTag, Properties::EnablePolymer>();
            static constexpr bool has_energy_ = getPropValue<TypeTag, Properties::EnableEnergy>();

            // TODO: where we should put these types, WellInterface or Well Model?
            // or there is some other strategy, like TypeTag
            typedef Dune::FieldVector<Scalar, numEq    > VectorBlockType;
            typedef Dune::BlockVector<VectorBlockType> BVector;

            typedef Dune::FieldMatrix<Scalar, numEq, numEq > MatrixBlockType;

            typedef BlackOilPolymerModule<TypeTag> PolymerModule;

            // For the conversion between the surface volume rate and resrevoir voidage rate
            using RateConverterType = RateConverter::
                SurfaceToReservoirVoidage<FluidSystem, std::vector<int> >;

            BlackoilWellModel(Simulator& ebosSimulator);

            void init();

            /////////////
            // <eWoms auxiliary module stuff>
            /////////////
            unsigned numDofs() const
            // No extra dofs are inserted for wells. (we use a Schur complement.)
            { return 0; }

            void addNeighbors(std::vector<NeighborSet>& neighbors) const;

            void applyInitial()
            {}

            void linearize(SparseMatrixAdapter& jacobian, GlobalEqVector& res);

            void postSolve(GlobalEqVector& deltaX)
            {
                recoverWellSolutionAndUpdateWellState(deltaX);
            }

            /////////////
            // </ eWoms auxiliary module stuff>
            /////////////

            template <class Restarter>
            void deserialize(Restarter& /* res */)
            {
                // TODO (?)
            }

            /*!
             * \brief This method writes the complete state of the well
             *        to the harddisk.
             */
            template <class Restarter>
            void serialize(Restarter& /* res*/)
            {
                // TODO (?)
            }

            void beginEpisode()
            {
                beginReportStep(ebosSimulator_.episodeIndex());
            }

            void beginTimeStep();

            void beginIteration()
            {
                assemble(ebosSimulator_.model().newtonMethod().numIterations(),
                         ebosSimulator_.timeStepSize());
            }

            void endIteration()
            { }

            void endTimeStep()
            {
                timeStepSucceeded(ebosSimulator_.time(), ebosSimulator_.timeStepSize());
            }

            void endEpisode()
            {
                endReportStep();
            }

            template <class Context>
            void computeTotalRatesForDof(RateVector& rate,
                                         const Context& context,
                                         unsigned spaceIdx,
                                         unsigned timeIdx) const;


            using WellInterfacePtr = std::shared_ptr<WellInterface<TypeTag> >;
            WellInterfacePtr well(const std::string& wellName) const;

            void initFromRestartFile(const RestartValue& restartValues);

            data::GroupAndNetworkValues
            groupAndNetworkData(const int reportStepIdx, const Schedule& sched) const
            {
                auto grp_nwrk_values = ::Opm::data::GroupAndNetworkValues{};

                this->assignGroupValues(reportStepIdx, sched, grp_nwrk_values.groupData);
                this->assignNodeValues(grp_nwrk_values.nodeData);

                return grp_nwrk_values;
            }


            /*
              The dynamic state of the well model is maintained with an instance
              of the WellState class. Currently we have
              three different wellstate instances:

               1. The currently active wellstate is in the active_well_state_
                  member. That is the state which is mutated by the simulator.

               2. In the case timestep fails to converge and we must go back and
                  try again with a smaller timestep we need to recover the last
                  valid wellstate. This is maintained with the
                  last_valid_well_state_ member and the functions
                  commitWellState() and resetWellState().

                3. For the NUPCOL functionality we should either use the
                   currently active wellstate or a wellstate frozen at max
                   nupcol iterations. This is handled with the member
                   nupcol_well_state_ and the initNupcolWellState() function.
            */


            /*
              Immutable version of the currently active wellstate.
            */
            const WellState& wellState() const
            {
                return this->active_wgstate_.well_state;
            }

            /*
              Mutable version of the currently active wellstate.
            */
            WellState& wellState()
            {
                return this->active_wgstate_.well_state;
            }

            /*
              Will return the last good wellstate. This is typcially used when
              initializing a new report step where the Schedule object might
              have introduced new wells. The wellstate returned by
              prevWellState() must have been stored with the commitWellState()
              function first.
            */
            const WellState& prevWellState() const
            {
                return this->last_valid_wgstate_.well_state;
            }

            const WGState& prevWGState() const
            {
                return this->last_valid_wgstate_;
            }
            /*
              Will return the currently active nupcolWellState; must initialize
              the internal nupcol wellstate with initNupcolWellState() first.
            */
            const WellState& nupcolWellState() const
            {
                return this->nupcol_wgstate_.well_state;
            }

            /*
              Will assign the internal member last_valid_well_state_ to the
              current value of the this->active_well_state_. The state stored
              with storeWellState() can then subsequently be recovered with the
              resetWellState() method.
            */
            void commitWGState()
            {
                this->last_valid_wgstate_ = this->active_wgstate_;
            }

            /*
              Will store a copy of the input argument well_state in the
              last_valid_well_state_ member, that state can then be recovered
              with a subsequent call to resetWellState().
            */
            void commitWGState(WGState wgstate)
            {
                this->last_valid_wgstate_ = std::move(wgstate);
            }

            /*
              Will update the internal variable active_well_state_ to whatever
              was stored in the last_valid_well_state_ member. This function
              works in pair with commitWellState() which should be called first.
            */
            void resetWGState()
            {
                this->active_wgstate_ = this->last_valid_wgstate_;
            }

            /*
              Will store the current active wellstate in the nupcol_well_state_
              member. This can then be subsequently retrieved with accessor
              nupcolWellState().
            */
            void updateNupcolWGState()
            {
                this->nupcol_wgstate_ = this->active_wgstate_;
            }

            const GroupState& groupState() const
            {
                return this->active_wgstate_.group_state;
            }

            data::Wells wellData() const
            {
                auto wsrpt = this->wellState().report(UgGridHelpers::globalCell(grid()),
                                                      [this](const int well_ndex) -> bool
                                                      {
                                                          return this->wasDynamicallyShutThisTimeStep(well_ndex);
                                                      });

                this->assignWellGuideRates(wsrpt);
                this->assignShutConnections(wsrpt);

                return wsrpt;
            }

            // substract Binv(D)rw from r;
            void apply( BVector& r) const;

            // subtract B*inv(D)*C * x from A*x
            void apply(const BVector& x, BVector& Ax) const;

#if HAVE_CUDA || HAVE_OPENCL
            // accumulate the contributions of all Wells in the WellContributions object
            void getWellContributions(WellContributions& x) const;
#endif

            // apply well model with scaling of alpha
            void applyScaleAdd(const Scalar alpha, const BVector& x, BVector& Ax) const;

            // Check if well equations is converged.
            ConvergenceReport getWellConvergence(const std::vector<Scalar>& B_avg, const bool checkGroupConvergence = false) const;

            const PhaseUsage& phaseUsage() const { return phase_usage_; }

            const SimulatorReportSingle& lastReport() const;

            void addWellContributions(SparseMatrixAdapter& jacobian) const
            {
                for ( const auto& well: well_container_ ) {
                    well->addWellContributions(jacobian);
                }
            }

            // called at the beginning of a report step
            void beginReportStep(const int time_step);

            /// Return true if any well has a THP constraint.
            bool hasTHPConstraints() const;

            /// Shut down any single well, but only if it is in prediction mode.
            /// Returns true if the well was actually found and shut.
            bool forceShutWellByNameIfPredictionMode(const std::string& wellname, const double simulation_time);

            void updateEclWells(const int timeStepIdx, const std::unordered_set<std::string>& wells);
            bool hasWell(const std::string& wname);
            double wellPI(const int well_index) const;
            double wellPI(const std::string& well_name) const;

            void updatePerforationIntensiveQuantities();
            // it should be able to go to prepareTimeStep(), however, the updateWellControls() and initPrimaryVariablesEvaluation()
            // makes it a little more difficult. unless we introduce if (iterationIdx != 0) to avoid doing the above functions
            // twice at the beginning of the time step
            /// Calculating the explict quantities used in the well calculation. By explicit, we mean they are cacluated
            /// at the beginning of the time step and no derivatives are included in these quantities
            void calculateExplicitQuantities(DeferredLogger& deferred_logger) const;
            // some preparation work, mostly related to group control and RESV,
            // at the beginning of each time step (Not report step)
            void prepareTimeStep(DeferredLogger& deferred_logger);
            void initPrimaryVariablesEvaluation() const;
            void updateWellControls(DeferredLogger& deferred_logger, const bool checkGroupControls);
            WellInterfacePtr getWell(const std::string& well_name) const;
        protected:
            Simulator& ebosSimulator_;

            std::vector< Well > wells_ecl_{};
            std::vector< std::vector<PerforationData> > well_perf_data_{};
            std::vector< WellProdIndexCalculator > prod_index_calc_{};
            std::vector<int> local_shut_wells_{};

            std::vector< ParallelWellInfo > parallel_well_info_;
            std::vector< ParallelWellInfo* > local_parallel_well_info_;

            bool wells_active_{false};

            // a vector of all the wells.
            std::vector<WellInterfacePtr > well_container_{};

            // Map from logically cartesian cell indices to compressed ones.
            // Cells not in the interior are not mapped. This deactivates
            // these for distributed wells and makes the distribution non-overlapping.
            std::vector<int> cartesian_to_compressed_{};

            std::vector<bool> is_cell_perforated_{};

            std::function<bool(const Well&)> not_on_process_{};

            void initializeWellProdIndCalculators();
            void initializeWellPerfData();
            void initializeWellState(const int           timeStepIdx,
                                     const SummaryState& summaryState);

            // create the well container
            std::vector<WellInterfacePtr > createWellContainer(const int time_step);

            void inferLocalShutWells();

            WellInterfacePtr
            createWellPointer(const int wellID,
                              const int time_step) const;

            template <typename WellType>
            std::unique_ptr<WellType>
            createTypedWellPointer(const int wellID,
                                   const int time_step) const;

            WellInterfacePtr createWellForWellTest(const std::string& well_name, const int report_step, DeferredLogger& deferred_logger) const;


            const ModelParameters param_;
            bool terminal_output_{false};
            std::vector<int> pvt_region_idx_{};
            PhaseUsage phase_usage_;
            size_t global_num_cells_{};
            // the number of the cells in the local grid
            size_t local_num_cells_{};
            double gravity_{};
            std::vector<double> depth_{};
            bool initial_step_{};
            bool report_step_starts_{};
            bool glift_debug = false;
            bool alternative_well_rate_init_{};

            std::optional<int> last_run_wellpi_{};

            std::unique_ptr<RateConverterType> rateConverter_{};
            std::unique_ptr<VFPProperties> vfp_properties_{};

            SimulatorReportSingle last_report_{};

            WellTestState wellTestState_{};
            std::unique_ptr<GuideRate> guideRate_{};

            std::map<std::string, double> node_pressures_{}; // Storing network pressures for output.
            mutable std::unordered_set<std::string> closed_this_step_{};

            // used to better efficiency of calcuation
            mutable BVector scaleAddRes_{};

            std::vector<Scalar> B_avg_{};

            const Grid& grid() const
            { return ebosSimulator_.vanguard().grid(); }

            const EclipseState& eclState() const
            { return ebosSimulator_.vanguard().eclState(); }

            const Schedule& schedule() const
            { return ebosSimulator_.vanguard().schedule(); }

            void gliftDebug(
                const std::string &msg,
                DeferredLogger& deferred_logger) const;

            /// \brief Get the wells of our partition that are not shut.
            /// \param timeStepIdx The index of the time step.
            /// \param[out] globalNumWells the number of wells globally.
            std::vector< Well > getLocalWells(const int timeStepId) const;

            /// \brief Create the parallel well information
            /// \param localWells The local wells from ECL schedule
            std::vector< ParallelWellInfo* >
            createLocalParallelWellInfo(const std::vector<Well>& localWells);

            // compute the well fluxes and assemble them in to the reservoir equations as source terms
            // and in the well equations.
            void assemble(const int iterationIdx,
                          const double dt);

            // called at the end of a time step
            void timeStepSucceeded(const double& simulationTime, const double dt);

            // called at the end of a report step
            void endReportStep();

            // using the solution x to recover the solution xw for wells and applying
            // xw to update Well State
            void recoverWellSolutionAndUpdateWellState(const BVector& x);

            void updateAndCommunicateGroupData();
            void updateNetworkPressures();

            // setting the well_solutions_ based on well_state.
            void updatePrimaryVariables(DeferredLogger& deferred_logger);

            void setupCartesianToCompressed_(const int* global_cell, int local_num__cells);

            void setRepRadiusPerfLength();


            void updateAverageFormationFactor();

            // Calculating well potentials for each well
            void updateWellPotentials(const int reportStepIdx, const bool onlyAfterEvent, DeferredLogger& deferred_logger);

            const std::vector<double>& wellPerfEfficiencyFactors() const;

            void calculateEfficiencyFactors(const int reportStepIdx);

            void calculateProductivityIndexValuesShutWells(const int reportStepIdx, DeferredLogger& deferred_logger);
            void calculateProductivityIndexValues(DeferredLogger& deferred_logger);
            void calculateProductivityIndexValues(const WellInterface<TypeTag>* wellPtr,
                                                  DeferredLogger& deferred_logger);

            // The number of components in the model.
            int numComponents() const;

            int numLocalWells() const;

            int numPhases() const;

            int reportStepIndex() const;

            void assembleWellEq(const double dt, DeferredLogger& deferred_logger);

            void maybeDoGasLiftOptimize(DeferredLogger& deferred_logger);

            void gliftDebugShowALQ(DeferredLogger& deferred_logger);

            void gasLiftOptimizationStage2(DeferredLogger& deferred_logger,
                GLiftProdWells &prod_wells, GLiftOptWells &glift_wells,
                GLiftWellStateMap &map);

            void extractLegacyCellPvtRegionIndex_();

            void extractLegacyDepth_();

            /// return true if wells are available in the reservoir
            bool wellsActive() const;

            void setWellsActive(const bool wells_active);

            /// return true if wells are available on this process
            bool localWellsActive() const;

            /// upate the wellTestState related to economic limits
            void updateWellTestState(const double& simulationTime, WellTestState& wellTestState) const;

            void wellTesting(const int timeStepIdx, const double simulationTime, DeferredLogger& deferred_logger);

            // convert well data from opm-common to well state from opm-core
            void loadRestartData( const data::Wells& wells,
                               const data::GroupAndNetworkValues& grpNwrkValues,
                               const PhaseUsage& phases,
                               const bool handle_ms_well,
                               WellState& state );

            // whether there exists any multisegment well open on this process
            bool anyMSWellOpenLocal() const;

            const Well& getWellEcl(const std::string& well_name) const;

            void updateGroupIndividualControls(DeferredLogger& deferred_logger, std::set<std::string>& switched_groups);
            void updateGroupIndividualControl(const Group& group, DeferredLogger& deferred_logger, std::set<std::string>& switched_groups);
            bool checkGroupConstraints(const Group& group, DeferredLogger& deferred_logger) const;
            Group::ProductionCMode checkGroupProductionConstraints(const Group& group, DeferredLogger& deferred_logger) const;
            Group::InjectionCMode checkGroupInjectionConstraints(const Group& group, const Phase& phase) const;
            void checkGconsaleLimits(const Group& group, WellState& well_state, DeferredLogger& deferred_logger );

            void updateGroupHigherControls(DeferredLogger& deferred_logger, std::set<std::string>& switched_groups);
            void checkGroupHigherConstraints(const Group& group, DeferredLogger& deferred_logger, std::set<std::string>& switched_groups);

            void actionOnBrokenConstraints(const Group& group, const Group::ExceedAction& exceed_action, const Group::ProductionCMode& newControl, DeferredLogger& deferred_logger);

            void actionOnBrokenConstraints(const Group& group, const Group::InjectionCMode& newControl, const Phase& topUpPhase, DeferredLogger& deferred_logger);

            void updateWsolvent(const Group& group, const Schedule& schedule, const int reportStepIdx, const WellState& wellState);

            void setWsolvent(const Group& group, const Schedule& schedule, const int reportStepIdx, double wsolvent);

            void runWellPIScaling(const int timeStepIdx, DeferredLogger& local_deferredLogger);

            bool wasDynamicallyShutThisTimeStep(const int well_index) const;

            void assignWellGuideRates(data::Wells& wsrpt) const;
            void assignShutConnections(data::Wells& wsrpt) const;
            void assignGroupValues(const int                               reportStepIdx,
                                   const Schedule&                         sched,
                                   std::map<std::string, data::GroupData>& gvalues) const;

            void assignNodeValues(std::map<std::string, data::NodeData>& gvalues) const;

            std::unordered_map<std::string, data::GroupGuideRates>
            calculateAllGroupGuiderates(const int reportStepIdx, const Schedule& sched) const;

            void assignGroupControl(const Group& group, data::GroupData& gdata) const;
            data::GuideRateValue getGuideRateValues(const Well& well) const;
            data::GuideRateValue getGuideRateValues(const Group& group) const;
            data::GuideRateValue getGuideRateInjectionGroupValues(const Group& group) const;
            void getGuideRateValues(const GuideRate::RateVector& qs,
                                    const bool                   is_inj,
                                    const std::string&           wgname,
                                    data::GuideRateValue&        grval) const;

            void assignGroupGuideRates(const Group& group,
                                       const std::unordered_map<std::string, data::GroupGuideRates>& groupGuideRates,
                                       data::GroupData& gdata) const;

             void computeWellTemperature();                       
        private:
            GroupState& groupState() { return this->active_wgstate_.group_state; }
            BlackoilWellModel(Simulator& ebosSimulator, const PhaseUsage& pu);
            /*
              The various wellState members should be accessed and modified
              through the accessor functions wellState(), prevWellState(),
              commitWellState(), resetWellState(), nupcolWellState() and
              updateNupcolWellState().
            */
            WGState active_wgstate_;
            WGState last_valid_wgstate_;
            WGState nupcol_wgstate_;

        };


} // namespace Opm

#include "BlackoilWellModel_impl.hpp"
#endif
