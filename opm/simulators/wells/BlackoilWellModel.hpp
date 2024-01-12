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

#include <cassert>
#include <cstddef>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRate.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestState.hpp>

#include <opm/simulators/flow/countGlobalCells.hpp>
#include <opm/simulators/flow/SubDomain.hpp>

#include <opm/simulators/wells/BlackoilWellModelGeneric.hpp>
#include <opm/simulators/wells/BlackoilWellModelGuideRates.hpp>
#include <opm/simulators/wells/GasLiftGroupInfo.hpp>
#include <opm/simulators/wells/GasLiftSingleWell.hpp>
#include <opm/simulators/wells/GasLiftSingleWellGeneric.hpp>
#include <opm/simulators/wells/GasLiftWellState.hpp>
#include <opm/simulators/wells/MultisegmentWell.hpp>
#include <opm/simulators/wells/ParallelWBPCalculation.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>
#include <opm/simulators/wells/PerforationData.hpp>
#include <opm/simulators/wells/RateConverter.hpp>
#include <opm/simulators/wells/RegionAverageCalculator.hpp>
#include <opm/simulators/wells/StandardWell.hpp>
#include <opm/simulators/wells/VFPInjProperties.hpp>
#include <opm/simulators/wells/VFPProdProperties.hpp>
#include <opm/simulators/wells/WGState.hpp>
#include <opm/simulators/wells/WellGroupHelpers.hpp>
#include <opm/simulators/wells/WellInterface.hpp>
#include <opm/simulators/wells/WellProdIndexCalculator.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <opm/simulators/timestepping/SimulatorReport.hpp>
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
                                , public BlackoilWellModelGeneric
        {
        public:
            // ---------      Types      ---------
            typedef BlackoilModelParametersEbos<TypeTag> ModelParameters;

            using Grid = GetPropType<TypeTag, Properties::Grid>;
            using EquilGrid = GetPropType<TypeTag, Properties::EquilGrid>;
            using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
            using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
            using Indices = GetPropType<TypeTag, Properties::Indices>;
            using Simulator = GetPropType<TypeTag, Properties::Simulator>;
            using Scalar = GetPropType<TypeTag, Properties::Scalar>;
            using RateVector = GetPropType<TypeTag, Properties::RateVector>;
            using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
            using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
            using GasLiftSingleWell = typename WellInterface<TypeTag>::GasLiftSingleWell;
            using GLiftOptWells = typename BlackoilWellModelGeneric::GLiftOptWells;
            using GLiftProdWells = typename BlackoilWellModelGeneric::GLiftProdWells;
            using GLiftWellStateMap =
                typename BlackoilWellModelGeneric::GLiftWellStateMap;
            using GLiftEclWells = typename GasLiftGroupInfo::GLiftEclWells;
            using GLiftSyncGroups = typename GasLiftSingleWellGeneric::GLiftSyncGroups;
            constexpr static std::size_t pressureVarIndex = GetPropType<TypeTag, Properties::Indices>::pressureSwitchIdx;
            typedef typename BaseAuxiliaryModule<TypeTag>::NeighborSet NeighborSet;

            static const int numEq = Indices::numEq;
            static const int solventSaturationIdx = Indices::solventSaturationIdx;
            static constexpr bool has_solvent_ = getPropValue<TypeTag, Properties::EnableSolvent>();
            static constexpr bool has_polymer_ = getPropValue<TypeTag, Properties::EnablePolymer>();
            static constexpr bool has_energy_ = getPropValue<TypeTag, Properties::EnableEnergy>();
            static constexpr bool has_micp_ = getPropValue<TypeTag, Properties::EnableMICP>();

            // TODO: where we should put these types, WellInterface or Well Model?
            // or there is some other strategy, like TypeTag
            typedef Dune::FieldVector<Scalar, numEq    > VectorBlockType;
            typedef Dune::BlockVector<VectorBlockType> BVector;

            typedef BlackOilPolymerModule<TypeTag> PolymerModule;
            typedef BlackOilMICPModule<TypeTag> MICPModule;

            // For the conversion between the surface volume rate and resrevoir voidage rate
            using RateConverterType = RateConverter::
                SurfaceToReservoirVoidage<FluidSystem, std::vector<int> >;

            // For computing average pressured used by gpmaint
            using AverageRegionalPressureType = RegionAverageCalculator::
                AverageRegionalPressure<FluidSystem, std::vector<int> >;

            using Domain = SubDomain<Grid>;

            BlackoilWellModel(Simulator& ebosSimulator);

            void init();
            void initWellContainer(const int reportStepIdx) override;

            /////////////
            // <eWoms auxiliary module stuff>
            /////////////
            unsigned numDofs() const override
            // No extra dofs are inserted for wells. (we use a Schur complement.)
            { return 0; }

            void addNeighbors(std::vector<NeighborSet>& neighbors) const override;

            void applyInitial() override
            {}

            void linearize(SparseMatrixAdapter& jacobian, GlobalEqVector& res) override;
            void linearizeDomain(const Domain& domain, SparseMatrixAdapter& jacobian, GlobalEqVector& res);

            void postSolve(GlobalEqVector& deltaX) override
            {
                recoverWellSolutionAndUpdateWellState(deltaX);
            }

            void postSolveDomain(GlobalEqVector& deltaX, const Domain& domain)
            {
                recoverWellSolutionAndUpdateWellStateDomain(deltaX, domain);
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
                OPM_TIMEBLOCK(beginEpsiode);
                beginReportStep(ebosSimulator_.episodeIndex());
            }

            void beginTimeStep();

            void beginIteration()
            {
                OPM_TIMEBLOCK(beginIteration);
                assemble(ebosSimulator_.model().newtonMethod().numIterations(),
                         ebosSimulator_.timeStepSize());
            }

            void endIteration()
            { }

            void endTimeStep()
            {
                OPM_TIMEBLOCK(endTimeStep);
                timeStepSucceeded(ebosSimulator_.time(), ebosSimulator_.timeStepSize());
            }

            void endEpisode()
            {
                endReportStep();
            }

            void computeTotalRatesForDof(RateVector& rate,
                                         unsigned globalIdx) const;

            template <class Context>
            void computeTotalRatesForDof(RateVector& rate,
                                         const Context& context,
                                         unsigned spaceIdx,
                                         unsigned timeIdx) const;


            using WellInterfacePtr = std::shared_ptr<WellInterface<TypeTag> >;

            using BlackoilWellModelGeneric::initFromRestartFile;
            void initFromRestartFile(const RestartValue& restartValues)
            {
                initFromRestartFile(restartValues,
                                    this->ebosSimulator_.vanguard().transferWTestState(),
                                    grid().size(0),
                                    param_.use_multisegment_well_);
            }

            using BlackoilWellModelGeneric::prepareDeserialize;
            void prepareDeserialize(const int report_step)
            {
                prepareDeserialize(report_step, grid().size(0),
                                   param_.use_multisegment_well_);
            }

            data::Wells wellData() const
            {
                auto wsrpt = this->wellState()
                    .report(ebosSimulator_.vanguard().globalCell().data(),
                            [this](const int well_index) -> bool
                {
                    return this->wasDynamicallyShutThisTimeStep(well_index);
                });

                const auto& tracerRates = ebosSimulator_.problem().tracerModel().getWellTracerRates();
                this->assignWellTracerRates(wsrpt, tracerRates);


                BlackoilWellModelGuideRates(*this).assignWellGuideRates(wsrpt, this->reportStepIndex());
                this->assignShutConnections(wsrpt, this->reportStepIndex());

                return wsrpt;
            }

            data::WellBlockAveragePressures wellBlockAveragePressures() const
            {
                return this->computeWellBlockAveragePressures();
            }

            // subtract Binv(D)rw from r;
            void apply( BVector& r) const;

            // subtract B*inv(D)*C * x from A*x
            void apply(const BVector& x, BVector& Ax) const;

            // accumulate the contributions of all Wells in the WellContributions object
            void getWellContributions(WellContributions& x) const;

            // apply well model with scaling of alpha
            void applyScaleAdd(const Scalar alpha, const BVector& x, BVector& Ax) const;

            // Check if well equations is converged.
            ConvergenceReport getWellConvergence(const std::vector<Scalar>& B_avg, const bool checkWellGroupControls = false) const;

            // Check if well equations are converged locally.
            ConvergenceReport getDomainWellConvergence(const Domain& domain,
                                                       const std::vector<Scalar>& B_avg,
                                                       DeferredLogger& local_deferredLogger) const;

            const SimulatorReportSingle& lastReport() const;

            void addWellContributions(SparseMatrixAdapter& jacobian) const;

            // add source from wells to the reservoir matrix
            void addReservoirSourceTerms(GlobalEqVector& residual,
                                         std::vector<typename SparseMatrixAdapter::MatrixBlock*>& diagMatAddress) const;

            // called at the beginning of a report step
            void beginReportStep(const int time_step);

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
            void initPrimaryVariablesEvaluationDomain(const Domain& domain) const;

            std::pair<bool, bool>
            updateWellControls(const bool mandatory_network_balance, DeferredLogger& deferred_logger, const bool relax_network_tolerance = false);

            void updateAndCommunicate(const int reportStepIdx,
                                      const int iterationIdx,
                                      DeferredLogger& deferred_logger);

            bool updateGroupControls(const Group& group,
                                    DeferredLogger& deferred_logger,
                                    const int reportStepIdx,
                                    const int iterationIdx);

            WellInterfacePtr getWell(const std::string& well_name) const;
            bool hasWell(const std::string& well_name) const;

            using PressureMatrix = Dune::BCRSMatrix<Opm::MatrixBlock<double, 1, 1>>;

            void addWellPressureEquations(PressureMatrix& jacobian, const BVector& weights,const bool use_well_weights) const;

            void addWellPressureEquationsStruct(PressureMatrix& jacobian) const;

            void initGliftEclWellMap(GLiftEclWells &ecl_well_map);

            /// \brief Get list of local nonshut wells
            const std::vector<WellInterfacePtr>& localNonshutWells() const
            {
                return well_container_;
            }

            // prototype for assemble function for ASPIN solveLocal()
            // will try to merge back to assemble() when done prototyping
            void assembleDomain(const int iterationIdx,
                                const double dt,
                                const Domain& domain);
            void updateWellControlsDomain(DeferredLogger& deferred_logger, const Domain& domain);

            void logPrimaryVars() const;
            std::vector<double> getPrimaryVarsDomain(const Domain& domain) const;
            void setPrimaryVarsDomain(const Domain& domain, const std::vector<double>& vars);

            void setupDomains(const std::vector<Domain>& domains);

        protected:
            Simulator& ebosSimulator_;

            // a vector of all the wells.
            std::vector<WellInterfacePtr> well_container_{};

            std::vector<bool> is_cell_perforated_{};

            void initializeWellState(const int timeStepIdx);

            // create the well container
            void createWellContainer(const int report_step) override;

            WellInterfacePtr
            createWellPointer(const int wellID,
                              const int report_step) const;

            template <typename WellType>
            std::unique_ptr<WellType>
            createTypedWellPointer(const int wellID,
                                   const int time_step) const;

            WellInterfacePtr createWellForWellTest(const std::string& well_name, const int report_step, DeferredLogger& deferred_logger) const;


            const ModelParameters param_;
            std::size_t global_num_cells_{};
            // the number of the cells in the local grid
            std::size_t local_num_cells_{};
            double gravity_{};
            std::vector<double> depth_{};
            bool alternative_well_rate_init_{};

            std::unique_ptr<RateConverterType> rateConverter_{};
            std::map<std::string, std::unique_ptr<AverageRegionalPressureType>> regionalAveragePressureCalculator_{};

            struct WBPCalcID
            {
                std::optional<typename std::vector<WellInterfacePtr>::size_type> openWellIdx_{};
                std::size_t wbpCalcIdx_{};
            };

            std::vector<WBPCalcID> wbpCalcMap_{};

            SimulatorReportSingle last_report_{};

            // Pre-step network solve at static reservoir conditions (group and well states might be updated)
            void doPreStepNetworkRebalance(DeferredLogger& deferred_logger);

            // used to better efficiency of calcuation
            mutable BVector scaleAddRes_{};

            std::vector<Scalar> B_avg_{};

            // Keep track of the domain of each well, if using subdomains.
            std::map<std::string, int> well_domain_;

            const Grid& grid() const
            { return ebosSimulator_.vanguard().grid(); }

            const EquilGrid& equilGrid() const
            { return ebosSimulator_.vanguard().equilGrid(); }

            const EclipseState& eclState() const
            { return ebosSimulator_.vanguard().eclState(); }

            // compute the well fluxes and assemble them in to the reservoir equations as source terms
            // and in the well equations.
            void assemble(const int iterationIdx,
                          const double dt);

            // well controls and network pressures affect each other and are solved in an iterative manner.
            // the function handles one iteration of updating well controls and network pressures.
            // it is possible to decouple the update of well controls and network pressures further.
            // the returned two booleans are {continue_due_to_network, well_group_control_changed}, respectively
            std::pair<bool, bool> updateWellControlsAndNetworkIteration(const bool mandatory_network_balance,
                                                                        const bool relax_network_tolerance,
                                                                        const double dt,
                                                                        DeferredLogger& local_deferredLogger);

            bool updateWellControlsAndNetwork(const bool mandatory_network_balance,
                                              const double dt,
                                              DeferredLogger& local_deferredLogger);

            /// Update rank's notion of intersecting wells and their
            /// associate solution variables.
            ///
            /// \param[in] reportStepIdx Report step.
            ///
            /// \param[in] enableWellPIScaling Whether or not to enable WELPI
            ///   scaling.  Typically enabled (i.e., true) only at the start
            ///   of a report step.
            void initializeLocalWellStructure(const int  reportStepIdx,
                                              const bool enableWellPIScaling);

            /// Initialize group control modes/constraints and group solution state.
            ///
            /// \param[in] reportStepIdx Report step.
            void initializeGroupStructure(const int reportStepIdx);

            // called at the end of a time step
            void timeStepSucceeded(const double simulationTime, const double dt);

            // called at the end of a report step
            void endReportStep();

            // using the solution x to recover the solution xw for wells and applying
            // xw to update Well State
            void recoverWellSolutionAndUpdateWellState(const BVector& x);

            // using the solution x to recover the solution xw for wells and applying
            // xw to update Well State
            void recoverWellSolutionAndUpdateWellStateDomain(const BVector& x, const Domain& domain);

            // setting the well_solutions_ based on well_state.
            void updatePrimaryVariables(DeferredLogger& deferred_logger);

            void initializeWBPCalculationService();

            data::WellBlockAveragePressures
            computeWellBlockAveragePressures() const;

            ParallelWBPCalculation::EvaluatorFactory
            makeWellSourceEvaluatorFactory(const std::vector<Well>::size_type wellIdx) const;

            void registerOpenWellsForWBPCalculation();

            void updateAverageFormationFactor();

            void computePotentials(const std::size_t widx,
                                   const WellState& well_state_copy,
                                   std::string& exc_msg,
                                   ExceptionType::ExcEnum& exc_type,
                                   DeferredLogger& deferred_logger) override;

            const std::vector<double>& wellPerfEfficiencyFactors() const;

            void calculateProductivityIndexValuesShutWells(const int reportStepIdx, DeferredLogger& deferred_logger) override;
            void calculateProductivityIndexValues(DeferredLogger& deferred_logger) override;
            void calculateProductivityIndexValues(const WellInterface<TypeTag>* wellPtr,
                                                  DeferredLogger& deferred_logger);

            // The number of components in the model.
            int numComponents() const;

            int reportStepIndex() const;

            void assembleWellEq(const double dt, DeferredLogger& deferred_logger);
            void assembleWellEqDomain(const double dt, const Domain& domain, DeferredLogger& deferred_logger);

            void prepareWellsBeforeAssembling(const double dt, DeferredLogger& deferred_logger);

            // TODO: finding a better naming
            void assembleWellEqWithoutIteration(const double dt, DeferredLogger& deferred_logger);

            bool maybeDoGasLiftOptimize(DeferredLogger& deferred_logger);

            void gasLiftOptimizationStage1(DeferredLogger& deferred_logger,
                GLiftProdWells &prod_wells, GLiftOptWells &glift_wells,
                GasLiftGroupInfo &group_info, GLiftWellStateMap &state_map);

            // cannot be const since it accesses the non-const WellState
            void gasLiftOptimizationStage1SingleWell(WellInterface<TypeTag> *well,
                DeferredLogger& deferred_logger,
                GLiftProdWells &prod_wells, GLiftOptWells &glift_wells,
                GasLiftGroupInfo &group_info, GLiftWellStateMap &state_map,
                GLiftSyncGroups& groups_to_sync);

            void extractLegacyCellPvtRegionIndex_();

            void extractLegacyDepth_();

            /// upate the wellTestState related to economic limits
            void updateWellTestState(const double& simulationTime, WellTestState& wellTestState) const;

            void wellTesting(const int timeStepIdx, const double simulationTime, DeferredLogger& deferred_logger);

            void calcRates(const int fipnum,
                           const int pvtreg,
                           const std::vector<double>& production_rates,
                           std::vector<double>& resv_coeff) override;

            void calcInjRates(const int fipnum,
                           const int pvtreg,
                           std::vector<double>& resv_coeff) override;

            void computeWellTemperature();

            int compressedIndexForInterior(int cartesian_cell_idx) const override {
                return ebosSimulator_.vanguard().compressedIndexForInterior(cartesian_cell_idx);
            }

        private:
            BlackoilWellModel(Simulator& ebosSimulator, const PhaseUsage& pu);
        };


} // namespace Opm

#include "BlackoilWellModel_impl.hpp"
#endif
