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

#if HAVE_MPI
#define RESERVOIR_COUPLING_ENABLED
#endif
#ifdef RESERVOIR_COUPLING_ENABLED
#include <opm/simulators/flow/ReservoirCouplingMaster.hpp>
#include <opm/simulators/flow/ReservoirCouplingSlave.hpp>
#endif

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixmatrix.hh>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRate.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestState.hpp>

#include <opm/material/densead/Math.hpp>

#include <opm/simulators/flow/countGlobalCells.hpp>
#include <opm/simulators/flow/FlowBaseVanguard.hpp>

#include <opm/simulators/linalg/matrixblock.hh>

#include <opm/simulators/timestepping/SimulatorReport.hpp>
#include <opm/simulators/timestepping/gatherConvergenceReport.hpp>

#include <opm/simulators/utils/DeferredLogger.hpp>

#include <opm/simulators/wells/BlackoilWellModelGasLift.hpp>
#include <opm/simulators/wells/BlackoilWellModelGeneric.hpp>
#include <opm/simulators/wells/BlackoilWellModelGuideRates.hpp>
#include <opm/simulators/wells/GasLiftGroupInfo.hpp>
#include <opm/simulators/wells/GasLiftSingleWell.hpp>
#include <opm/simulators/wells/GasLiftSingleWellGeneric.hpp>
#include <opm/simulators/wells/GasLiftWellState.hpp>
#include <opm/simulators/wells/GuideRateHandler.hpp>
#include <opm/simulators/wells/MultisegmentWell.hpp>
#include <opm/simulators/wells/ParallelWBPCalculation.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>
#include <opm/simulators/wells/PerforationData.hpp>
#include <opm/simulators/wells/RateConverter.hpp>
#include <opm/simulators/wells/RegionAverageCalculator.hpp>
#include <opm/simulators/wells/StandardWell.hpp>
#include <opm/simulators/wells/VFPInjProperties.hpp>
#include <opm/simulators/wells/VFPProdProperties.hpp>
#include <opm/simulators/wells/WellConnectionAuxiliaryModule.hpp>
#include <opm/simulators/wells/WellGroupHelpers.hpp>
#include <opm/simulators/wells/WellInterface.hpp>
#include <opm/simulators/wells/WellProdIndexCalculator.hpp>
#include <opm/simulators/wells/WellState.hpp>
#include <opm/simulators/wells/WGState.hpp>

#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

namespace Opm {

template<class Scalar> class BlackoilWellModelNldd;
template<class T> class SparseTable;

#if COMPILE_GPU_BRIDGE
template<class Scalar> class WellContributions;
#endif

        /// Class for handling the blackoil well model.
        template<typename TypeTag>
        class BlackoilWellModel : public WellConnectionAuxiliaryModule<TypeTag, BlackoilWellModel<TypeTag>>
                                , public BlackoilWellModelGeneric<GetPropType<TypeTag,
                                                                              Properties::Scalar>>
        {
        public:
            // ---------      Types      ---------
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
            using ModelParameters = BlackoilModelParameters<Scalar>;

            using WellConnectionModule = WellConnectionAuxiliaryModule<TypeTag, BlackoilWellModel<TypeTag>>;

            constexpr static std::size_t pressureVarIndex = GetPropType<TypeTag, Properties::Indices>::pressureSwitchIdx;

            static const int numEq = Indices::numEq;
            static const int solventSaturationIdx = Indices::solventSaturationIdx;
            static constexpr bool has_solvent_ = getPropValue<TypeTag, Properties::EnableSolvent>();
            static constexpr bool has_polymer_ = getPropValue<TypeTag, Properties::EnablePolymer>();
            static constexpr bool has_energy_ = getPropValue<TypeTag, Properties::EnableEnergy>();
            static constexpr bool has_micp_ = getPropValue<TypeTag, Properties::EnableMICP>();

            // TODO: where we should put these types, WellInterface or Well Model?
            // or there is some other strategy, like TypeTag
            using VectorBlockType = Dune::FieldVector<Scalar, numEq>;
            using BVector = Dune::BlockVector<VectorBlockType>;

            using PolymerModule = BlackOilPolymerModule<TypeTag>;
            using MICPModule = BlackOilMICPModule<TypeTag>;

            // For the conversion between the surface volume rate and reservoir voidage rate
            using RateConverterType = RateConverter::
                SurfaceToReservoirVoidage<FluidSystem, std::vector<int> >;

            // For computing average pressured used by gpmaint
            using AverageRegionalPressureType = RegionAverageCalculator::
                AverageRegionalPressure<FluidSystem, std::vector<int> >;

            explicit BlackoilWellModel(Simulator& simulator);

            void init();
            void initWellContainer(const int reportStepIdx) override;

            void beginEpisode()
            {
                OPM_TIMEBLOCK(beginEpsiode);
                beginReportStep(simulator_.episodeIndex());
            }

            void beginTimeStep();

            void beginIteration()
            {
                OPM_TIMEBLOCK(beginIteration);
                assemble(simulator_.model().newtonMethod().numIterations(),
                         simulator_.timeStepSize());
            }

            void endIteration()
            { }

            void endTimeStep()
            {
                OPM_TIMEBLOCK(endTimeStep);
                timeStepSucceeded(simulator_.time(), simulator_.timeStepSize());
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

            using BlackoilWellModelGeneric<Scalar>::initFromRestartFile;
            void initFromRestartFile(const RestartValue& restartValues)
            {
                initFromRestartFile(restartValues,
                                    this->simulator_.vanguard().transferWTestState(),
                                    grid().size(0),
                                    param_.use_multisegment_well_,
                                    this->simulator_.vanguard().enableDistributedWells());
            }

            using BlackoilWellModelGeneric<Scalar>::prepareDeserialize;
            void prepareDeserialize(const int report_step)
            {
                prepareDeserialize(report_step, grid().size(0),
                                   param_.use_multisegment_well_,
                                   this->simulator_.vanguard().enableDistributedWells());
            }

            data::Wells wellData() const
            {
                auto wsrpt = this->wellState()
                    .report(this->simulator_.vanguard().globalCell().data(),
                            [this](const int well_index)
                {
                    return this->wasDynamicallyShutThisTimeStep(well_index);
                });

                BlackoilWellModelGuideRates(*this)
                    .assignWellGuideRates(wsrpt, this->reportStepIndex());

                this->assignWellTracerRates(wsrpt);

                if (const auto& rspec = eclState().runspec();
                    rspec.co2Storage() || rspec.h2Storage())
                {
                    // The gas reference density (surface condition) is the
                    // same for all PVT regions in CO2STORE/H2STORE runs so,
                    // for simplicity, we use region zero (0) here.

                    this->assignMassGasRate(wsrpt, FluidSystem::referenceDensity(FluidSystem::gasPhaseIdx, 0));
                }

                this->assignWellTargets(wsrpt);

                this->assignDynamicWellStatus(wsrpt, this->reportStepIndex());

                // Assigning (a subset of the) property values in shut
                // connections should be the last step of wellData().
                this->assignShutConnections(wsrpt, this->reportStepIndex());

                return wsrpt;
            }

            data::WellBlockAveragePressures wellBlockAveragePressures() const
            {
                return this->wbp_.computeWellBlockAveragePressures(this->gravity_);
            }

#if COMPILE_GPU_BRIDGE
            // accumulate the contributions of all Wells in the WellContributions object
            void getWellContributions(WellContributions<Scalar>& x) const;
#endif

            // Check if well equations is converged.
            ConvergenceReport getWellConvergence(const std::vector<Scalar>& B_avg, const bool checkWellGroupControls = false) const;

            const SimulatorReportSingle& lastReport() const;

            void addWellContributions(SparseMatrixAdapter& jacobian) const;

            // add source from wells to the reservoir matrix
            void addReservoirSourceTerms(GlobalEqVector& residual,
                                         const std::vector<typename SparseMatrixAdapter::MatrixBlock*>& diagMatAddress) const;

            // called at the beginning of a report step
            void beginReportStep(const int time_step);

            // it should be able to go to prepareTimeStep(), however, the updateWellControls()
            // makes it a little more difficult. unless we introduce if (iterationIdx != 0) to avoid doing the above function
            // twice at the beginning of the time step
            /// Calculating the explict quantities used in the well calculation. By explicit, we mean they are cacluated
            /// at the beginning of the time step and no derivatives are included in these quantities
            void calculateExplicitQuantities(DeferredLogger& deferred_logger) const;
            // some preparation work, mostly related to group control and RESV,
            // at the beginning of each time step (Not report step)
            void prepareTimeStep(DeferredLogger& deferred_logger);

            bool
            updateWellControls(DeferredLogger& deferred_logger);

            std::tuple<bool, Scalar>
            updateNetworks(const bool mandatory_network_balance, DeferredLogger& deferred_logger, const bool relax_network_tolerance = false);


            void updateAndCommunicate(const int reportStepIdx,
                                      const int iterationIdx,
                                      DeferredLogger& deferred_logger);

            bool updateGroupControls(const Group& group,
                                    DeferredLogger& deferred_logger,
                                    const int reportStepIdx,
                                    const int iterationIdx);

            WellInterfacePtr getWell(const std::string& well_name) const;

            using PressureMatrix = Dune::BCRSMatrix<Opm::MatrixBlock<Scalar, 1, 1>>;

            void addWellPressureEquations(PressureMatrix& jacobian,
                                          const BVector& weights,
                                          const bool use_well_weights) const;
            void addWellPressureEquationsStruct(PressureMatrix& jacobian) const;
            void addWellPressureEquationsDomain(PressureMatrix& jacobian,
                                                const BVector& weights,
                                                const bool use_well_weights,
                                                const int domainIndex) const
            {
                if (!nldd_) {
                    OPM_THROW(std::logic_error, "Attempt to access NLDD data without a NLDD solver");
                }
                return nldd_->addWellPressureEquations(jacobian,
                                                       weights,
                                                       use_well_weights,
                                                       domainIndex);
            }

            /// \brief Get list of local nonshut wells
            const std::vector<WellInterfacePtr>& localNonshutWells() const
            {
                return well_container_;
            }

            const SparseTable<int>& well_local_cells() const
            {
                if (!nldd_) {
                    OPM_THROW(std::logic_error, "Attempt to access NLDD data without a NLDD solver");
                }
                return nldd_->well_local_cells();
            }

            const std::map<std::string, int>& well_domain() const
            {
                if (!nldd_) {
                    OPM_THROW(std::logic_error, "Attempt to access NLDD data without a NLDD solver");
                }

                return nldd_->well_domain();
            }

            auto begin() const { return well_container_.begin(); }
            auto end() const { return well_container_.end(); }
            bool empty() const { return well_container_.empty(); }

            bool addMatrixContributions() const
            { return param_.matrix_add_well_contributions_; }

            int numStrictIterations() const
            { return param_.strict_outer_iter_wells_; }

            int compressedIndexForInterior(int cartesian_cell_idx) const override
            {
                return simulator_.vanguard().compressedIndexForInterior(cartesian_cell_idx);
            }

            int compressedIndexForInteriorLGR(const std::string& lgr_tag, const Connection& conn) const override
            {
                return simulator_.vanguard().compressedIndexForInteriorLGR(lgr_tag, conn);
            }            

            // using the solution x to recover the solution xw for wells and applying
            // xw to update Well State
            void recoverWellSolutionAndUpdateWellState(const BVector& x);

            // using the solution x to recover the solution xw for wells and applying
            // xw to update Well State
            void recoverWellSolutionAndUpdateWellStateDomain(const BVector& x,
                                                             const int domainIdx);

            const Grid& grid() const
            { return simulator_.vanguard().grid(); }

            const Simulator& simulator() const
            { return simulator_; }

            void setNlddAdapter(BlackoilWellModelNldd<TypeTag>* mod)
            { nldd_ = mod; }

#ifdef RESERVOIR_COUPLING_ENABLED
            ReservoirCouplingMaster& reservoirCouplingMaster() {
                return *(this->simulator_.reservoirCouplingMaster());
            }
            ReservoirCouplingSlave& reservoirCouplingSlave() {
                return *(this->simulator_.reservoirCouplingSlave());
            }
            bool isReservoirCouplingMaster() const {
                return this->simulator_.reservoirCouplingMaster() != nullptr;
            }
            bool isReservoirCouplingSlave() const {
                return this->simulator_.reservoirCouplingSlave() != nullptr;
            }
            void setReservoirCouplingMaster(ReservoirCouplingMaster *master)
            {
                this->guide_rate_handler_.setReservoirCouplingMaster(master);
            }
            void setReservoirCouplingSlave(ReservoirCouplingSlave *slave)
            {
                this->guide_rate_handler_.setReservoirCouplingSlave(slave);
            }
        #endif
        protected:
            Simulator& simulator_;

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

            WellInterfacePtr createWellForWellTest(const std::string& well_name,
                                                   const int report_step,
                                                   DeferredLogger& deferred_logger) const;

            const ModelParameters param_;
            std::size_t global_num_cells_{};
            // the number of the cells in the local grid
            std::size_t local_num_cells_{};
            Scalar gravity_{};
            std::vector<Scalar> depth_{};
            bool alternative_well_rate_init_{};
            std::map<std::string, Scalar> well_group_thp_calc_;
            std::unique_ptr<RateConverterType> rateConverter_{};
            std::map<std::string, std::unique_ptr<AverageRegionalPressureType>> regionalAveragePressureCalculator_{};

            SimulatorReportSingle last_report_{};
            GuideRateHandler<Scalar> guide_rate_handler_{};

            // Pre-step network solve at static reservoir conditions (group and well states might be updated)
            void doPreStepNetworkRebalance(DeferredLogger& deferred_logger);

            std::vector<Scalar> B_avg_{};

            const EquilGrid& equilGrid() const
            { return simulator_.vanguard().equilGrid(); }

            const EclipseState& eclState() const
            { return simulator_.vanguard().eclState(); }

            // compute the well fluxes and assemble them in to the reservoir equations as source terms
            // and in the well equations.
            void assemble(const int iterationIdx,
                          const double dt);

            // well controls and network pressures affect each other and are solved in an iterative manner.
            // the function handles one iteration of updating well controls and network pressures.
            // it is possible to decouple the update of well controls and network pressures further.
            // the returned two booleans are {continue_due_to_network, well_group_control_changed}, respectively
            std::tuple<bool, bool, Scalar> updateWellControlsAndNetworkIteration(const bool mandatory_network_balance,
                                                                        const bool relax_network_tolerance,
                                                                        const bool optimize_gas_lift,
                                                                        const double dt,
                                                                        DeferredLogger& local_deferredLogger);

            bool updateWellControlsAndNetwork(const bool mandatory_network_balance,
                                              const double dt,
                                              DeferredLogger& local_deferredLogger);

            bool computeWellGroupThp(const double dt, DeferredLogger& local_deferredLogger);

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

            // setting the well_solutions_ based on well_state.
            void updatePrimaryVariables(DeferredLogger& deferred_logger);

            void updateAverageFormationFactor();

            void computePotentials(const std::size_t widx,
                                   const WellState<Scalar>& well_state_copy,
                                   std::string& exc_msg,
                                   ExceptionType::ExcEnum& exc_type,
                                   DeferredLogger& deferred_logger) override;

            const std::vector<Scalar>& wellPerfEfficiencyFactors() const;

            void calculateProductivityIndexValuesShutWells(const int reportStepIdx, DeferredLogger& deferred_logger) override;
            void calculateProductivityIndexValues(DeferredLogger& deferred_logger) override;
            void calculateProductivityIndexValues(const WellInterface<TypeTag>* wellPtr,
                                                  DeferredLogger& deferred_logger);

            // The number of components in the model.
            int numComponents() const;

            int reportStepIndex() const;

            void assembleWellEq(const double dt, DeferredLogger& deferred_logger);

            void prepareWellsBeforeAssembling(const double dt, DeferredLogger& deferred_logger);

            // TODO: finding a better naming
            void assembleWellEqWithoutIteration(const double dt, DeferredLogger& deferred_logger);

            void extractLegacyCellPvtRegionIndex_();

            void extractLegacyDepth_();

            /// upate the wellTestState related to economic limits
            void updateWellTestState(const double& simulationTime, WellTestState& wellTestState) const;

            void wellTesting(const int timeStepIdx, const double simulationTime, DeferredLogger& deferred_logger);

            void calcResvCoeff(const int fipnum,
                               const int pvtreg,
                               const std::vector<Scalar>& production_rates,
                               std::vector<Scalar>& resv_coeff) override;

            void calcInjResvCoeff(const int fipnum,
                                  const int pvtreg,
                                  std::vector<Scalar>& resv_coeff) override;

            void computeWellTemperature();

        private:
            BlackoilWellModel(Simulator& simulator, const PhaseUsage& pu);

            BlackoilWellModelGasLift<TypeTag> gaslift_;
            BlackoilWellModelNldd<TypeTag>* nldd_ = nullptr; //!< NLDD well model adapter (not owned)

            // These members are used to avoid reallocation in specific functions
            // instead of using local variables.
            // Their state is not relevant between function calls, so they can
            // (and must) be mutable, as the functions using them are const.
            mutable BVector x_local_;

            void assignWellTracerRates(data::Wells& wsrpt) const;
        };

} // namespace Opm

#include "BlackoilWellModel_impl.hpp"

#endif // OPM_BLACKOILWELLMODEL_HEADER_INCLUDED
