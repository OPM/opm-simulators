/*
  Copyright 2016 - 2019 SINTEF Digital, Mathematics & Cybernetics.
  Copyright 2016 - 2018 Equinor ASA.
  Copyright 2017 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2016 - 2018 Norce AS

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

#ifndef OPM_BLACKOILWELLMODEL_IMPL_HEADER_INCLUDED
#define OPM_BLACKOILWELLMODEL_IMPL_HEADER_INCLUDED

// Improve IDE experience
#ifndef OPM_BLACKOILWELLMODEL_HEADER_INCLUDED
#include <config.h>
#include <opm/simulators/wells/BlackoilWellModel.hpp>
#endif

#include <opm/grid/utility/cartesianToCompressed.hpp>
#include <opm/common/utility/numeric/RootFinders.hpp>

#include <opm/input/eclipse/Schedule/Network/Balance.hpp>
#include <opm/input/eclipse/Schedule/Network/ExtNetwork.hpp>
#include <opm/input/eclipse/Schedule/Well/PAvgDynamicSourceData.hpp>
#include <opm/input/eclipse/Schedule/Well/WellMatcher.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestConfig.hpp>

#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <opm/simulators/wells/BlackoilWellModelConstraints.hpp>
#include <opm/simulators/wells/BlackoilWellModelNldd.hpp>
#include <opm/simulators/wells/ParallelPAvgDynamicSourceData.hpp>
#include <opm/simulators/wells/ParallelWBPCalculation.hpp>
#include <opm/simulators/wells/VFPProperties.hpp>
#include <opm/simulators/wells/WellBhpThpCalculator.hpp>
#include <opm/simulators/wells/WellGroupControls.hpp>
#include <opm/simulators/wells/WellGroupHelpers.hpp>
#include <opm/simulators/wells/TargetCalculator.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/utils/MPIPacker.hpp>
#include <opm/simulators/utils/phaseUsageFromDeck.hpp>

#if COMPILE_GPU_BRIDGE
#include <opm/simulators/linalg/gpubridge/WellContributions.hpp>
#endif

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <utility>
#include <optional>

#include <fmt/format.h>

namespace Opm {
    template<typename TypeTag>
    BlackoilWellModel<TypeTag>::
    BlackoilWellModel(Simulator& simulator, const PhaseUsage& phase_usage)
        : WellConnectionModule(*this, simulator.gridView().comm())
        , BlackoilWellModelGeneric<Scalar>(simulator.vanguard().schedule(),
                                           gaslift_,
                                           simulator.vanguard().summaryState(),
                                           simulator.vanguard().eclState(),
                                           phase_usage,
                                           simulator.gridView().comm())
        , simulator_(simulator)
        , gaslift_(this->terminal_output_, this->phase_usage_)
    {
        local_num_cells_ = simulator_.gridView().size(0);

        // Number of cells the global grid view
        global_num_cells_ = simulator_.vanguard().globalNumCells();

        {
            auto& parallel_wells = simulator.vanguard().parallelWells();

            this->parallel_well_info_.reserve(parallel_wells.size());
            for( const auto& name_bool : parallel_wells) {
                this->parallel_well_info_.emplace_back(name_bool, grid().comm());
            }
        }

        this->alternative_well_rate_init_ =
            Parameters::Get<Parameters::AlternativeWellRateInit>();

        using SourceDataSpan =
            typename PAvgDynamicSourceData<Scalar>::template SourceDataSpan<Scalar>;

        this->wbp_.initializeSources(
            [this](const std::size_t globalIndex)
            { return this->compressedIndexForInterior(globalIndex); },
            [this](const int localCell, SourceDataSpan sourceTerms)
            {
                using Item = typename SourceDataSpan::Item;

                const auto* intQuants = this->simulator_.model()
                    .cachedIntensiveQuantities(localCell, /*timeIndex = */0);
                const auto& fs = intQuants->fluidState();

                sourceTerms
                    .set(Item::PoreVol, intQuants->porosity().value() *
                         this->simulator_.model().dofTotalVolume(localCell))
                    .set(Item::Depth, this->depth_[localCell]);

                constexpr auto io = FluidSystem::oilPhaseIdx;
                constexpr auto ig = FluidSystem::gasPhaseIdx;
                constexpr auto iw = FluidSystem::waterPhaseIdx;

                // Ideally, these would be 'constexpr'.
                const auto haveOil = FluidSystem::phaseIsActive(io);
                const auto haveGas = FluidSystem::phaseIsActive(ig);
                const auto haveWat = FluidSystem::phaseIsActive(iw);

                auto weightedPhaseDensity = [&fs](const auto ip)
                {
                    return fs.saturation(ip).value() * fs.density(ip).value();
                };

                if (haveOil)      { sourceTerms.set(Item::Pressure, fs.pressure(io).value()); }
                else if (haveGas) { sourceTerms.set(Item::Pressure, fs.pressure(ig).value()); }
                else              { sourceTerms.set(Item::Pressure, fs.pressure(iw).value()); }

                // Strictly speaking, assumes SUM(s[p]) == 1.
                auto rho = 0.0;
                if (haveOil) { rho += weightedPhaseDensity(io); }
                if (haveGas) { rho += weightedPhaseDensity(ig); }
                if (haveWat) { rho += weightedPhaseDensity(iw); }

                sourceTerms.set(Item::MixtureDensity, rho);
            }
        );
    }

    template<typename TypeTag>
    BlackoilWellModel<TypeTag>::
    BlackoilWellModel(Simulator& simulator) :
        BlackoilWellModel(simulator, phaseUsageFromDeck(simulator.vanguard().eclState()))
    {}


    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    init()
    {
        extractLegacyCellPvtRegionIndex_();
        extractLegacyDepth_();

        gravity_ = simulator_.problem().gravity()[2];

        this->initial_step_ = true;

        // add the eWoms auxiliary module for the wells to the list
        simulator_.model().addAuxiliaryModule(this);

        is_cell_perforated_.resize(local_num_cells_, false);
    }


    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    initWellContainer(const int reportStepIdx)
    {
        const uint64_t effective_events_mask = ScheduleEvents::WELL_STATUS_CHANGE
                        + ScheduleEvents::NEW_WELL;
        const auto& events = this->schedule()[reportStepIdx].wellgroup_events();
        for (auto& wellPtr : this->well_container_) {
            const bool well_opened_this_step = this->report_step_starts_ &&
                                               events.hasEvent(wellPtr->name(),
                                                               effective_events_mask);
            wellPtr->init(&this->phase_usage_, this->depth_, this->gravity_,
                          this->B_avg_, well_opened_this_step);
        }
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    beginReportStep(const int timeStepIdx)
    {
        DeferredLogger local_deferredLogger{};

        this->report_step_starts_ = true;
        this->report_step_start_events_ = this->schedule()[timeStepIdx].wellgroup_events();

        this->rateConverter_ = std::make_unique<RateConverterType>
            (this->phase_usage_, std::vector<int>(this->local_num_cells_, 0));

        {
            // WELPI scaling runs at start of report step.
            const auto enableWellPIScaling = true;
            this->initializeLocalWellStructure(timeStepIdx, enableWellPIScaling);
        }

        this->initializeGroupStructure(timeStepIdx);

        const auto& comm = this->simulator_.vanguard().grid().comm();

        OPM_BEGIN_PARALLEL_TRY_CATCH()
        {
            // Create facility for calculating reservoir voidage volumes for
            // purpose of RESV controls.
            this->rateConverter_->template defineState<ElementContext>(this->simulator_);

            // Update VFP properties.
            {
                const auto& sched_state = this->schedule()[timeStepIdx];

                this->vfp_properties_ = std::make_unique<VFPProperties<Scalar>>
                    (sched_state.vfpinj(), sched_state.vfpprod(), this->wellState());
            }
        }
        OPM_END_PARALLEL_TRY_CATCH_LOG(local_deferredLogger,
                                       "beginReportStep() failed: ",
                                       this->terminal_output_, comm)

        // Store the current well and group states in order to recover in
        // the case of failed iterations
        this->commitWGState();

        this->wellStructureChangedDynamically_ = false;
    }





    template <typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    initializeLocalWellStructure(const int  reportStepIdx,
                                 const bool enableWellPIScaling)
    {
        DeferredLogger local_deferredLogger{};

        const auto& comm = this->simulator_.vanguard().grid().comm();

        // Wells_ecl_ holds this rank's wells, both open and stopped/shut.
        this->wells_ecl_ = this->getLocalWells(reportStepIdx);
        this->local_parallel_well_info_ =
            this->createLocalParallelWellInfo(this->wells_ecl_);

        // At least initializeWellState() might be throw an exception in
        // UniformTabulated2DFunction.  Playing it safe by extending the
        // scope a bit.
        OPM_BEGIN_PARALLEL_TRY_CATCH()
        {
            this->initializeWellPerfData();
            this->initializeWellState(reportStepIdx);
            this->wbp_.initializeWBPCalculationService();

            if (this->param_.use_multisegment_well_ && this->anyMSWellOpenLocal()) {
                this->wellState().initWellStateMSWell(this->wells_ecl_, &this->prevWellState());
            }

            this->initializeWellProdIndCalculators();

            if (enableWellPIScaling && this->schedule()[reportStepIdx].events()
                .hasEvent(ScheduleEvents::Events::WELL_PRODUCTIVITY_INDEX))
            {
                this->runWellPIScaling(reportStepIdx, local_deferredLogger);
            }
        }
        OPM_END_PARALLEL_TRY_CATCH_LOG(local_deferredLogger,
                                       "Failed to initialize local well structure: ",
                                       this->terminal_output_, comm)
    }





    template <typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    initializeGroupStructure(const int reportStepIdx)
    {
        DeferredLogger local_deferredLogger{};

        const auto& comm = this->simulator_.vanguard().grid().comm();

        OPM_BEGIN_PARALLEL_TRY_CATCH()
        {
            const auto& fieldGroup =
                this->schedule().getGroup("FIELD", reportStepIdx);

            WellGroupHelpers<Scalar>::setCmodeGroup(fieldGroup,
                                                    this->schedule(),
                                                    this->summaryState(),
                                                    reportStepIdx,
                                                    this->groupState());

            // Define per region average pressure calculators for use by
            // pressure maintenance groups (GPMAINT keyword).
            if (this->schedule()[reportStepIdx].has_gpmaint()) {
                WellGroupHelpers<Scalar>::setRegionAveragePressureCalculator
                    (fieldGroup,
                     this->schedule(),
                     reportStepIdx,
                     this->eclState_.fieldProps(),
                     this->phase_usage_,
                     this->regionalAveragePressureCalculator_);
            }
        }
        OPM_END_PARALLEL_TRY_CATCH_LOG(local_deferredLogger,
                                       "Failed to initialize group structure: ",
                                       this->terminal_output_, comm)
    }





    // called at the beginning of a time step
    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    beginTimeStep()
    {
        OPM_TIMEBLOCK(beginTimeStep);

        this->updateAverageFormationFactor();

        DeferredLogger local_deferredLogger;

        this->switched_prod_groups_.clear();
        this->switched_inj_groups_.clear();

        if (this->wellStructureChangedDynamically_) {
            // Something altered the well structure/topology.  Possibly
            // WELSPECS/COMPDAT and/or WELOPEN run from an ACTIONX block.
            // Reconstruct the local wells to account for the new well
            // structure.
            const auto reportStepIdx =
                this->simulator_.episodeIndex();

            // Disable WELPI scaling when well structure is updated in the
            // middle of a report step.
            const auto enableWellPIScaling = false;

            this->initializeLocalWellStructure(reportStepIdx, enableWellPIScaling);
            this->initializeGroupStructure(reportStepIdx);

            this->commitWGState();

            // Reset topology flag to signal that we've handled this
            // structure change.  That way we don't end up here in
            // subsequent calls to beginTimeStep() unless there's a new
            // dynamic change to the well structure during a report step.
            this->wellStructureChangedDynamically_ = false;
        }

        this->resetWGState();

        const int reportStepIdx = simulator_.episodeIndex();
        this->updateAndCommunicateGroupData(reportStepIdx,
                                            simulator_.model().newtonMethod().numIterations(),
                                            param_.nupcol_group_rate_tolerance_,
                                            local_deferredLogger);

        this->wellState().updateWellsDefaultALQ(this->schedule(), reportStepIdx, this->summaryState());
        this->wellState().gliftTimeStepInit();

        const double simulationTime = simulator_.time();
        OPM_BEGIN_PARALLEL_TRY_CATCH();
        {
            // test wells
            wellTesting(reportStepIdx, simulationTime, local_deferredLogger);

            // create the well container
            createWellContainer(reportStepIdx);

            // Wells are active if they are active wells on at least one process.
            const Grid& grid = simulator_.vanguard().grid();
            this->wells_active_ = grid.comm().max(!this->well_container_.empty());

            // do the initialization for all the wells
            // TODO: to see whether we can postpone of the intialization of the well containers to
            // optimize the usage of the following several member variables
            this->initWellContainer(reportStepIdx);

            // update the updated cell flag
            std::fill(is_cell_perforated_.begin(), is_cell_perforated_.end(), false);
            for (auto& well : well_container_) {
                well->updatePerforatedCell(is_cell_perforated_);
            }

            // calculate the efficiency factors for each well
            this->calculateEfficiencyFactors(reportStepIdx);

            if constexpr (has_polymer_)
            {
                if (PolymerModule::hasPlyshlog() || getPropValue<TypeTag, Properties::EnablePolymerMW>() ) {
                    this->setRepRadiusPerfLength();
                }
            }

        }
        OPM_END_PARALLEL_TRY_CATCH_LOG(local_deferredLogger, "beginTimeStep() failed: ",
                                        this->terminal_output_, simulator_.vanguard().grid().comm());

        for (auto& well : well_container_) {
            well->setVFPProperties(this->vfp_properties_.get());
            well->setGuideRate(&this->guideRate_);
        }

        this->updateFiltrationModelsPreStep(local_deferredLogger);

        // Close completions due to economic reasons
        for (auto& well : well_container_) {
            well->closeCompletions(this->wellTestState());
        }

        // we need the inj_multiplier from the previous time step
        this->initInjMult();

        const auto& summaryState = this->summaryState();
        if (alternative_well_rate_init_) {
            // Update the well rates of well_state_, if only single-phase rates, to
            // have proper multi-phase rates proportional to rates at bhp zero.
            // This is done only for producers, as injectors will only have a single
            // nonzero phase anyway.
            for (const auto& well : well_container_) {
                //const bool zero_target = well->stoppedOrZeroRateTarget(simulator_, this->wellState(), local_deferredLogger);
                if (well->isProducer()){// && !zero_target) {
                    well->initializeWellState(simulator_, this->groupState(), summaryState, this->wellState(), local_deferredLogger);
                }
            }
        }

        for (const auto& well : well_container_) {
            if (well->isVFPActive(local_deferredLogger)){
                well->setPrevSurfaceRates(this->wellState(), this->prevWellState());
            }
        }

        // calculate the well potentials
        try {
            this->updateWellPotentials(reportStepIdx,
                                       /*onlyAfterEvent*/true,
                                       simulator_.vanguard().summaryConfig(),
                                       local_deferredLogger);
        } catch ( std::runtime_error& e ) {
            const std::string msg = "A zero well potential is returned for output purposes. ";
            local_deferredLogger.warning("WELL_POTENTIAL_CALCULATION_FAILED", msg);
        }

        //update guide rates
        const auto& comm = simulator_.vanguard().grid().comm();
        std::vector<Scalar> pot(this->numPhases(), 0.0);
        const Group& fieldGroup = this->schedule().getGroup("FIELD", reportStepIdx);
        WellGroupHelpers<Scalar>::updateGuideRates(fieldGroup,
                                                   this->schedule(),
                                                   summaryState,
                                                   this->phase_usage_,
                                                   reportStepIdx,
                                                   simulationTime,
                                                   this->wellState(),
                                                   this->groupState(),
                                                   comm,
                                                   &this->guideRate_,
                                                   pot,
                                                   local_deferredLogger);
        std::string exc_msg;
        auto exc_type = ExceptionType::NONE;
        // update gpmaint targets
        if (this->schedule_[reportStepIdx].has_gpmaint()) {
            for (const auto& calculator : regionalAveragePressureCalculator_) {
                calculator.second->template defineState<ElementContext>(simulator_);
            }
            const double dt = simulator_.timeStepSize();
            WellGroupHelpers<Scalar>::updateGpMaintTargetForGroups(fieldGroup,
                                                                   this->schedule_,
                                                                   regionalAveragePressureCalculator_,
                                                                   reportStepIdx,
                                                                   dt,
                                                                   this->wellState(),
                                                                   this->groupState());
        }
        try {
            // Compute initial well solution for new wells and injectors that change injection type i.e. WAG.
            for (auto& well : well_container_) {
                const uint64_t effective_events_mask = ScheduleEvents::WELL_STATUS_CHANGE
                        + ScheduleEvents::INJECTION_TYPE_CHANGED
                        + ScheduleEvents::WELL_SWITCHED_INJECTOR_PRODUCER
                        + ScheduleEvents::NEW_WELL;

                const auto& events = this->schedule()[reportStepIdx].wellgroup_events();
                const bool event = this->report_step_starts_ && events.hasEvent(well->name(), effective_events_mask);
                const bool dyn_status_change = this->wellState().well(well->name()).status
                        != this->prevWellState().well(well->name()).status;

                if (event || dyn_status_change) {
                    try {
                        well->updateWellStateWithTarget(simulator_, this->groupState(), this->wellState(), local_deferredLogger);
                        well->calculateExplicitQuantities(simulator_, this->wellState(), local_deferredLogger);
                        well->solveWellEquation(simulator_, this->wellState(), this->groupState(), local_deferredLogger);
                    } catch (const std::exception& e) {
                        const std::string msg = "Compute initial well solution for new well " + well->name() + " failed. Continue with zero initial rates";
                        local_deferredLogger.warning("WELL_INITIAL_SOLVE_FAILED", msg);
                    }
                }
            }
        }
        // Catch clauses for all errors setting exc_type and exc_msg
        OPM_PARALLEL_CATCH_CLAUSE(exc_type, exc_msg);

        if (exc_type != ExceptionType::NONE) {
            const std::string msg = "Compute initial well solution for new wells failed. Continue with zero initial rates";
            local_deferredLogger.warning("WELL_INITIAL_SOLVE_FAILED", msg);
        }

        logAndCheckForExceptionsAndThrow(local_deferredLogger,
                                         exc_type, "beginTimeStep() failed: " + exc_msg, this->terminal_output_, comm);

    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::wellTesting(const int timeStepIdx,
                                            const double simulationTime,
                                            DeferredLogger& deferred_logger)
    {
        for (const std::string& well_name : this->getWellsForTesting(timeStepIdx, simulationTime)) {
            const Well& wellEcl = this->schedule().getWell(well_name, timeStepIdx);
            if (wellEcl.getStatus() == Well::Status::SHUT)
                continue;

            WellInterfacePtr well = createWellForWellTest(well_name, timeStepIdx, deferred_logger);
            // some preparation before the well can be used
            well->init(&this->phase_usage_, depth_, gravity_, B_avg_, true);

            Scalar well_efficiency_factor = wellEcl.getEfficiencyFactor() *
                                            this->wellState().getGlobalEfficiencyScalingFactor(well_name);
            WellGroupHelpers<Scalar>::accumulateGroupEfficiencyFactor(this->schedule().getGroup(wellEcl.groupName(),
                                                                                                timeStepIdx),
                                                                      this->schedule(),
                                                                      timeStepIdx,
                                                                      well_efficiency_factor);

            well->setWellEfficiencyFactor(well_efficiency_factor);
            well->setVFPProperties(this->vfp_properties_.get());
            well->setGuideRate(&this->guideRate_);

            // initialize rates/previous rates to prevent zero fractions in vfp-interpolation
            if (well->isProducer()) {
                well->initializeWellState(simulator_, this->groupState(), this->summaryState(), this->wellState(), deferred_logger);
            }
            if (well->isVFPActive(deferred_logger)) {
                well->setPrevSurfaceRates(this->wellState(), this->prevWellState());
            }

            const auto& network = this->schedule()[timeStepIdx].network();
            if (network.active() && !this->node_pressures_.empty()) {
                if (well->isProducer()) {
                    const auto it = this->node_pressures_.find(well->wellEcl().groupName());
                    if (it != this->node_pressures_.end()) {
                        // The well belongs to a group which has a network nodal pressure,
                        // set the dynamic THP constraint based on the network nodal pressure
                        const Scalar nodal_pressure = it->second;
                        well->setDynamicThpLimit(nodal_pressure);
                    }
                }
            }
            try {
                using GLiftEclWells = typename GasLiftGroupInfo<Scalar>::GLiftEclWells;
                GLiftEclWells ecl_well_map;
                gaslift_.initGliftEclWellMap(well_container_, ecl_well_map);
                well->wellTesting(simulator_,
                                  simulationTime,
                                  this->wellState(),
                                  this->groupState(),
                                  this->wellTestState(),
                                  this->phase_usage_,
                                  ecl_well_map,
                                  this->well_open_times_,
                                  deferred_logger);
            } catch (const std::exception& e) {
                const std::string msg = fmt::format("Exception during testing of well: {}. The well will not open.\n Exception message: {}", wellEcl.name(), e.what());
                deferred_logger.warning("WELL_TESTING_FAILED", msg);
            }
        }
    }





    // called at the end of a report step
    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    endReportStep()
    {
        // Clear the communication data structures for above values.
        for (auto&& pinfo : this->local_parallel_well_info_)
        {
            pinfo.get().clear();
        }
    }





    // called at the end of a report step
    template<typename TypeTag>
    const SimulatorReportSingle&
    BlackoilWellModel<TypeTag>::
    lastReport() const {return last_report_; }





    // called at the end of a time step
    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    timeStepSucceeded(const double simulationTime, const double dt)
    {
        this->closed_this_step_.clear();

        // time step is finished and we are not any more at the beginning of an report step
        this->report_step_starts_ = false;
        const int reportStepIdx = simulator_.episodeIndex();

        DeferredLogger local_deferredLogger;
        for (const auto& well : well_container_) {
            if (getPropValue<TypeTag, Properties::EnablePolymerMW>() && well->isInjector()) {
                well->updateWaterThroughput(dt, this->wellState());
            }
        }
        // update connection transmissibility factor and d factor (if applicable) in the wellstate
        for (const auto& well : well_container_) {
            well->updateConnectionTransmissibilityFactor(simulator_, this->wellState().well(well->indexOfWell()));
            well->updateConnectionDFactor(simulator_, this->wellState().well(well->indexOfWell()));
        }

        if (Indices::waterEnabled) {
            this->updateFiltrationModelsPostStep(dt, FluidSystem::waterPhaseIdx, local_deferredLogger);
        }

        // WINJMULT: At the end of the time step, update the inj_multiplier saved in WellState for later use
        this->updateInjMult(local_deferredLogger);

        // report well switching
        for (const auto& well : well_container_) {
            well->reportWellSwitching(this->wellState().well(well->indexOfWell()), local_deferredLogger);
        }
        // report group switching
        if (this->terminal_output_) {
            this->reportGroupSwitching(local_deferredLogger);
        }

        // update the rate converter with current averages pressures etc in
        rateConverter_->template defineState<ElementContext>(simulator_);

        // calculate the well potentials
        try {
            this->updateWellPotentials(reportStepIdx,
                                       /*onlyAfterEvent*/false,
                                       simulator_.vanguard().summaryConfig(),
                                       local_deferredLogger);
        } catch ( std::runtime_error& e ) {
            const std::string msg = "A zero well potential is returned for output purposes. ";
            local_deferredLogger.warning("WELL_POTENTIAL_CALCULATION_FAILED", msg);
        }

        updateWellTestState(simulationTime, this->wellTestState());

        // check group sales limits at the end of the timestep
        const Group& fieldGroup = this->schedule_.getGroup("FIELD", reportStepIdx);
        this->checkGEconLimits(fieldGroup, simulationTime,
                               simulator_.episodeIndex(), local_deferredLogger);
        this->checkGconsaleLimits(fieldGroup, this->wellState(),
                                  simulator_.episodeIndex(), local_deferredLogger);

        this->calculateProductivityIndexValues(local_deferredLogger);

        this->commitWGState();

        const Opm::Parallel::Communication& comm = grid().comm();
        DeferredLogger global_deferredLogger = gatherDeferredLogger(local_deferredLogger, comm);
        if (this->terminal_output_) {
            global_deferredLogger.logMessages();
        }

        //reporting output temperatures
        this->computeWellTemperature();
    }


    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    computeTotalRatesForDof(RateVector& rate,
                            unsigned elemIdx) const
    {
        rate = 0;

        if (!is_cell_perforated_[elemIdx]) {
            return;
        }

        for (const auto& well : well_container_)
            well->addCellRates(rate, elemIdx);
    }


    template<typename TypeTag>
    template <class Context>
    void
    BlackoilWellModel<TypeTag>::
    computeTotalRatesForDof(RateVector& rate,
                            const Context& context,
                            unsigned spaceIdx,
                            unsigned timeIdx) const
    {
        rate = 0;
        int elemIdx = context.globalSpaceIndex(spaceIdx, timeIdx);

        if (!is_cell_perforated_[elemIdx]) {
            return;
        }

        for (const auto& well : well_container_)
            well->addCellRates(rate, elemIdx);
    }



    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    initializeWellState(const int timeStepIdx)
    {
        const auto pressIx = []()
        {
            if (Indices::oilEnabled)   { return FluidSystem::oilPhaseIdx;   }
            if (Indices::waterEnabled) { return FluidSystem::waterPhaseIdx; }

            return FluidSystem::gasPhaseIdx;
        }();

        auto cellPressures = std::vector<Scalar>(this->local_num_cells_, Scalar{0});
        auto cellTemperatures = std::vector<Scalar>(this->local_num_cells_, Scalar{0});

        auto elemCtx = ElementContext { this->simulator_ };
        const auto& gridView = this->simulator_.vanguard().gridView();

        OPM_BEGIN_PARALLEL_TRY_CATCH();
        for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

            const auto ix = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& fs = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0).fluidState();

            cellPressures[ix] = fs.pressure(pressIx).value();
            cellTemperatures[ix] = fs.temperature(0).value();
        }
        OPM_END_PARALLEL_TRY_CATCH("BlackoilWellModel::initializeWellState() failed: ",
                                   this->simulator_.vanguard().grid().comm());

        this->wellState().init(cellPressures, cellTemperatures, this->schedule(), this->wells_ecl_,
                               this->local_parallel_well_info_, timeStepIdx,
                               &this->prevWellState(), this->well_perf_data_,
                               this->summaryState(), simulator_.vanguard().enableDistributedWells());
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    createWellContainer(const int report_step)
    {
        DeferredLogger local_deferredLogger;

        const int nw = this->numLocalWells();

        well_container_.clear();

        if (nw > 0) {
            well_container_.reserve(nw);

            const auto& wmatcher = this->schedule().wellMatcher(report_step);
            const auto& wcycle = this->schedule()[report_step].wcycle.get();

            // First loop and check for status changes. This is necessary
            // as wcycle needs the updated open/close times.
            std::for_each(this->wells_ecl_.begin(), this->wells_ecl_.end(),
                          [this, &wg_events = this->report_step_start_events_](const auto& well_ecl)
                          {
                              if (!well_ecl.hasConnections()) {
                                  // No connections in this well.  Nothing to do.
                                  return;
                              }

                              constexpr auto events_mask = ScheduleEvents::WELL_STATUS_CHANGE |
                                                           ScheduleEvents::REQUEST_OPEN_WELL;
                              const bool well_status_change =
                                  this->report_step_starts_ &&
                                  wg_events.hasEvent(well_ecl.name(), events_mask);
                              if (well_status_change) {
                                  if (well_ecl.getStatus() == WellStatus::OPEN) {
                                      this->well_open_times_.insert_or_assign(well_ecl.name(),
                                                                              this->simulator_.time());
                                      this->well_close_times_.erase(well_ecl.name());
                                  } else if (well_ecl.getStatus() == WellStatus::SHUT) {
                                      this->well_close_times_.insert_or_assign(well_ecl.name(),
                                                                               this->simulator_.time());
                                      this->well_open_times_.erase(well_ecl.name());
                                  }
                              }
                          });

            // Grab wcycle states. This needs to run before the schedule gets processed
            const auto cycle_states = wcycle.wellStatus(this->simulator_.time(),
                                                         wmatcher,
                                                         this->well_open_times_,
                                                         this->well_close_times_);

            for (int w = 0; w < nw; ++w) {
                const Well& well_ecl = this->wells_ecl_[w];

                if (!well_ecl.hasConnections()) {
                    // No connections in this well.  Nothing to do.
                    continue;
                }

                const std::string& well_name = well_ecl.name();
                const auto well_status = this->schedule()
                    .getWell(well_name, report_step).getStatus();

                if ((well_ecl.getStatus() == Well::Status::SHUT) ||
                    (well_status          == Well::Status::SHUT))
                {
                    // Due to ACTIONX the well might have been closed behind our back.
                    if (well_ecl.getStatus() != Well::Status::SHUT) {
                        this->closed_this_step_.insert(well_name);
                        this->wellState().shutWell(w);
                    }

                    this->well_open_times_.erase(well_name);
                    this->well_close_times_.erase(well_name);
                    continue;
                }

                // A new WCON keywords can re-open a well that was closed/shut due to Physical limit
                if (this->wellTestState().well_is_closed(well_name)) {
                    // The well was shut this timestep, we are most likely retrying
                    // a timestep without the well in question, after it caused
                    // repeated timestep cuts. It should therefore not be opened,
                    // even if it was new or received new targets this report step.
                    const bool closed_this_step = (this->wellTestState().lastTestTime(well_name) == simulator_.time());
                    // TODO: more checking here, to make sure this standard more specific and complete
                    // maybe there is some WCON keywords will not open the well
                    auto& events = this->wellState().well(w).events;
                    if (events.hasEvent(ScheduleEvents::REQUEST_OPEN_WELL)) {
                        if (!closed_this_step) {
                            this->wellTestState().open_well(well_name);
                            this->wellTestState().open_completions(well_name);
                            this->well_open_times_.insert_or_assign(well_name,
                                                                    this->simulator_.time());
                            this->well_close_times_.erase(well_name);
                        }
                        events.clearEvent(ScheduleEvents::REQUEST_OPEN_WELL);
                    }
                }

                // TODO: should we do this for all kinds of closing reasons?
                // something like wellTestState().hasWell(well_name)?
                bool wellIsStopped = false;
                if (this->wellTestState().well_is_closed(well_name))
                {
                    if (well_ecl.getAutomaticShutIn()) {
                        // shut wells are not added to the well container
                        this->wellState().shutWell(w);
                        this->well_close_times_.erase(well_name);
                        this->well_open_times_.erase(well_name);
                        continue;
                    } else {
                        if (!well_ecl.getAllowCrossFlow()) {
                            // stopped wells where cross flow is not allowed
                            // are not added to the well container
                            this->wellState().shutWell(w);
                            this->well_close_times_.erase(well_name);
                            this->well_open_times_.erase(well_name);
                            continue;
                        }
                        // stopped wells are added to the container but marked as stopped
                        this->wellState().stopWell(w);
                        wellIsStopped = true;
                    }
                }

                // shut wells with zero rante constraints and disallowing
                if (!well_ecl.getAllowCrossFlow()) {
                    const bool any_zero_rate_constraint = well_ecl.isProducer()
                        ? well_ecl.productionControls(this->summaryState_).anyZeroRateConstraint()
                        : well_ecl.injectionControls(this->summaryState_).anyZeroRateConstraint();
                    if (any_zero_rate_constraint) {
                        // Treat as shut, do not add to container.
                        local_deferredLogger.debug(fmt::format("  Well {} gets shut due to having zero rate constraint and disallowing crossflow ", well_ecl.name()) );
                        this->wellState().shutWell(w);
                        this->well_close_times_.erase(well_name);
                        this->well_open_times_.erase(well_name);
                        continue;
                    }
                }

                if (well_status == Well::Status::STOP) {
                    this->wellState().stopWell(w);
                    this->well_close_times_.erase(well_name);
                    this->well_open_times_.erase(well_name);
                    wellIsStopped = true;
                }

                if (!wcycle.empty()) {
                    const auto it = cycle_states.find(well_name);
                    if (it != cycle_states.end()) {
                        if (!it->second) {
                            this->wellState().shutWell(w);
                            continue;
                        } else {
                            this->wellState().openWell(w);
                        }
                    }
                }

                well_container_.emplace_back(this->createWellPointer(w, report_step));

                if (wellIsStopped) {
                    well_container_.back()->stopWell();
                    this->well_close_times_.erase(well_name);
                    this->well_open_times_.erase(well_name);
                }
            }

            if (!wcycle.empty()) {
                const auto schedule_open =
                    [&wg_events = this->report_step_start_events_](const std::string& name)
                    {
                        return wg_events.hasEvent(name, ScheduleEvents::REQUEST_OPEN_WELL);
                    };
                for (const auto& [wname, wscale] : wcycle.efficiencyScale(this->simulator_.time(),
                                                                          this->simulator_.timeStepSize(),
                                                                          wmatcher,
                                                                          this->well_open_times_,
                                                                          schedule_open))
                {
                    this->wellState().updateEfficiencyScalingFactor(wname, wscale);
                    this->schedule_.add_event(ScheduleEvents::WELLGROUP_EFFICIENCY_UPDATE, report_step);
                }
            }
        }

        // Collect log messages and print.

        const Opm::Parallel::Communication& comm = grid().comm();
        DeferredLogger global_deferredLogger = gatherDeferredLogger(local_deferredLogger, comm);
        if (this->terminal_output_) {
            global_deferredLogger.logMessages();
        }

        this->well_container_generic_.clear();
        for (auto& w : well_container_)
            this->well_container_generic_.push_back(w.get());

        const auto& network = this->schedule()[report_step].network();
        if (network.active() && !this->node_pressures_.empty()) {
            for (auto& well: this->well_container_generic_) {
                // Producers only, since we so far only support the
                // "extended" network model (properties defined by
                // BRANPROP and NODEPROP) which only applies to producers.
                if (well->isProducer()) {
                    const auto it = this->node_pressures_.find(well->wellEcl().groupName());
                    if (it != this->node_pressures_.end()) {
                        // The well belongs to a group which has a network nodal pressure,
                        // set the dynamic THP constraint based on the network nodal pressure
                        const Scalar nodal_pressure = it->second;
                        well->setDynamicThpLimit(nodal_pressure);
                    }
                }
            }
        }

        this->wbp_.registerOpenWellsForWBPCalculation();
    }





    template <typename TypeTag>
    typename BlackoilWellModel<TypeTag>::WellInterfacePtr
    BlackoilWellModel<TypeTag>::
    createWellPointer(const int wellID, const int report_step) const
    {
        const auto is_multiseg = this->wells_ecl_[wellID].isMultiSegment();

        if (! (this->param_.use_multisegment_well_ && is_multiseg)) {
            return this->template createTypedWellPointer<StandardWell<TypeTag>>(wellID, report_step);
        }
        else {
            return this->template createTypedWellPointer<MultisegmentWell<TypeTag>>(wellID, report_step);
        }
    }





    template <typename TypeTag>
    template <typename WellType>
    std::unique_ptr<WellType>
    BlackoilWellModel<TypeTag>::
    createTypedWellPointer(const int wellID, const int time_step) const
    {
        // Use the pvtRegionIdx from the top cell
        const auto& perf_data = this->well_perf_data_[wellID];

        // Cater for case where local part might have no perforations.
        const auto pvtreg = perf_data.empty()
            ? 0 : this->pvt_region_idx_[perf_data.front().cell_index];

        const auto& parallel_well_info = this->local_parallel_well_info_[wellID].get();
        const auto global_pvtreg = parallel_well_info.broadcastFirstPerforationValue(pvtreg);

        return std::make_unique<WellType>(this->wells_ecl_[wellID],
                                          parallel_well_info,
                                          time_step,
                                          this->param_,
                                          *this->rateConverter_,
                                          global_pvtreg,
                                          this->numComponents(),
                                          this->numPhases(),
                                          wellID,
                                          perf_data);
    }





    template<typename TypeTag>
    typename BlackoilWellModel<TypeTag>::WellInterfacePtr
    BlackoilWellModel<TypeTag>::
    createWellForWellTest(const std::string& well_name,
                          const int report_step,
                          DeferredLogger& deferred_logger) const
    {
        // Finding the location of the well in wells_ecl
        const auto it = std::find_if(this->wells_ecl_.begin(),
                                     this->wells_ecl_.end(),
                                     [&well_name](const auto& w)
                                     { return well_name == w.name(); });
        // It should be able to find in wells_ecl.
        if (it == this->wells_ecl_.end()) {
            OPM_DEFLOG_THROW(std::logic_error,
                             fmt::format("Could not find well {} in wells_ecl ", well_name),
                             deferred_logger);
        }

        const int pos = static_cast<int>(std::distance(this->wells_ecl_.begin(), it));
        return this->createWellPointer(pos, report_step);
    }



    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    doPreStepNetworkRebalance(DeferredLogger& deferred_logger)
    {
        OPM_TIMEFUNCTION();
        const double dt = this->simulator_.timeStepSize();
        // TODO: should we also have the group and network backed-up here in case the solution did not get converged?
        auto& well_state = this->wellState();

        const bool changed_well_group = updateWellControlsAndNetwork(true, dt, deferred_logger);
        assembleWellEqWithoutIteration(dt, deferred_logger);
        const bool converged = this->getWellConvergence(this->B_avg_, true).converged() && !changed_well_group;

        OPM_BEGIN_PARALLEL_TRY_CATCH();
        for (auto& well : this->well_container_) {
            well->solveEqAndUpdateWellState(simulator_, well_state, deferred_logger);
        }
        OPM_END_PARALLEL_TRY_CATCH("BlackoilWellModel::doPreStepNetworkRebalance() failed: ",
                                    this->simulator_.vanguard().grid().comm());

        if (!converged) {
            const std::string msg = fmt::format("Initial (pre-step) network balance did not converge.");
            deferred_logger.warning(msg);
        }
    }




    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    assemble(const int iterationIdx,
             const double dt)
    {
        OPM_TIMEFUNCTION();
        DeferredLogger local_deferredLogger;

        if constexpr (BlackoilWellModelGasLift<TypeTag>::glift_debug) {
            if (gaslift_.terminalOutput()) {
                const std::string msg =
                    fmt::format("assemble() : iteration {}" , iterationIdx);
                gaslift_.gliftDebug(msg, local_deferredLogger);
            }
        }
        last_report_ = SimulatorReportSingle();
        Dune::Timer perfTimer;
        perfTimer.start();
        this->closed_offending_wells_.clear();

        {
            const int episodeIdx = simulator_.episodeIndex();
            const auto& network = this->schedule()[episodeIdx].network();
            if (!this->wellsActive() && !network.active()) {
                return;
            }
        }

        if (iterationIdx == 0 && this->wellsActive()) {
	    OPM_TIMEBLOCK(firstIterationAssmble);
            // try-catch is needed here as updateWellControls
            // contains global communication and has either to
            // be reached by all processes or all need to abort
            // before.
            OPM_BEGIN_PARALLEL_TRY_CATCH();
            {
                calculateExplicitQuantities(local_deferredLogger);
                prepareTimeStep(local_deferredLogger);
            }
            OPM_END_PARALLEL_TRY_CATCH_LOG(local_deferredLogger,
                                           "assemble() failed (It=0): ",
                                           this->terminal_output_, grid().comm());
        }

        const bool well_group_control_changed = updateWellControlsAndNetwork(false, dt, local_deferredLogger);

        // even when there is no wells active, the network nodal pressure still need to be updated through updateWellControlsAndNetwork()
        // but there is no need to assemble the well equations
        if ( ! this->wellsActive() ) {
            return;
        }

        assembleWellEqWithoutIteration(dt, local_deferredLogger);

        // if group or well control changes we don't consider the
        // case converged
        last_report_.well_group_control_changed = well_group_control_changed;
        last_report_.assemble_time_well += perfTimer.stop();
    }




    template<typename TypeTag>
    bool
    BlackoilWellModel<TypeTag>::
    updateWellControlsAndNetwork(const bool mandatory_network_balance,
                                 const double dt,
                                 DeferredLogger& local_deferredLogger)
    {
        OPM_TIMEFUNCTION();
        // not necessarily that we always need to update once of the network solutions
        bool do_network_update = true;
        bool well_group_control_changed = false;
        Scalar network_imbalance = 0.0;
        // after certain number of the iterations, we use relaxed tolerance for the network update
        const std::size_t iteration_to_relax = param_.network_max_strict_outer_iterations_;
        // after certain number of the iterations, we terminate
        const std::size_t max_iteration = param_.network_max_outer_iterations_;
        std::size_t network_update_iteration = 0;
        while (do_network_update) {
            if (network_update_iteration >= max_iteration ) {
                // only output to terminal if we at the last newton iterations where we try to balance the network.
                const int episodeIdx = simulator_.episodeIndex();
                const int iterationIdx = simulator_.model().newtonMethod().numIterations();
                if (this->shouldBalanceNetwork(episodeIdx, iterationIdx + 1)) {
                    const std::string msg = fmt::format("Maximum of {:d} network iterations has been used and we stop the update, \n"
                        "and try again after the next Newton iteration (imbalance = {:.2e} bar, ctrl_change = {})",
                        max_iteration, network_imbalance*1.0e-5, well_group_control_changed);
                    local_deferredLogger.debug(msg);
                } else {
                    if (this->terminal_output_) {
                        const std::string msg = fmt::format("Maximum of {:d} network iterations has been used and we stop the update. \n"
                            "The simulator will continue with unconverged network results (imbalance = {:.2e} bar, ctrl_change = {})",
                            max_iteration, network_imbalance*1.0e-5, well_group_control_changed);
                        local_deferredLogger.info(msg);
                    }
                }
                break;
            }
            if (this->terminal_output_ && (network_update_iteration == iteration_to_relax) ) {
                local_deferredLogger.debug("We begin using relaxed tolerance for network update now after " + std::to_string(iteration_to_relax) + " iterations ");
            }
            const bool relax_network_balance = network_update_iteration >= iteration_to_relax;
            // Never optimize gas lift in last iteration, to allow network convergence (unless max_iter < 2)
            const bool optimize_gas_lift = ( (network_update_iteration + 1) < std::max(max_iteration, static_cast<std::size_t>(2)) );
            std::tie(well_group_control_changed, do_network_update, network_imbalance) =
                    updateWellControlsAndNetworkIteration(mandatory_network_balance, relax_network_balance, optimize_gas_lift, dt,local_deferredLogger);
            ++network_update_iteration;
        }
        return well_group_control_changed;
    }




    template<typename TypeTag>
    std::tuple<bool, bool, typename BlackoilWellModel<TypeTag>::Scalar>
    BlackoilWellModel<TypeTag>::
    updateWellControlsAndNetworkIteration(const bool mandatory_network_balance,
                                          const bool relax_network_tolerance,
                                          const bool optimize_gas_lift,
                                          const double dt,
                                          DeferredLogger& local_deferredLogger)
    {
        OPM_TIMEFUNCTION();
        auto [well_group_control_changed, more_inner_network_update, network_imbalance] =
                updateWellControls(mandatory_network_balance,
                                   local_deferredLogger,
                                   relax_network_tolerance);

        bool alq_updated = false;
        OPM_BEGIN_PARALLEL_TRY_CATCH();
        {
            if (optimize_gas_lift) alq_updated = gaslift_.maybeDoGasLiftOptimize(simulator_,
                                                          well_container_,
                                                          this->wellState(),
                                                          this->groupState(),
                                                          local_deferredLogger);

            prepareWellsBeforeAssembling(dt, local_deferredLogger);
        }
        OPM_END_PARALLEL_TRY_CATCH_LOG(local_deferredLogger,
                                       "updateWellControlsAndNetworkIteration() failed: ",
                                       this->terminal_output_, grid().comm());

        // update guide rates
        const int reportStepIdx = simulator_.episodeIndex();
        if (alq_updated || BlackoilWellModelGuideRates(*this).
                              guideRateUpdateIsNeeded(reportStepIdx)) {
            const double simulationTime = simulator_.time();
            const auto& comm = simulator_.vanguard().grid().comm();
            const auto& summaryState = simulator_.vanguard().summaryState();
            std::vector<Scalar> pot(this->numPhases(), 0.0);
            const Group& fieldGroup = this->schedule().getGroup("FIELD", reportStepIdx);
            WellGroupHelpers<Scalar>::updateGuideRates(fieldGroup,
                                                       this->schedule(),
                                                       summaryState,
                                                       this->phase_usage_,
                                                       reportStepIdx,
                                                       simulationTime,
                                                       this->wellState(),
                                                       this->groupState(),
                                                       comm,
                                                       &this->guideRate_,
                                                       pot,
                                                       local_deferredLogger);
        }
        // we need to re-iterate the network when the well group controls changed or gaslift/alq is changed or
        // the inner iterations are did not converge
        const int iterationIdx = simulator_.model().newtonMethod().numIterations();
        const bool more_network_update = this->shouldBalanceNetwork(reportStepIdx, iterationIdx) &&
                    (more_inner_network_update || well_group_control_changed || alq_updated);
        return {well_group_control_changed, more_network_update, network_imbalance};
    }

    // This function is to be used for well groups in an extended network that act as a subsea manifold
    // The wells of such group should have a common THP and total phase rate(s) obeying (if possible)
    // the well group constraint set by GCONPROD
    template <typename TypeTag>
    bool
    BlackoilWellModel<TypeTag>::
    computeWellGroupThp(const double dt, DeferredLogger& local_deferredLogger)
    {
        OPM_TIMEFUNCTION();
        const int reportStepIdx = this->simulator_.episodeIndex();
        const auto& network = this->schedule()[reportStepIdx].network();
        const auto& balance = this->schedule()[reportStepIdx].network_balance();
        const Scalar thp_tolerance = balance.thp_tolerance();

        if (!network.active()) {
            return false;
        }

        auto& well_state = this->wellState();
        auto& group_state = this->groupState();

        bool well_group_thp_updated = false;
        for (const std::string& nodeName : network.node_names()) {
            const bool has_choke = network.node(nodeName).as_choke();
            if (has_choke) {
                const auto& summary_state = this->simulator_.vanguard().summaryState();
                const Group& group = this->schedule().getGroup(nodeName, reportStepIdx);

                const auto pu = this->phase_usage_;
                //TODO: Auto choke combined with RESV control is not supported
                std::vector<Scalar> resv_coeff(pu.num_phases, 1.0);
                Scalar gratTargetFromSales = 0.0;
                if (group_state.has_grat_sales_target(group.name()))
                    gratTargetFromSales = group_state.grat_sales_target(group.name());

                const auto ctrl = group.productionControls(summary_state);
                auto cmode_tmp = ctrl.cmode;
                Scalar target_tmp{0.0};
                bool fld_none = false;
                if (cmode_tmp == Group::ProductionCMode::FLD || cmode_tmp == Group::ProductionCMode::NONE) {
                    fld_none = true;
                    // Target is set for an ancestor group. Target for autochoke group to be 
                    // derived via group guide rates
                    const Scalar efficiencyFactor = 1.0;
                    const Group& parentGroup = this->schedule().getGroup(group.parent(), reportStepIdx);
                    auto target = WellGroupControls<Scalar>::getAutoChokeGroupProductionTargetRate(
                                                            group.name(),
                                                            parentGroup,
                                                            well_state,
                                                            group_state,
                                                            this->schedule(),
                                                            summary_state,
                                                            resv_coeff,
                                                            efficiencyFactor,
                                                            reportStepIdx,
                                                            pu,
                                                            &this->guideRate_,
                                                            local_deferredLogger);
                    target_tmp = target.first;
                    cmode_tmp = target.second;
                }
                const auto cmode = cmode_tmp;
                WGHelpers::TargetCalculator tcalc(cmode, pu, resv_coeff,
                                                  gratTargetFromSales, nodeName, group_state,
                                                  group.has_gpmaint_control(cmode));
                if (!fld_none)
                {
                    // Target is set for the autochoke group itself
                    target_tmp = tcalc.groupTarget(ctrl, local_deferredLogger);
                }

                const Scalar orig_target = target_tmp;

                auto mismatch = [&] (auto group_thp) {
                    Scalar group_rate(0.0);
                    Scalar rate(0.0);
                    for (auto& well : this->well_container_) {
                        std::string well_name = well->name();
                        auto& ws = well_state.well(well_name);
                        if (group.hasWell(well_name)) {
                            well->setDynamicThpLimit(group_thp);
                            const Well& well_ecl = this->wells_ecl_[well->indexOfWell()];
                            const auto inj_controls = Well::InjectionControls(0);
                            const auto prod_controls = well_ecl.productionControls(summary_state);
                            well->iterateWellEqWithSwitching(this->simulator_, dt, inj_controls, prod_controls, well_state, group_state, local_deferredLogger,  false, false);
                            rate = -tcalc.calcModeRateFromRates(ws.surface_rates);
                            group_rate += rate;
                        }
                    }
                    return (group_rate - orig_target)/orig_target;
                };

                const auto upbranch = network.uptree_branch(nodeName);
                const auto it = this->node_pressures_.find((*upbranch).uptree_node());
                const Scalar nodal_pressure = it->second;
                Scalar well_group_thp = nodal_pressure;

                std::optional<Scalar> autochoke_thp;
                if (auto iter = this->well_group_thp_calc_.find(nodeName); iter != this->well_group_thp_calc_.end()) {
                    autochoke_thp = this->well_group_thp_calc_.at(nodeName);
                }

                //Find an initial bracket
                std::array<Scalar, 2> range_initial;
                if (!autochoke_thp.has_value()){
                    Scalar min_thp, max_thp;
                    // Retrieve the terminal pressure of the associated root of the manifold group
                    std::string node_name =  nodeName;
                    while (!network.node(node_name).terminal_pressure().has_value()) {
                        auto branch = network.uptree_branch(node_name).value();
                        node_name = branch.uptree_node();
                    }
                    min_thp = network.node(node_name).terminal_pressure().value();
                    WellBhpThpCalculator<Scalar>::bruteForceBracketCommonTHP(mismatch, min_thp, max_thp);
                    // Narrow down the bracket
                    Scalar low1, high1;
                    std::array<Scalar, 2> range = {Scalar{0.9}*min_thp, Scalar{1.1}*max_thp};
                    std::optional<Scalar> appr_sol;
                    WellBhpThpCalculator<Scalar>::bruteForceBracketCommonTHP(mismatch, range, low1, high1, appr_sol, 0.0, local_deferredLogger);
                    min_thp = low1;
                    max_thp = high1;
                    range_initial = {min_thp, max_thp};
                }

                if (!autochoke_thp.has_value() || autochoke_thp.value() > nodal_pressure) {
                    // The bracket is based on the initial bracket or on a range based on a previous calculated group thp
                    std::array<Scalar, 2> range = autochoke_thp.has_value() ?
                        std::array<Scalar, 2>{Scalar{0.9} * autochoke_thp.value(),
                                              Scalar{1.1} * autochoke_thp.value()} : range_initial;
                    Scalar low, high;
                    std::optional<Scalar> approximate_solution;
                    const Scalar tolerance1 = thp_tolerance;
                    local_deferredLogger.debug("Using brute force search to bracket the group THP");
                    const bool finding_bracket = WellBhpThpCalculator<Scalar>::bruteForceBracketCommonTHP(mismatch, range, low, high, approximate_solution, tolerance1, local_deferredLogger);

                    if (approximate_solution.has_value()) {
                        autochoke_thp = *approximate_solution;
                        local_deferredLogger.debug("Approximate group THP value found: "  + std::to_string(autochoke_thp.value()));
                    } else if (finding_bracket) {
                        const Scalar tolerance2 = thp_tolerance;
                        const int max_iteration_solve = 100;
                        int iteration = 0;
                        autochoke_thp = RegulaFalsiBisection<ThrowOnError>::
                                         solve(mismatch, low, high, max_iteration_solve, tolerance2, iteration);
                        local_deferredLogger.debug(" bracket = [" + std::to_string(low) + ", " + std::to_string(high) + "], " +
                                                   "iteration = " + std::to_string(iteration));
                        local_deferredLogger.debug("Group THP value = " + std::to_string(autochoke_thp.value()));
                    } else {
                        autochoke_thp.reset();
                        local_deferredLogger.debug("Group THP solve failed due to bracketing failure");
                    }
                }
                 if (autochoke_thp.has_value()) {
                    well_group_thp_calc_[nodeName] = autochoke_thp.value();
                    // Note: The node pressure of the auto-choke node is set to well_group_thp in computeNetworkPressures()
                    // and must be larger or equal to the pressure of the uptree node of its branch.
                    well_group_thp = std::max(autochoke_thp.value(), nodal_pressure);
                }

                for (auto& well : this->well_container_) {
                    std::string well_name = well->name();

                    if (well->isInjector() || !well->wellEcl().predictionMode())
                        continue;

                    if (group.hasWell(well_name)) {
                        well->setDynamicThpLimit(well_group_thp);
                    }
                    const auto& ws = this->wellState().well(well->indexOfWell());
                    const bool thp_is_limit = ws.production_cmode == Well::ProducerCMode::THP;
                    if (thp_is_limit) {
                        well->prepareWellBeforeAssembling(this->simulator_, dt, this->wellState(), this->groupState(), local_deferredLogger);
                    }
                }

                // Use the group THP in computeNetworkPressures().
                const auto& current_well_group_thp = group_state.is_autochoke_group(nodeName) ? group_state.well_group_thp(nodeName) : 1e30;
                if (std::abs(current_well_group_thp - well_group_thp) > balance.pressure_tolerance()) {
                    well_group_thp_updated = true;
                    group_state.update_well_group_thp(nodeName, well_group_thp);
                }
            }
        }
        return well_group_thp_updated;
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    assembleWellEq(const double dt, DeferredLogger& deferred_logger)
    {
        OPM_TIMEFUNCTION();
        for (auto& well : well_container_) {
            well->assembleWellEq(simulator_, dt, this->wellState(), this->groupState(), deferred_logger);
        }
    }


    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    prepareWellsBeforeAssembling(const double dt, DeferredLogger& deferred_logger)
    {
        OPM_TIMEFUNCTION();
        for (auto& well : well_container_) {
            well->prepareWellBeforeAssembling(simulator_, dt, this->wellState(), this->groupState(), deferred_logger);
        }
    }


    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    assembleWellEqWithoutIteration(const double dt, DeferredLogger& deferred_logger)
    {
        OPM_TIMEFUNCTION();
        // We make sure that all processes throw in case there is an exception
        // on one of them (WetGasPvt::saturationPressure might throw if not converged)
        OPM_BEGIN_PARALLEL_TRY_CATCH();

        for (auto& well: well_container_) {
            well->assembleWellEqWithoutIteration(simulator_, dt, this->wellState(), this->groupState(),
                                                 deferred_logger);
        }
        OPM_END_PARALLEL_TRY_CATCH_LOG(deferred_logger, "BlackoilWellModel::assembleWellEqWithoutIteration failed: ",
                                       this->terminal_output_, grid().comm());

    }

#if COMPILE_GPU_BRIDGE
    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    getWellContributions(WellContributions<Scalar>& wellContribs) const
    {
        // prepare for StandardWells
        wellContribs.setBlockSize(StandardWell<TypeTag>::Indices::numEq, StandardWell<TypeTag>::numStaticWellEq);

        for(unsigned int i = 0; i < well_container_.size(); i++){
            auto& well = well_container_[i];
            std::shared_ptr<StandardWell<TypeTag> > derived = std::dynamic_pointer_cast<StandardWell<TypeTag> >(well);
            if (derived) {
                wellContribs.addNumBlocks(derived->linSys().getNumBlocks());
            }
        }

        // allocate memory for data from StandardWells
        wellContribs.alloc();

        for(unsigned int i = 0; i < well_container_.size(); i++){
            auto& well = well_container_[i];
            // maybe WellInterface could implement addWellContribution()
            auto derived_std = std::dynamic_pointer_cast<StandardWell<TypeTag>>(well);
            if (derived_std) {
                derived_std->linSys().extract(derived_std->numStaticWellEq, wellContribs);
            } else {
                auto derived_ms = std::dynamic_pointer_cast<MultisegmentWell<TypeTag> >(well);
                if (derived_ms) {
                    derived_ms->linSys().extract(wellContribs);
                } else {
                    OpmLog::warning("Warning unknown type of well");
                }
            }
        }
    }
#endif

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    addWellContributions(SparseMatrixAdapter& jacobian) const
    {
        for ( const auto& well: well_container_ ) {
            well->addWellContributions(jacobian);
        }
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    addWellPressureEquations(PressureMatrix& jacobian,
                             const BVector& weights,
                             const bool use_well_weights) const
    {
        int nw = this->numLocalWellsEnd();
        int rdofs = local_num_cells_;
        for ( int i = 0; i < nw; i++ ) {
            int wdof = rdofs + i;
            jacobian[wdof][wdof] = 1.0;// better scaling ?
        }

        for (const auto& well : well_container_) {
            well->addWellPressureEquations(jacobian,
                                           weights,
                                           pressureVarIndex,
                                           use_well_weights,
                                           this->wellState());
        }
    }

    template <typename TypeTag>
    void BlackoilWellModel<TypeTag>::
    addReservoirSourceTerms(GlobalEqVector& residual,
                            const std::vector<typename SparseMatrixAdapter::MatrixBlock*>& diagMatAddress) const
    {
        // NB this loop may write multiple times to the same element
        // if a cell is perforated by more than one well, so it should
        // not be OpenMP-parallelized.
        for (const auto& well : well_container_) {
            if (!well->isOperableAndSolvable() && !well->wellIsStopped()) {
                continue;
            }
            const auto& cells = well->cells();
            const auto& rates = well->connectionRates();
            for (unsigned perfIdx = 0; perfIdx < rates.size(); ++perfIdx) {
                unsigned cellIdx = cells[perfIdx];
                auto rate = rates[perfIdx];
                rate *= -1.0;
                VectorBlockType res(0.0);
                using MatrixBlockType = typename SparseMatrixAdapter::MatrixBlock;
                MatrixBlockType bMat(0.0);
                simulator_.model().linearizer().setResAndJacobi(res, bMat, rate);
                residual[cellIdx] += res;
                *diagMatAddress[cellIdx] += bMat;
            }
        }
    }


    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    addWellPressureEquationsStruct(PressureMatrix& jacobian) const
    {
        int nw =  this->numLocalWellsEnd();
        int rdofs = local_num_cells_;
        for (int i = 0; i < nw; ++i) {
            int wdof = rdofs + i;
            jacobian.entry(wdof,wdof) = 1.0;// better scaling ?
        }
        const auto wellconnections = this->getMaxWellConnections();
        for (int i = 0; i < nw; ++i) {
            const auto& perfcells = wellconnections[i];
            for (int perfcell : perfcells) {
                int wdof = rdofs + i;
                jacobian.entry(wdof, perfcell) = 0.0;
                jacobian.entry(perfcell, wdof) = 0.0;
            }
        }
    }


    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    recoverWellSolutionAndUpdateWellState(const BVector& x)
    {
        DeferredLogger local_deferredLogger;
        OPM_BEGIN_PARALLEL_TRY_CATCH();
        {
            for (const auto& well : well_container_) {
                const auto& cells = well->cells();
                x_local_.resize(cells.size());

                for (size_t i = 0; i < cells.size(); ++i) {
                    x_local_[i] = x[cells[i]];
                }
                well->recoverWellSolutionAndUpdateWellState(simulator_, x_local_,
                                                            this->wellState(), local_deferredLogger);
            }
        }
        OPM_END_PARALLEL_TRY_CATCH_LOG(local_deferredLogger,
                                       "recoverWellSolutionAndUpdateWellState() failed: ",
                                       this->terminal_output_, simulator_.vanguard().grid().comm());
    }


    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    recoverWellSolutionAndUpdateWellStateDomain(const BVector& x, const int domainIdx)
    {
        if (!nldd_) {
            OPM_THROW(std::logic_error, "Attempt to call NLDD method without a NLDD solver");
        }

        return nldd_->recoverWellSolutionAndUpdateWellState(x, domainIdx);
    }


    template<typename TypeTag>
    ConvergenceReport
    BlackoilWellModel<TypeTag>::
    getWellConvergence(const std::vector<Scalar>& B_avg, bool checkWellGroupControls) const
    {

        DeferredLogger local_deferredLogger;
        // Get global (from all processes) convergence report.
        ConvergenceReport local_report;
        const int iterationIdx = simulator_.model().newtonMethod().numIterations();
        for (const auto& well : well_container_) {
            if (well->isOperableAndSolvable() || well->wellIsStopped()) {
                local_report += well->getWellConvergence(
                        simulator_, this->wellState(), B_avg, local_deferredLogger,
                        iterationIdx > param_.strict_outer_iter_wells_);
            } else {
                ConvergenceReport report;
                using CR = ConvergenceReport;
                report.setWellFailed({CR::WellFailure::Type::Unsolvable, CR::Severity::Normal, -1, well->name()});
                local_report += report;
            }
        }

        const Opm::Parallel::Communication comm = grid().comm();
        DeferredLogger global_deferredLogger = gatherDeferredLogger(local_deferredLogger, comm);
        ConvergenceReport report = gatherConvergenceReport(local_report, comm);

        // the well_group_control_changed info is already communicated
        if (checkWellGroupControls) {
            report.setWellGroupTargetsViolated(this->lastReport().well_group_control_changed);
        }

        if (this->terminal_output_) {
            global_deferredLogger.logMessages();

            // Log debug messages for NaN or too large residuals.
            for (const auto& f : report.wellFailures()) {
                if (f.severity() == ConvergenceReport::Severity::NotANumber) {
                        OpmLog::debug("NaN residual found with phase " + std::to_string(f.phase()) + " for well " + f.wellName());
                } else if (f.severity() == ConvergenceReport::Severity::TooLarge) {
                        OpmLog::debug("Too large residual found with phase " + std::to_string(f.phase()) + " for well " + f.wellName());
                }
            }
        }
        return report;
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    calculateExplicitQuantities(DeferredLogger& deferred_logger) const
    {
        // TODO: checking isOperableAndSolvable() ?
        for (auto& well : well_container_) {
            well->calculateExplicitQuantities(simulator_, this->wellState(), deferred_logger);
        }
    }





    template<typename TypeTag>
    std::tuple<bool, bool, typename BlackoilWellModel<TypeTag>::Scalar>
    BlackoilWellModel<TypeTag>::
    updateWellControls(const bool mandatory_network_balance,
                       DeferredLogger& deferred_logger,
                       const bool relax_network_tolerance)
    {
        OPM_TIMEFUNCTION();
        const int episodeIdx = simulator_.episodeIndex();
        const auto& network = this->schedule()[episodeIdx].network();
        if (!this->wellsActive() && !network.active()) {
            return {false, false, 0.0};
        }

        const int iterationIdx = simulator_.model().newtonMethod().numIterations();
        const auto& comm = simulator_.vanguard().grid().comm();
        this->updateAndCommunicateGroupData(episodeIdx, iterationIdx, param_.nupcol_group_rate_tolerance_, deferred_logger);

        // network related
        Scalar network_imbalance = 0.0;
        bool more_network_update = false;
        if (this->shouldBalanceNetwork(episodeIdx, iterationIdx) || mandatory_network_balance) {
            OPM_TIMEBLOCK(BalanceNetowork);
            const double dt = this->simulator_.timeStepSize();
            // Calculate common THP for subsea manifold well group (item 3 of NODEPROP set to YES)
            bool well_group_thp_updated = computeWellGroupThp(dt, deferred_logger);
            const int max_number_of_sub_iterations = param_.network_max_sub_iterations_;
            const Scalar network_pressure_update_damping_factor = param_.network_pressure_update_damping_factor_;
            const Scalar network_max_pressure_update = param_.network_max_pressure_update_in_bars_ * unit::barsa;
            bool more_network_sub_update = false;
            for (int i = 0; i < max_number_of_sub_iterations; i++) {
                const auto local_network_imbalance = this->updateNetworkPressures(episodeIdx, network_pressure_update_damping_factor, network_max_pressure_update);
                network_imbalance = comm.max(local_network_imbalance);
                const auto& balance = this->schedule()[episodeIdx].network_balance();
                constexpr Scalar relaxation_factor = 10.0;
                const Scalar tolerance = relax_network_tolerance ? relaxation_factor * balance.pressure_tolerance() : balance.pressure_tolerance();
                more_network_sub_update = this->networkActive() && network_imbalance > tolerance;
                if (!more_network_sub_update)
                    break;

                for (const auto& well : well_container_) {
                    if (well->isInjector() || !well->wellEcl().predictionMode())
                         continue;

                    const auto it = this->node_pressures_.find(well->wellEcl().groupName());
                    if (it != this->node_pressures_.end()) {
                        const auto& ws = this->wellState().well(well->indexOfWell());
                        const bool thp_is_limit = ws.production_cmode == Well::ProducerCMode::THP;
                        if (thp_is_limit) {
                            well->prepareWellBeforeAssembling(this->simulator_, dt, this->wellState(), this->groupState(), deferred_logger);
                        }
                    }
                }
                this->updateAndCommunicateGroupData(episodeIdx, iterationIdx, param_.nupcol_group_rate_tolerance_, deferred_logger);
            }
            more_network_update = more_network_sub_update || well_group_thp_updated;
        }

        bool changed_well_group = false;
        // Check group individual constraints.
        const Group& fieldGroup = this->schedule().getGroup("FIELD", episodeIdx);
        changed_well_group = updateGroupControls(fieldGroup, deferred_logger, episodeIdx, iterationIdx);

        // Check wells' group constraints and communicate.
        bool changed_well_to_group = false;
        {
            OPM_TIMEBLOCK(UpdateWellControls);
            // For MS Wells a linear solve is performed below and the matrix might be singular.
            // We need to communicate the exception thrown to the others and rethrow.
            OPM_BEGIN_PARALLEL_TRY_CATCH()
                for (const auto& well : well_container_) {
                    const auto mode = WellInterface<TypeTag>::IndividualOrGroup::Group;
                    const bool changed_well = well->updateWellControl(simulator_, mode, this->wellState(), this->groupState(), deferred_logger);
                    if (changed_well) {
                        changed_well_to_group = changed_well || changed_well_to_group;
                    }
                }
            OPM_END_PARALLEL_TRY_CATCH("BlackoilWellModel: updating well controls failed: ",
                                       simulator_.gridView().comm());
        }

        changed_well_to_group = comm.sum(static_cast<int>(changed_well_to_group));
        if (changed_well_to_group) {
            updateAndCommunicate(episodeIdx, iterationIdx, deferred_logger);
            changed_well_group = true;
        }

        // Check individual well constraints and communicate.
        bool changed_well_individual = false;
        {
            // For MS Wells a linear solve is performed below and the matrix might be singular.
            // We need to communicate the exception thrown to the others and rethrow.
            OPM_BEGIN_PARALLEL_TRY_CATCH()
                for (const auto& well : well_container_) {
                    const auto mode = WellInterface<TypeTag>::IndividualOrGroup::Individual;
                    const bool changed_well = well->updateWellControl(simulator_, mode, this->wellState(), this->groupState(), deferred_logger);
                    if (changed_well) {
                        changed_well_individual = changed_well || changed_well_individual;
                    }
                }
            OPM_END_PARALLEL_TRY_CATCH("BlackoilWellModel: updating well controls failed: ",
                                       simulator_.gridView().comm());
        }

        changed_well_individual = comm.sum(static_cast<int>(changed_well_individual));
        if (changed_well_individual) {
            updateAndCommunicate(episodeIdx, iterationIdx, deferred_logger);
            changed_well_group = true;
        }

        // update wsolvent fraction for REIN wells
        this->updateWsolvent(fieldGroup, episodeIdx,  this->nupcolWellState());

        return { changed_well_group, more_network_update, network_imbalance };
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    updateAndCommunicate(const int reportStepIdx,
                         const int iterationIdx,
                         DeferredLogger& deferred_logger)
    {
        this->updateAndCommunicateGroupData(reportStepIdx,
                                            iterationIdx,
                                            param_.nupcol_group_rate_tolerance_,
                                            deferred_logger);

        // updateWellStateWithTarget might throw for multisegment wells hence we
        // have a parallel try catch here to thrown on all processes.
        OPM_BEGIN_PARALLEL_TRY_CATCH()
        // if a well or group change control it affects all wells that are under the same group
        for (const auto& well : well_container_) {
            // We only want to update wells under group-control here
            const auto& ws = this->wellState().well(well->indexOfWell());
            if (ws.production_cmode ==  Well::ProducerCMode::GRUP ||
                ws.injection_cmode == Well::InjectorCMode::GRUP)
            {
                well->updateWellStateWithTarget(simulator_, this->groupState(),
                                                this->wellState(), deferred_logger);
            }
        }
        OPM_END_PARALLEL_TRY_CATCH("BlackoilWellModel::updateAndCommunicate failed: ",
                                   simulator_.gridView().comm())
        this->updateAndCommunicateGroupData(reportStepIdx,
                                            iterationIdx,
                                            param_.nupcol_group_rate_tolerance_,
                                            deferred_logger);
    }

    template<typename TypeTag>
    bool
    BlackoilWellModel<TypeTag>::
    updateGroupControls(const Group& group,
                        DeferredLogger& deferred_logger,
                        const int reportStepIdx,
                        const int iterationIdx)
    {
        OPM_TIMEFUNCTION();
        bool changed = false;
        // restrict the number of group switches but only after nupcol iterations.
        const int nupcol = this->schedule()[reportStepIdx].nupcol();
        const int max_number_of_group_switches = iterationIdx < nupcol ? 9999 : param_.max_number_of_group_switches_;
        bool changed_hc = this->checkGroupHigherConstraints( group, deferred_logger, reportStepIdx, max_number_of_group_switches);
        if (changed_hc) {
            changed = true;
            updateAndCommunicate(reportStepIdx, iterationIdx, deferred_logger);
        }

        bool changed_individual =
            BlackoilWellModelConstraints(*this).
                updateGroupIndividualControl(group,
                                             reportStepIdx,
                                             max_number_of_group_switches,
                                             this->switched_inj_groups_,
                                             this->switched_prod_groups_,
                                             this->closed_offending_wells_,
                                             this->groupState(),
                                             this->wellState(),
                                             deferred_logger);

        if (changed_individual) {
            changed = true;
            updateAndCommunicate(reportStepIdx, iterationIdx, deferred_logger);
        }
        // call recursively down the group hierarchy
        for (const std::string& groupName : group.groups()) {
            bool changed_this = updateGroupControls( this->schedule().getGroup(groupName, reportStepIdx), deferred_logger, reportStepIdx,iterationIdx);
            changed = changed || changed_this;
        }
        return changed;
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    updateWellTestState(const double& simulationTime, WellTestState& wellTestState) const
    {
        OPM_TIMEFUNCTION();
        DeferredLogger local_deferredLogger;
        for (const auto& well : well_container_) {
            const auto& wname = well->name();
            const auto wasClosed = wellTestState.well_is_closed(wname);
            well->checkWellOperability(simulator_,
                                       this->wellState(),
                                       local_deferredLogger);
            const bool under_zero_target =
                well->wellUnderZeroGroupRateTarget(this->simulator_,
                                                   this->wellState(),
                                                   local_deferredLogger);
            well->updateWellTestState(this->wellState().well(wname),
                                      simulationTime,
                                      /*writeMessageToOPMLog=*/ true,
                                      under_zero_target,
                                      wellTestState,
                                      local_deferredLogger);

            if (!wasClosed && wellTestState.well_is_closed(wname)) {
                this->closed_this_step_.insert(wname);
            }
        }

        const Opm::Parallel::Communication comm = grid().comm();
        DeferredLogger global_deferredLogger =
            gatherDeferredLogger(local_deferredLogger, comm);

        for (const auto& [group_name, to] : this->closed_offending_wells_) {
            if (this->hasOpenLocalWell(to.second) &&
                !this->wasDynamicallyShutThisTimeStep(to.second))
            {
                wellTestState.close_well(to.second,
                                         WellTestConfig::Reason::GROUP,
                                         simulationTime);
                this->updateClosedWellsThisStep(to.second);
                const std::string msg =
                    fmt::format("Procedure on exceeding {} limit is WELL for group {}. "
                                "Well {} is {}.",
                                to.first,
                                group_name,
                                to.second,
                                "shut");
                global_deferredLogger.info(msg);
            }
        }

        if (this->terminal_output_) {
            global_deferredLogger.logMessages();
        }
    }


    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::computePotentials(const std::size_t widx,
                                                  const WellState<Scalar>& well_state_copy,
                                                  std::string& exc_msg,
                                                  ExceptionType::ExcEnum& exc_type,
                                                  DeferredLogger& deferred_logger)
    {
        OPM_TIMEFUNCTION();
        const int np = this->numPhases();
        std::vector<Scalar> potentials;
        const auto& well = well_container_[widx];
        std::string cur_exc_msg;
        auto cur_exc_type = ExceptionType::NONE;
        try {
            well->computeWellPotentials(simulator_, well_state_copy, potentials, deferred_logger);
        }
        // catch all possible exception and store type and message.
        OPM_PARALLEL_CATCH_CLAUSE(cur_exc_type, cur_exc_msg);
        if (cur_exc_type != ExceptionType::NONE) {
            exc_msg += fmt::format("\nFor well {}: {}", well->name(), cur_exc_msg);
        }
        exc_type = std::max(exc_type, cur_exc_type);
        // Store it in the well state
        // potentials is resized and set to zero in the beginning of well->ComputeWellPotentials
        // and updated only if sucessfull. i.e. the potentials are zero for exceptions
        auto& ws = this->wellState().well(well->indexOfWell());
        for (int p = 0; p < np; ++p) {
            // make sure the potentials are positive
            ws.well_potentials[p] = std::max(Scalar{0.0}, potentials[p]);
        }
    }



    template <typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    calculateProductivityIndexValues(DeferredLogger& deferred_logger)
    {
        for (const auto& wellPtr : this->well_container_) {
            this->calculateProductivityIndexValues(wellPtr.get(), deferred_logger);
        }
    }





    template <typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    calculateProductivityIndexValuesShutWells(const int reportStepIdx,
                                              DeferredLogger& deferred_logger)
    {
        // For the purpose of computing PI/II values, it is sufficient to
        // construct StandardWell instances only.  We don't need to form
        // well objects that honour the 'isMultisegment()' flag of the
        // corresponding "this->wells_ecl_[shutWell]".

        for (const auto& shutWell : this->local_shut_wells_) {
            if (!this->wells_ecl_[shutWell].hasConnections()) {
                // No connections in this well.  Nothing to do.
                continue;
            }

            auto wellPtr = this->template createTypedWellPointer
                <StandardWell<TypeTag>>(shutWell, reportStepIdx);

            wellPtr->init(&this->phase_usage_, this->depth_, this->gravity_, this->B_avg_, true);

            this->calculateProductivityIndexValues(wellPtr.get(), deferred_logger);
        }
    }





    template <typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    calculateProductivityIndexValues(const WellInterface<TypeTag>* wellPtr,
                                     DeferredLogger& deferred_logger)
    {
        wellPtr->updateProductivityIndex(this->simulator_,
                                         this->prod_index_calc_[wellPtr->indexOfWell()],
                                         this->wellState(),
                                         deferred_logger);
    }



    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    prepareTimeStep(DeferredLogger& deferred_logger)
    {
        // Check if there is a network with active prediction wells at this time step.
        const auto episodeIdx = simulator_.episodeIndex();
        this->updateNetworkActiveState(episodeIdx);

        // Rebalance the network initially if any wells in the network have status changes
        // (Need to check this before clearing events)
        const bool do_prestep_network_rebalance = param_.pre_solve_network_ && this->needPreStepNetworkRebalance(episodeIdx);

        for (const auto& well : well_container_) {
            auto& events = this->wellState().well(well->indexOfWell()).events;
            if (events.hasEvent(WellState<Scalar>::event_mask)) {
                well->updateWellStateWithTarget(simulator_, this->groupState(), this->wellState(), deferred_logger);
                well->updatePrimaryVariables(simulator_, this->wellState(), deferred_logger);
                // There is no new well control change input within a report step,
                // so next time step, the well does not consider to have effective events anymore.
                events.clearEvent(WellState<Scalar>::event_mask);
            }
            // these events only work for the first time step within the report step
            if (events.hasEvent(ScheduleEvents::REQUEST_OPEN_WELL)) {
                events.clearEvent(ScheduleEvents::REQUEST_OPEN_WELL);
            }
            // solve the well equation initially to improve the initial solution of the well model
            if (param_.solve_welleq_initially_ && well->isOperableAndSolvable()) {
                try {
                    well->solveWellEquation(simulator_, this->wellState(), this->groupState(), deferred_logger);
                } catch (const std::exception& e) {
                    const std::string msg = "Compute initial well solution for " + well->name() + " initially failed. Continue with the previous rates";
                    deferred_logger.warning("WELL_INITIAL_SOLVE_FAILED", msg);
                }
            }
            // If we're using local well solves that include control switches, they also update
            // operability, so reset before main iterations begin
            well->resetWellOperability();
        }
        updatePrimaryVariables(deferred_logger);

        // Actually do the pre-step network rebalance, using the updated well states and initial solutions
        if (do_prestep_network_rebalance) doPreStepNetworkRebalance(deferred_logger);
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    updateAverageFormationFactor()
    {
        std::vector< Scalar > B_avg(numComponents(), Scalar() );
        const auto& grid = simulator_.vanguard().grid();
        const auto& gridView = grid.leafGridView();
        ElementContext elemCtx(simulator_);

        OPM_BEGIN_PARALLEL_TRY_CATCH();
        for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

            const auto& intQuants = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& fs = intQuants.fluidState();

            for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx)
            {
                if (!FluidSystem::phaseIsActive(phaseIdx)) {
                    continue;
                }

                const unsigned compIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
                auto& B  = B_avg[ compIdx ];

                B += 1 / fs.invB(phaseIdx).value();
            }
            if constexpr (has_solvent_) {
                auto& B  = B_avg[solventSaturationIdx];
                B += 1 / intQuants.solventInverseFormationVolumeFactor().value();
            }
        }
        OPM_END_PARALLEL_TRY_CATCH("BlackoilWellModel::updateAverageFormationFactor() failed: ", grid.comm())

        // compute global average
        grid.comm().sum(B_avg.data(), B_avg.size());
        for (auto& bval : B_avg)
        {
            bval /= global_num_cells_;
        }
        B_avg_ = B_avg;
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    updatePrimaryVariables(DeferredLogger& deferred_logger)
    {
        for (const auto& well : well_container_) {
            well->updatePrimaryVariables(simulator_, this->wellState(), deferred_logger);
        }
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::extractLegacyCellPvtRegionIndex_()
    {
        const auto& grid = simulator_.vanguard().grid();
        const auto& eclProblem = simulator_.problem();
        const unsigned numCells = grid.size(/*codim=*/0);

        this->pvt_region_idx_.resize(numCells);
        for (unsigned cellIdx = 0; cellIdx < numCells; ++cellIdx) {
            this->pvt_region_idx_[cellIdx] =
                eclProblem.pvtRegionIndex(cellIdx);
        }
    }

    // The number of components in the model.
    template<typename TypeTag>
    int
    BlackoilWellModel<TypeTag>::numComponents() const
    {
        // The numComponents here does not reflect the actual number of the components in the system.
        // It more or less reflects the number of mass conservation equations for the well equations.
        // For example, in the current formulation, we do not have the polymer conservation equation
        // in the well equations. As a result, for an oil-water-polymer system, this function will return 2.
        // In some way, it makes this function appear to be confusing from its name, and we need
        // to revisit/revise this function again when extending the variants of system that flow can simulate.
        int numComp = this->numPhases() < 3 ? this->numPhases() : FluidSystem::numComponents;
        if constexpr (has_solvent_) {
            numComp++;
        }
        return numComp;
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::extractLegacyDepth_()
    {
        const auto& eclProblem = simulator_.problem();
        depth_.resize(local_num_cells_);
        for (unsigned cellIdx = 0; cellIdx < local_num_cells_; ++cellIdx) {
            depth_[cellIdx] = eclProblem.dofCenterDepth(cellIdx);
        }
    }

    template<typename TypeTag>
    typename BlackoilWellModel<TypeTag>::WellInterfacePtr
    BlackoilWellModel<TypeTag>::
    getWell(const std::string& well_name) const
    {
        // finding the iterator of the well in wells_ecl
        auto well = std::find_if(well_container_.begin(),
                                     well_container_.end(),
                                     [&well_name](const WellInterfacePtr& elem)->bool {
                                         return elem->name() == well_name;
                                     });

        assert(well != well_container_.end());

        return *well;
    }

    template <typename TypeTag>
    int
    BlackoilWellModel<TypeTag>::
    reportStepIndex() const
    {
        return std::max(this->simulator_.episodeIndex(), 0);
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    calcResvCoeff(const int fipnum,
                  const int pvtreg,
                  const std::vector<Scalar>& production_rates,
                  std::vector<Scalar>& resv_coeff)
    {
        rateConverter_->calcCoeff(fipnum, pvtreg, production_rates, resv_coeff);
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    calcInjResvCoeff(const int fipnum,
                     const int pvtreg,
                     std::vector<Scalar>& resv_coeff)
    {
        rateConverter_->calcInjCoeff(fipnum, pvtreg, resv_coeff);
    }


    template <typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    computeWellTemperature()
    {
        if constexpr (has_energy_) {
            int np = this->numPhases();
            Scalar cellInternalEnergy;
            Scalar cellBinv;
            Scalar cellDensity;
            Scalar perfPhaseRate;
            const int nw = this->numLocalWells();
            for (auto wellID = 0*nw; wellID < nw; ++wellID) {
                const Well& well = this->wells_ecl_[wellID];
                auto& ws = this->wellState().well(wellID);
                if (well.isInjector()) {
                    if (ws.status != WellStatus::STOP) {
                        this->wellState().well(wellID).temperature = well.inj_temperature();
                        continue;
                    }
                }

                std::array<Scalar,2> weighted{0.0,0.0};
                auto& [weighted_temperature, total_weight] = weighted;

                auto& well_info = this->local_parallel_well_info_[wellID].get();
                auto& perf_data = ws.perf_data;
                auto& perf_phase_rate = perf_data.phase_rates;

                using int_type = decltype(this->well_perf_data_[wellID].size());
                for (int_type perf = 0, end_perf = this->well_perf_data_[wellID].size(); perf < end_perf; ++perf) {
                    const int cell_idx = this->well_perf_data_[wellID][perf].cell_index;
                    const auto& intQuants = simulator_.model().intensiveQuantities(cell_idx, /*timeIdx=*/0);
                    const auto& fs = intQuants.fluidState();

                    // we on only have one temperature pr cell any phaseIdx will do
                    Scalar cellTemperatures = fs.temperature(/*phaseIdx*/0).value();

                    Scalar weight_factor = 0.0;
                    for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
                        if (!FluidSystem::phaseIsActive(phaseIdx)) {
                            continue;
                        }
                        cellInternalEnergy = fs.enthalpy(phaseIdx).value() -
                                             fs.pressure(phaseIdx).value() / fs.density(phaseIdx).value();
                        cellBinv = fs.invB(phaseIdx).value();
                        cellDensity = fs.density(phaseIdx).value();
                        perfPhaseRate = perf_phase_rate[perf*np + phaseIdx];
                        weight_factor += cellDensity * perfPhaseRate / cellBinv * cellInternalEnergy / cellTemperatures;
                    }
                    weight_factor = std::abs(weight_factor) + 1e-13;
                    total_weight += weight_factor;
                    weighted_temperature += weight_factor * cellTemperatures;
                }
                well_info.communication().sum(weighted.data(), 2);
                this->wellState().well(wellID).temperature = weighted_temperature / total_weight;
            }
        }
    }
} // namespace Opm

#endif
