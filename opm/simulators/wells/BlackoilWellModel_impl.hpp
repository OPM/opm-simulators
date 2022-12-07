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

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/core/props/phaseUsageFromDeck.hpp>
#include <opm/grid/utility/cartesianToCompressed.hpp>
#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <opm/simulators/wells/BlackoilWellModelConstraints.hpp>
#include <opm/simulators/wells/VFPProperties.hpp>
#include <opm/simulators/utils/MPIPacker.hpp>
#include <opm/simulators/linalg/bda/WellContributions.hpp>

#if HAVE_MPI
#include <ebos/eclmpiserializer.hh>
#endif

#include <algorithm>
#include <utility>

#include <fmt/format.h>

namespace Opm {
    template<typename TypeTag>
    BlackoilWellModel<TypeTag>::
    BlackoilWellModel(Simulator& ebosSimulator, const PhaseUsage& phase_usage)
        : BlackoilWellModelGeneric(ebosSimulator.vanguard().schedule(),
                                   ebosSimulator.vanguard().summaryState(),
                                   ebosSimulator.vanguard().eclState(),
                                   phase_usage,
                                   ebosSimulator.gridView().comm())
        , ebosSimulator_(ebosSimulator)
    {
        terminal_output_ = ((ebosSimulator.gridView().comm().rank() == 0) &&
                           EWOMS_GET_PARAM(TypeTag, bool, EnableTerminalOutput));

        local_num_cells_ = ebosSimulator_.gridView().size(0);

        // Number of cells the global grid view
        global_num_cells_ = ebosSimulator_.vanguard().globalNumCells();

        {
            auto& parallel_wells = ebosSimulator.vanguard().parallelWells();

            this->parallel_well_info_.reserve(parallel_wells.size());
            for( const auto& name_bool : parallel_wells) {
                this->parallel_well_info_.emplace_back(name_bool, grid().comm());
            }
        }

        this->alternative_well_rate_init_ =
            EWOMS_GET_PARAM(TypeTag, bool, AlternativeWellRateInit);
    }

    template<typename TypeTag>
    BlackoilWellModel<TypeTag>::
    BlackoilWellModel(Simulator& ebosSimulator) :
        BlackoilWellModel(ebosSimulator, phaseUsageFromDeck(ebosSimulator.vanguard().eclState()))
    {}


    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    init()
    {
        extractLegacyCellPvtRegionIndex_();
        extractLegacyDepth_();

        gravity_ = ebosSimulator_.problem().gravity()[2];

        initial_step_ = true;

        // add the eWoms auxiliary module for the wells to the list
        ebosSimulator_.model().addAuxiliaryModule(this);

        is_cell_perforated_.resize(local_num_cells_, false);
    }


    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    initWellContainer(const int reportStepIdx)
    {
        const uint64_t effective_events_mask = ScheduleEvents::WELL_STATUS_CHANGE
                        + ScheduleEvents::NEW_WELL;
        const auto& events = schedule()[reportStepIdx].wellgroup_events();
        for (auto& wellPtr : this->well_container_) {
            const bool well_opened_this_step = report_step_starts_ && events.hasEvent(wellPtr->name(), effective_events_mask);
            wellPtr->init(&this->phase_usage_, this->depth_, this->gravity_,
                          this->local_num_cells_, this->B_avg_, well_opened_this_step);
        }
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    addNeighbors(std::vector<NeighborSet>& neighbors) const
    {
        if (!param_.matrix_add_well_contributions_) {
            return;
        }

        // Create cartesian to compressed mapping
        const auto& schedule_wells = schedule().getWellsatEnd();

        // initialize the additional cell connections introduced by wells.
        for (const auto& well : schedule_wells)
        {
            std::vector<int> wellCells;
            // All possible connections of the well
            const auto& connectionSet = well.getConnections();
            wellCells.reserve(connectionSet.size());

            for ( size_t c=0; c < connectionSet.size(); c++ )
            {
                const auto& connection = connectionSet.get(c);
                int compressed_idx = compressedIndexForInterior(connection.global_index());
                if ( compressed_idx >= 0 ) { // Ignore connections in inactive/remote cells.
                    wellCells.push_back(compressed_idx);
                }
            }

            for (int cellIdx : wellCells) {
                neighbors[cellIdx].insert(wellCells.begin(),
                                          wellCells.end());
            }
        }
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    linearize(SparseMatrixAdapter& jacobian, GlobalEqVector& res)
    {
        if (!param_.matrix_add_well_contributions_)
        {
            OPM_BEGIN_PARALLEL_TRY_CATCH();
            {
                // if the well contributions are not supposed to be included explicitly in
                // the matrix, we only apply the vector part of the Schur complement here.
                for (const auto& well: well_container_) {
                    // r = r - duneC_^T * invDuneD_ * resWell_
                    well->apply(res);
                }
            }
            OPM_END_PARALLEL_TRY_CATCH("BlackoilWellModel::linearize failed: ",
                                       ebosSimulator_.gridView().comm());
            return;
        }

        for (const auto& well: well_container_) {
            well->addWellContributions(jacobian);

            // applying the well residual to reservoir residuals
            // r = r - duneC_^T * invDuneD_ * resWell_
            well->apply(res);
        }
    }


    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    beginReportStep(const int timeStepIdx)
    {
        DeferredLogger local_deferredLogger;

        report_step_starts_ = true;

        const Grid& grid = ebosSimulator_.vanguard().grid();
        const auto& summaryState = ebosSimulator_.vanguard().summaryState();
        // Make wells_ecl_ contain only this partition's wells.
        wells_ecl_ = getLocalWells(timeStepIdx);
        this->local_parallel_well_info_ = createLocalParallelWellInfo(wells_ecl_);

        // at least initializeWellState might be throw
        // exception in opm-material (UniformTabulated2DFunction.hpp)
        // playing it safe by extending the scope a bit.
        OPM_BEGIN_PARALLEL_TRY_CATCH();
        {

            // The well state initialize bhp with the cell pressure in the top cell.
            // We must therefore provide it with updated cell pressures
            this->initializeWellPerfData();
            this->initializeWellState(timeStepIdx, summaryState);

            // handling MS well related
            if (param_.use_multisegment_well_&& anyMSWellOpenLocal()) { // if we use MultisegmentWell model
                this->wellState().initWellStateMSWell(wells_ecl_, &this->prevWellState());
            }

            const Group& fieldGroup = schedule().getGroup("FIELD", timeStepIdx);
            WellGroupHelpers::setCmodeGroup(fieldGroup, schedule(), summaryState, timeStepIdx, this->wellState(), this->groupState());

            // Compute reservoir volumes for RESV controls.
            rateConverter_.reset(new RateConverterType (phase_usage_,
                                                        std::vector<int>(local_num_cells_, 0)));
            rateConverter_->template defineState<ElementContext>(ebosSimulator_);

            // Compute regional average pressures used by gpmaint
            if (schedule_[timeStepIdx].has_gpmaint()) {
                const auto& fp = this->eclState_.fieldProps();
                const auto& fipnum = fp.get_int("FIPNUM");
                regionalAveragePressureCalculator_.reset(new AverageRegionalPressureType (phase_usage_,fipnum));
            }

            {
                const auto& sched_state = this->schedule()[timeStepIdx];
                // update VFP properties
                vfp_properties_.reset(new VFPProperties( sched_state.vfpinj(), sched_state.vfpprod(), this->prevWellState()));
                this->initializeWellProdIndCalculators();
                if (sched_state.events().hasEvent(ScheduleEvents::Events::WELL_PRODUCTIVITY_INDEX)) {
                    this->runWellPIScaling(timeStepIdx, local_deferredLogger);
                }
            }
        }
        OPM_END_PARALLEL_TRY_CATCH_LOG(local_deferredLogger, "beginReportStep() failed: ",
                                       terminal_output_, grid.comm());
        // Store the current well state, to be able to recover in the case of failed iterations
        this->commitWGState();
    }


    // called at the beginning of a time step
    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    beginTimeStep()
    {
        updatePerforationIntensiveQuantities();
        updateAverageFormationFactor();
        DeferredLogger local_deferredLogger;
        switched_prod_groups_.clear();
        switched_inj_groups_.clear();

        this->resetWGState();
        const int reportStepIdx = ebosSimulator_.episodeIndex();
        updateAndCommunicateGroupData(reportStepIdx,
                                      ebosSimulator_.model().newtonMethod().numIterations());
        this->wellState().gliftTimeStepInit();
        const double simulationTime = ebosSimulator_.time();
        OPM_BEGIN_PARALLEL_TRY_CATCH();
        {
            // test wells
            wellTesting(reportStepIdx, simulationTime, local_deferredLogger);

            // create the well container
            createWellContainer(reportStepIdx);

            // Wells are active if they are active wells on at least one process.
            const Grid& grid = ebosSimulator_.vanguard().grid();
            wells_active_ = !this->well_container_.empty();
            wells_active_ = grid.comm().max(wells_active_);

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
            calculateEfficiencyFactors(reportStepIdx);

            if constexpr (has_polymer_)
            {
                if (PolymerModule::hasPlyshlog() || getPropValue<TypeTag, Properties::EnablePolymerMW>() ) {
                    setRepRadiusPerfLength();
                }
            }

        }
        OPM_END_PARALLEL_TRY_CATCH_LOG(local_deferredLogger, "beginTimeStep() failed: ",
                                        terminal_output_, ebosSimulator_.vanguard().grid().comm());

        for (auto& well : well_container_) {
            well->setVFPProperties(vfp_properties_.get());
            well->setGuideRate(&guideRate_);
        }

        // Close completions due to economical reasons
        for (auto& well : well_container_) {
            well->closeCompletions(wellTestState());
        }

        if (alternative_well_rate_init_) {
            // Update the well rates of well_state_, if only single-phase rates, to
            // have proper multi-phase rates proportional to rates at bhp zero.
            // This is done only for producers, as injectors will only have a single
            // nonzero phase anyway.
            for (auto& well : well_container_) {
                if (well->isProducer()) {
                    well->updateWellStateRates(ebosSimulator_, this->wellState(), local_deferredLogger);
                }
            }
        }

        // calculate the well potentials
        try {
            updateWellPotentials(reportStepIdx,
                                 /*onlyAfterEvent*/true,
                                 ebosSimulator_.vanguard().summaryConfig(),
                                 local_deferredLogger);
        } catch ( std::runtime_error& e ) {
            const std::string msg = "A zero well potential is returned for output purposes. ";
            local_deferredLogger.warning("WELL_POTENTIAL_CALCULATION_FAILED", msg);
        }

        //update guide rates
        const auto& comm = ebosSimulator_.vanguard().grid().comm();
        const auto& summaryState = ebosSimulator_.vanguard().summaryState();
        std::vector<double> pot(numPhases(), 0.0);
        const Group& fieldGroup = schedule().getGroup("FIELD", reportStepIdx);
        WellGroupHelpers::updateGuideRates(fieldGroup, schedule(), summaryState, this->phase_usage_, reportStepIdx, simulationTime,
                                           this->wellState(), this->groupState(), comm, &this->guideRate_, pot, local_deferredLogger);
        std::string exc_msg;
        auto exc_type = ExceptionType::NONE;
        // update gpmaint targets
        if (schedule_[reportStepIdx].has_gpmaint()) {
            regionalAveragePressureCalculator_->template defineState<ElementContext>(ebosSimulator_);
            const double dt = ebosSimulator_.timeStepSize();
            WellGroupHelpers::updateGpMaintTargetForGroups(fieldGroup,
                                                           schedule_, *regionalAveragePressureCalculator_, reportStepIdx, dt, this->wellState(), this->groupState());
        }
        try {
            // Compute initial well solution for new wells and injectors that change injection type i.e. WAG.
            for (auto& well : well_container_) {
                const uint64_t effective_events_mask = ScheduleEvents::WELL_STATUS_CHANGE
                        + ScheduleEvents::INJECTION_TYPE_CHANGED
                        + ScheduleEvents::WELL_SWITCHED_INJECTOR_PRODUCER
                        + ScheduleEvents::NEW_WELL;

                const auto& events = schedule()[reportStepIdx].wellgroup_events();
                const bool event = report_step_starts_ && events.hasEvent(well->name(), effective_events_mask);
                const bool dyn_status_change = this->wellState().well(well->name()).status
                        != this->prevWellState().well(well->name()).status;

                if (event || dyn_status_change) {
                    try {
                        well->updateWellStateWithTarget(ebosSimulator_, this->groupState(), this->wellState(), local_deferredLogger);
                        well->calculateExplicitQuantities(ebosSimulator_, this->wellState(), local_deferredLogger);
                        well->solveWellEquation(ebosSimulator_, this->wellState(), this->groupState(), local_deferredLogger);
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
            exc_type, "beginTimeStep() failed: " + exc_msg, terminal_output_, comm);

    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::wellTesting(const int timeStepIdx,
                                            const double simulationTime,
                                            DeferredLogger& deferred_logger)
    {
        const auto& wtest_config = schedule()[timeStepIdx].wtest_config();
        if (!wtest_config.empty()) { // there is a WTEST request
            const std::vector<std::string> wellsForTesting = wellTestState()
                .test_wells(wtest_config, simulationTime);
            for (const std::string& well_name : wellsForTesting) {

                const Well& wellEcl = schedule().getWell(well_name, timeStepIdx);
                if (wellEcl.getStatus() == Well::Status::SHUT)
                    continue;

                WellInterfacePtr well = createWellForWellTest(well_name, timeStepIdx, deferred_logger);
                // some preparation before the well can be used
                well->init(&phase_usage_, depth_, gravity_, local_num_cells_, B_avg_, true);

                double well_efficiency_factor = wellEcl.getEfficiencyFactor();
                WellGroupHelpers::accumulateGroupEfficiencyFactor(schedule().getGroup(wellEcl.groupName(), timeStepIdx),
                                                                  schedule(), timeStepIdx, well_efficiency_factor);

                well->setWellEfficiencyFactor(well_efficiency_factor);
                well->setVFPProperties(vfp_properties_.get());
                well->setGuideRate(&guideRate_);

                well->wellTesting(ebosSimulator_, simulationTime, this->wellState(), this->groupState(), wellTestState(), deferred_logger);
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
    timeStepSucceeded(const double& simulationTime, const double dt)
    {
        this->closed_this_step_.clear();

        // time step is finished and we are not any more at the beginning of an report step
        report_step_starts_ = false;
        const int reportStepIdx = ebosSimulator_.episodeIndex();

        DeferredLogger local_deferredLogger;
        for (const auto& well : well_container_) {
            if (getPropValue<TypeTag, Properties::EnablePolymerMW>() && well->isInjector()) {
                well->updateWaterThroughput(dt, this->wellState());
            }
        }
        // report well switching
        for (const auto& well : well_container_) {
            well->reportWellSwitching(this->wellState().well(well->indexOfWell()), local_deferredLogger);
        }
        // report group switching
        if (terminal_output_) {

            for (const auto& [name, to] : switched_prod_groups_) {
                const Group::ProductionCMode& oldControl = this->prevWGState().group_state.production_control(name);
                std::string from = Group::ProductionCMode2String(oldControl);
                if (to != from) {
                    std::string msg = "    Production Group " + name
                    + " control mode changed from ";
                    msg += from;
                    msg += " to " + to;
                    local_deferredLogger.info(msg);
                }
            }
            for (const auto& [key, to] : switched_inj_groups_) {
                const std::string& name = key.first;
                const Opm::Phase& phase = key.second;

                const Group::InjectionCMode& oldControl = this->prevWGState().group_state.injection_control(name, phase);
                std::string from = Group::InjectionCMode2String(oldControl);
                if (to != from) {
                    std::string msg = "    Injection Group " + name
                    + " control mode changed from ";
                    msg += from;
                    msg += " to " + to;
                    local_deferredLogger.info(msg);
                }
            }
        }

        // update the rate converter with current averages pressures etc in
        rateConverter_->template defineState<ElementContext>(ebosSimulator_);

        // calculate the well potentials
        try {
            updateWellPotentials(reportStepIdx,
                                 /*onlyAfterEvent*/false,
                                 ebosSimulator_.vanguard().summaryConfig(),
                                 local_deferredLogger);
        } catch ( std::runtime_error& e ) {
            const std::string msg = "A zero well potential is returned for output purposes. ";
            local_deferredLogger.warning("WELL_POTENTIAL_CALCULATION_FAILED", msg);
        }

        updateWellTestState(simulationTime, wellTestState());

        // check group sales limits at the end of the timestep
        const Group& fieldGroup = schedule_.getGroup("FIELD", reportStepIdx);
        checkGconsaleLimits(fieldGroup, this->wellState(),
                            ebosSimulator_.episodeIndex(), local_deferredLogger);

        this->calculateProductivityIndexValues(local_deferredLogger);

        this->commitWGState();
 
        const Opm::Parallel::Communication& comm = grid().comm();
        DeferredLogger global_deferredLogger = gatherDeferredLogger(local_deferredLogger, comm);
        if (terminal_output_) {
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

        if (!is_cell_perforated_[elemIdx])
            return;

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

        if (!is_cell_perforated_[elemIdx])
            return;

        for (const auto& well : well_container_)
            well->addCellRates(rate, elemIdx);
    }



    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    initializeWellState(const int           timeStepIdx,
                        const SummaryState& summaryState)
    {
        std::vector<double> cellPressures(this->local_num_cells_, 0.0);
        ElementContext elemCtx(ebosSimulator_);

        const auto& gridView = ebosSimulator_.vanguard().gridView();

        OPM_BEGIN_PARALLEL_TRY_CATCH();
        for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

            const auto& fs = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0).fluidState();
            // copy of get perfpressure in Standard well except for value
            double& perf_pressure = cellPressures[elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0)];
            if (Indices::oilEnabled) {
                perf_pressure = fs.pressure(FluidSystem::oilPhaseIdx).value();
            } else if (Indices::waterEnabled) {
                perf_pressure = fs.pressure(FluidSystem::waterPhaseIdx).value();
            } else {
                perf_pressure = fs.pressure(FluidSystem::gasPhaseIdx).value();
            }
        }
        OPM_END_PARALLEL_TRY_CATCH("BlackoilWellModel::initializeWellState() failed: ", ebosSimulator_.vanguard().grid().comm());

        this->wellState().init(cellPressures, schedule(), wells_ecl_, local_parallel_well_info_, timeStepIdx,
                               &this->prevWellState(), well_perf_data_,
                               summaryState);
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    createWellContainer(const int time_step)
    {
        DeferredLogger local_deferredLogger;

        const int nw = numLocalWells();

        well_container_.clear();

        if (nw > 0) {
            well_container_.reserve(nw);

            for (int w = 0; w < nw; ++w) {
                const Well& well_ecl = wells_ecl_[w];

                if (well_ecl.getConnections().empty()) {
                    // No connections in this well.  Nothing to do.
                    continue;
                }

                const std::string& well_name = well_ecl.name();
                const auto well_status = this->schedule()
                    .getWell(well_name, time_step).getStatus();

                if ((well_ecl.getStatus() == Well::Status::SHUT) ||
                    (well_status          == Well::Status::SHUT))
                {
                    // Due to ACTIONX the well might have been closed behind our back.
                    if (well_ecl.getStatus() != Well::Status::SHUT) {
                        this->closed_this_step_.insert(well_name);
                        this->wellState().shutWell(w);
                    }

                    continue;
                }

                // A new WCON keywords can re-open a well that was closed/shut due to Physical limit
                if (this->wellTestState().well_is_closed(well_name)) {
                    // TODO: more checking here, to make sure this standard more specific and complete
                    // maybe there is some WCON keywords will not open the well
                    auto& events = this->wellState().well(w).events;
                    if (events.hasEvent(WellState::event_mask)) {
                        if (wellTestState().lastTestTime(well_name) == ebosSimulator_.time()) {
                            // The well was shut this timestep, we are most likely retrying
                            // a timestep without the well in question, after it caused
                            // repeated timestep cuts. It should therefore not be opened,
                            // even if it was new or received new targets this report step.
                            events.clearEvent(WellState::event_mask);
                        } else {
                            wellTestState().open_well(well_name);
                            wellTestState().open_completions(well_name);
                        }
                    }
                }

                // TODO: should we do this for all kinds of closing reasons?
                // something like wellTestState().hasWell(well_name)?
                bool wellIsStopped = false;
                if (wellTestState().well_is_closed(well_name))
                {
                    if (well_ecl.getAutomaticShutIn()) {
                        // shut wells are not added to the well container
                        this->wellState().shutWell(w);
                        continue;
                    } else {
                        if (!well_ecl.getAllowCrossFlow()) {
                            // stopped wells where cross flow is not allowed
                            // are not added to the well container
                            this->wellState().shutWell(w);
                            continue;
                        }
                        // stopped wells are added to the container but marked as stopped
                        this->wellState().stopWell(w);
                        wellIsStopped = true;
                    }
                }

                // If a production well disallows crossflow and its
                // (prediction type) rate control is zero, then it is effectively shut.
                if (!well_ecl.getAllowCrossFlow() && well_ecl.isProducer() && well_ecl.predictionMode()) {
                    const auto& summaryState = ebosSimulator_.vanguard().summaryState();
                    const auto prod_controls = well_ecl.productionControls(summaryState);

                    auto is_zero = [](const double x)
                    {
                        return std::isfinite(x) && !std::isnormal(x);
                    };

                    bool zero_rate_control = false;
                    switch (prod_controls.cmode) {
                    case Well::ProducerCMode::ORAT:
                        zero_rate_control = is_zero(prod_controls.oil_rate);
                        break;

                    case Well::ProducerCMode::WRAT:
                        zero_rate_control = is_zero(prod_controls.water_rate);
                        break;

                    case Well::ProducerCMode::GRAT:
                        zero_rate_control = is_zero(prod_controls.gas_rate);
                        break;

                    case Well::ProducerCMode::LRAT:
                        zero_rate_control = is_zero(prod_controls.liquid_rate);
                        break;

                    case Well::ProducerCMode::RESV:
                        zero_rate_control = is_zero(prod_controls.resv_rate);
                        break;
                    default:
                        // Might still have zero rate controls, but is pressure controlled.
                        zero_rate_control = false;
                        break;
                    }

                    if (zero_rate_control) {
                        // Treat as shut, do not add to container.
                        local_deferredLogger.info("  Well shut due to zero rate control and disallowing crossflow: " + well_ecl.name());
                        this->wellState().shutWell(w);
                        continue;
                    }
                }

                if (well_status == Well::Status::STOP) {
                    this->wellState().stopWell(w);
                    wellIsStopped = true;
                }

                well_container_.emplace_back(this->createWellPointer(w, time_step));

                if (wellIsStopped)
                    well_container_.back()->stopWell();
            }
        }

        // Collect log messages and print.
        
        const Opm::Parallel::Communication& comm = grid().comm();
        DeferredLogger global_deferredLogger = gatherDeferredLogger(local_deferredLogger, comm);
        if (terminal_output_) {
            global_deferredLogger.logMessages();
        }

        well_container_generic_.clear();
        for (auto& w : well_container_)
          well_container_generic_.push_back(w.get());
    }





    template <typename TypeTag>
    typename BlackoilWellModel<TypeTag>::WellInterfacePtr
    BlackoilWellModel<TypeTag>::
    createWellPointer(const int wellID, const int time_step) const
    {
        const auto is_multiseg = this->wells_ecl_[wellID].isMultiSegment();

        if (! (this->param_.use_multisegment_well_ && is_multiseg)) {
            return this->template createTypedWellPointer<StandardWell<TypeTag>>(wellID, time_step);
        }
        else {
            return this->template createTypedWellPointer<MultisegmentWell<TypeTag>>(wellID, time_step);
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
            ? 0 : pvt_region_idx_[perf_data.front().cell_index];

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
        const int nw_wells_ecl = wells_ecl_.size();
        int index_well_ecl = 0;
        for (; index_well_ecl < nw_wells_ecl; ++index_well_ecl) {
          if (well_name == wells_ecl_[index_well_ecl].name()) {
                break;
            }
        }
        // It should be able to find in wells_ecl.
        if (index_well_ecl == nw_wells_ecl) {
            OPM_DEFLOG_THROW(std::logic_error, "Could not find well " << well_name << " in wells_ecl ", deferred_logger);
        }

        return this->createWellPointer(index_well_ecl, report_step);
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    assemble(const int iterationIdx,
             const double dt)
    {

        DeferredLogger local_deferredLogger;
        if (this->glift_debug) {
            const std::string msg = fmt::format(
                "assemble() : iteration {}" , iterationIdx);
            gliftDebug(msg, local_deferredLogger);
        }
        last_report_ = SimulatorReportSingle();
        Dune::Timer perfTimer;
        perfTimer.start();

        if ( ! wellsActive() ) {
            return;
        }

        updatePerforationIntensiveQuantities();

        if (iterationIdx == 0) {
            // try-catch is needed here as updateWellControls
            // contains global communication and has either to
            // be reached by all processes or all need to abort
            // before.
            OPM_BEGIN_PARALLEL_TRY_CATCH();
            {
                calculateExplicitQuantities(local_deferredLogger);
                prepareTimeStep(local_deferredLogger);
            }
            OPM_END_PARALLEL_TRY_CATCH_LOG(local_deferredLogger, "assemble() failed (It=0): ",
                                               terminal_output_, grid().comm());
        }

        const bool well_group_control_changed = assembleImpl(iterationIdx, dt, 0, local_deferredLogger);

        // if group or well control changes we don't consider the
        // case converged
        last_report_.well_group_control_changed = well_group_control_changed;
        last_report_.assemble_time_well += perfTimer.stop();
    }





    template<typename TypeTag>
    bool
    BlackoilWellModel<TypeTag>::
    assembleImpl(const int iterationIdx,
                 const double dt,
                 const std::size_t recursion_level,
                 DeferredLogger& local_deferredLogger)
    {

        auto [well_group_control_changed, network_changed, network_imbalance] = updateWellControls(local_deferredLogger);

        bool alq_updated = false;
        OPM_BEGIN_PARALLEL_TRY_CATCH();
        {
            // Set the well primary variables based on the value of well solutions
            initPrimaryVariablesEvaluation();

            alq_updated = maybeDoGasLiftOptimize(local_deferredLogger);
            assembleWellEq(dt, local_deferredLogger);
        }
        OPM_END_PARALLEL_TRY_CATCH_LOG(local_deferredLogger, "assemble() failed: ",
                                       terminal_output_, grid().comm());

        //update guide rates
        const int reportStepIdx = ebosSimulator_.episodeIndex();
        if (alq_updated || BlackoilWellModelGuideRates(*this).
                           guideRateUpdateIsNeeded(reportStepIdx)) {
            const double simulationTime = ebosSimulator_.time();
            const auto& comm = ebosSimulator_.vanguard().grid().comm();
            const auto& summaryState = ebosSimulator_.vanguard().summaryState();
            std::vector<double> pot(numPhases(), 0.0);
            const Group& fieldGroup = schedule().getGroup("FIELD", reportStepIdx);
            WellGroupHelpers::updateGuideRates(fieldGroup, schedule(), summaryState, this->phase_usage_, reportStepIdx, simulationTime,
                                        this->wellState(), this->groupState(), comm, &this->guideRate_, pot, local_deferredLogger);
        }

        // Maybe do a recursive call to iterate network and well controls.
        if (network_changed) {
            if (shouldBalanceNetwork(reportStepIdx, iterationIdx)) {
                const auto& balance = schedule()[reportStepIdx].network_balance();
                // Iterate if not converged, and number of iterations is not yet max (NETBALAN item 3).
                if (recursion_level < balance.pressure_max_iter() && network_imbalance > balance.pressure_tolerance()) {
                    well_group_control_changed = assembleImpl(iterationIdx, dt, recursion_level + 1, local_deferredLogger);
                }
            }
        }
        return well_group_control_changed;
    }





    template<typename TypeTag>
    bool
    BlackoilWellModel<TypeTag>::
    maybeDoGasLiftOptimize(DeferredLogger& deferred_logger)
    {
        bool do_glift_optimization = false;
        int num_wells_changed = 0;
        const double simulation_time = ebosSimulator_.time();
        const double min_wait = ebosSimulator_.vanguard().schedule().glo(ebosSimulator_.episodeIndex()).min_wait();
        // We only optimize if a min_wait time has past.
        // If all_newton is true we still want to optimize several times pr timestep
        // i.e. we also optimize if check simulation_time == last_glift_opt_time_
        // that is when the last_glift_opt_time is already updated with the current time step
        if ( simulation_time == last_glift_opt_time_  || simulation_time >= (last_glift_opt_time_ + min_wait)) {
            do_glift_optimization = true;
            last_glift_opt_time_ = simulation_time;
        }

        if (do_glift_optimization) {
            GLiftOptWells glift_wells;
            GLiftProdWells prod_wells;
            GLiftWellStateMap state_map;
            // NOTE: To make GasLiftGroupInfo (see below) independent of the TypeTag
            //  associated with *this (i.e. BlackoilWellModel<TypeTag>) we observe
            //  that GasLiftGroupInfo's only dependence on *this is that it needs to
            //  access the eclipse Wells in the well container (the eclipse Wells
            //  themselves are independent of the TypeTag).
            //  Hence, we extract them from the well container such that we can pass
            //  them to the GasLiftGroupInfo constructor.
            GLiftEclWells ecl_well_map;
            initGliftEclWellMap(ecl_well_map);
            GasLiftGroupInfo group_info {
                ecl_well_map,
                ebosSimulator_.vanguard().schedule(),
                ebosSimulator_.vanguard().summaryState(),
                ebosSimulator_.episodeIndex(),
                ebosSimulator_.model().newtonMethod().numIterations(),
                phase_usage_,
                deferred_logger,
                this->wellState(),
                this->groupState(),
                ebosSimulator_.vanguard().grid().comm(),
                this->glift_debug
            };
            group_info.initialize();
            gasLiftOptimizationStage1(
                deferred_logger, prod_wells, glift_wells, group_info, state_map);
            gasLiftOptimizationStage2(
                deferred_logger, prod_wells, glift_wells, group_info, state_map,
                ebosSimulator_.episodeIndex());
            if (this->glift_debug) gliftDebugShowALQ(deferred_logger);
            num_wells_changed = glift_wells.size();
        }
        num_wells_changed = this->comm_.sum(num_wells_changed);
        return num_wells_changed > 0;
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    gasLiftOptimizationStage1(DeferredLogger& deferred_logger,
        GLiftProdWells &prod_wells, GLiftOptWells &glift_wells,
        GasLiftGroupInfo &group_info, GLiftWellStateMap &state_map)
    {
        auto comm = ebosSimulator_.vanguard().grid().comm();
        int num_procs = comm.size();
        // NOTE: Gas lift optimization stage 1 seems to be difficult
        //  to do in parallel since the wells are optimized on different
        //  processes and each process needs to know the current ALQ allocated
        //  to each group it is a memeber of in order to check group limits and avoid
        //  allocating more ALQ than necessary.  (Surplus ALQ is removed in
        //  stage 2). In stage1, as each well is adding ALQ, the current group ALQ needs
        //  to be communicated to the other processes.  But there is no common
        //  synchronization point that all process will reach in the
        //  runOptimizeLoop_() in GasLiftSingleWell.cpp.
        //
        //  TODO: Maybe a better solution could be invented by distributing
        //    wells according to certain parent groups. Then updated group rates
        //    might not have to be communicated to the other processors.

        //  Currently, the best option seems to be to run this part sequentially
        //    (not in parallel).
        //
        //  TODO: The simplest approach seems to be if a) one process could take
        //    ownership of all the wells (the union of all the wells in the
        //    well_container_ of each process) then this process could do the
        //    optimization, while the other processes could wait for it to
        //    finish (e.g. comm.barrier()), or alternatively, b) if all
        //    processes could take ownership of all the wells.  Then there
        //    would be no need for synchronization here..
        //
        for (int i = 0; i< num_procs; i++) {
            int num_rates_to_sync = 0;  // communication variable
            GLiftSyncGroups groups_to_sync;
            if (comm.rank() ==  i) {
                // Run stage1: Optimize single wells while also checking group limits
                for (const auto& well : well_container_) {
                    // NOTE: Only the wells in "group_info" needs to be optimized
                    if (group_info.hasWell(well->name())) {
                        gasLiftOptimizationStage1SingleWell(
                            well.get(), deferred_logger, prod_wells, glift_wells,
                            group_info, state_map, groups_to_sync
                        );
                    }
                }
                num_rates_to_sync = groups_to_sync.size();
            }
            num_rates_to_sync = comm.sum(num_rates_to_sync);
            if (num_rates_to_sync > 0) {
                std::vector<int> group_indexes;
                group_indexes.reserve(num_rates_to_sync);
                std::vector<double> group_alq_rates;
                group_alq_rates.reserve(num_rates_to_sync);
                std::vector<double> group_oil_rates;
                group_oil_rates.reserve(num_rates_to_sync);
                std::vector<double> group_gas_rates;
                group_gas_rates.reserve(num_rates_to_sync);
                std::vector<double> group_water_rates;
                group_water_rates.reserve(num_rates_to_sync);
                if (comm.rank() == i) {
                    for (auto idx : groups_to_sync) {
                        auto [oil_rate, gas_rate, water_rate, alq] = group_info.getRates(idx);
                        group_indexes.push_back(idx);
                        group_oil_rates.push_back(oil_rate);
                        group_gas_rates.push_back(gas_rate);
                        group_water_rates.push_back(water_rate);
                        group_alq_rates.push_back(alq);
                    }
                } else {
                    group_indexes.resize(num_rates_to_sync);
                    group_oil_rates.resize(num_rates_to_sync);
                    group_gas_rates.resize(num_rates_to_sync);
                    group_water_rates.resize(num_rates_to_sync);
                    group_alq_rates.resize(num_rates_to_sync);
                }
#if HAVE_MPI
                EclMpiSerializer ser(comm);
                ser.broadcast(i, group_indexes, group_oil_rates,
                              group_gas_rates, group_water_rates, group_alq_rates);
#endif
                if (comm.rank() != i) {
                    for (int j=0; j<num_rates_to_sync; j++) {
                        group_info.updateRate(group_indexes[j],
                            group_oil_rates[j], group_gas_rates[j], group_water_rates[j], group_alq_rates[j]);
                    }
                }
                if (this->glift_debug) {
                    int counter = 0;
                    if (comm.rank() == i) {
                        counter = this->wellState().gliftGetDebugCounter();
                    }
                    counter = comm.sum(counter);
                    if (comm.rank() != i) {
                        this->wellState().gliftSetDebugCounter(counter);
                    }
                }
            }
        }
    }

    // NOTE: this method cannot be const since it passes this->wellState()
    //   (see below) to the GasLiftSingleWell constructor which accepts WellState
    //   as a non-const reference..
    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    gasLiftOptimizationStage1SingleWell(WellInterface<TypeTag> *well,
        DeferredLogger& deferred_logger,
        GLiftProdWells &prod_wells, GLiftOptWells &glift_wells,
        GasLiftGroupInfo &group_info, GLiftWellStateMap &state_map,
        GLiftSyncGroups& sync_groups)
    {
        const auto& summary_state = ebosSimulator_.vanguard().summaryState();
        std::unique_ptr<GasLiftSingleWell> glift
            = std::make_unique<GasLiftSingleWell>(
                *well, ebosSimulator_, summary_state,
                deferred_logger, this->wellState(), this->groupState(),
                group_info, sync_groups, this->comm_, this->glift_debug);
        auto state = glift->runOptimize(
            ebosSimulator_.model().newtonMethod().numIterations());
        if (state) {
            state_map.insert({well->name(), std::move(state)});
            glift_wells.insert({well->name(), std::move(glift)});
            return;
        }
        prod_wells.insert({well->name(), well});
    }


    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    initGliftEclWellMap(GLiftEclWells &ecl_well_map)
    {
        for ( const auto& well: well_container_ ) {
            ecl_well_map.try_emplace(
                well->name(), &(well->wellEcl()), well->indexOfWell());
        }
    }


    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    assembleWellEq(const double dt, DeferredLogger& deferred_logger)
    {
        for (auto& well : well_container_) {
            well->assembleWellEq(ebosSimulator_, dt, this->wellState(), this->groupState(), deferred_logger);
        }
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    apply( BVector& r) const
    {
        for (auto& well : well_container_) {
            well->apply(r);
        }
    }


    // Ax = A x - C D^-1 B x
    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    apply(const BVector& x, BVector& Ax) const
    {
        for (auto& well : well_container_) {
            well->apply(x, Ax);
        }
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    getWellContributions(WellContributions& wellContribs) const
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

    // Ax = Ax - alpha * C D^-1 B x
    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    applyScaleAdd(const Scalar alpha, const BVector& x, BVector& Ax) const
    {
        if (this->well_container_.empty()) {
            return;
        }

        if( scaleAddRes_.size() != Ax.size() ) {
            scaleAddRes_.resize( Ax.size() );
        }

        scaleAddRes_ = 0.0;
        // scaleAddRes_  = - C D^-1 B x
        apply( x, scaleAddRes_ );
        // Ax = Ax + alpha * scaleAddRes_
        Ax.axpy( alpha, scaleAddRes_ );
    }

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
    addWellPressureEquations(PressureMatrix& jacobian, const BVector& weights,const bool use_well_weights) const
    {
        int nw =  this->numLocalWellsEnd();
        int rdofs = local_num_cells_;
        for ( int i = 0; i < nw; i++ ){
            int wdof = rdofs + i; 
            jacobian[wdof][wdof] = 1.0;// better scaling ?
        }

        for ( const auto& well : well_container_ ) {
            well->addWellPressureEquations(jacobian, weights, pressureVarIndex, use_well_weights, this->wellState());
        }
    }

    template<typename TypeTag>
    int
    BlackoilWellModel<TypeTag>::
    numLocalWellsEnd() const
    {
        auto w = schedule().getWellsatEnd();
        w.erase(std::remove_if(w.begin(), w.end(), not_on_process_), w.end());
        return w.size();
    }

    template<typename TypeTag>
    std::vector<std::vector<int>>
    BlackoilWellModel<TypeTag>::
    getMaxWellConnections() const
    {
        std::vector<std::vector<int>> wells;

        auto schedule_wells = schedule().getWellsatEnd();
        schedule_wells.erase(std::remove_if(schedule_wells.begin(), schedule_wells.end(), not_on_process_), schedule_wells.end());
        wells.reserve(schedule_wells.size());

        // initialize the additional cell connections introduced by wells.
        for ( const auto& well : schedule_wells )
        {
            std::vector<int> compressed_well_perforations;
            // All possible completions of the well
            const auto& completionSet = well.getConnections();
            compressed_well_perforations.reserve(completionSet.size());

            for (const auto& connection: well.getConnections())
            {
                const int compressed_idx = compressedIndexForInterior(connection.global_index());
                if ( compressed_idx >= 0 ) // Ignore completions in inactive/remote cells.
                {
                    compressed_well_perforations.push_back(compressed_idx);
                }
            }

            // also include wells with no perforations in case
            std::sort(compressed_well_perforations.begin(),
                      compressed_well_perforations.end());

            wells.push_back(compressed_well_perforations);
        }
        return wells;
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    addWellPressureEquationsStruct(PressureMatrix& jacobian) const
    {
        int nw =  this->numLocalWellsEnd();
        int rdofs = local_num_cells_;
        for(int i=0; i < nw; i++){
            int wdof = rdofs + i; 
            jacobian.entry(wdof,wdof) = 1.0;// better scaling ?
        }
        std::vector<std::vector<int>> wellconnections = getMaxWellConnections();
        for(int i=0; i < nw; i++){
            const auto& perfcells = wellconnections[i];
            for(int perfcell : perfcells){
                int wdof = rdofs + i; 
                jacobian.entry(wdof,perfcell) = 0.0;
                jacobian.entry(perfcell, wdof) = 0.0;
            }
        }
    }

    template<typename TypeTag>
    int
    BlackoilWellModel<TypeTag>::
    numLocalNonshutWells() const
    {
        return well_container_.size();
    }
    
    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    recoverWellSolutionAndUpdateWellState(const BVector& x)
    {
         
        DeferredLogger local_deferredLogger;
        OPM_BEGIN_PARALLEL_TRY_CATCH();
        {
            for (auto& well : well_container_) {
                well->recoverWellSolutionAndUpdateWellState(x, this->wellState(), local_deferredLogger);
            }

        }
        OPM_END_PARALLEL_TRY_CATCH_LOG(local_deferredLogger,
                                       "recoverWellSolutionAndUpdateWellState() failed: ",
                                       terminal_output_, ebosSimulator_.vanguard().grid().comm());

    }




    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    initPrimaryVariablesEvaluation() const
    {
        for (auto& well : well_container_) {
            well->initPrimaryVariablesEvaluation();
        }
    }






    template<typename TypeTag>
    ConvergenceReport
    BlackoilWellModel<TypeTag>::
    getWellConvergence(const std::vector<Scalar>& B_avg, bool checkWellGroupControls) const
    {

        DeferredLogger local_deferredLogger;
        // Get global (from all processes) convergence report.
        ConvergenceReport local_report;
        const int iterationIdx = ebosSimulator_.model().newtonMethod().numIterations();
        for (const auto& well : well_container_) {
            if (well->isOperableAndSolvable() || well->wellIsStopped()) {
                local_report += well->getWellConvergence(this->wellState(), B_avg, local_deferredLogger, iterationIdx > param_.strict_outer_iter_wells_ );
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

        if (terminal_output_) {
            global_deferredLogger.logMessages();
        }
        // Log debug messages for NaN or too large residuals.
        if (terminal_output_) {
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
            well->calculateExplicitQuantities(ebosSimulator_, this->wellState(), deferred_logger);
        }
    }





    template<typename TypeTag>
    bool
    BlackoilWellModel<TypeTag>::
    shouldBalanceNetwork(const int reportStepIdx, const int iterationIdx) const
    {
        const auto& balance = schedule()[reportStepIdx].network_balance();
        if (balance.mode() == Network::Balance::CalcMode::TimeStepStart) {
            return iterationIdx == 0;
        } else if (balance.mode() == Network::Balance::CalcMode::NUPCOL) {
            const int nupcol = schedule()[reportStepIdx].nupcol();
            return iterationIdx < nupcol;
        } else {
            // We do not support any other rebalancing modes,
            // i.e. TimeInterval based rebalancing is not available.
            // This should be warned about elsewhere, so we choose to
            // avoid spamming with a warning here.
            return false;
        }
    }





    template<typename TypeTag>
    std::tuple<bool, bool, double>
    BlackoilWellModel<TypeTag>::
    updateWellControls(DeferredLogger& deferred_logger)
    {
        // Even if there are no wells active locally, we cannot
        // return as the DeferredLogger uses global communication.
        // For no well active globally we simply return.
        if( !wellsActive() ) return { false, false, 0.0 };

        const int episodeIdx = ebosSimulator_.episodeIndex();
        const int iterationIdx = ebosSimulator_.model().newtonMethod().numIterations();
        const auto& comm = ebosSimulator_.vanguard().grid().comm();
        updateAndCommunicateGroupData(episodeIdx, iterationIdx);

        const auto [local_network_changed, local_network_imbalance]
            = shouldBalanceNetwork(episodeIdx, iterationIdx) ?
            updateNetworkPressures(episodeIdx) : std::make_pair(false, 0.0);
        const bool network_changed = comm.sum(local_network_changed);
        const double network_imbalance = comm.max(local_network_imbalance);

        bool changed_well_group = false;
        // Check group individual constraints.
        const int nupcol = schedule()[episodeIdx].nupcol();
        // don't switch group control when iterationIdx > nupcol
        // to avoid oscilations between group controls
        if (iterationIdx <= nupcol) {
            const Group& fieldGroup = schedule().getGroup("FIELD", episodeIdx);
            changed_well_group = updateGroupControls(fieldGroup, deferred_logger, episodeIdx, iterationIdx);
        }
        // Check wells' group constraints and communicate.
        bool changed_well_to_group = false;
        for (const auto& well : well_container_) {
            const auto mode = WellInterface<TypeTag>::IndividualOrGroup::Group;
            const bool changed_well = well->updateWellControl(ebosSimulator_, mode, this->wellState(), this->groupState(), deferred_logger);
            if (changed_well) {
                changed_well_to_group = changed_well || changed_well_to_group;
            }
        }

        changed_well_to_group = comm.sum(changed_well_to_group);
        if (changed_well_to_group) {
            updateAndCommunicate(episodeIdx, iterationIdx, deferred_logger);
            changed_well_group = true;
        }

        // Check individual well constraints and communicate.
        bool changed_well_individual = false;
        for (const auto& well : well_container_) {
            const auto mode = WellInterface<TypeTag>::IndividualOrGroup::Individual;
            const bool changed_well = well->updateWellControl(ebosSimulator_, mode, this->wellState(), this->groupState(), deferred_logger);
            if (changed_well) {
                changed_well_individual = changed_well || changed_well_individual;
            }
        }
        changed_well_individual = comm.sum(changed_well_individual);
        if (changed_well_individual) {
            updateAndCommunicate(episodeIdx, iterationIdx, deferred_logger);
            changed_well_group = true;
        }

        // update wsolvent fraction for REIN wells
        const Group& fieldGroup = schedule().getGroup("FIELD", episodeIdx);
        updateWsolvent(fieldGroup, episodeIdx,  this->nupcolWellState());

        return { changed_well_group, network_changed, network_imbalance };
    }


    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    updateAndCommunicate(const int reportStepIdx,
                         const int iterationIdx,
                         DeferredLogger& deferred_logger)
    {
        updateAndCommunicateGroupData(reportStepIdx, iterationIdx);
        // if a well or group change control it affects all wells that are under the same group
        for (const auto& well : well_container_) {
            well->updateWellStateWithTarget(ebosSimulator_, this->groupState(), this->wellState(), deferred_logger);
        }
        updateAndCommunicateGroupData(reportStepIdx, iterationIdx);
    }

    template<typename TypeTag>
    bool
    BlackoilWellModel<TypeTag>::
    updateGroupControls(const Group& group,
                        DeferredLogger& deferred_logger,
                        const int reportStepIdx,
                        const int iterationIdx)
    {
        bool changed = false;
        bool changed_hc = checkGroupHigherConstraints( group, deferred_logger, reportStepIdx);
        if (changed_hc) {
            changed = true;
            updateAndCommunicate(reportStepIdx, iterationIdx, deferred_logger);
        }
        bool changed_individual =
            BlackoilWellModelConstraints(*this).
                updateGroupIndividualControl(group,
                                             reportStepIdx,
                                             this->switched_inj_groups_,
                                             this->switched_prod_groups_,
                                             this->groupState(),
                                             this->wellState(),
                                             deferred_logger);
        if (changed_individual) {
            changed = true;
            updateAndCommunicate(reportStepIdx, iterationIdx, deferred_logger);
        }
        // call recursively down the group hierarchy
        for (const std::string& groupName : group.groups()) {
            bool changed_this = updateGroupControls( schedule().getGroup(groupName, reportStepIdx), deferred_logger, reportStepIdx,iterationIdx);
            changed = changed || changed_this;
        }
        return changed;
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    updateWellTestState(const double& simulationTime, WellTestState& wellTestState) const
    {
        DeferredLogger local_deferredLogger;
        for (const auto& well : well_container_) {
            const auto& wname = well->name();
            const auto wasClosed = wellTestState.well_is_closed(wname);
            well->checkWellOperability(ebosSimulator_, this->wellState(), local_deferredLogger);
            well->updateWellTestState(this->wellState().well(wname), simulationTime, /*writeMessageToOPMLog=*/ true, wellTestState, local_deferredLogger);

            if (!wasClosed && wellTestState.well_is_closed(wname)) {
                this->closed_this_step_.insert(wname);
            }
        }
               
        const Opm::Parallel::Communication comm = grid().comm();
        DeferredLogger global_deferredLogger = gatherDeferredLogger(local_deferredLogger, comm);
        if (terminal_output_) {
            global_deferredLogger.logMessages();
        }
    }


    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::computePotentials(const std::size_t widx,
                                                  const WellState& well_state_copy,
                                                  std::string& exc_msg,
                                                  ExceptionType::ExcEnum& exc_type,
                                                  DeferredLogger& deferred_logger)
    {
        const int np = numPhases();
        std::vector<double> potentials;
        const auto& well= well_container_[widx];
        try {
            well->computeWellPotentials(ebosSimulator_, well_state_copy, potentials, deferred_logger);
        }
        // catch all possible exception and store type and message.
        OPM_PARALLEL_CATCH_CLAUSE(exc_type, exc_msg);
        // Store it in the well state
        // potentials is resized and set to zero in the beginning of well->ComputeWellPotentials
        // and updated only if sucessfull. i.e. the potentials are zero for exceptions
        auto& ws = this->wellState().well(well->indexOfWell());
        for (int p = 0; p < np; ++p) {
            // make sure the potentials are positive
            ws.well_potentials[p] = std::max(0.0, potentials[p]);
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
            if (this->wells_ecl_[shutWell].getConnections().empty()) {
                // No connections in this well.  Nothing to do.
                continue;
            }

            auto wellPtr = this->template createTypedWellPointer
                <StandardWell<TypeTag>>(shutWell, reportStepIdx);

            wellPtr->init(&this->phase_usage_, this->depth_, this->gravity_,
                          this->local_num_cells_, this->B_avg_, true);

            this->calculateProductivityIndexValues(wellPtr.get(), deferred_logger);
        }
    }





    template <typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    calculateProductivityIndexValues(const WellInterface<TypeTag>* wellPtr,
                                     DeferredLogger& deferred_logger)
    {
        wellPtr->updateProductivityIndex(this->ebosSimulator_,
                                         this->prod_index_calc_[wellPtr->indexOfWell()],
                                         this->wellState(),
                                         deferred_logger);
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    prepareTimeStep(DeferredLogger& deferred_logger)
    {
        for (const auto& well : well_container_) {
            auto& events = this->wellState().well(well->indexOfWell()).events;
            if (events.hasEvent(WellState::event_mask)) {
                well->updateWellStateWithTarget(ebosSimulator_, this->groupState(), this->wellState(), deferred_logger);
                well->updatePrimaryVariables(this->wellState(), deferred_logger);
                well->initPrimaryVariablesEvaluation();
                // There is no new well control change input within a report step,
                // so next time step, the well does not consider to have effective events anymore.
                events.clearEvent(WellState::event_mask);
            }
            // solve the well equation initially to improve the initial solution of the well model
            if (param_.solve_welleq_initially_ && well->isOperableAndSolvable()) {
                try {
                    well->solveWellEquation(ebosSimulator_, this->wellState(), this->groupState(), deferred_logger);
                } catch (const std::exception& e) {
                    const std::string msg = "Compute initial well solution for " + well->name() + " initially failed. Continue with the privious rates";
                    deferred_logger.warning("WELL_INITIAL_SOLVE_FAILED", msg);
                }
            }
        }
        updatePrimaryVariables(deferred_logger);
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    updateAverageFormationFactor()
    {
        std::vector< Scalar > B_avg(numComponents(), Scalar() );
        const auto& grid = ebosSimulator_.vanguard().grid();
        const auto& gridView = grid.leafGridView();
        ElementContext elemCtx(ebosSimulator_);

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
            well->updatePrimaryVariables(this->wellState(), deferred_logger);
        }
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::extractLegacyCellPvtRegionIndex_()
    {
        const auto& grid = ebosSimulator_.vanguard().grid();
        const auto& eclProblem = ebosSimulator_.problem();
        const unsigned numCells = grid.size(/*codim=*/0);

        pvt_region_idx_.resize(numCells);
        for (unsigned cellIdx = 0; cellIdx < numCells; ++cellIdx) {
            pvt_region_idx_[cellIdx] =
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
        if (numPhases() < 3) {
            return numPhases();
        }
        int numComp = FluidSystem::numComponents;
        if constexpr (has_solvent_) {
            numComp ++;
        }

        return numComp;
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::extractLegacyDepth_()
    {
        const auto& eclProblem = ebosSimulator_.problem();
        depth_.resize(local_num_cells_);
        for (unsigned cellIdx = 0; cellIdx < local_num_cells_; ++cellIdx) {
            depth_[cellIdx] = eclProblem.dofCenterDepth(cellIdx);
        }
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    updatePerforationIntensiveQuantities()
    {
        ElementContext elemCtx(ebosSimulator_);
        const auto& gridView = ebosSimulator_.gridView();

        OPM_BEGIN_PARALLEL_TRY_CATCH();
        for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
            elemCtx.updatePrimaryStencil(elem);
            int elemIdx = elemCtx.globalSpaceIndex(0, 0);

            if (!is_cell_perforated_[elemIdx]) {
                continue;
            }
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
        }
        OPM_END_PARALLEL_TRY_CATCH("BlackoilWellModel::updatePerforationIntensiveQuantities() failed: ", ebosSimulator_.vanguard().grid().comm());
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

    template<typename TypeTag>
    bool
    BlackoilWellModel<TypeTag>::
    hasWell(const std::string& well_name) const
    {
        return std::any_of(well_container_.begin(), well_container_.end(),
            [&well_name](const WellInterfacePtr& elem) -> bool
        {
            return elem->name() == well_name;
        });
    }




    template <typename TypeTag>
    int
    BlackoilWellModel<TypeTag>::
    reportStepIndex() const
    {
        return std::max(this->ebosSimulator_.episodeIndex(), 0);
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    calcRates(const int fipnum,
              const int pvtreg,
              std::vector<double>& resv_coeff)
    {
        rateConverter_->calcCoeff(fipnum, pvtreg, resv_coeff);
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    calcInjRates(const int fipnum,
              const int pvtreg,
              std::vector<double>& resv_coeff)
    {
        rateConverter_->calcInjCoeff(fipnum, pvtreg, resv_coeff);
    }


    template <typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    computeWellTemperature()
    {
        if (!has_energy_)
            return;

        int np = numPhases();
        double cellInternalEnergy;
        double cellBinv;
        double cellDensity;
        double perfPhaseRate;
        const int nw = numLocalWells();
        for (auto wellID = 0*nw; wellID < nw; ++wellID) {
            const Well& well = wells_ecl_[wellID];
            if (well.isInjector())
                continue;

            int connpos = 0;
            for (int i = 0; i < wellID; ++i) {
                connpos += well_perf_data_[i].size();
            }
            connpos *= np;
            std::array<double,2> weighted{0.0,0.0};
            auto& [weighted_temperature, total_weight] = weighted;

            auto& well_info = local_parallel_well_info_[wellID].get();
            const int num_perf_this_well = well_info.communication().sum(well_perf_data_[wellID].size());
            auto& ws = this->wellState().well(wellID);
            auto& perf_data = ws.perf_data;
            auto& perf_phase_rate = perf_data.phase_rates;

            for (int perf = 0; perf < num_perf_this_well; ++perf) {
                const int cell_idx = well_perf_data_[wellID][perf].cell_index;
                const auto& intQuants = *(ebosSimulator_.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
                const auto& fs = intQuants.fluidState();

                // we on only have one temperature pr cell any phaseIdx will do
                double cellTemperatures = fs.temperature(/*phaseIdx*/0).value();

                double weight_factor = 0.0;
                for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx)
                {
                    if (!FluidSystem::phaseIsActive(phaseIdx)) {
                        continue;
                    }
                    cellInternalEnergy = fs.enthalpy(phaseIdx).value() - fs.pressure(phaseIdx).value() / fs.density(phaseIdx).value();
                    cellBinv = fs.invB(phaseIdx).value();
                    cellDensity = fs.density(phaseIdx).value();
                    perfPhaseRate = perf_phase_rate[ perf*np + phaseIdx ];
                    weight_factor += cellDensity  * perfPhaseRate/cellBinv * cellInternalEnergy/cellTemperatures;
                }
                total_weight += weight_factor;
                weighted_temperature += weight_factor * cellTemperatures;
            }
            well_info.communication().sum(weighted.data(), 2);
            this->wellState().well(wellID).temperature = weighted_temperature/total_weight;
        }
    }



    template <typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    assignWellTracerRates(data::Wells& wsrpt) const
    {
        const auto & wellTracerRates = ebosSimulator_.problem().tracerModel().getWellTracerRates();

        if (wellTracerRates.empty())
            return; // no tracers

        for (const auto& wTR : wellTracerRates) {
            std::string wellName = wTR.first.first;
            auto xwPos = wsrpt.find(wellName);
            if (xwPos == wsrpt.end()) { // No well results.
                continue;
            }
            std::string tracerName = wTR.first.second;
            double rate = wTR.second;
            xwPos->second.rates.set(data::Rates::opt::tracer, rate, tracerName);
        }
    }





} // namespace Opm
