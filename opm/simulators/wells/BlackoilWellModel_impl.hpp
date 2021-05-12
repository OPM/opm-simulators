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

#include <opm/parser/eclipse/Units/UnitSystem.hpp>

#include <algorithm>
#include <utility>

#include <fmt/format.h>

namespace Opm {
    template<typename TypeTag>
    BlackoilWellModel<TypeTag>::
    BlackoilWellModel(Simulator& ebosSimulator, const PhaseUsage& phase_usage)
        : ebosSimulator_(ebosSimulator)
        , terminal_output_((ebosSimulator.gridView().comm().rank() == 0) &&
                           EWOMS_GET_PARAM(TypeTag, bool, EnableTerminalOutput))
        , phase_usage_(phase_usage)
        , active_wgstate_(phase_usage)
        , last_valid_wgstate_(phase_usage)
        , nupcol_wgstate_(phase_usage)
    {
        // Create the guide rate container.
        this->guideRate_ =
            std::make_unique<GuideRate>(ebosSimulator_.vanguard().schedule());

        local_num_cells_ = ebosSimulator_.gridView().size(0);

        // Number of cells the global grid view
        global_num_cells_ = ebosSimulator_.vanguard().globalNumCells();

        // Set up cartesian mapping.
        {
            const auto& grid = this->ebosSimulator_.vanguard().grid();
            const auto& cartDims = UgGridHelpers::cartDims(grid);
            setupCartesianToCompressed_(UgGridHelpers::globalCell(grid),
                                        cartDims[0] * cartDims[1] * cartDims[2]);

            auto& parallel_wells = ebosSimulator.vanguard().parallelWells();
            this->parallel_well_info_.assign(parallel_wells.begin(),
                                             parallel_wells.end());
        }

        const auto numProcs = ebosSimulator.gridView().comm().size();
        this->not_on_process_ = [this, numProcs](const Well& well) {
            if (numProcs == decltype(numProcs){1})
                return false;

            // Recall: false indicates NOT active!
            const auto value = std::make_pair(well.name(), true);
            auto candidate = std::lower_bound(this->parallel_well_info_.begin(),
                                              this->parallel_well_info_.end(),
                                              value);

            return (candidate == this->parallel_well_info_.end())
                || (*candidate != value);
        };

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
                int compressed_idx = cartesian_to_compressed_
                    .at(connection.global_index());

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
        if (!localWellsActive())
            return;

        if (!param_.matrix_add_well_contributions_) {
            // if the well contributions are not supposed to be included explicitly in
            // the matrix, we only apply the vector part of the Schur complement here.
            for (const auto& well: well_container_) {
                // r = r - duneC_^T * invDuneD_ * resWell_
                well->apply(res);
            }
            return;
        }

        for (const auto& well: well_container_) {
            well->addWellContributions(jacobian);

            // applying the well residual to reservoir residuals
            // r = r - duneC_^T * invDuneD_ * resWell_
            well->apply(res);
        }
    }


    /// Return true if any well has a THP constraint.
    template<typename TypeTag>
    bool
    BlackoilWellModel<TypeTag>::
    hasTHPConstraints() const
    {
        int local_result = false;
        const auto& summaryState = ebosSimulator_.vanguard().summaryState();
        for (const auto& well : well_container_) {
            if (well->wellHasTHPConstraints(summaryState)) {
                local_result=true;
            }
        }
        return grid().comm().max(local_result);
    }




    /// Return true if the well was found and shut.
    template<typename TypeTag>
    bool
    BlackoilWellModel<TypeTag>::
    forceShutWellByNameIfPredictionMode(const std::string& wellname,
                                        const double simulation_time)
    {
        // Only add the well to the closed list on the
        // process that owns it.
        int well_was_shut = 0;
        for (const auto& well : well_container_) {
            if (well->name() == wellname && !well->wellIsStopped()) {
                if (well->underPredictionMode()) {
                    wellTestState_.closeWell(wellname, WellTestConfig::Reason::PHYSICAL, simulation_time);
                    well_was_shut = 1;
                }
                break;
            }
        }

        // Communicate across processes if a well was shut.
        well_was_shut = ebosSimulator_.vanguard().grid().comm().max(well_was_shut);

        // Only log a message on the output rank.
        if (terminal_output_ && well_was_shut) {
            const std::string msg = "Well " + wellname
                + " will be shut because it cannot get converged.";
            OpmLog::info(msg);
        }

        return (well_was_shut == 1);
    }


    template<typename TypeTag>
    std::vector< Well >
    BlackoilWellModel<TypeTag>::
    getLocalWells(const int timeStepIdx) const
    {
        auto w = schedule().getWells(timeStepIdx);
        w.erase(std::remove_if(w.begin(), w.end(), not_on_process_), w.end());
        return w;
    }

    template<typename TypeTag>
    std::vector< ParallelWellInfo* >
    BlackoilWellModel<TypeTag>::createLocalParallelWellInfo(const std::vector<Well>& wells)
    {
        std::vector< ParallelWellInfo* > local_parallel_well_info;
        local_parallel_well_info.reserve(wells.size());
        for (const auto& well : wells)
        {
            auto wellPair = std::make_pair(well.name(), true);
            auto pwell = std::lower_bound(parallel_well_info_.begin(),
                                          parallel_well_info_.end(),
                                          wellPair);
            assert(pwell != parallel_well_info_.end() &&
                   *pwell == wellPair);
            local_parallel_well_info.push_back(&(*pwell));
        }
        return local_parallel_well_info;
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
        local_parallel_well_info_ = createLocalParallelWellInfo(wells_ecl_);

        // The well state initialize bhp with the cell pressure in the top cell.
        // We must therefore provide it with updated cell pressures
        this->initializeWellPerfData();
        this->initializeWellState(timeStepIdx, summaryState);

        // Wells are active if they are active wells on at least
        // one process.
        wells_active_ = localWellsActive() ? 1 : 0;
        wells_active_ = grid.comm().max(wells_active_);

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

        {
            const auto& sched_state = this->schedule()[timeStepIdx];
            // update VFP properties
            vfp_properties_.reset(new VFPProperties( sched_state.vfpinj(), sched_state.vfpprod()) );
            this->initializeWellProdIndCalculators();
            if (sched_state.events().hasEvent(ScheduleEvents::Events::WELL_PRODUCTIVITY_INDEX)) {
                this->runWellPIScaling(timeStepIdx, local_deferredLogger);
            }
        }

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

        this->resetWGState();
        updateAndCommunicateGroupData();
        this->wellState().gliftTimeStepInit();
        const int reportStepIdx = ebosSimulator_.episodeIndex();
        const double simulationTime = ebosSimulator_.time();
        std::string exc_msg;
        auto exc_type = ExceptionType::NONE;
        try {
            // test wells
            wellTesting(reportStepIdx, simulationTime, local_deferredLogger);

            // create the well container
            well_container_ = createWellContainer(reportStepIdx);

            // do the initialization for all the wells
            // TODO: to see whether we can postpone of the intialization of the well containers to
            // optimize the usage of the following several member variables
            for (auto& well : well_container_) {
                well->init(&phase_usage_, depth_, gravity_, local_num_cells_, B_avg_);
            }

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
        } catch (const std::runtime_error& e) {
            exc_type = ExceptionType::RUNTIME_ERROR;
            exc_msg = e.what();
        } catch (const std::invalid_argument& e) {
            exc_type = ExceptionType::INVALID_ARGUMENT;
            exc_msg = e.what();
        } catch (const std::logic_error& e) {
            exc_type = ExceptionType::LOGIC_ERROR;
            exc_msg = e.what();
        } catch (const std::exception& e) {
            exc_type = ExceptionType::DEFAULT;
            exc_msg = e.what();
        }

        logAndCheckForExceptionsAndThrow(local_deferredLogger, exc_type, "beginTimeStep() failed: " + exc_msg, terminal_output_);

        for (auto& well : well_container_) {
            well->setVFPProperties(vfp_properties_.get());
            well->setGuideRate(guideRate_.get());
        }

        // Close completions due to economical reasons
        for (auto& well : well_container_) {
            well->closeCompletions(wellTestState_);
        }

        // calculate the well potentials
        try {
            std::vector<double> well_potentials;
            computeWellPotentials(well_potentials, reportStepIdx, local_deferredLogger);
        } catch ( std::runtime_error& e ) {
            const std::string msg = "A zero well potential is returned for output purposes. ";
            local_deferredLogger.warning("WELL_POTENTIAL_CALCULATION_FAILED", msg);
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

        //update guide rates
        const auto& comm = ebosSimulator_.vanguard().grid().comm();
        const auto& summaryState = ebosSimulator_.vanguard().summaryState();
        std::vector<double> pot(numPhases(), 0.0);
        const Group& fieldGroup = schedule().getGroup("FIELD", reportStepIdx);
        WellGroupHelpers::updateGuideRateForProductionGroups(fieldGroup, schedule(), phase_usage_, reportStepIdx, simulationTime, this->wellState(), this->groupState(), comm, guideRate_.get(), pot);
        WellGroupHelpers::updateGuideRatesForInjectionGroups(fieldGroup, schedule(), summaryState, phase_usage_, reportStepIdx, this->wellState(), this->groupState(), guideRate_.get(), local_deferredLogger);
        WellGroupHelpers::updateGuideRatesForWells(schedule(), phase_usage_, reportStepIdx, simulationTime, this->wellState(), comm, guideRate_.get());

        try {
            // Compute initial well solution for new wells and injectors that change injection type i.e. WAG.
            for (auto& well : well_container_) {
                const uint64_t effective_events_mask = ScheduleEvents::WELL_STATUS_CHANGE
                        + ScheduleEvents::INJECTION_TYPE_CHANGED
                        + ScheduleEvents::WELL_SWITCHED_INJECTOR_PRODUCER
                        + ScheduleEvents::NEW_WELL;

                const auto& events = schedule()[reportStepIdx].wellgroup_events();
                const bool event = report_step_starts_ && events.hasEvent(well->name(), effective_events_mask);
                if (event) {
                    try {
                        well->updateWellStateWithTarget(ebosSimulator_, this->wellState(), local_deferredLogger);
                        well->calculateExplicitQuantities(ebosSimulator_, this->wellState(), local_deferredLogger);
                        well->solveWellEquation(ebosSimulator_, this->wellState(), this->groupState(), local_deferredLogger);
                    } catch (const std::exception& e) {
                        const std::string msg = "Compute initial well solution for new well " + well->name() + " failed. Continue with zero initial rates";
                        local_deferredLogger.warning("WELL_INITIAL_SOLVE_FAILED", msg);
                    }
                }
            }
        } catch (const std::runtime_error& e) {
            exc_type = ExceptionType::RUNTIME_ERROR;
            exc_msg = e.what();
        } catch (const std::invalid_argument& e) {
            exc_type = ExceptionType::INVALID_ARGUMENT;
            exc_msg = e.what();
        } catch (const std::logic_error& e) {
            exc_type = ExceptionType::LOGIC_ERROR;
            exc_msg = e.what();
        } catch (const std::exception& e) {
            exc_type = ExceptionType::DEFAULT;
            exc_msg = e.what();
        }

        if (exc_type != ExceptionType::NONE) {
            const std::string msg = "Compute initial well solution for new wells failed. Continue with zero initial rates";
            local_deferredLogger.warning("WELL_INITIAL_SOLVE_FAILED", msg);
        }

        logAndCheckForExceptionsAndThrow(local_deferredLogger,
            exc_type, "beginTimeStep() failed: " + exc_msg, terminal_output_);

    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::gliftDebug(
        const std::string &msg, DeferredLogger &deferred_logger) const
    {
        if (this->glift_debug) {
            const std::string message = fmt::format(
                "  GLIFT (DEBUG) : BlackoilWellModel : {}", msg);
            deferred_logger.info(message);
        }
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::wellTesting(const int timeStepIdx,
                                            const double simulationTime,
                                            DeferredLogger& deferred_logger)
    {
        const auto& wtest_config = schedule()[timeStepIdx].wtest_config();
        if (wtest_config.size() != 0) { // there is a WTEST request
            const auto wellsForTesting = wellTestState_
                .updateWells(wtest_config, wells_ecl_, simulationTime);

            for (const auto& testWell : wellsForTesting) {
                const std::string& well_name = testWell.first;

                // this is the well we will test
                WellInterfacePtr well = createWellForWellTest(well_name, timeStepIdx, deferred_logger);

                // some preparation before the well can be used
                well->init(&phase_usage_, depth_, gravity_, local_num_cells_, B_avg_);
                const Well& wellEcl = schedule().getWell(well_name, timeStepIdx);
                double well_efficiency_factor = wellEcl.getEfficiencyFactor();
                WellGroupHelpers::accumulateGroupEfficiencyFactor(schedule().getGroup(wellEcl.groupName(), timeStepIdx),
                                                                  schedule(), timeStepIdx, well_efficiency_factor);

                well->setWellEfficiencyFactor(well_efficiency_factor);
                well->setVFPProperties(vfp_properties_.get());
                well->setGuideRate(guideRate_.get());

                const WellTestConfig::Reason testing_reason = testWell.second;

                well->wellTesting(ebosSimulator_, simulationTime, timeStepIdx,
                                  testing_reason, this->wellState(), this->groupState(), wellTestState_, deferred_logger);
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
        for (auto&& pinfo : local_parallel_well_info_)
        {
            pinfo->clear();
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
        updateWellTestState(simulationTime, wellTestState_);

        // update the rate converter with current averages pressures etc in
        rateConverter_->template defineState<ElementContext>(ebosSimulator_);

        // calculate the well potentials
        try {
            std::vector<double> well_potentials;

            computeWellPotentials(well_potentials, reportStepIdx, local_deferredLogger);
        } catch ( std::runtime_error& e ) {
            const std::string msg = "A zero well potential is returned for output purposes. ";
            local_deferredLogger.warning("WELL_POTENTIAL_CALCULATION_FAILED", msg);
        }

        // check group sales limits at the end of the timestep
        const Group& fieldGroup = schedule().getGroup("FIELD", reportStepIdx);
        checkGconsaleLimits(fieldGroup, this->wellState(), local_deferredLogger);

        this->calculateProductivityIndexValues(local_deferredLogger);

        this->commitWGState();

        DeferredLogger global_deferredLogger = gatherDeferredLogger(local_deferredLogger);
        if (terminal_output_) {
            global_deferredLogger.logMessages();
        }

        //reporting output temperatures
        this->computeWellTemperature();
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
    typename BlackoilWellModel<TypeTag>::WellInterfacePtr
    BlackoilWellModel<TypeTag>::
    well(const std::string& wellName) const
    {
        for (const auto& well : well_container_) {
            if (well->name() == wellName) {
                return well;
            }
        }
        OPM_THROW(std::invalid_argument, "The well with name " + wellName + " is not in the well Container");
        return nullptr;
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    initFromRestartFile(const RestartValue& restartValues)
    {
        // The restart step value is used to identify wells present at the given
        // time step. Wells that are added at the same time step as RESTART is initiated
        // will not be present in a restart file. Use the previous time step to retrieve
        // wells that have information written to the restart file.
        const int report_step = std::max(eclState().getInitConfig().getRestartStep() - 1, 0);
        const auto& summaryState = ebosSimulator_.vanguard().summaryState();
        // wells_ecl_ should only contain wells on this processor.
        wells_ecl_ = getLocalWells(report_step);
        local_parallel_well_info_ = createLocalParallelWellInfo(wells_ecl_);

        this->initializeWellProdIndCalculators();
        initializeWellPerfData();

        const int nw = wells_ecl_.size();
        if (nw > 0) {
            const auto phaseUsage = phaseUsageFromDeck(eclState());
            const size_t numCells = UgGridHelpers::numCells(grid());
            const bool handle_ms_well = (param_.use_multisegment_well_ && anyMSWellOpenLocal());
            this->wellState().resize(wells_ecl_, local_parallel_well_info_, schedule(), handle_ms_well, numCells, well_perf_data_, summaryState); // Resize for restart step
            loadRestartData(restartValues.wells, restartValues.grp_nwrk, phaseUsage, handle_ms_well, this->wellState());
        }

        this->commitWGState();
        initial_step_ = false;
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    initializeWellProdIndCalculators()
    {
        this->prod_index_calc_.clear();
        this->prod_index_calc_.reserve(this->wells_ecl_.size());
        for (const auto& well : this->wells_ecl_) {
            this->prod_index_calc_.emplace_back(well);
        }
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    initializeWellPerfData()
    {
        well_perf_data_.resize(wells_ecl_.size());
        int well_index = 0;
        for (const auto& well : wells_ecl_) {
            int completion_index = 0;
            // INVALID_ECL_INDEX marks no above perf available
            int completion_index_above = ParallelWellInfo::INVALID_ECL_INDEX;
            well_perf_data_[well_index].clear();
            well_perf_data_[well_index].reserve(well.getConnections().size());
            CheckDistributedWellConnections checker(well, *local_parallel_well_info_[well_index]);
            bool hasFirstPerforation = false;
            bool firstOpenCompletion = true;
            auto& parallelWellInfo = *local_parallel_well_info_[well_index];
            parallelWellInfo.beginReset();

            for (const auto& completion : well.getConnections()) {
                const int active_index =
                    cartesian_to_compressed_[completion.global_index()];
                if (completion.state() == Connection::State::OPEN) {
                    if (active_index >= 0) {
                        if (firstOpenCompletion)
                        {
                            hasFirstPerforation = true;
                        }
                        checker.connectionFound(completion_index);
                        PerforationData pd;
                        pd.cell_index = active_index;
                        pd.connection_transmissibility_factor = completion.CF();
                        pd.satnum_id = completion.satTableId();
                        pd.ecl_index = completion_index;
                        well_perf_data_[well_index].push_back(pd);
                        parallelWellInfo.pushBackEclIndex(completion_index_above,
                                                          completion_index);
                    }
                    firstOpenCompletion = false;
                    // Next time this index is the one above as each open completion is
                    // is stored somehwere.
                    completion_index_above = completion_index;
                } else {
                    checker.connectionFound(completion_index);
                    if (completion.state() != Connection::State::SHUT) {
                        OPM_THROW(std::runtime_error,
                                  "Completion state: " << Connection::State2String(completion.state()) << " not handled");
                    }
                }
                // Note: we rely on the connections being filtered! I.e. there are only connections
                // to active cells in the global grid.
                ++completion_index;
            }
            parallelWellInfo.endReset();
            checker.checkAllConnectionsFound();
            parallelWellInfo.communicateFirstPerforation(hasFirstPerforation);
            ++well_index;
        }
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
        const auto& elemEndIt = gridView.template end</*codim=*/0>();
        for (auto elemIt = gridView.template begin</*codim=*/0>();
             elemIt != elemEndIt;
             ++elemIt)
        {
            if (elemIt->partitionType() != Dune::InteriorEntity) {
                continue;
            }

            elemCtx.updatePrimaryStencil(*elemIt);
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

        this->wellState().init(cellPressures, schedule(), wells_ecl_, local_parallel_well_info_, timeStepIdx,
                               &this->prevWellState(), well_perf_data_,
                               summaryState);
    }





    template<typename TypeTag>
    std::vector<typename BlackoilWellModel<TypeTag>::WellInterfacePtr >
    BlackoilWellModel<TypeTag>::
    createWellContainer(const int time_step)
    {
        std::vector<WellInterfacePtr> well_container;

        DeferredLogger local_deferredLogger;

        const int nw = numLocalWells();

        if (nw > 0) {
            well_container.reserve(nw);

            for (int w = 0; w < nw; ++w) {
                const Well& well_ecl = wells_ecl_[w];
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
                if (this->wellTestState_.hasWellClosed(well_name)) {
                    // TODO: more checking here, to make sure this standard more specific and complete
                    // maybe there is some WCON keywords will not open the well
                    auto& events = this->wellState().events(w);
                    if (events.hasEvent(WellStateFullyImplicitBlackoil::event_mask)) {
                        if (wellTestState_.lastTestTime(well_name) == ebosSimulator_.time()) {
                            // The well was shut this timestep, we are most likely retrying
                            // a timestep without the well in question, after it caused
                            // repeated timestep cuts. It should therefore not be opened,
                            // even if it was new or received new targets this report step.
                            events.clearEvent(WellStateFullyImplicitBlackoil::event_mask);
                        } else {
                            wellTestState_.openWell(well_name);
                        }
                    }
                }

                // TODO: should we do this for all kinds of closing reasons?
                // something like wellTestState_.hasWell(well_name)?
                bool wellIsStopped = false;
                if (wellTestState_.hasWellClosed(well_name, WellTestConfig::Reason::ECONOMIC) ||
                    wellTestState_.hasWellClosed(well_name, WellTestConfig::Reason::PHYSICAL))
                {
                    if (well_ecl.getAutomaticShutIn()) {
                        // shut wells are not added to the well container
                        this->wellState().shutWell(w);
                        continue;
                    } else {
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

                well_container.emplace_back(this->createWellPointer(w, time_step));

                if (wellIsStopped)
                    well_container.back()->stopWell();
            }
        }

        // Collect log messages and print.
        DeferredLogger global_deferredLogger = gatherDeferredLogger(local_deferredLogger);
        if (terminal_output_) {
            global_deferredLogger.logMessages();
        }

        return well_container;
    }





    template <typename TypeTag>
    void BlackoilWellModel<TypeTag>::
    inferLocalShutWells()
    {
        this->local_shut_wells_.clear();

        const auto nw = this->numLocalWells();

        auto used = std::vector<bool>(nw, false);
        for (const auto& wellPtr : this->well_container_) {
            used[wellPtr->indexOfWell()] = true;
        }

        for (auto wellID = 0; wellID < nw; ++wellID) {
            if (! used[wellID]) {
                this->local_shut_wells_.push_back(wellID);
            }
        }
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

        const auto& parallel_well_info = *local_parallel_well_info_[wellID];
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
                                          this->wellState().firstPerfIndex()[wellID],
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

        auto exc_type = ExceptionType::NONE;
        std::string exc_msg;
        try {
            if (iterationIdx == 0) {
                calculateExplicitQuantities(local_deferredLogger);
                prepareTimeStep(local_deferredLogger);
            }
            updateWellControls(local_deferredLogger, /* check group controls */ true);

            // Set the well primary variables based on the value of well solutions
            initPrimaryVariablesEvaluation();

            maybeDoGasLiftOptimize(local_deferredLogger);
            assembleWellEq(dt, local_deferredLogger);
        } catch (const std::runtime_error& e) {
            exc_type = ExceptionType::RUNTIME_ERROR;
            exc_msg = e.what();
        } catch (const std::invalid_argument& e) {
            exc_type = ExceptionType::INVALID_ARGUMENT;
            exc_msg = e.what();
        } catch (const std::logic_error& e) {
            exc_type = ExceptionType::LOGIC_ERROR;
            exc_msg = e.what();
        } catch (const std::exception& e) {
            exc_type = ExceptionType::DEFAULT;
            exc_msg = e.what();
        }
        logAndCheckForExceptionsAndThrow(local_deferredLogger, exc_type, "assemble() failed: " + exc_msg, terminal_output_);
        last_report_.converged = true;
        last_report_.assemble_time_well += perfTimer.stop();
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    maybeDoGasLiftOptimize(DeferredLogger& deferred_logger)
    {
        this->wellState().enableGliftOptimization();
        GLiftOptWells glift_wells;
        GLiftProdWells prod_wells;
        GLiftWellStateMap state_map;
        // Stage1: Optimize single wells not checking any group limits
        for (auto& well : well_container_) {
            well->gasLiftOptimizationStage1(
                this->wellState(), ebosSimulator_, deferred_logger,
                prod_wells, glift_wells, state_map);
        }
        gasLiftOptimizationStage2(deferred_logger, prod_wells, glift_wells, state_map);
        if (this->glift_debug) gliftDebugShowALQ(deferred_logger);
        this->wellState().disableGliftOptimization();
    }

    // If a group has any production rate constraints, and/or a limit
    // on its total rate of lift gas supply,  allocate lift gas
    // preferentially to the wells that gain the most benefit from
    // it. Lift gas increments are allocated in turn to the well that
    // currently has the largest weighted incremental gradient. The
    // procedure takes account of any limits on the group production
    // rate or lift gas supply applied to any level of group.
    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    gasLiftOptimizationStage2(DeferredLogger& deferred_logger,
        GLiftProdWells &prod_wells, GLiftOptWells &glift_wells,
        GLiftWellStateMap &glift_well_state_map)
    {

        GasLiftStage2 glift {*this, ebosSimulator_, deferred_logger, this->wellState(),
                             prod_wells, glift_wells, glift_well_state_map};
        glift.runOptimize();
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    gliftDebugShowALQ(DeferredLogger& deferred_logger)
    {
        for (auto& well : this->well_container_) {
            if (well->isProducer()) {
                auto alq = this->wellState().getALQ(well->name());
                const std::string msg = fmt::format("ALQ_REPORT : {} : {}",
                    well->name(), alq);
                gliftDebug(msg, deferred_logger);
            }
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
        if ( ! localWellsActive() ) {
            return;
        }

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
        // TODO: do we still need localWellsActive()?
        if ( ! localWellsActive() ) {
            return;
        }

        for (auto& well : well_container_) {
            well->apply(x, Ax);
        }
    }

#if HAVE_CUDA || HAVE_OPENCL
    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    getWellContributions(WellContributions& wellContribs) const
    {
        // prepare for StandardWells
        wellContribs.setBlockSize(StandardWell<TypeTag>::numEq, StandardWell<TypeTag>::numStaticWellEq);

        for(unsigned int i = 0; i < well_container_.size(); i++){
            auto& well = well_container_[i];
            std::shared_ptr<StandardWell<TypeTag> > derived = std::dynamic_pointer_cast<StandardWell<TypeTag> >(well);
            if (derived) {
                unsigned int numBlocks;
                derived->getNumBlocks(numBlocks);
                wellContribs.addNumBlocks(numBlocks);
            }
        }

        // allocate memory for data from StandardWells
        wellContribs.alloc();

        for(unsigned int i = 0; i < well_container_.size(); i++){
            auto& well = well_container_[i];
            // maybe WellInterface could implement addWellContribution()
            auto derived_std = std::dynamic_pointer_cast<StandardWell<TypeTag> >(well);
            if (derived_std) {
                derived_std->addWellContribution(wellContribs);
            } else {
                auto derived_ms = std::dynamic_pointer_cast<MultisegmentWell<TypeTag> >(well);
                if (derived_ms) {
                    derived_ms->addWellContribution(wellContribs);
                } else {
                    OpmLog::warning("Warning unknown type of well");
                }
            }
        }
    }
#endif

    // Ax = Ax - alpha * C D^-1 B x
    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    applyScaleAdd(const Scalar alpha, const BVector& x, BVector& Ax) const
    {
        if ( ! localWellsActive() ) {
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
    recoverWellSolutionAndUpdateWellState(const BVector& x)
    {
        DeferredLogger local_deferredLogger;
        auto exc_type = ExceptionType::NONE;
        std::string exc_msg;
        try {
            if (localWellsActive()) {
                for (auto& well : well_container_) {
                    well->recoverWellSolutionAndUpdateWellState(x, this->wellState(), local_deferredLogger);
                }
            }
        } catch (const std::runtime_error& e) {
            exc_type = ExceptionType::RUNTIME_ERROR;
            exc_msg = e.what();
        } catch (const std::invalid_argument& e) {
            exc_type = ExceptionType::INVALID_ARGUMENT;
            exc_msg = e.what();
        } catch (const std::logic_error& e) {
            exc_type = ExceptionType::LOGIC_ERROR;
            exc_msg = e.what();
        } catch (const std::exception& e) {
            exc_type = ExceptionType::DEFAULT;
            exc_msg = e.what();
        }
        logAndCheckForExceptionsAndThrow(local_deferredLogger, exc_type, "recoverWellSolutionAndUpdateWellState() failed: " + exc_msg, terminal_output_);
    }




    template<typename TypeTag>
    bool
    BlackoilWellModel<TypeTag>::
    wellsActive() const
    {
        return wells_active_;
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    setWellsActive(const bool wells_active)
    {
        wells_active_ = wells_active;
    }





    template<typename TypeTag>
    bool
    BlackoilWellModel<TypeTag>::
    localWellsActive() const
    {
        return numLocalWells() > 0;
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
    getWellConvergence(const std::vector<Scalar>& B_avg, bool checkGroupConvergence) const
    {

        DeferredLogger local_deferredLogger;
        // Get global (from all processes) convergence report.
        ConvergenceReport local_report;
        for (const auto& well : well_container_) {
            if (well->isOperable() ) {
                local_report += well->getWellConvergence(this->wellState(), B_avg, local_deferredLogger);
            }
        }
        DeferredLogger global_deferredLogger = gatherDeferredLogger(local_deferredLogger);
        if (terminal_output_) {
            global_deferredLogger.logMessages();
        }

        ConvergenceReport report = gatherConvergenceReport(local_report);

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

        if (checkGroupConvergence) {
            const int reportStepIdx = ebosSimulator_.episodeIndex();
            const Group& fieldGroup = schedule().getGroup("FIELD", reportStepIdx);
            bool violated = checkGroupConstraints(fieldGroup, global_deferredLogger);
            report.setGroupConverged(!violated);
        }
        return report;
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    calculateExplicitQuantities(DeferredLogger& deferred_logger) const
    {
        // TODO: checking isOperable() ?
        for (auto& well : well_container_) {
            well->calculateExplicitQuantities(ebosSimulator_, this->wellState(), deferred_logger);
        }
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    updateWellControls(DeferredLogger& deferred_logger, const bool checkGroupControls)
    {
        // Even if there are no wells active locally, we cannot
        // return as the DeferredLogger uses global communication.
        // For no well active globally we simply return.
        if( !wellsActive() ) return ;

        updateAndCommunicateGroupData();

        updateNetworkPressures();

        std::set<std::string> switched_wells;
        std::set<std::string> switched_groups;

        if (checkGroupControls) {
            // Check group individual constraints.
            updateGroupIndividualControls(deferred_logger, switched_groups);

            // Check group's constraints from higher levels.
            updateGroupHigherControls(deferred_logger, switched_groups);

            updateAndCommunicateGroupData();

            // Check wells' group constraints and communicate.
            for (const auto& well : well_container_) {
                const auto mode = WellInterface<TypeTag>::IndividualOrGroup::Group;
                const bool changed = well->updateWellControl(ebosSimulator_, mode, this->wellState(), this->groupState(), deferred_logger);
                if (changed) {
                    switched_wells.insert(well->name());
                }
            }
            updateAndCommunicateGroupData();
        }

        // Check individual well constraints and communicate.
        for (const auto& well : well_container_) {
            if (switched_wells.count(well->name())) {
                continue;
            }
            const auto mode = WellInterface<TypeTag>::IndividualOrGroup::Individual;
            well->updateWellControl(ebosSimulator_, mode, this->wellState(), this->groupState(), deferred_logger);
        }
        updateAndCommunicateGroupData();

    }




    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    updateNetworkPressures()
    {
        // Get the network and return if inactive.
        const int reportStepIdx = ebosSimulator_.episodeIndex();
        const auto& network = schedule()[reportStepIdx].network();
        if (!network.active()) {
            return;
        }
        node_pressures_ = WellGroupHelpers::computeNetworkPressures(network, this->wellState(), this->groupState(), *(vfp_properties_->getProd()), schedule(), reportStepIdx);

        // Set the thp limits of wells
        for (auto& well : well_container_) {
            // Producers only, since we so far only support the
            // "extended" network model (properties defined by
            // BRANPROP and NODEPROP) which only applies to producers.
            if (well->isProducer()) {
                const auto it = node_pressures_.find(well->wellEcl().groupName());
                if (it != node_pressures_.end()) {
                    // The well belongs to a group with has a network pressure constraint,
                    // set the dynamic THP constraint of the well accordingly.
                    well->setDynamicThpLimit(it->second);
                }
            }
        }
    }




    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    updateAndCommunicateGroupData()
    {
        const int reportStepIdx = ebosSimulator_.episodeIndex();
        const Group& fieldGroup = schedule().getGroup("FIELD", reportStepIdx);
        const int nupcol = schedule()[reportStepIdx].nupcol();
        const int iterationIdx = ebosSimulator_.model().newtonMethod().numIterations();

        // This builds some necessary lookup structures, so it must be called
        // before we copy to well_state_nupcol_.
        const auto& comm = ebosSimulator_.vanguard().grid().comm();
        this->wellState().updateGlobalIsGrup(comm);

        if (iterationIdx < nupcol) {
            this->updateNupcolWGState();
        }

        auto& well_state = this->wellState();
        const auto& well_state_nupcol = this->nupcolWellState();
        const auto& summaryState = ebosSimulator_.vanguard().summaryState();
        // the group target reduction rates needs to be update since wells may have swicthed to/from GRUP control
        // Currently the group target reduction does not honor NUPCOL. TODO: is that true?
        std::vector<double> groupTargetReduction(numPhases(), 0.0);
        WellGroupHelpers::updateGroupTargetReduction(fieldGroup, schedule(), reportStepIdx, /*isInjector*/ false, phase_usage_, *guideRate_, well_state_nupcol, well_state, this->groupState(), groupTargetReduction);
        std::vector<double> groupTargetReductionInj(numPhases(), 0.0);
        WellGroupHelpers::updateGroupTargetReduction(fieldGroup, schedule(), reportStepIdx, /*isInjector*/ true, phase_usage_, *guideRate_, well_state_nupcol, well_state, this->groupState(), groupTargetReductionInj);

        WellGroupHelpers::updateREINForGroups(fieldGroup, schedule(), reportStepIdx, phase_usage_, summaryState, well_state_nupcol, well_state, this->groupState());
        WellGroupHelpers::updateVREPForGroups(fieldGroup, schedule(), reportStepIdx, well_state_nupcol, well_state, this->groupState());

        WellGroupHelpers::updateReservoirRatesInjectionGroups(fieldGroup, schedule(), reportStepIdx, well_state_nupcol, well_state, this->groupState());
        WellGroupHelpers::updateGroupProductionRates(fieldGroup, schedule(), reportStepIdx, well_state_nupcol, well_state, this->groupState());

        // We use the rates from the privious time-step to reduce oscilations
        WellGroupHelpers::updateWellRates(fieldGroup, schedule(), reportStepIdx, this->prevWellState(), well_state);

        // Set ALQ for off-process wells to zero
        for (const auto& wname : schedule().wellNames(reportStepIdx)) {
            const bool is_producer = schedule().getWell(wname, reportStepIdx).isProducer();
            const bool not_on_this_process = well_state.wellMap().count(wname) == 0;
            if (is_producer && not_on_this_process) {
                well_state.setALQ(wname, 0.0);
            }
        }

        well_state.communicateGroupRates(comm);
        this->groupState().communicate_rates(comm);
        // compute wsolvent fraction for REIN wells
        updateWsolvent(fieldGroup, schedule(), reportStepIdx,  well_state_nupcol);

    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    updateWellTestState(const double& simulationTime, WellTestState& wellTestState) const
    {
        DeferredLogger local_deferredLogger;
        for (const auto& well : well_container_) {
            const auto wasClosed = wellTestState.hasWellClosed(well->name());

            well->updateWellTestState(this->wellState(), simulationTime, /*writeMessageToOPMLog=*/ true, wellTestState, local_deferredLogger);

            if (!wasClosed && wellTestState.hasWellClosed(well->name())) {
                this->closed_this_step_.insert(well->name());
            }
        }

        DeferredLogger global_deferredLogger = gatherDeferredLogger(local_deferredLogger);
        if (terminal_output_) {
            global_deferredLogger.logMessages();
        }
    }



    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    computeWellPotentials(std::vector<double>& well_potentials, const int reportStepIdx, DeferredLogger& deferred_logger)
    {
        // number of wells and phases
        const int nw = numLocalWells();
        const int np = numPhases();
        well_potentials.resize(nw * np, 0.0);

        auto well_state = this->wellState();

        const SummaryConfig& summaryConfig = ebosSimulator_.vanguard().summaryConfig();
        const bool write_restart_file = ebosSimulator_.vanguard().schedule().write_rst_file(reportStepIdx);
        auto exc_type = ExceptionType::NONE;
        std::string exc_msg;
        for (const auto& well : well_container_) {
            const bool needed_for_summary = ((summaryConfig.hasSummaryKey( "WWPI:" + well->name()) ||
                                              summaryConfig.hasSummaryKey( "WOPI:" + well->name()) ||
                                              summaryConfig.hasSummaryKey( "WGPI:" + well->name())) && well->isInjector()) ||
                    ((summaryConfig.hasSummaryKey( "WWPP:" + well->name()) ||
                      summaryConfig.hasSummaryKey( "WOPP:" + well->name()) ||
                      summaryConfig.hasSummaryKey( "WGPP:" + well->name())) && well->isProducer());

            bool needPotentialsForGuideRate = true;//eclWell.getGuideRatePhase() == Well::GuideRateTarget::UNDEFINED;
            if (write_restart_file || needed_for_summary || needPotentialsForGuideRate)
            {
                try {
                    std::vector<double> potentials;
                    well->computeWellPotentials(ebosSimulator_, well_state, potentials, deferred_logger);
                    // putting the sucessfully calculated potentials to the well_potentials
                    for (int p = 0; p < np; ++p) {
                        well_potentials[well->indexOfWell() * np + p] = std::abs(potentials[p]);
                    }
                } catch (const std::runtime_error& e) {
                    exc_type = ExceptionType::RUNTIME_ERROR;
                    exc_msg = e.what();
                } catch (const std::invalid_argument& e) {
                    exc_type = ExceptionType::INVALID_ARGUMENT;
                    exc_msg = e.what();
                } catch (const std::logic_error& e) {
                    exc_type = ExceptionType::LOGIC_ERROR;
                    exc_msg = e.what();
                } catch (const std::exception& e) {
                    exc_type = ExceptionType::DEFAULT;
                    exc_msg = e.what();
                }
            }
        }
        logAndCheckForExceptionsAndThrow(deferred_logger, exc_type,
                                         "computeWellPotentials() failed: " + exc_msg,
                                         terminal_output_);

        // Store it in the well state
        this->wellState().wellPotentials() = well_potentials;
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
            auto wellPtr = this->template createTypedWellPointer
                <StandardWell<TypeTag>>(shutWell, reportStepIdx);

            wellPtr->init(&this->phase_usage_, this->depth_, this->gravity_,
                          this->local_num_cells_, this->B_avg_);

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
        auto exc_type = ExceptionType::NONE;
        std::string exc_msg;
        try {
            for (const auto& well : well_container_) {
                well->checkWellOperability(ebosSimulator_, this->wellState(), deferred_logger);

                if (!well->isOperable() ) continue;

                auto& events = this->wellState().events(well->indexOfWell());
                if (events.hasEvent(WellStateFullyImplicitBlackoil::event_mask)) {
                    well->updateWellStateWithTarget(ebosSimulator_, this->wellState(), deferred_logger);
                    // There is no new well control change input within a report step,
                    // so next time step, the well does not consider to have effective events anymore.
                    events.clearEvent(WellStateFullyImplicitBlackoil::event_mask);
                }

                // solve the well equation initially to improve the initial solution of the well model
                if (param_.solve_welleq_initially_) {
                    well->solveWellEquation(ebosSimulator_, this->wellState(), this->groupState(), deferred_logger);
                }

             }  // end of for (const auto& well : well_container_)
            updatePrimaryVariables(deferred_logger);
        } catch (const std::runtime_error& e) {
            exc_type = ExceptionType::RUNTIME_ERROR;
            exc_msg = e.what();
        } catch (const std::invalid_argument& e) {
            exc_type = ExceptionType::INVALID_ARGUMENT;
            exc_msg = e.what();
        } catch (const std::logic_error& e) {
            exc_type = ExceptionType::LOGIC_ERROR;
            exc_msg = e.what();
        } catch (const std::exception& e) {
            exc_type = ExceptionType::DEFAULT;
            exc_msg = e.what();
        }
        logAndCheckForExceptionsAndThrow(deferred_logger, exc_type, "prepareTimestep() failed: " + exc_msg, terminal_output_);
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    calculateEfficiencyFactors(const int reportStepIdx)
    {
        if ( !localWellsActive() ) {
            return;
        }

        for (auto& well : well_container_) {
            const Well& wellEcl = well->wellEcl();
            double well_efficiency_factor = wellEcl.getEfficiencyFactor();
            WellGroupHelpers::accumulateGroupEfficiencyFactor(schedule().getGroup(wellEcl.groupName(), reportStepIdx), schedule(), reportStepIdx, well_efficiency_factor);
            well->setWellEfficiencyFactor(well_efficiency_factor);
        }
    }



    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    setupCartesianToCompressed_(const int* global_cell, int number_of_cartesian_cells)
    {
        cartesian_to_compressed_.resize(number_of_cartesian_cells, -1);
        if (global_cell) {
            auto elemIt = ebosSimulator_.gridView().template begin</*codim=*/ 0>();
            for (unsigned i = 0; i < local_num_cells_; ++i) {
                // Skip perforations in the overlap/ghost for distributed wells.
                if (elemIt->partitionType() == Dune::InteriorEntity)
                {
                    assert(ebosSimulator_.gridView().indexSet().index(*elemIt) == static_cast<int>(i));
                    cartesian_to_compressed_[global_cell[i]] = i;
                }
                ++elemIt;
            }
        }
        else {
            for (unsigned i = 0; i < local_num_cells_; ++i) {
                cartesian_to_compressed_[i] = i;
            }
        }

    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    setRepRadiusPerfLength()
    {
        for (const auto& well : well_container_) {
            well->setRepRadiusPerfLength(cartesian_to_compressed_);
        }
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
        const auto& elemEndIt = gridView.template end</*codim=*/0, Dune::Interior_Partition>();

        for (auto elemIt = gridView.template begin</*codim=*/0, Dune::Interior_Partition>();
             elemIt != elemEndIt; ++elemIt)
        {
            elemCtx.updatePrimaryStencil(*elemIt);
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

        // compute global average
        grid.comm().sum(B_avg.data(), B_avg.size());
        for(auto& bval: B_avg)
        {
            bval/=global_num_cells_;
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
        if (wellsActive()  && numPhases() < 3) {
            return numPhases();
        }
        int numComp = FluidSystem::numComponents;
        if constexpr (has_solvent_) {
            numComp ++;
        }

        return numComp;
    }

    template<typename TypeTag>
    int
    BlackoilWellModel<TypeTag>:: numLocalWells() const
    {
        return wells_ecl_.size();
    }

    template<typename TypeTag>
    int
    BlackoilWellModel<TypeTag>::numPhases() const
    {
        return phase_usage_.num_phases;
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::extractLegacyDepth_()
    {
        const auto& grid = ebosSimulator_.vanguard().grid();
        const unsigned numCells = grid.size(/*codim=*/0);

        depth_.resize(numCells);
        for (unsigned cellIdx = 0; cellIdx < numCells; ++cellIdx) {
            depth_[cellIdx] = UgGridHelpers::cellCenterDepth( grid, cellIdx );
        }
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    updatePerforationIntensiveQuantities() {
        ElementContext elemCtx(ebosSimulator_);
        const auto& gridView = ebosSimulator_.gridView();
        const auto& elemEndIt = gridView.template end</*codim=*/0, Dune::Interior_Partition>();
        for (auto elemIt = gridView.template begin</*codim=*/0, Dune::Interior_Partition>();
             elemIt != elemEndIt;
             ++elemIt)
        {

            elemCtx.updatePrimaryStencil(*elemIt);
            int elemIdx = elemCtx.globalSpaceIndex(0, 0);

            if (!is_cell_perforated_[elemIdx]) {
                continue;
            }
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
        }
    }


    template<typename TypeTag>
    bool
    BlackoilWellModel<TypeTag>::
    hasWell(const std::string& wname) {
        auto iter = std::find_if( this->wells_ecl_.begin(), this->wells_ecl_.end(), [&wname](const Well& well) { return well.name() == wname; });
        return (iter != this->wells_ecl_.end());
    }


    // convert well data from opm-common to well state from opm-core
    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    loadRestartData( const data::Wells& rst_wells,
                     const data::GroupAndNetworkValues& grpNwrkValues,
                     const PhaseUsage& phases,
                     const bool handle_ms_well,
                     WellStateFullyImplicitBlackoil& well_state)
    {
        using GPMode = Group::ProductionCMode;
        using GIMode = Group::InjectionCMode;

        using rt = data::Rates::opt;
        const auto np = phases.num_phases;

        std::vector< rt > phs( np );
        if( phases.phase_used[BlackoilPhases::Aqua] ) {
            phs.at( phases.phase_pos[BlackoilPhases::Aqua] ) = rt::wat;
        }

        if( phases.phase_used[BlackoilPhases::Liquid] ) {
            phs.at( phases.phase_pos[BlackoilPhases::Liquid] ) = rt::oil;
        }

        if( phases.phase_used[BlackoilPhases::Vapour] ) {
            phs.at( phases.phase_pos[BlackoilPhases::Vapour] ) = rt::gas;
        }

        for( const auto& wm : well_state.wellMap() ) {
            const auto well_index = wm.second[ 0 ];
            const auto& rst_well = rst_wells.at( wm.first );
            well_state.bhp()[ well_index ] = rst_well.bhp;
            well_state.temperature()[ well_index ] = rst_well.temperature;

            if (rst_well.current_control.isProducer) {
                well_state.currentProductionControl( well_index, rst_well.current_control.prod);
            }
            else {
                well_state.currentInjectionControl( well_index, rst_well.current_control.inj);
            }

            const auto wellrate_index = well_index * np;
            for( size_t i = 0; i < phs.size(); ++i ) {
                assert( rst_well.rates.has( phs[ i ] ) );
                well_state.wellRates()[ wellrate_index + i ] = rst_well.rates.get( phs[ i ] );
            }

            auto& perf_pressure = well_state.perfPress(well_index);
            auto& perf_rates = well_state.perfRates(well_index);
            auto * perf_phase_rates = well_state.mutable_perfPhaseRates().data() + wm.second[1]*np;
            const auto& perf_data = this->well_perf_data_[well_index];

            for (std::size_t perf_index = 0; perf_index < perf_data.size(); perf_index++) {
                const auto& pd = perf_data[perf_index];
                const auto& rst_connection = rst_well.connections[pd.ecl_index];
                perf_pressure[perf_index] = rst_connection.pressure;
                perf_rates[perf_index] = rst_connection.reservoir_rate;
                for (int phase_index = 0; phase_index < np; ++phase_index)
                    perf_phase_rates[perf_index*np + phase_index] = rst_connection.rates.get(phs[phase_index]);
            }

            if (handle_ms_well && !rst_well.segments.empty()) {
                // we need the well_ecl_ information
                const std::string& well_name = wm.first;
                const Well& well_ecl = getWellEcl(well_name);

                const WellSegments& segment_set = well_ecl.getSegments();

                const int top_segment_index = well_state.topSegmentIndex(well_index);
                const auto& segments = rst_well.segments;

                // \Note: eventually we need to hanlde the situations that some segments are shut
                assert(0u + segment_set.size() == segments.size());

                for (const auto& segment : segments) {
                    const int segment_index = segment_set.segmentNumberToIndex(segment.first);

                    // recovering segment rates and pressure from the restart values
                    const auto pres_idx = data::SegmentPressures::Value::Pressure;
                    well_state.segPress()[top_segment_index + segment_index] = segment.second.pressures[pres_idx];

                    const auto& segment_rates = segment.second.rates;
                    for (int p = 0; p < np; ++p) {
                        well_state.segRates()[(top_segment_index + segment_index) * np + p] = segment_rates.get(phs[p]);
                    }
                }
            }
        }

        for (const auto& [group, value] : grpNwrkValues.groupData) {
            const auto cpc = value.currentControl.currentProdConstraint;
            const auto cgi = value.currentControl.currentGasInjectionConstraint;
            const auto cwi = value.currentControl.currentWaterInjectionConstraint;

            if (cpc != GPMode::NONE) {
                this->groupState().production_control(group, cpc);
            }

            if (cgi != GIMode::NONE) {
                this->groupState().injection_control(group, Phase::GAS, cgi);
            }

            if (cwi != GIMode::NONE) {
                this->groupState().injection_control(group, Phase::WATER, cwi);
            }
        }
    }





    template<typename TypeTag>
    bool
    BlackoilWellModel<TypeTag>::
    anyMSWellOpenLocal() const
    {
        for (const auto& well : wells_ecl_) {
            if (well.isMultiSegment()) {
                return true;
            }
        }
        return false;
    }





    template<typename TypeTag>
    const Well&
    BlackoilWellModel<TypeTag>::
    getWellEcl(const std::string& well_name) const
    {
        // finding the iterator of the well in wells_ecl
        auto well_ecl = std::find_if(wells_ecl_.begin(),
                                     wells_ecl_.end(),
                                     [&well_name](const Well& elem)->bool {
                                         return elem.name() == well_name;
                                     });

        assert(well_ecl != wells_ecl_.end());

        return *well_ecl;
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
    void
    BlackoilWellModel<TypeTag>::
    updateGroupIndividualControls(DeferredLogger& deferred_logger, std::set<std::string>& switched_groups)
    {
        const int reportStepIdx = ebosSimulator_.episodeIndex();

        const int nupcol = schedule()[reportStepIdx].nupcol();
        const int iterationIdx = ebosSimulator_.model().newtonMethod().numIterations();
        // don't switch group control when iterationIdx > nupcol
        // to avoid oscilations between group controls
        if (iterationIdx > nupcol)
            return;

        const Group& fieldGroup = schedule().getGroup("FIELD", reportStepIdx);
        updateGroupIndividualControl(fieldGroup, deferred_logger, switched_groups);
    }



    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    updateGroupIndividualControl(const Group& group, DeferredLogger& deferred_logger, std::set<std::string>& switched_groups) {

        const int reportStepIdx = ebosSimulator_.episodeIndex();
        const bool skip = switched_groups.count(group.name());
        if (!skip && group.isInjectionGroup())
        {
            const Phase all[] = {Phase::WATER, Phase::OIL, Phase::GAS};
            for (Phase phase : all) {
                if (!group.hasInjectionControl(phase)) {
                    continue;
                }
                Group::InjectionCMode newControl = checkGroupInjectionConstraints(group, phase);
                if (newControl != Group::InjectionCMode::NONE)
                {
                    switched_groups.insert(group.name());
                    actionOnBrokenConstraints(group, newControl, phase, deferred_logger);
                }
            }
        }
        if (!skip && group.isProductionGroup()) {
            Group::ProductionCMode newControl = checkGroupProductionConstraints(group, deferred_logger);
            const auto& summaryState = ebosSimulator_.vanguard().summaryState();
            const auto controls = group.productionControls(summaryState);
            if (newControl != Group::ProductionCMode::NONE)
            {
                switched_groups.insert(group.name());
                actionOnBrokenConstraints(group, controls.exceed_action, newControl, deferred_logger);
            }
        }

        // call recursively down the group hiearchy
        for (const std::string& groupName : group.groups()) {
            updateGroupIndividualControl( schedule().getGroup(groupName, reportStepIdx), deferred_logger, switched_groups);
        }
    }

    template<typename TypeTag>
    bool
    BlackoilWellModel<TypeTag>::
    checkGroupConstraints(const Group& group, DeferredLogger& deferred_logger) const {

        const int reportStepIdx = ebosSimulator_.episodeIndex();
        if (group.isInjectionGroup()) {
            const Phase all[] = {Phase::WATER, Phase::OIL, Phase::GAS};
            for (Phase phase : all) {
                if (!group.hasInjectionControl(phase)) {
                    continue;
                }
                Group::InjectionCMode newControl = checkGroupInjectionConstraints(group, phase);
                if (newControl != Group::InjectionCMode::NONE) {
                    return true;
                }
            }
        }
        if (group.isProductionGroup()) {
            Group::ProductionCMode newControl = checkGroupProductionConstraints(group, deferred_logger);
            if (newControl != Group::ProductionCMode::NONE)
            {
                return true;
            }
        }

        // call recursively down the group hiearchy
        bool violated = false;
        for (const std::string& groupName : group.groups()) {
            violated = violated || checkGroupConstraints( schedule().getGroup(groupName, reportStepIdx), deferred_logger);
        }
        return violated;
    }



    template<typename TypeTag>
    Group::ProductionCMode
    BlackoilWellModel<TypeTag>::
    checkGroupProductionConstraints(const Group& group, DeferredLogger& deferred_logger) const {

        const int reportStepIdx = ebosSimulator_.episodeIndex();
        const auto& summaryState = ebosSimulator_.vanguard().summaryState();
        const auto& comm = ebosSimulator_.vanguard().grid().comm();
        const auto& well_state = this->wellState();

        const auto controls = group.productionControls(summaryState);
        const Group::ProductionCMode& currentControl = this->groupState().production_control(group.name());

        if (group.has_control(Group::ProductionCMode::ORAT))
        {
            if (currentControl != Group::ProductionCMode::ORAT)
            {
                double current_rate = 0.0;
                current_rate += WellGroupHelpers::sumWellRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Liquid], false);

                // sum over all nodes
                current_rate = comm.sum(current_rate);

                if (controls.oil_target < current_rate  ) {
                    return Group::ProductionCMode::ORAT;
                }
            }
        }

        if (group.has_control(Group::ProductionCMode::WRAT))
        {
            if (currentControl != Group::ProductionCMode::WRAT)
            {

                double current_rate = 0.0;
                current_rate += WellGroupHelpers::sumWellRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Aqua], false);

                // sum over all nodes
                current_rate = comm.sum(current_rate);

                if (controls.water_target < current_rate  ) {
                    return Group::ProductionCMode::WRAT;
                }
            }
        }
        if (group.has_control(Group::ProductionCMode::GRAT))
        {
            if (currentControl != Group::ProductionCMode::GRAT)
            {
                double current_rate = 0.0;
                current_rate += WellGroupHelpers::sumWellRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Vapour], false);

                // sum over all nodes
                current_rate = comm.sum(current_rate);
                if (controls.gas_target < current_rate  ) {
                    return Group::ProductionCMode::GRAT;
                }
            }
        }
        if (group.has_control(Group::ProductionCMode::LRAT))
        {
            if (currentControl != Group::ProductionCMode::LRAT)
            {
                double current_rate = 0.0;
                current_rate += WellGroupHelpers::sumWellRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Liquid], false);
                current_rate += WellGroupHelpers::sumWellRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Aqua], false);

                // sum over all nodes
                current_rate = comm.sum(current_rate);

                if (controls.liquid_target < current_rate  ) {
                     return Group::ProductionCMode::LRAT;
                }
            }
        }

        if (group.has_control(Group::ProductionCMode::CRAT))
        {
            OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "CRAT control for production groups not implemented" , deferred_logger);
        }
        if (group.has_control(Group::ProductionCMode::RESV))
        {
            if (currentControl != Group::ProductionCMode::RESV)
            {
                double current_rate = 0.0;
                current_rate += WellGroupHelpers::sumWellResRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Aqua], true);
                current_rate += WellGroupHelpers::sumWellResRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Liquid], true);
                current_rate += WellGroupHelpers::sumWellResRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Vapour], true);

                // sum over all nodes
                current_rate = comm.sum(current_rate);

                if (controls.resv_target < current_rate  ) {
                    return Group::ProductionCMode::RESV;
                }
            }

        }
        if (group.has_control(Group::ProductionCMode::PRBL))
        {
            OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "PRBL control for production groups not implemented", deferred_logger);
        }
        return Group::ProductionCMode::NONE;
    }


    template<typename TypeTag>
    Group::InjectionCMode
    BlackoilWellModel<TypeTag>::
    checkGroupInjectionConstraints(const Group& group, const Phase& phase) const {

        const int reportStepIdx = ebosSimulator_.episodeIndex();
        const auto& summaryState = ebosSimulator_.vanguard().summaryState();
        const auto& comm = ebosSimulator_.vanguard().grid().comm();
        const auto& well_state = this->wellState();

        int phasePos;
        if (phase == Phase::GAS && phase_usage_.phase_used[BlackoilPhases::Vapour] )
            phasePos = phase_usage_.phase_pos[BlackoilPhases::Vapour];
        else if (phase == Phase::OIL && phase_usage_.phase_used[BlackoilPhases::Liquid])
            phasePos = phase_usage_.phase_pos[BlackoilPhases::Liquid];
        else if (phase == Phase::WATER && phase_usage_.phase_used[BlackoilPhases::Aqua] )
            phasePos = phase_usage_.phase_pos[BlackoilPhases::Aqua];
        else
            OPM_THROW(std::runtime_error, "Unknown phase" );

        const auto& controls = group.injectionControls(phase, summaryState);
        auto currentControl = this->groupState().injection_control(group.name(), phase);

        if (controls.has_control(Group::InjectionCMode::RATE))
        {
            if (currentControl != Group::InjectionCMode::RATE)
            {
                double current_rate = 0.0;
                current_rate += WellGroupHelpers::sumWellRates(group, schedule(), well_state, reportStepIdx, phasePos, /*isInjector*/true);

                // sum over all nodes
                current_rate = comm.sum(current_rate);

                if (controls.surface_max_rate < current_rate) {
                    return Group::InjectionCMode::RATE;
                }
            }
        }
        if (controls.has_control(Group::InjectionCMode::RESV))
        {
            if (currentControl != Group::InjectionCMode::RESV)
            {
                double current_rate = 0.0;
                current_rate += WellGroupHelpers::sumWellResRates(group, schedule(), well_state, reportStepIdx, phasePos, /*isInjector*/true);
                // sum over all nodes
                current_rate = comm.sum(current_rate);

                if (controls.resv_max_rate < current_rate) {
                    return Group::InjectionCMode::RESV;
                }
            }
        }
        if (controls.has_control(Group::InjectionCMode::REIN))
        {
            if (currentControl != Group::InjectionCMode::REIN)
            {
                double production_Rate = 0.0;
                const Group& groupRein = schedule().getGroup(controls.reinj_group, reportStepIdx);
                production_Rate += WellGroupHelpers::sumWellRates(groupRein, schedule(), well_state, reportStepIdx, phasePos, /*isInjector*/false);

                // sum over all nodes
                production_Rate = comm.sum(production_Rate);

                double current_rate = 0.0;
                current_rate += WellGroupHelpers::sumWellRates(group, schedule(), well_state, reportStepIdx, phasePos, /*isInjector*/true);

                // sum over all nodes
                current_rate = comm.sum(current_rate);

                if (controls.target_reinj_fraction*production_Rate < current_rate) {
                    return Group::InjectionCMode::REIN;
                }
            }
        }
        if (controls.has_control(Group::InjectionCMode::VREP))
        {
            if (currentControl != Group::InjectionCMode::VREP)
            {
                double voidage_rate = 0.0;
                const Group& groupVoidage = schedule().getGroup(controls.voidage_group, reportStepIdx);
                voidage_rate += WellGroupHelpers::sumWellResRates(groupVoidage, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Aqua], false);
                voidage_rate += WellGroupHelpers::sumWellResRates(groupVoidage, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Liquid], false);
                voidage_rate += WellGroupHelpers::sumWellResRates(groupVoidage, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Vapour], false);

                // sum over all nodes
                voidage_rate = comm.sum(voidage_rate);

                double total_rate = 0.0;
                total_rate += WellGroupHelpers::sumWellResRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Aqua], true);
                total_rate += WellGroupHelpers::sumWellResRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Liquid], true);
                total_rate += WellGroupHelpers::sumWellResRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Vapour], true);

                // sum over all nodes
                total_rate = comm.sum(total_rate);

                if (controls.target_void_fraction*voidage_rate < total_rate) {
                    return Group::InjectionCMode::VREP;
                }
            }
        }
        return Group::InjectionCMode::NONE;
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    checkGconsaleLimits(const Group& group, WellState& well_state, DeferredLogger& deferred_logger)
    {
        const int reportStepIdx = ebosSimulator_.episodeIndex();
         // call recursively down the group hiearchy
        for (const std::string& groupName : group.groups()) {
            checkGconsaleLimits( schedule().getGroup(groupName, reportStepIdx), well_state, deferred_logger);
        }

        // only for groups with gas injection controls
        if (!group.hasInjectionControl(Phase::GAS)) {
            return;
        }

        // check if gconsale is used for this group
        if (!schedule()[reportStepIdx].gconsale().has(group.name()))
            return;

        std::ostringstream ss;

        const auto& summaryState = ebosSimulator_.vanguard().summaryState();
        const auto& comm = ebosSimulator_.vanguard().grid().comm();

        const auto& gconsale = schedule()[reportStepIdx].gconsale().get(group.name(), summaryState);
        const Group::ProductionCMode& oldProductionControl = this->groupState().production_control(group.name());


        int gasPos = phase_usage_.phase_pos[BlackoilPhases::Vapour];
        double production_rate = WellGroupHelpers::sumWellRates(group, schedule(), well_state, reportStepIdx, gasPos, /*isInjector*/false);
        double injection_rate = WellGroupHelpers::sumWellRates(group, schedule(), well_state, reportStepIdx, gasPos, /*isInjector*/true);

        // sum over all nodes
        injection_rate = comm.sum(injection_rate);
        production_rate = comm.sum(production_rate);

        double sales_rate = production_rate - injection_rate;
        double production_target = gconsale.sales_target + injection_rate;

        // add import rate and substract consumption rate for group for gas
        if (schedule()[reportStepIdx].gconsump().has(group.name())) {
            const auto& gconsump = schedule()[reportStepIdx].gconsump().get(group.name(), summaryState);
            if (phase_usage_.phase_used[BlackoilPhases::Vapour]) {
                sales_rate += gconsump.import_rate;
                sales_rate -= gconsump.consumption_rate;
                production_target -= gconsump.import_rate;
                production_target += gconsump.consumption_rate;
            }
        }
        if (sales_rate > gconsale.max_sales_rate) {

            switch(gconsale.max_proc) {
            case GConSale::MaxProcedure::NONE: {
                if (oldProductionControl != Group::ProductionCMode::GRAT && oldProductionControl != Group::ProductionCMode::NONE) {
                    ss << "Group sales exceed maximum limit, but the action is NONE for " + group.name() + ". Nothing happens";
                }
                break;
                }
            case GConSale::MaxProcedure::CON: {
                OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "GCONSALE exceed limit CON not implemented", deferred_logger);
                break;
            }
            case GConSale::MaxProcedure::CON_P: {
                OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "GCONSALE exceed limit CON_P not implemented", deferred_logger);
                break;
            }
            case GConSale::MaxProcedure::WELL: {
                OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "GCONSALE exceed limit WELL not implemented", deferred_logger);
                break;
            }
            case GConSale::MaxProcedure::PLUG: {
                OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "GCONSALE exceed limit PLUG not implemented", deferred_logger);
                break;
            }
            case GConSale::MaxProcedure::MAXR: {
                OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "GCONSALE exceed limit MAXR not implemented", deferred_logger);
                break;
            }
            case GConSale::MaxProcedure::END: {
                OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "GCONSALE exceed limit END not implemented", deferred_logger);
                break;
            }
            case GConSale::MaxProcedure::RATE: {
                this->groupState().production_control(group.name(), Group::ProductionCMode::GRAT);
                ss << "Maximum GCONSALE limit violated for " << group.name() << ". The group is switched from ";
                ss << Group::ProductionCMode2String(oldProductionControl) << " to " << Group::ProductionCMode2String(Group::ProductionCMode::GRAT);
                ss << " and limited by the maximum sales rate after consumption and import are considered" ;
                this->groupState().update_grat_sales_target(group.name(), production_target);
                break;
            }
            default:
                throw("Invalid procedure for maximum rate limit selected for group" + group.name());
            }
        }
        if (sales_rate < gconsale.min_sales_rate) {
            const Group::ProductionCMode& currentProductionControl = this->groupState().production_control(group.name());
            if ( currentProductionControl == Group::ProductionCMode::GRAT ) {
                ss << "Group " + group.name() + " has sale rate less then minimum permitted value and is under GRAT control. \n";
                ss << "The GRAT is increased to meet the sales minimum rate. \n";
                this->groupState().update_grat_sales_target(group.name(), production_target);
            //} else if () {//TODO add action for WGASPROD
            //} else if () {//TODO add action for drilling queue
            } else {
                ss << "Group " + group.name() + " has sale rate less then minimum permitted value but cannot increase the group production rate \n";
                ss << "or adjust gas production using WGASPROD or drill new wells to meet the sales target. \n";
                ss << "Note that WGASPROD and drilling queues are not implemented in Flow. No action is taken. \n ";
            }
        }
        if (gconsale.sales_target < 0.0) {
            OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + " has sale rate target less then zero. Not implemented in Flow" , deferred_logger);
        }

        auto cc = Dune::MPIHelper::getCollectiveCommunication();
        if (!ss.str().empty() && cc.rank() == 0)
            deferred_logger.info(ss.str());

    }


    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    actionOnBrokenConstraints(const Group& group, const Group::ExceedAction& exceed_action, const Group::ProductionCMode& newControl, DeferredLogger& deferred_logger) {

        const Group::ProductionCMode oldControl = this->groupState().production_control(group.name());

        std::ostringstream ss;

        switch(exceed_action) {
        case Group::ExceedAction::NONE: {
            if (oldControl != newControl && oldControl != Group::ProductionCMode::NONE) {
                ss << "Group production exceed action is NONE for group " + group.name() + ". Nothing happens.";
            }
            break;
        }
        case Group::ExceedAction::CON: {
            OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "GroupProductionExceedLimit CON not implemented", deferred_logger);
            break;
        }
        case Group::ExceedAction::CON_PLUS: {
            OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "GroupProductionExceedLimit CON_PLUS not implemented", deferred_logger);
            break;
        }
        case Group::ExceedAction::WELL: {
            OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "GroupProductionExceedLimit WELL not implemented", deferred_logger);
            break;
        }
        case Group::ExceedAction::PLUG: {
            OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "GroupProductionExceedLimit PLUG not implemented", deferred_logger);
            break;
        }
        case Group::ExceedAction::RATE: {
            if (oldControl != newControl) {
                this->groupState().production_control(group.name(), newControl);
                ss << "Switching production control mode for group "<< group.name()
                   << " from " << Group::ProductionCMode2String(oldControl)
                   << " to " << Group::ProductionCMode2String(newControl);
            }
            break;
        }
        default:
            throw("Invalid procedure for maximum rate limit selected for group" + group.name());
        }

        auto cc = Dune::MPIHelper::getCollectiveCommunication();
        if (!ss.str().empty() && cc.rank() == 0)
            deferred_logger.info(ss.str());


    }


    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    actionOnBrokenConstraints(const Group& group, const Group::InjectionCMode& newControl, const Phase& controlPhase, DeferredLogger& deferred_logger) {
        auto oldControl = this->groupState().injection_control(group.name(), controlPhase);

        std::ostringstream ss;
        if (oldControl != newControl) {
            const std::string from = Group::InjectionCMode2String(oldControl);
            ss << "Switching injection control mode for group "<< group.name()
               << " from " << Group::InjectionCMode2String(oldControl)
               << " to " << Group::InjectionCMode2String(newControl);
            this->groupState().injection_control(group.name(), controlPhase, newControl);
        }
        auto cc = Dune::MPIHelper::getCollectiveCommunication();
        if (!ss.str().empty() && cc.rank() == 0)
            deferred_logger.info(ss.str());

    }




    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    updateGroupHigherControls(DeferredLogger& deferred_logger, std::set<std::string>& switched_groups)
    {
        const int reportStepIdx = ebosSimulator_.episodeIndex();
        const Group& fieldGroup = schedule().getGroup("FIELD", reportStepIdx);
        checkGroupHigherConstraints(fieldGroup, deferred_logger, switched_groups);
    }


    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    checkGroupHigherConstraints(const Group& group, DeferredLogger& deferred_logger, std::set<std::string>& switched_groups)
    {
        // Set up coefficients for RESV <-> surface rate conversion.
        // Use the pvtRegionIdx from the top cell of the first well.
        // TODO fix this!
        // This is only used for converting RESV rates.
        // What is the proper approach?
        const auto& comm = ebosSimulator_.vanguard().grid().comm();
        const int fipnum = 0;
        int pvtreg = well_perf_data_.empty() || well_perf_data_[0].empty()
            ? pvt_region_idx_[0]
            : pvt_region_idx_[well_perf_data_[0][0].cell_index];

        if ( comm.size() > 1)
        {
            // Just like in the sequential case the pvtregion is determined
            // by the first cell of the first well. What is the first well
            // is decided by the order in the Schedule using Well::seqIndex()
            int firstWellIndex = well_perf_data_.empty() ?
                std::numeric_limits<int>::max() : wells_ecl_[0].seqIndex();
            auto regIndexPair = std::make_pair(pvtreg, firstWellIndex);
            std::vector<decltype(regIndexPair)> pairs(comm.size());
            comm.allgather(&regIndexPair, 1, pairs.data());
            pvtreg = std::min_element(pairs.begin(), pairs.end(),
                                      [](const auto& p1, const auto& p2){ return p1.second < p2.second;})
                ->first;
        }

        std::vector<double> resv_coeff(phase_usage_.num_phases, 0.0);
        rateConverter_->calcCoeff(fipnum, pvtreg, resv_coeff);

        const int reportStepIdx = ebosSimulator_.episodeIndex();
        const auto& summaryState = ebosSimulator_.vanguard().summaryState();

        std::vector<double> rates(phase_usage_.num_phases, 0.0);

        const bool skip = switched_groups.count(group.name()) || group.name() == "FIELD";

        if (!skip && group.isInjectionGroup()) {
            // Obtain rates for group.
            for (int phasePos = 0; phasePos < phase_usage_.num_phases; ++phasePos) {
                const double local_current_rate = WellGroupHelpers::sumWellRates(group, schedule(), this->wellState(), reportStepIdx, phasePos, /* isInjector */ true);
                // Sum over all processes
                rates[phasePos] = comm.sum(local_current_rate);
            }
            const Phase all[] = { Phase::WATER, Phase::OIL, Phase::GAS };
            for (Phase phase : all) {
                // Check higher up only if under individual (not FLD) control.
                auto currentControl = this->groupState().injection_control(group.name(), phase);
                if (currentControl != Group::InjectionCMode::FLD && group.injectionGroupControlAvailable(phase)) {
                    const Group& parentGroup = schedule().getGroup(group.parent(), reportStepIdx);
                    const std::pair<bool, double> changed = WellGroupHelpers::checkGroupConstraintsInj(
                        group.name(),
                        group.parent(),
                        parentGroup,
                        this->wellState(),
                        this->groupState(),
                        reportStepIdx,
                        guideRate_.get(),
                        rates.data(),
                        phase,
                        phase_usage_,
                        group.getGroupEfficiencyFactor(),
                        schedule(),
                        summaryState,
                        resv_coeff,
                        deferred_logger);
                    if (changed.first) {
                        switched_groups.insert(group.name());
                        actionOnBrokenConstraints(group, Group::InjectionCMode::FLD, phase, deferred_logger);
                    }
                }
            }
        }

        if (!skip && group.isProductionGroup()) {
            // Obtain rates for group.
            for (int phasePos = 0; phasePos < phase_usage_.num_phases; ++phasePos) {
                const double local_current_rate = WellGroupHelpers::sumWellRates(group, schedule(), this->wellState(), reportStepIdx, phasePos, /* isInjector */ false);
                // Sum over all processes
                rates[phasePos] = -comm.sum(local_current_rate);
            }
            // Check higher up only if under individual (not FLD) control.
            const Group::ProductionCMode& currentControl = this->groupState().production_control(group.name());
                if (currentControl != Group::ProductionCMode::FLD && group.productionGroupControlAvailable()) {
                    const Group& parentGroup = schedule().getGroup(group.parent(), reportStepIdx);
                    const std::pair<bool, double> changed = WellGroupHelpers::checkGroupConstraintsProd(
                        group.name(),
                        group.parent(),
                        parentGroup,
                        this->wellState(),
                        this->groupState(),
                        reportStepIdx,
                        guideRate_.get(),
                        rates.data(),
                        phase_usage_,
                        group.getGroupEfficiencyFactor(),
                        schedule(),
                        summaryState,
                        resv_coeff,
                        deferred_logger);
                    if (changed.first) {
                        switched_groups.insert(group.name());
                        const auto exceed_action = group.productionControls(summaryState).exceed_action;
                        actionOnBrokenConstraints(group, exceed_action, Group::ProductionCMode::FLD, deferred_logger);
                    }
                }
        }

        // call recursively down the group hiearchy
        for (const std::string& groupName : group.groups()) {
            checkGroupHigherConstraints( schedule().getGroup(groupName, reportStepIdx), deferred_logger, switched_groups);
         }
    }








    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    updateEclWells(const int timeStepIdx, const std::unordered_set<std::string>& wells) {
        const auto& schedule = this->ebosSimulator_.vanguard().schedule();
        for (const auto& wname : wells) {
            auto well_iter = std::find_if( this->wells_ecl_.begin(), this->wells_ecl_.end(), [wname] (const auto& well) -> bool { return well.name() == wname;});
            if (well_iter != this->wells_ecl_.end()) {
                auto well_index = std::distance( this->wells_ecl_.begin(), well_iter );
                this->wells_ecl_[well_index] = schedule.getWell(wname, timeStepIdx);

                const auto& well = this->wells_ecl_[well_index];
                auto& pd     = this->well_perf_data_[well_index];
                auto  pdIter = pd.begin();
                for (const auto& conn : well.getConnections()) {
                    if (conn.state() != Connection::State::SHUT) {
                        pdIter->connection_transmissibility_factor = conn.CF();
                        ++pdIter;
                    }
                }
                this->wellState().updateStatus(well_index, well.getStatus());
                this->wellState().resetConnectionTransFactors(well_index, pd);
                this->prod_index_calc_[well_index].reInit(well);
            }
        }
    }



    template<typename TypeTag>
    double
    BlackoilWellModel<TypeTag>::
    wellPI(const int well_index) const
    {
        const auto& pu = this->phase_usage_;
        const auto  np = this->numPhases();
        const auto* pi = &this->wellState().productivityIndex()[np*well_index + 0];

        const auto preferred = this->wells_ecl_[well_index].getPreferredPhase();
        switch (preferred) { // Should really have LIQUID = OIL + WATER here too...
        case Phase::WATER:
            return pu.phase_used[BlackoilPhases::PhaseIndex::Aqua]
                ? pi[pu.phase_pos[BlackoilPhases::PhaseIndex::Aqua]]
                : 0.0;

        case Phase::OIL:
            return pu.phase_used[BlackoilPhases::PhaseIndex::Liquid]
                ? pi[pu.phase_pos[BlackoilPhases::PhaseIndex::Liquid]]
                : 0.0;

        case Phase::GAS:
            return pu.phase_used[BlackoilPhases::PhaseIndex::Vapour]
                ? pi[pu.phase_pos[BlackoilPhases::PhaseIndex::Vapour]]
                : 0.0;

        default:
            throw std::invalid_argument {
                "Unsupported preferred phase " +
                std::to_string(static_cast<int>(preferred))
            };
        }
    }





    template<typename TypeTag>
    double
    BlackoilWellModel<TypeTag>::
    wellPI(const std::string& well_name) const
    {
        auto well_iter = std::find_if(this->wells_ecl_.begin(), this->wells_ecl_.end(),
            [&well_name](const Well& well)
        {
            return well.name() == well_name;
        });

        if (well_iter == this->wells_ecl_.end()) {
            throw std::logic_error { "Could not find well: " + well_name };
        }

        auto well_index = std::distance(this->wells_ecl_.begin(), well_iter);
        return this->wellPI(well_index);
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
    runWellPIScaling(const int timeStepIdx, DeferredLogger& local_deferredLogger)
    {
        if (this->last_run_wellpi_.has_value() && (*this->last_run_wellpi_ == timeStepIdx)) {
            // We've already run WELPI scaling for this report step.  Most
            // common for the very first report step.  Don't redo WELPI scaling.
            return;
        }

        auto hasWellPIEvent = [this, timeStepIdx](const int well_index) -> bool
        {
            return this->schedule()[timeStepIdx].wellgroup_events()
                .hasEvent(this->wells_ecl_[well_index].name(),
                          ScheduleEvents::Events::WELL_PRODUCTIVITY_INDEX);
        };

        auto updateEclWell = [this, timeStepIdx](const int well_index) -> void
        {
            const auto& schedule = this->schedule();
            const auto& wname = this->wells_ecl_[well_index].name();
            this->wells_ecl_[well_index] = schedule.getWell(wname, timeStepIdx);

            const auto& well = this->wells_ecl_[well_index];
            auto& pd     = this->well_perf_data_[well_index];
            auto  pdIter = pd.begin();
            for (const auto& conn : well.getConnections()) {
                if (conn.state() != Connection::State::SHUT) {
                    pdIter->connection_transmissibility_factor = conn.CF();
                    ++pdIter;
                }
            }
            this->wellState().resetConnectionTransFactors(well_index, pd);
            this->prod_index_calc_[well_index].reInit(well);
        };


        auto rescaleWellPI =
            [this, timeStepIdx](const int    well_index,
                                const double newWellPI) -> void
        {
            const auto& wname = this->wells_ecl_[well_index].name();

            auto& schedule = this->ebosSimulator_.vanguard().schedule(); // Mutable
            schedule.applyWellProdIndexScaling(wname, timeStepIdx, newWellPI);
        };

        // Minimal well setup to compute PI/II values
        {
            auto saved_previous_wgstate = this->prevWGState();
            this->commitWGState();

            this->well_container_ = this->createWellContainer(timeStepIdx);
            this->inferLocalShutWells();

            for (auto& wellPtr : this->well_container_) {
                wellPtr->init(&this->phase_usage_, this->depth_, this->gravity_,
                              this->local_num_cells_, this->B_avg_);
            }

            this->calculateProductivityIndexValues(local_deferredLogger);
            this->calculateProductivityIndexValuesShutWells(timeStepIdx, local_deferredLogger);

            this->commitWGState(std::move(saved_previous_wgstate));
        }

        const auto nw = this->numLocalWells();
        for (auto wellID = 0*nw; wellID < nw; ++wellID) {
            if (hasWellPIEvent(wellID)) {
                rescaleWellPI(wellID, this->wellPI(wellID));
                updateEclWell(wellID);
            }
        }

        this->last_run_wellpi_ = timeStepIdx;
    }





    template <typename TypeTag>
    bool
    BlackoilWellModel<TypeTag>::
    wasDynamicallyShutThisTimeStep(const int well_index) const
    {
        return this->closed_this_step_.find(this->wells_ecl_[well_index].name())
            != this->closed_this_step_.end();
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    updateWsolvent(const Group& group, const Schedule& schedule, const int reportStepIdx, const WellStateFullyImplicitBlackoil& wellState) {
        for (const std::string& groupName : group.groups()) {
            const Group& groupTmp = schedule.getGroup(groupName, reportStepIdx);
            updateWsolvent(groupTmp, schedule, reportStepIdx, wellState);
        }

        if (group.isProductionGroup())
            return;

        auto currentGroupControl = this->groupState().injection_control(group.name(), Phase::GAS);
        if( currentGroupControl == Group::InjectionCMode::REIN ) {
            int gasPos = phase_usage_.phase_pos[BlackoilPhases::Vapour];
            const auto& summaryState = ebosSimulator_.vanguard().summaryState();
            const auto& controls = group.injectionControls(Phase::GAS, summaryState);
            const Group& groupRein = schedule.getGroup(controls.reinj_group, reportStepIdx);
            double gasProductionRate = WellGroupHelpers::sumWellRates(groupRein, schedule, wellState, reportStepIdx, gasPos, /*isInjector*/false);
            double solventProductionRate = WellGroupHelpers::sumSolventRates(groupRein, schedule, wellState, reportStepIdx, /*isInjector*/false);

            const auto& comm = ebosSimulator_.vanguard().grid().comm();
            solventProductionRate = comm.sum(solventProductionRate);
            gasProductionRate = comm.sum(gasProductionRate);

            double wsolvent = 0.0;
            if (std::abs(gasProductionRate) > 1e-6)
                wsolvent = solventProductionRate / gasProductionRate;

            setWsolvent(group, schedule, reportStepIdx, wsolvent);
        }
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    setWsolvent(const Group& group, const Schedule& schedule, const int reportStepIdx, double wsolvent) {
        for (const std::string& groupName : group.groups()) {
            const Group& groupTmp = schedule.getGroup(groupName, reportStepIdx);
            setWsolvent(groupTmp, schedule, reportStepIdx, wsolvent);
        }

        for (const std::string& wellName : group.wells()) {
            const auto& wellTmp = schedule.getWell(wellName, reportStepIdx);
            if (wellTmp.getStatus() == Well::Status::SHUT)
                continue;

            auto well = getWell(wellName);
            well->setWsolvent(wsolvent);
        }
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    assignWellGuideRates(data::Wells& wsrpt) const
    {
        for (const auto& well : this->wells_ecl_) {
            auto xwPos = wsrpt.find(well.name());
            if (xwPos == wsrpt.end()) { // No well results.  Unexpected.
                continue;
            }

            xwPos->second.guide_rates = this->getGuideRateValues(well);
        }
    }





    template <typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    assignShutConnections(data::Wells& wsrpt) const
    {
        auto wellID = 0;

        for (const auto& well : this->wells_ecl_) {
            auto& xwel = wsrpt[well.name()]; // data::Wells is a std::map<>

            xwel.dynamicStatus = this->schedule()
                .getWell(well.name(), this->reportStepIndex()).getStatus();

            const auto wellIsOpen = xwel.dynamicStatus == Well::Status::OPEN;
            auto skip = [wellIsOpen](const Connection& conn)
            {
                return wellIsOpen && (conn.state() != Connection::State::SHUT);
            };

            if (this->wellTestState_.hasWellClosed(well.name()) &&
                !this->wasDynamicallyShutThisTimeStep(wellID))
            {
                xwel.dynamicStatus = well.getAutomaticShutIn()
                    ? Well::Status::SHUT : Well::Status::STOP;
            }

            auto& xcon = xwel.connections;
            for (const auto& conn : well.getConnections()) {
                if (skip(conn)) {
                    continue;
                }

                auto& xc = xcon.emplace_back();
                xc.index = conn.global_index();
                xc.pressure = xc.reservoir_rate = 0.0;

                xc.effective_Kh = conn.Kh();
                xc.trans_factor = conn.CF();
            }

            ++wellID;
        }
    }





    template <typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    assignGroupValues(const int                               reportStepIdx,
                      const Schedule&                         sched,
                      std::map<std::string, data::GroupData>& gvalues) const
    {
        const auto groupGuideRates =
            this->calculateAllGroupGuiderates(reportStepIdx, sched);

        for (const auto& gname : sched.groupNames(reportStepIdx)) {
            const auto& grup = sched.getGroup(gname, reportStepIdx);

            auto& gdata = gvalues[gname];
            this->assignGroupControl(grup, gdata);
            this->assignGroupGuideRates(grup, groupGuideRates, gdata);
        }
    }

    template <typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    assignNodeValues(std::map<std::string, data::NodeData>& nodevalues) const
    {
        nodevalues.clear();
        for (const auto& [node, pressure] : node_pressures_) {
            nodevalues.emplace(node, data::NodeData{pressure});
        }
    }

    template<typename TypeTag>
    std::unordered_map<std::string, data::GroupGuideRates>
    BlackoilWellModel<TypeTag>::
    calculateAllGroupGuiderates(const int reportStepIdx, const Schedule& sched) const
    {
        auto gr = std::unordered_map<std::string, data::GroupGuideRates>{};
        auto up = std::vector<std::string>{};

        // Start at well level, accumulate contributions towards root of
        // group tree (FIELD group).

        for (const auto& wname : sched.wellNames(reportStepIdx)) {
            if (! (this->wellState().hasWellRates(wname) &&
                   this->guideRate_->has(wname)))
            {
                continue;
            }

            const auto& well   = sched.getWell(wname, reportStepIdx);
            const auto& parent = well.groupName();

            if (parent == "FIELD") {
                // Well parented directly to "FIELD".  Inadvisable and
                // unexpected, but nothing to do about that here.  Just skip
                // this guide rate contribution.
                continue;
            }

            auto& grval = well.isInjector()
                ? gr[parent].injection
                : gr[parent].production;

            grval += this->getGuideRateValues(well);
            up.push_back(parent);
        }

        // Propagate accumulated guide rates up towards root of group tree.
        // Override accumulation if there is a GUIDERAT specification that
        // applies to a group.
        std::sort(up.begin(), up.end());
        auto start = 0*up.size();
        auto u     = std::unique(up.begin(), up.end());
        auto nu    = std::distance(up.begin(), u);
        while (nu > 0) {
            const auto ntot = up.size();

            for (auto gi = 0*nu; gi < nu; ++gi) {
                const auto& gname = up[start + gi];
                const auto& group = sched.getGroup(gname, reportStepIdx);

                if (this->guideRate_->has(gname)) {
                    gr[gname].production = this->getGuideRateValues(group);
                }

                if (this->guideRate_->has(gname, Phase::WATER)
                        || this->guideRate_->has(gname, Phase::GAS)) {
                    gr[gname].injection = this->getGuideRateInjectionGroupValues(group);
                }

                const auto parent = group.parent();
                if (parent == "FIELD") { continue; }

                gr[parent].injection  += gr[gname].injection;
                gr[parent].production += gr[gname].production;
                up.push_back(parent);
            }

            start = ntot;

            auto begin = up.begin() + ntot;
            std::sort(begin, up.end());
            u  = std::unique(begin, up.end());
            nu = std::distance(begin, u);
        }

        return gr;
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    assignGroupControl(const Group& group, data::GroupData& gdata) const
    {
        const auto& gname     = group.name();
        const auto  grup_type = group.getGroupType();
        auto&       cgc       = gdata.currentControl;

        cgc.currentProdConstraint =
            ::Opm::Group::ProductionCMode::NONE;

        cgc.currentGasInjectionConstraint =
        cgc.currentWaterInjectionConstraint =
            ::Opm::Group::InjectionCMode::NONE;

        if (this->groupState().has_production_control(gname)) {
            cgc.currentProdConstraint = this->groupState().production_control(gname);
        }

        if ((grup_type == ::Opm::Group::GroupType::INJECTION) ||
            (grup_type == ::Opm::Group::GroupType::MIXED))
        {
            if (this->groupState().has_injection_control(gname, ::Opm::Phase::WATER)) {
                cgc.currentWaterInjectionConstraint = this->groupState().injection_control(gname, Phase::WATER);
            }

            if (this->groupState().has_injection_control(gname, ::Opm::Phase::GAS)) {
                cgc.currentGasInjectionConstraint = this->groupState().injection_control(gname, Phase::GAS);
            }
        }
    }

    template <typename TypeTag>
    data::GuideRateValue
    BlackoilWellModel<TypeTag>::
    getGuideRateValues(const Well& well) const
    {
        auto grval = data::GuideRateValue{};

        assert (this->guideRate_ != nullptr);

        const auto& wname = well.name();
        if (! this->wellState().hasWellRates(wname)) {
            // No flow rates for 'wname' -- might be before well comes
            // online (e.g., for the initial condition before simulation
            // starts).
            return grval;
        }

        if (! this->guideRate_->has(wname)) {
            // No guiderates exist for 'wname'.
            return grval;
        }

        const auto qs = WellGroupHelpers::
            getWellRateVector(this->wellState(), this->phase_usage_, wname);

        this->getGuideRateValues(qs, well.isInjector(), wname, grval);

        return grval;
    }
    template <typename TypeTag>
    data::GuideRateValue
    BlackoilWellModel<TypeTag>::
    getGuideRateInjectionGroupValues(const Group& group) const
    {
        auto grval = data::GuideRateValue{};

        assert (this->guideRate_ != nullptr);

        const auto& gname = group.name();
        if (this->guideRate_->has(gname, Phase::GAS)) {
            grval.set(data::GuideRateValue::Item::Gas,
                      this->guideRate_->get(gname, Phase::GAS));
        }
        if (this->guideRate_->has(gname, Phase::WATER)) {
            grval.set(data::GuideRateValue::Item::Water,
                      this->guideRate_->get(gname, Phase::WATER));
        }
        return grval;
    }

    template <typename TypeTag>
    data::GuideRateValue
    BlackoilWellModel<TypeTag>::
    getGuideRateValues(const Group& group) const
    {
        auto grval = data::GuideRateValue{};

        assert (this->guideRate_ != nullptr);

        const auto& gname = group.name();

        if ( ! this->groupState().has_production_rates(gname)) {
            // No flow rates for production group 'gname' -- might be before group comes
            // online (e.g., for the initial condition before simulation
            // starts).
            return grval;
        }

        if (! this->guideRate_->has(gname)) {
            // No guiderates exist for 'gname'.
            return grval;
        }

        const auto qs = WellGroupHelpers::getProductionGroupRateVector(this->groupState(), this->phase_usage_, gname);

        const auto is_inj = false; // This procedure only applies to G*PGR.
        this->getGuideRateValues(qs, is_inj, gname, grval);

        return grval;
    }

    template <typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    getGuideRateValues(const GuideRate::RateVector& qs,
                       const bool                   is_inj,
                       const std::string&           wgname,
                       data::GuideRateValue&        grval) const
    {
        auto getGR = [this, &wgname, &qs](const GuideRateModel::Target t)
        {
            return this->guideRate_->get(wgname, t, qs);
        };

        // Note: GuideRate does currently (2020-07-20) not support Target::RES.
        grval.set(data::GuideRateValue::Item::Gas,
                  getGR(GuideRateModel::Target::GAS));

        grval.set(data::GuideRateValue::Item::Water,
                  getGR(GuideRateModel::Target::WAT));

        if (! is_inj) {
            // Producer.  Extract "all" guiderate values.
            grval.set(data::GuideRateValue::Item::Oil,
                      getGR(GuideRateModel::Target::OIL));
        }
    }

    template <typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    assignGroupGuideRates(const Group& group,
                          const std::unordered_map<std::string, data::GroupGuideRates>& groupGuideRates,
                          data::GroupData& gdata) const
    {
        auto& prod = gdata.guideRates.production;  prod.clear();
        auto& inj  = gdata.guideRates.injection;   inj .clear();

        auto xgrPos = groupGuideRates.find(group.name());
        if ((xgrPos == groupGuideRates.end()) ||
            !this->guideRate_->has(group.name()))
        {
            // No guiderates defined for this group.
            return;
        }

        const auto& xgr = xgrPos->second;
        prod = xgr.production;
        inj  = xgr.injection;
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
            double weighted_temperature = 0.0;
            double total_weight = 0.0;

            auto& well_info = *local_parallel_well_info_[wellID];
            const int num_perf_this_well = well_info.communication().sum(well_perf_data_[wellID].size());

            for (int perf = 0; perf < num_perf_this_well; ++perf) {
                const int cell_idx = well_perf_data_[wellID][perf].cell_index;
                const auto& intQuants = *(ebosSimulator_.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
                const auto& fs = intQuants.fluidState();

                double cellTemperatures = fs.temperature(/*phaseIdx*/0).value();
                double weight_factor = 0.0;
                for (int phaseIdx = 0; phaseIdx  < np; ++phaseIdx) {
                    cellInternalEnergy = fs.enthalpy(phaseIdx).value() - fs.pressure(phaseIdx).value() / fs.density(phaseIdx).value();
                    cellBinv = fs.invB(phaseIdx).value();
                    cellDensity = fs.density(phaseIdx).value();
                    perfPhaseRate = this->wellState().perfPhaseRates()[connpos + perf*np + phaseIdx ];
                    weight_factor += cellDensity  * perfPhaseRate/cellBinv * cellInternalEnergy/cellTemperatures;
                }
                total_weight += weight_factor;
                weighted_temperature += weight_factor * cellTemperatures;
            }
            weighted_temperature = well_info.communication().sum(weighted_temperature);
            total_weight = well_info.communication().sum(total_weight);
            this->wellState().temperature()[wellID] =  weighted_temperature/total_weight;
        }
    }
} // namespace Opm
