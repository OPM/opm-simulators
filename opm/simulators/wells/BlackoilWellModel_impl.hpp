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
#include <opm/simulators/wells/SimFIBODetails.hpp>
#include <opm/core/props/phaseUsageFromDeck.hpp>

namespace Opm {
    template<typename TypeTag>
    BlackoilWellModel<TypeTag>::
    BlackoilWellModel(Simulator& ebosSimulator)
        : ebosSimulator_(ebosSimulator)
        , has_solvent_(GET_PROP_VALUE(TypeTag, EnableSolvent))
        , has_polymer_(GET_PROP_VALUE(TypeTag, EnablePolymer))
    {
        terminal_output_ = false;
        if (ebosSimulator.gridView().comm().rank() == 0)
            terminal_output_ = EWOMS_GET_PARAM(TypeTag, bool, EnableTerminalOutput);

        // Create the guide rate container.
        guideRate_.reset(new GuideRate (ebosSimulator_.vanguard().schedule()));

        // calculate the number of elements of the compressed sequential grid. this needs
        // to be done in two steps because the dune communicator expects a reference as
        // argument for sum()
        const auto& gridView = ebosSimulator_.gridView();
        number_of_cells_ = gridView.size(/*codim=*/0);
        global_nc_ = gridView.comm().sum(number_of_cells_);

        // Set up cartesian mapping.
        const auto& grid = ebosSimulator_.vanguard().grid();
        const auto& cartDims = Opm::UgGridHelpers::cartDims(grid);
        setupCartesianToCompressed_(Opm::UgGridHelpers::globalCell(grid),
                                    cartDims[0]*cartDims[1]*cartDims[2]);
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    init()
    {
        const Opm::EclipseState& eclState = ebosSimulator_.vanguard().eclState();

        extractLegacyCellPvtRegionIndex_();
        extractLegacyDepth_();

        phase_usage_ = phaseUsageFromDeck(eclState);

        gravity_ = ebosSimulator_.problem().gravity()[2];

        initial_step_ = true;

        // add the eWoms auxiliary module for the wells to the list
        ebosSimulator_.model().addAuxiliaryModule(this);

        is_cell_perforated_.resize(number_of_cells_, false);
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
        const auto& cartesianSize = Opm::UgGridHelpers::cartDims(grid());

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
                int i = connection.getI();
                int j = connection.getJ();
                int k = connection.getK();
                int cart_grid_idx = i + cartesianSize[0]*(j + cartesianSize[1]*k);
                int compressed_idx = cartesian_to_compressed_.at(cart_grid_idx);

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
    void
    BlackoilWellModel<TypeTag>::
    beginReportStep(const int timeStepIdx)
    {
        Opm::DeferredLogger local_deferredLogger;
        report_step_starts_ = true;

        const Grid& grid = ebosSimulator_.vanguard().grid();
        const auto& summaryState = ebosSimulator_.vanguard().summaryState();
        int globalNumWells = 0;
        // Make wells_ecl_ contain only this partition's non-shut wells.
        {
            const auto& defunct_well_names = ebosSimulator_.vanguard().defunctWellNames();
            auto is_shut_or_defunct = [&defunct_well_names](const Well& well) {
                return (well.getStatus() == Well::Status::SHUT) || (defunct_well_names.find(well.name()) != defunct_well_names.end());
            };
            auto w = schedule().getWells(timeStepIdx);
            globalNumWells = w.size();
            w.erase(std::remove_if(w.begin(), w.end(), is_shut_or_defunct), w.end());
            wells_ecl_.swap(w);
        }
        initializeWellPerfData();

        // Wells are active if they are active wells on at least
        // one process.
        wells_active_ = localWellsActive() ? 1 : 0;
        wells_active_ = grid.comm().max(wells_active_);

        // The well state initialize bhp with the cell pressure in the top cell.
        // We must therefore provide it with updated cell pressures
        size_t nc = number_of_cells_;
        std::vector<double> cellPressures(nc, 0.0);
        ElementContext elemCtx(ebosSimulator_);
        const auto& gridView = ebosSimulator_.vanguard().gridView();
        const auto& elemEndIt = gridView.template end</*codim=*/0>();
        for (auto elemIt = gridView.template begin</*codim=*/0>();
             elemIt != elemEndIt;
             ++elemIt)
        {
            const auto& elem = *elemIt;
            if (elem.partitionType() != Dune::InteriorEntity) {
                continue;
            }
            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

            const unsigned cellIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& intQuants = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& fs = intQuants.fluidState();
            // copy of get perfpressure in Standard well
            // exept for value
            double perf_pressure = 0.0;
            if (Indices::oilEnabled) {
                perf_pressure = fs.pressure(FluidSystem::oilPhaseIdx).value();
            } else {
                if (Indices::waterEnabled) {
                    perf_pressure = fs.pressure(FluidSystem::waterPhaseIdx).value();
                } else {
                    perf_pressure = fs.pressure(FluidSystem::gasPhaseIdx).value();
                }
            }
            cellPressures[cellIdx] = perf_pressure;
        }
        well_state_.init(cellPressures, schedule(), wells_ecl_, timeStepIdx, &previous_well_state_, phase_usage_, well_perf_data_, summaryState, globalNumWells);

        // handling MS well related
        if (param_.use_multisegment_well_&& anyMSWellOpenLocal()) { // if we use MultisegmentWell model
            well_state_.initWellStateMSWell(wells_ecl_, phase_usage_, &previous_well_state_);
        }

        const int nw = wells_ecl_.size();
        for (int w = 0; w <nw; ++w) {
            const auto& well = wells_ecl_[w];
            const uint64_t effective_events_mask = ScheduleEvents::WELL_STATUS_CHANGE
                                                 + ScheduleEvents::PRODUCTION_UPDATE
                                                 + ScheduleEvents::INJECTION_UPDATE
                                                 + ScheduleEvents::NEW_WELL;

            if(!schedule().hasWellGroupEvent(well.name(), effective_events_mask, timeStepIdx))
                continue;

            if (well.isProducer()) {
                const auto controls = well.productionControls(summaryState);
                well_state_.currentProductionControls()[w] = controls.cmode;
            }
            else {
                const auto controls = well.injectionControls(summaryState);
                well_state_.currentInjectionControls()[w] = controls.cmode;
            }
        }
        const Group& fieldGroup = schedule().getGroup("FIELD", timeStepIdx);
        WellGroupHelpers::setCmodeGroup(fieldGroup, schedule(), summaryState, timeStepIdx, well_state_);

        // Compute reservoir volumes for RESV controls.
        rateConverter_.reset(new RateConverterType (phase_usage_,
                                                    std::vector<int>(number_of_cells_, 0)));
        rateConverter_->template defineState<ElementContext>(ebosSimulator_);

        // update VFP properties
        vfp_properties_.reset (new VFPProperties<VFPInjProperties,VFPProdProperties> (
                                   schedule().getVFPInjTables(timeStepIdx),
                                   schedule().getVFPProdTables(timeStepIdx)) );

        // update the previous well state. This is used to restart failed steps.
        previous_well_state_ = well_state_;

    }


    // called at the beginning of a time step
    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    beginTimeStep() {

        Opm::DeferredLogger local_deferredLogger;

        well_state_ = previous_well_state_;

        const int reportStepIdx = ebosSimulator_.episodeIndex();
        const double simulationTime = ebosSimulator_.time();

        int exception_thrown = 0;
        try {
            // test wells
            wellTesting(reportStepIdx, simulationTime, local_deferredLogger);

            // create the well container
            well_container_ = createWellContainer(reportStepIdx);

            // do the initialization for all the wells
            // TODO: to see whether we can postpone of the intialization of the well containers to
            // optimize the usage of the following several member variables
            for (auto& well : well_container_) {
                well->init(&phase_usage_, depth_, gravity_, number_of_cells_);
            }

            // update the updated cell flag
            std::fill(is_cell_perforated_.begin(), is_cell_perforated_.end(), false);
            for (auto& well : well_container_) {
                well->updatePerforatedCell(is_cell_perforated_);
            }

            // calculate the efficiency factors for each well
            calculateEfficiencyFactors(reportStepIdx);

            if (has_polymer_)
            {
                const Grid& grid = ebosSimulator_.vanguard().grid();
                if (PolymerModule::hasPlyshlog() || GET_PROP_VALUE(TypeTag, EnablePolymerMW) ) {
                        computeRepRadiusPerfLength(grid, local_deferredLogger);
                }
            }
        } catch (std::exception& e) {
            exception_thrown = 1;
        }

        logAndCheckForExceptionsAndThrow(local_deferredLogger, exception_thrown, "beginTimeStep() failed.", terminal_output_);

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

        //compute well guideRates
        const auto& comm = ebosSimulator_.vanguard().grid().comm();
        WellGroupHelpers::updateGuideRatesForWells(schedule(), phase_usage_, reportStepIdx, simulationTime, well_state_, comm, guideRate_.get());
    }


    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::wellTesting(const int timeStepIdx, const double simulationTime, Opm::DeferredLogger& deferred_logger) {
        const auto& wtest_config = schedule().wtestConfig(timeStepIdx);
        if (wtest_config.size() != 0) { // there is a WTEST request

            // average B factors are required for the convergence checking of well equations
            // Note: this must be done on all processes, even those with
            // no wells needing testing, otherwise we will have locking.
            std::vector< Scalar > B_avg(numComponents(), Scalar() );
            computeAverageFormationFactor(B_avg);

            const auto& wellsForTesting = wellTestState_.updateWells(wtest_config, wells_ecl_, simulationTime);
            for (const auto& testWell : wellsForTesting) {
                const std::string& well_name = testWell.first;

                // this is the well we will test
                WellInterfacePtr well = createWellForWellTest(well_name, timeStepIdx, deferred_logger);

                // some preparation before the well can be used
                well->init(&phase_usage_, depth_, gravity_, number_of_cells_);
                const Well& wellEcl = schedule().getWell(well_name, timeStepIdx);
                double well_efficiency_factor = wellEcl.getEfficiencyFactor();
                WellGroupHelpers::accumulateGroupEfficiencyFactor(schedule().getGroup(wellEcl.groupName(), timeStepIdx), schedule(), timeStepIdx, well_efficiency_factor);
                well->setWellEfficiencyFactor(well_efficiency_factor);
                well->setVFPProperties(vfp_properties_.get());
                well->setGuideRate(guideRate_.get());

                const WellTestConfig::Reason testing_reason = testWell.second;

                well->wellTesting(ebosSimulator_, B_avg, simulationTime, timeStepIdx,
                                  testing_reason, well_state_, wellTestState_, deferred_logger);
            }
        }
    }

    // called at the end of a report step
    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    endReportStep() {
    }

    // called at the end of a report step
    template<typename TypeTag>
    const SimulatorReport&
    BlackoilWellModel<TypeTag>::
    lastReport() const {return last_report_; }

    // called at the end of a time step
    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    timeStepSucceeded(const double& simulationTime, const double dt) {

        // time step is finished and we are not any more at the beginning of an report step
        report_step_starts_ = false;

        Opm::DeferredLogger local_deferredLogger;
        for (const auto& well : well_container_) {
            if (GET_PROP_VALUE(TypeTag, EnablePolymerMW) && well->isInjector()) {
                well->updateWaterThroughput(dt, well_state_);
            }
        }
        updateWellTestState(simulationTime, wellTestState_);

        // update the rate converter with current averages pressures etc in
        rateConverter_->template defineState<ElementContext>(ebosSimulator_);

        // calculate the well potentials
        try {
            std::vector<double> well_potentials;
            const int reportStepIdx = ebosSimulator_.episodeIndex();
            computeWellPotentials(well_potentials, reportStepIdx, local_deferredLogger);
        } catch ( std::runtime_error& e ) {
            const std::string msg = "A zero well potential is returned for output purposes. ";
            local_deferredLogger.warning("WELL_POTENTIAL_CALCULATION_FAILED", msg);
        }
        previous_well_state_ = well_state_;

        Opm::DeferredLogger global_deferredLogger = gatherDeferredLogger(local_deferredLogger);
        if (terminal_output_) {
            global_deferredLogger.logMessages();
        }

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
        int globalNumWells = 0;
        // Make wells_ecl_ contain only this partition's non-shut wells.
        {
            const auto& defunct_well_names = ebosSimulator_.vanguard().defunctWellNames();
            auto is_shut_or_defunct = [&defunct_well_names](const Well& well) {
                return (well.getStatus() == Well::Status::SHUT) || (defunct_well_names.find(well.name()) != defunct_well_names.end());
            };
            auto w = schedule().getWells(report_step);
            globalNumWells = w.size();
            w.erase(std::remove_if(w.begin(), w.end(), is_shut_or_defunct), w.end());
            wells_ecl_.swap(w);
        }

        initializeWellPerfData();

        const int nw = wells_ecl_.size();
        if (nw > 0) {
            const auto phaseUsage = phaseUsageFromDeck(eclState());
            const size_t numCells = Opm::UgGridHelpers::numCells(grid());
            const bool handle_ms_well = (param_.use_multisegment_well_ && anyMSWellOpenLocal());
            well_state_.resize(wells_ecl_, schedule(), handle_ms_well, numCells, phaseUsage, well_perf_data_, summaryState, globalNumWells); // Resize for restart step
            wellsToState(restartValues.wells, phaseUsage, handle_ms_well, well_state_);
        }

        previous_well_state_ = well_state_;

        initial_step_ = false;
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    initializeWellPerfData()
    {
        const auto& grid = ebosSimulator_.vanguard().grid();
        const auto& cartDims = Opm::UgGridHelpers::cartDims(grid);
        well_perf_data_.resize(wells_ecl_.size());
        first_perf_index_.clear();
        first_perf_index_.resize(wells_ecl_.size() + 1, 0);
        int well_index = 0;
        for (const auto& well : wells_ecl_) {
            well_perf_data_[well_index].clear();
            well_perf_data_[well_index].reserve(well.getConnections().size());
            for (const auto& completion : well.getConnections()) {
                if (completion.state() == Connection::State::OPEN) {
                    const int i = completion.getI();
                    const int j = completion.getJ();
                    const int k = completion.getK();
                    const int cart_grid_indx = i + cartDims[0] * (j + cartDims[1] * k);
                    const int active_index = cartesian_to_compressed_[cart_grid_indx];
                    if (active_index < 0) {
                        const std::string msg
                            = ("Cell with i,j,k indices " + std::to_string(i) + " " + std::to_string(j) + " "
                               + std::to_string(k) + " not found in grid (well = " + well.name() + ").");
                        OPM_THROW(std::runtime_error, msg);
                    } else {
                        PerforationData pd;
                        pd.cell_index = active_index;
                        pd.connection_transmissibility_factor = completion.CF();
                        pd.satnum_id = completion.satTableId();
                        well_perf_data_[well_index].push_back(pd);
                    }
                } else {
                    if (completion.state() != Connection::State::SHUT) {
                        OPM_THROW(std::runtime_error,
                                  "Completion state: " << Connection::State2String(completion.state()) << " not handled");
                    }
                }
            }
            first_perf_index_[well_index + 1] = first_perf_index_[well_index] + well_perf_data_[well_index].size();
            ++well_index;
        }
    }





    template<typename TypeTag>
    std::vector<typename BlackoilWellModel<TypeTag>::WellInterfacePtr >
    BlackoilWellModel<TypeTag>::
    createWellContainer(const int time_step)
    {
        std::vector<WellInterfacePtr> well_container;

        Opm::DeferredLogger local_deferredLogger;

        const int nw = numLocalWells();

        if (nw > 0) {
            well_container.reserve(nw);
            for (int w = 0; w < nw; ++w) {
                const Well& well_ecl = wells_ecl_[w];
                const std::string& well_name = well_ecl.name();

                // A new WCON keywords can re-open a well that was closed/shut due to Physical limit
                if ( wellTestState_.hasWellClosed(well_name)) {
                    // TODO: more checking here, to make sure this standard more specific and complete
                    // maybe there is some WCON keywords will not open the well
                    if (well_state_.effectiveEventsOccurred(w)) {
                        if (wellTestState_.lastTestTime(well_name) == ebosSimulator_.time()) {
                            // The well was shut this timestep, we are most likely retrying
                            // a timestep without the well in question, after it caused
                            // repeated timestep cuts. It should therefore not be opened,
                            // even if it was new or received new targets this report step.
                            well_state_.setEffectiveEventsOccurred(w, false);
                        } else {
                            wellTestState_.openWell(well_name);
                        }
                    }
                }

                // TODO: should we do this for all kinds of closing reasons?
                // something like wellTestState_.hasWell(well_name)?
                bool wellIsStopped = false;
                if ( wellTestState_.hasWellClosed(well_name, WellTestConfig::Reason::ECONOMIC) ||
                    wellTestState_.hasWellClosed(well_name, WellTestConfig::Reason::PHYSICAL) ) {
                    if( well_ecl.getAutomaticShutIn() ) {
                        // shut wells are not added to the well container
                        well_state_.shutWell(w);
                        continue;
                    } else {
                        // stopped wells are added to the container but marked as stopped
                        well_state_.thp()[w] = 0.;
                        wellIsStopped = true;
                    }
                }

                // Due to ACTIONX the well might have been closed 'behind our back'.
                const auto well_status = schedule().getWell(well_name, time_step).getStatus();
                if (well_status == Well::Status::SHUT) {
                    well_state_.shutWell(w);
                    continue;
                }

                // If a production well disallows crossflow and its
                // (prediction type) rate control is zero, then it is effectively shut.
                if (!well_ecl.getAllowCrossFlow() && well_ecl.isProducer() && well_ecl.predictionMode()) {
                    const auto& summaryState = ebosSimulator_.vanguard().summaryState();
                    auto prod_controls = well_ecl.productionControls(summaryState);
                    bool zero_rate_control = false;
                    switch (prod_controls.cmode) {
                    case Well::ProducerCMode::ORAT:
                        zero_rate_control = (prod_controls.oil_rate == 0.0);
                        break;
                    case Well::ProducerCMode::WRAT:
                        zero_rate_control = (prod_controls.water_rate == 0.0);
                        break;
                    case Well::ProducerCMode::GRAT:
                        zero_rate_control = (prod_controls.gas_rate == 0.0);
                        break;
                    case Well::ProducerCMode::LRAT:
                        zero_rate_control = (prod_controls.liquid_rate == 0.0);
                        break;
                    case Well::ProducerCMode::RESV:
                        zero_rate_control = (prod_controls.resv_rate == 0.0);
                        break;
                    default:
                        // Might still have zero rate controls, but is pressure controlled.
                        zero_rate_control = false;
                    }
                    if (zero_rate_control) {
                        // Treat as shut, do not add to container.
                        local_deferredLogger.info("  Well shut due to zero rate control and disallowing crossflow: " + well_ecl.name());
                        well_state_.shutWell(w);
                        continue;
                    }
                }

                if (well_status == Well::Status::STOP) {
                    well_state_.thp()[w] = 0.;
                    wellIsStopped = true;
                }

                // Use the pvtRegionIdx from the top cell
                const int well_cell_top = well_perf_data_[w][0].cell_index;
                const int pvtreg = pvt_region_idx_[well_cell_top];

                if (!well_ecl.isMultiSegment() || !param_.use_multisegment_well_) {
                    well_container.emplace_back(new StandardWell<TypeTag>(well_ecl,
                                                                          time_step,
                                                                          param_,
                                                                          *rateConverter_,
                                                                          pvtreg,
                                                                          numComponents(),
                                                                          numPhases(),
                                                                          w,
                                                                          first_perf_index_[w],
                                                                          well_perf_data_[w]));
                } else {
                    well_container.emplace_back(new MultisegmentWell<TypeTag>(well_ecl,
                                                                              time_step,
                                                                              param_,
                                                                              *rateConverter_,
                                                                              pvtreg,
                                                                              numComponents(),
                                                                              numPhases(),
                                                                              w,
                                                                              first_perf_index_[w],
                                                                              well_perf_data_[w]));
                }
                if (wellIsStopped)
                    well_container.back()->stopWell();
            }
        }

        // Collect log messages and print.
        Opm::DeferredLogger global_deferredLogger = gatherDeferredLogger(local_deferredLogger);
        if (terminal_output_) {
            global_deferredLogger.logMessages();
        }

        return well_container;
    }





    template<typename TypeTag>
    typename BlackoilWellModel<TypeTag>::WellInterfacePtr
    BlackoilWellModel<TypeTag>::
    createWellForWellTest(const std::string& well_name,
                          const int report_step,
                          Opm::DeferredLogger& deferred_logger) const
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

        const Well& well_ecl = wells_ecl_[index_well_ecl];

        // Use the pvtRegionIdx from the top cell
        const int well_cell_top = well_perf_data_[index_well_ecl][0].cell_index;
        const int pvtreg = pvt_region_idx_[well_cell_top];

        if (!well_ecl.isMultiSegment() || !param_.use_multisegment_well_) {
            return WellInterfacePtr(new StandardWell<TypeTag>(well_ecl,
                                                              report_step,
                                                              param_,
                                                              *rateConverter_,
                                                              pvtreg,
                                                              numComponents(),
                                                              numPhases(),
                                                              index_well_ecl,
                                                              first_perf_index_[index_well_ecl],
                                                              well_perf_data_[index_well_ecl]));
        } else {
            return WellInterfacePtr(new MultisegmentWell<TypeTag>(well_ecl,
                                                                  report_step,
                                                                  param_,
                                                                  *rateConverter_,
                                                                  pvtreg,
                                                                  numComponents(),
                                                                  numPhases(),
                                                                  index_well_ecl,
                                                                  first_perf_index_[index_well_ecl],
                                                                  well_perf_data_[index_well_ecl]));
        }
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    assemble(const int iterationIdx,
             const double dt)
    {

        last_report_ = SimulatorReport();

        if ( ! wellsActive() ) {
            return;
        }

        Opm::DeferredLogger local_deferredLogger;

        updatePerforationIntensiveQuantities();

        int exception_thrown = 0;
        try {
            if (iterationIdx == 0) {
                calculateExplicitQuantities(local_deferredLogger);
                prepareTimeStep(local_deferredLogger);
            }
            updateWellControls(local_deferredLogger, /* check group controls */ true);

            // Set the well primary variables based on the value of well solutions
            initPrimaryVariablesEvaluation();

            std::vector< Scalar > B_avg(numComponents(), Scalar() );
            computeAverageFormationFactor(B_avg);

            if (param_.solve_welleq_initially_ && iterationIdx == 0) {
                // solve the well equations as a pre-processing step
                last_report_ = solveWellEq(B_avg, dt, local_deferredLogger);


                if (initial_step_) {
                    // update the explicit quantities to get the initial fluid distribution in the well correct.
                    calculateExplicitQuantities(local_deferredLogger);
                    prepareTimeStep(local_deferredLogger);
                    last_report_ = solveWellEq(B_avg, dt, local_deferredLogger);
                    initial_step_ = false;
                }
                // TODO: should we update the explicit related here again, or even prepareTimeStep().
                // basically, this is a more updated state from the solveWellEq based on fixed
                // reservoir state, will tihs be a better place to inialize the explict information?
            }

            assembleWellEq(B_avg, dt, local_deferredLogger);

        } catch (std::exception& e) {
            exception_thrown = 1;
        }
        logAndCheckForExceptionsAndThrow(local_deferredLogger, exception_thrown, "assemble() failed.", terminal_output_);

        last_report_.converged = true;
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    assembleWellEq(const std::vector<Scalar>& B_avg, const double dt, Opm::DeferredLogger& deferred_logger)
    {
        for (auto& well : well_container_) {
            well->assembleWellEq(ebosSimulator_, B_avg, dt, well_state_, deferred_logger);
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

#if HAVE_CUDA
    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    getWellContributions(WellContributions& wellContribs) const
    {
        wellContribs.setBlockSize(StandardWell<TypeTag>::numEq, StandardWell<TypeTag>::numStaticWellEq);
        for(unsigned int i = 0; i < well_container_.size(); i++){
            auto& well = well_container_[i];
            std::shared_ptr<StandardWell<TypeTag> > derived = std::dynamic_pointer_cast<StandardWell<TypeTag> >(well);
            unsigned int numBlocks;
            derived->getNumBlocks(numBlocks);
            wellContribs.addNumBlocks(numBlocks);
        }
        wellContribs.alloc();
        for(unsigned int i = 0; i < well_container_.size(); i++){
            auto& well = well_container_[i];
            std::shared_ptr<StandardWell<TypeTag> > derived = std::dynamic_pointer_cast<StandardWell<TypeTag> >(well);
            if (derived) {
                derived->addWellContribution(wellContribs);
            } else {
                OpmLog::warning("Warning only StandardWell is supported by WellContributions for GPU");
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
        Opm::DeferredLogger local_deferredLogger;

        int exception_thrown = 0;
        try {
            if (localWellsActive()) {
                for (auto& well : well_container_) {
                    well->recoverWellSolutionAndUpdateWellState(x, well_state_, local_deferredLogger);
                }
            }
        } catch (std::exception& e) {
            exception_thrown = 1;
        }
        logAndCheckForExceptionsAndThrow(local_deferredLogger, exception_thrown, "recoverWellSolutionAndUpdateWellState() failed.", terminal_output_);
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
    SimulatorReport
    BlackoilWellModel<TypeTag>::
    solveWellEq(const std::vector<Scalar>& B_avg, const double dt, Opm::DeferredLogger& deferred_logger)
    {
        WellState well_state0 = well_state_;

        const int max_iter = param_.max_welleq_iter_;

        int it  = 0;
        bool converged;
        int exception_thrown = 0;
        do {
            try {
                assembleWellEq(B_avg, dt, deferred_logger);
            } catch (std::exception& e) {
                exception_thrown = 1;
            }
            // We need to check on all processes, as getWellConvergence() below communicates on all processes.
            logAndCheckForExceptionsAndThrow(deferred_logger, exception_thrown, "solveWellEq() failed.", terminal_output_);

            const auto report = getWellConvergence(B_avg);
            converged = report.converged();

            if (converged) {
                break;
            }

            try {
                if( localWellsActive() )
                {
                    for (auto& well : well_container_) {
                        well->solveEqAndUpdateWellState(well_state_, deferred_logger);
                    }
                }
                // updateWellControls uses communication
                // Therefore the following is executed if there
                // are active wells anywhere in the global domain.
                if( wellsActive() )
                {
                    updateWellControls(deferred_logger, /*don't switch group controls*/false);
                    initPrimaryVariablesEvaluation();
                }
            } catch (std::exception& e) {
                exception_thrown = 1;
            }

            logAndCheckForExceptionsAndThrow(deferred_logger, exception_thrown, "solveWellEq() failed.", terminal_output_);
            ++it;
        } while (it < max_iter);

        try {
            if (converged) {
                if (terminal_output_) {
                    deferred_logger.debug("Well equation solution gets converged with " + std::to_string(it) + " iterations");
                }
            } else {
                if (terminal_output_) {
                    deferred_logger.debug("Well equation solution failed in getting converged with " + std::to_string(it) + " iterations");
                }
                well_state_ = well_state0;
                updatePrimaryVariables(deferred_logger);
            }
        } catch (std::exception& e) {
            exception_thrown = 1;
        }

        logAndCheckForExceptionsAndThrow(deferred_logger, exception_thrown, "solveWellEq() failed.", terminal_output_);

        SimulatorReport report;
        report.converged = converged;
        report.total_well_iterations = it;
        return report;
    }





    template<typename TypeTag>
    ConvergenceReport
    BlackoilWellModel<TypeTag>::
    getWellConvergence(const std::vector<Scalar>& B_avg, bool checkGroupConvergence) const
    {

        Opm::DeferredLogger local_deferredLogger;
        // Get global (from all processes) convergence report.
        ConvergenceReport local_report;
        for (const auto& well : well_container_) {
            if (well->isOperable() ) {
                local_report += well->getWellConvergence(well_state_, B_avg, local_deferredLogger);
            }
        }
        Opm::DeferredLogger global_deferredLogger = gatherDeferredLogger(local_deferredLogger);
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
    calculateExplicitQuantities(Opm::DeferredLogger& deferred_logger) const
    {
        // TODO: checking isOperable() ?
        for (auto& well : well_container_) {
            well->calculateExplicitQuantities(ebosSimulator_, well_state_, deferred_logger);
        }
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    updateWellControls(Opm::DeferredLogger& deferred_logger, const bool checkGroupControls)
    {
        // Even if there are no wells active locally, we cannot
        // return as the DeferredLogger uses global communication.
        // For no well active globally we simply return.
        if( !wellsActive() ) return ;

        updateAndCommunicateGroupData();

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
                const bool changed = well->updateWellControl(ebosSimulator_, mode, well_state_, deferred_logger);
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
            well->updateWellControl(ebosSimulator_, mode, well_state_, deferred_logger);
        }
        updateAndCommunicateGroupData();

    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    updateAndCommunicateGroupData()
    {
        const int reportStepIdx = ebosSimulator_.episodeIndex();
        const Group& fieldGroup = schedule().getGroup("FIELD", reportStepIdx);
        const int nupcol = schedule().getNupcol(reportStepIdx);
        const int iterationIdx = ebosSimulator_.model().newtonMethod().numIterations();

        // This builds some necessary lookup structures, so it must be called
        // before we copy to well_state_nupcol_.
        const auto& comm = ebosSimulator_.vanguard().grid().comm();
        well_state_.updateGlobalIsGrup(schedule(), reportStepIdx, comm);

        if (iterationIdx < nupcol) {
            well_state_nupcol_ = well_state_;
        }

        // the group target reduction rates needs to be update since wells may have swicthed to/from GRUP control
        // Currently the group target reduction does not honor NUPCOL. TODO: is that true?
        std::vector<double> groupTargetReduction(numPhases(), 0.0);
        WellGroupHelpers::updateGroupTargetReduction(fieldGroup, schedule(), reportStepIdx, /*isInjector*/ false, phase_usage_, *guideRate_, well_state_nupcol_, well_state_, groupTargetReduction);
        std::vector<double> groupTargetReductionInj(numPhases(), 0.0);
        WellGroupHelpers::updateGroupTargetReduction(fieldGroup, schedule(), reportStepIdx, /*isInjector*/ true, phase_usage_, *guideRate_, well_state_nupcol_, well_state_, groupTargetReductionInj);

        const double simulationTime = ebosSimulator_.time();
        std::vector<double> pot(numPhases(), 0.0);
        WellGroupHelpers::updateGuideRateForGroups(fieldGroup, schedule(), phase_usage_, reportStepIdx, simulationTime, /*isInjector*/ false, well_state_, comm, guideRate_.get(), pot);
        std::vector<double> potInj(numPhases(), 0.0);
        WellGroupHelpers::updateGuideRateForGroups(fieldGroup, schedule(), phase_usage_, reportStepIdx, simulationTime, /*isInjector*/ true, well_state_, comm, guideRate_.get(), potInj);

        const auto& summaryState = ebosSimulator_.vanguard().summaryState();
        WellGroupHelpers::updateREINForGroups(fieldGroup, schedule(), reportStepIdx, phase_usage_, summaryState, well_state_nupcol_, well_state_);
        WellGroupHelpers::updateVREPForGroups(fieldGroup, schedule(), reportStepIdx, well_state_nupcol_, well_state_);

        WellGroupHelpers::updateReservoirRatesInjectionGroups(fieldGroup, schedule(), reportStepIdx, well_state_nupcol_, well_state_);
        WellGroupHelpers::updateGroupProductionRates(fieldGroup, schedule(), reportStepIdx, well_state_nupcol_, well_state_);
        WellGroupHelpers::updateWellRates(fieldGroup, schedule(), reportStepIdx, well_state_nupcol_, well_state_);
        well_state_.communicateGroupRates(comm);

        // compute wsolvent fraction for REIN wells
        updateWsolvent(fieldGroup, schedule(), reportStepIdx,  well_state_nupcol_);

    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    updateWellTestState(const double& simulationTime, WellTestState& wellTestState) const
    {
        Opm::DeferredLogger local_deferredLogger;
        for (const auto& well : well_container_) {
            well->updateWellTestState(well_state_, simulationTime, /*writeMessageToOPMLog=*/ true, wellTestState, local_deferredLogger);
        }
        Opm::DeferredLogger global_deferredLogger = gatherDeferredLogger(local_deferredLogger);
        if (terminal_output_) {
            global_deferredLogger.logMessages();
        }
    }



    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    computeWellPotentials(std::vector<double>& well_potentials, const int reportStepIdx, Opm::DeferredLogger& deferred_logger)
    {
        // number of wells and phases
        const int nw = numLocalWells();
        const int np = numPhases();
        well_potentials.resize(nw * np, 0.0);

        auto well_state_copy = well_state_;

        // average B factors are required for the convergence checking of well equations
        // Note: this must be done on all processes, even those with
        // no wells needing testing, otherwise we will have locking.
        std::vector< Scalar > B_avg(numComponents(), Scalar() );
        computeAverageFormationFactor(B_avg);

        const Opm::SummaryConfig& summaryConfig = ebosSimulator_.vanguard().summaryConfig();
        const bool write_restart_file = ebosSimulator_.vanguard().schedule().restart().getWriteRestartFile(reportStepIdx);
        int exception_thrown = 0;
        try {
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
                    std::vector<double> potentials;
                    well->computeWellPotentials(ebosSimulator_, B_avg, well_state_copy, potentials, deferred_logger);
                    // putting the sucessfully calculated potentials to the well_potentials
                    for (int p = 0; p < np; ++p) {
                        well_potentials[well->indexOfWell() * np + p] = std::abs(potentials[p]);
                    }
                }
            } // end of for (int w = 0; w < nw; ++w)
        } catch (std::exception& e) {
            exception_thrown = 1;
        }

        logAndCheckForExceptionsAndThrow(deferred_logger, exception_thrown, "computeWellPotentials() failed.", terminal_output_);

        // Store it in the well state
        well_state_.wellPotentials() = well_potentials;

    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    prepareTimeStep(Opm::DeferredLogger& deferred_logger)
    {
        int exception_thrown = 0;
        try {
            for (const auto& well : well_container_) {
                well->checkWellOperability(ebosSimulator_, well_state_, deferred_logger);
            }
            // since the controls are all updated, we should update well_state accordingly
            for (const auto& well : well_container_) {
                const int w = well->indexOfWell();
                if (!well->isOperable() ) continue;

                if (well_state_.effectiveEventsOccurred(w) ) {
                    well->updateWellStateWithTarget(ebosSimulator_, well_state_, deferred_logger);
                }

                // there is no new well control change input within a report step,
                // so next time step, the well does not consider to have effective events anymore
                // TODO: if we can know whether this is the first time step within the report step,
                // we do not need to set it to false
                // TODO: we should do this at the end of the time step in case we will need it within
                // this time step somewhere
                if (well_state_.effectiveEventsOccurred(w) ) {
                    well_state_.setEffectiveEventsOccurred(w, false);
                }
            }  // end of for (const auto& well : well_container_)
            updatePrimaryVariables(deferred_logger);
        } catch (std::exception& e) {
            exception_thrown = 1;
        }
        logAndCheckForExceptionsAndThrow(deferred_logger, exception_thrown, "prepareTimestep() failed.", terminal_output_);
    }





    template<typename TypeTag>
    const typename BlackoilWellModel<TypeTag>::WellState&
    BlackoilWellModel<TypeTag>::
    wellState() const { return well_state_; }

    template<typename TypeTag>
    const typename BlackoilWellModel<TypeTag>::WellState&
    BlackoilWellModel<TypeTag>::
    wellState(const WellState& well_state OPM_UNUSED) const { return wellState(); }


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
            for (unsigned i = 0; i < number_of_cells_; ++i) {
                cartesian_to_compressed_[global_cell[i]] = i;
            }
        }
        else {
            for (unsigned i = 0; i < number_of_cells_; ++i) {
                cartesian_to_compressed_[i] = i;
            }
        }

    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    computeRepRadiusPerfLength(const Grid& grid, Opm::DeferredLogger& deferred_logger)
    {
        for (const auto& well : well_container_) {
            well->computeRepRadiusPerfLength(grid, cartesian_to_compressed_, deferred_logger);
        }
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    computeAverageFormationFactor(std::vector<Scalar>& B_avg) const
    {
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
            if (has_solvent_) {
                auto& B  = B_avg[solventSaturationIdx];
                B += 1 / intQuants.solventInverseFormationVolumeFactor().value();
            }
        }

        // compute global average
        grid.comm().sum(B_avg.data(), B_avg.size());
        for(auto& bval: B_avg)
        {
            bval/=global_nc_;
        }
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    updatePrimaryVariables(Opm::DeferredLogger& deferred_logger)
    {
        for (const auto& well : well_container_) {
            well->updatePrimaryVariables(well_state_, deferred_logger);
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
        if (has_solvent_) {
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
            depth_[cellIdx] = Opm::UgGridHelpers::cellCenterDepth( grid, cellIdx );
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


    // convert well data from opm-common to well state from opm-core
    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    wellsToState( const data::Wells& wells,
                  const PhaseUsage& phases,
                  const bool handle_ms_well,
                  WellStateFullyImplicitBlackoil& state) const
    {

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

        for( const auto& wm : state.wellMap() ) {
            const auto well_index = wm.second[ 0 ];
            const auto& well = wells.at( wm.first );
            state.bhp()[ well_index ] = well.bhp;
            state.temperature()[ well_index ] = well.temperature;

            if (well.current_control.isProducer) {
                state.currentProductionControls()[ well_index ] = well.current_control.prod;
            }
            else {
                state.currentInjectionControls()[ well_index ] = well.current_control.inj;
            }

            const auto wellrate_index = well_index * np;
            for( size_t i = 0; i < phs.size(); ++i ) {
                assert( well.rates.has( phs[ i ] ) );
                state.wellRates()[ wellrate_index + i ] = well.rates.get( phs[ i ] );
            }

            const auto perforation_pressure = []( const data::Connection& comp ) {
                return comp.pressure;
            };

            const auto perforation_reservoir_rate = []( const data::Connection& comp ) {
                return comp.reservoir_rate;
            };
            std::transform( well.connections.begin(),
                            well.connections.end(),
                            state.perfPress().begin() + wm.second[ 1 ],
                    perforation_pressure );

            std::transform( well.connections.begin(),
                            well.connections.end(),
                            state.perfRates().begin() + wm.second[ 1 ],
                    perforation_reservoir_rate );

            int local_comp_index = 0;
            for (const data::Connection& comp : well.connections) {
                const int global_comp_index = wm.second[1] + local_comp_index;
                for (int phase_index = 0; phase_index < np; ++phase_index) {
                    state.perfPhaseRates()[global_comp_index*np + phase_index] = comp.rates.get(phs[phase_index]);
                }
                ++local_comp_index;
            }

            if (handle_ms_well && !well.segments.empty()) {
                // we need the well_ecl_ information
                const std::string& well_name = wm.first;
                const Well& well_ecl = getWellEcl(well_name);

                const WellSegments& segment_set = well_ecl.getSegments();

                const int top_segment_index = state.topSegmentIndex(well_index);
                const auto& segments = well.segments;

                // \Note: eventually we need to hanlde the situations that some segments are shut
                assert(0u + segment_set.size() == segments.size());

                for (const auto& segment : segments) {
                    const int segment_index = segment_set.segmentNumberToIndex(segment.first);

                    // recovering segment rates and pressure from the restart values
                    const auto pres_idx = Opm::data::SegmentPressures::Value::Pressure;
                    state.segPress()[top_segment_index + segment_index] = segment.second.pressures[pres_idx];

                    const auto& segment_rates = segment.second.rates;
                    for (int p = 0; p < np; ++p) {
                        state.segRates()[(top_segment_index + segment_index) * np + p] = segment_rates.get(phs[p]);
                    }
                }
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
    updateGroupIndividualControls(Opm::DeferredLogger& deferred_logger, std::set<std::string>& switched_groups)
    {
        const int reportStepIdx = ebosSimulator_.episodeIndex();
        const Group& fieldGroup = schedule().getGroup("FIELD", reportStepIdx);
        updateGroupIndividualControl(fieldGroup, deferred_logger, switched_groups);
    }



    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    updateGroupIndividualControl(const Group& group, Opm::DeferredLogger& deferred_logger, std::set<std::string>& switched_groups) {

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
    checkGroupConstraints(const Group& group, Opm::DeferredLogger& deferred_logger) const {

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
    checkGroupProductionConstraints(const Group& group, Opm::DeferredLogger& deferred_logger) const {

        const int reportStepIdx = ebosSimulator_.episodeIndex();
        const auto& summaryState = ebosSimulator_.vanguard().summaryState();
        const auto& comm = ebosSimulator_.vanguard().grid().comm();
        const auto& well_state = well_state_;

        const auto controls = group.productionControls(summaryState);
        const Group::ProductionCMode& currentControl = well_state.currentProductionGroupControl(group.name());

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
        const auto& well_state = well_state_;

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
        const Group::InjectionCMode& currentControl = well_state.currentInjectionGroupControl(phase, group.name());

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
        // Handle GCONSALE
        if (schedule().gConSale(reportStepIdx).has(group.name())) {

            if (controls.phase != Phase::GAS)
                OPM_THROW(std::runtime_error, "Group " + group.name() + " has GCONSALE control but is not a GAS group" );

            const auto& gconsale = schedule().gConSale(reportStepIdx).get(group.name(), summaryState);

            double sales_rate = 0.0;
            int gasPos = phase_usage_.phase_pos[BlackoilPhases::Vapour];
            sales_rate += WellGroupHelpers::sumWellRates(group, schedule(), well_state, reportStepIdx, gasPos, /*isInjector*/false);
            sales_rate -= WellGroupHelpers::sumWellRates(group, schedule(), well_state, reportStepIdx, gasPos, /*isInjector*/true);

            // sum over all nodes
            sales_rate = comm.sum(sales_rate);

            // add import rate and substract consumption rate for group for gas
            if (schedule().gConSump(reportStepIdx).has(group.name())) {
                const auto& gconsump = schedule().gConSump(reportStepIdx).get(group.name(), summaryState);
                if (phase_usage_.phase_used[BlackoilPhases::Vapour]) {
                    sales_rate += gconsump.import_rate;
                    sales_rate -= gconsump.consumption_rate;
                }
            }
            if (sales_rate > gconsale.max_sales_rate) {
                OPM_THROW(std::runtime_error, "Group " + group.name() + " has sale rate more then the maximum permitted value. Not implemented in Flow" );
            }
            if (sales_rate < gconsale.min_sales_rate) {
                OPM_THROW(std::runtime_error, "Group " + group.name() + " has sale rate less then minimum permitted value. Not implemented in Flow" );
            }
            if (gconsale.sales_target < 0.0) {
                OPM_THROW(std::runtime_error, "Group " + group.name() + " has sale rate target less then zero. Not implemented in Flow" );
            }

        }
        return Group::InjectionCMode::NONE;
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    actionOnBrokenConstraints(const Group& group, const Group::ExceedAction& exceed_action, const Group::ProductionCMode& newControl, Opm::DeferredLogger& deferred_logger) {

        auto& well_state = well_state_;
        const Group::ProductionCMode& oldControl = well_state.currentProductionGroupControl(group.name());

        std::ostringstream ss;

        switch(exceed_action) {
        case Group::ExceedAction::NONE: {
            if (oldControl != newControl && oldControl != Group::ProductionCMode::NONE) {
                ss << "Group production exceed limit is NONE for group " + group.name() + ". Nothing happens";
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
                well_state.setCurrentProductionGroupControl(group.name(), newControl);
                ss << "Switching control mode for group "<< group.name() << " to " << Group::ProductionCMode2String(newControl);
            }
            break;
        }
        default:
            throw("Invalid procedure for maximum rate limit selected for group" + group.name());
        }

        auto cc = Dune::MPIHelper::getCollectiveCommunication();
        if (cc.size() > 1) {
            ss << " on rank " << cc.rank();
        }
        if (!ss.str().empty())
            deferred_logger.info(ss.str());


    }


    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    actionOnBrokenConstraints(const Group& group, const Group::InjectionCMode& newControl, const Phase& controlPhase, Opm::DeferredLogger& deferred_logger) {
        auto& well_state = well_state_;
        const Group::InjectionCMode& oldControl = well_state.currentInjectionGroupControl(controlPhase, group.name());

        std::ostringstream ss;
        if (oldControl != newControl) {
            const std::string from = Group::InjectionCMode2String(oldControl);
            ss << "Group " << group.name() << " exceeding "
               << from << " limit \n";
            ss << "Switching control mode for group "<< group.name() << " to " << Group::InjectionCMode2String(newControl);
            auto cc = Dune::MPIHelper::getCollectiveCommunication();
            if (cc.size() > 1) {
                ss << " on rank " << cc.rank();
            }
            well_state.setCurrentInjectionGroupControl(controlPhase, group.name(), newControl);
        }

        if (!ss.str().empty())
            deferred_logger.info(ss.str());

    }




    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    updateGroupHigherControls(Opm::DeferredLogger& deferred_logger, std::set<std::string>& switched_groups)
    {
        const int reportStepIdx = ebosSimulator_.episodeIndex();
        const Group& fieldGroup = schedule().getGroup("FIELD", reportStepIdx);
        checkGroupHigherConstraints(fieldGroup, deferred_logger, switched_groups);
    }


    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    checkGroupHigherConstraints(const Group& group, Opm::DeferredLogger& deferred_logger, std::set<std::string>& switched_groups)
    {
        // Set up coefficients for RESV <-> surface rate conversion.
        // Use the pvtRegionIdx from the top cell of the first well.
        // TODO fix this!
        // This is only used for converting RESV rates.
        // What is the proper approach?
        const int fipnum = 0;
        const int pvtreg = well_perf_data_.empty()
            ? pvt_region_idx_[0]
            : pvt_region_idx_[well_perf_data_[0][0].cell_index];
        std::vector<double> resv_coeff(phase_usage_.num_phases, 0.0);
        rateConverter_->calcCoeff(fipnum, pvtreg, resv_coeff);

        const int reportStepIdx = ebosSimulator_.episodeIndex();
        const auto& summaryState = ebosSimulator_.vanguard().summaryState();

        std::vector<double> rates(phase_usage_.num_phases, 0.0);
        const auto& comm = ebosSimulator_.vanguard().grid().comm();

        const bool skip = switched_groups.count(group.name()) || group.name() == "FIELD";

        if (!skip && group.isInjectionGroup()) {
            // Obtain rates for group.
            for (int phasePos = 0; phasePos < phase_usage_.num_phases; ++phasePos) {
                const double local_current_rate = WellGroupHelpers::sumWellRates(
                    group, schedule(), well_state_, reportStepIdx, phasePos, /* isInjector */ true);
                // Sum over all processes
                rates[phasePos] = comm.sum(local_current_rate);
            }
            const Phase all[] = { Phase::WATER, Phase::OIL, Phase::GAS };
            for (Phase phase : all) {
                // Check higher up only if under individual (not FLD) control.
                const Group::InjectionCMode& currentControl = well_state_.currentInjectionGroupControl(phase, group.name());
                if (currentControl != Group::InjectionCMode::FLD) {
                    const Group& parentGroup = schedule().getGroup(group.parent(), reportStepIdx);
                    const std::pair<bool, double> changed = WellGroupHelpers::checkGroupConstraintsInj(
                        group.name(),
                        group.parent(),
                        parentGroup,
                        well_state_,
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
                const double local_current_rate = WellGroupHelpers::sumWellRates(
                    group, schedule(), well_state_, reportStepIdx, phasePos, /* isInjector */ false);
                // Sum over all processes
                rates[phasePos] = -comm.sum(local_current_rate);
            }
            // Check higher up only if under individual (not FLD) control.
            const Group::ProductionCMode& currentControl = well_state_.currentProductionGroupControl(group.name());
                if (currentControl != Group::ProductionCMode::FLD) {
                    const Group& parentGroup = schedule().getGroup(group.parent(), reportStepIdx);
                    const std::pair<bool, double> changed = WellGroupHelpers::checkGroupConstraintsProd(
                        group.name(),
                        group.parent(),
                        parentGroup,
                        well_state_,
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
    updateWsolvent(const Group& group, const Schedule& schedule, const int reportStepIdx, const WellStateFullyImplicitBlackoil& wellState) {
        for (const std::string& groupName : group.groups()) {
            const Group& groupTmp = schedule.getGroup(groupName, reportStepIdx);
            updateWsolvent(groupTmp, schedule, reportStepIdx, wellState);
        }

        if (group.isProductionGroup())
            return;

        const Group::InjectionCMode& currentGroupControl = wellState.currentInjectionGroupControl(Phase::GAS, group.name());
        if( currentGroupControl == Group::InjectionCMode::REIN ) {
            int gasPos = phase_usage_.phase_pos[BlackoilPhases::Vapour];
            double gasProductionRate = WellGroupHelpers::sumWellRates(group, schedule, wellState, reportStepIdx, gasPos, /*isInjector*/false);
            double solventProductionRate = WellGroupHelpers::sumSolventRates(group, schedule, wellState, reportStepIdx, /*isInjector*/false);

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


} // namespace Opm
