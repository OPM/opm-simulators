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
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    init()
    {
        const Opm::EclipseState& eclState = ebosSimulator_.vanguard().eclState();

        gravity_ = ebosSimulator_.problem().gravity()[2];

        extractLegacyCellPvtRegionIndex_();
        extractLegacyDepth_();

        phase_usage_ = phaseUsageFromDeck(eclState);

        const auto& gridView = ebosSimulator_.gridView();

        // calculate the number of elements of the compressed sequential grid. this needs
        // to be done in two steps because the dune communicator expects a reference as
        // argument for sum()
        number_of_cells_ = gridView.size(/*codim=*/0);
        global_nc_ = gridView.comm().sum(number_of_cells_);
        gravity_ = ebosSimulator_.problem().gravity()[2];

        extractLegacyCellPvtRegionIndex_();
        extractLegacyDepth_();
        initial_step_ = true;

        const auto& grid = ebosSimulator_.vanguard().grid();
        const auto& cartDims = Opm::UgGridHelpers::cartDims(grid);
        setupCartesianToCompressed_(Opm::UgGridHelpers::globalCell(grid),
                                    cartDims[0]*cartDims[1]*cartDims[2]);

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
        const auto& schedule_wells = schedule().getWells2atEnd();
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

        const Grid& grid = ebosSimulator_.vanguard().grid();
        const auto& defunct_well_names = ebosSimulator_.vanguard().defunctWellNames();
        const auto& eclState = ebosSimulator_.vanguard().eclState();
        const auto& summaryState = ebosSimulator_.vanguard().summaryState();
        wells_ecl_ = schedule().getWells2(timeStepIdx);

        // Create wells and well state.
        wells_manager_.reset( new WellsManager (eclState,
                                                schedule(),
                                                summaryState,
                                                timeStepIdx,
                                                Opm::UgGridHelpers::numCells(grid),
                                                Opm::UgGridHelpers::globalCell(grid),
                                                Opm::UgGridHelpers::cartDims(grid),
                                                Opm::UgGridHelpers::dimensions(grid),
                                                Opm::UgGridHelpers::cell2Faces(grid),
                                                Opm::UgGridHelpers::beginFaceCentroids(grid),
                                                grid.comm().size() > 1,
                                                defunct_well_names) );

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
        well_state_.init(wells(), cellPressures, schedule(), wells_ecl_, timeStepIdx, &previous_well_state_, phase_usage_);

        // handling MS well related
        if (param_.use_multisegment_well_&& anyMSWellOpenLocal(wells())) { // if we use MultisegmentWell model
            well_state_.initWellStateMSWell(wells(), wells_ecl_, phase_usage_, &previous_well_state_);
        }

        const int nw = wells()->number_of_wells;
        for (int w = 0; w <nw; ++w) {
            const int nw_wells_ecl = wells_ecl_.size();
            int index_well_ecl = 0;
            const std::string well_name(wells()->name[w]);
            for (; index_well_ecl < nw_wells_ecl; ++index_well_ecl) {
                if (well_name == wells_ecl_[index_well_ecl].name()) {
                    break;
                }
            }

            // It should be able to find in wells_ecl.
            if (index_well_ecl == nw_wells_ecl) {
                OPM_THROW(std::logic_error, "Could not find well " << well_name << " in wells_ecl ");
            }
            const auto well = wells_ecl_[index_well_ecl];
            const uint64_t effective_events_mask = ScheduleEvents::WELL_STATUS_CHANGE
                                                 + ScheduleEvents::PRODUCTION_UPDATE
                                                 + ScheduleEvents::INJECTION_UPDATE
                                                 + ScheduleEvents::NEW_WELL;

            if(!schedule().hasWellEvent(well_name, effective_events_mask, timeStepIdx))
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
        const Group2& fieldGroup = schedule().getGroup2("FIELD", timeStepIdx);
        wellGroupHelpers::setCmodeGroup(fieldGroup, schedule(), summaryState, timeStepIdx, well_state_);

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
            well_container_ = createWellContainer(reportStepIdx, local_deferredLogger);

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
        const int np = numPhases();
        for (const auto& well : well_container_) {
            const double& oilpot = well_state_.wellPotentials()[well->indexOfWell() * np + phase_usage_.phase_pos[BlackoilPhases::Liquid]];
            const double& gaspot = well_state_.wellPotentials()[well->indexOfWell() * np + phase_usage_.phase_pos[BlackoilPhases::Vapour]];
            const double& waterpot = well_state_.wellPotentials()[well->indexOfWell() * np + phase_usage_.phase_pos[BlackoilPhases::Aqua]];
            guideRate_->compute(well->name(), reportStepIdx, simulationTime, oilpot, gaspot, waterpot);
        }
        const Group2& fieldGroup = schedule().getGroup2("FIELD", reportStepIdx);
        wellGroupHelpers::updateGuideRateForGroups(fieldGroup, schedule(), phase_usage_, reportStepIdx, simulationTime, guideRate_.get(), well_state_);
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
                const Well2& wellEcl = schedule().getWell2(well_name, timeStepIdx);
                double well_efficiency_factor = wellEcl.getEfficiencyFactor();
                wellGroupHelpers::accumulateGroupEfficiencyFactor(schedule().getGroup2(wellEcl.groupName(), timeStepIdx), schedule(), timeStepIdx, well_efficiency_factor);
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

        Opm::DeferredLogger local_deferredLogger;
        for (const auto& well : well_container_) {
            if (GET_PROP_VALUE(TypeTag, EnablePolymerMW) && well->wellType() == INJECTOR) {
                well->updateWaterThroughput(dt, well_state_);
            }
        }
        updateWellTestState(simulationTime, wellTestState_);

        // update the rate converter with current averages pressures etc in
        rateConverter_->template defineState<ElementContext>(ebosSimulator_);

        const int reportStepIdx = ebosSimulator_.episodeIndex();
        // calculate the well potentials
        try {
            std::vector<double> well_potentials;
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
        const auto& defunctWellNames = ebosSimulator_.vanguard().defunctWellNames();

        // The restart step value is used to identify wells present at the given
        // time step. Wells that are added at the same time step as RESTART is initiated
        // will not be present in a restart file. Use the previous time step to retrieve
        // wells that have information written to the restart file.
        const int report_step = std::max(eclState().getInitConfig().getRestartStep() - 1, 0);
        const auto& summaryState = ebosSimulator_.vanguard().summaryState();

        WellsManager wellsmanager(eclState(),
                                  schedule(),
                                  summaryState,
                                  report_step,
                                  Opm::UgGridHelpers::numCells(grid()),
                                  Opm::UgGridHelpers::globalCell(grid()),
                                  Opm::UgGridHelpers::cartDims(grid()),
                                  Opm::UgGridHelpers::dimensions(grid()),
                                  Opm::UgGridHelpers::cell2Faces(grid()),
                                  Opm::UgGridHelpers::beginFaceCentroids(grid()),
                                  grid().comm().size() > 1,
                                  defunctWellNames);

        const Wells* wells = wellsmanager.c_wells();

        wells_ecl_ = schedule().getWells2(report_step);

        const int nw = wells->number_of_wells;
        if (nw > 0) {
            const auto phaseUsage = phaseUsageFromDeck(eclState());
            const size_t numCells = Opm::UgGridHelpers::numCells(grid());
            const bool handle_ms_well = (param_.use_multisegment_well_ && anyMSWellOpenLocal(wells));
            well_state_.resize(wells, wells_ecl_, schedule(), handle_ms_well, numCells, phaseUsage); // Resize for restart step
            wellsToState(restartValues.wells, phaseUsage, handle_ms_well, well_state_);
        }

        // for ecl compatible restart the current controls are not written
        const auto& ioCfg = eclState().getIOConfig();
        const auto ecl_compatible_rst = ioCfg.getEclCompatibleRST();        
        if (true || ecl_compatible_rst) { // always set the control from the schedule
            for (int w = 0; w <nw; ++w) {
                const int nw_wells_ecl = wells_ecl_.size();
                int index_well_ecl = 0;
                const std::string well_name(wells->name[w]);
                for (; index_well_ecl < nw_wells_ecl; ++index_well_ecl) {
                    if (well_name == wells_ecl_[index_well_ecl].name()) {
                        break;
                    }
                }

                // It should be able to find in wells_ecl.
                if (index_well_ecl == nw_wells_ecl) {
                    OPM_THROW(std::logic_error, "Could not find well " << well_name << " in wells_ecl ");
                }
                const auto& well = wells_ecl_[index_well_ecl];

                if (well.isProducer()) {
                    const auto controls = well.productionControls(summaryState);
                    well_state_.currentProductionControls()[w] = controls.cmode;
                }
                else {
                    const auto controls = well.injectionControls(summaryState);
                    well_state_.currentInjectionControls()[w] = controls.cmode;
                }
            }
        }

        previous_well_state_ = well_state_;

        initial_step_ = false;
    }





    template<typename TypeTag>
    std::vector<typename BlackoilWellModel<TypeTag>::WellInterfacePtr >
    BlackoilWellModel<TypeTag>::
    createWellContainer(const int time_step, Opm::DeferredLogger& deferred_logger)
    {
        std::vector<WellInterfacePtr> well_container;

        const int nw = numWells();

        if (nw > 0) {
            well_container.reserve(nw);

            // With the following way, it will have the same order with wells struct
            // Hopefully, it can generate the same residual history with master branch
            for (int w = 0; w < nw; ++w) {
                const std::string well_name = std::string(wells()->name[w]);

                // finding the location of the well in wells_ecl
                const int nw_wells_ecl = wells_ecl_.size();
                int index_well = 0;
                for (; index_well < nw_wells_ecl; ++index_well) {
                    if (well_name == wells_ecl_[index_well].name()) {
                        break;
                    }
                }

                // It should be able to find in wells_ecl.
                if (index_well == nw_wells_ecl) {
                    OPM_DEFLOG_THROW(std::runtime_error, "Could not find well " + well_name + " in wells_ecl ", deferred_logger);
                }

                const Well2& well_ecl = wells_ecl_[index_well];

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
                        // TODO: make a function from well_state side to handle the following
                        well_state_.thp()[w] = 0.;
                        well_state_.bhp()[w] = 0.;
                        const int np = numPhases();
                        for (int p = 0; p < np; ++p) {
                            well_state_.wellRates()[np * w + p] = 0.;
                        }
                        continue;
                    } else {
                        // stopped wells are added to the container but marked as stopped
                        wellIsStopped = true;
                    }
                }

                // Use the pvtRegionIdx from the top cell
                const int well_cell_top = wells()->well_cells[wells()->well_connpos[w]];
                const int pvtreg = pvt_region_idx_[well_cell_top];

                if ( !well_ecl.isMultiSegment() || !param_.use_multisegment_well_) {
                    well_container.emplace_back(new StandardWell<TypeTag>(well_ecl, time_step, wells(),
                                                param_, *rateConverter_, pvtreg, numComponents() ) );
                } else {
                    well_container.emplace_back(new MultisegmentWell<TypeTag>(well_ecl, time_step, wells(),
                                                param_, *rateConverter_, pvtreg, numComponents() ) );                    
                }
                if (wellIsStopped)
                    well_container.back()->stopWell();
            }
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

        const Well2& well_ecl = wells_ecl_[index_well_ecl];

        // Finding the location of the well in wells struct.
        const int nw = numWells();
        int well_index_wells = -999;
        for (int w = 0; w < nw; ++w) {
            if (well_name == std::string(wells()->name[w])) {
                well_index_wells = w;
                break;
            }
        }

        if (well_index_wells < 0) {
            OPM_DEFLOG_THROW(std::logic_error, "Could not find the well  " << well_name << " in the well struct ", deferred_logger);
        }

        // Use the pvtRegionIdx from the top cell
        const int well_cell_top = wells()->well_cells[wells()->well_connpos[well_index_wells]];
        const int pvtreg = pvt_region_idx_[well_cell_top];

        if ( !well_ecl.isMultiSegment() || !param_.use_multisegment_well_) {
             return WellInterfacePtr(new StandardWell<TypeTag>(well_ecl, report_step, wells(),
                                                 param_, *rateConverter_, pvtreg, numComponents() ) );
        } else {
             return WellInterfacePtr(new MultisegmentWell<TypeTag>(well_ecl, report_step, wells(),
                                                 param_, *rateConverter_, pvtreg, numComponents() ) );
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
            // only check group controls for iterationIdx smaller then nupcol
            const int reportStepIdx = ebosSimulator_.episodeIndex();
            const int nupcol = schedule().getNupcol(reportStepIdx);
            updateWellControls(local_deferredLogger, iterationIdx < nupcol);
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
        return numWells() > 0;
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
    getWellConvergence(const std::vector<Scalar>& B_avg) const
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
    updateWellControls(Opm::DeferredLogger& deferred_logger, const bool checkGroupControl)
    {
        // Even if there are no wells active locally, we cannot
        // return as the DeferredLogger uses global communication.
        // For no well active globally we simply return.
        if( !wellsActive() ) return ;

        // update group controls
        if (checkGroupControl) {
            const int reportStepIdx = ebosSimulator_.episodeIndex();
            const Group2& fieldGroup = schedule().getGroup2("FIELD", reportStepIdx);
            checkGroupConstraints(fieldGroup, deferred_logger);
        }

        for (const auto& well : well_container_) {
            well->updateWellControl(ebosSimulator_, well_state_, deferred_logger);
        }
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
        const int nw = numWells();
        const int np = numPhases();
        well_potentials.resize(nw * np, 0.0);

        auto well_state_copy = well_state_;

        // average B factors are required for the convergence checking of well equations
        // Note: this must be done on all processes, even those with
        // no wells needing testing, otherwise we will have locking.
        std::vector< Scalar > B_avg(numComponents(), Scalar() );
        computeAverageFormationFactor(B_avg);

        const Opm::SummaryConfig& summaryConfig = ebosSimulator_.vanguard().summaryConfig();
        const bool write_restart_file = ebosSimulator_.vanguard().eclState().getRestartConfig().getWriteRestartFile(reportStepIdx);
        int exception_thrown = 0;
        try {
            for (const auto& well : well_container_) {
                const bool needed_for_summary = ((summaryConfig.hasSummaryKey( "WWPI:" + well->name()) ||
                                                  summaryConfig.hasSummaryKey( "WOPI:" + well->name()) ||
                                                  summaryConfig.hasSummaryKey( "WGPI:" + well->name())) && well->wellType() == INJECTOR) ||
                                                ((summaryConfig.hasSummaryKey( "WWPP:" + well->name()) ||
                                                  summaryConfig.hasSummaryKey( "WOPP:" + well->name()) ||
                                                  summaryConfig.hasSummaryKey( "WGPP:" + well->name())) && well->wellType() == PRODUCER);

                const Well2& eclWell = well->wellEcl();
                bool needPotentialsForGuideRate = eclWell.getGuideRatePhase() == Well2::GuideRateTarget::UNDEFINED;
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
            const Well2& wellEcl = well->wellEcl();
            double well_efficiency_factor = wellEcl.getEfficiencyFactor();
            wellGroupHelpers::accumulateGroupEfficiencyFactor(schedule().getGroup2(wellEcl.groupName(), reportStepIdx), schedule(), reportStepIdx, well_efficiency_factor);
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
        if (numWells() > 0 && numPhases() < 3) {
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
    BlackoilWellModel<TypeTag>:: numWells() const
    {
        return wells() ? wells()->number_of_wells : 0;
    }

    template<typename TypeTag>
    int
    BlackoilWellModel<TypeTag>:: numPhases() const
    {
        return wells() ? wells()->number_of_phases : 1;
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::extractLegacyDepth_()
    {
        const auto& grid = ebosSimulator_.vanguard().grid();
        const unsigned numCells = grid.size(/*codim=*/0);

        depth_.resize(numCells);
        for (unsigned cellIdx = 0; cellIdx < numCells; ++cellIdx) {
            depth_[cellIdx] =
                grid.cellCenterDepth(cellIdx);
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

            //state.currentInjectionControls()[ well_index ] = static_cast<Opm::Well2::InjectorCMode>(well.injectionControl);
            //state.currentProductionControls()[ well_index ] = static_cast<Well2::ProducerCMode>(well.productionControl);

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
                const Well2& well_ecl = getWellEcl(well_name);

                const WellSegments& segment_set = well_ecl.getSegments();

                const int top_segment_index = state.topSegmentIndex(well_index);
                const auto& segments = well.segments;

                // \Note: eventually we need to hanlde the situations that some segments are shut
                assert(0u + segment_set.size() == segments.size());

                for (const auto& segment : segments) {
                    const int segment_index = segment_set.segmentNumberToIndex(segment.first);

                    // recovering segment rates and pressure from the restart values
                    state.segPress()[top_segment_index + segment_index] = segment.second.pressure;

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
    anyMSWellOpenLocal(const Wells* wells) const
    {
        bool any_ms_well_open = false;

        const int nw = wells->number_of_wells;
        for (int w = 0; w < nw; ++w) {
            const std::string well_name = std::string(wells->name[w]);

            const Well2& well_ecl = getWellEcl(well_name);

            if (well_ecl.isMultiSegment() ) {
                any_ms_well_open = true;
                break;
            }
        }
        return any_ms_well_open;
    }





    template<typename TypeTag>
    const Well2&
    BlackoilWellModel<TypeTag>::
    getWellEcl(const std::string& well_name) const
    {
        // finding the iterator of the well in wells_ecl
        auto well_ecl = std::find_if(wells_ecl_.begin(),
                                     wells_ecl_.end(),
                                     [&well_name](const Well2& elem)->bool {
                                         return elem.name() == well_name;
                                     });

        assert(well_ecl != wells_ecl_.end());

        return *well_ecl;
    }
    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    checkGroupConstraints(const Group2& group, Opm::DeferredLogger& deferred_logger) {

        // call recursively
        const int reportStepIdx = ebosSimulator_.episodeIndex();
        for (const std::string& groupName : group.groups()) {
            checkGroupConstraints( schedule().getGroup2(groupName, reportStepIdx), deferred_logger);
        }

        const auto& summaryState = ebosSimulator_.vanguard().summaryState();
        auto& well_state = well_state_;

        if (group.isInjectionGroup())
        {
            const auto controls = group.injectionControls(summaryState);
            int phasePos;
            switch (controls.phase) {
            case Phase::WATER:
            {
                phasePos = phase_usage_.phase_pos[BlackoilPhases::Aqua];
                break;
            }
            case Phase::OIL:
            {
                phasePos = phase_usage_.phase_pos[BlackoilPhases::Liquid];
                break;
            }
            case Phase::GAS:
            {
                phasePos = phase_usage_.phase_pos[BlackoilPhases::Vapour];
                break;
            }
            default:
                throw("Expected WATER, OIL or GAS as type for group injector: " + group.name());
            }

            if (group.has_control(Group2::InjectionCMode::NONE))
            {
                // do nothing??
            }

            if (group.has_control(Group2::InjectionCMode::RATE))
            {
                double current_rate = 0.0;
                current_rate += wellGroupHelpers::sumWellRates(group, schedule(), well_state, reportStepIdx, phasePos, /*isInjector*/true);
                if (controls.surface_max_rate < current_rate) {
                    actionOnBrokenConstraints(group, Group2::InjectionCMode::RATE, reportStepIdx, deferred_logger);
                }
            }
            if (group.has_control(Group2::InjectionCMode::RESV))
            {
                double current_rate = 0.0;
                current_rate += wellGroupHelpers::sumWellResRates(group, schedule(), well_state, reportStepIdx, phasePos, /*isInjector*/true);
                if (controls.resv_max_rate < current_rate) {
                    actionOnBrokenConstraints(group, Group2::InjectionCMode::RESV, reportStepIdx, deferred_logger);
                }                    }
            if (group.has_control(Group2::InjectionCMode::REIN))
            {
                double production_Rate = 0.0;
                production_Rate += wellGroupHelpers::sumWellRates(group, schedule(), well_state, reportStepIdx, phasePos, /*isInjector*/false);

                double current_rate = 0.0;
                current_rate += wellGroupHelpers::sumWellRates(group, schedule(), well_state, reportStepIdx, phasePos, /*isInjector*/true);

                if (controls.target_reinj_fraction*production_Rate < current_rate) {
                    actionOnBrokenConstraints(group, Group2::InjectionCMode::REIN, reportStepIdx, deferred_logger);
                }                    }
            if (group.has_control(Group2::InjectionCMode::VREP))
            {
                double voidage_Rate = 0.0;
                voidage_Rate += wellGroupHelpers::sumWellResRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Aqua], false);
                voidage_Rate += wellGroupHelpers::sumWellResRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Liquid], false);
                voidage_Rate += wellGroupHelpers::sumWellResRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Vapour], false);


                double current_rate = 0.0;
                current_rate += wellGroupHelpers::sumWellResRates(group, schedule(), well_state, reportStepIdx, phasePos, /*isInjector*/true);

                if (controls.target_void_fraction*voidage_Rate < current_rate) {
                    actionOnBrokenConstraints(group, Group2::InjectionCMode::VREP, reportStepIdx, deferred_logger);
                }
            }
            if (group.has_control(Group2::InjectionCMode::FLD))
            {
                // do nothing???
                //OPM_THROW(std::runtime_error, "Group " + group.name() + "FLD control for injecting groups not implemented" );
            }

        } else if (group.isProductionGroup())
        {
            const auto controls = group.productionControls(summaryState);

            if (group.has_control(Group2::ProductionCMode::NONE))
            {

            }
            if (group.has_control(Group2::ProductionCMode::ORAT))
            {
                double current_rate = 0.0;
                current_rate += wellGroupHelpers::sumWellRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Liquid], false);
                if (controls.oil_target < current_rate  ) {
                    actionOnBrokenConstraints(group, controls.exceed_action, Group2::ProductionCMode::ORAT, reportStepIdx, deferred_logger);
                }
            }

            if (group.has_control(Group2::ProductionCMode::WRAT))
            {
                double current_rate = 0.0;
                current_rate += wellGroupHelpers::sumWellRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Aqua], false);

                if (controls.water_target < current_rate  ) {
                    actionOnBrokenConstraints(group, controls.exceed_action, Group2::ProductionCMode::WRAT, reportStepIdx, deferred_logger);
                }
            }
            if (group.has_control(Group2::ProductionCMode::GRAT))
            {
                double current_rate = 0.0;
                current_rate += wellGroupHelpers::sumWellRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Vapour], false);
                if (controls.gas_target < current_rate  ) {
                    actionOnBrokenConstraints(group, controls.exceed_action, Group2::ProductionCMode::GRAT, reportStepIdx, deferred_logger);
                }
            }
            if (group.has_control(Group2::ProductionCMode::LRAT))
            {
                double current_rate = 0.0;
                current_rate += wellGroupHelpers::sumWellRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Liquid], false);
                current_rate += wellGroupHelpers::sumWellRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Aqua], false);
                if (controls.liquid_target < current_rate  ) {
                    actionOnBrokenConstraints(group, controls.exceed_action, Group2::ProductionCMode::LRAT, reportStepIdx, deferred_logger);
                }
            }

            if (group.has_control(Group2::ProductionCMode::CRAT))
            {
                OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "CRAT control for production groups not implemented" , deferred_logger);

            }
            if (group.has_control(Group2::ProductionCMode::RESV))
            {
                double current_rate = 0.0;
                current_rate += wellGroupHelpers::sumWellResRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Aqua], true);
                current_rate += wellGroupHelpers::sumWellResRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Liquid], true);
                current_rate += wellGroupHelpers::sumWellResRates(group, schedule(), well_state, reportStepIdx, phase_usage_.phase_pos[BlackoilPhases::Vapour], true);

                if (controls.resv_target < current_rate  ) {
                    actionOnBrokenConstraints(group, controls.exceed_action, Group2::ProductionCMode::RESV, reportStepIdx, deferred_logger);
                }

            }
            if (group.has_control(Group2::ProductionCMode::PRBL))
            {
                OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "PRBL control for production groups not implemented", deferred_logger);
            }
            if (group.has_control(Group2::ProductionCMode::FLD))
            {
                // do nothing???
            }
        } else {

            //neither production or injecting group FIELD?
        }



    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    actionOnBrokenConstraints(const Group2& group, const Group2::ExceedAction& exceed_action, const Group2::ProductionCMode& newControl, const int reportStepIdx, Opm::DeferredLogger& deferred_logger) {

        auto& well_state = well_state_;       
        const Group2::ProductionCMode& oldControl = well_state.currentProductionGroupControl(group.name());
        const std::string from = Group2::ProductionCMode2String(oldControl);

        std::ostringstream ss;
        ss << "Group " << group.name() << " exceeding "
           << from << " limit  \n";
        switch(exceed_action) {
        case Group2::ExceedAction::NONE: {
            OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "GroupProductionExceedLimit NONE not implemented", deferred_logger);
            break;
        }
        case Group2::ExceedAction::CON: {
            OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "GroupProductionExceedLimit CON not implemented", deferred_logger);
            break;
        }
        case Group2::ExceedAction::CON_PLUS: {
            OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "GroupProductionExceedLimit CON_PLUS not implemented", deferred_logger);
            break;
        }
        case Group2::ExceedAction::WELL: {
            OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "GroupProductionExceedLimit WELL not implemented", deferred_logger);
            break;
        }
        case Group2::ExceedAction::PLUG: {
            OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "GroupProductionExceedLimit PLUG not implemented", deferred_logger);
            break;
        }
        case Group2::ExceedAction::RATE: {
            well_state.setCurrentProductionGroupControl(group.name(), newControl);
            ss << "Switching control mode for group to " << Group2::ProductionCMode2String(newControl)
               << " \n Wells in group " + group.name() + " switches to GRUP control limit";
            wellGroupHelpers::setGroupControl(group, schedule(), reportStepIdx, false, well_state);
            break;
        }
        default:
            throw("Invalid procedure for maximum rate limit selected for group" + group.name());
        }

        auto cc = Dune::MPIHelper::getCollectiveCommunication();
        if (cc.size() > 1) {
            ss << " on rank " << cc.rank();
        }
        deferred_logger.info(ss.str());


    }


    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    actionOnBrokenConstraints(const Group2& group, const Group2::InjectionCMode& newControl, const int reportStepIdx, Opm::DeferredLogger& deferred_logger) {
        auto& well_state = well_state_;
        const Group2::InjectionCMode& oldControl = well_state.currentInjectionGroupControl(group.name());
        const std::string from = Group2::InjectionCMode2String(oldControl);
        std::ostringstream ss;
        ss << "Group " << group.name() << " exceeding "
           << from << " limit \n";
        ss << "Switching control mode for group to " << Group2::InjectionCMode2String(newControl)
           << " \n Wells in group " + group.name() + " switches to GRUP control limit";
        auto cc = Dune::MPIHelper::getCollectiveCommunication();
        if (cc.size() > 1) {
            ss << " on rank " << cc.rank();
        }
        deferred_logger.info(ss.str());
        well_state.setCurrentInjectionGroupControl(group.name(), newControl);
        wellGroupHelpers::setGroupControl(group, schedule(), reportStepIdx, /*isInjector*/true, well_state);
    }


} // namespace Opm
