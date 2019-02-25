


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
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    init(const Opm::EclipseState& eclState, const Opm::Schedule& schedule)
    {
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
        int last_time_step = schedule().getTimeMap().size() - 1;
        const auto& schedule_wells = schedule().getWells();
        const auto& cartesianSize = Opm::UgGridHelpers::cartDims(grid());

        // initialize the additional cell connections introduced by wells.
        for (const auto well : schedule_wells)
        {
            std::vector<int> wellCells;
            // All possible connections of the well
            const auto& connectionSet = well->getConnections(last_time_step);
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
    linearize(SparseMatrixAdapter& mat , GlobalEqVector& res)
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
            well->addWellContributions(mat.istlMatrix());

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
        for (const auto& well : well_container_) {
            if (well->wellHasTHPConstraints()) {
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
            if (well->name() == wellname) {
                if (well->underPredictionMode()) {
                    wellTestState_.addClosedWell(wellname, WellTestConfig::Reason::PHYSICAL, simulation_time);
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
        const Grid& grid = ebosSimulator_.vanguard().grid();
        const auto& defunct_well_names = ebosSimulator_.vanguard().defunctWellNames();
        const auto& eclState = ebosSimulator_.vanguard().eclState();
        wells_ecl_ = schedule().getWells(timeStepIdx);

        // Create wells and well state.
        wells_manager_.reset( new WellsManager (eclState,
                                                schedule(),
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

            const double p = fs.pressure(FluidSystem::oilPhaseIdx).value();
            cellPressures[cellIdx] = p;
        }
        well_state_.init(wells(), cellPressures, wells_ecl_, timeStepIdx, &previous_well_state_, phase_usage_);

        // handling MS well related
        if (param_.use_multisegment_well_) { // if we use MultisegmentWell model
            for (const auto& well : wells_ecl_) {
                // TODO: this is acutally not very accurate, because sometimes a deck just claims a MS well
                // while keep the well shut. More accurately, we should check if the well exisits in the Wells
                // structure here
                if (well->isMultiSegment(timeStepIdx) ) { // there is one well is MS well
                    well_state_.initWellStateMSWell(wells(), wells_ecl_, timeStepIdx, phase_usage_, previous_well_state_);
                    break;
                }
            }
        }

        // update the previous well state. This is used to restart failed steps.
        previous_well_state_ = well_state_;

        // Compute reservoir volumes for RESV controls.
        rateConverter_.reset(new RateConverterType (phase_usage_,
                                                    std::vector<int>(number_of_cells_, 0)));
        computeRESV(timeStepIdx);

        // update VFP properties
        vfp_properties_.reset (new VFPProperties<VFPInjProperties,VFPProdProperties> (
                                   schedule().getVFPInjTables(timeStepIdx),
                                   schedule().getVFPProdTables(timeStepIdx)) );



    }


    // called at the beginning of a time step
    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    beginTimeStep() {
        well_state_ = previous_well_state_;

        const int reportStepIdx = ebosSimulator_.episodeIndex();
        const double simulationTime = ebosSimulator_.time();

        // test wells
        wellTesting(reportStepIdx, simulationTime);

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
        calculateEfficiencyFactors();

        if (has_polymer_)
        {
            const Grid& grid = ebosSimulator_.vanguard().grid();
            if (PolymerModule::hasPlyshlog() || GET_PROP_VALUE(TypeTag, EnablePolymerMW) ) {
                computeRepRadiusPerfLength(grid);
            }
        }

        for (auto& well : well_container_) {
            well->setVFPProperties(vfp_properties_.get());
        }

        // Close completions due to economical reasons
        for (auto& well : well_container_) {
            well->closeCompletions(wellTestState_);
        }

    }


    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::wellTesting(const int timeStepIdx, const double simulationTime) {
        Opm::DeferredLogger local_deferredLogger;
        const auto& wtest_config = schedule().wtestConfig(timeStepIdx);
        if (wtest_config.size() != 0) { // there is a WTEST request

            // average B factors are required for the convergence checking of well equations
            // Note: this must be done on all processes, even those with
            // no wells needing testing, otherwise we will have locking.
            std::vector< Scalar > B_avg(numComponents(), Scalar() );
            computeAverageFormationFactor(B_avg);

            const auto& wellsForTesting = wellTestState_.updateWell(wtest_config, simulationTime);
            for (const auto& testWell : wellsForTesting) {
                const std::string& well_name = testWell.first;

                // this is the well we will test
                WellInterfacePtr well = createWellForWellTest(well_name, timeStepIdx);

                // some preparation before the well can be used
                well->init(&phase_usage_, depth_, gravity_, number_of_cells_);
                const WellNode& well_node = wellCollection().findWellNode(well_name);
                const double well_efficiency_factor = well_node.getAccumulativeEfficiencyFactor();
                well->setWellEfficiencyFactor(well_efficiency_factor);
                well->setVFPProperties(vfp_properties_.get());

                const WellTestConfig::Reason testing_reason = testWell.second;

                well->wellTesting(ebosSimulator_, B_avg, simulationTime, timeStepIdx,
                                  testing_reason, well_state_, wellTestState_, local_deferredLogger);
            }
        }
        Opm::DeferredLogger global_deferredLogger = gatherDeferredLogger(local_deferredLogger);
        if (terminal_output_) {
            global_deferredLogger.logMessages();
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
        // TODO: when necessary
        rateConverter_->template defineState<ElementContext>(ebosSimulator_);
        for (const auto& well : well_container_) {
            well->calculateReservoirRates(well_state_);
            if (GET_PROP_VALUE(TypeTag, EnablePolymerMW) && well->wellType() == INJECTOR) {
                well->updateWaterThroughput(dt, well_state_);
            }
        }
        updateWellTestState(simulationTime, wellTestState_);

        // calculate the well potentials for output
        // TODO: when necessary
        try
        {
            std::vector<double> well_potentials;
            computeWellPotentials(well_potentials);
        }
        catch ( std::runtime_error& e )
        {
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
        WellsManager wellsmanager(eclState(),
                                  schedule(),
                                  // The restart step value is used to identify wells present at the given
                                  // time step. Wells that are added at the same time step as RESTART is initiated
                                  // will not be present in a restart file. Use the previous time step to retrieve
                                  // wells that have information written to the restart file.
                                  std::max(eclState().getInitConfig().getRestartStep() - 1, 0),
                                  Opm::UgGridHelpers::numCells(grid()),
                                  Opm::UgGridHelpers::globalCell(grid()),
                                  Opm::UgGridHelpers::cartDims(grid()),
                                  Opm::UgGridHelpers::dimensions(grid()),
                                  Opm::UgGridHelpers::cell2Faces(grid()),
                                  Opm::UgGridHelpers::beginFaceCentroids(grid()),
                                  grid().comm().size() > 1,
                                  defunctWellNames);

        const Wells* wells = wellsmanager.c_wells();

        const int nw = wells->number_of_wells;
        if (nw > 0) {
            const auto phaseUsage = phaseUsageFromDeck(eclState());
            const size_t numCells = Opm::UgGridHelpers::numCells(grid());
            well_state_.resize(wells, numCells, phaseUsage); // Resize for restart step
            wellsToState(restartValues.wells, phaseUsage, well_state_);
            previous_well_state_ = well_state_;
        }
        initial_step_ = false;
    }





    template<typename TypeTag>
    std::vector<typename BlackoilWellModel<TypeTag>::WellInterfacePtr >
    BlackoilWellModel<TypeTag>::
    createWellContainer(const int time_step)
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
                    if (well_name == wells_ecl_[index_well]->name()) {
                        break;
                    }
                }

                // It should be able to find in wells_ecl.
                if (index_well == nw_wells_ecl) {
                    OPM_THROW(std::logic_error, "Could not find well " << well_name << " in wells_ecl ");
                }

                const Well* well_ecl = wells_ecl_[index_well];

                // A new WCON keywords can re-open a well that was closed/shut due to Physical limit
                if ( wellTestState_.hasWell(well_name, WellTestConfig::Reason::PHYSICAL ) ) {
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
                if ( wellTestState_.hasWell(well_name, WellTestConfig::Reason::ECONOMIC) ||
                     wellTestState_.hasWell(well_name, WellTestConfig::Reason::PHYSICAL) ) {
                    if( well_ecl->getAutomaticShutIn() ) {
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
                        // close wells are added to the container but marked as closed
                        struct WellControls* well_controls = wells()->ctrls[w];
                        well_controls_stop_well(well_controls);
                    }
                }

                // Use the pvtRegionIdx from the top cell
                const int well_cell_top = wells()->well_cells[wells()->well_connpos[w]];
                const int pvtreg = pvt_region_idx_[well_cell_top];

                if ( !well_ecl->isMultiSegment(time_step) || !param_.use_multisegment_well_) {
                    if ( GET_PROP_VALUE(TypeTag, EnablePolymerMW) && well_ecl->isInjector(time_step) ) {
                        well_container.emplace_back(new StandardWellV<TypeTag>(well_ecl, time_step, wells(),
                                                    param_, *rateConverter_, pvtreg, numComponents() ) );
                    } else {
                        well_container.emplace_back(new StandardWell<TypeTag>(well_ecl, time_step, wells(),
                                                    param_, *rateConverter_, pvtreg, numComponents() ) );
                    }
                } else {
                    well_container.emplace_back(new MultisegmentWell<TypeTag>(well_ecl, time_step, wells(),
                                                param_, *rateConverter_, pvtreg, numComponents() ) );
                }
            }
        }
        return well_container;
    }





    template<typename TypeTag>
    typename BlackoilWellModel<TypeTag>::WellInterfacePtr
    BlackoilWellModel<TypeTag>::
    createWellForWellTest(const std::string& well_name,
                          const int report_step) const
    {
        // Finding the location of the well in wells_ecl
        const int nw_wells_ecl = wells_ecl_.size();
        int index_well_ecl = 0;
        for (; index_well_ecl < nw_wells_ecl; ++index_well_ecl) {
            if (well_name == wells_ecl_[index_well_ecl]->name()) {
                break;
            }
        }
        // It should be able to find in wells_ecl.
        if (index_well_ecl == nw_wells_ecl) {
            OPM_THROW(std::logic_error, "Could not find well " << well_name << " in wells_ecl ");
        }

        const Well* well_ecl = wells_ecl_[index_well_ecl];

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
            OPM_THROW(std::logic_error, "Could not find the well  " << well_name << " in the well struct ");
        }

        // Use the pvtRegionIdx from the top cell
        const int well_cell_top = wells()->well_cells[wells()->well_connpos[well_index_wells]];
        const int pvtreg = pvt_region_idx_[well_cell_top];

        if ( !well_ecl->isMultiSegment(report_step) || !param_.use_multisegment_well_) {
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

        if (iterationIdx == 0) {
            calculateExplicitQuantities();
            prepareTimeStep(local_deferredLogger);
        }

        updateWellControls();
        // Set the well primary variables based on the value of well solutions
        initPrimaryVariablesEvaluation();

        if (param_.solve_welleq_initially_ && iterationIdx == 0) {
            // solve the well equations as a pre-processing step
            last_report_ = solveWellEq(dt, local_deferredLogger);

            if (initial_step_) {
                // update the explicit quantities to get the initial fluid distribution in the well correct.
                calculateExplicitQuantities();
                prepareTimeStep(local_deferredLogger);
                last_report_ = solveWellEq(dt, local_deferredLogger);
                initial_step_ = false;
            }
            // TODO: should we update the explicit related here again, or even prepareTimeStep().
            // basically, this is a more updated state from the solveWellEq based on fixed
            // reservoir state, will tihs be a better place to inialize the explict information?
        }
        assembleWellEq(dt, local_deferredLogger);

        Opm::DeferredLogger global_deferredLogger = gatherDeferredLogger(local_deferredLogger);
        if (terminal_output_) {
            global_deferredLogger.logMessages();
        }


        last_report_.converged = true;
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    assembleWellEq(const double dt, Opm::DeferredLogger& deferred_logger)
    {
        for (auto& well : well_container_) {
            well->assembleWellEq(ebosSimulator_, dt, well_state_, deferred_logger);
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
        if (!localWellsActive())
            return;

        for (auto& well : well_container_) {
            well->recoverWellSolutionAndUpdateWellState(x, well_state_);
        }
    }




    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    resetWellControlFromState() const
    {
        for (auto& well : well_container_) {
            WellControls* wc = well->wellControls();
            well_controls_set_current( wc, well_state_.currentControls()[well->indexOfWell()]);
        }
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
    solveWellEq(const double dt, Opm::DeferredLogger& deferred_logger)
    {
        WellState well_state0 = well_state_;

        const int numComp = numComponents();
        std::vector< Scalar > B_avg( numComp, Scalar() );
        computeAverageFormationFactor(B_avg);

        const int max_iter = param_.max_welleq_iter_;

        int it  = 0;
        bool converged;
        do {
            assembleWellEq(dt, deferred_logger);

            const auto report = getWellConvergence(B_avg);
            converged = report.converged();

            // checking whether the group targets are converged
            if (wellCollection().groupControlActive()) {
                converged = converged && wellCollection().groupTargetConverged(well_state_.wellRates());
            }

            if (converged) {
                break;
            }

            ++it;
            if( localWellsActive() )
            {
                for (auto& well : well_container_) {
                    well->solveEqAndUpdateWellState(well_state_);
                }
            }
            // updateWellControls uses communication
            // Therefore the following is executed if there
            // are active wells anywhere in the global domain.
            if( wellsActive() )
            {
                updateWellControls();
                initPrimaryVariablesEvaluation();
            }
        } while (it < max_iter);

        if (converged) {
            if (terminal_output_) {
                deferred_logger.debug("Well equation solution gets converged with " + std::to_string(it) + " iterations");
            }
        } else {
            if (terminal_output_) {
                deferred_logger.debug("Well equation solution failed in getting converged with " + std::to_string(it) + " iterations");
            }

            well_state_ = well_state0;
            updatePrimaryVariables();
            // also recover the old well controls
            for (const auto& well : well_container_) {
                const int index_of_well = well->indexOfWell();
                WellControls* wc = well->wellControls();
                well_controls_set_current(wc, well_state_.currentControls()[index_of_well]);
            }
        }

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
        // Get global (from all processes) convergence report.
        ConvergenceReport local_report;
        for (const auto& well : well_container_) {
            if (well->isOperable() ) {
                local_report += well->getWellConvergence(B_avg);
            }
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
    calculateExplicitQuantities() const
    {
        // TODO: checking isOperable() ?
        for (auto& well : well_container_) {
            well->calculateExplicitQuantities(ebosSimulator_, well_state_);
        }
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    updateWellControls()
    {
        // Even if there are no wells active locally, we cannot
        // return as the DeferredLogger uses global communication.
        // For no well active globally we simply return.
        if( !wellsActive() ) return ;

        Opm::DeferredLogger local_deferredLogger;

        for (const auto& well : well_container_) {
            well->updateWellControl(ebosSimulator_, well_state_, local_deferredLogger);
        }

        updateGroupControls(local_deferredLogger);

        Opm::DeferredLogger global_deferredLogger = gatherDeferredLogger(local_deferredLogger);
        if (terminal_output_) {
            global_deferredLogger.logMessages();
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
    computeWellPotentials(std::vector<double>& well_potentials)
    {
        Opm::DeferredLogger local_deferredLogger;
        // number of wells and phases
        const int nw = numWells();
        const int np = numPhases();
        well_potentials.resize(nw * np, 0.0);

        const Opm::SummaryConfig& summaryConfig = ebosSimulator_.vanguard().summaryConfig();
        for (const auto& well : well_container_) {
            // Only compute the well potential when asked for
            bool needed_for_output = ((summaryConfig.hasSummaryKey( "WWPI:" + well->name()) ||
                                       summaryConfig.hasSummaryKey( "WOPI:" + well->name()) ||
                                       summaryConfig.hasSummaryKey( "WGPI:" + well->name())) && well->wellType() == INJECTOR) ||
                                    ((summaryConfig.hasSummaryKey( "WWPP:" + well->name()) ||
                                                       summaryConfig.hasSummaryKey( "WOPP:" + well->name()) ||
                                                       summaryConfig.hasSummaryKey( "WGPP:" + well->name())) && well->wellType() == PRODUCER);

            if (needed_for_output || wellCollection().requireWellPotentials())
            {
                std::vector<double> potentials;
                well->computeWellPotentials(ebosSimulator_, well_state_, potentials, local_deferredLogger);

                // putting the sucessfully calculated potentials to the well_potentials
                for (int p = 0; p < np; ++p) {
                    well_potentials[well->indexOfWell() * np + p] = std::abs(potentials[p]);
                }
            }
        } // end of for (int w = 0; w < nw; ++w)

        // Store it in the well state
        well_state_.wellPotentials() = well_potentials;

        Opm::DeferredLogger global_deferredLogger = gatherDeferredLogger(local_deferredLogger);
        if (terminal_output_) {
            global_deferredLogger.logMessages();
        }

    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    prepareTimeStep(Opm::DeferredLogger& deferred_logger)
    {

        if ( wellCollection().havingVREPGroups() ) {
            rateConverter_->template defineState<ElementContext>(ebosSimulator_);
        }

        // after restarting, the well_controls can be modified while
        // the well_state still uses the old control index
        // we need to synchronize these two.
        // keep in mind that we set the control index of well_state to be the same with
        // with the wellControls from the deck when we create well_state at the beginning of the report step
        resetWellControlFromState();

        // process group control related
        prepareGroupControl();

        for (const auto& well : well_container_) {
            well->checkWellOperability(ebosSimulator_, well_state_, deferred_logger);
        }

        // since the controls are all updated, we should update well_state accordingly
        for (const auto& well : well_container_) {
            const int w = well->indexOfWell();
            WellControls* wc = well->wellControls();
            const int control = well_controls_get_current(wc);
            well_state_.currentControls()[w] = control;

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

        updatePrimaryVariables();
    }








    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    prepareGroupControl()
    {
        // group control related processing
        if (wellCollection().groupControlActive()) {
            for (const auto& well : well_container_) {
                WellControls* wc = well->wellControls();
                WellNode& well_node = wellCollection().findWellNode(well->name());

                // handling the situation that wells do not have a valid control
                // it happens the well specified with GRUP and restarting due to non-convergencing
                // putting the well under group control for this situation
                int ctrl_index = well_controls_get_current(wc);

                const int group_control_index = well_node.groupControlIndex();
                if (group_control_index >= 0 && ctrl_index < 0) {
                    // put well under group control
                    well_controls_set_current(wc, group_control_index);
                    well_state_.currentControls()[well->indexOfWell()] = group_control_index;
                }

                // Final step, update whehter the well is under group control or individual control
                // updated ctrl_index from the well control
                ctrl_index = well_controls_get_current(wc);
                if (well_node.groupControlIndex() >= 0 && ctrl_index == well_node.groupControlIndex()) {
                    // under group control
                    well_node.setIndividualControl(false);
                } else {
                    // individual control
                    well_node.setIndividualControl(true);
                }
            }

            if (wellCollection().requireWellPotentials()) {

                // calculate the well potentials
                std::vector<double> well_potentials;
                computeWellPotentials(well_potentials);

                // update/setup guide rates for each well based on the well_potentials
                // TODO: this is one of two places that still need Wells struct. In this function, only the well names
                // well types are used, probably the order of the wells to locate the correct values in well_potentials.
                wellCollection().setGuideRatesWithPotentials(wells(), phase_usage_, well_potentials);
            }

            applyVREPGroupControl();

            if (!wellCollection().groupControlApplied()) {
                wellCollection().applyGroupControls();
            } else {
                wellCollection().updateWellTargets(well_state_.wellRates());
            }
        }
    }

    template<typename TypeTag>
    const WellCollection&
    BlackoilWellModel<TypeTag>::
    wellCollection() const
    {
        return wells_manager_->wellCollection();
    }

    template<typename TypeTag>
    WellCollection&
    BlackoilWellModel<TypeTag>::
    wellCollection()
    {
        return wells_manager_->wellCollection();
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
    calculateEfficiencyFactors()
    {
        if ( !localWellsActive() ) {
            return;
        }

        for (auto& well : well_container_) {
            const std::string& well_name = well->name();
            const WellNode& well_node = wellCollection().findWellNode(well_name);

            const double well_efficiency_factor = well_node.getAccumulativeEfficiencyFactor();

            well->setWellEfficiencyFactor(well_efficiency_factor);
        }
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    computeWellVoidageRates(std::vector<double>& well_voidage_rates,
                            std::vector<double>& voidage_conversion_coeffs) const
    {
        if ( !localWellsActive() ) {
            return;
        }
        // TODO: for now, we store the voidage rates for all the production wells.
        // For injection wells, the rates are stored as zero.
        // We only store the conversion coefficients for all the injection wells.
        // Later, more delicate model will be implemented here.
        // And for the moment, group control can only work for serial running.
        const int nw = numWells();

        const int np = numPhases();

        // we calculate the voidage rate for each well, that means the sum of all the phases.
        well_voidage_rates.resize(nw, 0);
        // store the conversion coefficients, while only for the use of injection wells.
        voidage_conversion_coeffs.resize(nw * np, 1.0);

        std::vector<double> well_rates(np, 0.0);
        std::vector<double> convert_coeff(np, 1.0);

         for (auto& well : well_container_) {
            const bool is_producer = well->wellType() == PRODUCER;
            const int well_cell_top =well->cells()[0];
            const int w = well->indexOfWell();
            const int pvtRegionIdx = pvt_region_idx_[well_cell_top];

            // not sure necessary to change all the value to be positive
            if (is_producer) {
                std::transform(well_state_.wellRates().begin() + np * w,
                               well_state_.wellRates().begin() + np * (w + 1),
                               well_rates.begin(), std::negate<double>());

                // the average hydrocarbon conditions of the whole field will be used
                const int fipreg = 0; // Not considering FIP for the moment.

                rateConverter_->calcCoeff(fipreg, pvtRegionIdx, convert_coeff);
                well_voidage_rates[w] = std::inner_product(well_rates.begin(), well_rates.end(),
                                                           convert_coeff.begin(), 0.0);
            } else {
                // TODO: Not sure whether will encounter situation with all zero rates
                // and whether it will cause problem here.
                std::copy(well_state_.wellRates().begin() + np * w,
                          well_state_.wellRates().begin() + np * (w + 1),
                          well_rates.begin());
                // the average hydrocarbon conditions of the whole field will be used
                const int fipreg = 0; // Not considering FIP for the moment.
                rateConverter_->calcCoeff(fipreg, pvtRegionIdx, convert_coeff);
                std::copy(convert_coeff.begin(), convert_coeff.end(),
                          voidage_conversion_coeffs.begin() + np * w);
            }
        }
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    applyVREPGroupControl()
    {
        if ( wellCollection().havingVREPGroups() ) {
            std::vector<double> well_voidage_rates;
            std::vector<double> voidage_conversion_coeffs;
            computeWellVoidageRates(well_voidage_rates, voidage_conversion_coeffs);
            wellCollection().applyVREPGroupControls(well_voidage_rates, voidage_conversion_coeffs);

            // for the wells under group control, update the control index for the well_state_ and well_controls
            for (const WellNode* well_node : wellCollection().getLeafNodes()) {
                if (well_node->isInjector() && !well_node->individualControl()) {
                    const int well_index = well_node->selfIndex();
                    well_state_.currentControls()[well_index] = well_node->groupControlIndex();

                    WellControls* wc = well_container_[well_index]->wellControls();
                    well_controls_set_current(wc, well_node->groupControlIndex());
                }
            }
        }
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    updateGroupControls(Opm::DeferredLogger& deferred_logger)
    {

        if (wellCollection().groupControlActive()) {
           for (auto& well : well_container_) {
                // update whether well is under group control
                // get well node in the well collection
                WellNode& well_node = wellCollection().findWellNode(well->name());

                // update whehter the well is under group control or individual control
                const int current = well_state_.currentControls()[well->indexOfWell()];
                if (well_node.groupControlIndex() >= 0 && current == well_node.groupControlIndex()) {
                    // under group control
                    well_node.setIndividualControl(false);
                } else {
                    // individual control
                    well_node.setIndividualControl(true);
                }
            }

            applyVREPGroupControl();
            // upate the well targets following group controls
            // it will not change the control mode, only update the targets
            wellCollection().updateWellTargets(well_state_.wellRates());

            // TODO: we should only do the well is involved in the update group targets
            for (auto& well : well_container_) {
                well->updateWellStateWithTarget(ebosSimulator_, well_state_, deferred_logger);
                well->updatePrimaryVariables(well_state_);
            }
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
    computeRepRadiusPerfLength(const Grid& grid)
    {
        for (const auto& well : well_container_) {
            well->computeRepRadiusPerfLength(grid, cartesian_to_compressed_);
        }
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    computeAverageFormationFactor(std::vector<double>& B_avg) const
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
    updatePrimaryVariables()
    {
        for (const auto& well : well_container_) {
            well->updatePrimaryVariables(well_state_);
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
        if (numPhases() == 2) {
            return 2;
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


    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    computeRESV(const std::size_t step)
    {

        const std::vector<int>& resv_wells = SimFIBODetails::resvWells(wells());

        int global_number_resv_wells = resv_wells.size();
        global_number_resv_wells = ebosSimulator_.gridView().comm().sum(global_number_resv_wells);
        if ( global_number_resv_wells > 0 )
        {
            rateConverter_->template defineState<ElementContext>(ebosSimulator_);
        }

        if (! resv_wells.empty()) {
            typedef SimFIBODetails::WellMap WellMap;
            const WellMap& wmap = SimFIBODetails::mapWells(wells_ecl_);

            for (std::vector<int>::const_iterator
                     rp = resv_wells.begin(), e = resv_wells.end();
                 rp != e; ++rp)
            {
                WellControls* ctrl = wells()->ctrls[*rp];
                const bool is_producer = wells()->type[*rp] == PRODUCER;
                const int well_cell_top = wells()->well_cells[wells()->well_connpos[*rp]];
                const int pvtreg = pvt_region_idx_[well_cell_top];

                // RESV control mode, all wells
                {
                    const int rctrl = SimFIBODetails::resv_control(ctrl);
                    const int np = numPhases();
                    std::vector<double> distr (np);

                    if (0 <= rctrl) {
                        const int fipreg = 0; // Hack.  Ignore FIP regions.
                        rateConverter_->calcCoeff(fipreg, pvtreg, distr);

                        if (!is_producer) { // injectors
                            well_controls_assert_number_of_phases(ctrl, np);

                            // original distr contains 0 and 1 to indicate phases under control
                            const double* old_distr = well_controls_get_current_distr(ctrl);

                            for (int p = 0; p < np; ++p) {
                                distr[p] *= old_distr[p];
                            }
                        }

                        well_controls_iset_distr(ctrl, rctrl, & distr[0]);

                        // for the WCONHIST wells, we need to calculate the RESV rates since it can not be specified directly
                        if (is_producer) {
                            const WellMap::const_iterator i = wmap.find(wells()->name[*rp]);

                            if (i == wmap.end()) {
                                OPM_THROW(std::runtime_error, "Failed to find the well " << wells()->name[*rp] << " in wmap.");
                            }
                            const auto* wp = i->second;
                            const WellProductionProperties& production_properties = wp->getProductionProperties(step);
                            // historical phase rates
                            std::vector<double> hrates(np);
                            SimFIBODetails::historyRates(phase_usage_, production_properties, hrates);

                            std::vector<double> hrates_resv(np);
                            rateConverter_->calcReservoirVoidageRates(fipreg, pvtreg, hrates, hrates_resv);

                            const double target = -std::accumulate(hrates_resv.begin(), hrates_resv.end(), 0.0);

                            well_controls_iset_target(ctrl, rctrl, target);
                        }
                    }
                }
            }
        }
    }


    // convert well data from opm-common to well state from opm-core
    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    wellsToState( const data::Wells& wells,
                       PhaseUsage phases,
                       WellStateFullyImplicitBlackoil& state ) {

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
            state.currentControls()[ well_index ] = well.control;
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
        }
    }

} // namespace Opm
