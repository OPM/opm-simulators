


namespace Opm {


    template<typename TypeTag>
    BlackoilWellModel<TypeTag>::
    BlackoilWellModel(Simulator& ebosSimulator,
                      const ModelParameters& param,
                      const bool terminal_output)
        : ebosSimulator_(ebosSimulator)
        , param_(param)
        , terminal_output_(terminal_output)
        , has_solvent_(GET_PROP_VALUE(TypeTag, EnableSolvent))
        , has_polymer_(GET_PROP_VALUE(TypeTag, EnablePolymer))
    {
        const auto& eclState = ebosSimulator_.gridManager().eclState();
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
    }



    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    beginReportStep(const int timeStepIdx)
    {
        const Grid& grid = ebosSimulator_.gridManager().grid();
        const auto& defunct_well_names = ebosSimulator_.gridManager().defunctWellNames();
        const auto& eclState = ebosSimulator_.gridManager().eclState();
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
                                                dynamic_list_econ_limited_,
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
        const auto& gridView = ebosSimulator_.gridManager().gridView();
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
        well_state_.init(wells(), cellPressures, previous_well_state_, phase_usage_);

        std::map<std::string, std::vector<int> > perforation_mapping;
        // handling MS well related
        if (param_.use_multisegment_well_) { // if we use MultisegmentWell model
            for (const auto& well : wells_ecl_) {
                // there is one active well is MS well
                if (well->isMultiSegment(timeStepIdx) && well->getStatus(timeStepIdx) != WellCommon::SHUT) {
                    const int* global_cell = Opm::UgGridHelpers::globalCell(grid);

                    std::map<int,int> cartesian_to_compressed;
                    setupCompressedToCartesian(global_cell, number_of_cells_,
                                               cartesian_to_compressed);

                    const int* cart_dims = Opm::UgGridHelpers::cartDims(grid);

                    well_state_.initWellStateMSWell(wells(), wells_ecl_, timeStepIdx, phase_usage_,
                                                    previous_well_state_, cartesian_to_compressed, cart_dims,
                                                    perforation_mapping);
                    break;
                }
            }
        }

        // Compute reservoir volumes for RESV controls.
        rateConverter_.reset(new RateConverterType (phase_usage_,
                                         std::vector<int>(number_of_cells_, 0)));
        computeRESV(timeStepIdx);

        // create the well container
        well_container_ = createWellContainer(timeStepIdx, perforation_mapping);

        // do the initialization for all the wells
        // TODO: to see whether we can postpone of the intialization of the well containers to
        // optimize the usage of the following several member variables
        for (auto& well : well_container_) {
            well->init(&phase_usage_, depth_, gravity_, number_of_cells_);
        }

        // calculate the efficiency factors for each well
        calculateEfficiencyFactors();

        if (has_polymer_)
        {
            if (PolymerModule::hasPlyshlog()) {
                computeRepRadiusPerfLength(grid);
            }
        }

        // compute VFP properties
        vfp_properties_.reset (new VFPProperties (
                                   eclState.getTableManager().getVFPInjTables(),
                                   eclState.getTableManager().getVFPProdTables()) );

        // update the previous well state. This is used to restart failed steps.
        previous_well_state_ = well_state_;


    }


    // called at the beginning of a time step
    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    beginTimeStep() {
        well_state_ = previous_well_state_;

        if (wellCollection().havingVREPGroups() ) {
            rateConverter_->template defineState<ElementContext>(ebosSimulator_);

        }
    }

    // only use this for restart.
    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    setRestartWellState(const WellState& well_state) { previous_well_state_ = well_state; }

    // called at the end of a report step
    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    endReportStep() {
        // update the list contanining information of closed wells
        // and connections due to economical limits
        // Used by the wellManager
        updateListEconLimited(dynamic_list_econ_limited_);
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
    timeStepSucceeded() {
        previous_well_state_ = well_state_;
    }

    template<typename TypeTag>
    std::vector<typename BlackoilWellModel<TypeTag>::WellInterfacePtr >
    BlackoilWellModel<TypeTag>::
    createWellContainer(const int time_step,
                        const std::map<std::string, std::vector<int> >& perforation_mapping) const
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

                // Use the pvtRegionIdx from the top cell
                const int well_cell_top = wells()->well_cells[wells()->well_connpos[w]];
                const int pvtreg = pvt_region_idx_[well_cell_top];

                if ( !well_ecl->isMultiSegment(time_step) || !param_.use_multisegment_well_) {
                    well_container.emplace_back(new StandardWell<TypeTag>(well_ecl, time_step, wells(),
                                                param_, *rateConverter_, pvtreg, numComponents() ) );
                } else {
                    const auto it = perforation_mapping.find(well_name);
                    if (it == perforation_mapping.end()) { // not found
                        OPM_THROW(std::runtime_error, "Could not find well " + well_name + " in perforation_mapping");
                    }
                    const std::vector<int>& perforation_mapping_well = it->second;
                    well_container.emplace_back(new MultisegmentWell<TypeTag>(well_ecl, time_step, wells(),
                                                param_, *rateConverter_, pvtreg, numComponents(), perforation_mapping_well) );
                }
            }
        }
        return well_container;
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

        updatePerforationIntensiveQuantities();

        if (iterationIdx == 0) {
            prepareTimeStep();
        }

        updateWellControls();
        // Set the well primary variables based on the value of well solutions
        initPrimaryVariablesEvaluation();

        if (iterationIdx == 0) {
            calculateExplicitQuantities();
        }

        if (param_.solve_welleq_initially_ && iterationIdx == 0) {
            // solve the well equations as a pre-processing step
            last_report_ = solveWellEq(dt);
        }
        assembleWellEq(dt, false);

        last_report_.converged = true;
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    assembleWellEq(const double dt,
                   bool only_wells)
    {
        for (int w = 0; w < numWells(); ++w) {
            well_container_[w]->assembleWellEq(ebosSimulator_, dt, well_state_, only_wells);
        }
    }





    // applying the well residual to reservoir residuals
    // r = r - duneC_^T * invDuneD_ * resWell_
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
        const int        nw   = numWells();

        assert(nw == int(well_container_.size()) );

        for (int w = 0; w < nw; ++w) {
            WellControls* wc = well_container_[w]->wellControls();
            well_controls_set_current( wc, well_state_.currentControls()[w]);
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
    solveWellEq(const double dt)
    {
        const int nw = numWells();
        WellState well_state0 = well_state_;

        const int numComp = numComponents();
        std::vector< Scalar > B_avg( numComp, Scalar() );
        computeAverageFormationFactor(B_avg);

        const int max_iter = param_.max_welleq_iter_;

        int it  = 0;
        bool converged;
        do {
            assembleWellEq(dt, true);

            converged = getWellConvergence(B_avg);

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
            if ( terminal_output_ ) {
                OpmLog::debug("Well equation solution gets converged with " + std::to_string(it) + " iterations");
            }
        } else {
            if ( terminal_output_ ) {
                OpmLog::debug("Well equation solution failed in getting converged with " + std::to_string(it) + " iterations");
            }

            well_state_ = well_state0;
            updatePrimaryVariables();
            // also recover the old well controls
            for (int w = 0; w < nw; ++w) {
                WellControls* wc = well_container_[w]->wellControls();
                well_controls_set_current(wc, well_state_.currentControls()[w]);
            }
        }

        SimulatorReport report;
        report.converged = converged;
        report.total_well_iterations = it;
        return report;
    }





    template<typename TypeTag>
    bool
    BlackoilWellModel<TypeTag>::
    getWellConvergence(const std::vector<Scalar>& B_avg) const
    {
        ConvergenceReport report;

        for (const auto& well : well_container_) {
            report += well->getWellConvergence(B_avg);
        }

        // checking NaN residuals
        {
            bool nan_residual_found = report.nan_residual_found;
            const auto& grid = ebosSimulator_.gridManager().grid();
            int value = nan_residual_found ? 1 : 0;

            nan_residual_found = grid.comm().max(value);

            if (nan_residual_found) {
                for (const auto& well : report.nan_residual_wells) {
                    OpmLog::debug("NaN residual found with phase " + well.phase_name + " for well " + well.well_name);
                }
                OPM_THROW(Opm::NumericalProblem, "NaN residual found!");
            }
        }

        // checking too large residuals
        {
            bool too_large_residual_found = report.too_large_residual_found;
            const auto& grid = ebosSimulator_.gridManager().grid();
            int value = too_large_residual_found ? 1 : 0;

            too_large_residual_found = grid.comm().max(value);
            if (too_large_residual_found) {
                for (const auto& well : report.too_large_residual_wells) {
                    OpmLog::debug("Too large residual found with phase " + well.phase_name + " fow well " + well.well_name);
                }
                OPM_THROW(Opm::NumericalProblem, "Too large residual found!");
            }
        }

        // checking convergence
        bool converged_well = report.converged;
        {
            const auto& grid = ebosSimulator_.gridManager().grid();
            int value = converged_well ? 1 : 0;

            converged_well = grid.comm().min(value);
        }

        return converged_well;
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    calculateExplicitQuantities() const
    {
         for (auto& well : well_container_) {
             well->calculateExplicitQuantities(ebosSimulator_, well_state_);
         }
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    updateWellControls()
    {
        // Even if there no wells active locally, we cannot
        // return as the Destructor of the WellSwitchingLogger
        // uses global communication. For no well active globally
        // we simply return.
        if( !wellsActive() ) return ;

#if HAVE_OPENMP
#endif // HAVE_OPENMP
        wellhelpers::WellSwitchingLogger logger;

        for (const auto& well : well_container_) {
            well->updateWellControl(well_state_, logger);
        }

        updateGroupControls();
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    updateListEconLimited(DynamicListEconLimited& list_econ_limited) const
    {
        for (const auto& well : well_container_) {
            well->updateListEconLimited(well_state_, list_econ_limited);
        }
    }



    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    computeWellPotentials(std::vector<double>& well_potentials)
    {
        // number of wells and phases
        const int nw = numWells();
        const int np = numPhases();
        well_potentials.resize(nw * np, 0.0);

        for (int w = 0; w < nw; ++w) {
            std::vector<double> potentials;
            well_container_[w]->computeWellPotentials(ebosSimulator_, well_state_, potentials);

            // putting the sucessfully calculated potentials to the well_potentials
            for (int p = 0; p < np; ++p) {
                well_potentials[w * np + p] = std::abs(potentials[p]);
            }
        } // end of for (int w = 0; w < nw; ++w)
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    prepareTimeStep()
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

        // since the controls are all updated, we should update well_state accordingly
        for (int w = 0; w < numWells(); ++w) {
            WellControls* wc = well_container_[w]->wellControls();
            const int control = well_controls_get_current(wc);
            well_state_.currentControls()[w] = control;
            // TODO: for VFP control, the perf_densities are still zero here, investigate better
            // way to handle it later.
            well_container_[w]->updateWellStateWithTarget(well_state_);

            // The wells are not considered to be newly added
            // for next time step
            if (well_state_.isNewWell(w) ) {
                well_state_.setNewWell(w, false);
            }
        }  // end of for (int w = 0; w < nw; ++w)

    }








    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    prepareGroupControl()
    {
        // group control related processing
        if (wellCollection().groupControlActive()) {
            for (int w = 0; w < numWells(); ++w) {
                WellControls* wc = well_container_[w]->wellControls();
                WellNode& well_node = wellCollection().findWellNode(well_container_[w]->name());

                // handling the situation that wells do not have a valid control
                // it happens the well specified with GRUP and restarting due to non-convergencing
                // putting the well under group control for this situation
                int ctrl_index = well_controls_get_current(wc);

                const int group_control_index = well_node.groupControlIndex();
                if (group_control_index >= 0 && ctrl_index < 0) {
                    // put well under group control
                    well_controls_set_current(wc, group_control_index);
                    well_state_.currentControls()[w] = group_control_index;
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

        const int nw = numWells();

        for (int w = 0; w < nw; ++w) {
            const std::string well_name = well_container_[w]->name();
            const WellNode& well_node = wellCollection().findWellNode(well_name);

            const double well_efficiency_factor = well_node.getAccumulativeEfficiencyFactor();

            well_container_[w]->setWellEfficiencyFactor(well_efficiency_factor);
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

        for (int w = 0; w < nw; ++w) {
            const bool is_producer = well_container_[w]->wellType() == PRODUCER;
            const int well_cell_top = well_container_[w]->cells()[0];
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
    updateGroupControls()
    {

        if (wellCollection().groupControlActive()) {
            for (int w = 0; w < numWells(); ++w) {
                // update whether well is under group control
                // get well node in the well collection
                WellNode& well_node = wellCollection().findWellNode(well_container_[w]->name());

                // update whehter the well is under group control or individual control
                const int current = well_state_.currentControls()[w];
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

            for (int w = 0; w < numWells(); ++w) {
                well_container_[w]->updateWellStateWithTarget(well_state_);
            }
        }
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    setupCompressedToCartesian(const int* global_cell, int number_of_cells, std::map<int,int>& cartesian_to_compressed ) const
    {
        if (global_cell) {
            for (int i = 0; i < number_of_cells; ++i) {
                cartesian_to_compressed.insert(std::make_pair(global_cell[i], i));
            }
        }
        else {
            for (int i = 0; i < number_of_cells; ++i) {
                cartesian_to_compressed.insert(std::make_pair(i, i));
            }
        }

    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    computeRepRadiusPerfLength(const Grid& grid)
    {
        // TODO, the function does not work for parallel running
        // to be fixed later.
        const int* global_cell = Opm::UgGridHelpers::globalCell(grid);

        std::map<int,int> cartesian_to_compressed;
        setupCompressedToCartesian(global_cell, number_of_cells_,
                                    cartesian_to_compressed);

        for (const auto& well : well_container_) {
            well->computeRepRadiusPerfLength(grid, cartesian_to_compressed);
        }
    }





    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    computeAverageFormationFactor(std::vector<double>& B_avg) const
    {
        const auto& grid = ebosSimulator_.gridManager().grid();
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
        const auto& grid = ebosSimulator_.gridManager().grid();
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
        return wells()->number_of_wells;
    }

    template<typename TypeTag>
    int
    BlackoilWellModel<TypeTag>:: numPhases() const
    {
        return wells()->number_of_phases;
    }

    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::extractLegacyDepth_()
    {
        const auto& grid = ebosSimulator_.gridManager().grid();
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
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
        }
    }


    template<typename TypeTag>
    void
    BlackoilWellModel<TypeTag>::
    computeRESV(const std::size_t step)
    {
        typedef SimFIBODetails::WellMap WellMap;

        const WellMap& wmap = SimFIBODetails::mapWells(wells_ecl_);

        const std::vector<int>& resv_wells = SimFIBODetails::resvWells(wells(), step, wmap);

        int global_number_resv_wells = resv_wells.size();
        global_number_resv_wells = ebosSimulator_.gridView().comm().sum(global_number_resv_wells);
        if ( global_number_resv_wells > 0 )
        {
            rateConverter_->template defineState<ElementContext>(ebosSimulator_);
        }

        if (! resv_wells.empty()) {
            const PhaseUsage&                    pu = phase_usage_;
            const std::vector<double>::size_type np = pu.num_phases;

            std::vector<double> distr (np);
            std::vector<double> hrates(np);

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

                    if (0 <= rctrl) {
                        const int fipreg = 0; // Hack.  Ignore FIP regions.
                        rateConverter_->calcCoeff(fipreg, pvtreg, distr);

                        if (!is_producer) { // injectors
                            well_controls_assert_number_of_phases(ctrl, np);

                            // original distr contains 0 and 1 to indicate phases under control
                            const double* old_distr = well_controls_get_current_distr(ctrl);

                            for (size_t p = 0; p < np; ++p) {
                                distr[p] *= old_distr[p];
                            }
                        }

                        well_controls_iset_distr(ctrl, rctrl, & distr[0]);
                    }
                }

                // RESV control, WCONHIST wells.  A bit of duplicate
                // work, regrettably.
                if (is_producer && wells()->name[*rp] != 0) {
                    WellMap::const_iterator i = wmap.find(wells()->name[*rp]);

                    if (i != wmap.end()) {
                        const auto* wp = i->second;

                        const WellProductionProperties& p =
                            wp->getProductionProperties(step);

                        if (! p.predictionMode) {
                            // History matching (WCONHIST/RESV)
                            SimFIBODetails::historyRates(pu, p, hrates);

                            const int fipreg = 0; // Hack.  Ignore FIP regions.
                            rateConverter_->calcCoeff(fipreg, pvtreg, distr);

                            // WCONHIST/RESV target is sum of all
                            // observed phase rates translated to
                            // reservoir conditions.  Recall sign
                            // convention: Negative for producers.
                            const double target =
                                - std::inner_product(distr.begin(), distr.end(),
                                                     hrates.begin(), 0.0);

                            well_controls_clear(ctrl);
                            well_controls_assert_number_of_phases(ctrl, int(np));

                            static const double invalid_alq = -std::numeric_limits<double>::max();
                            static const int invalid_vfp = -std::numeric_limits<int>::max();

                            const int ok_resv =
                                well_controls_add_new(RESERVOIR_RATE, target,
                                                      invalid_alq, invalid_vfp,
                                                      & distr[0], ctrl);

                            // For WCONHIST the BHP limit is set to 1 atm.
                            // or a value specified using WELTARG
                            double bhp_limit = (p.BHPLimit > 0) ? p.BHPLimit : unit::convert::from(1.0, unit::atm);
                            const int ok_bhp =
                                well_controls_add_new(BHP, bhp_limit,
                                                      invalid_alq, invalid_vfp,
                                                      NULL, ctrl);

                            if (ok_resv != 0 && ok_bhp != 0) {
                                well_state_.currentControls()[*rp] = 0;
                                well_controls_set_current(ctrl, 0);
                            }
                        }
                    }
                }
            }
        }

        if( wells() )
        {
            for (int w = 0, nw = numWells(); w < nw; ++w) {
                WellControls* ctrl = wells()->ctrls[w];
                const bool is_producer = wells()->type[w] == PRODUCER;
                if (!is_producer && wells()->name[w] != 0) {
                    WellMap::const_iterator i = wmap.find(wells()->name[w]);
                    if (i != wmap.end()) {
                        const auto* wp = i->second;
                        const WellInjectionProperties& injector = wp->getInjectionProperties(step);
                        if (!injector.predictionMode) {
                            //History matching WCONINJEH
                            static const double invalid_alq = -std::numeric_limits<double>::max();
                            static const int invalid_vfp = -std::numeric_limits<int>::max();
                            // For WCONINJEH the BHP limit is set to a large number
                            // or a value specified using WELTARG
                            double bhp_limit = (injector.BHPLimit > 0) ? injector.BHPLimit : std::numeric_limits<double>::max();
                            const int ok_bhp =
                                well_controls_add_new(BHP, bhp_limit,
                                                      invalid_alq, invalid_vfp,
                                                      NULL, ctrl);
                            if (!ok_bhp) {
                                OPM_THROW(std::runtime_error, "Failed to add well control.");
                            }
                        }
                    }
                }
            }
        }
    }

} // namespace Opm
