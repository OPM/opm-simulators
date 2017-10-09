


namespace Opm {

    template<typename TypeTag>
    StandardWellsDense<TypeTag>::
    StandardWellsDense(const Wells* wells_arg,
                       WellCollection* well_collection,
                       const std::vector< const Well* >& wells_ecl,
                       const ModelParameters& param,
                       const RateConverterType& rate_converter,
                       const bool terminal_output,
                       const int current_timeIdx,
                       std::vector<int>& pvt_region_idx)
       : wells_active_(wells_arg!=nullptr)
       , wells_(wells_arg)
       , wells_ecl_(wells_ecl)
       , number_of_wells_(wells_arg ? (wells_arg->number_of_wells) : 0)
       , number_of_phases_(wells_arg ? (wells_arg->number_of_phases) : 0) // TODO: not sure if it is proper for this way
       , well_container_(createWellContainer(wells_arg, wells_ecl, current_timeIdx) )
       , well_collection_(well_collection)
       , param_(param)
       , terminal_output_(terminal_output)
       , has_solvent_(GET_PROP_VALUE(TypeTag, EnableSolvent))
       , has_polymer_(GET_PROP_VALUE(TypeTag, EnablePolymer))
       , current_timeIdx_(current_timeIdx)
       , rate_converter_(rate_converter)
       , pvt_region_idx_(pvt_region_idx)
    {
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    init(const PhaseUsage phase_usage_arg,
         const std::vector<bool>& active_arg,
         const double gravity_arg,
         const std::vector<double>& depth_arg,
         long int global_nc,
         const Grid& grid)
    {
        // has to be set always for the convergence check!
        global_nc_   = global_nc;

        phase_usage_ = phase_usage_arg;
        active_ = active_arg;

        if ( ! localWellsActive() ) {
            return;
        }

        calculateEfficiencyFactors();

#ifndef NDEBUG
        const auto& pu = phase_usage_;
        const int np = pu.num_phases;

        // assumes the gas fractions are stored after water fractions
        // WellVariablePositions needs to be changed for 2p runs
        assert (np == 3 || (np == 2 && !pu.phase_used[Gas]) );
#endif

        if (has_polymer_)
        {
            if (PolymerModule::hasPlyshlog()) {
                computeRepRadiusPerfLength(grid);
            }
        }

        number_of_cells_ = Opm::UgGridHelpers::numCells(grid);
        // do the initialization for all the wells
        // TODO: to see whether we can postpone of the intialization of the well containers to
        // optimize the usage of the following several member variables
        for (auto& well : well_container_) {
            well->init(&phase_usage_, &active_, depth_arg, gravity_arg, number_of_cells_);
        }
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    setVFPProperties(const VFPProperties*  vfp_properties_arg)
    {
        for (auto& well : well_container_) {
            well->setVFPProperties(vfp_properties_arg);
        }
    }






    template<typename TypeTag>
    int
    StandardWellsDense<TypeTag>::
    numWells() const
    {
        return number_of_wells_;
    }





    template<typename TypeTag>
    std::vector<typename StandardWellsDense<TypeTag>::WellInterfacePtr >
    StandardWellsDense<TypeTag>::
    createWellContainer(const Wells* wells,
                        const std::vector< const Well* >& wells_ecl,
                        const int time_step)
    {
        std::vector<WellInterfacePtr> well_container;

        const int nw = wells ? (wells->number_of_wells) : 0;

        if (nw > 0) {
            well_container.reserve(nw);

            // With the following way, it will have the same order with wells struct
            // Hopefully, it can generate the same residual history with master branch
            for (int w = 0; w < nw; ++w) {
                const std::string well_name = std::string(wells->name[w]);

                // finding the location of the well in wells_ecl
                const int nw_wells_ecl = wells_ecl.size();
                int index_well = 0;
                for (; index_well < nw_wells_ecl; ++index_well) {
                    if (well_name == wells_ecl[index_well]->name()) {
                        break;
                    }
                }

                // It should be able to find in wells_ecl.
                if (index_well == nw_wells_ecl) {
                    OPM_THROW(std::logic_error, "Could not find well " << well_name << " in wells_ecl ");
                }

                const Well* well_ecl = wells_ecl[index_well];
                // TODO: stopping throwing when encoutnering MS wells for now.
                /* if (well_ecl->isMultiSegment(time_step)) {
                    OPM_THROW(Opm::NumericalProblem, "Not handling Multisegment Wells for now");
                } */

                // Basically, we are handling all the wells as StandardWell for the moment
                well_container.emplace_back(new StandardWell<TypeTag>(well_ecl, time_step, wells) );
            }
        }
        return well_container;
    }





    template<typename TypeTag>
    SimulatorReport
    StandardWellsDense<TypeTag>::
    assemble(Simulator& ebosSimulator,
             const int iterationIdx,
             const double dt,
             WellState& well_state)
    {

        if (iterationIdx == 0) {
            prepareTimeStep(ebosSimulator, well_state);
        }

        SimulatorReport report;
        if ( ! wellsActive() ) {
            return report;
        }

        updateWellControls(well_state);
        // Set the well primary variables based on the value of well solutions
        initPrimaryVariablesEvaluation();

        if (iterationIdx == 0) {
            computeWellConnectionPressures(ebosSimulator, well_state);
            computeAccumWells();
        }

        if (param_.solve_welleq_initially_ && iterationIdx == 0) {
            // solve the well equations as a pre-processing step
            report = solveWellEq(ebosSimulator, dt, well_state);
        }
        assembleWellEq(ebosSimulator, dt, well_state, false);

        report.converged = true;
        return report;
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    assembleWellEq(Simulator& ebosSimulator,
                   const double dt,
                   WellState& well_state,
                   bool only_wells) const
    {
        for (int w = 0; w < number_of_wells_; ++w) {
            well_container_[w]->assembleWellEq(ebosSimulator, dt, well_state, only_wells);
        }
    }





    // applying the well residual to reservoir residuals
    // r = r - duneC_^T * invDuneD_ * resWell_
    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
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
    StandardWellsDense<TypeTag>::
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
    StandardWellsDense<TypeTag>::
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
    StandardWellsDense<TypeTag>::
    recoverWellSolutionAndUpdateWellState(const BVector& x, WellState& well_state) const
    {
        for (auto& well : well_container_) {
            well->recoverWellSolutionAndUpdateWellState(x, param_, well_state);
        }
    }





    template<typename TypeTag>
    int
    StandardWellsDense<TypeTag>::
    flowPhaseToEbosPhaseIdx( const int phaseIdx ) const
    {
        const auto& pu = phase_usage_;
        if (active_[Water] && pu.phase_pos[Water] == phaseIdx)
            return FluidSystem::waterPhaseIdx;
        if (active_[Oil] && pu.phase_pos[Oil] == phaseIdx)
            return FluidSystem::oilPhaseIdx;
        if (active_[Gas] && pu.phase_pos[Gas] == phaseIdx)
            return FluidSystem::gasPhaseIdx;

        assert(phaseIdx < 3);
        // for other phases return the index
        return phaseIdx;
    }





    template<typename TypeTag>
    int
    StandardWellsDense<TypeTag>::
    numPhases() const
    {
        return number_of_phases_;
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    resetWellControlFromState(const WellState& xw) const
    {
        const int        nw   = wells_->number_of_wells;
        for (int w = 0; w < nw; ++w) {
            WellControls* wc = wells_->ctrls[w];
            well_controls_set_current( wc, xw.currentControls()[w]);
        }
    }





    template<typename TypeTag>
    bool
    StandardWellsDense<TypeTag>::
    wellsActive() const
    {
        return wells_active_;
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    setWellsActive(const bool wells_active)
    {
        wells_active_ = wells_active;
    }





    template<typename TypeTag>
    bool
    StandardWellsDense<TypeTag>::
    localWellsActive() const
    {
        return number_of_wells_ > 0;
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    initPrimaryVariablesEvaluation() const
    {
        for (auto& well : well_container_) {
            well->initPrimaryVariablesEvaluation();
        }
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    computeAccumWells() const
    {
        for (auto& well : well_container_) {
            well->computeAccumWell();
        }
    }





    template<typename TypeTag>
    SimulatorReport
    StandardWellsDense<TypeTag>::
    solveWellEq(Simulator& ebosSimulator,
                const double dt,
                WellState& well_state) const
    {
        const int nw = number_of_wells_;
        WellState well_state0 = well_state;

        const int numComp = numComponents();
        std::vector< Scalar > B_avg( numComp, Scalar() );
        computeAverageFormationFactor(ebosSimulator, B_avg);

        int it  = 0;
        bool converged;
        do {
            assembleWellEq(ebosSimulator, dt, well_state, true);

            converged = getWellConvergence(ebosSimulator, B_avg);

            // checking whether the group targets are converged
            if (wellCollection()->groupControlActive()) {
                converged = converged && wellCollection()->groupTargetConverged(well_state.wellRates());
            }

            if (converged) {
                break;
            }

            ++it;
            if( localWellsActive() )
            {
                for (auto& well : well_container_) {
                    well->solveEqAndUpdateWellState(param_, well_state);
                }
            }
            // updateWellControls uses communication
            // Therefore the following is executed if there
            // are active wells anywhere in the global domain.
            if( wellsActive() )
            {
                updateWellControls(well_state);
                initPrimaryVariablesEvaluation();
            }
        } while (it < 15);

        if (converged) {
            if ( terminal_output_ ) {
                OpmLog::debug("Well equation solution gets converged with " + std::to_string(it) + " iterations");
            }
        } else {
            if ( terminal_output_ ) {
                OpmLog::debug("Well equation solution failed in getting converged with " + std::to_string(it) + " iterations");
            }

            well_state = well_state0;
            updatePrimaryVariables(well_state);
            // also recover the old well controls
            for (int w = 0; w < nw; ++w) {
                WellControls* wc = well_container_[w]->wellControls();
                well_controls_set_current(wc, well_state.currentControls()[w]);
            }
        }

        SimulatorReport report;
        report.converged = converged;
        report.total_well_iterations = it;
        return report;
    }





    template<typename TypeTag>
    bool
    StandardWellsDense<TypeTag>::
    getWellConvergence(Simulator& ebosSimulator,
                       const std::vector<Scalar>& B_avg) const
    {
        ConvergenceReport report;

        for (const auto& well : well_container_) {
            report += well->getWellConvergence(ebosSimulator, B_avg, param_);
        }

        // checking NaN residuals
        {
            bool nan_residual_found = report.nan_residual_found;
            const auto& grid = ebosSimulator.gridManager().grid();
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
            const auto& grid = ebosSimulator.gridManager().grid();
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
            const auto& grid = ebosSimulator.gridManager().grid();
            int value = converged_well ? 1 : 0;

            converged_well = grid.comm().min(value);
        }

        return converged_well;
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    computeWellConnectionPressures(const Simulator& ebosSimulator,
                                   const WellState& xw) const
    {
         if( ! localWellsActive() ) return ;

         for (auto& well : well_container_) {
             well->computeWellConnectionPressures(ebosSimulator, xw);
         }
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    updateWellControls(WellState& xw) const
    {
        // Even if there no wells active locally, we cannot
        // return as the Destructor of the WellSwitchingLogger
        // uses global communication. For no well active globally
        // we simply return.
        if( !wellsActive() ) return ;

        wellhelpers::WellSwitchingLogger logger;

        for (const auto& well : well_container_) {
            well->updateWellControl(xw, logger);
        }

        updateGroupControls(xw);
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    updateListEconLimited(const Schedule& schedule,
                          const int current_step,
                          const Wells* wells_struct,
                          const WellState& well_state,
                          DynamicListEconLimited& list_econ_limited) const
    {
        for (const auto& well : well_container_) {
            well->updateListEconLimited(well_state, list_econ_limited);
        }
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    computeWellPotentials(const Simulator& ebosSimulator,
                          const WellState& well_state,
                          std::vector<double>& well_potentials) const
    {
        updatePrimaryVariables(well_state);
        computeWellConnectionPressures(ebosSimulator, well_state);

        // initialize the primary variables in Evaluation, which is used in computePerfRate for computeWellPotentials
        // TODO: for computeWellPotentials, no derivative is required actually
        initPrimaryVariablesEvaluation();

        // number of wells and phases
        const int nw = number_of_wells_;
        const int np = number_of_phases_;
        well_potentials.resize(nw * np, 0.0);

        for (int w = 0; w < nw; ++w) {
            std::vector<double> potentials;
            well_container_[w]->computeWellPotentials(ebosSimulator, well_state, potentials);

            // putting the sucessfully calculated potentials to the well_potentials
            for (int p = 0; p < np; ++p) {
                well_potentials[w * np + p] = std::abs(potentials[p]);
            }
        } // end of for (int w = 0; w < nw; ++w)
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    prepareTimeStep(const Simulator& ebos_simulator,
                    WellState& well_state)
    {
        const int nw = number_of_wells_;
        for (int w = 0; w < nw; ++w) {
            // after restarting, the well_controls can be modified while
            // the well_state still uses the old control index
            // we need to synchronize these two.
            resetWellControlFromState(well_state);

            if (wellCollection()->groupControlActive()) {
                WellControls* wc = well_container_[w]->wellControls();
                WellNode& well_node = well_collection_->findWellNode(well_container_[w]->name());

                // handling the situation that wells do not have a valid control
                // it happens the well specified with GRUP and restarting due to non-convergencing
                // putting the well under group control for this situation
                int ctrl_index = well_controls_get_current(wc);

                const int group_control_index = well_node.groupControlIndex();
                if (group_control_index >= 0 && ctrl_index < 0) {
                    // put well under group control
                    well_controls_set_current(wc, group_control_index);
                    well_state.currentControls()[w] = group_control_index;
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
        }

        if (well_collection_->groupControlActive()) {
            if (well_collection_->requireWellPotentials()) {

                // calculate the well potentials
                std::vector<double> well_potentials;
                computeWellPotentials(ebos_simulator, well_state, well_potentials);

                // update/setup guide rates for each well based on the well_potentials
                // TODO: this is one of two places that still need Wells struct. In this function, only the well names
                // well types are used, probably the order of the wells to locate the correct values in well_potentials.
                well_collection_->setGuideRatesWithPotentials(wells_, phase_usage_, well_potentials);
            }

            applyVREPGroupControl(well_state);

            if (!wellCollection()->groupControlApplied()) {
                wellCollection()->applyGroupControls();
            } else {
                wellCollection()->updateWellTargets(well_state.wellRates());
            }
        }

        // since the controls are all updated, we should update well_state accordingly
        for (int w = 0; w < nw; ++w) {
            WellControls* wc = well_container_[w]->wellControls();
            const int control = well_controls_get_current(wc);
            well_state.currentControls()[w] = control;
            // TODO: for VFP control, the perf_densities are still zero here, investigate better
            // way to handle it later.
            well_container_[w]->updateWellStateWithTarget(control, well_state);

            // The wells are not considered to be newly added
            // for next time step
            if (well_state.isNewWell(w) ) {
                well_state.setNewWell(w, false);
            }
        }  // end of for (int w = 0; w < nw; ++w)
    }





    template<typename TypeTag>
    WellCollection*
    StandardWellsDense<TypeTag>::
    wellCollection() const
    {
        return well_collection_;
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    calculateEfficiencyFactors()
    {
        if ( !localWellsActive() ) {
            return;
        }

        const int nw = number_of_wells_;

        for (int w = 0; w < nw; ++w) {
            const std::string well_name = well_container_[w]->name();
            const WellNode& well_node = wellCollection()->findWellNode(well_name);

            const double well_efficiency_factor = well_node.getAccumulativeEfficiencyFactor();

            well_container_[w]->setWellEfficiencyFactor(well_efficiency_factor);
        }
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    computeWellVoidageRates(const WellState& well_state,
                            std::vector<double>& well_voidage_rates,
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
                std::transform(well_state.wellRates().begin() + np * w,
                               well_state.wellRates().begin() + np * (w + 1),
                               well_rates.begin(), std::negate<double>());

                // the average hydrocarbon conditions of the whole field will be used
                const int fipreg = 0; // Not considering FIP for the moment.

                rate_converter_.calcCoeff(fipreg, pvtRegionIdx, convert_coeff);
                well_voidage_rates[w] = std::inner_product(well_rates.begin(), well_rates.end(),
                                                           convert_coeff.begin(), 0.0);
            } else {
                // TODO: Not sure whether will encounter situation with all zero rates
                // and whether it will cause problem here.
                std::copy(well_state.wellRates().begin() + np * w,
                          well_state.wellRates().begin() + np * (w + 1),
                          well_rates.begin());
                // the average hydrocarbon conditions of the whole field will be used
                const int fipreg = 0; // Not considering FIP for the moment.
                rate_converter_.calcCoeff(fipreg, pvtRegionIdx, convert_coeff);
                std::copy(convert_coeff.begin(), convert_coeff.end(),
                          voidage_conversion_coeffs.begin() + np * w);
            }
        }
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    applyVREPGroupControl(WellState& well_state) const
    {
        if ( wellCollection()->havingVREPGroups() ) {
            std::vector<double> well_voidage_rates;
            std::vector<double> voidage_conversion_coeffs;
            computeWellVoidageRates(well_state, well_voidage_rates, voidage_conversion_coeffs);
            wellCollection()->applyVREPGroupControls(well_voidage_rates, voidage_conversion_coeffs);

            // for the wells under group control, update the control index for the well_state and well_controls
            for (const WellNode* well_node : wellCollection()->getLeafNodes()) {
                if (well_node->isInjector() && !well_node->individualControl()) {
                    const int well_index = well_node->selfIndex();
                    well_state.currentControls()[well_index] = well_node->groupControlIndex();

                    WellControls* wc = well_container_[well_index]->wellControls();
                    well_controls_set_current(wc, well_node->groupControlIndex());
                }
            }
        }
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
    updateGroupControls(WellState& well_state) const
    {

        if (wellCollection()->groupControlActive()) {
            for (int w = 0; w < number_of_wells_; ++w) {
                // update whether well is under group control
                // get well node in the well collection
                WellNode& well_node = well_collection_->findWellNode(well_container_[w]->name());

                // update whehter the well is under group control or individual control
                const int current = well_state.currentControls()[w];
                if (well_node.groupControlIndex() >= 0 && current == well_node.groupControlIndex()) {
                    // under group control
                    well_node.setIndividualControl(false);
                } else {
                    // individual control
                    well_node.setIndividualControl(true);
                }
            }

            applyVREPGroupControl(well_state);
            // upate the well targets following group controls
            // it will not change the control mode, only update the targets
            wellCollection()->updateWellTargets(well_state.wellRates());

            for (int w = 0; w < number_of_wells_; ++w) {
                // TODO: check whether we need current argument in updateWellStateWithTarget
                // maybe there is some circumstances that the current is different from the one
                // in the WellState.
                // while probalby, the current argument can be removed
                const int current = well_state.currentControls()[w];
                well_container_[w]->updateWellStateWithTarget(current, well_state);
            }
        }
    }





    template<typename TypeTag>
    void
    StandardWellsDense<TypeTag>::
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
    StandardWellsDense<TypeTag>::
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
    StandardWellsDense<TypeTag>::
    computeAverageFormationFactor(Simulator& ebosSimulator,
                                  std::vector<double>& B_avg) const
    {
        const int np = numPhases();

        const auto& grid = ebosSimulator.gridManager().grid();
        const auto& gridView = grid.leafGridView();
        ElementContext elemCtx(ebosSimulator);
        const auto& elemEndIt = gridView.template end</*codim=*/0, Dune::Interior_Partition>();

        for (auto elemIt = gridView.template begin</*codim=*/0, Dune::Interior_Partition>();
             elemIt != elemEndIt; ++elemIt)
        {
            elemCtx.updatePrimaryStencil(*elemIt);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

            const auto& intQuants = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& fs = intQuants.fluidState();

            for ( int phaseIdx = 0; phaseIdx < np; ++phaseIdx )
            {
                auto& B  = B_avg[ phaseIdx ];
                const int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(phaseIdx);

                B += 1 / fs.invB(ebosPhaseIdx).value();
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
    StandardWellsDense<TypeTag>::
    updatePrimaryVariables(const WellState& well_state) const
    {
        for (const auto& well : well_container_) {
            well->updatePrimaryVariables(well_state);
        }
    }

} // namespace Opm
