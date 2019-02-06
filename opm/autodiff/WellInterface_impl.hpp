/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2018 IRIS

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


namespace Opm
{


    template<typename TypeTag>
    WellInterface<TypeTag>::
    WellInterface(const Well* well, const int time_step, const Wells* wells,
                  const ModelParameters& param,
                  const RateConverterType& rate_converter,
                  const int pvtRegionIdx,
                  const int num_components)
    : well_ecl_(well)
    , current_step_(time_step)
    , param_(param)
    , rateConverter_(rate_converter)
    , pvtRegionIdx_(pvtRegionIdx)
    , num_components_(num_components)
    {
        if (!well) {
            OPM_THROW(std::invalid_argument, "Null pointer of Well is used to construct WellInterface");
        }

        if (time_step < 0) {
            OPM_THROW(std::invalid_argument, "Negtive time step is used to construct WellInterface");
        }

        if (!wells) {
            OPM_THROW(std::invalid_argument, "Null pointer of Wells is used to construct WellInterface");
        }

        const std::string& well_name = well->name();

        // looking for the location of the well in the wells struct
        int index_well;
        for (index_well = 0; index_well < wells->number_of_wells; ++index_well) {
            if (well_name == std::string(wells->name[index_well])) {
                break;
            }
        }

        // should not enter the constructor if the well does not exist in the wells struct
        // here, just another assertion.
        assert(index_well != wells->number_of_wells);

        index_of_well_ = index_well;
        well_type_ = wells->type[index_well];
        number_of_phases_ = wells->number_of_phases;

        // copying the comp_frac
        {
            comp_frac_.resize(number_of_phases_);
            const int index_begin = index_well * number_of_phases_;
            std::copy(wells->comp_frac + index_begin,
                      wells->comp_frac + index_begin + number_of_phases_, comp_frac_.begin() );
        }

        well_controls_ = wells->ctrls[index_well];

        ref_depth_ = wells->depth_ref[index_well];

        // perforations related
        {
            const int perf_index_begin = wells->well_connpos[index_well];
            const int perf_index_end = wells->well_connpos[index_well + 1];
            number_of_perforations_ = perf_index_end - perf_index_begin;
            first_perf_ = perf_index_begin;

            well_cells_.resize(number_of_perforations_);
            std::copy(wells->well_cells + perf_index_begin,
                      wells->well_cells + perf_index_end,
                      well_cells_.begin() );

            well_index_.resize(number_of_perforations_);
            std::copy(wells->WI + perf_index_begin,
                      wells->WI + perf_index_end,
                      well_index_.begin() );

            saturation_table_number_.resize(number_of_perforations_);
            std::copy(wells->sat_table_id + perf_index_begin,
                      wells->sat_table_id + perf_index_end,
                      saturation_table_number_.begin() );
        }

        well_efficiency_factor_ = 1.0;

        connectionRates_.resize(number_of_perforations_);

        well_productivity_index_logger_counter_ = 0;

    }

    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    updatePerforatedCell(std::vector<bool>& is_cell_perforated)
    {

        for (int perf_idx = 0; perf_idx<number_of_perforations_; ++perf_idx) {
            is_cell_perforated[well_cells_[perf_idx]] = true;
        }
    }



    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    init(const PhaseUsage* phase_usage_arg,
         const std::vector<double>& /* depth_arg */,
         const double gravity_arg,
         const int /* num_cells */)
    {
        phase_usage_ = phase_usage_arg;
        gravity_ = gravity_arg;
    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    setVFPProperties(const VFPProperties<VFPInjProperties,VFPProdProperties>* vfp_properties_arg)
    {
        vfp_properties_ = vfp_properties_arg;
    }





    template<typename TypeTag>
    const std::string&
    WellInterface<TypeTag>::
    name() const
    {
        return well_ecl_->name();
    }





    template<typename TypeTag>
    WellType
    WellInterface<TypeTag>::
    wellType() const
    {
        return well_type_;
    }





    template<typename TypeTag>
    WellControls*
    WellInterface<TypeTag>::
    wellControls() const
    {
        return well_controls_;
    }

    template<typename TypeTag>
    const int
    WellInterface<TypeTag>::
    indexOfWell() const
    {
        return index_of_well_;
    }






    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    getAllowCrossFlow() const
    {
        return well_ecl_->getAllowCrossFlow();
    }




    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    setWellEfficiencyFactor(const double efficiency_factor)
    {
        well_efficiency_factor_ = efficiency_factor;
    }



    template<typename TypeTag>
    const Well*
    WellInterface<TypeTag>::
    wellEcl() const
    {
        return well_ecl_;
    }



    template<typename TypeTag>
    const PhaseUsage&
    WellInterface<TypeTag>::
    phaseUsage() const
    {
        assert(phase_usage_);

        return *phase_usage_;
    }





    template<typename TypeTag>
    int
    WellInterface<TypeTag>::
    flowPhaseToEbosCompIdx( const int phaseIdx ) const
    {
        const auto& pu = phaseUsage();
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) && pu.phase_pos[Water] == phaseIdx)
            return Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && pu.phase_pos[Oil] == phaseIdx)
            return Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && pu.phase_pos[Gas] == phaseIdx)
            return Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);

        // for other phases return the index
        return phaseIdx;
    }

    template<typename TypeTag>
    int
    WellInterface<TypeTag>::
    ebosCompIdxToFlowCompIdx( const unsigned compIdx ) const
    {
        const auto& pu = phaseUsage();
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) && Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx) == compIdx)
            return pu.phase_pos[Water];
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx) == compIdx)
            return pu.phase_pos[Oil];
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx) == compIdx)
            return pu.phase_pos[Gas];

        // for other phases return the index
        return compIdx;
    }




    template<typename TypeTag>
    double
    WellInterface<TypeTag>::
    wsolvent() const
    {
        if (!has_solvent) {
            return 0.0;
        }

        WellInjectionProperties injection = well_ecl_->getInjectionProperties(current_step_);
        if (injection.injectorType == WellInjector::GAS) {
            double solvent_fraction = well_ecl_->getSolventFraction(current_step_);
            return solvent_fraction;
        } else {
            // Not a gas injection well => no solvent.
            return 0.0;
        }
    }





    template<typename TypeTag>
    double
    WellInterface<TypeTag>::
    wpolymer() const
    {
        if (!has_polymer) {
            return 0.0;
        }

        WellInjectionProperties injection = well_ecl_->getInjectionProperties(current_step_);
        WellPolymerProperties polymer = well_ecl_->getPolymerProperties(current_step_);

        if (injection.injectorType == WellInjector::WATER) {
            const double polymer_injection_concentration = polymer.m_polymerConcentration;
            return polymer_injection_concentration;
        } else {
            // Not a water injection well => no polymer.
            return 0.0;
        }
    }





    template<typename TypeTag>
    double
    WellInterface<TypeTag>::
    mostStrictBhpFromBhpLimits() const
    {
        double bhp;

        // initial bhp value, making the value not usable
        switch( well_type_ ) {
        case INJECTOR:
            bhp = std::numeric_limits<double>::max();
            break;
        case PRODUCER:
            bhp = -std::numeric_limits<double>::max();
            break;
        default:
            OPM_THROW(std::logic_error, "Expected PRODUCER or INJECTOR type for well " << name());
        }

        // The number of the well controls/constraints
        const int nwc = well_controls_get_num(well_controls_);

        for (int ctrl_index = 0; ctrl_index < nwc; ++ctrl_index) {
            // finding a BHP constraint
            if (well_controls_iget_type(well_controls_, ctrl_index) == BHP) {
                // get the bhp constraint value, it should always be postive assummingly
                const double bhp_target = well_controls_iget_target(well_controls_, ctrl_index);

                switch(well_type_) {
                case INJECTOR: // using the lower bhp contraint from Injectors
                    if (bhp_target < bhp) {
                        bhp = bhp_target;
                    }
                    break;
                case PRODUCER:
                    if (bhp_target > bhp) {
                        bhp = bhp_target;
                    }
                    break;
                default:
                    OPM_THROW(std::logic_error, "Expected PRODUCER or INJECTOR type for well " << name());
                } // end of switch
            }
        }

        return bhp;
    }





    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    wellHasTHPConstraints() const
    {
        return getTHPControlIndex() >= 0;
    }




    template<typename TypeTag>
    double
    WellInterface<TypeTag>::
    getTHPConstraint() const
    {
        const int thp_control_index = getTHPControlIndex();

        if (thp_control_index < 0) {
            OPM_THROW(std::runtime_error, " there is no THP constraint/limit for well " << name()
                                          << ", while we are requesting it ");
        }

        return well_controls_iget_target(well_controls_, thp_control_index);
    }




    template<typename TypeTag>
    int
    WellInterface<TypeTag>::
    getTHPControlIndex() const
    {
        const int nwc = well_controls_get_num(well_controls_);
        for (int ctrl_index = 0; ctrl_index < nwc; ++ctrl_index) {
            if (well_controls_iget_type(well_controls_, ctrl_index) == THP) {
                return ctrl_index;
            }
        }
        return -1;
    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    updateWellControl(/* const */ Simulator& ebos_simulator,
                      WellState& well_state,
                      Opm::DeferredLogger& deferred_logger) /* const */
    {
        auto cc = Dune::MPIHelper::getCollectiveCommunication();
        const int np = number_of_phases_;
        const int w = index_of_well_;

        const int old_control_index = well_state.currentControls()[w];

        // Find, for each well, if any constraints are broken. If so,
        // switch control to first broken constraint.
        WellControls* wc = well_controls_;

        // Loop over all controls except the current one, and also
        // skip any RESERVOIR_RATE controls, since we cannot
        // handle those.
        const int nwc = well_controls_get_num(wc);
        // the current control index
        int current = well_state.currentControls()[w];
        int ctrl_index = 0;
        for (; ctrl_index < nwc; ++ctrl_index) {
            if (ctrl_index == current) {
                // This is the currently used control, so it is
                // used as an equation. So this is not used as an
                // inequality constraint, and therefore skipped.
                continue;
            }
            if (wellhelpers::constraintBroken(
                    well_state.bhp(), well_state.thp(), well_state.wellRates(),
                    w, np, well_type_, wc, ctrl_index)) {

                // if the well can not work under THP / BHP control, we should not switch to THP / BHP control
                const bool cannot_switch_to_bhp = well_controls_iget_type(wc, ctrl_index) == BHP && !operability_status_.isOperableUnderBHPLimit();
                const bool cannot_switch_to_thp = well_controls_iget_type(wc, ctrl_index) == THP && !operability_status_.isOperableUnderTHPLimit();
                const bool cannot_switch = cannot_switch_to_bhp || cannot_switch_to_thp;
                if ( !cannot_switch ) {

                    // ctrl_index will be the index of the broken constraint after the loop.
                    break;
                } else {
                    // before we figure out to handle it, we give some debug information here
                    if ( well_controls_iget_type(wc, ctrl_index) == BHP && !operability_status_.isOperableUnderBHPLimit() ) {
                        deferred_logger.debug("well " + name() + " breaks the BHP limit, while it is not operable under BHP limit");
                    }

                    if ( well_controls_iget_type(wc, ctrl_index) == THP && !operability_status_.isOperableUnderTHPLimit() ) {
                        deferred_logger.debug("well " + name() + " breaks the THP limit, while it is not operable under THP limit");
                    }
                }
            }
        }

        if (ctrl_index != nwc) {
            // Constraint number ctrl_index was broken, switch to it.
            well_state.currentControls()[w] = ctrl_index;
            current = well_state.currentControls()[w];
            well_controls_set_current( wc, current);
        }

        // the new well control indices after all the related updates,
        const int updated_control_index = well_state.currentControls()[w];

        // checking whether control changed
        if (updated_control_index != old_control_index) {

            auto from = well_controls_iget_type(wc, old_control_index);
            auto to  = well_controls_iget_type(wc, updated_control_index);
            std::ostringstream ss;
            ss << "    Switching control mode for well " << name()
               << " from " << modestring[from]
               << " to " <<  modestring[to];
            if (cc.size() > 1) {
               ss << " on rank " << cc.rank();
            }
            deferred_logger.info(ss.str());
        }

        if (updated_control_index != old_control_index) { //  || well_collection_->groupControlActive()) {
            updateWellStateWithTarget(ebos_simulator, well_state, deferred_logger);
            updatePrimaryVariables(well_state);
        }
    }





    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    underPredictionMode() const
    {
        bool under_prediction_mode = false;

        switch( well_type_ ) {
        case PRODUCER:
            under_prediction_mode = well_ecl_->getProductionProperties(current_step_).predictionMode;
            break;
        case INJECTOR:
            under_prediction_mode = well_ecl_->getInjectionProperties(current_step_).predictionMode;
            break;
        default:
            OPM_THROW(std::logic_error, "Expected PRODUCER or INJECTOR type for well " << name());
        }

        return under_prediction_mode;
    }





    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    checkRateEconLimits(const WellEconProductionLimits& econ_production_limits,
                        const WellState& well_state,
                        Opm::DeferredLogger& deferred_logger) const
    {
        const Opm::PhaseUsage& pu = phaseUsage();
        const int np = number_of_phases_;

        if (econ_production_limits.onMinOilRate()) {
            assert(FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx));
            const double oil_rate = well_state.wellRates()[index_of_well_ * np + pu.phase_pos[ Oil ] ];
            const double min_oil_rate = econ_production_limits.minOilRate();
            if (std::abs(oil_rate) < min_oil_rate) {
                return true;
            }
        }

        if (econ_production_limits.onMinGasRate() ) {
            assert(FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx));
            const double gas_rate = well_state.wellRates()[index_of_well_ * np + pu.phase_pos[ Gas ] ];
            const double min_gas_rate = econ_production_limits.minGasRate();
            if (std::abs(gas_rate) < min_gas_rate) {
                return true;
            }
        }

        if (econ_production_limits.onMinLiquidRate() ) {
            assert(FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx));
            assert(FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx));
            const double oil_rate = well_state.wellRates()[index_of_well_ * np + pu.phase_pos[ Oil ] ];
            const double water_rate = well_state.wellRates()[index_of_well_ * np + pu.phase_pos[ Water ] ];
            const double liquid_rate = oil_rate + water_rate;
            const double min_liquid_rate = econ_production_limits.minLiquidRate();
            if (std::abs(liquid_rate) < min_liquid_rate) {
                return true;
            }
        }

        if (econ_production_limits.onMinReservoirFluidRate()) {
            deferred_logger.warning("NOT_SUPPORTING_MIN_RESERVOIR_FLUID_RATE", "Minimum reservoir fluid production rate limit is not supported yet");
        }

        return false;
    }






    template<typename TypeTag>
    typename WellInterface<TypeTag>::RatioCheckTuple
    WellInterface<TypeTag>::
    checkMaxWaterCutLimit(const WellEconProductionLimits& econ_production_limits,
                          const WellState& well_state) const
    {
        bool water_cut_limit_violated = false;
        int worst_offending_completion = INVALIDCOMPLETION;
        double violation_extent = -1.0;

        const int np = number_of_phases_;
        const Opm::PhaseUsage& pu = phaseUsage();
        const int well_number = index_of_well_;

        assert(FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx));
        assert(FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx));

        const double oil_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Oil ] ];
        const double water_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Water ] ];
        const double liquid_rate = oil_rate + water_rate;
        double water_cut;
        if (std::abs(liquid_rate) != 0.) {
            water_cut = water_rate / liquid_rate;
        } else {
            water_cut = 0.0;
        }

        const double max_water_cut_limit = econ_production_limits.maxWaterCut();
        if (water_cut > max_water_cut_limit) {
            water_cut_limit_violated = true;
        }

        if (water_cut_limit_violated) {
            // need to handle the worst_offending_connection
            const int perf_start = first_perf_;
            const int perf_number = number_of_perforations_;

            std::vector<double> water_cut_perf(perf_number);
            for (int perf = 0; perf < perf_number; ++perf) {
                const int i_perf = perf_start + perf;
                const double oil_perf_rate = well_state.perfPhaseRates()[i_perf * np + pu.phase_pos[ Oil ] ];
                const double water_perf_rate = well_state.perfPhaseRates()[i_perf * np + pu.phase_pos[ Water ] ];
                const double liquid_perf_rate = oil_perf_rate + water_perf_rate;
                if (std::abs(liquid_perf_rate) != 0.) {
                    water_cut_perf[perf] = water_perf_rate / liquid_perf_rate;
                } else {
                    water_cut_perf[perf] = 0.;
                }
            }
            const auto& completions = well_ecl_->getCompletions(current_step_);
            const auto& connections = well_ecl_->getConnections(current_step_);

            int complnumIdx = 0;
            std::vector<double> water_cut_in_completions(completions.size(), 0.0);
            for (const auto& completion : completions) {
                int complnum = completion.first;
                for (int perf = 0; perf < perf_number; ++perf) {
                    if (complnum == connections.get ( perf ).complnum()) {
                        water_cut_in_completions[complnumIdx] +=  water_cut_perf[perf];
                    }
                }
                complnumIdx++;
            }

            double max_water_cut_perf = 0.;
            complnumIdx = 0;
            for (const auto& completion : completions) {
                if (water_cut_in_completions[complnumIdx] > max_water_cut_perf) {
                    worst_offending_completion = completion.first;
                    max_water_cut_perf = water_cut_in_completions[complnumIdx];
                }
                complnumIdx++;
            }

            assert(max_water_cut_limit != 0.);
            assert(worst_offending_completion != INVALIDCOMPLETION);
            violation_extent = max_water_cut_perf / max_water_cut_limit;
        }

        return std::make_tuple(water_cut_limit_violated, worst_offending_completion, violation_extent);
    }





    template<typename TypeTag>
    typename WellInterface<TypeTag>::RatioCheckTuple
    WellInterface<TypeTag>::
    checkRatioEconLimits(const WellEconProductionLimits& econ_production_limits,
                         const WellState& well_state,
                         Opm::DeferredLogger& deferred_logger) const
    {
        // TODO: not sure how to define the worst-offending completion when more than one
        //       ratio related limit is violated.
        //       The defintion used here is that we define the violation extent based on the
        //       ratio between the value and the corresponding limit.
        //       For each violated limit, we decide the worst-offending completion separately.
        //       Among the worst-offending completions, we use the one has the biggest violation
        //       extent.

        bool any_limit_violated = false;
        int worst_offending_completion = INVALIDCOMPLETION;
        double violation_extent = -1.0;

        if (econ_production_limits.onMaxWaterCut()) {
            const RatioCheckTuple water_cut_return = checkMaxWaterCutLimit(econ_production_limits, well_state);
            bool water_cut_violated = std::get<0>(water_cut_return);
            if (water_cut_violated) {
                any_limit_violated = true;
                const double violation_extent_water_cut = std::get<2>(water_cut_return);
                if (violation_extent_water_cut > violation_extent) {
                    violation_extent = violation_extent_water_cut;
                    worst_offending_completion = std::get<1>(water_cut_return);
                }
            }
        }

        if (econ_production_limits.onMaxGasOilRatio()) {
            deferred_logger.warning("NOT_SUPPORTING_MAX_GOR", "the support for max Gas-Oil ratio is not implemented yet!");
        }

        if (econ_production_limits.onMaxWaterGasRatio()) {
            deferred_logger.warning("NOT_SUPPORTING_MAX_WGR", "the support for max Water-Gas ratio is not implemented yet!");
        }

        if (econ_production_limits.onMaxGasLiquidRatio()) {
            deferred_logger.warning("NOT_SUPPORTING_MAX_GLR", "the support for max Gas-Liquid ratio is not implemented yet!");
        }

        if (any_limit_violated) {
            assert(worst_offending_completion != INVALIDCOMPLETION);
            assert(violation_extent > 1.);
        }

        return std::make_tuple(any_limit_violated, worst_offending_completion, violation_extent);
    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    updateWellTestState(const WellState& well_state,
                        const double& simulationTime,
                        const bool& writeMessageToOPMLog,
                        WellTestState& wellTestState,
                        Opm::DeferredLogger& deferred_logger) const
    {

        // currently, we only updateWellTestState for producers
        if (wellType() != PRODUCER) {
            return;
        }

        // Based on current understanding, only under prediction mode, we need to shut well due to various
        // reasons or limits. With more knowlage or testing cases later, this might need to be corrected.
        if (!underPredictionMode() ) {
            return;
        }

        // updating well test state based on physical (THP/BHP) limits.
        updateWellTestStatePhysical(well_state, simulationTime, writeMessageToOPMLog, wellTestState, deferred_logger);

        // updating well test state based on Economic limits.
        updateWellTestStateEconomic(well_state, simulationTime, writeMessageToOPMLog, wellTestState, deferred_logger);

        // TODO: well can be shut/closed due to other reasons
    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    updateWellTestStatePhysical(const WellState& well_state,
                                const double simulation_time,
                                const bool write_message_to_opmlog,
                                WellTestState& well_test_state,
                                Opm::DeferredLogger& deferred_logger) const
    {
        if (!isOperable()) {
            well_test_state.addClosedWell(name(), WellTestConfig::Reason::PHYSICAL, simulation_time);
            if (write_message_to_opmlog) {
                // TODO: considering auto shut in?
                const std::string msg = "well " + name()
                             + std::string(" will be shut as it can not operate under current reservoir condition");
                deferred_logger.info(msg);
            }
        }

    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    updateWellTestStateEconomic(const WellState& well_state,
                                const double simulation_time,
                                const bool write_message_to_opmlog,
                                WellTestState& well_test_state,
                                Opm::DeferredLogger& deferred_logger) const
    {
        const WellEconProductionLimits& econ_production_limits = well_ecl_->getEconProductionLimits(current_step_);

        // if no limit is effective here, then continue to the next well
        if ( !econ_production_limits.onAnyEffectiveLimit() ) {
            return;
        }

        // flag to check if the mim oil/gas rate limit is violated
        bool rate_limit_violated = false;

        // for the moment, we only handle rate limits, not handling potential limits
        // the potential limits should not be difficult to add
        const WellEcon::QuantityLimitEnum& quantity_limit = econ_production_limits.quantityLimit();
        if (quantity_limit == WellEcon::POTN) {
            const std::string msg = std::string("POTN limit for well ") + name() + std::string(" is not supported for the moment. \n")
                                  + std::string("All the limits will be evaluated based on RATE. ");
            deferred_logger.warning("NOT_SUPPORTING_POTN", msg);
        }

        if (econ_production_limits.onAnyRateLimit()) {
            rate_limit_violated = checkRateEconLimits(econ_production_limits, well_state, deferred_logger);
        }

        if (rate_limit_violated) {
            if (econ_production_limits.endRun()) {
                const std::string warning_message = std::string("ending run after well closed due to economic limits")
                                                  + std::string("is not supported yet \n")
                                                  + std::string("the program will keep running after ") + name()
                                                  + std::string(" is closed");
                deferred_logger.warning("NOT_SUPPORTING_ENDRUN", warning_message);
            }

            if (econ_production_limits.validFollowonWell()) {
                deferred_logger.warning("NOT_SUPPORTING_FOLLOWONWELL", "opening following on well after well closed is not supported yet");
            }

            well_test_state.addClosedWell(name(), WellTestConfig::Reason::ECONOMIC, simulation_time);
            if (write_message_to_opmlog) {
                if (well_ecl_->getAutomaticShutIn()) {
                    const std::string msg = std::string("well ") + name() + std::string(" will be shut due to rate economic limit");
                    deferred_logger.info(msg);
                } else {
                    const std::string msg = std::string("well ") + name() + std::string(" will be stopped due to rate economic limit");
                    deferred_logger.info(msg);
                }
            }
            // the well is closed, not need to check other limits
            return;
        }

        // checking for ratio related limits, mostly all kinds of ratio.
        bool ratio_limits_violated = false;
        RatioCheckTuple ratio_check_return;

        if (econ_production_limits.onAnyRatioLimit()) {
            ratio_check_return = checkRatioEconLimits(econ_production_limits, well_state, deferred_logger);
            ratio_limits_violated = std::get<0>(ratio_check_return);
        }

        if (ratio_limits_violated) {
            const WellEcon::WorkoverEnum workover = econ_production_limits.workover();
            switch (workover) {
                case WellEcon::CON:
                {
                    const int worst_offending_completion = std::get<1>(ratio_check_return);

                    well_test_state.addClosedCompletion(name(), worst_offending_completion, simulation_time);
                    if (write_message_to_opmlog) {
                        if (worst_offending_completion < 0) {
                            const std::string msg = std::string("Connection ") + std::to_string(- worst_offending_completion)
                                    + std::string(" for well ") + name() + std::string(" will be closed due to economic limit");
                            deferred_logger.info(msg);
                        } else {
                            const std::string msg = std::string("Completion ") + std::to_string(worst_offending_completion)
                                    + std::string(" for well ") + name() + std::string(" will be closed due to economic limit");
                            deferred_logger.info(msg);
                        }
                    }

                    bool allCompletionsClosed = true;
                    const auto& connections = well_ecl_->getConnections(current_step_);
                    for (const auto& connection : connections) {
                        if (!well_test_state.hasCompletion(name(), connection.complnum())) {
                            allCompletionsClosed = false;
                        }
                    }

                    if (allCompletionsClosed) {
                        well_test_state.addClosedWell(name(), WellTestConfig::Reason::ECONOMIC, simulation_time);
                        if (write_message_to_opmlog) {
                            if (well_ecl_->getAutomaticShutIn()) {
                                const std::string msg = name() + std::string(" will be shut due to last completion closed");
                            	deferred_logger.info(msg);
                            } else {
                                const std::string msg = name() + std::string(" will be stopped due to last completion closed");
                                deferred_logger.info(msg);
                            }
                        }
                    }
                    break;
                }
                case WellEcon::WELL:
                {
                well_test_state.addClosedWell(name(), WellTestConfig::Reason::ECONOMIC, 0);
                if (write_message_to_opmlog) {
                    if (well_ecl_->getAutomaticShutIn()) {
                        // tell the controll that the well is closed
                        const std::string msg = name() + std::string(" will be shut due to ratio economic limit");
                        deferred_logger.info(msg);
                    } else {
                        const std::string msg = name() + std::string(" will be stopped due to ratio economic limit");
                        deferred_logger.info(msg);
                    }
                }
                    break;
                }
                case WellEcon::NONE:
                    break;
                default:
                {
                    deferred_logger.warning("NOT_SUPPORTED_WORKOVER_TYPE",
                                    "not supporting workover type " + WellEcon::WorkoverEnumToString(workover) );
                }
            }
        }
    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    wellTesting(Simulator& simulator, const std::vector<double>& B_avg,
                const double simulation_time, const int report_step,
                const WellTestConfig::Reason testing_reason,
                /* const */ WellState& well_state,
                WellTestState& well_test_state,
                Opm::DeferredLogger& deferred_logger)
    {
        if (testing_reason == WellTestConfig::Reason::PHYSICAL) {
            wellTestingPhysical(simulator, B_avg, simulation_time, report_step,
                                well_state, well_test_state, deferred_logger);
        }

        if (testing_reason == WellTestConfig::Reason::ECONOMIC) {
            wellTestingEconomic(simulator, B_avg, simulation_time, report_step,
                                well_state, well_test_state, deferred_logger);
        }
    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    wellTestingEconomic(Simulator& simulator, const std::vector<double>& B_avg,
                        const double simulation_time, const int report_step,
                        const WellState& well_state, WellTestState& welltest_state, Opm::DeferredLogger& deferred_logger)
    {
        deferred_logger.debug(" well " + name() + " is being tested for economic limits");

        WellState well_state_copy = well_state;

        updatePrimaryVariables(well_state_copy);
        initPrimaryVariablesEvaluation();

        // create a well
        WellTestState welltest_state_temp;

        bool testWell = true;
        // if a well is closed because all completions are closed, we need to check each completion
        // individually. We first open all completions, then we close one by one by calling updateWellTestState
        // untill the number of closed completions do not increase anymore.
        while (testWell) {
            const size_t original_number_closed_completions = welltest_state_temp.sizeCompletions();
            solveWellForTesting(simulator, well_state_copy, B_avg, deferred_logger);
            updateWellTestState(well_state_copy, simulation_time, /*writeMessageToOPMLog=*/ false, welltest_state_temp, deferred_logger);
            closeCompletions(welltest_state_temp);

            // Stop testing if the well is closed or shut due to all completions shut
            // Also check if number of completions has increased. If the number of closed completions do not increased
            // we stop the testing.
            // TODO: it can be tricky here, if the well is shut/closed due to other reasons
            if ( welltest_state_temp.sizeWells() > 0 ||
                (original_number_closed_completions == welltest_state_temp.sizeCompletions()) ) {
                    testWell = false; // this terminates the while loop
            }
        }

        // update wellTestState if the well test succeeds
        if (!welltest_state_temp.hasWell(name(), WellTestConfig::Reason::ECONOMIC)) {
            welltest_state.openWell(name());
            const std::string msg = std::string("well ") + name() + std::string(" is re-opened");
            deferred_logger.info(msg);

            // also reopen completions
            for (auto& completion : well_ecl_->getCompletions(report_step)) {
                if (!welltest_state_temp.hasCompletion(name(), completion.first)) {
                    welltest_state.dropCompletion(name(), completion.first);
                }
            }
        }
    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    computeRepRadiusPerfLength(const Grid& grid,
                               const std::vector<int>& cartesian_to_compressed)
    {
        const int* cart_dims = Opm::UgGridHelpers::cartDims(grid);
        auto cell_to_faces = Opm::UgGridHelpers::cell2Faces(grid);
        auto begin_face_centroids = Opm::UgGridHelpers::beginFaceCentroids(grid);

        const int nperf = number_of_perforations_;

        perf_rep_radius_.clear();
        perf_length_.clear();
        bore_diameters_.clear();

        perf_rep_radius_.reserve(nperf);
        perf_length_.reserve(nperf);
        bore_diameters_.reserve(nperf);

        // COMPDAT handling
        const auto& connectionSet = well_ecl_->getConnections(current_step_);
        for (size_t c=0; c<connectionSet.size(); c++) {
            const auto& connection = connectionSet.get(c);
            if (connection.state() == WellCompletion::OPEN) {
                const int i = connection.getI();
                const int j = connection.getJ();
                const int k = connection.getK();

                const int* cpgdim = cart_dims;
                const int cart_grid_indx = i + cpgdim[0]*(j + cpgdim[1]*k);
                const int cell = cartesian_to_compressed[cart_grid_indx];

                if (cell < 0) {
                    OPM_THROW(std::runtime_error, "Cell with i,j,k indices " << i << ' ' << j << ' '
                              << k << " not found in grid (well = " << name() << ')');
                }

                {
                    double radius = connection.rw();
                    const std::array<double, 3> cubical =
                    WellsManagerDetail::getCubeDim<3>(cell_to_faces, begin_face_centroids, cell);

                    double re; // area equivalent radius of the grid block
                    double perf_length; // the length of the well perforation

                    switch (connection.dir()) {
                        case Opm::WellCompletion::DirectionEnum::X:
                            re = std::sqrt(cubical[1] * cubical[2] / M_PI);
                            perf_length = cubical[0];
                            break;
                        case Opm::WellCompletion::DirectionEnum::Y:
                            re = std::sqrt(cubical[0] * cubical[2] / M_PI);
                            perf_length = cubical[1];
                            break;
                        case Opm::WellCompletion::DirectionEnum::Z:
                            re = std::sqrt(cubical[0] * cubical[1] / M_PI);
                            perf_length = cubical[2];
                            break;
                        default:
                            OPM_THROW(std::runtime_error, " Dirtecion of well is not supported ");
                    }

                    const double repR = std::sqrt(re * radius);
                    perf_rep_radius_.push_back(repR);
                    perf_length_.push_back(perf_length);
                    bore_diameters_.push_back(2. * radius);
                }
            }
        }
    }

    template<typename TypeTag>
    double
    WellInterface<TypeTag>::scalingFactor(const int phaseIdx) const
    {
        const WellControls* wc = well_controls_;
        const double* distr = well_controls_get_current_distr(wc);

        if (well_controls_get_current_type(wc) == RESERVOIR_RATE) {
            if (has_solvent && phaseIdx == contiSolventEqIdx ) {
                typedef Ewoms::BlackOilSolventModule<TypeTag> SolventModule;
                double coeff = 0;
                rateConverter_.template calcCoeffSolvent<SolventModule>(0, pvtRegionIdx_, coeff);
                return coeff;
            }
            // TODO: use the rateConverter here as well.
            return distr[phaseIdx];
        }
        const auto& pu = phaseUsage();
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) && pu.phase_pos[Water] == phaseIdx)
            return 1.0;
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && pu.phase_pos[Oil] == phaseIdx)
            return 1.0;
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && pu.phase_pos[Gas] == phaseIdx)
            return 0.01;
        if (has_solvent && phaseIdx == contiSolventEqIdx )
            return 0.01;

        // we should not come this far
        assert(false);
        return 1.0;
    }





    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::isVFPActive() const
    {
        // since the well_controls only handles the VFP number when THP constraint/target is there.
        // we need to get the table number through the parser, in case THP constraint/target is not there.
        // When THP control/limit is not active, if available VFP table is provided, we will still need to
        // update THP value. However, it will only used for output purpose.

        if (well_type_ == PRODUCER) { // producer
            const int table_id = well_ecl_->getProductionProperties(current_step_).VFPTableNumber;
            if (table_id <= 0) {
                return false;
            } else {
                if (vfp_properties_->getProd()->hasTable(table_id)) {
                    return true;
                } else {
                    OPM_THROW(std::runtime_error, "VFPPROD table " << std::to_string(table_id) << " is specfied,"
                              << " for well " << name() << ", while we could not access it during simulation");
                }
            }

        } else { // injector
            const int table_id = well_ecl_->getInjectionProperties(current_step_).VFPTableNumber;
            if (table_id <= 0) {
                return false;
            } else {
                if (vfp_properties_->getInj()->hasTable(table_id)) {
                    return true;
                } else {
                    OPM_THROW(std::runtime_error, "VFPINJ table " << std::to_string(table_id) << " is specfied,"
                              << " for well " << name() << ", while we could not access it during simulation");
                }
            }
        }
    }





    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    solveWellEqUntilConverged(Simulator& ebosSimulator,
                              const std::vector<double>& B_avg,
                              WellState& well_state,
                              Opm::DeferredLogger& deferred_logger)
    {
        const int max_iter = param_.max_welleq_iter_;
        int it = 0;
        const double dt = 1.0; //not used for the well tests
        bool converged;
        WellState well_state0 = well_state;
        do {
            assembleWellEq(ebosSimulator, dt, well_state, deferred_logger);

            auto report = getWellConvergence(B_avg);
            converged = report.converged();
            if (converged) {
                break;
            }

            ++it;
            solveEqAndUpdateWellState(well_state);

            updateWellControl(ebosSimulator, well_state, deferred_logger);
            initPrimaryVariablesEvaluation();
        } while (it < max_iter);

        return converged;
    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::calculateReservoirRates(WellState& well_state) const
    {
        const int fipreg = 0; // not considering the region for now
        const int np = number_of_phases_;

        std::vector<double> surface_rates(np, 0.0);
        const int well_rate_index = np * index_of_well_;
        for (int p = 0; p < np; ++p) {
            surface_rates[p] = well_state.wellRates()[well_rate_index + p];
        }

        std::vector<double> voidage_rates(np, 0.0);
        rateConverter_.calcReservoirVoidageRates(fipreg, pvtRegionIdx_, surface_rates, voidage_rates);

        for (int p = 0; p < np; ++p) {
            well_state.wellReservoirRates()[well_rate_index + p] = voidage_rates[p];
        }
    }

    template<typename TypeTag>
    void
    WellInterface<TypeTag>::closeCompletions(WellTestState& wellTestState)
    {
        const auto& connections = well_ecl_->getConnections(current_step_);
        int perfIdx = 0;
        for (const auto& connection : connections) {
            if (wellTestState.hasCompletion(name(), connection.complnum())) {
                well_index_[perfIdx] = 0.0;
            }
            perfIdx++;
        }
    }

    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    solveWellForTesting(Simulator& ebosSimulator, WellState& well_state,
                        const std::vector<double>& B_avg,
                        Opm::DeferredLogger& deferred_logger)
    {
        // keep a copy of the original well state
        const WellState well_state0 = well_state;
        const bool converged = solveWellEqUntilConverged(ebosSimulator, B_avg, well_state, deferred_logger);
        if (converged) {
            deferred_logger.debug("WellTest: Well equation for well " + name() +  " converged");
        } else {
            const int max_iter = param_.max_welleq_iter_;
            deferred_logger.debug("WellTest: Well equation for well " +name() + " failed converging in "
                          + std::to_string(max_iter) + " iterations");
            well_state = well_state0;
        }
    }

    template<typename TypeTag>
    void
    WellInterface<TypeTag>::scaleProductivityIndex(const int perfIdx, double& productivity_index, Opm::DeferredLogger& deferred_logger)
    {

        const auto& connection = well_ecl_->getConnections(current_step_)[perfIdx];

        const bool new_well = well_ecl_->hasEvent(ScheduleEvents::NEW_WELL , current_step_);

        if (well_ecl_->getDrainageRadius(current_step_) < 0) {
            if (new_well && perfIdx == 0) {
                deferred_logger.warning("PRODUCTIVITY_INDEX_WARNING", "Negative drainage radius not supported. The productivity index is set to zero");
            }
            productivity_index = 0.0;
            return;
        }

        if (connection.r0() > well_ecl_->getDrainageRadius(current_step_)) {
            if (new_well && well_productivity_index_logger_counter_ < 1) {
                deferred_logger.info("PRODUCTIVITY_INDEX_INFO", "The effective radius is larger than the well drainage radius for well " + name() +
                             " They are set to equal in the well productivity index calculations");
                well_productivity_index_logger_counter_++;
            }
            return;
        }

        // For zero drainage radius the productivity index is just the transmissibility times the mobility
        if (well_ecl_->getDrainageRadius(current_step_) == 0) {
            return;
        }

        // Scale the productivity index to account for the drainage radius.
        // Assumes steady radial flow only valied for horizontal wells
        productivity_index *=
        (std::log(connection.r0() / connection.rw()) + connection.skinFactor()) /
        (std::log(well_ecl_->getDrainageRadius(current_step_) / connection.rw()) + connection.skinFactor());
    }

    template<typename TypeTag>
    void
    WellInterface<TypeTag>::addCellRates(RateVector& rates, int cellIdx) const
    {
        for (int perfIdx = 0; perfIdx < number_of_perforations_; ++perfIdx) {
            if (cells()[perfIdx] == cellIdx) {
                for (int i = 0; i < RateVector::dimension; ++i) {
                    rates[i] += connectionRates_[perfIdx][i];
                }
            }
        }
    }

    template<typename TypeTag>
    typename WellInterface<TypeTag>::Scalar
    WellInterface<TypeTag>::volumetricSurfaceRateForConnection(int cellIdx, int phaseIdx) const {
        for (int perfIdx = 0; perfIdx < number_of_perforations_; ++perfIdx) {
            if (cells()[perfIdx] == cellIdx) {
                const unsigned activeCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
                return connectionRates_[perfIdx][activeCompIdx].value();
            }
        }
        OPM_THROW(std::invalid_argument, "The well with name " + name()
                  + " does not perforate cell " + std::to_string(cellIdx));
        return 0.0;
    }





    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    isOperable() const {
        return operability_status_.isOperable();
    }


}
