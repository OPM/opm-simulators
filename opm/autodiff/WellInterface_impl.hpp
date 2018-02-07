/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.

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
    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    init(const PhaseUsage* phase_usage_arg,
         const std::vector<double>& /* depth_arg */,
         const double gravity_arg,
         const size_t /* num_cells */)
    {
        phase_usage_ = phase_usage_arg;
        gravity_ = gravity_arg;
    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    setVFPProperties(const VFPProperties* vfp_properties_arg)
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
        const int nwc = well_controls_get_num(well_controls_);
        for (int ctrl_index = 0; ctrl_index < nwc; ++ctrl_index) {
            if (well_controls_iget_type(well_controls_, ctrl_index) == THP) {
                return true;
            }
        }
        return false;
    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    updateWellControl(WellState& well_state,
                      wellhelpers::WellSwitchingLogger& logger) const
    {
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
                // ctrl_index will be the index of the broken constraint after the loop.
                break;
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
            logger.wellSwitched(name(),
                                well_controls_iget_type(wc, old_control_index),
                                well_controls_iget_type(wc, updated_control_index));
        }

        if (updated_control_index != old_control_index) { //  || well_collection_->groupControlActive()) {
            updateWellStateWithTarget(well_state);
        }
    }





    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    checkRateEconLimits(const WellEconProductionLimits& econ_production_limits,
                        const WellState& well_state) const
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
            OpmLog::warning("NOT_SUPPORTING_MIN_RESERVOIR_FLUID_RATE", "Minimum reservoir fluid production rate limit is not supported yet");
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
        int worst_offending_connection = INVALIDCONNECTION;
        bool last_connection = false;
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

            last_connection = (perf_number == 1);
            if (last_connection) {
                worst_offending_connection = 0;
                violation_extent = water_cut_perf[0] / max_water_cut_limit;
                return std::make_tuple(water_cut_limit_violated, last_connection, worst_offending_connection, violation_extent);
            }

            double max_water_cut_perf = 0.;
            for (int perf = 0; perf < perf_number; ++perf) {
                if (water_cut_perf[perf] > max_water_cut_perf) {
                    worst_offending_connection = perf;
                    max_water_cut_perf = water_cut_perf[perf];
                }
            }

            assert(max_water_cut_perf != 0.);
            assert((worst_offending_connection >= 0) && (worst_offending_connection < perf_number));

            violation_extent = max_water_cut_perf / max_water_cut_limit;
        }

        return std::make_tuple(water_cut_limit_violated, last_connection, worst_offending_connection, violation_extent);
    }





    template<typename TypeTag>
    typename WellInterface<TypeTag>::RatioCheckTuple
    WellInterface<TypeTag>::
    checkRatioEconLimits(const WellEconProductionLimits& econ_production_limits,
                         const WellState& well_state) const
    {
        // TODO: not sure how to define the worst-offending connection when more than one
        //       ratio related limit is violated.
        //       The defintion used here is that we define the violation extent based on the
        //       ratio between the value and the corresponding limit.
        //       For each violated limit, we decide the worst-offending connection separately.
        //       Among the worst-offending connections, we use the one has the biggest violation
        //       extent.

        bool any_limit_violated = false;
        bool last_connection = false;
        int worst_offending_connection = INVALIDCONNECTION;
        double violation_extent = -1.0;

        if (econ_production_limits.onMaxWaterCut()) {
            const RatioCheckTuple water_cut_return = checkMaxWaterCutLimit(econ_production_limits, well_state);
            bool water_cut_violated = std::get<0>(water_cut_return);
            if (water_cut_violated) {
                any_limit_violated = true;
                const double violation_extent_water_cut = std::get<3>(water_cut_return);
                if (violation_extent_water_cut > violation_extent) {
                    violation_extent = violation_extent_water_cut;
                    worst_offending_connection = std::get<2>(water_cut_return);
                    last_connection = std::get<1>(water_cut_return);
                }
            }
        }

        if (econ_production_limits.onMaxGasOilRatio()) {
            OpmLog::warning("NOT_SUPPORTING_MAX_GOR", "the support for max Gas-Oil ratio is not implemented yet!");
        }

        if (econ_production_limits.onMaxWaterGasRatio()) {
            OpmLog::warning("NOT_SUPPORTING_MAX_WGR", "the support for max Water-Gas ratio is not implemented yet!");
        }

        if (econ_production_limits.onMaxGasLiquidRatio()) {
            OpmLog::warning("NOT_SUPPORTING_MAX_GLR", "the support for max Gas-Liquid ratio is not implemented yet!");
        }

        if (any_limit_violated) {
            assert(worst_offending_connection >=0);
            assert(violation_extent > 1.);
        }

        return std::make_tuple(any_limit_violated, last_connection, worst_offending_connection, violation_extent);
    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    updateListEconLimited(const WellState& well_state,
                          DynamicListEconLimited& list_econ_limited) const
    {
        // economic limits only apply for production wells.
        if (wellType() != PRODUCER) {
            return;
        }

        // flag to check if the mim oil/gas rate limit is violated
        bool rate_limit_violated = false;
        const WellEconProductionLimits& econ_production_limits = well_ecl_->getEconProductionLimits(current_step_);

        // if no limit is effective here, then continue to the next well
        if ( !econ_production_limits.onAnyEffectiveLimit() ) {
            return;
        }

        const std::string well_name = name();

        // for the moment, we only handle rate limits, not handling potential limits
        // the potential limits should not be difficult to add
        const WellEcon::QuantityLimitEnum& quantity_limit = econ_production_limits.quantityLimit();
        if (quantity_limit == WellEcon::POTN) {
            const std::string msg = std::string("POTN limit for well ") + well_name + std::string(" is not supported for the moment. \n")
                                  + std::string("All the limits will be evaluated based on RATE. ");
            OpmLog::warning("NOT_SUPPORTING_POTN", msg);
        }

        if (econ_production_limits.onAnyRateLimit()) {
            rate_limit_violated = checkRateEconLimits(econ_production_limits, well_state);
        }

        if (rate_limit_violated) {
            if (econ_production_limits.endRun()) {
                const std::string warning_message = std::string("ending run after well closed due to economic limits is not supported yet \n")
                                                  + std::string("the program will keep running after ") + well_name + std::string(" is closed");
                OpmLog::warning("NOT_SUPPORTING_ENDRUN", warning_message);
            }

            if (econ_production_limits.validFollowonWell()) {
                OpmLog::warning("NOT_SUPPORTING_FOLLOWONWELL", "opening following on well after well closed is not supported yet");
            }

            if (well_ecl_->getAutomaticShutIn()) {
                list_econ_limited.addShutWell(well_name);
                const std::string msg = std::string("well ") + well_name + std::string(" will be shut in due to rate economic limit");
                    OpmLog::info(msg);
            } else {
                list_econ_limited.addStoppedWell(well_name);
                const std::string msg = std::string("well ") + well_name + std::string(" will be stopped due to rate economic limit");
                OpmLog::info(msg);
            }
            // the well is closed, not need to check other limits
            return;
        }

        // checking for ratio related limits, mostly all kinds of ratio.
        bool ratio_limits_violated = false;
        RatioCheckTuple ratio_check_return;

        if (econ_production_limits.onAnyRatioLimit()) {
            ratio_check_return = checkRatioEconLimits(econ_production_limits, well_state);
            ratio_limits_violated = std::get<0>(ratio_check_return);
        }

        if (ratio_limits_violated) {
            const WellEcon::WorkoverEnum workover = econ_production_limits.workover();
            switch (workover) {
                case WellEcon::CON:
                {
                    const bool last_connection = std::get<1>(ratio_check_return);
                    const int worst_offending_connection = std::get<2>(ratio_check_return);

                    assert((worst_offending_connection >= 0) && (worst_offending_connection < number_of_perforations_));

                    const int cell_worst_offending_connection = well_cells_[worst_offending_connection];
                    list_econ_limited.addClosedConnectionsForWell(well_name, cell_worst_offending_connection);
                    const std::string msg = std::string("Connection ") + std::to_string(worst_offending_connection) + std::string(" for well ")
                                            + well_name + std::string(" will be closed due to economic limit");
                    OpmLog::info(msg);

                    if (last_connection) {
                        // TODO: there is more things to check here
                        list_econ_limited.addShutWell(well_name);
                        const std::string msg2 = well_name + std::string(" will be shut due to the last connection closed");
                        OpmLog::info(msg2);
                    }
                    break;
                }
                case WellEcon::WELL:
                {
                    if (well_ecl_->getAutomaticShutIn()) {
                        list_econ_limited.addShutWell(well_name);
                        const std::string msg = well_name + std::string(" will be shut due to ratio economic limit");
                        OpmLog::info(msg);
                    } else {
                        list_econ_limited.addStoppedWell(well_name);
                        const std::string msg = well_name + std::string(" will be stopped due to ratio economic limit");
                        OpmLog::info(msg);
                    }
                    break;
                }
                case WellEcon::NONE:
                    break;
                default:
                {
                    OpmLog::warning("NOT_SUPPORTED_WORKOVER_TYPE", "not supporting workover type " + WellEcon::WorkoverEnumToString(workover) );
                }
            }
        }
    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    computeRepRadiusPerfLength(const Grid& grid,
                               const std::map<int, int>& cartesian_to_compressed)
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
        const auto& completionSet = well_ecl_->getCompletions(current_step_);
        for (size_t c=0; c<completionSet.size(); c++) {
            const auto& completion = completionSet.get(c);
            if (completion.getState() == WellCompletion::OPEN) {
                const int i = completion.getI();
                const int j = completion.getJ();
                const int k = completion.getK();

                const int* cpgdim = cart_dims;
                const int cart_grid_indx = i + cpgdim[0]*(j + cpgdim[1]*k);
                const std::map<int, int>::const_iterator cgit = cartesian_to_compressed.find(cart_grid_indx);
                if (cgit == cartesian_to_compressed.end()) {
                    OPM_THROW(std::runtime_error, "Cell with i,j,k indices " << i << ' ' << j << ' '
                              << k << " not found in grid (well = " << name() << ')');
                }
                const int cell = cgit->second;

                {
                    double radius = 0.5*completion.getDiameter();
                    if (radius <= 0.0) {
                        radius = 0.5*unit::feet;
                        OPM_MESSAGE("**** Warning: Well bore internal radius set to " << radius);
                    }

                    const std::array<double, 3> cubical =
                    WellsManagerDetail::getCubeDim<3>(cell_to_faces, begin_face_centroids, cell);

                    WellCompletion::DirectionEnum direction = completion.getDirection();

                    double re; // area equivalent radius of the grid block
                    double perf_length; // the length of the well perforation

                    switch (direction) {
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

}
