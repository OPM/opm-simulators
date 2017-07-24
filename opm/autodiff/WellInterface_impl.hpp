/*
  Copyright 2017 SINTEF ICT, Applied Mathematics.
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
    WellInterface(const Well* well, const int time_step, const Wells* wells)
    {

        // TODO: trying to use wells struct as little as possible here, be prepared to
        // remove the wells struct in future
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

        name_ = well_name;
        index_of_well_ = index_well;
        well_type_ = wells->type[index_well];
        allow_cf_ = wells->allow_cf[index_well];
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

            well_cell_.resize(number_of_perforations_);
            std::copy(wells->well_cells + perf_index_begin,
                      wells->well_cells + perf_index_end,
                      well_cell_.begin() );

            well_index_.resize(number_of_perforations_);
            std::copy(wells->WI + perf_index_begin,
                      wells->WI + perf_index_end,
                      well_index_.begin() );

            saturation_table_number_.resize(number_of_perforations_);
            std::copy(wells->sat_table_id + perf_index_begin,
                      wells->sat_table_id + perf_index_end,
                      saturation_table_number_.begin() );


            // TODO: not sure about the processing of depth for perforations here
            // Will revisit here later. There are different ways and the definition for different wells
            // can be different, it is possible that we need to remove this from the WellInterface
            perf_depth_.resize(number_of_perforations_, 0.);
            const auto& completion_set = well->getCompletions(time_step);
            for (int i = 0; i < number_of_perforations_; ++i) {
                perf_depth_[i] = completion_set.get(i).getCenterDepth();
            }
        }

        well_efficiency_factor_ = 1.0;
        // TODO: need to calculate based on wellCollections, or it should happen in the Well Model side.
    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    init(const PhaseUsage* phase_usage_arg,
         const std::vector<bool>* active_arg,
         const VFPProperties* vfp_properties_arg,
         const double gravity_arg,
         const int /* num_cells */)
    {
        phase_usage_ = phase_usage_arg;
        active_ = active_arg;
        vfp_properties_ = vfp_properties_arg;
        gravity_ = gravity_arg;
    }





    template<typename TypeTag>
    const std::string&
    WellInterface<TypeTag>::
    name() const
    {
        return name_;
    }





    template<typename TypeTag>
    int
    WellInterface<TypeTag>::
    indexOfWell() const
    {
        return index_of_well_;
    }





    template<typename TypeTag>
    WellType
    WellInterface<TypeTag>::
    wellType() const
    {
        return well_type_;
    }





    template<typename TypeTag>
    int
    WellInterface<TypeTag>::
    numberOfPhases() const
    {
        return number_of_phases_;
    }




    template<typename TypeTag>
    const std::vector<double>&
    WellInterface<TypeTag>::
    compFrac() const
    {
        return comp_frac_;
    }





    template<typename TypeTag>
    WellControls*
    WellInterface<TypeTag>::
    wellControls() const
    {
        return well_controls_;
    }





    template<typename TypeTag>
    const std::vector<int>&
    WellInterface<TypeTag>::
    saturationTableNumber() const
    {
        return saturation_table_number_;
    }





    template<typename TypeTag>
    int
    WellInterface<TypeTag>::
    numberOfPerforations() const
    {
        return number_of_perforations_;
    }





    template<typename TypeTag>
    const std::vector<double>&
    WellInterface<TypeTag>::
    wellIndex() const
    {
        return well_index_;
    }





    template<typename TypeTag>
    const std::vector<double>&
    WellInterface<TypeTag>::
    perfDepth() const
    {
        return perf_depth_;
    }





    template<typename TypeTag>
    const std::vector<int>&
    WellInterface<TypeTag>::
    wellCells() const
    {
        return well_cell_;
    }





    template<typename TypeTag>
    const std::vector<bool>&
    WellInterface<TypeTag>::
    active() const
    {
        assert(active_);

        return *active_;
    }





    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    allowCrossFlow() const
    {
        return allow_cf_;
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
        const int phaseToComp[ 3 ] = { FluidSystem::waterCompIdx, FluidSystem::oilCompIdx, FluidSystem::gasCompIdx};
        if (phaseIdx > 2 )
            return phaseIdx;
        return phaseToComp[ phaseIdx ];
    }





    template<typename TypeTag>
    int
    WellInterface<TypeTag>::
    flowToEbosPvIdx( const int flowPv ) const
    {
        const int flowToEbos[ 3 ] = {
            BlackoilIndices::pressureSwitchIdx,
            BlackoilIndices::waterSaturationIdx,
            BlackoilIndices::compositionSwitchIdx
        };

        if (flowPv > 2 )
            return flowPv;

        return flowToEbos[ flowPv ];
    }





    template<typename TypeTag>
    int
    WellInterface<TypeTag>::
    flowPhaseToEbosPhaseIdx( const int phaseIdx ) const
    {
        assert(phaseIdx < 3);
        const int flowToEbos[ 3 ] = { FluidSystem::waterPhaseIdx, FluidSystem::oilPhaseIdx, FluidSystem::gasPhaseIdx };
        return flowToEbos[ phaseIdx ];
    }





    template<typename TypeTag>
    int
    WellInterface<TypeTag>::
    numPhases() const
    {
        return number_of_phases_;
    }





    template<typename TypeTag>
    int
    WellInterface<TypeTag>::
    numComponents() const
    {
        if (numPhases() == 2) {
                return 2;
        }

        int numComp = FluidSystem::numComponents;

        if (has_solvent) {
                    numComp ++;
        }
        return numComp;
    }




    template<typename TypeTag>
    double
    WellInterface<TypeTag>::
    wsolvent() const
    {
        // TODO: not handling it for the moment
        // TODO: it needs information from the well_ecl
        // TODO: will decide on well_ecl role later.
        // It can be just one member variable and no need to deal with well_ecl at all
        return 0.0;
    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    computeWellPotentials(const Simulator& ebosSimulator,
                          const WellState& well_state,
                          std::vector<double>& well_potentials) const
    {
        const int np = numberOfPhases();

        well_potentials.resize(np, 0.0);

        // get the bhp value based on the bhp constraints
        const double bhp = mostStrictBhpFromBhpLimits();

        // does the well have a THP related constraint?
        if ( !wellHasTHPConstraints() ) {
            assert(std::abs(bhp) != std::numeric_limits<double>::max());

            computeWellRatesWithBhp(ebosSimulator, bhp, well_potentials);
        } else {
            // the well has a THP related constraint
            // checking whether a well is newly added, it only happens at the beginning of the report step
            if ( !well_state.isNewWell(index_of_well_) ) {
                for (int p = 0; p < np; ++p) {
                    // This is dangerous for new added well
                    // since we are not handling the initialization correctly for now
                    well_potentials[p] = well_state.wellRates()[index_of_well_ * np + p];
                }
            } else {
                // We need to generate a reasonable rates to start the iteration process
                computeWellRatesWithBhp(ebosSimulator, bhp, well_potentials);
                for (double& value : well_potentials) {
                    // make the value a little safer in case the BHP limits are default ones
                    // TODO: a better way should be a better rescaling based on the investigation of the VFP table.
                    const double rate_safety_scaling_factor = 0.00001;
                    value *= rate_safety_scaling_factor;
                }
            }

            // potentials = computeWellPotentialWithTHP(ebosSimulator, w, bhp, potentials);
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
}
