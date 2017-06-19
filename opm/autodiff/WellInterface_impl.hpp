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
        number_of_phases_ = wells->number_of_phases;

        // copying the comp_frac
        {
            comp_frac_.resize(number_of_phases_);
            const int index_begin = index_well * number_of_phases_;
            std::copy(wells->comp_frac + index_begin,
                      wells->comp_frac + index_begin + number_of_phases_, comp_frac_.begin() );
        }

        well_controls_ = wells->ctrls[index_well];

        // perforations related
        {
            const int perf_index_begin = wells->well_connpos[index_well];
            const int perf_index_end = wells->well_connpos[index_well + 1];
            number_of_perforations_ = perf_index_end - perf_index_begin;

            well_cell_.resize(number_of_perforations_);
            std::copy(wells->well_cells + perf_index_begin,
                      wells->well_cells + perf_index_end,
                      well_cell_.begin() );

            well_index_.resize(number_of_perforations_);
            std::copy(wells->WI + perf_index_begin,
                      wells->WI + perf_index_end,
                      well_index_.begin() );

            // TODO: not sure about the processing of depth for perforations here
            // Will revisit here later. There are different ways and the definition for different wells
            // can be different, it is possible that we need to remove this from the WellInterface
            perf_depth_.resize(number_of_perforations_, 0.);
            const auto& completion_set = well->getCompletions(time_step);
            for (int i = 0; i < number_of_perforations_; ++i) {
                perf_depth_[i] = completion_set.get(i).getCenterDepth();
            }
        }
    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    init(const PhaseUsage* phase_usage_arg,
         const std::vector<bool>* active_arg,
         const VFPProperties* vfp_properties_arg,
         const double gravity_arg)
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
    const PhaseUsage&
    WellInterface<TypeTag>::
    phaseUsage() const
    {
        assert(phase_usage_);

        return *phase_usage_;
    }


}
