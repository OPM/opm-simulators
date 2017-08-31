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



#ifndef OPM_WELLSTATEMSWELL_HEADER_INCLUDED
#define OPM_WELLSTATEMSWELL_HEADER_INCLUDED

namespace Opm
{
    // TODO: eventually, there should be a base class for WellState
    // TODO: the copy from old WellState should be trival
    class WellStateMSWell
    {
    public:
        // TODO: looks like will need ResarvoirState.
        // TODO: It is possible, we can not use unique_ptr here if we need dynamic_cast, we should
        // cast it as a reference.
        WellStateMSWell(const Well* well_ecl, const int np)
         : well_name_()
         , nseg_(ms_well.numberOfSegments())
         , nperf_(ms_well.numberOfPerforations())
         , num_phases_(ms_well.numPhases())
         , wellrate_(num_phases_)
         , segphaserate_(nseg_ * num_phases_)
         , segpress_(nseg_)
        {
            const WellControls* ctrl = ms_well.wellControls();
            current_control_ = well_controls_get_current(ctrl);

            if (!well_controls_well_is_stopped(ctrl)) {
                // Open Wells
                // 1. Rates: initialize well rates to match controls if type is SURFACE_RATE. Otherwise, we
                // cannot set the correct value here, so we aasign a small rate with the correct sign so that any
                // logic depending on that sign will work as expected.
                if (well_controls_get_current_type(ctrl) == SURFACE_RATE) {
                    const double rate_target = well_controls_get_current_target(ctrl);
                    const double* distr = well_controls_get_current_distr( ctrl );
                    for (int p = 0; p < num_phases_; ++p) {
                        wellrate_[p] = rate_target * distr[p];
                    }
                } else {
                    const double small_rate = 1e-14;
                    const double sign = (ms_well.wellType() == INJECTOR) ? 1.0 : -1.0;
                    for (int p = 0; p < np; ++p) {
                        wellrate_[p] = small_rate * sign;
                    }
               }

               // 2. Bhp:
               if (well_controls_get_current_type(ctrl) == BHP) {
                    bhp_ = well_controls_get_current_target(ctrl);
               } else {
                    const int first_cell = ms_well.wellCells()[0];
                    const double safety_factor = (ms_well.wellType() == INJECTOR) ? 1.01 : 0.99;
                    bhp_ = safety_factor* state.pressure()[first_cell];
               }
               // 3. Thp:
               if (well_controls_get_current_type(ctrl) == THP) {
                    thp_ = well_controls_get_current_target(ctrl);
               } else {
                    thp_ = bhp_;
               }

               // 4. Perf rates and pressures
                for (int i = 0; i < nperf_; ++i) {
                    for (int p = 0; p < num_phases_; ++p) {
                        perfphaserates_ [np * i + p] = wellrate_[p] / double(nperf_);
                    }
                    const double safety_factor = (ms_well.wellType() == INJECTOR) ? 1.01 : 0.99;
                    const int cell_index = ms_well.wellCells()[i];
                    perfpress_[i] = safety_factor * state.pressure()[cell_index];
                }

                /* // 5. Segment pressures and rates
                // top segment
                segpress_[0] = bhp_;
                for (int i = 1; i < nseg_; ++i) {
                    // the first perforation of the segment
                    const int first_perforation = ms_well. [0];
                    segpress_[i] = perfpress_[first_perforation];
                }
                // segment rates
                for (int i = 1; i < nseg_; ++i) {
                    // basically summ up the perforations
                } */
            } else {
                // Stopped well
                // 1. WellRates: 0
                // 2. Bhp: assign bhp equal to bhp control, if applicable, otherwise
                // assign equal to first perforation cell pressure.
                if (well_controls_get_current_type(ctrl) == BHP) {
                    bhp_[w] = well_controls_get_current_target(ctrl);
                } else {
                    const int first_cell = ms_well.wellCells()[0];
                    bhp_[w] = state.pressure()[first_cell];
                }
               // 3. Thp: assign thp equal to thp control, if applicable,
               // otherwise assign equal to bhp value.
               if (well_controls_get_current_type(ctrl) == THP) {
                    thp()[w] = well_controls_get_current_target( ctrl );
               } else {
                    thp()[w] = bhp()[w];
               }
               // 4. Perforation pressures and phase rates
               // 5. Segment pressures and phase rates
            }
        }


        const std::string& name() const { return well_name_; }

        int numSegments() const { return nseg_; }

        int numPhases() const { return num_phases_; }

        int numPerforations() const { return nperf_; }

    private:
        // well name is the id of the well
        const std::string well_name_;
        // number of segments
        const int nseg_;
        // number of perforations
        const int nperf_;
        // number of phases
        const int num_phases_;
        // bhp
        double bhp_;
        // thp
        double thp_;
        // well rates
        std::vector<double> wellrate_;
        // segment phase rates
        std::vector<double> segphaserate_;
        // segment pressure
        std::vector<double> segpress_;
        // perforation phase rates // TODO: it can go to the base class, since STDWells also have that
        std::vector<double> perfphaserates_;
        // current control
        int current_control_;
    };

    // TODO: only when the segment and perforation infromation, and well type are the same,
    // we can copy the Well State. Otherwise, we can just use the initialized one.
    // we should have a function here bool consistent() here to decide if we can copy
    // it is not a very safe judgement. well types are also needed to be considered, while
    // we do not have the information from the previous WellModel. Not sure whether to add a WellType here.
    // we can also make an experiment to see if we need to copy the old Well State.
    bool copyable(const WellStateMSWell& well_state, const WellStateMSWell& prev_state) {
        return ( well_state.name() == prev_state.name() &&
                 well_state.numSegments() == prev_state.numSegments() &&
                 well_state.numPerforations() == prev_state.numPerforations() &&
                 well_state.numPhases() == prev_state.numPhases() );
    }
}

#endif
