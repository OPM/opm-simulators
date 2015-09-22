/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_WELLSTATEMULTISEGMENT_HEADER_INCLUDED
#define OPM_WELLSTATEMULTISEGMENT_HEADER_INCLUDED


#include <opm/core/wells.h>
#include <opm/core/well_controls.h>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/WellMultiSegment.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>
#include <vector>
#include <cassert>
#include <string>
#include <utility>
#include <map>
#include <algorithm>
#include <array>

namespace Opm
{

    /// The state of a set of multi-sgemnet wells
    /// since we are avoiding to use the old wells structure
    /// it makes it might be a good idea not to relate this State to the WellState
    class WellStateMultiSegment
         : public WellStateFullyImplicitBlackoil
    {
    public:

        typedef WellMultiSegment::V V;

        // typedef std::array< int, 3 >  mapentry_t;
        // typedef std::map< std::string, mapentry_t > WellMapType;
        // this map needs to change a little bit?
        typedef struct {
            int well_number;
            int start_segment;
            int number_of_segments;
            int start_perforation;
            int number_of_perforations;
            std::vector<int> start_perforation_segment; // the starting position of perforation inside the segment
            std::vector<int> number_of_perforations_segment; // the numbers for perforations for the segments
        } MapentryType;

        typedef std::map<std::string, MapentryType> WellMapType;
        // MAYNOT NEED THIS

        /// Allocate and initialize if wells is non-null.  Also tries
        /// to give useful initial values to the bhp(), wellRates()
        /// and perfPhaseRates() fields, depending on controls
        /// the PrevState here must be the same with State
        template <class State, class PrevState>
        void init(const std::vector<WellMutliSegmentConstPtr>& wells, const State& state, const PrevState& prevState)
        {
            const int nw = wells.size();
            if (nw == 0) {
                return;
            }

            const int np = wells[0]->numberOfPhases(); // number of the phases

            int nperf = 0; // the number of the perforations
            int nseg = 0; // the nubmer of the segments

            for (int iw = 0; iw < nw; ++iw) {
                nperf += wells[iw]->numberOfPerforations();
                nseg += wells[iw]->numberOfSegments();
            }

            bhp_.resize(nw);
            thp_.resize(nw);

            // deciding to add the following variables temporarily
            // TODO: making it better later
            np_ = np;
            nseg_ = nseg;
            nperf_ = nperf_;
            nwells_ = nw;

            wellrates_.resize(nw * np, 0.0);

            current_controls_.resize(nw);
            for(int iw = 0; iw < nw; ++iw) {
                current_controls_[iw] = well_controls_get_current(wells[iw]->wellControls());
            }

            for (int iw = 0; iw < nw; ++iw) {
                assert((wells[i]->wellType() == INJECTOR) || (wells[i]->wellType() == PRODUCER));
                const WellControls* ctrl = wells[iw]->wellControls();
            }

            int start_segment = 0;
            int start_perforation = 0;

            perfphaserates_.clear();
            perfphaserates_.resize(nperf * np, 0.0);

            perfpress_.clear();
            perfpress_.resize(nperf, -1.0e100);

            segphaserates_.clear();
            segphaserates_.resize(nseg * np, 0.0);

            segpress_.clear();
            segpress_.resize(nseg, -1.0e100);


            for (int w = 0; w < nw; ++w) {
                assert((wells[w]->wellType() == INJECTOR) || (wells[w]->wellType() == PRODUCER));
                const WellControls* ctrl = wells[w]->wellControls();

                std::string well_name(wells[w]->name());
                // Initialize the wellMap_ here
                MapentryType& wellMapEntry = wellMap_[well_name];
                wellMapEntry.well_number = w;
                wellMapEntry.start_segment = start_segment;
                wellMapEntry.number_of_segments = wells[w]->numberOfSegments();
                wellMapEntry.start_perforation = start_perforation;
                wellMapEntry.number_of_perforations = wells[w]->numberOfPerforations();

                int start_perforation_segment = 0;
                wellMapEntry.start_perforation_segment.resize(wellMapEntry.number_of_segments);
                wellMapEntry.number_of_perforations_segment.resize(wellMapEntry.number_of_segments);

                for (int i = 0; i < wellMapEntry.number_of_segments; ++i) {
                    wellMapEntry.start_perforation_segment[i] = start_perforation_segment;
                    wellMapEntry.number_of_perforations_segment[i] = wells[w]->segmentPerforations()[i].size();
                    start_perforation_segment += wellMapEntry.number_of_perforations_segment[i];
                }
                assert(start_perforation_segment == wellMapEntry.number_of_perforations);


                if (well_controls_well_is_stopped(ctrl)) {
                    // 1. WellRates: 0
                    // 2. Bhp: assign bhp equal to bhp control, if applicable, otherwise
                    // assign equal to first perforation cell pressure.
                    if (well_controls_get_current_type(ctrl) == BHP) {
                        bhp_[w] = well_controls_get_current_target(ctrl);
                    } else {
                        const int first_cell = wells[0]->wellCells()[0];
                        bhp_[w] = state.pressure()[first_cell];
                    }
                    // 3. Thp: assign thp equal to thp control, if applicable,
                    // otherwise assign equal to bhp value.
                    if (well_controls_get_current_type(ctrl) == THP) {
                        thp_[w] = well_controls_get_current_target( ctrl );
                    } else {
                        thp_[w] = bhp_[w];
                    }
                    // 4. Perforation pressures and phase rates
                    // 5. Segment pressures and phase rates

                } else {
                    // Open Wells
                    // 1. Rates: initialize well rates to match controls if type is SURFACE_RATE. Otherwise, we
                    // cannot set the correct value here, so we aasign a small rate with the correct sign so that any
                    // logic depending on that sign will work as expected.
                    if (well_controls_get_current_type(ctrl) == SURFACE_RATE) {
                        const double rate_target = well_controls_get_current_target(ctrl);
                        const double * distr = well_controls_get_current_distr( ctrl );
                        for (int p = 0; p < np; ++p) {
                            wellrates_[np * w + p] = rate_target * distr[p];
                        }
                    } else {
                        const double small_rate = 1e-14;
                        const double sign = (wells[w]->wellType() == INJECTOR) ? 1.0 : -1.0;
                        for (int p = 0; p < np; ++p) {
                            wellrates_[np * w + p] = small_rate * sign;
                        }
                    }

                    // 2. Bhp:
                    if (well_controls_get_current_type(ctrl) == BHP) {
                        bhp_[w] = well_controls_get_current_target(ctrl);
                    } else {
                        const int first_cell = wells[0]->wellCells()[0];
                        const double safety_factor = (wells[w]->wellType() == INJECTOR) ? 1.01 : 0.99;
                        bhp_[w] = safety_factor* state.pressure()[first_cell];
                    }
                    // 3. Thp:
                    if (well_controls_get_current_type(ctrl) == THP) {
                        thp_[w] = well_controls_get_current_target(ctrl);
                    } else {
                        thp_[w] = bhp_[w];
                    }

                    // 4. Perf rates and pressures
                    int number_of_perforations = wellMapEntry.number_of_perforations;
                    for (int i = 0; i < number_of_perforations; ++i) {
                        for (int p = 0; p < np; ++p) {
                            perfphaserates_[np * (i + start_perforation) + p] = wellrates_[np * w + p] / double(number_of_perforations);
                        }
                        perfpress_[i + start_perforation] = state.pressure()[wells[w]->wellCells()[i + start_perforation]];
                    }

                    // 5. Segment rates and pressures
                    int number_of_segments = wellMapEntry.number_of_segments;
                    // the seg_pressure is the same with the first pref_pressure. For the top segment, it is the same with bhp,
                    // when under bhp control.
                    // the seg_rates will related to the sum of the perforation rates, and also trying to keep consistent with the
                    // well rates. Most importantly, the segment rates of the top segment is the same with the well rates
                    for (int i = 0; i < number_of_segments; ++i) {
                        /* for (int p = 0; p < np; ++p) {
                            segphaserates_[np * (i + start_segment) + p] = 0.;
                        } */
                        int first_perforation_segment = start_perforation + wellMapEntry.start_perforation_segment[i];
                        segpress_[i + start_segment] = perfpress_[first_perforation_segment];
                        // the segmnent pressure of the top segment should be the bhp
                    }

                    for (int p = 0; p < np; ++p) {
                        // std::vector<double> v_phase_rates(number_of_perforations);
                        // V v_perf_rates = V::Zero(number_of_perforations);
                        Eigen::VectorXd v_perf_rates(number_of_perforations);
                        for (int i = 0; i < number_of_perforations; ++i) {
                            v_perf_rates[i] = perfphaserates_[np * (i + start_perforation) + p];
                        }

                        Eigen::VectorXd v_segment_rates = wells[w]->wellOps().p2s_gather * v_perf_rates;

                        for (int i = 0; i < number_of_segments; ++i) {
                            segphaserates_[np * (i + start_segment) + p] = v_segment_rates[i];
                        }
                    }
                    // initialize the segmnet rates.
                    // it works in the analog way with the usual wells.

                    // How to initialize the perforation rates and the segment rates.?
                    // Perforation pressures can be set to the pressure of the corresponding grid cells?
                    // deviding the well rates by the number of the perforations
                    // then calculating the segment rate based on the rates for perforations and
                    // make sure the flow rate for the top segment is consistent with the well flow rates
                    // for pressure it is not that trival
                }

                start_segment += wellMapEntry.number_of_segments;
                start_perforation += wellMapEntry.number_of_perforations;

            }

            // Initialize current_controls_.
            // The controls set in the Wells object are treated as defaults,
            // and also used for initial values.
            current_controls_.resize(nw);
            for (int w = 0; w < nw; ++w) {
                current_controls_[w] = well_controls_get_current(wells[w]->wellControls());
            }

            // initialize wells that have been there before
            // order can change so the mapping is based on the well names
            if ( !(prevState.wellMap().empty()) )
            {
                typedef typename WellMapType::const_iterator const_iterator;
                const_iterator end_old = prevState.wellMap().end();
                const_iterator end_this = wellMap().end();

                for (int w = 0; w < nw; ++w) {
                    std::string well_name(wells[w]->name());
                    const_iterator it_old = prevState.wellMap().find(well_name);
                    const_iterator it_this = wellMap().find(well_name);

                    assert(it_this != end_this); // the current well must be present in the current well map

                    if (it_old != end_old) {
                        const int oldIndex = (*it_old).second.well_number;
                        const int newIndex = w;

                        // bhp
                        bhp()[newIndex] = prevState.bhp()[oldIndex];

                        // well rates
                        for( int i=0, idx=newIndex*np, oldidx=oldIndex*np; i<np; ++i, ++idx, ++oldidx )
                        {
                            wellRates()[ idx ] = prevState.wellRates()[ oldidx ];
                        }

                        const int num_perf_old_well = (*it_old).second.number_of_perforations;
                        const int num_seg_old_well = (*it_old).second.number_of_segments;

                        const int num_seg_this_well = (*it_this).second.number_of_segments;
                        const int num_perf_this_well = (*it_this).second.number_of_perforations;

                        // determing if the structure of the wells has been changed by comparing the number of segments and perforations
                        // may not be very safe.
                        // The strategy HAS to be changed later with experiments and analysis
                        if ((num_perf_old_well == num_seg_old_well) && (num_seg_old_well == num_seg_this_well)) {
                            const int old_start_perforation = (*it_old).second.start_perforation;
                            const int old_start_segment = (*it_old).second.start_segment;

                            const int this_start_perforation = (*it_this).second.start_perforation;
                            const int this_start_segment = (*it_this).second.start_segment;

                            // this is not good when the the well rates changed dramatically
                            for (int i = 0; i < num_seg_this_well * np; ++i) {
                                segphaserates_[this_start_segment * np + i] = prevState.segPhaseRates()[old_start_segment * np + i];
                            }

                            for (int i = 0; i < num_perf_this_well * np; ++i) {
                                perfphaserates_[this_start_perforation * np + i] = prevState.perfPhaseRates()[old_start_perforation * np + i];
                            }

                            // perf_pressures_
                            for (int i = 0; i < num_perf_this_well; ++i) {
                                // p
                                perfpress_[this_start_perforation + i] = prevState.perfPress()[old_start_perforation + i];
                            }

                            // segpress_
                            for (int i = 0; i < num_seg_this_well; ++i) {
                                // p
                                segpress_[this_start_segment + i] = prevState.segPress()[old_start_segment + i];
                            }
                            // current controls
                        }

                        // else {
                            // deviding the well rates by the number of the perforations
                            // then calculating the segment rate based on the rates for perforations and
                            // make sure the flow rate for the top segment is consistent with the well flow rates
                            // for pressure it is not that trival
                        // }

/*         typedef struct {
            int well_number;
            int start_segment;
            int number_of_segments;
            int start_perforation;
            int number_of_perforations;
            std::vector<int> start_perforation_segment; // the starting position of perforation inside the segment
            std::vector<int> number_of_perforations_segment; // the numbers for perforations for the segments
        } MapentryType; */

                        // peforation rates
                        // segment rates
                        // It really depends on if the structures on the segments and perforations are changed.
                        // TODO: if it will be reasonable to assume that if the numbers of segments and perforations are same, then
                        // the structures of the wells are not changed.
                        // Obviously it is not true.

                        // for the perforation rates, it is Okay to calculate by deviding the well rates by the perforation numbers.
                        // Then the segment rates are calculated based on the perforation rates and the well rates.
                        // The segment rates of top segments should be the same with the well rates.


                    }
                }


            }
            // TODO: maybe we should store the values of np, nw, nseg, nperf for the states for later use.
            // assert(start_perforation == total_perforation);
            // assert(start_segment == total_segment);
               /*  if (well_controls_well_is_stopped(ctrl)) {
                    // shut well: all the rates are zero.
                } else {
                    // Initialize the phase rates for each perforation by deviding the well rates by the number of perofrations
                    // Then using the perf rates to initialize the rates for the segments
                } */
        }

        std::vector<double>& segPhaseRates() { return segphaserates_; }
        const std::vector<double>& segPhaseRates() const { return segphaserates_; }

        std::vector<double>& segPress() { return segpress_; }
        const std::vector<double>& segPress() const { return segpress_; }

        std::vector<double>& perfPress() { return perfpress_; }
        const std::vector<double>& perfPress() const { return perfpress_; }

        std::vector<double>& perfPhaseRates() { return perfphaserates_; }
        const std::vector<double>& perfPhaseRates() const { return perfphaserates_; }

        std::vector<double>& bhp() { return bhp_; }
        const std::vector<double>& bhp() const { return bhp_; }

        std::vector<double>& thp() { return thp_; }
        const std::vector<double>& thp() const { return thp_; }

        std::vector<double>& wellRates() { return wellrates_; }
        const std::vector<double>& wellRates() const { return wellrates_; }

        std::vector<int>& currentControls() { return current_controls_; }
        const std::vector<int>& currentControls() const { return current_controls_; }

        // wellrate should be the out segment rates for the top segments

        const WellMapType& wellMap() const { return wellMap_; }
        WellMapType& wellMap() { return wellMap_; }

        const int numberOfPhases() const { return np_; }
        const int numberOfSegments() const { return nseg_; }
        const int numberOfPerforations() const { return nperf_; }
        const int numberOfWells() const { return nwells_; }
    private:
        std::vector<double> bhp_;
        std::vector<double> thp_;
        std::vector<double> wellrates_;
        // std::vector<double> temperature_;
        // pressure for the segment nodes
        std::vector<double> segpress_;
        // phase rates for the segments
        std::vector<double> segphaserates_;
        // phase rates for the completions
        std::vector<double> perfphaserates_;
        // pressure for the perforatins
        std::vector<double> perfpress_;
        // TODO: MIGHT NOT USE THE FOLLOWING VARIABLES AT THE
        // fractions for each segments (W, O, G)
        std::vector<double> segphasefrac_;
        // total flow rates for each segments, G_T
        std::vector<double> segtotalrate_;
        std::vector<int> current_controls_;

        WellMapType wellMap_;

        int nseg_;
        int np_;
        int nperf_;
        int nwells_;

    };

} // namespace Opm


#endif // OPM_WELLSTATEMULTISEGMENT_HEADER_INCLUDE
