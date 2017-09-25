/*
  Copyright 2016 IRIS

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

#ifndef OPM_WELLSTATEFULLYIMPLICITBLACKOILDENSE_HEADER_INCLUDED
#define OPM_WELLSTATEFULLYIMPLICITBLACKOILDENSE_HEADER_INCLUDED


#include <opm/core/wells.h>
#include <opm/core/well_controls.h>
#include <opm/core/simulator/WellState.hpp>
#include <opm/autodiff/BlackoilModelEnums.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <vector>
#include <cassert>
#include <string>
#include <utility>
#include <map>
#include <algorithm>
#include <array>
#include <cmath>

namespace Opm
{

    /// The state of a set of wells, tailored for use by the fully
    /// implicit blackoil simulator.
    class WellStateFullyImplicitBlackoilDense
        : public WellStateFullyImplicitBlackoil
    {
        typedef WellStateFullyImplicitBlackoil  BaseType;
    public:
        typedef BaseType :: WellMapType WellMapType;

        using BaseType :: wellRates;
        using BaseType :: bhp;
        using BaseType :: perfPress;
        using BaseType :: wellMap;
        using BaseType :: numWells;
        using BaseType :: numPhases;
        using BaseType :: perfPhaseRates;
        using BaseType :: currentControls;

        /// Allocate and initialize if wells is non-null.  Also tries
        /// to give useful initial values to the bhp(), wellRates()
        /// and perfPhaseRates() fields, depending on controls
        template <class PrevWellState>
        void init(const Wells* wells, const std::vector<double>& cellPressures, const PrevWellState& prevState, const PhaseUsage& pu)
        {

            // call init on base class
            BaseType :: init(wells, cellPressures, prevState);


            const int nw = wells->number_of_wells;
            if (nw == 0) {
                return;
            }
            const int nperf = wells->well_connpos[nw];
            perfRateSolvent_.clear();
            perfRateSolvent_.resize(nperf, 0.0);

            if (pu.has_solvent) {

                // intialize wells that have been there before
                // order may change so the mapping is based on the well name
                if( ! prevState.wellMap().empty() )
                {
                    typedef typename WellMapType :: const_iterator const_iterator;
                    const_iterator end = prevState.wellMap().end();
                    for (int w = 0; w < nw; ++w) {
                        std::string name( wells->name[ w ] );
                        const_iterator it = prevState.wellMap().find( name );
                        if( it != end )
                        {
                            const int newIndex = w;

                            // perfSolventRates
                            int oldPerf_idx = (*it).second[ 1 ];
                            const int num_perf_old_well = (*it).second[ 2 ];
                            const int num_perf_this_well = wells->well_connpos[newIndex + 1] - wells->well_connpos[newIndex];
                            if( num_perf_old_well == num_perf_this_well )
                            {
                                for (int perf = wells->well_connpos[ newIndex ];
                                     perf < wells->well_connpos[ newIndex + 1]; ++perf, ++oldPerf_idx )
                                {
                                    perfRateSolvent()[ perf ] = prevState.perfRateSolvent()[ oldPerf_idx ];
                                }
                            }
                        }
                    }
                }
            }
        }

        /// init the MS well related.
        template <typename PrevWellState>
        void initWellStateMSWell(const Wells* wells, const std::vector<const Well*>& wells_ecl,
                                 const int time_step, const PhaseUsage& pu, const PrevWellState& prev_well_state)
        {
            // still using the order in wells
            const int nw = wells->number_of_wells;
            if (nw == 0) {
                return;
            }

            nseg_ = 0.;
            // in the init function, the well rates and perforation rates have been initialized or copied from prevState
            // what we do here, is to set the segment rates and perforation rates
            for (int w = 0; w < nw; ++w) {
                const int nw_wells_ecl = wells_ecl.size();
                int index_well_ecl = 0;
                const std::string well_name(wells->name[w]);
                for (; index_well_ecl < nw_wells_ecl; ++index_well_ecl) {
                    if (well_name == wells_ecl[index_well_ecl]->name()) {
                        break;
                    }
                }

                // It should be able to find in wells_ecl.
                if (index_well_ecl == nw_wells_ecl) {
                    OPM_THROW(std::logic_error, "Could not find well " << well_name << " in wells_ecl ");
                }

                const Well* well_ecl = wells_ecl[index_well_ecl];
                top_segment_loc_.push_back(nseg_);
                if ( !well_ecl->isMultiSegment(time_step) ) { // not multi-segment well
                    nseg_ += 1;
                    segpress_.push_back(bhp()[w]);
                    const int np = numPhases();
                    for (int p = 0; p < np; ++p) {
                        segrates_.push_back(wellRates()[np * w + p]);
                    }
                } else { // it is a multi-segment well
                    const SegmentSet& segment_set = well_ecl->getSegmentSet(time_step);
                    // assuming the oder of the perforations in well_ecl is the same with Wells
                    const CompletionSet& completion_set = well_ecl->getCompletions(time_step);
                    const int nseg = segment_set.numberSegment();
                    const int nperf = completion_set.size();
                    nseg_ += nseg;
                    // we need to know for each segment, how many perforation it has and how many segments using it as outlet_segment
                    // that is why I think we should use a well model to initialize the WellState here
                    std::vector<std::vector<int>> segment_perforations(nseg);
                    for (int perf = 0; perf < nperf; ++perf) {
                        const Completion& completion = completion_set.get(perf);
                        const int segment_number = completion.getSegmentNumber();
                        const int segment_location = segment_set.numberToLocation(segment_number);
                        segment_perforations[segment_location].push_back(perf);
                    }

                    std::vector<std::vector<int>> segment_inlets(nseg);
                    for (int seg = 0; seg < nseg; ++seg) {
                        const Segment& segment = segment_set[seg];
                        const int segment_number = segment.segmentNumber();
                        const int outlet_segment_number = segment.outletSegment();
                        if (outlet_segment_number > 0) {
                            const int segment_location = segment_set.numberToLocation(segment_number);
                            const int outlet_segment_location = segment_set.numberToLocation(outlet_segment_number);
                            segment_inlets[outlet_segment_location].push_back(segment_location);
                        }
                    }


                    // for the segrates_, now it becomes a recursive solution procedure.
                    {
                        const int np = numPhases();
                        const int start_perf = wells->well_connpos[w];
                        const int start_perf_next_well = wells->well_connpos[w + 1];
                        assert(nperf == (start_perf_next_well - start_perf)); // make sure the information from wells_ecl consistent with wells
                        if (pu.phase_used[Gas]) {
                            const int gaspos = pu.phase_pos[Gas];
                            // scale the phase rates for Gas to avoid too bad initial guess for gas fraction
                            // it will probably benefit the standard well too, while it needs to be justified
                            // TODO: to see if this strategy can benefit StandardWell too
                            // TODO: it might cause big problem for gas rate control or if there is gas rate limit
                            for (int perf = 0; perf < nperf; perf++) {
                                const int perf_pos = start_perf + perf;
                                perfPhaseRates()[np * perf_pos + gaspos] *= 100.;
                            }
                        }

                        const std::vector<double> perforation_rates(perfPhaseRates().begin() + np * start_perf,
                                                                    perfPhaseRates().begin() + np * start_perf_next_well); // the perforation rates for this well
                        std::vector<double> segment_rates;
                        calculateSegmentRates(segment_inlets, segment_perforations, perforation_rates, np, 0 /* top segment */, segment_rates);
                        std::copy(segment_rates.begin(), segment_rates.end(), std::back_inserter(segrates_));
                    }

                    // for the segment pressure, the segment pressure is the same with the first perforation
                    // if there is no perforation associated with this segment, it uses the pressure from the outlet segment
                    // which requres the ordering is successful
                    // TODO: maybe also check the relation with the cell?
                    // Not sure what is the best way to handle the initialization, hopefully, the bad initialization can be
                    // improved during the solveWellEq process
                    {
                        // top segment is always the first one, and its pressure is the well bhp
                        segpress_.push_back(bhp()[w]);
                        const int top_segment = top_segment_loc_[w];
                        const int start_perf = wells->well_connpos[w];
                        for (int seg = 1; seg < nseg; ++seg) {
                            if ( !segment_perforations[seg].empty() ) {
                                const int first_perf = segment_perforations[seg][0];
                                segpress_.push_back(perfPress()[start_perf + first_perf]);
                            } else {
                                // segpress_.push_back(bhp); // may not be a good decision
                                // using the outlet segment pressure // it needs the ordering is correct
                                const int outlet_seg = segment_set[seg].outletSegment();
                                segpress_.push_back(segpress_[top_segment + segment_set.numberToLocation(outlet_seg)]);
                            }
                        }
                    }
                }
            }
            assert(int(segpress_.size()) == nseg_);
            assert(int(segrates_.size()) == nseg_ * numPhases() );

            if (!prev_well_state.wellMap().empty()) {
                // copying MS well related
                const auto& end = prev_well_state.wellMap().end();
                const int np = numPhases();
                for (int w = 0; w < nw; ++w) {
                    const std::string name( wells->name[w] );
                    const auto& it = prev_well_state.wellMap().find( name );

                    if (it != end) { // the well is found in the prev_well_state
                        // TODO: the well with same name can change a lot, like they might not have same number of segments
                        // we need to handle that later.
                        // for now, we just copy them.
                        const int old_index_well = (*it).second[0];
                        const int new_index_well = w;
                        const int old_top_segment_location = prev_well_state.topSegmentLocation(old_index_well);
                        const int new_top_segmnet_location = topSegmentLocation(new_index_well);
                        int number_of_segment = 0;
                        if (new_index_well == top_segment_loc_.size() - 1) {
                            number_of_segment = nseg_ - new_top_segmnet_location;
                        } else {
                            number_of_segment = topSegmentLocation(new_index_well + 1) - new_top_segmnet_location;
                        }

                        for (int i = 0; i < number_of_segment * np; ++i) {
                            segrates_[new_top_segmnet_location * np + i] = prev_well_state.segRates()[old_top_segment_location * np + i];
                        }

                        for (int i = 0; i < number_of_segment; ++i) {
                            segpress_[new_top_segmnet_location + i] = prev_well_state.segPress()[old_top_segment_location + i];
                        }
                    }
                }
            }
        }

        static void calculateSegmentRates(const std::vector<std::vector<int>>& segment_inlets, const std::vector<std::vector<int>>&segment_perforations,
                                          const std::vector<double>& perforation_rates, const int np, const int segment, std::vector<double>& segment_rates)
        {
            // the rate of the segment equals to the sum of the contribution from the perforations and inlet segment rates.
            // the first segment is always the top segment, its rates should be equal to the well rates.
            assert(segment_inlets.size() == segment_perforations.size());
            const int nseg = segment_inlets.size();
            if (segment == 0) { // beginning the calculation
                segment_rates.resize(np * nseg, 0.0);
            }
            // contributions from the perforations belong to this segment
            for (const int& perf : segment_perforations[segment]) {
                for (int p = 0; p < np; ++p) {
                    segment_rates[np * segment + p] += perforation_rates[np * perf + p];
                }
            }
            for (const int& inlet_seg : segment_inlets[segment]) {
                calculateSegmentRates(segment_inlets, segment_perforations, perforation_rates, np, inlet_seg, segment_rates);
                for (int p = 0; p < np; ++p) {
                    segment_rates[np * segment + p] += segment_rates[np * inlet_seg + p];
                }
            }
        }

        template <class State>
        void resize(const Wells* wells, const State& state, const PhaseUsage& pu ) {
            const WellStateFullyImplicitBlackoilDense dummy_state{}; // Init with an empty previous state only resizes
            init(wells, state.pressure(), dummy_state, pu) ;
        }

        /// One rate pr well connection.
        std::vector<double>& perfRateSolvent() { return perfRateSolvent_; }
        const std::vector<double>& perfRateSolvent() const { return perfRateSolvent_; }

        /// One rate pr well
        double solventWellRate(const int w) const {
            double solvent_well_rate = 0.0;
            for (int perf = wells_->well_connpos[w]; perf < wells_->well_connpos[w+1]; ++perf ) {
                solvent_well_rate += perfRateSolvent_[perf];
            }
            return solvent_well_rate;
        }


        data::Wells report(const PhaseUsage& pu) const override {
            data::Wells res = BaseType::report(pu);
            const int nw = WellState::numWells();
            // If there are now wells numPhases throws a floating point
            // exception.
            if (nw == 0) {
                return res;
            }
            if (pu.has_solvent) {
                // add solvent component
                for( int w = 0; w < nw; ++w ) {
                    using rt = data::Rates::opt;
                    res.at( wells_->name[ w ]).rates.set( rt::solvent, solventWellRate(w) );
                }
            }

            return res;
        }

        const std::vector<double>& segRates() const
        {
            return segrates_;
        }

        std::vector<double>& segRates()
        {
            return segrates_;
        }

        const std::vector<double>& segPress() const
        {
            return segpress_;
        }

        std::vector<double>& segPress()
        {
            return segpress_;
        }

        int numSegment() const
        {
            return nseg_;
        }

        int topSegmentLocation(const int w) const
        {
            assert(w < int(top_segment_loc_.size()) );
            return top_segment_loc_[w];
        }


    private:
        std::vector<double> perfRateSolvent_;

        // MS well related
        // for StandardWell, the segment number will be one
        std::vector<double> segrates_;
        std::vector<double> segpress_;
        std::vector<int> top_segment_loc_; // the index of the top segments
        int nseg_; // number of the segments

    };

} // namespace Opm


#endif // OPM_WELLSTATEFULLYIMPLICITBLACKOILDENSE_HEADER_INCLUDED
