
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

#include <opm/autodiff/WellMultiSegment.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/SegmentSet.hpp>


namespace Opm
{

    // how to find the order of the well?
    void WellMultiSegment::init(WellConstPtr well, size_t time_step, const Wells* wells) {
        CompletionSetConstPtr completion_set = well->getCompletions(time_step);

        if (well->isMultiSegment()) {
            m_is_multi_segment_ = true;
            SegmentSetConstPtr segment_set = well->getSegmentSet(time_step);
            m_number_of_segments_ = segment_set->numberSegment();

            m_comp_pressure_drop_ = segment_set->compPressureDrop();
            m_multiphase_model_ = segment_set->multiPhaseModel();

            // m_number_of_perforations_ from wells
            // m_well_index_ from wells
            m_outlet_segment_.resize(m_number_of_segments_);
            m_segment_length_.resize(m_number_of_segments_);
            m_segment_internal_diameter_.resize(m_number_of_segments_);
            m_segment_roughness_.resize(m_number_of_segments_);
            m_segment_volume_.resize(m_number_of_segments_);
            // what is the point to do this?

            // we change the ID to location now for easier use later.
            for (size_t i = 0; i < m_number_of_segments_; ++i) {
                // now the segment for top segment is 0, then its outlet segment will be -1
                // it is also the flag to indicate the top segment
                m_outlet_segment_[i] = segment_set->numberToLocation((*segment_set)[i]->outletSegment());
                m_segment_length_[i] = (*segment_set)[i]->length();
                m_segment_depth_[i] = (*segment_set)[i]->depth();
                m_segment_internal_diameter_[i] = (*segment_set)[i]->internalDiameter();
                m_segment_roughness_[i] = (*segment_set)[i]->roughness();
                m_segment_volume_[i] = (*segment_set)[i]->volume();
            }
            // update the completion related information
            // find the location of the well in wells first
            // through well names?
            for (size_t i = 0; i < completion_set->size(); ++i) {

            }

            // how to build the relation between segments and completions
            // the completion should be in a global array (all the wells) or local scope (this well)?
            // we handle with in the local scope first
            // need to change the related segments to the new segment ID (location in the arrary)

        } else {
            // building a segment structure for the non-segmented well
            // basically, for each well, only one segment, which means the top segment
            // and all the completions will contribute to the top segment.
            // only hydrostatic pressure drop calculation (H--) and homogenous multiphase
            // flow model (HO) will be used.
            // so each well will have four primary variables
            // G_t, F_w, F_g, bhp
            // This is NOT something we are totally sure about at the moment.
            // The major point will be if the pressure difference between the
            // location of bhp reference point is exactly same with the current
            // way to compute the connection pressure difference.

        }
    }

    size_t WellMultiSegment::numberOfPerforations() const {
        return m_number_of_perforations_;
    }

    size_t WellMultiSegment::numberOfSegment() const {
        return m_number_of_segments_;
    }

    const std::vector<double>& WellMultiSegment::wellIndex() const {
        return m_well_index_;
    }

    const std::vector<double>& WellMultiSegment::perfDepth() const {
        return m_perf_depth_;
    }

    const std::vector<int>& WellMultiSegment::wellCell() const {
        return m_well_cell_;
    }

    const std::vector<int>& WellMultiSegment::outletSegment() const {
        return m_outlet_segment_;
    }

    const std::vector<double>& WellMultiSegment::segmentLength() const {
        return m_segment_length_;
    }

    const std::vector<double>& WellMultiSegment::segmentDepth() const {
        return m_segment_depth_;
    }

    const std::vector<double>& WellMultiSegment::segmentCrossArea() const {
        return m_segment_cross_area_;
    }

    const std::vector<double>& WellMultiSegment::segmentRoughness() const {
        return m_segment_roughness_;
    }

    const std::vector<double>& WellMultiSegment::segmentVolume() const {
        return m_segment_volume_;
    }

    const std::vector<std::vector<int>>& WellMultiSegment::segmentPerforatioins() const {
        return m_segment_perforations_;
    }
}
