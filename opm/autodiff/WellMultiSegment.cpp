
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

    WellMultiSegment::WellMultiSegment(WellConstPtr well, size_t time_step, const Wells* wells) {
        init(well, time_step, wells);
    }


    // how to find the order of the well?
    void WellMultiSegment::init(WellConstPtr well, size_t time_step, const Wells* wells) {
        m_well_name_ = well->name();
        CompletionSetConstPtr completion_set = well->getCompletions(time_step);

        if (well->isMultiSegment()) {
            m_is_multi_segment_ = true;
            SegmentSetConstPtr segment_set = well->getSegmentSet(time_step);
            m_number_of_segments_ = segment_set->numberSegment();
            m_number_of_perforations_ = completion_set->size();
            m_comp_pressure_drop_ = segment_set->compPressureDrop();
            m_multiphase_model_ = segment_set->multiPhaseModel();

            // m_number_of_perforations_ from wells
            // m_well_index_ from wells
            m_outlet_segment_.resize(m_number_of_segments_);
            m_segment_length_.resize(m_number_of_segments_);
            m_segment_depth_.resize(m_number_of_segments_);
            m_segment_internal_diameter_.resize(m_number_of_segments_);
            m_segment_roughness_.resize(m_number_of_segments_);
            m_segment_volume_.resize(m_number_of_segments_);
            m_segment_perforations_.resize(m_number_of_segments_);
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
            int index_well;
            for (index_well = 0; index_well < wells->number_of_wells; ++index_well) {
                if (m_well_name_ == std::string(wells->name[index_well])) {
                    break;
                }
                 // std::cout << std::string(wells->name[i]) << "1" << std::endl;
            }

            std::vector<int> temp_well_cell;
            std::vector<double> temp_well_index;

            if (index_well == wells->number_of_wells) {
                throw std::runtime_error(" did not find the well  " + m_well_name_ + "\n");
            } else {

                m_well_type_ = wells->type[index_well];
                m_well_controls_ = wells->ctrls[index_well];
                m_number_of_phases_ = wells->number_of_phases;
                m_comp_frac_.resize(m_number_of_phases_);
                std::copy(wells->comp_frac + index_well * m_number_of_phases_,
                        wells->comp_frac + (index_well + 1) * m_number_of_phases_, m_comp_frac_.begin());

                int index_begin = wells->well_connpos[index_well];
                int index_end = wells->well_connpos[index_well + 1];
                // m_number_of_perforations_ = index_end - index_begin;


                for(int i = index_begin; i < index_end; ++i) {
                    temp_well_cell.push_back(wells->well_cells[i]);
                    temp_well_index.push_back(wells->WI[i]);
                // copy the WI and  well_cell_ informatin to m_well_index_ and m_well_cell_
                // maybe also the depth of the perforations.
                }
            }

            std::vector<double> temp_perf_depth;
            temp_perf_depth.resize(m_number_of_perforations_);

            for (size_t i = 0; i < completion_set->size(); ++i) {
                int i_segment = completion_set->get(i)->getSegmentNumber();
                // convert the original segment id to be the current segment id,
                // which is the location in the array.
                i_segment = segment_set->numberToLocation(i_segment);
                m_segment_perforations_[i_segment].push_back(i);
                temp_perf_depth[i] = completion_set->get(i)->getCenterDepth();
                // TODO: how to handle the center depth of the perforation from the COMPSEGS
            }

            // reordering the perforation related informations
            // so that the perforations belong to the same segment will be continuous
            m_well_cell_.resize(m_number_of_perforations_);
            m_well_index_.resize(m_number_of_perforations_);
            m_perf_depth_.resize(m_number_of_perforations_);

            int perf_count = 0;
            for (size_t is = 0; is < m_number_of_segments_; ++is) {
                for (size_t iperf = 0; iperf < m_segment_perforations_[is].size(); ++iperf) {
                    int perf_number = m_segment_perforations_[is][iperf];
                    m_well_cell_[perf_count] = temp_well_cell[perf_number];
                    m_well_index_[perf_count] = temp_well_index[perf_number];
                    m_perf_depth_[perf_count] = temp_perf_depth[perf_number];
                    m_segment_perforations_[is][iperf] = perf_count;
                    ++perf_count;
                }
            }

            assert(perf_count == m_number_of_perforations_);
            // std::cin.ignore();

            // how to build the relation between segments and completions
            // the completion should be in a global array (all the wells) or local scope (this well)?
            // As the first version, probably it should be OK to handle the wells seperately.
            // Even it will be easier to genrate the matrix in an interleaved way.
            // For each segment, a 4X4 block.
            // we handle with in the local scope first
            // need to change the related segments to the new segment ID (location in the arrary)

        } else {
            m_is_multi_segment_ = false;
            m_number_of_segments_ = 1;
            m_comp_pressure_drop_ = WellSegment::H__;
            m_multiphase_model_ = WellSegment::HO;

            m_outlet_segment_.resize(m_number_of_segments_);
            m_segment_length_.resize(m_number_of_segments_);
            m_segment_depth_.resize(m_number_of_segments_);
            m_segment_internal_diameter_.resize(m_number_of_segments_);
            m_segment_roughness_.resize(m_number_of_segments_);
            m_segment_volume_.resize(m_number_of_segments_);
            m_segment_perforations_.resize(m_number_of_segments_);

                // now the segment for top segment is 0, then its outlet segment will be -1
                // it is also the flag to indicate the top segment
            m_outlet_segment_[0] = -1;
            // TODO: decide the following quantities later.
            // m_segment_length_[i] = (*segment_set)[i]->length();
            // m_segment_depth_[i] = (*segment_set)[i]->depth();
            // m_segment_internal_diameter_[i] = (*segment_set)[i]->internalDiameter();
            // m_segment_roughness_[i] = (*segment_set)[i]->roughness();
            // m_segment_volume_[i] = (*segment_set)[i]->volume();

            // update the completion related information
            // find the location of the well in wells first
            // through well names?
            int index_well;
            for (index_well = 0; index_well < wells->number_of_wells; ++index_well) {
                if (m_well_name_ == std::string(wells->name[index_well])) {
                    break;
                }
            }

            if (index_well == wells->number_of_wells) {
                throw std::runtime_error(" did not find the well  " + m_well_name_ + "\n");
            } else {
                m_well_type_ = wells->type[index_well];
                m_well_controls_ = wells->ctrls[index_well];
                m_number_of_phases_ = wells->number_of_phases;
                m_comp_frac_.resize(m_number_of_phases_);
                std::copy(wells->comp_frac + index_well * m_number_of_phases_,
                        wells->comp_frac + (index_well + 1) * m_number_of_phases_, m_comp_frac_.begin());

                int index_begin = wells->well_connpos[index_well];
                int index_end = wells->well_connpos[index_well + 1];
                m_number_of_perforations_ = index_end - index_begin;

                for(int i = index_begin; i < index_end; ++i) {
                    m_well_cell_.push_back(wells->well_cells[i]);
                    m_well_index_.push_back(wells->WI[i]);
                }
            }

            // TODO: not sure if we need the perf_depth_.
            // m_perf_depth_.resize(m_number_of_perforations_);
            m_segment_perforations_[0].resize(m_number_of_perforations_);

            for (size_t i = 0; i < m_number_of_perforations_; ++i) {
                m_segment_perforations_[0][i] = i;
            }

            // std::cin.ignore();

            // how to build the relation between segments and completions
            // the completion should be in a global array (all the wells) or local scope (this well)?
            // As the first version, probably it should be OK to handle the wells seperately.
            // Even it will be easier to genrate the matrix in an interleaved way.
            // For each segment, a 4X4 block.
            // we handle with in the local scope first
            // need to change the related segments to the new segment ID (location in the arrary)


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
            // And also, with Multisegment well, the way to calculate the density of
            // mixture for wellbore hydrostatic head is always AVG.
            // while for usual wells, it is typically SEG.
        }

    }

    const enum WellType WellMultiSegment::wellType() const {
        return m_well_type_;
    }

    const struct WellControls* WellMultiSegment::wellControls() const {
        return m_well_controls_;
    }

    const size_t WellMultiSegment::numberOfPerforations() const {
        return m_number_of_perforations_;
    }

    const size_t WellMultiSegment::numberOfSegment() const {
        return m_number_of_segments_;
    }

    const std::vector<double>& WellMultiSegment::compFrac() const {
        return m_comp_frac_;
    }

    const size_t WellMultiSegment::numberOfPhases() const {
        return m_number_of_phases_;
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
