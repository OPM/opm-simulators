
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
#include <opm/parser/eclipse/EclipseState/Schedule/CompletionSet.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well.hpp>


namespace Opm
{

    WellMultiSegment::WellMultiSegment(const Well* well, size_t time_step, const Wells* wells) {
        m_well_name_ = well->name();
        if (well->isMultiSegment(time_step)) {
            initMultiSegmentWell(well, time_step, wells);
        } else {
            initNonMultiSegmentWell(well, time_step, wells);
        }
        updateWellOps();
    }

    void WellMultiSegment::initMultiSegmentWell(const Well* well, size_t time_step, const Wells* wells) {

        const auto& completion_set = well->getCompletions(time_step);

        m_is_multi_segment_ = true;
        const auto& segment_set = well->getSegmentSet(time_step);
        m_number_of_segments_ = segment_set.numberSegment();
        m_number_of_perforations_ = completion_set.size();
        m_comp_pressure_drop_ = segment_set.compPressureDrop();
        m_multiphase_model_ = segment_set.multiPhaseModel();

        m_outlet_segment_.resize(m_number_of_segments_);
        m_inlet_segments_.resize(m_number_of_segments_);
        m_segment_length_.resize(m_number_of_segments_);
        m_segment_depth_.resize(m_number_of_segments_);
        m_segment_internal_diameter_.resize(m_number_of_segments_);
        m_segment_roughness_.resize(m_number_of_segments_);
        m_segment_cross_area_.resize(m_number_of_segments_, 0.);
        m_segment_volume_.resize(m_number_of_segments_);
        m_segment_perforations_.resize(m_number_of_segments_);

        // we change the ID to location now for easier use later.
        for (int i = 0; i < m_number_of_segments_; ++i) {
            // The segment number for top segment is 0, the segment number of its outlet segment will be -1
            m_outlet_segment_[i] = segment_set.segmentNumberToIndex(segment_set[i].outletSegment());
            m_segment_length_[i] = segment_set[i].totalLength();
            m_segment_depth_[i] = segment_set[i].depth();
            m_segment_internal_diameter_[i] = segment_set[i].internalDiameter();
            m_segment_roughness_[i] = segment_set[i].roughness();
            m_segment_cross_area_[i] = segment_set[i].crossArea();
            m_segment_volume_[i] = segment_set[i].volume();
        }

        // update the completion related information
        // find the location of the well in wells
        int index_well;
        for (index_well = 0; index_well < wells->number_of_wells; ++index_well) {
            if (m_well_name_ == std::string(wells->name[index_well])) {
                break;
            }
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


            for(int i = index_begin; i < index_end; ++i) {
                // copy the WI and  well_cell_ informatin to m_well_index_ and m_well_cell_
                // maybe also the depth of the perforations.
                temp_well_cell.push_back(wells->well_cells[i]);
                temp_well_index.push_back(wells->WI[i]);
            }
        }

        std::vector<double> temp_perf_depth;
        temp_perf_depth.resize(m_number_of_perforations_);

        for (int i = 0; i < (int)completion_set.size(); ++i) {
            int i_segment = completion_set.get(i).getSegmentNumber();
            // using the location of the segment in the array as the segment number/id.
            // TODO: it can be helpful for output or postprocessing if we can keep the original number.
            i_segment = segment_set.segmentNumberToIndex(i_segment);
            m_segment_perforations_[i_segment].push_back(i);
            temp_perf_depth[i] = completion_set.get(i).getCenterDepth();
        }

        // reordering the perforation related informations
        // so that the perforations belong to the same segment will be continuous
        m_well_cell_.resize(m_number_of_perforations_);
        m_well_index_.resize(m_number_of_perforations_);
        m_perf_depth_.resize(m_number_of_perforations_);
        m_segment_cell_.resize(m_number_of_segments_, -1);

        int perf_count = 0;
            for (int is = 0; is < m_number_of_segments_; ++is) {
            // TODO: the grid cell related to a segment should be calculated based on the location
            //       of the segment node.
            //       As the current temporary solution, the grid cell related to a segment determined by the
            //       first perforation cell related to the segment.
            //       when no perforation is related to the segment, use it outlet segment's cell.
            const int nperf = m_segment_perforations_[is].size();
            if (nperf > 0) {
                const int first_perf_number = m_segment_perforations_[is][0];
                m_segment_cell_[is] = temp_well_cell[first_perf_number];
                for (int iperf = 0; iperf < nperf; ++iperf) {
                    const int perf_number = m_segment_perforations_[is][iperf];
                    m_well_cell_[perf_count] = temp_well_cell[perf_number];
                    m_well_index_[perf_count] = temp_well_index[perf_number];
                    m_perf_depth_[perf_count] = temp_perf_depth[perf_number];
                    m_segment_perforations_[is][iperf] = perf_count;
                    ++perf_count;
                }
            } else {
                // using the cell of its outlet segment
                const int i_outlet_segment = m_outlet_segment_[is];
                if (i_outlet_segment < 0) {
                    assert(is == 0); // it must be the top segment
                    OPM_THROW(std::logic_error, "Top segment is not related to any perforation, its related cell must be calculated based the location of its segment node, which is not implemented yet \n");
                } else {
                    if (m_well_cell_[i_outlet_segment] < 0) {
                    OPM_THROW(std::logic_error, "The segment cell of its outlet segment is not determined yet, the current implementation does not support this \n");
                    } else {
                        m_segment_cell_[is] = m_segment_cell_[i_outlet_segment];
                    }
                }
            }
        }

        assert(perf_count == m_number_of_perforations_);

        // update m_inlet_segments_
        // top segment does not have a outlet segment
        for (int is = 1; is < m_number_of_segments_; ++is) {
            const int index_outlet = m_outlet_segment_[is];
            m_inlet_segments_[index_outlet].push_back(is);
        }

    }

    void WellMultiSegment::initNonMultiSegmentWell(const Well* well, size_t time_step, const Wells* wells) {

        const auto& completion_set = well->getCompletions(time_step);

        m_is_multi_segment_ = false;
        m_number_of_segments_ = 1;
        m_comp_pressure_drop_ = WellSegment::H__;
        m_multiphase_model_ = WellSegment::HO;

        m_outlet_segment_.resize(m_number_of_segments_, -1);
        m_segment_length_.resize(m_number_of_segments_, 0.);
        m_segment_depth_.resize(m_number_of_segments_, 0.);
        m_segment_internal_diameter_.resize(m_number_of_segments_, 0.);
        m_segment_roughness_.resize(m_number_of_segments_, 0.);
        m_segment_cross_area_.resize(m_number_of_segments_, 0.);
        m_segment_volume_.resize(m_number_of_segments_, 0.);
        m_segment_perforations_.resize(m_number_of_segments_);

        // update the completion related information
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
            // set the segment depth to be the bhp reference depth
            m_segment_depth_[0] = wells->depth_ref[index_well];
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
            m_segment_cell_.resize(1, -1);
            m_segment_cell_[0] = m_well_cell_[0];
        }

        // TODO: not sure if we need the perf_depth_.
        m_perf_depth_.resize(m_number_of_perforations_, 0.);
        m_segment_perforations_[0].resize(m_number_of_perforations_);

        for (int i = 0; i < m_number_of_perforations_; ++i) {
            m_segment_perforations_[0][i] = i;
            m_perf_depth_[i] = completion_set.get(i).getCenterDepth();
        }

        m_inlet_segments_.resize(m_number_of_segments_);
    }

    void WellMultiSegment::updateWellOps() {
        m_wops_.s2p = Matrix(m_number_of_perforations_, m_number_of_segments_);
        m_wops_.p2s = Matrix(m_number_of_segments_, m_number_of_perforations_);

        typedef Eigen::Triplet<double> Tri;

        std::vector<Tri> s2p;
        std::vector<Tri> p2s;

        s2p.reserve(m_number_of_perforations_);
        p2s.reserve(m_number_of_perforations_);

        for(int s = 0; s < (int)m_number_of_segments_; ++s) {
            int temp_nperf = m_segment_perforations_[s].size();
            // some segment may not have any perforation
            assert(temp_nperf >= 0);
            for (int perf = 0; perf < temp_nperf; ++perf) {
                const int index_perf = m_segment_perforations_[s][perf];
                s2p.push_back(Tri(index_perf, s, 1.0));
                p2s.push_back(Tri(s, index_perf, 1.0));
            }
        }
        m_wops_.s2p.setFromTriplets(s2p.begin(), s2p.end());
        m_wops_.p2s.setFromTriplets(p2s.begin(), p2s.end());

        m_wops_.s2s_gather = Matrix(m_number_of_segments_, m_number_of_segments_);
        std::vector<Tri> s2s_gather;

        s2s_gather.reserve(m_number_of_segments_ * m_number_of_segments_);

        std::vector<Tri> s2s_inlets;
        s2s_inlets.reserve(m_number_of_segments_);
        std::vector<Tri> s2s_outlet;
        s2s_outlet.reserve(m_number_of_segments_);

        for (int s = 0; s < (int)m_number_of_segments_; ++s) {
            s2s_gather.push_back(Tri(s, s, 1.0));
            int s_outlet = m_outlet_segment_[s];
            if (s_outlet >=0) {
                s2s_inlets.push_back(Tri(s_outlet, s, 1.0));
                s2s_outlet.push_back(Tri(s, s_outlet, 1.0));
            }
            int temp_s = s;
            while (m_outlet_segment_[temp_s] >=0) {
                s2s_gather.push_back(Tri(m_outlet_segment_[temp_s], s, 1.0));
                temp_s = m_outlet_segment_[temp_s];
            }
        }

        m_wops_.s2s_gather.setFromTriplets(s2s_gather.begin(), s2s_gather.end());

        m_wops_.p2s_gather = Matrix(m_number_of_segments_, m_number_of_perforations_);
        m_wops_.p2s_gather = m_wops_.s2s_gather * m_wops_.p2s;

        m_wops_.s2s_inlets = Matrix(m_number_of_segments_, m_number_of_segments_);
        m_wops_.s2s_inlets.setFromTriplets(s2s_inlets.begin(), s2s_inlets.end());

        m_wops_.s2s_outlet = Matrix(m_number_of_segments_, m_number_of_segments_);
        m_wops_.s2s_outlet.setFromTriplets(s2s_outlet.begin(), s2s_outlet.end());

        m_wops_.p2s_average = Matrix(m_number_of_segments_, m_number_of_perforations_);
        std::vector<Tri> p2s_average;
        p2s_average.reserve(m_number_of_segments_);

        for (int s = 0; s < (int)m_number_of_segments_; ++s) {
            const int nperf = m_segment_perforations_[s].size();
            if (nperf > 0) {
                p2s_average.push_back(Tri(s, s, 1.0/nperf));
            }
        }

        // constructing the diagonal matrix to do the averaging for p2s
        Matrix temp_averaging_p2s = Matrix(m_number_of_segments_, m_number_of_segments_);
        temp_averaging_p2s.setFromTriplets(p2s_average.begin(), p2s_average.end());
        m_wops_.p2s_average = temp_averaging_p2s * m_wops_.p2s;
    }

    const std::string& WellMultiSegment::name() const {
        return m_well_name_;
    }

    bool WellMultiSegment::isMultiSegmented() const {
        return m_is_multi_segment_;
    }

    WellType WellMultiSegment::wellType() const {
        return m_well_type_;
    }

    const WellControls* WellMultiSegment::wellControls() const {
        return m_well_controls_;
    }

    int WellMultiSegment::numberOfPerforations() const {
        return m_number_of_perforations_;
    }

    int WellMultiSegment::numberOfSegments() const {
        return m_number_of_segments_;
    }

    std::string WellMultiSegment::compPressureDrop() const {
        return WellSegment::CompPressureDropEnumToString(m_comp_pressure_drop_);
    }

    const std::vector<double>& WellMultiSegment::compFrac() const {
        return m_comp_frac_;
    }

    int WellMultiSegment::numberOfPhases() const {
        return m_number_of_phases_;
    }

    const std::vector<double>& WellMultiSegment::wellIndex() const {
        return m_well_index_;
    }

    const std::vector<double>& WellMultiSegment::perfDepth() const {
        return m_perf_depth_;
    }

    const std::vector<int>& WellMultiSegment::wellCells() const {
        return m_well_cell_;
    }

    const std::vector<int>& WellMultiSegment::segmentCells() const {
        return m_segment_cell_;
    }

    const std::vector<int>& WellMultiSegment::outletSegment() const {
        return m_outlet_segment_;
    }

    const std::vector<std::vector<int>>& WellMultiSegment::inletSegments() const {
        return m_inlet_segments_;
    }

    const std::vector<double>& WellMultiSegment::segmentLength() const {
        return m_segment_length_;
    }

    const std::vector<double>& WellMultiSegment::segmentDepth() const {
        return m_segment_depth_;
    }

    const std::vector<double>& WellMultiSegment::segmentDiameter() const {
        return m_segment_internal_diameter_;
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

    const std::vector<std::vector<int>>& WellMultiSegment::segmentPerforations() const {
        return m_segment_perforations_;
    }

    const WellMultiSegment::WellOps& WellMultiSegment::wellOps() const {
        return m_wops_;
    }
}
