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

#ifndef OPM_WELLMULTISEGMENT_HEADER_INCLUDED
#define OPM_WELLMULTISEGMENT_HEADER_INCLUDED


#include <opm/core/utility/platform_dependent/disable_warnings.h>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <opm/core/utility/platform_dependent/reenable_warnings.h>

#include <opm/core/wells.h>
#include <opm/core/well_controls.h>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/ScheduleEnums.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/SegmentSet.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well.hpp>
#include <vector>
#include <cassert>
#include <string>
#include <utility>
#include <map>
#include <algorithm>
#include <array>

namespace Opm
{

    class WellMultiSegment
    {
    public:

        typedef Eigen::Array<double, Eigen::Dynamic, 1> V;
        typedef Eigen::SparseMatrix<double> M;

        WellMultiSegment(WellConstPtr well, size_t time_step, const Wells* wells);

        const std::string& name() const;
        const bool isMultiSegmented() const;
        const size_t numberOfPerforations() const;
        const size_t numberOfSegments() const;

        const struct WellControls* wellControls() const;
        const std::vector<double>& compFrac() const;

        const size_t numberOfPhases() const;

        const enum WellType wellType() const;
        const std::vector<double>& wellIndex() const;
        const std::vector<double>& perfDepth() const;
        const std::vector<int>& wellCells() const;
        const std::vector<int>& outletSegment() const;
        const std::vector<std::vector<int>>& inletSegments() const;
        const std::vector<double>& segmentLength() const;
        const std::vector<double>& segmentDepth() const;
        const std::vector<double>& segmentCrossArea() const;
        const std::vector<double>& segmentRoughness() const;
        const std::vector<double>& segmentVolume() const;
        const std::vector<std::vector<int>>& segmentPerforations() const;

        struct WellOps {
            M s2p;              // segment -> perf (scatter)
            M p2s;              // perf -> segment (gather)
            M p2s_average;      // perf -> segment (avarage)
            // M w2p;              // well -> perf (scatter)
            // M p2w;              // perf - > well (gather)
                                // but since only one well, so it is just an arrary
                                // not needed now.
            // M w2s;              // well -> segment (scatter)
            M s2s_gather;       // segment -> segment (in an accumlative way)
                                // means the outlet segments will gather all the contribution
                                // from all the inlet segments in a recurisive way
            M p2s_gather;       // perforation -> segment (in an accumative way)
            M s2s_outlet;       // segment -> their outlet segments
        };

        const WellOps& wellOps() const;

    private:
        // for the moment, we use the information from wells.
        // TODO: remove the dependency on wells from opm-core.
        void init(WellConstPtr well, size_t time_step, const Wells* wells);

    private:
        // name of the well
        std::string m_well_name_;
        // flag to indicate if this well is a
        // multi-segmented well
        bool m_is_multi_segment_;
        // well type
        // INJECTOR or PRODUCER
        enum WellType m_well_type_;
        // number of phases
        size_t m_number_of_phases_;
        // component fractions for each well
        std::vector<double> m_comp_frac_;
        // controls for this well
        // using pointer for temporary
        // changing it when figuring out how to using it
        struct WellControls *m_well_controls_;
        // components of the pressure drop to be included
        WellSegment::CompPresureDropEnum m_comp_pressure_drop_;
        // multi-phase flow model
        WellSegment::MultiPhaseModelEnum m_multiphase_model_;
        // number of perforation for this well
        size_t m_number_of_perforations_;
        // number of segments for this well
        size_t m_number_of_segments_;
        // well index for each completion
        std::vector<double> m_well_index_;
        // depth for each completion // form the keyword COMPSEGS?
        std::vector<double> m_perf_depth_;
        // well cell for each completion
        std::vector<int> m_well_cell_;
        // how to organize the segment structure here?
        // indicate the outlet segment for each segment
        // maybe here we can use the location in the vector
        // at the moment, we still use the ID number
        // then a mapping from the ID number to the actual location will be required
        // The ID is already changed to the location now.
        std::vector<int> m_outlet_segment_;
        // for convinience, we store the inlet segments for each segment
        std::vector<std::vector<int>> m_inlet_segments_;
        // this one is not necessary any more, since the segment number is its location.
        // std::map<int, int> m_number_to_location_;
        // has not decided to use the absolute length from the well head
        // or the length of this single segment
        // using the absolute length first
        std::vector<double> m_segment_length_;
        // the depth of the segmnet node
        std::vector<double> m_segment_depth_;
        // the internal diameter of the segment
        std::vector<double> m_segment_internal_diameter_;
        // the roughness of the segment
        std::vector<double> m_segment_roughness_;
        // the cross area of the segment
        std::vector<double> m_segment_cross_area_;
        // the volume of the segment
        std::vector<double> m_segment_volume_;
        // the completions that is related to each segment
        // the completions's ids are their location in the vector m_well_index_
        // m_well_cell_
        // This is also assuming the order of the completions in Well is the same with
        // the order of the completions in wells.
        std::vector<std::vector<int>> m_segment_perforations_;

        WellOps m_wops_;
    };

    typedef std::shared_ptr<WellMultiSegment> WellMultiSegmentPtr;
    typedef std::shared_ptr<const WellMultiSegment> WellMultiSegmentConstPtr;

} // namespace Opm


#endif // OPM_WELLMULTISEGMENT_HEADER_INCLUDE
