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


#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <opm/core/wells.h>
#include <opm/core/well_controls.h>
#include <opm/core/simulator/WellState.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/ScheduleEnums.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/MSW/SegmentSet.hpp>
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

        typedef Eigen::SparseMatrix<double> Matrix;

        /// Constructor of WellMultiSegment
        /// \param[in]       well information from EclipseState
        /// \param[in]       current time step
        /// \param[in[       pointer to Wells structure, to be removed eventually
        WellMultiSegment(WellConstPtr well, size_t time_step, const Wells* wells);

        /// Well name.
        const std::string& name() const;

        /// Flag indicating if the well is a multi-segment well.
        bool isMultiSegmented() const;

        /// Number of the perforations.
        int numberOfPerforations() const;

        /// Number of the segments.
        int numberOfSegments() const;

        /// Components of the pressure drop invloved.
        /// HFA Hydrostatic + friction + acceleration
        /// HF- Hydrostatic + friction
        /// H-- Hydrostatic only.
        std::string compPressureDrop() const;

        /// Well control.
        const WellControls* wellControls() const;

        /// Component fractions for each well.
        const std::vector<double>& compFrac() const;

        /// Number of phases.
        int numberOfPhases() const;

        /// Well type.
        WellType wellType() const;

        /// Well productivity index.
        const std::vector<double>& wellIndex() const;

        /// Depth of the perforations.
        const std::vector<double>& perfDepth() const;

        /// Indices of the grid blocks that perforations are completed in.
        const std::vector<int>& wellCells() const;

        /// Indices of the gird blocks that segments locate at.
        const std::vector<int>& segmentCells() const;

        /// Outlet segments, a segment (except top segment) can only have one outlet segment.
        /// For top segment, its outlet segments is -1 always, which means no outlet segment for top segment.
        const std::vector<int>& outletSegment() const;

        /// Inlet segments. a segment can have more than one inlet segments.
        const std::vector<std::vector<int>>& inletSegments() const;

        /// The length of the segment nodes down the wellbore from the reference point.
        const std::vector<double>& segmentLength() const;

        /// The depth of the segment nodes.
        const std::vector<double>& segmentDepth() const;

        /// Tubing internal diameter.
        const std::vector<double>& segmentDiameter() const;

        /// Cross-sectional area.
        const std::vector<double>& segmentCrossArea() const;

        /// Effective absolute roughness of the tube.
        const std::vector<double>& segmentRoughness() const;

        /// Volume of segment.
        const std::vector<double>& segmentVolume() const;

        /// Perforations related to each segment.
        const std::vector<std::vector<int>>& segmentPerforations() const;

        /// Struct for the well operator matrices.
        /// All the operator matrics only apply to the one specifi well.
        struct WellOps {
            Matrix s2p;              // segment -> perf (scatter)
            Matrix p2s;              // perf -> segment (gather)
            Matrix p2s_average;      // perf -> segment (avarage)
            Matrix s2s_gather;       // segment -> segment (in an accumlative way)
                                // means the outlet segments will gather all the contribution
                                // from all the inlet segments in a recurisive way
            Matrix p2s_gather;       // perforation -> segment (in an accumative way)
            Matrix s2s_inlets;       // segment -> its inlet segments
            Matrix s2s_outlet;      // segment -> its outlet segment
        };

        /// Well operator matrics
        const WellOps& wellOps() const;

    private:
        // for the moment, we use the information from wells.
        // TODO: remove the dependency on wells from opm-core.
        void initMultiSegmentWell(WellConstPtr well, size_t time_step, const Wells* wells);
        void initNonMultiSegmentWell(WellConstPtr well, size_t time_step, const Wells* wells);
        void updateWellOps();

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
        int m_number_of_phases_;
        // component fractions for each well
        std::vector<double> m_comp_frac_;
        // controls for this well
        // using pointer for temporary
        // changing it when figuring out how to using it
        struct WellControls *m_well_controls_;
        // components of the pressure drop to be included
        WellSegment::CompPressureDropEnum m_comp_pressure_drop_;
        // multi-phase flow model
        WellSegment::MultiPhaseModelEnum m_multiphase_model_;
        // number of perforation for this well
        int m_number_of_perforations_;
        // number of segments for this well
        int m_number_of_segments_;
        // well index for each completion
        std::vector<double> m_well_index_;
        // depth for each completion // form the keyword COMPSEGS?
        std::vector<double> m_perf_depth_;
        // well cell for each completion
        std::vector<int> m_well_cell_;
        // cell for each segment
        std::vector<int> m_segment_cell_;
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
