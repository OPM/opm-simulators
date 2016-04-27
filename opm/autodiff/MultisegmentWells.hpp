/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 Statoil ASA.

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


#ifndef OPM_MULTISEGMENTWELLS_HEADER_INCLUDED
#define OPM_MULTISEGMENTWELLS_HEADER_INCLUDED

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <cassert>

#include <opm/core/props/BlackoilPhases.hpp>

#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/autodiff/BlackoilModelEnums.hpp>
#include <opm/autodiff/BlackoilPropsAdInterface.hpp>

#include <opm/autodiff/WellMultiSegment.hpp>



namespace Opm {


    /// Class for handling the multi-segment well model
    class MultisegmentWells {
    public:
            // ---------      Types      ---------
            using ADB = AutoDiffBlock<double>;
            using Vector = ADB::V;

            // Well operations and data needed.
            struct MultisegmentWellOps {
                explicit MultisegmentWellOps(const std::vector<WellMultiSegmentConstPtr>& wells_ms);
                Eigen::SparseMatrix<double> w2p;              // well -> perf (scatter)
                Eigen::SparseMatrix<double> p2w;              // perf -> well (gather)
                Eigen::SparseMatrix<double> w2s;              // well -> segment (scatter)
                Eigen::SparseMatrix<double> s2w;              // segment -> well (gather)
                Eigen::SparseMatrix<double> s2p;              // segment -> perf (scatter)
                Eigen::SparseMatrix<double> p2s;              // perf -> segment (gather)
                Eigen::SparseMatrix<double> s2s_inlets;       // segment -> its inlet segments
                Eigen::SparseMatrix<double> s2s_outlet;       // segment -> its outlet segment
                Eigen::SparseMatrix<double> topseg2w;         // top segment -> well
                AutoDiffMatrix eliminate_topseg;              // change the top segment related to be zero
                std::vector<int> well_cells;                  // the set of perforated cells
                Vector conn_trans_factors;                         // connection transmissibility factors
                bool has_multisegment_wells;                  // flag indicating whether there is any muli-segment well
            };

            // copied from BlackoilModelBase
            // should put to somewhere better
            using DataBlock =  Eigen::Array<double,
                                            Eigen::Dynamic,
                                            Eigen::Dynamic,
                                            Eigen::RowMajor>;

            // ---------  Public methods  ---------
            // TODO: using a vector of WellMultiSegmentConstPtr for now
            // TODO: it should use const Wells or something else later.
            MultisegmentWells(const std::vector<WellMultiSegmentConstPtr>& wells_multisegment, const int np);

            const std::vector<WellMultiSegmentConstPtr>& wells() const;
            const MultisegmentWellOps& wellOps() const;

            int numSegment() const { return nseg_total_; };
            int numPerf() const { return nperf_total_; };

            const Vector& wellPerforationCellPressureDiffs() const { return well_perforation_cell_pressure_diffs_; };
            Vector& wellPerforationCellPressureDiffs() { return well_perforation_cell_pressure_diffs_; };

            const ADB& wellSegmentPerforationPressureDiffs() const { return well_segment_perforation_pressure_diffs_; };
            ADB& wellSegmentPerforationPressureDiffs() { return well_segment_perforation_pressure_diffs_; };

            const Vector& wellSegmentPerforationDepthDiffs() const { return well_segment_perforation_depth_diffs_; };
            Vector& wellSegmentPerforationDepthDiffs() { return well_segment_perforation_depth_diffs_; };

            const Vector& wellPerforationCellDensities() const { return well_perforation_cell_densities_; };
            Vector& wellPerforationCellDensities() { return well_perforation_cell_densities_; };

            const ADB& wellSegmentDensities() const { return well_segment_densities_; };
            ADB& wellSegmentDensities() { return well_segment_densities_; };

            const ADB& wellSegmentPressureDelta() const { return well_segment_pressures_delta_; };
            ADB& wellSegmentPressureDelta() { return well_segment_pressures_delta_; };

            const std::vector<Vector>& segmentCompSurfVolumeInitial() const { return segment_comp_surf_volume_initial_; };
            std::vector<Vector>& segmentCompSurfVolumeInitial() { return segment_comp_surf_volume_initial_; };

            const std::vector<ADB>& segmentCompSurfVolumeCurrent() const { return segment_comp_surf_volume_current_; };
            std::vector<ADB>& segmentCompSurfVolumeCurrent() { return segment_comp_surf_volume_current_; };

            const ADB& segmentMassFlowRates() const { return segment_mass_flow_rates_; };
            ADB& segmentMassFlowRates() { return segment_mass_flow_rates_; };

            const ADB& segmentViscosities() const { return segment_viscosities_; };
            ADB& segmentViscosities() { return segment_viscosities_; };

            const std::vector<int>& topWellSegments() const { return top_well_segments_; };
            std::vector<int>& topWellSegments() { return top_well_segments_; };

            const Vector& segVDt() const { return segvdt_; };
            Vector& segVDt() { return segvdt_; };



            template <class WellState>
            void
            updateWellState(const Vector& dwells,
                            const int np,
                            const double dpmaxrel,
                            WellState& well_state) const;

            // TODO: some arguments can be removed later
            // TODO: compi will be required in the multisegment wells
            template <class SolutionState>
            void
            computeWellFlux(const SolutionState& state,
                            const Opm::PhaseUsage& pu,
                            const std::vector<bool>& active,
                            const Vector& well_perforation_pressure_diffs,
                            const DataBlock& compi,
                            const std::vector<ADB>& mob_perfcells,
                            const std::vector<ADB>& b_perfcells,
                            const int np,
                            Vector& aliveWells,
                            std::vector<ADB>& cq_s) const;


            // Calculate the density of the mixture in the segments
            // And the surface volume of the components in the segments by dt
            template <class SolutionState>
            void
            computeSegmentFluidProperties(const SolutionState& state,
                                          const std::vector<PhasePresence>& pc,
                                          const std::vector<bool>& active,
                                          const BlackoilPropsAdInterface& fluid,
                                          const int np);

            void
            computeSegmentPressuresDelta(const double grav);


    protected:
        // TODO: probably a wells_active_ will be required here.
        const std::vector<WellMultiSegmentConstPtr> wells_multisegment_;
        const MultisegmentWellOps wops_ms_;
        int nseg_total_;
        int nperf_total_;

        // Pressure correction due to the different depth of the perforation
        // and the cell center of the grid block
        // For the non-segmented wells, since the perforation are forced to be
        // at the center of the grid cell, it should be ZERO.
        // It only applies to the mutli-segmented wells.
        Vector well_perforation_cell_pressure_diffs_;

        // Pressure correction due to the depth differennce between segment depth and perforation depth.
        ADB well_segment_perforation_pressure_diffs_;

        // The depth difference between segment nodes and perforations
        Vector well_segment_perforation_depth_diffs_;

        // the average of the fluid densities in the grid block
        // which is used to calculate the hydrostatic head correction due to the depth difference of the perforation
        // and the cell center of the grid block
        Vector well_perforation_cell_densities_;

        // the density of the fluid mixture in the segments
        // which is calculated in an implicit way
        ADB well_segment_densities_;

        // the hydrostatic pressure drop between segment nodes
        // calculated with the above density of fluid mixtures
        // for the top segment, they should always be zero for the moment.
        ADB well_segment_pressures_delta_;

        // the surface volume of components in the segments
        // the initial value at the beginning of the time step
        std::vector<Vector>   segment_comp_surf_volume_initial_;

        // the value within the current iteration.
        std::vector<ADB> segment_comp_surf_volume_current_;

        // the mass flow rate in the segments
        ADB segment_mass_flow_rates_;

        // the viscosity of the fluid mixture in the segments
        // TODO: it is only used to calculate the Reynolds number as we know
        //       maybe it is not better just to store the Reynolds number?
        ADB segment_viscosities_;

        std::vector<int> top_well_segments_;

        // segment volume by dt (time step)
        // to handle the volume effects of the segment
        Vector segvdt_;

    };

} // namespace Opm

#include "MultisegmentWells_impl.hpp"

#endif // OPM_MULTISEGMENTWELLS_HEADER_INCLUDED
