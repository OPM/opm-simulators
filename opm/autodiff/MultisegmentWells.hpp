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

#include <dune/common/parallel/mpihelper.hh>

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <cassert>

#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/wells/WellCollection.hpp>

#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/autodiff/BlackoilModelEnums.hpp>
#include <opm/autodiff/BlackoilPropsAdFromDeck.hpp>
#include <opm/autodiff/LinearisedBlackoilResidual.hpp>
#include <opm/autodiff/WellHelpers.hpp>
#include <opm/autodiff/VFPPropertiesAdb.hpp>

#include <opm/autodiff/WellMultiSegment.hpp>
#include <opm/autodiff/WellDensitySegmented.hpp>
#include "WellSwitchingLogger.hpp"



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
            using Communication =
                Dune::CollectiveCommunication<typename Dune::MPIHelper
                                              ::MPICommunicator>;

            // ---------  Public methods  ---------
            // TODO: using a vector of WellMultiSegmentConstPtr for now
            // TODO: it should use const Wells or something else later.
            MultisegmentWells(const Wells* wells_arg,
                              WellCollection* well_collection,
                              const std::vector< const Well* >& wells_ecl,
                              const int time_step);

            std::vector<WellMultiSegmentConstPtr> createMSWellVector(const Wells* wells_arg,
                                                                     const std::vector< const Well* >& wells_ecl,
                                                                     const int time_step);

            void init(const BlackoilPropsAdFromDeck* fluid_arg,
                      const std::vector<bool>* active_arg,
                      const std::vector<PhasePresence>* pc_arg,
                      const VFPPropertiesAdb*  vfp_properties_arg,
                      const double gravity_arg,
                      const Vector& depth_arg);

            const std::vector<WellMultiSegmentConstPtr>& msWells() const;
            const MultisegmentWellOps& wellOps() const;

            const Wells& wells() const;

            const Wells* wellsPointer() const;

            int numPhases() const { return num_phases_; };

            int numWells() const { return wells_multisegment_.size(); }

            int numSegment() const { return nseg_total_; };

            int numPerf() const { return nperf_total_; };

            bool wellsActive() const { return wells_active_; };

            void setWellsActive(const bool wells_active) { wells_active_ = wells_active; };

            bool localWellsActive() const { return ! wells_multisegment_.empty(); };

            int numWellVars() const { return (num_phases_ + 1) * nseg_total_; };

            template <class ReservoirResidualQuant, class SolutionState>
            void extractWellPerfProperties(const SolutionState& state,
                                           const std::vector<ReservoirResidualQuant>& rq,
                                           std::vector<ADB>& mob_perfcells,
                                           std::vector<ADB>& b_perfcells) const;

            Vector& wellPerforationCellPressureDiffs() { return well_perforation_cell_pressure_diffs_; };

            Vector& wellSegmentPerforationDepthDiffs() { return well_segment_perforation_depth_diffs_; };

            const Vector& wellPerforationCellDensities() const { return well_perforation_cell_densities_; };
            Vector& wellPerforationCellDensities() { return well_perforation_cell_densities_; };

            const std::vector<Vector>& segmentCompSurfVolumeInitial() const { return segment_comp_surf_volume_initial_; };
            std::vector<Vector>& segmentCompSurfVolumeInitial() { return segment_comp_surf_volume_initial_; };

            const std::vector<ADB>& segmentCompSurfVolumeCurrent() const { return segment_comp_surf_volume_current_; };

            const std::vector<int>& topWellSegments() const { return top_well_segments_; };
            std::vector<int>& topWellSegments() { return top_well_segments_; };

            Vector& segVDt() { return segvdt_; };

            const Vector& wellPerforationDensities() const { return well_perforation_densities_; };
            Vector& wellPerforationDensities() { return well_perforation_densities_; };

            const Vector& wellPerforationPressureDiffs() const { return well_perforation_pressure_diffs_; };
            Vector& wellPerforationPressureDiffs() { return well_perforation_pressure_diffs_; };

            template <class WellState>
            void
            updateWellState(const Vector& dwells,
                            const double dpmaxrel,
                            WellState& well_state) const;

            // TODO: some arguments can be removed later
            // TODO: compi will be required in the multisegment wells
            template <class SolutionState>
            void
            computeWellFlux(const SolutionState& state,
                            const std::vector<ADB>& mob_perfcells,
                            const std::vector<ADB>& b_perfcells,
                            Vector& aliveWells,
                            std::vector<ADB>& cq_s) const;

            template <class SolutionState, class WellState>
            void updatePerfPhaseRatesAndPressures(const std::vector<ADB>& cq_s,
                                                  const SolutionState& state,
                                                  WellState& xw) const;


            // Calculate the density of the mixture in the segments
            // And the surface volume of the components in the segments by dt
            template <class SolutionState>
            void
            computeSegmentFluidProperties(const SolutionState& state);

            void
            computeSegmentPressuresDelta(const double grav);

            template <class SolutionState>
            void
            addWellFluxEq(const std::vector<ADB>& cq_s,
                          const SolutionState& state,
                          LinearisedBlackoilResidual& residual);

            template <class SolutionState, class WellState>
            void
            addWellControlEq(const SolutionState& state,
                             const WellState& xw,
                             const Vector& aliveWells,
                             LinearisedBlackoilResidual& residual);

            template <class WellState>
            void
            updateWellControls(WellState& xw) const;

            // TODO: these code are same with the StandardWells
            // to find a better solution later.
            void
            variableStateWellIndices(std::vector<int>& indices,
                                     int& next) const;

            template <class SolutionState>
            void
            variableStateExtractWellsVars(const std::vector<int>& indices,
                                          std::vector<ADB>& vars,
                                          SolutionState& state) const;

            std::vector<int>
            variableWellStateIndices() const;

            template <class WellState>
            void
            variableWellStateInitials(const WellState& xw,
                                      std::vector<Vector>& vars0) const;

            template <class SolutionState, class WellState>
            void computeWellConnectionPressures(const SolutionState& state,
                                                const WellState& xw,
                                                const std::vector<ADB>& kr_adb,
                                                const std::vector<ADB>& fluid_density);

            WellCollection* wellCollection() const;

            void calculateEfficiencyFactors();

            const Vector& wellPerfEfficiencyFactors() const;

            // just return the passed well state
            template<class WellState>
            const WellState& wellState(const WellState& well_state) const { return well_state; }


    protected:
        // TODO: probably a wells_active_ will be required here.
        bool wells_active_;
        std::vector<WellMultiSegmentConstPtr> wells_multisegment_;
        MultisegmentWellOps wops_ms_;
        // It will probably need to be updated during running time.
        WellCollection* well_collection_;

        // The efficiency factor for each connection
        // It is specified based on wells and groups
        // By default, they should all be one.
        Vector well_perforation_efficiency_factors_;

        const int num_phases_;
        int nseg_total_;
        int nperf_total_;

        // TODO: put the Wells here to simplify the interface
        // TODO: at the moment, both wells_ and wells_multisegment_
        // TODO: acutally contain all the wells
        // TODO: they should be split eventually.
        const Wells* wells_;

        const BlackoilPropsAdFromDeck* fluid_;
        const std::vector<bool>*  active_;
        const std::vector<PhasePresence>*  phase_condition_;
        const VFPPropertiesAdb* vfp_properties_;
        double gravity_;
        // The depth of the all the cell centers
        // It is different from the perforation depth in MultisegmentWells
        Vector perf_cell_depth_;

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

        // The depth difference between the perforations and the perforation cells.
        Vector perf_cell_depth_diffs_;

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

        // technically, they are only useful for standard wells
        // since at the moment, we are handling both the standard
        // wells and the multi-segment wells under the MultisegmentWells
        // we need them to remove the dependency on StandardWells
        Vector well_perforation_densities_;
        Vector well_perforation_pressure_diffs_;

    };

} // namespace Opm

#include "MultisegmentWells_impl.hpp"

#endif // OPM_MULTISEGMENTWELLS_HEADER_INCLUDED
