/*
  Copyright 2013, 2015 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_BLACKOIMULTISEGMENTLMODEL_IMPL_HEADER_INCLUDED
#define OPM_BLACKOIMULTISEGMENTLMODEL_IMPL_HEADER_INCLUDED

#include <opm/autodiff/BlackoilMultiSegmentModel.hpp>

#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/autodiff/GridHelpers.hpp>
#include <opm/autodiff/BlackoilPropsAdInterface.hpp>
#include <opm/autodiff/GeoProps.hpp>
#include <opm/autodiff/WellDensitySegmented.hpp>
#include <opm/autodiff/VFPProperties.hpp>
#include <opm/autodiff/VFPProdProperties.hpp>
#include <opm/autodiff/VFPInjProperties.hpp>

#include <opm/core/grid.h>
#include <opm/core/linalg/LinearSolverInterface.hpp>
#include <opm/core/linalg/ParallelIstlInformation.hpp>
#include <opm/core/props/rock/RockCompressibility.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/well_controls.h>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <vector>
//#include <fstream>


namespace Opm {


    namespace detail
    {
        ADB onlyWellDerivs(const ADB& x)
        {
            V val = x.value();
            const int nb = x.numBlocks();
            if (nb < 2) {
                OPM_THROW(std::logic_error, "Called onlyWellDerivs() with argument that has " << nb << " blocks.");
            }
            std::vector<M> derivs = { x.derivative()[nb - 2], x.derivative()[nb - 1] };
            return ADB::function(std::move(val), std::move(derivs));
        }
    } // namespace detail




    template <class Grid>
    BlackoilMultiSegmentModel<Grid>::
    BlackoilMultiSegmentModel(const typename Base::ModelParameters&  param,
                  const Grid&                     grid ,
                  const BlackoilPropsAdInterface& fluid,
                  const DerivedGeology&           geo  ,
                  const RockCompressibility*      rock_comp_props,
                  const Wells*                    wells_arg,
                  const NewtonIterationBlackoilInterface&    linsolver,
                  Opm::EclipseStateConstPtr eclState,
                  const bool has_disgas,
                  const bool has_vapoil,
                  const bool terminal_output,
                  const std::vector<WellMultiSegmentConstPtr>& wells_multisegment)
        : Base(param, grid, fluid, geo, rock_comp_props, wells_arg, linsolver,
               eclState, has_disgas, has_vapoil, terminal_output)
        , well_segment_perforation_pressure_diffs_(ADB::null())
        , well_segment_densities_(ADB::null())
        , well_segment_pressures_delta_(ADB::null())
        , segment_comp_surf_volume_initial_(fluid.numPhases())
        , segment_comp_surf_volume_current_(fluid.numPhases(), ADB::null())
        , segment_mass_flow_rates_(ADB::null())
        , segment_viscosities_(ADB::null())
        , wells_multisegment_(wells_multisegment)
        , wops_ms_(wells_multisegment)
    {
    }






    template <class Grid>
    BlackoilMultiSegmentModel<Grid>::
    MultiSegmentWellOps::MultiSegmentWellOps(const std::vector<WellMultiSegmentConstPtr>& wells_ms)
    {
        if (wells_ms.empty()) {
            return;
        }

        // Count the total number of perforations and segments.
        const int nw = wells_ms.size();
        int total_nperf = 0;
        int total_nseg = 0;
        for (int w = 0; w < nw; ++w) {
            total_nperf += wells_ms[w]->numberOfPerforations();
            total_nseg += wells_ms[w]->numberOfSegments();
        }

        // Create well_cells and conn_trans_factors.
        well_cells.reserve(total_nperf);
        conn_trans_factors.resize(total_nperf);
        int well_perf_start = 0;
        for (int w = 0; w < nw; ++w) {
            WellMultiSegmentConstPtr well = wells_ms[w];
            well_cells.insert(well_cells.end(), well->wellCells().begin(), well->wellCells().end());
            const std::vector<double>& perf_trans = well->wellIndex();
            std::copy(perf_trans.begin(), perf_trans.end(), conn_trans_factors.data() + well_perf_start);
            well_perf_start += well->numberOfPerforations();
        }
        assert(well_perf_start == total_nperf);
        assert(int(well_cells.size()) == total_nperf);

        // Create all the operator matrices,
        // using the setFromTriplets() method.
        s2s_inlets = Eigen::SparseMatrix<double>(total_nseg, total_nseg);
        s2s_outlet = Eigen::SparseMatrix<double>(total_nseg, total_nseg);
        s2w        = Eigen::SparseMatrix<double>(nw, total_nseg);
        w2s        = Eigen::SparseMatrix<double>(total_nseg, nw);
        topseg2w   = Eigen::SparseMatrix<double>(nw, total_nseg);
        s2p = Eigen::SparseMatrix<double>(total_nperf, total_nseg);
        p2s = Eigen::SparseMatrix<double>(total_nseg, total_nperf);
        typedef Eigen::Triplet<double> Tri;
        std::vector<Tri> s2s_inlets_vector;
        std::vector<Tri> s2s_outlet_vector;
        std::vector<Tri> s2w_vector;
        std::vector<Tri> w2s_vector;
        std::vector<Tri> topseg2w_vector;
        std::vector<Tri> s2p_vector;
        std::vector<Tri> p2s_vector;
        V topseg_zero = V::Ones(total_nseg);
        s2s_inlets_vector.reserve(total_nseg);
        s2s_outlet_vector.reserve(total_nseg);
        s2w_vector.reserve(total_nseg);
        w2s_vector.reserve(total_nseg);
        topseg2w_vector.reserve(nw);
        s2p_vector.reserve(total_nperf);
        p2s_vector.reserve(total_nperf);
        int seg_start = 0;
        int perf_start = 0;
        for (int w = 0; w < nw; ++w) {
            const int ns = wells_ms[w]->numberOfSegments();
            const int np = wells_ms[w]->numberOfPerforations();
            for (int seg = 0; seg < ns; ++seg) {
                const int seg_ind = seg_start + seg;
                w2s_vector.push_back(Tri(seg_ind, w, 1.0));
                s2w_vector.push_back(Tri(w, seg_ind, 1.0));
                if (seg == 0) {
                    topseg2w_vector.push_back(Tri(w, seg_ind, 1.0));
                    topseg_zero(seg_ind) = 0.0;
                }
                int seg_outlet = wells_ms[w]->outletSegment()[seg];
                if (seg_outlet >= 0) {
                    const int outlet_ind = seg_start + seg_outlet;
                    s2s_inlets_vector.push_back(Tri(outlet_ind, seg_ind, 1.0));
                    s2s_outlet_vector.push_back(Tri(seg_ind, outlet_ind, 1.0));
                }

                const auto& seg_perf = wells_ms[w]->segmentPerforations()[seg];
                // the number of perforations related to this segment
                const int npseg = seg_perf.size();
                for (int perf = 0; perf < npseg; ++perf) {
                    const int perf_ind = perf_start + seg_perf[perf];
                    s2p_vector.push_back(Tri(perf_ind, seg_ind, 1.0));
                    p2s_vector.push_back(Tri(seg_ind, perf_ind, 1.0));
                }
            }
            seg_start += ns;
            perf_start += np;
        }

        s2s_inlets.setFromTriplets(s2s_inlets_vector.begin(), s2s_inlets_vector.end());
        s2s_outlet.setFromTriplets(s2s_outlet_vector.begin(), s2s_outlet_vector.end());
        w2s.setFromTriplets(w2s_vector.begin(), w2s_vector.end());
        s2w.setFromTriplets(s2w_vector.begin(), s2w_vector.end());
        topseg2w.setFromTriplets(topseg2w_vector.begin(), topseg2w_vector.end());
        s2p.setFromTriplets(s2p_vector.begin(), s2p_vector.end());
        p2s.setFromTriplets(p2s_vector.begin(), p2s_vector.end());

        w2p = Eigen::SparseMatrix<double>(total_nperf, nw);
        p2w = Eigen::SparseMatrix<double>(nw, total_nperf);
        w2p = s2p * w2s;
        p2w = s2w * p2s;

        eliminate_topseg = AutoDiffMatrix(topseg_zero.matrix().asDiagonal());
    }






    template <class Grid>
    void
    BlackoilMultiSegmentModel<Grid>::
    prepareStep(const double dt,
                ReservoirState& reservoir_state,
                WellState& well_state)
    {
        pvdt_ = geo_.poreVolume() / dt;
        if (active_[Gas]) {
            updatePrimalVariableFromState(reservoir_state);
        }

        top_well_segments_ = well_state.topSegmentLoc();

        const int nw = wellsMultiSegment().size();
        const int nseg_total = well_state.numSegments();
        std::vector<double> segment_volume;
        segment_volume.reserve(nseg_total);
        for (int w = 0; w < nw; ++w) {
            WellMultiSegmentConstPtr well = wellsMultiSegment()[w];
            const std::vector<double>& segment_volume_well = well->segmentVolume();
            segment_volume.insert(segment_volume.end(), segment_volume_well.begin(), segment_volume_well.end());
        }
        assert(int(segment_volume.size()) == nseg_total);
        segvdt_ = Eigen::Map<V>(segment_volume.data(), nseg_total) / dt;
    }





    template <class Grid>
    int
    BlackoilMultiSegmentModel<Grid>::numWellVars() const
    {
        // For each segment, we have a pressure variable, and one flux per phase.
        const int nseg = wops_ms_.p2s.rows();
        return (numPhases() + 1) * nseg;
    }





    template <class Grid>
    void
    BlackoilMultiSegmentModel<Grid>::makeConstantState(SolutionState& state) const
    {
        Base::makeConstantState(state);
        state.segp  = ADB::constant(state.segp.value());
        state.segqs = ADB::constant(state.segqs.value());
    }





    template <class Grid>
    void
    BlackoilMultiSegmentModel<Grid>::variableWellStateInitials(const WellState& xw, std::vector<V>& vars0) const
    {
        // Initial well rates
        if ( wellsMultiSegment().size() > 0 )
        {
            // Need to reshuffle well segment rates, from phase running fastest
            const int nseg = xw.numSegments();
            const int np = xw.numPhases();

            // The transpose() below switches the ordering of the segment rates
            const DataBlock segrates = Eigen::Map<const DataBlock>(& xw.segPhaseRates()[0], nseg, np).transpose();
            // segment phase rates in surface volume
            const V segqs = Eigen::Map<const V>(segrates.data(), nseg * np);
            vars0.push_back(segqs);

            // for the pressure of the segments
            const V segp = Eigen::Map<const V>(& xw.segPress()[0], xw.segPress().size());
            vars0.push_back(segp);
        }
        else
        {
            // push null sates for segqs and segp
            vars0.push_back(V());
            vars0.push_back(V());
        }
    }





    template <class Grid>
    void
    BlackoilMultiSegmentModel<Grid>::variableStateExtractWellsVars(const std::vector<int>& indices,
                                                                   std::vector<ADB>& vars,
                                                                   SolutionState& state) const
    {
        // TODO: using the original Qs for the segment rates for now, to be fixed eventually.
        // TODO: using the original Bhp for the segment pressures for now, to be fixed eventually.

        // segment phase rates in surface volume
        state.segqs = std::move(vars[indices[Qs]]);

        // segment pressures
        state.segp = std::move(vars[indices[Bhp]]);

        // The qs and bhp are no longer primary variables, but could
        // still be used in computations. They are identical to the
        // pressures and flows of the top segments.
        const int np = numPhases();
        const int ns = state.segp.size();
        const int nw = top_well_segments_.size();
        state.qs = ADB::constant(ADB::V::Zero(np*nw));
        for (int phase = 0; phase < np; ++phase) {
            // Extract segment fluxes for this phase (ns consecutive elements).
            ADB segqs_phase = subset(state.segqs, Span(ns, 1, ns*phase));
            // Extract top segment fluxes (= well fluxes)
            ADB wellqs_phase = subset(segqs_phase, top_well_segments_);
            // Expand to full size of qs (which contains all phases) and add.
            state.qs += superset(wellqs_phase, Span(nw, 1, nw*phase), nw*np);
        }
        state.bhp = subset(state.segp, top_well_segments_);
    }




    // TODO: This is just a preliminary version, remains to be improved later when we decide a better way
    // TODO: to intergrate the usual wells and multi-segment wells.
    template <class Grid>
    void BlackoilMultiSegmentModel<Grid>::computeWellConnectionPressures(const SolutionState& state,
                                                                         const WellState& xw)
    {
        if( ! wellsActive() ) return ;

        using namespace Opm::AutoDiffGrid;
        // 1. Compute properties required by computeConnectionPressureDelta().
        //    Note that some of the complexity of this part is due to the function
        //    taking std::vector<double> arguments, and not Eigen objects.
        const int nperf_total = xw.numPerforations();
        const int nw = xw.numWells();

        // the well cells for multisegment wells and non-segmented wells should be counted seperatedly
        // indexing should be put in WellState
        std::vector<int> well_cells;
        std::vector<int> well_cells_segmented_idx;
        std::vector<int> well_cells_non_segmented_idx;
        std::vector<int> well_cells_segmented;
        std::vector<int> well_cells_non_segmented;

        well_cells_segmented_idx.reserve(nperf_total);
        well_cells_non_segmented_idx.reserve(nperf_total);
        well_cells_segmented.reserve(nperf_total);
        well_cells_non_segmented.reserve(nperf_total);
        well_cells.reserve(nperf_total);
        for (int i = 0; i < nw; ++i) {
            const std::vector<int>& temp_well_cells = wellsMultiSegment()[i]->wellCells();
            const int num_cells_this_well = temp_well_cells.size();
            const int n_current = well_cells.size();
            if (wellsMultiSegment()[i]->isMultiSegmented()) {
                for (int j = 0; j < num_cells_this_well; ++j) {
                    well_cells_segmented_idx.push_back(j + n_current);
                }
                well_cells_segmented.insert(well_cells_segmented.end(), temp_well_cells.begin(), temp_well_cells.end());
            } else {
                for (int j = 0; j < num_cells_this_well; ++j) {
                    well_cells_non_segmented_idx.push_back(j + n_current);
                }
                well_cells_non_segmented.insert(well_cells_non_segmented.end(), temp_well_cells.begin(), temp_well_cells.end());
            }
            well_cells.insert(well_cells.end(), temp_well_cells.begin(), temp_well_cells.end());
        }

        well_perforation_densities_ = V::Zero(nperf_total);

        // Compute the average pressure in each well block
        // The following code is necessary for non-segmented wells.
        // For the multi-segmented wells.
        // Two hydrostatic heads need to be computed.
        // One is the hydrostatic head between the cell and the completion depth
        // The density is the density mixture of that grid block
        // one is the hydrostatic head between the segment and the completion depth
        // The density is the density of the fluid mixture in the related segment

        // TODO: the following part should be modified to fit the plan for only the non-segmented wells
        // const int nperf_nonsegmented = well_cells_non_segmented.size();

        const V perf_press = Eigen::Map<const V>(xw.perfPress().data(), nperf_total);
        // const V perf_press_nonsegmented = subset(perf_press, well_cells_non_segmented_idx);

        V avg_press = perf_press * 0.0;

        // for the non-segmented wells, calculated the average pressures.
        // If it is the top perforation, then average with the bhp().
        // If it is not the top perforation, then average with the perforation above it().
        // TODO: Make sure the order of the perforation are not changed, comparing with the old wells structure.
        int start_segment = 0;
        for (int w = 0; w < nw; ++w) {
            const int nseg = wellsMultiSegment()[w]->numberOfSegments();
            if (wellsMultiSegment()[w]->isMultiSegmented()) {
                // maybe we should give some reasonable values to prevent the following calculations fail
                start_segment += nseg;
                continue;
            }

            std::string well_name(wellsMultiSegment()[w]->name());
            typedef typename WellStateMultiSegment::SegmentedWellMapType::const_iterator const_iterator;
            const_iterator it_well = xw.segmentedWellMap().find(well_name);
            assert(it_well != xw.segmentedWellMap().end());

            const int start_perforation = (*it_well).second.start_perforation;
            const int end_perforation = start_perforation + (*it_well).second.number_of_perforations;
            for (int perf = start_perforation; perf < end_perforation; ++perf) {
                const double p_above = perf == start_perforation ? state.segp.value()[start_segment] : perf_press[perf - 1];
                const double p_avg = (perf_press[perf] + p_above)/2;
                avg_press[perf] = p_avg;
            }
            start_segment += nseg;
        }
        assert(start_segment == xw.numSegments());

        // Use cell values for the temperature as the wells don't knows its temperature yet.
        const ADB perf_temp = subset(state.temperature, well_cells);

        // Compute b, rsmax, rvmax values for perforations.
        // Evaluate the properties using average well block pressures
        // and cell values for rs, rv, phase condition and temperature.
        const ADB avg_press_ad = ADB::constant(avg_press);
        std::vector<PhasePresence> perf_cond(nperf_total);
        const std::vector<PhasePresence>& pc = phaseCondition();
        for (int perf = 0; perf < nperf_total; ++perf) {
            perf_cond[perf] = pc[well_cells[perf]];
        }
        const PhaseUsage& pu = fluid_.phaseUsage();
        DataBlock b(nperf_total, pu.num_phases);
        std::vector<double> rsmax_perf(nperf_total, 0.0);
        std::vector<double> rvmax_perf(nperf_total, 0.0);
        if (pu.phase_used[BlackoilPhases::Aqua]) {
            const V bw = fluid_.bWat(avg_press_ad, perf_temp, well_cells).value();
            b.col(pu.phase_pos[BlackoilPhases::Aqua]) = bw;
        }
        assert(active_[Oil]);
        const V perf_so =  subset(state.saturation[pu.phase_pos[Oil]].value(), well_cells);
        if (pu.phase_used[BlackoilPhases::Liquid]) {
            const ADB perf_rs = subset(state.rs, well_cells);
            const V bo = fluid_.bOil(avg_press_ad, perf_temp, perf_rs, perf_cond, well_cells).value();
            b.col(pu.phase_pos[BlackoilPhases::Liquid]) = bo;
            const V rssat = fluidRsSat(avg_press, perf_so, well_cells);
            rsmax_perf.assign(rssat.data(), rssat.data() + nperf_total);
        }
        if (pu.phase_used[BlackoilPhases::Vapour]) {
            const ADB perf_rv = subset(state.rv, well_cells);
            const V bg = fluid_.bGas(avg_press_ad, perf_temp, perf_rv, perf_cond, well_cells).value();
            b.col(pu.phase_pos[BlackoilPhases::Vapour]) = bg;
            const V rvsat = fluidRvSat(avg_press, perf_so, well_cells);
            rvmax_perf.assign(rvsat.data(), rvsat.data() + nperf_total);
        }
        // b is row major, so can just copy data.
        std::vector<double> b_perf(b.data(), b.data() + nperf_total * pu.num_phases);
        // Extract well connection depths.
        const V depth = cellCentroidsZToEigen(grid_);
        const V perfcelldepth = subset(depth, well_cells);
        std::vector<double> perf_cell_depth(perfcelldepth.data(), perfcelldepth.data() + nperf_total);

        // Surface density.
        // The compute density segment wants the surface densities as
        // an np * number of wells cells array
        V rho = superset(fluid_.surfaceDensity(0 , well_cells), Span(nperf_total, pu.num_phases, 0), nperf_total * pu.num_phases);
        for (int phase = 1; phase < pu.num_phases; ++phase) {
            rho += superset(fluid_.surfaceDensity(phase , well_cells), Span(nperf_total, pu.num_phases, phase), nperf_total * pu.num_phases);
        }
        std::vector<double> surf_dens_perf(rho.data(), rho.data() + nperf_total * pu.num_phases);

        // Gravity
        double grav = detail::getGravity(geo_.gravity(), dimensions(grid_));

        // 2. Compute densities
        std::vector<double> cd =
                WellDensitySegmented::computeConnectionDensities(
                        wells(), xw, fluid_.phaseUsage(),
                        b_perf, rsmax_perf, rvmax_perf, surf_dens_perf);

        // 3. Compute pressure deltas
        std::vector<double> cdp =
                WellDensitySegmented::computeConnectionPressureDelta(
                        wells(), perf_cell_depth, cd, grav);

        // 4. Store the results
        well_perforation_densities_ = Eigen::Map<const V>(cd.data(), nperf_total); // This one is not useful for segmented wells at all
        well_perforation_pressure_diffs_ = Eigen::Map<const V>(cdp.data(), nperf_total);

        // compute the average of the fluid densites in the well blocks.
        // the average is weighted according to the fluid relative permeabilities.
        const std::vector<ADB> kr_adb = Base::computeRelPerm(state);
        size_t temp_size = kr_adb.size();
        std::vector<V> perf_kr;
        for(size_t i = 0; i < temp_size; ++i) {
            // const ADB kr_phase_adb = subset(kr_adb[i], well_cells);
            const V kr_phase = (subset(kr_adb[i], well_cells)).value();
            perf_kr.push_back(kr_phase);
        }


        // compute the averaged density for the well block
        // TODO: for the non-segmented wells, they should be set to zero
        // TODO: for the moment, they are still calculated, while not used later.
        for (int i = 0; i < nperf_total; ++i) {
            double sum_kr = 0.;
            int np = perf_kr.size(); // make sure it is 3
            for (int p = 0;  p < np; ++p) {
                sum_kr += perf_kr[p][i];
            }

            for (int p = 0; p < np; ++p) {
                perf_kr[p][i] /= sum_kr;
            }
        }

        V rho_avg_perf = V::Constant(nperf_total, 0.0);
        // TODO: make sure the order of the density and the order of the kr are the same.
        for (int phaseIdx = 0; phaseIdx < fluid_.numPhases(); ++phaseIdx) {
            const int canonicalPhaseIdx = canph_[phaseIdx];
            const ADB fluid_density = fluidDensity(canonicalPhaseIdx, rq_[phaseIdx].b, state.rs, state.rv);
            const V rho_perf = subset(fluid_density, well_cells).value();
            // TODO: phaseIdx or canonicalPhaseIdx ?
            rho_avg_perf += rho_perf * perf_kr[phaseIdx];
        }

        // TODO: rho_avg_perf is actually the well_perforation_cell_densities_;
        well_perforation_cell_densities_ = Eigen::Map<const V>(rho_avg_perf.data(), nperf_total);

        // We should put this in a global class
        std::vector<double> perf_depth_vec;
        perf_depth_vec.reserve(nperf_total);
        for (int w = 0; w < nw; ++w) {
            WellMultiSegmentConstPtr well = wellsMultiSegment()[w];
            const std::vector<double>& perf_depth_well = well->perfDepth();
            perf_depth_vec.insert(perf_depth_vec.end(), perf_depth_well.begin(), perf_depth_well.end());
        }
        assert(int(perf_depth_vec.size()) == nperf_total);
        const V perf_depth = Eigen::Map<V>(perf_depth_vec.data(), nperf_total);

        const V perf_cell_depth_diffs = perf_depth - perfcelldepth;

        well_perforation_cell_pressure_diffs_ = grav * well_perforation_cell_densities_ * perf_cell_depth_diffs;


        // Calculating the depth difference between segment nodes and perforations.
        // TODO: should be put somewhere else for better clarity later
        well_segment_perforation_depth_diffs_ = V::Constant(nperf_total, -1e100);

        int start_perforation = 0;
        for (int w = 0; w < nw; ++w) {
            WellMultiSegmentConstPtr well = wellsMultiSegment()[w];
            const int nseg = well->numberOfSegments();
            const int nperf = well->numberOfPerforations();
            const std::vector<std::vector<int>>& segment_perforations = well->segmentPerforations();
            for (int s = 0; s < nseg; ++s) {
                const int nperf_seg = segment_perforations[s].size();
                const double segment_depth = well->segmentDepth()[s];
                for (int perf = 0; perf < nperf_seg; ++perf) {
                    const int perf_number = segment_perforations[s][perf] + start_perforation;
                    well_segment_perforation_depth_diffs_[perf_number] = segment_depth - perf_depth[perf_number];
                }
            }
            start_perforation += nperf;
        }
        assert(start_perforation == nperf_total);
    }







    template <class Grid>
    void
    BlackoilMultiSegmentModel<Grid>::
    assemble(const ReservoirState& reservoir_state,
             WellState& well_state,
             const bool initial_assembly)
    {
        using namespace Opm::AutoDiffGrid;

        // TODO: include VFP effect.
        // If we have VFP tables, we need the well connection
        // pressures for the "simple" hydrostatic correction
        // between well depth and vfp table depth.
        //  if (isVFPActive()) {
        //     SolutionState state = asImpl().variableState(reservoir_state, well_state);
        //     SolutionState state0 = state;
        //     asImpl().makeConstantState(state0);
        //     asImpl().computeWellConnectionPressures(state0, well_state);
        // }

        // Possibly switch well controls and updating well state to
        // get reasonable initial conditions for the wells
        asImpl().updateWellControls(well_state);

        // Create the primary variables.
        SolutionState state = asImpl().variableState(reservoir_state, well_state);

        if (initial_assembly) {
            // Create the (constant, derivativeless) initial state.
            SolutionState state0 = state;
            asImpl().makeConstantState(state0);
            // Compute initial accumulation contributions
            // and well connection pressures.
            asImpl().computeAccum(state0, 0);
            asImpl().computeSegmentFluidProperties(state0);
            const int np = numPhases();
            assert(np == int(segment_comp_surf_volume_initial_.size()));
            for (int phase = 0; phase < np; ++phase) {
                segment_comp_surf_volume_initial_[phase] = segment_comp_surf_volume_current_[phase].value();
            }
            asImpl().computeWellConnectionPressures(state0, well_state);
        }

        // OPM_AD_DISKVAL(state.pressure);
        // OPM_AD_DISKVAL(state.saturation[0]);
        // OPM_AD_DISKVAL(state.saturation[1]);
        // OPM_AD_DISKVAL(state.saturation[2]);
        // OPM_AD_DISKVAL(state.rs);
        // OPM_AD_DISKVAL(state.rv);
        // OPM_AD_DISKVAL(state.qs);
        // OPM_AD_DISKVAL(state.bhp);

        // -------- Mass balance equations --------
        asImpl().assembleMassBalanceEq(state);

        // -------- Well equations ----------

        if ( ! wellsActive() ) {
            return;
        }

        asImpl().computeSegmentFluidProperties(state);
        asImpl().computeSegmentPressuresDelta(state);

        std::vector<ADB> mob_perfcells;
        std::vector<ADB> b_perfcells;
        asImpl().extractWellPerfProperties(mob_perfcells, b_perfcells);
        if (param_.solve_welleq_initially_ && initial_assembly) {
            // solve the well equations as a pre-processing step
            asImpl().solveWellEq(mob_perfcells, b_perfcells, state, well_state);
        }

        // the perforation flux here are different
        // it is related to the segment location
        V aliveWells;
        std::vector<ADB> cq_s;
        asImpl().computeWellFlux(state, mob_perfcells, b_perfcells, aliveWells, cq_s);
        asImpl().updatePerfPhaseRatesAndPressures(cq_s, state, well_state);
        asImpl().addWellFluxEq(cq_s, state);
        asImpl().addWellContributionToMassBalanceEq(cq_s, state, well_state);
        asImpl().addWellControlEq(state, well_state, aliveWells);
    }





    template <class Grid>
    void
    BlackoilMultiSegmentModel<Grid>::computeWellFlux(const SolutionState& state,
                                                     const std::vector<ADB>& mob_perfcells,
                                                     const std::vector<ADB>& b_perfcells,
                                                     V& aliveWells,
                                                     std::vector<ADB>& cq_s) const
    {
        // if( ! wellsActive() ) return ;
        if (wellsMultiSegment().size() == 0) return;

        const int nw = wellsMultiSegment().size();
        const Opm::PhaseUsage& pu = fluid_.phaseUsage();

        aliveWells = V::Constant(nw, 1.0);

        const int np = numPhases();
        const int nseg = wops_ms_.s2p.cols();
        const int nperf = wops_ms_.s2p.rows();

        cq_s.resize(np, ADB::null());

        {
            const V& Tw = wops_ms_.conn_trans_factors;
            const std::vector<int>& well_cells = wops_ms_.well_cells;

            // determining in-flow (towards well-bore) or out-flow (towards reservoir)
            // for mutli-segmented wells and non-segmented wells, the calculation of the drawdown are different.
            const ADB& p_perfcells = subset(state.pressure, well_cells);
            const ADB& rs_perfcells = subset(state.rs, well_cells);
            const ADB& rv_perfcells = subset(state.rv, well_cells);

            const ADB& seg_pressures = state.segp;

            const ADB seg_pressures_perf = wops_ms_.s2p * seg_pressures;

            // Create selector for perforations of multi-segment vs. regular wells.
            V is_multisegment_well(nw);
            for (int w = 0; w < nw; ++w) {
                is_multisegment_well[w] = double(wellsMultiSegment()[w]->isMultiSegmented());
            }
            // Take one flag per well and expand to one flag per perforation.
            V is_multisegment_perf = wops_ms_.w2p * is_multisegment_well.matrix();
            Selector<double> msperf_selector(is_multisegment_perf, Selector<double>::NotEqualZero);

            // Compute drawdown.
            ADB h_nc = msperf_selector.select(well_segment_perforation_pressure_diffs_,
                                                    ADB::constant(well_perforation_pressure_diffs_));
            const V h_cj = msperf_selector.select(well_perforation_cell_pressure_diffs_, V::Zero(nperf));

            // Special handling for when we are called from solveWellEq().
            // TODO: restructure to eliminate need for special treatmemt.
            if ((h_nc.numBlocks() != 0) && (h_nc.numBlocks() != seg_pressures_perf.numBlocks())) {
                assert(seg_pressures_perf.numBlocks() == 2);
                assert(h_nc.numBlocks() > 2);
                h_nc = detail::onlyWellDerivs(h_nc);
                assert(h_nc.numBlocks() == 2);
            }

            ADB drawdown = (p_perfcells + h_cj - seg_pressures_perf - h_nc);

            // selects injection perforations
            V selectInjectingPerforations = V::Zero(nperf);
            // selects producing perforations
            V selectProducingPerforations = V::Zero(nperf);
            for (int c = 0; c < nperf; ++c){
                if (drawdown.value()[c] < 0)
                    selectInjectingPerforations[c] = 1;
                else
                    selectProducingPerforations[c] = 1;
            }

            // handling flow into wellbore
            // maybe there are something to do there make the procedure easier.
            std::vector<ADB> cq_ps(np, ADB::null());
            for (int phase = 0; phase < np; ++phase) {
                const ADB cq_p = -(selectProducingPerforations * Tw) * (mob_perfcells[phase] * drawdown);
                cq_ps[phase] = b_perfcells[phase] * cq_p;
            }

            if (active_[Oil] && active_[Gas]) {
                const int oilpos = pu.phase_pos[Oil];
                const int gaspos = pu.phase_pos[Gas];
                const ADB cq_psOil = cq_ps[oilpos];
                const ADB cq_psGas = cq_ps[gaspos];
                cq_ps[gaspos] += rs_perfcells * cq_psOil;
                cq_ps[oilpos] += rv_perfcells * cq_psGas;
            }

            // hadling flow out from wellbore
            ADB total_mob = mob_perfcells[0];
            for (int phase = 1; phase < np; ++phase) {
                total_mob += mob_perfcells[phase];
            }

            // injection perforations total volume rates
            const ADB cqt_i = -(selectInjectingPerforations * Tw) * (total_mob * drawdown);

            // compute wellbore mixture for injecting perforations
            // The wellbore mixture depends on the inflow from the reservoir
            // and the well injection rates.
            // TODO: should this based on the segments?
            // TODO: for the usual wells, the well rates are the sum of the perforations.
            // TODO: for multi-segmented wells, the segment rates are not the sum of the perforations.

            // TODO: two options here
            // TODO: 1. for each segment, only the inflow from the perforations related to this segment are considered.
            // TODO: 2. for each segment, the inflow from the perforrations related to this segment and also all the inflow
            // TODO: from the upstreaming sgments and their perforations need to be considered.
            // TODO: This way can be the more consistent way, while let us begin with the first option. The second option
            // TODO: involves one operations that are not valid now. (i.e. how to transverse from the leaves to the root,
            // TODO: although we can begin from the brutal force way)

            // TODO: stop using wells() here.
            const DataBlock compi = Eigen::Map<const DataBlock>(wells().comp_frac, nw, np);
            std::vector<ADB> wbq(np, ADB::null());
            ADB wbqt = ADB::constant(V::Zero(nseg));

            for (int phase = 0; phase < np; ++phase) {
                const ADB& q_ps = wops_ms_.p2s * cq_ps[phase];
                const ADB& q_s = subset(state.segqs, Span(nseg, 1, phase * nseg));
                Selector<double> injectingPhase_selector(q_s.value(), Selector<double>::GreaterZero);

                const int pos = pu.phase_pos[phase];

                // this is per segment
                wbq[phase] = (wops_ms_.w2s * ADB::constant(compi.col(pos)) * injectingPhase_selector.select(q_s, ADB::constant(V::Zero(nseg)))) - q_ps;

                // TODO: it should be a single value for this certain well.
                // TODO: it need to be changed later to handle things more consistently
                // or there should be an earsier way to decide if the well is dead.
                wbqt += wbq[phase];
            }

            // Set aliveWells.
            // the first value of the wbqt is the one to decide if the well is dead
            // or there should be some dead segments?
            {
                int topseg = 0;
                for (int w = 0; w < nw; ++w) {
                    if (wbqt.value()[topseg] == 0.0) { // yes we really mean == here, no fuzzyness
                        aliveWells[w] = 0.0;
                    }
                    topseg += wells_multisegment_[w]->numberOfSegments();
                }
            }

            // compute wellbore mixture at standard conditions.
            // before, the determination of alive wells is based on wells.
            // now, will there be any dead segment? I think no.
            // TODO: it is not clear if the cmix_s should be based on segment or the well
            std::vector<ADB> cmix_s(np, ADB::null());
            Selector<double> aliveWells_selector(aliveWells, Selector<double>::NotEqualZero);
            for (int phase = 0; phase < np; ++phase) {
                const int pos = pu.phase_pos[phase];
                const ADB phase_fraction = wops_ms_.topseg2w * (wbq[phase] / wbqt);
                cmix_s[phase] = wops_ms_.w2p * aliveWells_selector.select(phase_fraction, ADB::constant(compi.col(pos)));
            }

            // compute volume ration between connection at standard conditions
            ADB volumeRatio = ADB::constant(V::Zero(nperf));
            const ADB d = V::Constant(nperf,1.0) -  rv_perfcells * rs_perfcells;

            for (int phase = 0; phase < np; ++phase) {
                ADB tmp = cmix_s[phase];
                if (phase == Oil && active_[Gas]) {
                    const int gaspos = pu.phase_pos[Gas];
                    tmp = tmp - rv_perfcells * cmix_s[gaspos] / d;
                }
                if (phase == Gas && active_[Oil]) {
                    const int oilpos = pu.phase_pos[Oil];
                    tmp = tmp - rs_perfcells * cmix_s[oilpos] / d;
                }
                volumeRatio += tmp / b_perfcells[phase];
            }

            // injecting connections total volumerates at standard conditions
            ADB cqt_is = cqt_i/volumeRatio;

            // connection phase volumerates at standard conditions
            for (int phase = 0; phase < np; ++phase) {
                cq_s[phase] = cq_ps[phase] + cmix_s[phase]*cqt_is;
            }
        }
    }





    template <class Grid>
    void BlackoilMultiSegmentModel<Grid>::updatePerfPhaseRatesAndPressures(const std::vector<ADB>& cq_s,
                                                                           const SolutionState& state,
                                                                           WellState& xw) const
    {
        // Update the perforation phase rates (used to calculate the pressure drop in the wellbore).
        const int np = numPhases();
        const int nw = wellsMultiSegment().size();
        const int nperf_total = xw.perfPress().size();

        V cq = superset(cq_s[0].value(), Span(nperf_total, np, 0), nperf_total * np);
        for (int phase = 1; phase < np; ++phase) {
            cq += superset(cq_s[phase].value(), Span(nperf_total, np, phase), nperf_total * np);
        }
        xw.perfPhaseRates().assign(cq.data(), cq.data() + nperf_total * np);

        // Update the perforation pressures for usual wells first to recover the resutls
        // without mutlti segment wells. For segment wells, it has not been decided if
        // we need th concept of preforation pressures
        xw.perfPress().resize(nperf_total, -1.e100);

        const V& cdp = well_perforation_pressure_diffs_;
        int start_segment = 0;
        int start_perforation = 0;
        for (int i = 0; i < nw; ++i) {
            WellMultiSegmentConstPtr well = wellsMultiSegment()[i];
            const int nperf = well->numberOfPerforations();
            const int nseg = well->numberOfSegments();
            if (well->isMultiSegmented()) {
                start_segment += nseg;
                start_perforation += nperf;
                continue;
            }
            const V cdp_well = subset(cdp, Span(nperf, 1, start_perforation));
            const ADB segp = subset(state.segp, Span(nseg, 1, start_segment));
            const V perfpressure = (well->wellOps().s2p * segp.value().matrix()).array() + cdp_well;
            std::copy(perfpressure.data(), perfpressure.data() + nperf, &xw.perfPress()[start_perforation]);

            start_segment += nseg;
            start_perforation += nperf;
        }
    }



    template <class Grid>
    void BlackoilMultiSegmentModel<Grid>::addWellFluxEq(const std::vector<ADB>& cq_s,
                                                        const SolutionState& state)
    {
        // the well flux equations are for each segment and each phase.
        //    /delta m_p_n / dt  - /sigma Q_pi - /sigma q_pj + Q_pn = 0
        // 1. It is the gain of the amount of the component p in the segment n during the
        //    current time step under stock-tank conditions.
        //    It is used to handle the volume storage effects of the wellbore.
        //    We need the information from the previous step and the crrent time step.
        // 2. for the second term, it is flow into the segment from the inlet segments,
        //    which are unknown and treated implictly.
        // 3. for the third term, it is the inflow through the perforations.
        // 4. for the last term, it is the outlet rates and also the segment rates,
        //    which are the primary variable.
        const int np = numPhases();
        const int nseg_total = state.segp.size();

        ADB segqs = state.segqs;

        std::vector<ADB> segment_volume_change_dt(np, ADB::null());
        for (int phase = 0; phase < np; ++phase) {
            // Gain of the surface volume of each component in the segment by dt
            segment_volume_change_dt[phase] = segment_comp_surf_volume_current_[phase] -
                                              segment_comp_surf_volume_initial_[phase];

            // Special handling for when we are called from solveWellEq().
            // TODO: restructure to eliminate need for special treatmemt.
            if (segment_volume_change_dt[phase].numBlocks() != segqs.numBlocks()) {
                assert(segment_volume_change_dt[phase].numBlocks() > 2);
                assert(segqs.numBlocks() == 2);
                segment_volume_change_dt[phase] = detail::onlyWellDerivs(segment_volume_change_dt[phase]);
                assert(segment_volume_change_dt[phase].numBlocks() == 2);
            }

            const ADB cq_s_seg = wops_ms_.p2s * cq_s[phase];
            const ADB segqs_phase = subset(segqs, Span(nseg_total, 1, phase * nseg_total));
            segqs -= superset(cq_s_seg + wops_ms_.s2s_inlets * segqs_phase + segment_volume_change_dt[phase],
                              Span(nseg_total, 1, phase * nseg_total), np * nseg_total);
        }

        residual_.well_flux_eq = segqs;
    }



    template <class Grid>
    void BlackoilMultiSegmentModel<Grid>::updateWellControls(WellState& xw) const
    {
        if( ! wellsActive() ) return ;

        std::string modestring[4] = { "BHP", "THP", "RESERVOIR_RATE", "SURFACE_RATE" };
        // Find, for each well, if any constraints are broken. If so,
        // switch control to first broken constraint.
        const int np = wellsMultiSegment()[0]->numberOfPhases();
        const int nw = wellsMultiSegment().size();
        for (int w = 0; w < nw; ++w) {
            const WellControls* wc = wellsMultiSegment()[w]->wellControls();
            // The current control in the well state overrides
            // the current control set in the Wells struct, which
            // is instead treated as a default.
            int current = xw.currentControls()[w];
            // Loop over all controls except the current one, and also
            // skip any RESERVOIR_RATE controls, since we cannot
            // handle those.
            const int nwc = well_controls_get_num(wc);
            int ctrl_index = 0;
            for (; ctrl_index < nwc; ++ctrl_index) {
                if (ctrl_index == current) {
                    // This is the currently used control, so it is
                    // used as an equation. So this is not used as an
                    // inequality constraint, and therefore skipped.
                    continue;
                }
                if (detail::constraintBroken(
                        xw.bhp(), xw.thp(), xw.wellRates(),
                        w, np, wellsMultiSegment()[w]->wellType(), wc, ctrl_index)) {
                    // ctrl_index will be the index of the broken constraint after the loop.
                    break;
                }
            }

            if (ctrl_index != nwc) {
                // Constraint number ctrl_index was broken, switch to it.
                if (terminal_output_)
                {
                    std::cout << "Switching control mode for well " << wellsMultiSegment()[w]->name()
                              << " from " << modestring[well_controls_iget_type(wc, current)]
                              << " to " << modestring[well_controls_iget_type(wc, ctrl_index)] << std::endl;
                }
                xw.currentControls()[w] = ctrl_index;
                current = xw.currentControls()[w];
            }

            // Get gravity for THP hydrostatic corrrection
            // const double gravity = detail::getGravity(geo_.gravity(), UgGridHelpers::dimensions(grid_));

            // Updating well state and primary variables.
            // Target values are used as initial conditions for BHP, THP, and SURFACE_RATE
            const double target = well_controls_iget_target(wc, current);
            const double* distr = well_controls_iget_distr(wc, current);
            switch (well_controls_iget_type(wc, current)) {
            case BHP:
                xw.bhp()[w] = target;
                xw.segPress()[xw.topSegmentLoc()[w]] = target;
                break;

            case THP: {
                OPM_THROW(std::runtime_error, "THP control is not implemented for multi-sgement wells yet!!");
            }

            case RESERVOIR_RATE:
                // No direct change to any observable quantity at
                // surface condition.  In this case, use existing
                // flow rates as initial conditions as reservoir
                // rate acts only in aggregate.
                break;

            case SURFACE_RATE:
                for (int phase = 0; phase < np; ++phase) {
                    if (distr[phase] > 0.0) {
                        xw.wellRates()[np * w + phase] = target * distr[phase];
                        // TODO: consider changing all (not just top) segment rates
                        // to make them consistent, it could possibly improve convergence.
                        xw.segPhaseRates()[np * xw.topSegmentLoc()[w] + phase] = target * distr[phase];
                    }
                }
                break;
            }

        }
    }





    template <class Grid>
    bool BlackoilMultiSegmentModel<Grid>::solveWellEq(const std::vector<ADB>& mob_perfcells,
                                                      const std::vector<ADB>& b_perfcells,
                                                      SolutionState& state,
                                                      WellState& well_state)
    {
        const bool converged = Base::solveWellEq(mob_perfcells, b_perfcells, state, well_state);

        if (converged) {
            // We must now update the state.segp and state.segqs members,
            // that the base version does not know about.
            const int np = numPhases();
            const int nseg_total =well_state.numSegments();
            {
                // We will set the segp primary variable to the new ones,
                // but we do not change the derivatives here.
                ADB::V new_segp = Eigen::Map<ADB::V>(well_state.segPress().data(), nseg_total);
                // Avoiding the copy below would require a value setter method
                // in AutoDiffBlock.
                std::vector<ADB::M> old_segp_derivs = state.segp.derivative();
                state.segp = ADB::function(std::move(new_segp), std::move(old_segp_derivs));
            }
            {
                // Need to reshuffle well rates, from phase running fastest
                // to wells running fastest.
                // The transpose() below switches the ordering.
                const DataBlock segrates = Eigen::Map<const DataBlock>(well_state.segPhaseRates().data(), nseg_total, np).transpose();
                ADB::V new_segqs = Eigen::Map<const V>(segrates.data(), nseg_total * np);
                std::vector<ADB::M> old_segqs_derivs = state.segqs.derivative();
                state.segqs = ADB::function(std::move(new_segqs), std::move(old_segqs_derivs));
            }

            // This is also called by the base version, but since we have updated
            // state.segp we must call it again.
            asImpl().computeWellConnectionPressures(state, well_state);
        }

        return converged;
    }





    template <class Grid>
    void BlackoilMultiSegmentModel<Grid>::addWellControlEq(const SolutionState& state,
                                                           const WellState& xw,
                                                           const V& aliveWells)
    {
        // the name of the function is a a little misleading.
        // Basically it is the function for the pressure equation.
        // And also, it work as the control equation when it is the segment
        if( wellsMultiSegment().empty() ) return;

        const int np = numPhases();
        const int nw = wellsMultiSegment().size();
        const int nseg_total = xw.numSegments();

        ADB aqua   = ADB::constant(ADB::V::Zero(nseg_total));
        ADB liquid = ADB::constant(ADB::V::Zero(nseg_total));
        ADB vapour = ADB::constant(ADB::V::Zero(nseg_total));

        if (active_[Water]) {
            aqua += subset(state.segqs, Span(nseg_total, 1, BlackoilPhases::Aqua * nseg_total));
        }
        if (active_[Oil]) {
            liquid += subset(state.segqs, Span(nseg_total, 1, BlackoilPhases::Liquid * nseg_total));
        }
        if (active_[Gas]) {
            vapour += subset(state.segqs, Span(nseg_total, 1, BlackoilPhases::Vapour * nseg_total));
        }

        // THP control is not implemented for the moment.

        // Hydrostatic correction variables
        ADB::V rho_v = ADB::V::Zero(nw);
        ADB::V vfp_ref_depth_v = ADB::V::Zero(nw);

        // Target vars
        ADB::V bhp_targets  = ADB::V::Zero(nw);
        ADB::V rate_targets = ADB::V::Zero(nw);
        Eigen::SparseMatrix<double>  rate_distr(nw, np*nw);

        // Selection variables
        // well selectors
        std::vector<int> bhp_well_elems;
        std::vector<int> rate_well_elems;
        // segment selectors
        std::vector<int> bhp_top_elems;
        std::vector<int> rate_top_elems;
        std::vector<int> rate_top_phase_elems;
        std::vector<int> others_elems;

        //Run through all wells to calculate BHP/RATE targets
        //and gather info about current control
        int start_segment = 0;
        for (int w = 0; w < nw; ++w) {
            const struct WellControls* wc = wellsMultiSegment()[w]->wellControls();

            // The current control in the well state overrides
            // the current control set in the Wells struct, which
            // is instead treated as a default.
            const int current = xw.currentControls()[w];

            const int nseg = wellsMultiSegment()[w]->numberOfSegments();

            switch (well_controls_iget_type(wc, current)) {
            case BHP:
            {
                bhp_well_elems.push_back(w);
                bhp_top_elems.push_back(start_segment);
                bhp_targets(w)  = well_controls_iget_target(wc, current);
                rate_targets(w) = -1e100;
                for (int p = 0; p < np; ++p) {
                    rate_top_phase_elems.push_back(np * start_segment + p);
                }
            }
            break;

            case THP:
            {
                OPM_THROW(std::runtime_error, "THP control is not implemented for multi-sgement wells yet!!");
            }
            break;

            case RESERVOIR_RATE: // Intentional fall-through
            case SURFACE_RATE:
            {
                rate_well_elems.push_back(w);
                rate_top_elems.push_back(start_segment);
                for (int p = 0; p < np; ++p) {
                    rate_top_phase_elems.push_back(np * start_segment + p);
                }
                // RESERVOIR and SURFACE rates look the same, from a
                // high-level point of view, in the system of
                // simultaneous linear equations.

                const double* const distr =
                    well_controls_iget_distr(wc, current);

                for (int p = 0; p < np; ++p) {
                    rate_distr.insert(w, p*nw + w) = distr[p];
                }

                bhp_targets(w)  = -1.0e100;
                rate_targets(w) = well_controls_iget_target(wc, current);
            }
            break;

            }

            for (int i = 1; i < nseg; ++i) {
                others_elems.push_back(i + start_segment);
            }
            start_segment += nseg;
        }

        // for each segment: 1, if the segment is the top segment, then control equation
        //                   2, if the segment is not the top segment, then the pressure equation
        const ADB bhp_residual = subset(state.segp, bhp_top_elems) - subset(bhp_targets, bhp_well_elems);
        const ADB rate_residual = subset(rate_distr * subset(state.segqs, rate_top_phase_elems) - rate_targets, rate_well_elems);

        ADB others_residual = ADB::constant(V::Zero(nseg_total));

        // Special handling for when we are called from solveWellEq().
        // TODO: restructure to eliminate need for special treatmemt.
        ADB wspd = (state.segp.numBlocks() == 2)
            ? detail::onlyWellDerivs(well_segment_pressures_delta_)
            : well_segment_pressures_delta_;

        others_residual = wops_ms_.eliminate_topseg * (state.segp - wops_ms_.s2s_outlet * state.segp + wspd);

        //       all the control equations
        // TODO: can be optimized better
        ADB well_eq_topsegment = subset(superset(bhp_residual, bhp_top_elems, nseg_total) +
                                        superset(rate_residual, rate_top_elems, nseg_total), xw.topSegmentLoc());

        // For wells that are dead (not flowing), and therefore not communicating
        // with the reservoir, we set the equation to be equal to the well's total
        // flow. This will be a solution only if the target rate is also zero.
        Eigen::SparseMatrix<double> rate_summer(nw, np*nw);
        for (int w = 0; w < nw; ++w) {
            for (int phase = 0; phase < np; ++phase) {
                rate_summer.insert(w, phase*nw + w) = 1.0;
            }
        }
        Selector<double> alive_selector(aliveWells, Selector<double>::NotEqualZero);
        // TODO: Here only handles the wells, or the top segments
        // should we also handle some non-alive non-top segments?
        // should we introduce the cocept of non-alive segments?
        // At the moment, we only handle the control equations
        well_eq_topsegment = alive_selector.select(well_eq_topsegment, rate_summer * subset(state.segqs, rate_top_phase_elems));

        /* residual_.well_eq = superset(bhp_residual, bhp_top_elems, nseg_total) +
                            superset(rate_residual, rate_top_elems, nseg_total) +
                            superset(others_residual, others_elems, nseg_total); */
        residual_.well_eq = superset(well_eq_topsegment, xw.topSegmentLoc(), nseg_total) +
                            others_residual;

    }





    template <class Grid>
    void
    BlackoilMultiSegmentModel<Grid>::updateWellState(const V& dwells,
                                                     WellState& well_state)
    {

        if (!wellsMultiSegment().empty())
        {
            const int np = numPhases();
            const int nw = wellsMultiSegment().size();
            const int nseg_total = well_state.numSegments();

            // Extract parts of dwells corresponding to each part.
            int varstart = 0;
            const V dsegqs = subset(dwells, Span(np * nseg_total, 1, varstart));
            varstart += dsegqs.size();
            const V dsegp = subset(dwells, Span(nseg_total, 1, varstart));
            varstart += dsegp.size();
            assert(varstart == dwells.size());
            const double dpmaxrel = dpMaxRel();


            // segment phase rates update
            // in dwells, the phase rates are ordered by phase.
            // while in WellStateMultiSegment, the phase rates are ordered by segments
            const DataBlock wsr = Eigen::Map<const DataBlock>(dsegqs.data(), np, nseg_total).transpose();
            const V dwsr = Eigen::Map<const V>(wsr.data(), nseg_total * np);
            const V wsr_old = Eigen::Map<const V>(&well_state.segPhaseRates()[0], nseg_total * np);
            const V sr = wsr_old - dwsr;
            std::copy(&sr[0], &sr[0] + sr.size(), well_state.segPhaseRates().begin());


            // segment pressure updates
            const V segp_old = Eigen::Map<const V>(&well_state.segPress()[0], nseg_total, 1);
            // TODO: applying the pressure change limiter to all the segments, not sure if it is the correct thing to do
            const V dsegp_limited = sign(dsegp) * dsegp.abs().min(segp_old.abs() * dpmaxrel);
            const V segp = segp_old - dsegp_limited;
            std::copy(&segp[0], &segp[0] + segp.size(), well_state.segPress().begin());

            // update the well rates and bhps, which are not anymore primary vabriables.
            // they are updated directly from the updated segment phase rates and segment pressures.

            // Bhp update.
            V bhp = V::Zero(nw);
            V wr = V::Zero(nw * np);
            // it is better to use subset

            int start_segment = 0;
            for (int w = 0; w < nw; ++w) {
                bhp[w] = well_state.segPress()[start_segment];
                // insert can be faster
                for (int p = 0; p < np; ++p) {
                    wr[p + np * w] = well_state.segPhaseRates()[p + np * start_segment];
                }

                const int nseg = wellsMultiSegment()[w]->numberOfSegments();
                start_segment += nseg;
            }

            assert(start_segment == nseg_total);
            std::copy(&bhp[0], &bhp[0] + bhp.size(), well_state.bhp().begin());
            std::copy(&wr[0], &wr[0] + wr.size(), well_state.wellRates().begin());

            // TODO: handling the THP control related.
        }
    }




    template <class Grid>
    void
    BlackoilMultiSegmentModel<Grid>::computeSegmentFluidProperties(const SolutionState& state)
    {
        const int nw = wellsMultiSegment().size();
        const int nseg_total = state.segp.size();
        const int np = numPhases();

        // although we will calculate segment density for non-segmented wells at the same time,
        // while under most of the cases, they will not be used,
        // since for most of the cases, the density calculation for non-segment wells are
        // set to be 'SEG' way, which is not a option for multi-segment wells.
        // When the density calcuation for non-segmented wells are set to 'AVG', then
        // the density calculation of the mixtures can be the same, while it remains to be verified.

        // The grid cells associated with segments.
        // TODO: shoud be computed once and stored in WellState or global Wells structure or class.
        std::vector<int> segment_cells;
        segment_cells.reserve(nseg_total);
        for (int w = 0; w < nw; ++w) {
            const std::vector<int>& segment_cells_well = wellsMultiSegment()[w]->segmentCells();
            segment_cells.insert(segment_cells.end(), segment_cells_well.begin(), segment_cells_well.end());
        }
        assert(int(segment_cells.size()) == nseg_total);

        const ADB segment_temp = subset(state.temperature, segment_cells);
        // using the segment pressure or the average pressure
        // using the segment pressure first
        const ADB& segment_press = state.segp;

        // Compute PVT properties for segments.
        std::vector<PhasePresence> segment_cond(nseg_total);
        const std::vector<PhasePresence>& pc = phaseCondition();
        for (int s = 0; s < nseg_total; ++s) {
            segment_cond[s] = pc[segment_cells[s]];
        }
        std::vector<ADB> b_seg(np, ADB::null());
        // Viscosities for different phases
        std::vector<ADB> mu_seg(np, ADB::null());
        ADB rsmax_seg = ADB::null();
        ADB rvmax_seg = ADB::null();
        const PhaseUsage& pu = fluid_.phaseUsage();
        if (pu.phase_used[Water]) {
            b_seg[pu.phase_pos[Water]] = fluid_.bWat(segment_press, segment_temp, segment_cells);
            mu_seg[pu.phase_pos[Water]] = fluid_.muWat(segment_press, segment_temp, segment_cells);
        }
        assert(active_[Oil]);
        const ADB segment_so = subset(state.saturation[pu.phase_pos[Oil]], segment_cells);
        if (pu.phase_used[Oil]) {
            const ADB segment_rs = subset(state.rs, segment_cells);
            b_seg[pu.phase_pos[Oil]] = fluid_.bOil(segment_press, segment_temp, segment_rs,
                                                   segment_cond, segment_cells);
            rsmax_seg = fluidRsSat(segment_press, segment_so, segment_cells);
            mu_seg[pu.phase_pos[Oil]] = fluid_.muOil(segment_press, segment_temp, segment_rs,
                                                     segment_cond, segment_cells);
        }
        assert(active_[Gas]);
        if (pu.phase_used[Gas]) {
            const ADB segment_rv = subset(state.rv, segment_cells);
            b_seg[pu.phase_pos[Gas]] = fluid_.bGas(segment_press, segment_temp, segment_rv,
                                                   segment_cond, segment_cells);
            rvmax_seg = fluidRvSat(segment_press, segment_so, segment_cells);
            mu_seg[pu.phase_pos[Gas]] = fluid_.muGas(segment_press, segment_temp, segment_rv,
                                                   segment_cond, segment_cells);
        }

        // Extract segment flow by phase (segqs) and compute total surface rate.
        ADB tot_surface_rate = ADB::constant(V::Zero(nseg_total));
        std::vector<ADB> segqs(np, ADB::null());
        for (int phase = 0; phase < np; ++phase) {
            segqs[phase] = subset(state.segqs, Span(nseg_total, 1, phase * nseg_total));
            tot_surface_rate += segqs[phase];
        }

        // TODO: later this will be implmented as a global mapping
        std::vector<std::vector<double>> comp_frac(np, std::vector<double>(nseg_total, 0.0));
        int start_segment = 0;
        for (int w = 0; w < nw; ++w) {
            WellMultiSegmentConstPtr well = wellsMultiSegment()[w];
            const int nseg = well->numberOfSegments();
            const std::vector<double>& comp_frac_well = well->compFrac();
            for (int phase = 0; phase < np; ++phase) {
                for (int s = 0; s < nseg; ++s) {
                    comp_frac[phase][s + start_segment] = comp_frac_well[phase];
                }
            }
            start_segment += nseg;
        }
        assert(start_segment == nseg_total);

        // Compute mix.
        // 'mix' contains the component fractions under surface conditions.
        std::vector<ADB> mix(np, ADB::null());
        for (int phase = 0; phase < np; ++phase) {
            // initialize to be the compFrac for each well,
            // then update only the one with non-zero total volume rate
            mix[phase] = ADB::constant(Eigen::Map<V>(comp_frac[phase].data(), nseg_total));
        }
        // There should be a better way to do this.
        Selector<double> non_zero_tot_rate(tot_surface_rate.value(), Selector<double>::NotEqualZero);
        for (int phase = 0; phase < np; ++phase) {
            mix[phase] = non_zero_tot_rate.select(segqs[phase] / tot_surface_rate, mix[phase]);
        }

        // Calculate rs and rv.
        ADB rs = ADB::constant(V::Zero(nseg_total));
        ADB rv = rs;
        const int gaspos = pu.phase_pos[Gas];
        const int oilpos = pu.phase_pos[Oil];
        Selector<double> non_zero_mix_oilpos(mix[oilpos].value(), Selector<double>::GreaterZero);
        Selector<double> non_zero_mix_gaspos(mix[gaspos].value(), Selector<double>::GreaterZero);
        // What is the better way to do this?
        // big values should not be necessary
        ADB big_values = ADB::constant(V::Constant(nseg_total, 1.e100));
        ADB mix_gas_oil = non_zero_mix_oilpos.select(mix[gaspos] / mix[oilpos], big_values);
        ADB mix_oil_gas = non_zero_mix_gaspos.select(mix[oilpos] / mix[gaspos], big_values);
        if (active_[Oil]) {
            V selectorUnderRsmax = V::Zero(nseg_total);
            V selectorAboveRsmax = V::Zero(nseg_total);
            for (int s = 0; s < nseg_total; ++s) {
                if (mix_gas_oil.value()[s] > rsmax_seg.value()[s]) {
                    selectorAboveRsmax[s] = 1.0;
                } else {
                    selectorUnderRsmax[s] = 1.0;
                }
            }
            rs = non_zero_mix_oilpos.select(selectorAboveRsmax * rsmax_seg + selectorUnderRsmax * mix_gas_oil, rs);
        }
        if (active_[Gas]) {
            V selectorUnderRvmax = V::Zero(nseg_total);
            V selectorAboveRvmax = V::Zero(nseg_total);
            for (int s = 0; s < nseg_total; ++s) {
                if (mix_oil_gas.value()[s] > rvmax_seg.value()[s]) {
                    selectorAboveRvmax[s] = 1.0;
                } else {
                    selectorUnderRvmax[s] = 1.0;
                }
            }
            rv = non_zero_mix_gaspos.select(selectorAboveRvmax * rvmax_seg + selectorUnderRvmax * mix_oil_gas, rv);
        }

        // Calculate the phase fraction under reservoir conditions.
        std::vector<ADB> x(np, ADB::null());
        for (int phase = 0; phase < np; ++phase) {
            x[phase] = mix[phase];
        }
        if (active_[Gas] && active_[Oil]) {
            x[gaspos] = (mix[gaspos] - mix[oilpos] * rs) / (V::Ones(nseg_total) - rs * rv);
            x[oilpos] = (mix[oilpos] - mix[gaspos] * rv) / (V::Ones(nseg_total) - rs * rv);
        }

        // Compute total reservoir volume to surface volume ratio.
        ADB volrat = ADB::constant(V::Zero(nseg_total));
        for (int phase = 0; phase < np; ++phase) {
            volrat += x[phase] / b_seg[phase];
        }

        // Compute segment densities.
        ADB dens = ADB::constant(V::Zero(nseg_total));
        for (int phase = 0; phase < np; ++phase) {
            const V surface_density = fluid_.surfaceDensity(phase, segment_cells);
            dens += surface_density * mix[phase];
        }
        well_segment_densities_ = dens / volrat;

        // Calculating the surface volume of each component in the segment
        assert(np == int(segment_comp_surf_volume_current_.size()));
        const ADB segment_surface_volume = segvdt_ / volrat;
        for (int phase = 0; phase < np; ++phase) {
            segment_comp_surf_volume_current_[phase] = segment_surface_volume * mix[phase];
        }

        // Mass flow rate of the segments
        segment_mass_flow_rates_ = ADB::constant(V::Zero(nseg_total));
        for (int phase = 0; phase < np; ++phase) {
            // TODO: how to remove one repeated surfaceDensity()
            const V surface_density = fluid_.surfaceDensity(phase, segment_cells);
            segment_mass_flow_rates_ += surface_density * segqs[phase];
        }

        // Viscosity of the fluid mixture in the segments
        segment_viscosities_ = ADB::constant(V::Zero(nseg_total));
        for (int phase = 0; phase < np; ++phase) {
            segment_viscosities_ += x[phase] * mu_seg[phase];
        }
    }


    template <class Grid>
    void
    BlackoilMultiSegmentModel<Grid>::computeSegmentPressuresDelta(const SolutionState& state)
    {
        const int nw = wellsMultiSegment().size();
        const int nseg_total = state.segp.size();

        // calculate the depth difference of the segments
        // TODO: we need to store the following values somewhere to avoid recomputation.
        V segment_depth_delta = V::Zero(nseg_total);
        int start_segment = 0;
        for (int w = 0; w < nw; ++w) {
            WellMultiSegmentConstPtr well = wellsMultiSegment()[w];
            const int nseg = well->numberOfSegments();
            for (int s = 1; s < nseg; ++s) {
                const int s_outlet = well->outletSegment()[s];
                assert(s_outlet >= 0 && s_outlet < nseg);
                segment_depth_delta[s + start_segment] = well->segmentDepth()[s_outlet] - well->segmentDepth()[s];
            }
            start_segment += nseg;
        }
        assert(start_segment == nseg_total);

        const double grav = detail::getGravity(geo_.gravity(), UgGridHelpers::dimensions(grid_));
        const ADB grav_adb = ADB::constant(V::Constant(nseg_total, grav));
        well_segment_pressures_delta_ = segment_depth_delta * grav_adb * well_segment_densities_;

        ADB well_segment_perforation_densities = wops_ms_.s2p * well_segment_densities_;
        well_segment_perforation_pressure_diffs_ = grav * well_segment_perforation_depth_diffs_ * well_segment_perforation_densities;

    }

} // namespace Opm

#endif // OPM_BLACKOILMODELBASE_IMPL_HEADER_INCLUDED
