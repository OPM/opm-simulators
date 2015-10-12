/*
  Copyright 2013, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014, 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2014, 2015 Statoil ASA.
  Copyright 2015 NTNU
  Copyright 2015 IRIS AS

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
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/Exceptions.hpp>
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
          // not all will be necessary eventually
        , well_perforation_densities_adb_(ADB::null())
        , well_perforation_pressure_diffs_adb_(ADB::null())
        , well_perforation_pressure_cell_diffs_adb_(ADB::null())
        , well_perforations_segment_pressure_diffs_(ADB::null())
        , well_perforation_cell_densities_adb_(ADB::null())
        , well_segment_densities_(ADB::null())
        , wells_multisegment_(wells_multisegment)
        {

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

        //TODO: handle the volume related.
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
            // push null sates for qs and pseg
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
        // TODO: using the original Qs for the segment rates for now.
        // TODO: using the original Bhp for the segment pressures for now.

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




    // TODO: This is just a WIP version
    // TODO: To be fxied to handle the situation with both usual wells and mutli-segment wells.
    template <class Grid>
    void BlackoilMultiSegmentModel<Grid>::computeWellConnectionPressures(const SolutionState& state,
                                                                         const WellState& xw)
    {
        if( ! wellsActive() ) return ;

        using namespace Opm::AutoDiffGrid;
        // 1. Compute properties required by computeConnectionPressureDelta().
        //    Note that some of the complexity of this part is due to the function
        //    taking std::vector<double> arguments, and not Eigen objects.
        const int nperf = xw.numPerforations();
        const int nw = xw.numWells();

        // the well cells for multisegment wells and non-segmented wells should be counted seperatedly
        // indexing should be put in WellState
        std::vector<int> well_cells;
        std::vector<int> well_cells_segmented_idx;
        std::vector<int> well_cells_non_segmented_idx;
        std::vector<int> well_cells_segmented;
        std::vector<int> well_cells_non_segmented;

        well_cells_segmented_idx.reserve(nperf);
        well_cells_non_segmented_idx.reserve(nperf);
        well_cells_segmented.reserve(nperf);
        well_cells_non_segmented.reserve(nperf);
        well_cells.reserve(nperf);
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

        // compute the average of the fluid densites in the well blocks.
        // the average is weighted according to the fluid relative permeabilities.
        const std::vector<ADB> kr_adb = Base::computeRelPerm(state);
        size_t temp_size = kr_adb.size();
        std::vector<V> perf_kr;
        for(size_t i = 0; i < temp_size; ++i) {
            const ADB kr_phase_adb = subset(kr_adb[i], well_cells);
            const V kr_phase = kr_phase_adb.value();
            perf_kr.push_back(kr_phase);
        }

        // compute the averaged density for the well block
        // TODO: for the non-segmented wells, they should be set to zero
        // TODO: for the moment, they are still calculated, while not used later.
        for (int i = 0; i < nperf; ++i) {
            double sum_kr = 0.;
            int np = perf_kr.size(); // make sure it is 3
            for (int p = 0;  p < np; ++p) {
                sum_kr += perf_kr[p][i];
            }

            for (int p = 0; p < np; ++p) {
                perf_kr[p][i] /= sum_kr;
            }
        }

        // ADB rho_avg_perf = ADB::const(V::Zero(nperf);
        V rho_avg_perf(nperf);
        // get the fluid densities to do the average?
        // TODO: make sure the order of the density and the order of the kr are the same.
        for (int phaseIdx = 0; phaseIdx < fluid_.numPhases(); ++phaseIdx) {
            const int canonicalPhaseIdx = canph_[phaseIdx];
            const ADB rho = fluidDensity(canonicalPhaseIdx, rq_[phaseIdx].b, state.rs, state.rv);
            const V rho_perf = subset(rho, well_cells).value();
            // TODO: phaseIdx or canonicalPhaseIdx ?
            rho_avg_perf = rho_perf * perf_kr[phaseIdx];
        }

        // TODO: rho_avg_perf is actually the well_perforation_cell_densities_;
        well_perforation_cell_densities_ = Eigen::Map<const V>(rho_avg_perf.data(), nperf);
        // TODO: just for the current cases
        well_perforation_densities_ = V::Zero(nperf);

        // For the non-segmented wells

        // const std::vector<ADB> perf_kr_adb = subset(kr_adb, well_cells);

        // const V perf_kr = perf_kr_adb.value();

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

        const V perf_press = Eigen::Map<const V>(xw.perfPress().data(), nperf);
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

            // for (int perf = wells().well_connpos[w]; perf < wells().well_connpos[w+1]; ++perf) {
            const int start_perforation = (*it_well).second.start_perforation;
            const int end_perforation = start_perforation + (*it_well).second.number_of_perforations;
            for (int perf = start_perforation; perf < end_perforation; ++perf) {
                const double p_above = perf == start_perforation ? state.segp.value()[start_segment] : perf_press[perf - 1];
                const double p_avg = (perf_press[perf] + p_above)/2;
                avg_press[perf] = p_avg;
            }
            // std::cout << " the bhp value is " << state.segp.value()[start_segment];
            start_segment += nseg;
        }

#if 0
        std::cout << " state.segp " << std::endl;
        std::cout << state.segp.value() << std::endl;
        std::cout << " perf_press " << std::endl;
        std::cout << perf_press << std::endl;
        std::cout << " avg_press " << std::endl;
        std::cout << avg_press << std::endl;
#endif


        // Use cell values for the temperature as the wells don't knows its temperature yet.
        const ADB perf_temp = subset(state.temperature, well_cells);

        // Compute b, rsmax, rvmax values for perforations.
        // Evaluate the properties using average well block pressures
        // and cell values for rs, rv, phase condition and temperature.
        const ADB avg_press_ad = ADB::constant(avg_press);
        std::vector<PhasePresence> perf_cond(nperf);
        const std::vector<PhasePresence>& pc = phaseCondition();
        for (int perf = 0; perf < nperf; ++perf) {
            perf_cond[perf] = pc[well_cells[perf]];
        }
        const PhaseUsage& pu = fluid_.phaseUsage();
        DataBlock b(nperf, pu.num_phases);
        std::vector<double> rsmax_perf(nperf, 0.0);
        std::vector<double> rvmax_perf(nperf, 0.0);
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
            rsmax_perf.assign(rssat.data(), rssat.data() + nperf);
        }
        if (pu.phase_used[BlackoilPhases::Vapour]) {
            const ADB perf_rv = subset(state.rv, well_cells);
            const V bg = fluid_.bGas(avg_press_ad, perf_temp, perf_rv, perf_cond, well_cells).value();
            b.col(pu.phase_pos[BlackoilPhases::Vapour]) = bg;
            const V rvsat = fluidRvSat(avg_press, perf_so, well_cells);
            rvmax_perf.assign(rvsat.data(), rvsat.data() + nperf);
        }
        // b is row major, so can just copy data.
        std::vector<double> b_perf(b.data(), b.data() + nperf * pu.num_phases);
        // Extract well connection depths.
        const V depth = cellCentroidsZToEigen(grid_);
        const V pdepth = subset(depth, well_cells);
        std::vector<double> perf_depth(pdepth.data(), pdepth.data() + nperf);

        // Surface density.
        DataBlock surf_dens(nperf, pu.num_phases);
        for (int phase = 0; phase < pu.num_phases; ++ phase) {
            surf_dens.col(phase) = V::Constant(nperf, fluid_.surfaceDensity()[pu.phase_pos[phase]]);
        }

        std::vector<double> surf_dens_perf(surf_dens.data(), surf_dens.data() + nperf * pu.num_phases);

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
                        wells(), perf_depth, cd, grav);

        // 4. Store the results
        well_perforation_densities_ = Eigen::Map<const V>(cd.data(), nperf); // This one is not useful for segmented wells at all
        well_perforation_pressure_diffs_ = Eigen::Map<const V>(cdp.data(), nperf);

        // TODO: A temporary approach. We calculate all the densities and pressure difference for all the perforations.
        // For the segmented wells, the h_nc;
        // Firstly, we need to compute the segmented densities first.
        // It must be implicit. So it must be ADB variable.
        // If we need to consider the rs and rv for all the segments,  the process will be similar with the above.
        // Are they actually zero for the current cases?
        // TODO
        well_perforations_segment_pressure_diffs_ = ADB::constant(V::Zero(xw.numPerforations()));
        well_perforation_pressure_cell_diffs_ = V::Zero(xw.numPerforations());
        well_perforatoin_cell_pressure_diffs_ = V::Zero(xw.numPerforations());
#if 0
        std::cout << "well_perforation_densities_ " << std::endl;
        std::cout << well_perforation_densities_ << std::endl;

        std::cout << "well_perforation_pressure_diffs_ " << std::endl;
        std::cout << well_perforation_pressure_diffs_ << std::endl;
#endif
    }







    template <class Grid>
    void
    BlackoilMultiSegmentModel<Grid>::
    assemble(const ReservoirState& reservoir_state,
             WellState& well_state,
             const bool initial_assembly)
    {
        using namespace Opm::AutoDiffGrid;

        // TODO include VFP effect.
        // If we have VFP tables, we need the well connection
        // pressures for the "simple" hydrostatic correction
        // between well depth and vfp table depth.
        //  if (isVFPActive()) {
        //     SolutionState state = asImpl().variableState(reservoir_state, well_state);
        //     SolutionState state0 = state;
        //     asImpl().makeConstantState(state0);
        //     computeWellConnectionPressures(state0, well_state);
        // }

        // Possibly switch well controls and updating well state to
        // get reasonable initial conditions for the wells
        updateWellControls(well_state);

        // Create the primary variables.
        SolutionState state = variableState(reservoir_state, well_state);

        if (initial_assembly) {
            // Create the (constant, derivativeless) initial state.
            SolutionState state0 = state;
            makeConstantState(state0);
            // Compute initial accumulation contributions
            // and well connection pressures.
            Base::computeAccum(state0, 0);
            computeWellConnectionPressures(state0, well_state);
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
        Base::assembleMassBalanceEq(state);

        // -------- Well equations ----------

        if ( ! wellsActive() ) {
            return;
        }

        V aliveWells;
        // const int np = wells().number_of_phases;
        const int np = well_state.numPhases();
        std::vector<ADB> cq_s(np, ADB::null());

        // const int nw = wellsMultiSegment().size();
        const int nw = well_state.numWells();
        const int nperf = well_state.numPerforations();
        std::vector<int> well_cells;
        well_cells.reserve(nperf);
        for (int i = 0; i < nw; ++i) {
            const std::vector<int>& temp_well_cells = wellsMultiSegment()[i]->wellCells();
            well_cells.insert(well_cells.end(), temp_well_cells.begin(), temp_well_cells.end());
        }

        assert(nperf == int(well_cells.size()));

        std::vector<ADB> mob_perfcells(np, ADB::null());
        std::vector<ADB> b_perfcells(np, ADB::null());
        for (int phase = 0; phase < np; ++phase) {
            mob_perfcells[phase] = subset(rq_[phase].mob, well_cells);
            b_perfcells[phase] = subset(rq_[phase].b, well_cells);
        }

        // TODO: it will be a good thing to try to solve welleq seperately every time
        // if (param_.solve_welleq_initially_ && initial_assembly) {
            // solve the well equations as a pre-processing step
        //     solveWellEq(mob_perfcells, b_perfcells, state, well_state);
        // }

        // the perforation flux here are different
        // it is related to the segment location
        computeWellFlux(state, mob_perfcells, b_perfcells, aliveWells, cq_s);
        computeSegmentDensities(state, cq_s, b_perfcells, well_state);
        updatePerfPhaseRatesAndPressures(cq_s, state, well_state);
        addWellFluxEq(cq_s, state);
        addWellContributionToMassBalanceEq(cq_s, state, well_state);
        addWellControlEq(state, well_state, aliveWells);
    }





    template <class Grid>
    void
    BlackoilMultiSegmentModel<Grid>::addWellContributionToMassBalanceEq(const std::vector<ADB>& cq_s,
                                                                const SolutionState&,
                                                                const WellState& xw)
    {
        // For the version at the moment, it has to be done one well by one well
        // later, we may need to develop a new wells class for optimization.
        const int nw = wellsMultiSegment().size();
        const int nc = Opm::AutoDiffGrid::numCells(grid_);
        const int np = numPhases();
        const int nperf = xw.numPerforations();

        std::vector<int> well_cells;

        for (int w = 0; w < nw; ++w) {
            WellMultiSegmentConstPtr well = wellsMultiSegment()[w];
            well_cells.insert(well_cells.end(), well->wellCells().begin(), well->wellCells().end());
        }

        assert(int(well_cells.size()) == nperf);

        for (int phase = 0; phase < np; ++phase) {
            residual_.material_balance_eq[phase] -= superset(cq_s[phase], well_cells, nc);
        }


        // TODO: it must be done one by one?
        // or we develop a new Wells class?
        static_cast<void>(nperf);
        // Add well contributions to mass balance equations
        // const int nc = Opm::AutoDiffGrid::numCells(grid_);
        // const int nw = wells().number_of_wells;
        // const int nperf = wells().well_connpos[nw];
        // const int np = wells().number_of_phases;
        // const std::vector<int> well_cells(wells().well_cells, wells().well_cells + nperf);
        // for (int phase = 0; phase < np; ++phase) {
        //     residual_.material_balance_eq[phase] -= superset(cq_s[phase], well_cells, nc);
        // }
    }





    template <class Grid>
    void
    BlackoilMultiSegmentModel<Grid>::computeWellFlux(const SolutionState& state,
                                                             const std::vector<ADB>& ,
                                                             const std::vector<ADB>& ,
                                                             V& aliveWells,
                                                             std::vector<ADB>& cq_s)
    {
        // if( ! wellsActive() ) return ;
        if (wellsMultiSegment().size() == 0) return;


        // TODO: AS the first version, we handle well one by one for two reasons
        // 1: trying to handle the segmented wells and non-segmented wells seperately,
        //    before we can handle the non-segmented wells in a segmented way
        // 2: currently, wells are stored in a vector.
        //    After we confirm that we can handle non-segmented wells in a segmented way,
        //    maybe we will have a wellsMultisegment class similar to the previous Wells structure,
        //    so that we can handle all the wells togeter.
        // At the moment, let us hanldle wells one by one.
        // For the moment, the major purpose of this function is to calculate all the perforation phase rates.

        const int nw = wellsMultiSegment().size();
        const Opm::PhaseUsage& pu = fluid_.phaseUsage();

        int start_perforation = 0;
        int start_segment = 0;

        aliveWells = V::Constant(nw, 1.0);


        // temporary, no place to store the information about total perforations and segments
        int total_nperf = 0;
        for (int w = 0; w < nw; ++w) {
            total_nperf += wellsMultiSegment()[w]->numberOfPerforations();
        }
        const int np = numPhases();

        for (int p = 0; p < np; ++p) {
            cq_s[p] = ADB::constant(V::Zero(total_nperf));
        }

        for (int w = 0; w < nw; ++w) {
            WellMultiSegmentConstPtr well = wellsMultiSegment()[w];
            const int nseg = well->numberOfSegments();
            const int nperf = well->numberOfPerforations();

            V Tw = Eigen::Map<const V>(well->wellIndex().data(), nperf);
            const std::vector<int>& well_cells = well->wellCells();

            // extract mob_perfcells and b_perfcells.
            std::vector<ADB> mob_perfcells(np, ADB::null());
            std::vector<ADB> b_perfcells(np, ADB::null());
            for (int phase = 0; phase < np; ++phase) {
                mob_perfcells[phase] = subset(rq_[phase].mob, well_cells);
                b_perfcells[phase] = subset(rq_[phase].b, well_cells);
            }

            // determining in-flow (towards well-bore) or out-flow (towards reservoir)
            // for mutli-segmented wells and non-segmented wells, the calculation of the drawdown are different.
            const ADB& p_perfcells = subset(state.pressure, well_cells);
            const ADB& rs_perfcells = subset(state.rs, well_cells);
            const ADB& rv_perfcells = subset(state.rv, well_cells);

            const ADB& seg_pressures = subset(state.segp, Span(nseg, 1, start_segment));

            ADB drawdown = ADB::null();

            const ADB seg_pressures_perf = well->wellOps().s2p * seg_pressures;

            if (well->isMultiSegmented())
            {
                // get H_nc
                const ADB& h_nc = subset(well_perforations_segment_pressure_diffs_, Span(nperf, 1, start_perforation));
                const V& h_cj = subset(well_perforatoin_cell_pressure_diffs_, Span(nperf, 1, start_perforation));

                // V seg_pressures_perf = V::Zero(nperf);
                // for (int i = 0; i < nseg; ++i) {
                //     int temp_nperf = well->segmentPerforations()[i].size();
                //     assert(temp_nperf > 0);
                //     for (int j = 0; j < temp_nperf; ++j) {
                //         int index_perf = well->segmentPerforations()[i][j];
                //         assert(index_perf <= nperf);
                        // set the perforation pressure to be the segment pressure
                        // similiar to a scatter operation
                //         seg_pressures_perf[index_perf] = seg_pressures.value()[i];
                //     }
                // }

                drawdown = (p_perfcells + h_cj - seg_pressures_perf - h_nc);
            }
            else
            // usual wells
            // only one segment each
            {
                const V& cdp = subset(well_perforation_pressure_diffs_, Span(nperf, 1, start_perforation));

                const ADB perf_pressures = well->wellOps().s2p * seg_pressures + cdp;

                drawdown = p_perfcells - perf_pressures;
            }

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
            // TODO: if the wells are producing wells.
            // TODO: if the wells are injection wells.
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
            // The wellbore mixture depends on the inflow from the reservoar
            // and the well injection rates.
            // compute avg. and total wellbore phase volumetric rates at standard conds

           // TODO: should this based on the segments?
           // TODO: for the usual wells, the well rates are the sum of the perforations.
           // TODO: for multi-segmented wells, the segment rates are not the sum of the perforations.

           // TODO: two options here
           // TODO: 1. for each segment, only the inflow from the perforations related to this segment are considered.
           // TODO: 2. for each segment, the inflow from the perforrations related to this segment and also all the inflow
           // TODO: from the upstreaming sgments and their perforations need to be considered.
           // TODO: This way can be the more consistent way, while let is begin with the first option. The second option
           // TODO: involves one operations that are not valid now. (i.e. how to transverse from the leaves to the root,
           // TODO: although we can begin from the brutal force way)

            const std::vector<double>& compi = well->compFrac();
            std::vector<ADB> wbq(np, ADB::null());
            ADB wbqt = ADB::constant(V::Zero(nseg));

            for (int phase = 0; phase < np; ++phase) {

                // const ADB& q_ps = well->wellOps().p2s * cq_ps[phase];
                const ADB& q_ps = well->wellOps().p2s_gather * cq_ps[phase];
                const int n_total_segments = state.segp.size();
                const ADB& q_s = subset(state.segqs, Span(nseg, 1, phase * n_total_segments + start_segment));
                Selector<double> injectingPhase_selector(q_s.value(), Selector<double>::GreaterZero);

                const int pos = pu.phase_pos[phase];

                // this is per segment
                wbq[phase] = (compi[pos] * injectingPhase_selector.select(q_s, ADB::constant(V::Zero(nseg)))) - q_ps;

                // TODO: it should be a single value for this certain well.
                // TODO: it need to be changed later to handle things more consistently
                // or there should be an earsier way to decide if the well is dead.
                wbqt += wbq[phase];
            }

            // the first value of the wbqt is the one to decide if the well is dead
            // or there should be some dead segments?
            // maybe not.
            if (wbqt.value()[0] == 0.) {
                aliveWells[w] = 0.;
            }

            // compute wellbore mixture at standard conditons
            // before, the determination of alive wells is based on wells.
            // now, will there be any dead segment? I think no.
            std::vector<ADB> cmix_s(np, ADB::constant(V::Zero(nperf)));
            if (aliveWells[w] > 0.) {
                for (int phase = 0; phase < np; ++phase) {
                    // const int pos = pu.phase_pos[phase];
                    // TODO: make sure wheter phase or pos be required here.
                    cmix_s[phase] = well->wellOps().s2p * (wbq[phase] / wbqt);
                }
            } else {
                for (int phase = 0; phase < np; ++phase) {
                    const int pos = pu.phase_pos[phase];
                    cmix_s[phase] += ADB::constant(V::Zero(nperf) + compi[pos]);
                }
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
                cq_s[phase] += superset(cq_ps[phase] + cmix_s[phase]*cqt_is, Span(nperf, 1, start_perforation), total_nperf);
            }

            start_perforation += nperf;
            start_segment += nseg;
        }

    }





    template <class Grid>
    void BlackoilMultiSegmentModel<Grid>::updatePerfPhaseRatesAndPressures(const std::vector<ADB>& cq_s,
                                                                           const SolutionState& state,
                                                                           WellState& xw)
    {
        // Update the perforation phase rates (used to calculate the pressure drop in the wellbore).
        // TODO: now it is so necesary to have a gobal wellsMultiSegment class to store some global information.
        const int np = numPhases();
        const int nw = wellsMultiSegment().size();
        const int nperf_total = xw.perfPress().size();

        V cq = superset(cq_s[0].value(), Span(nperf_total, np, 0), nperf_total * np);
        for (int phase = 1; phase < np; ++phase) {
            cq += superset(cq_s[phase].value(), Span(nperf_total, np, phase), nperf_total * np);
        }
        xw.perfPhaseRates().assign(cq.data(), cq.data() + nperf_total * np);

        // Update the perforation pressures for usual wells first to recover the resutls
        // without mutlti segment wells.
        // For segment wells, it has not been decided if we need th concept
        // of preforation pressures
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
        // TODO: how to update segment pressures and segment rates?
        // or we do not need here?

        // TODO: update the perforation pressures.
        // it should be based on the segment pressures
        // it should be done at least for the non-segmented wells
        // Then it makes it necessary to update the segment pressures and phase rates.
    }



    template <class Grid>
    void BlackoilMultiSegmentModel<Grid>::addWellFluxEq(const std::vector<ADB>& cq_s,
                                                        const SolutionState& state)
    {
        // the equations is for each segment
        const int np = numPhases();
        const int nw = wellsMultiSegment().size();
        const int nseg_total = state.segp.size();

        ADB segqs = state.segqs;

        for (int phase = 0; phase < np; ++phase) {
            int start_segment = 0;
            int start_perforation = 0;
            for (int w = 0; w < nw; ++w) {
                WellMultiSegmentConstPtr well = wellsMultiSegment()[w];
                // the equation is
                // /deta m_p_n - /sigma Q_pi - /sigma q_pj + Q_pn = 0
                // Q_pn + /deta m_p_n - /sigma Q_pi - /sigma q_pj = 0
                // 1. for the first term, we need the information from the previous step.
                //    in term of stock-tank conditions.
                //    it should the total volume of the well bore * volume ratio for each phase.
                //    it can be the phase rate / total rates (in surface units)
                //    For the first one.
                //    it will be (V_1 * S_1 - V_0 * S_0) / delta_T
                //    so we need the information from the previous step and the time step.
                // 2. for the second term, it is the inlet segments, which are also unknowns
                //    TODO:we can have a mapping for the inlet segments
                // 3. for the third term, it is the inflow.
                // 4. for the last term, it is the outlet rates, which are also unknowns

                // For this version, we will ignore the wellbore volume effects
                // In the DATA file, the well bore volumes are set to be zero already.
                // The equation is Q_pn - /sigma Q_pi - /sigma q_pj = 0;
                const int nperf = well->numberOfPerforations();
                const int nseg = well->numberOfSegments();

                // perforate rates for this well
                const ADB& cq_s_perf = subset(cq_s[phase], Span(nperf, 1, start_perforation));

                // sum of the perforate rates to its related segment
                const ADB& cq_s_seg = well->wellOps().p2s_gather * cq_s_perf;

                const int start_position = start_segment + phase * nseg_total;
                // the segment rates of this well
                const ADB& segqs_well = subset(segqs, Span(nseg, 1, start_position));

                if (well->isMultiSegmented()) {
                    segqs -= superset(cq_s_seg + well->wellOps().s2s_inlets * segqs_well, Span(nseg, 1, start_position), np * nseg_total);
                }
                else
                {
                    segqs -= superset(cq_s_seg, Span(1, 1, start_position), np * nseg_total);
                }
                start_segment += nseg;
                start_perforation += nperf;
            }
            assert(start_segment == nseg_total);
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
        // const Opm::PhaseUsage& pu = fluid_.phaseUsage();
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

#if 0
            // Debug output
            std::cout << " xw.bhp() is " << xw.bhp()[w] << std::endl;
            std::cout << " ctrl_index " << ctrl_index << " nwc " << nwc << std::endl;
#endif

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
                /* double aqua = 0.0;
                double liquid = 0.0;
                double vapour = 0.0;

                if (active_[ Water ]) {
                    aqua = xw.wellRates()[w*np + pu.phase_pos[ Water ] ];
                }
                if (active_[ Oil ]) {
                    liquid = xw.wellRates()[w*np + pu.phase_pos[ Oil ] ];
                }
                if (active_[ Gas ]) {
                    vapour = xw.wellRates()[w*np + pu.phase_pos[ Gas ] ];
                }

                const int vfp        = well_controls_iget_vfp(wc, current);
                const double& thp    = well_controls_iget_target(wc, current);
                const double& alq    = well_controls_iget_alq(wc, current);

                //Set *BHP* target by calculating bhp from THP
                const WellType& well_type = wellsMultiSegment()[w]->wellType;

                if (well_type == INJECTOR) {
                    // TODO: this needs to be updated
                    double dp = detail::computeHydrostaticCorrection(
                            wells(), w, vfp_properties_.getInj()->getTable(vfp)->getDatumDepth(),
                            well_perforation_densities_, gravity);

                    xw.bhp()[w] = vfp_properties_.getInj()->bhp(vfp, aqua, liquid, vapour, thp) - dp;
                }
                else if (well_type == PRODUCER) {
                    // TODO: this needs to be updated
                    double dp = detail::computeHydrostaticCorrection(
                            wells(), w, vfp_properties_.getProd()->getTable(vfp)->getDatumDepth(),
                            well_perforation_densities_, gravity);

                    xw.bhp()[w] = vfp_properties_.getProd()->bhp(vfp, aqua, liquid, vapour, thp, alq) - dp;
                }
                else {
                    OPM_THROW(std::logic_error, "Expected PRODUCER or INJECTOR type of well");
                }
                break; */
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
                        // to make them consistent, it could perhaps improve convergence.
                        xw.segPhaseRates()[np * xw.topSegmentLoc()[w] + phase] = target * distr[phase];
                    }
                }
                break;
            }

        }
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

#if 0
        // Debug output.
        std::cout << "output bhp_targets " << std::endl;
        for (int i = 0; i < nw; ++i) {
            std::cout << i << "     " <<  bhp_targets(i) << std::endl;
        }
        std::cout << " bhp_well_elems " << std::endl;
        std::cout << " the size of bhp_well_elems is " << bhp_well_elems.size() << std::endl;
        for (size_t i = 0; i < bhp_well_elems.size(); ++i) {
            std::cout << " bhp_well_elems " << i << "   is " << bhp_well_elems[i] << std::endl;
        }
        std::cout << " rate_well_elems " << std::endl;
        std::cout << " the size of rate_well_elems is " << rate_well_elems.size() << std::endl;
        for (size_t i = 0; i < rate_well_elems.size(); ++i) {
            std::cout << " rate_well_elems " << i << " is " << rate_well_elems[i] << std::endl;
        }
        std::cout << " bhp_top_elems " << std::endl;
        std::cout << " the size of bhp_top_elems is " << bhp_top_elems.size() << std::endl;
        for (size_t i = 0; i < bhp_top_elems.size(); ++i) {
            std::cout << " bhp_top_elems " << i << "   is " << bhp_top_elems[i] << std::endl;
        }
        std::cout << " rate_top_elems " << std::endl;
        std::cout << " the size of the rate_top_elems " << rate_top_elems.size() << std::endl;
        for (size_t i = 0; i < rate_top_elems.size(); ++i) {
            std::cout << " rate_top_elems " << i << "   is " << rate_top_elems[i] << std::endl;
        }
        std::cout << " the others_elems " << std::endl;
        std::cout << " the size of others_elems is " << others_elems.size() << std::endl;
        for(size_t i = 0; i < others_elems.size(); ++i) {
            std::cout << "others_elems " << i << "   is " << others_elems[i] << std::endl;
        }

        std::cin.ignore();
#endif


        // for each segment: 1, if the segment is the top segment, then control equation
        //                   2, if the segment is not the top segment, then the pressure equation
        const ADB bhp_residual = subset(state.segp, bhp_top_elems) - subset(bhp_targets, bhp_well_elems);
        const ADB rate_residual = subset(rate_distr * subset(state.segqs, rate_top_phase_elems) - rate_targets, rate_well_elems);
        // const ADB rate_residual = subset(rate_distr * state.segqs - rate_targets, rate_well_elems);

        ADB others_residual = ADB::constant(V::Zero(nseg_total));
        start_segment = 0;
        for (int w = 0; w < nw; ++w) {
            WellMultiSegmentConstPtr well = wellsMultiSegment()[w];
            const int nseg = well->numberOfSegments();
            ADB segp = subset(state.segqs, Span(nseg, 1, start_segment));
            ADB well_residual = segp - well->wellOps().s2s_outlet * segp;
            ADB others_well_residual = subset(well_residual, Span(nseg - 1, 1, 1));
            others_residual = others_residual +  superset(others_well_residual, Span(nseg - 1, 1, start_segment + 1), nseg_total);
            start_segment += nseg;
        }

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

        // OPM_AD_DUMP(residual_.well_eq);

        // Calculate residuals
        // for (int w = 0; w < nw; ++w) {
        //     const int nseg = wellsMultiSegment()[w]->numberOfSegments();
        //     for (int s = 0; s < nseg; ++s) {
                // assuming the top segment always be the first one.
                // Three types of the pressure loss calculation
                // hydrostatic term depends of th density of the fluid mixture withn the segment,
                // TODO: as the first version, wo do not consider the rs rv in the mass flow rate and
                // fluid mixture density cacluation.
                // the FVF is based on the segment pressure
                // then the depth difference of the well segments
                // surface unit
                // frictional pressure drop
                // acceleration pressure drop
          //   }
        // }

        // const ADB thp_inj_residual = state.bhp - bhp_from_thp_inj + dp_inj;
        // const ADB thp_prod_residual = state.bhp - bhp_from_thp_prod + dp_prod;
        // const ADB bhp_residual = state.bhp - bhp_targets;
        // const ADB rate_residual = rate_distr * state.qs - rate_targets;

        // Select the right residual for each well
        // residual_.well_eq = superset(subset(bhp_residual, bhp_elems), bhp_elems, nw) +
        //         superset(subset(thp_inj_residual, thp_inj_elems), thp_inj_elems, nw) +
        //         superset(subset(thp_prod_residual, thp_prod_elems), thp_prod_elems, nw) +
        //         superset(subset(rate_residual, rate_elems), rate_elems, nw);

        // For wells that are dead (not flowing), and therefore not communicating
        // with the reservoir, we set the equation to be equal to the well's total
        // flow. This will be a solution only if the target rate is also zero.
        // M rate_summer(nw, np*nw);
        // for (int w = 0; w < nw; ++w) {
        //     for (int phase = 0; phase < np; ++phase) {
        //         rate_summer.insert(w, phase*nw + w) = 1.0;
        //     }
        // }

        // Selector<double> alive_selector(aliveWells, Selector<double>::NotEqualZero);
        // residual_.well_eq = alive_selector.select(residual_.well_eq, rate_summer * state.qs);
        // OPM_AD_DUMP(residual_.well_eq);
    }





    template <class Grid>
    void BlackoilMultiSegmentModel<Grid>::updateState(const V& dx,
                                          ReservoirState& reservoir_state,
                                          WellState& well_state)
    {
        using namespace Opm::AutoDiffGrid;
        const int np = fluid_.numPhases();
        const int nc = numCells(grid_);
        const V null;
        assert(null.size() == 0);
        const V zero = V::Zero(nc);

        // Extract parts of dx corresponding to each part.
        const V dp = subset(dx, Span(nc));
        int varstart = nc;
        const V dsw = active_[Water] ? subset(dx, Span(nc, 1, varstart)) : null;
        varstart += dsw.size();

        const V dxvar = active_[Gas] ? subset(dx, Span(nc, 1, varstart)): null;
        varstart += dxvar.size();

        // Extract well parts np phase rates + bhp
        const int nseg_total = well_state.numSegments();
        const V dwells = subset(dx, Span((np+1)*nseg_total, 1, varstart));
        varstart += dwells.size();

        assert(varstart == dx.size());

        // Pressure update.
        const double dpmaxrel = dpMaxRel();
        const V p_old = Eigen::Map<const V>(&reservoir_state.pressure()[0], nc, 1);
        const V absdpmax = dpmaxrel*p_old.abs();
        const V dp_limited = sign(dp) * dp.abs().min(absdpmax);
        const V p = (p_old - dp_limited).max(zero);
        std::copy(&p[0], &p[0] + nc, reservoir_state.pressure().begin());


        // Saturation updates.
        const Opm::PhaseUsage& pu = fluid_.phaseUsage();
        const DataBlock s_old = Eigen::Map<const DataBlock>(& reservoir_state.saturation()[0], nc, np);
        const double dsmax = dsMax();

        V so;
        V sw;
        V sg;

        {
            V maxVal = zero;
            V dso = zero;
            if (active_[Water]){
                maxVal = dsw.abs().max(maxVal);
                dso = dso - dsw;
            }

            V dsg;
            if (active_[Gas]){
                dsg = isSg_ * dxvar - isRv_ * dsw;
                maxVal = dsg.abs().max(maxVal);
                dso = dso - dsg;
            }

            maxVal = dso.abs().max(maxVal);

            V step = dsmax/maxVal;
            step = step.min(1.);

            if (active_[Water]) {
                const int pos = pu.phase_pos[ Water ];
                const V sw_old = s_old.col(pos);
                sw = sw_old - step * dsw;
            }

            if (active_[Gas]) {
                const int pos = pu.phase_pos[ Gas ];
                const V sg_old = s_old.col(pos);
                sg = sg_old - step * dsg;
            }

            const int pos = pu.phase_pos[ Oil ];
            const V so_old = s_old.col(pos);
            so = so_old - step * dso;
        }

        // Appleyard chop process.
        auto ixg = sg < 0;
        for (int c = 0; c < nc; ++c) {
            if (ixg[c]) {
                sw[c] = sw[c] / (1-sg[c]);
                so[c] = so[c] / (1-sg[c]);
                sg[c] = 0;
            }
        }


        auto ixo = so < 0;
        for (int c = 0; c < nc; ++c) {
            if (ixo[c]) {
                sw[c] = sw[c] / (1-so[c]);
                sg[c] = sg[c] / (1-so[c]);
                so[c] = 0;
            }
        }

        auto ixw = sw < 0;
        for (int c = 0; c < nc; ++c) {
            if (ixw[c]) {
                so[c] = so[c] / (1-sw[c]);
                sg[c] = sg[c] / (1-sw[c]);
                sw[c] = 0;
            }
        }

        //const V sumSat = sw + so + sg;
        //sw = sw / sumSat;
        //so = so / sumSat;
        //sg = sg / sumSat;

        // Update the reservoir_state
        for (int c = 0; c < nc; ++c) {
            reservoir_state.saturation()[c*np + pu.phase_pos[ Water ]] = sw[c];
        }

        for (int c = 0; c < nc; ++c) {
            reservoir_state.saturation()[c*np + pu.phase_pos[ Gas ]] = sg[c];
        }

        if (active_[ Oil ]) {
            const int pos = pu.phase_pos[ Oil ];
            for (int c = 0; c < nc; ++c) {
                reservoir_state.saturation()[c*np + pos] = so[c];
            }
        }

        // Update rs and rv
        const double drmaxrel = drMaxRel();
        V rs;
        if (has_disgas_) {
            const V rs_old = Eigen::Map<const V>(&reservoir_state.gasoilratio()[0], nc);
            const V drs = isRs_ * dxvar;
            const V drs_limited = sign(drs) * drs.abs().min(rs_old.abs()*drmaxrel);
            rs = rs_old - drs_limited;
        }
        V rv;
        if (has_vapoil_) {
            const V rv_old = Eigen::Map<const V>(&reservoir_state.rv()[0], nc);
            const V drv = isRv_ * dxvar;
            const V drv_limited = sign(drv) * drv.abs().min(rv_old.abs()*drmaxrel);
            rv = rv_old - drv_limited;
        }

        // Sg is used as primal variable for water only cells.
        const double epsilon = std::sqrt(std::numeric_limits<double>::epsilon());
        auto watOnly = sw >  (1 - epsilon);

        // phase translation sg <-> rs
        std::fill(primalVariable_.begin(), primalVariable_.end(), PrimalVariables::Sg);

        if (has_disgas_) {
            const V rsSat0 = fluidRsSat(p_old, s_old.col(pu.phase_pos[Oil]), cells_);
            const V rsSat = fluidRsSat(p, so, cells_);
            // The obvious case
            auto hasGas = (sg > 0 && isRs_ == 0);

            // Set oil saturated if previous rs is sufficiently large
            const V rs_old = Eigen::Map<const V>(&reservoir_state.gasoilratio()[0], nc);
            auto gasVaporized =  ( (rs > rsSat * (1+epsilon) && isRs_ == 1 ) && (rs_old > rsSat0 * (1-epsilon)) );
            auto useSg = watOnly || hasGas || gasVaporized;
            for (int c = 0; c < nc; ++c) {
                if (useSg[c]) {
                    rs[c] = rsSat[c];
                } else {
                    primalVariable_[c] = PrimalVariables::RS;
                }
            }

        }

        // phase transitions so <-> rv
        if (has_vapoil_) {

            // The gas pressure is needed for the rvSat calculations
            const V gaspress_old = computeGasPressure(p_old, s_old.col(Water), s_old.col(Oil), s_old.col(Gas));
            const V gaspress = computeGasPressure(p, sw, so, sg);
            const V rvSat0 = fluidRvSat(gaspress_old, s_old.col(pu.phase_pos[Oil]), cells_);
            const V rvSat = fluidRvSat(gaspress, so, cells_);

            // The obvious case
            auto hasOil = (so > 0 && isRv_ == 0);

            // Set oil saturated if previous rv is sufficiently large
            const V rv_old = Eigen::Map<const V>(&reservoir_state.rv()[0], nc);
            auto oilCondensed = ( (rv > rvSat * (1+epsilon) && isRv_ == 1) && (rv_old > rvSat0 * (1-epsilon)) );
            auto useSg = watOnly || hasOil || oilCondensed;
            for (int c = 0; c < nc; ++c) {
                if (useSg[c]) {
                    rv[c] = rvSat[c];
                } else {
                    primalVariable_[c] = PrimalVariables::RV;
                }
            }

        }

        // Update the reservoir_state
        if (has_disgas_) {
            std::copy(&rs[0], &rs[0] + nc, reservoir_state.gasoilratio().begin());
        }

        if (has_vapoil_) {
            std::copy(&rv[0], &rv[0] + nc, reservoir_state.rv().begin());
        }


        updateWellState(dwells,well_state);

        // Update phase conditions used for property calculations.
        updatePhaseCondFromPrimalVariable();
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
#if 0
            std::cout << " output sr in updateWellState " << std::endl;
            std::cout << sr << std::endl;
            std::cin.ignore();

            std::cout << " output bhp is updateWellState " << std::endl;
            std::cout << segp << std::endl;
            std::cin.ignore();
#endif


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

#if 0
        // Debug output.
        std::cout << " output all the well state informations after updateWellState()" << std::endl;
        const int np = well_state.numPhases();
        const int nw = well_state.numWells();
        const int nperf_total = well_state.numPerforations();
        const int nseg_total = well_state.numSegments();

        std::cout << " number of wells : " << nw << " nubmer of segments : " << nseg_total << std::endl;
        std::cout << " number of phase : " << np << " nubmer of perforations " << nperf_total << std::endl;

        std::cout << " bhps : " << std::endl;
        for (int i = 0; i < nw; ++i) {
            std::cout << well_state.bhp()[i] << std::endl;
        }

        std::cout << " thps : " << std::endl;

        for (int i = 0; i < nw; ++i) {
            std::cout << well_state.thp()[i] << std::endl;
        }

        std::cout << " well rates " << std::endl;
        for (int i = 0; i < nw; ++i) {
            std::cout << i;
            for (int p = 0; p < np; ++p) {
                std::cout << " " << well_state.wellRates()[np * i + p];
            }
            std::cout << std::endl;
        }

        std::cout << " segment pressures and segment phase rates : " << std::endl;
        for (int i = 0; i < nseg_total; ++i) {
            std::cout << i << " " << well_state.segPress()[i];
            for (int p = 0; p < np; ++p) {
                std::cout << " " << well_state.segPhaseRates()[np * i + p];
            }
            std::cout << std::endl;
        }

        std::cout << " perf pressures and pref phase rates : " << std::endl;
        for (int i = 0; i < nperf_total; ++i) {
            std::cout << i << " " << well_state.perfPress()[i];
            for (int p = 0; p < np; ++p) {
                std::cout << " " << well_state.perfPhaseRates()[np * i + p];
            }
            std::cout << std::endl;
        }
        std::cout << " output all the well state informations after updateWellState() DONE!!!!" << std::endl;
#endif
    }







    template <class Grid>
    bool
    BlackoilMultiSegmentModel<Grid>::getConvergence(const double dt, const int iteration)
    {

        const double tol_mb    = param_.tolerance_mb_;
        const double tol_cnv   = param_.tolerance_cnv_;
        const double tol_wells = param_.tolerance_wells_;

        const int nc = Opm::AutoDiffGrid::numCells(grid_);
        // const int nw = localWellsActive() ? wells().number_of_wells : 0;
        const int nw = wellsMultiSegment().size();
        // no good way to store nseg?
        int nseg_total = 0;
        for (int w = 0; w < nw; ++w) {
            nseg_total += wellsMultiSegment()[w]->numberOfSegments();
        }
        const int np = numPhases();
        const int nm = numMaterials();
        assert(int(rq_.size()) == nm);

        const V pv = geo_.poreVolume();

        std::vector<double> R_sum(nm);
        std::vector<double> B_avg(nm);
        std::vector<double> maxCoeff(nm);
        std::vector<double> maxNormWell(np);
        Eigen::Array<V::Scalar, Eigen::Dynamic, Eigen::Dynamic> B(nc, nm);
        Eigen::Array<V::Scalar, Eigen::Dynamic, Eigen::Dynamic> R(nc, nm);
        Eigen::Array<V::Scalar, Eigen::Dynamic, Eigen::Dynamic> tempV(nc, nm);

        for ( int idx = 0; idx < nm; ++idx )
        {
            const ADB& tempB = rq_[idx].b;
            B.col(idx)       = 1./tempB.value();
            R.col(idx)       = residual_.material_balance_eq[idx].value();
            tempV.col(idx)   = R.col(idx).abs()/pv;
        }

        const double pvSum = convergenceReduction(B, tempV, R,
                                                  R_sum, maxCoeff, B_avg, maxNormWell,
                                                  nc, nseg_total);

        std::vector<double> CNV(nm);
        std::vector<double> mass_balance_residual(nm);
        std::vector<double> well_flux_residual(np);

        bool converged_MB = true;
        bool converged_CNV = true;
        bool converged_Well = true;
        // Finish computation
        for ( int idx = 0; idx < nm; ++idx )
        {
            CNV[idx]                    = B_avg[idx] * dt * maxCoeff[idx];
            mass_balance_residual[idx]  = std::abs(B_avg[idx]*R_sum[idx]) * dt / pvSum;
            converged_MB                = converged_MB && (mass_balance_residual[idx] < tol_mb);
            converged_CNV               = converged_CNV && (CNV[idx] < tol_cnv);
            // Well flux convergence is only for fluid phases, not other materials
            // in our current implementation.
            assert(nm >= np);
            if (idx < np) {
                well_flux_residual[idx] = B_avg[idx] * maxNormWell[idx];
                converged_Well = converged_Well && (well_flux_residual[idx] < tol_wells);
            }
        }

        const double residualWell     = detail::infinityNormWell(residual_.well_eq,
                                                                 linsolver_.parallelInformation());
        converged_Well = converged_Well && (residualWell < Opm::unit::barsa);
        const bool converged = converged_MB && converged_CNV && converged_Well;

        // Residual in Pascal can have high values and still be ok.
        const double maxWellResidualAllowed = 1000.0 * maxResidualAllowed();

        for (int idx = 0; idx < nm; ++idx) {
            if (std::isnan(mass_balance_residual[idx])
                || std::isnan(CNV[idx])
                || (idx < np && std::isnan(well_flux_residual[idx]))) {
                OPM_THROW(Opm::NumericalProblem, "NaN residual for phase " << materialName(idx));
            }
            if (mass_balance_residual[idx] > maxResidualAllowed()
                || CNV[idx] > maxResidualAllowed()
                || (idx < np && well_flux_residual[idx] > maxResidualAllowed())) {
                OPM_THROW(Opm::NumericalProblem, "Too large residual for phase " << materialName(idx));
            }
        }
        if (std::isnan(residualWell) || residualWell > maxWellResidualAllowed) {
            OPM_THROW(Opm::NumericalProblem, "NaN or too large residual for well control equation");
        }

        if ( terminal_output_ )
        {
            // Only rank 0 does print to std::cout
            if (iteration == 0) {
                std::cout << "\nIter";
                for (int idx = 0; idx < nm; ++idx) {
                    std::cout << "   MB(" << materialName(idx).substr(0, 3) << ") ";
                }
                for (int idx = 0; idx < nm; ++idx) {
                    std::cout << "    CNV(" << materialName(idx).substr(0, 1) << ") ";
                }
                for (int idx = 0; idx < np; ++idx) {
                    std::cout << "  W-FLUX(" << materialName(idx).substr(0, 1) << ")";
                }
                std::cout << '\n';
            }
            const std::streamsize oprec = std::cout.precision(3);
            const std::ios::fmtflags oflags = std::cout.setf(std::ios::scientific);
            std::cout << std::setw(4) << iteration;
            for (int idx = 0; idx < nm; ++idx) {
                std::cout << std::setw(11) << mass_balance_residual[idx];
            }
            for (int idx = 0; idx < nm; ++idx) {
                std::cout << std::setw(11) << CNV[idx];
            }
            for (int idx = 0; idx < np; ++idx) {
                std::cout << std::setw(11) << well_flux_residual[idx];
            }
            std::cout << std::endl;
            std::cout.precision(oprec);
            std::cout.flags(oflags);
        }
        return converged;
    }





    template <class Grid>
    void
    BlackoilMultiSegmentModel<Grid>::computeSegmentDensities(const SolutionState& /* state */,
                                                             const std::vector<ADB>& cq_s,
                                                             const std::vector<ADB>& b_perf,
                                                             const WellState& xw)
    {
        // how to calculate the formation valume factor for the segments is the most tricky part.
        // ave_rho = (\sigma(rho_surf * Q_surf)) / (\sigma(rho_surf * B))
        // how to calculate the rs and rv to obtain B
        // should we relate the segments with different cells
        // TODO: as the first solution, we calculate the rs and rv as the average rs and rv from the related perforations.
        // Based on the current framework, it has to be done one well by one well
        const int nw = xw.numWells();
        const int nseg_total = xw.numSegments();
        const int np = numPhases();

        assert(np == int(b_perf.size()));

        well_segment_densities_ = ADB::constant(V::Zero(nseg_total)); // initialize to be zero

        int start_segment = 0;
        int start_perforation = 0;

        // surface density
        std::vector<double> surf_dens(fluid_.surfaceDensity(), fluid_.surfaceDensity() + np);

        for (int w = 0; w < nw; ++w) {
            // TODO: to check to make sure b_perf is B or 1/B
            // TODO: now for non-segmented wells, we do not calculate the segment densities.
            //       If the WELSPECS set the AVG, the density for the wellbore can be same with the segment density.
            // TODO: which is one way to unify the implementation
            WellMultiSegmentConstPtr well = wellsMultiSegment()[w];
            const int nseg = well->numberOfSegments();
            const int nperf = well->numberOfPerforations();

            if (well->isMultiSegmented()) {
                // calculate the the formation volume factor.
                // calcuate the rs and rv first
                // at the moment, the rs and rv are the averaged rs and rv for the perforation cells related to the segment
                std::vector<ADB> b_perf_well(np, ADB::null());
                std::vector<ADB> cq_s_well(np, ADB::null());
                std::vector<ADB> b_segment(np, ADB::null());

                // accumlated phase rates for segments from the perforations
                std::vector<ADB> q_ps(np, ADB::null());
                std::vector<ADB> b_seg(np, ADB::null());
                for (int p = 0; p < np; ++p) {
                    cq_s_well[p] = subset(cq_s[p], Span(nperf, 1, start_perforation));
                    b_perf_well[p] = subset(b_perf[p], Span(nperf, 1, start_perforation));

                    q_ps[p] = well->wellOps().p2s_gather * cq_s_well[p];
                    b_seg[p] = well->wellOps().p2s_average * b_perf_well[p];
                }

                ADB density = ADB::constant(V::Zero(nseg));

                ADB temp_rho_q = ADB::constant(V::Zero(nseg));
                ADB temp_volume_rate = ADB::constant(V::Zero(nseg));

                for (int p = 0; p < np; ++p) {
                    temp_rho_q += q_ps[p] * surf_dens[p];
                    temp_volume_rate += q_ps[p] * b_seg[p];
                }

                density = temp_rho_q / temp_volume_rate;

                well_segment_densities_ += superset(density, Span(nseg, 1, start_segment), nseg_total);
                // formation factor volume B should be calculated based on the segment pressures
                // how to handle rs and rv to calculate B
                // ADB& seg_press = subset(state.segp, Span(nseg, 1, start_segment));

            }
            start_segment += nseg;
            start_perforation += nperf;
        }
    }



} // namespace Opm

#endif // OPM_BLACKOILMODELBASE_IMPL_HEADER_INCLUDED
