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
#include <opm/autodiff/BlackoilPropsAdFromDeck.hpp>
#include <opm/autodiff/GeoProps.hpp>
#include <opm/autodiff/WellDensitySegmented.hpp>
#include <opm/autodiff/VFPPropertiesAdb.hpp>
#include <opm/autodiff/VFPProdPropertiesAdb.hpp>
#include <opm/autodiff/VFPInjPropertiesAdb.hpp>

#include <opm/core/grid.h>
#include <opm/core/linalg/LinearSolverInterface.hpp>
#include <opm/core/linalg/ParallelIstlInformation.hpp>
#include <opm/core/props/rock/RockCompressibility.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>
#include <opm/parser/eclipse/Units/Units.hpp>
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
                  const BlackoilPropsAdFromDeck& fluid,
                  const DerivedGeology&           geo  ,
                  const RockCompressibility*      rock_comp_props,
                  const MultisegmentWells&        well_model,
                  const NewtonIterationBlackoilInterface&    linsolver,
                  std::shared_ptr< const EclipseState > eclState,
                  std::shared_ptr<const Schedule> schedule,
                  std::shared_ptr<const SummaryConfig> summary_config,
                  const bool has_disgas,
                  const bool has_vapoil,
                  const bool terminal_output)
        : Base(param, grid, fluid, geo, rock_comp_props, well_model, linsolver,
               eclState, schedule, summary_config, has_disgas, has_vapoil, terminal_output)
    {
    }





    template <class Grid>
    void
    BlackoilMultiSegmentModel<Grid>::
    prepareStep(const SimulatorTimerInterface& timer,
                const ReservoirState& reservoir_state,
                const WellState& well_state)
    {
        const double dt = timer.currentStepLength();
        pvdt_ = geo_.poreVolume() / dt;
        if (active_[Gas]) {
            updatePrimalVariableFromState(reservoir_state);
        }

        const int nw = wellsMultiSegment().size();

        if ( !msWellOps().has_multisegment_wells ) {
            wellModel().segVDt() = V::Zero(nw);
            return;
        }

        const int nseg_total = well_state.numSegments();
        std::vector<double> segment_volume;
        segment_volume.reserve(nseg_total);
        for (int w = 0; w < nw; ++w) {
            WellMultiSegmentConstPtr well = wellsMultiSegment()[w];
            const std::vector<double>& segment_volume_well = well->segmentVolume();
            segment_volume.insert(segment_volume.end(), segment_volume_well.begin(), segment_volume_well.end());
        }
        assert(int(segment_volume.size()) == nseg_total);
        wellModel().segVDt() = Eigen::Map<V>(segment_volume.data(), nseg_total) / dt;
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
    SimulatorReport
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
        wellModel().updateWellControls(well_state);

        // TODO: I do not think the multi_segment well can handle group control yet
        if (asImpl().wellModel().wellCollection()->groupControlActive()) {
            // enforce VREP control when necessary.
            Base::applyVREPGroupControl(reservoir_state, well_state);

            asImpl().wellModel().wellCollection()->updateWellTargets(well_state.wellRates());
        }

        // Create the primary variables.
        SolutionState state = asImpl().variableState(reservoir_state, well_state);

        if (initial_assembly) {
            // Create the (constant, derivativeless) initial state.
            SolutionState state0 = state;
            asImpl().makeConstantState(state0);
            // Compute initial accumulation contributions
            // and well connection pressures.
            asImpl().computeAccum(state0, 0);
            wellModel().computeSegmentFluidProperties(state0);
            const int np = numPhases();
            assert(np == int(wellModel().segmentCompSurfVolumeInitial().size()));
            for (int phase = 0; phase < np; ++phase) {
                wellModel().segmentCompSurfVolumeInitial()[phase] = wellModel().segmentCompSurfVolumeCurrent()[phase].value();
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
            SimulatorReport report;
            return report;
        }

        wellModel().computeSegmentFluidProperties(state);

        const double gravity = detail::getGravity(geo_.gravity(), UgGridHelpers::dimensions(grid_));
        wellModel().computeSegmentPressuresDelta(gravity);

        std::vector<ADB> mob_perfcells;
        std::vector<ADB> b_perfcells;
        SimulatorReport report;
        wellModel().extractWellPerfProperties(state, sd_.rq, mob_perfcells, b_perfcells);
        if (param_.solve_welleq_initially_ && initial_assembly) {
            // solve the well equations as a pre-processing step
            report = asImpl().solveWellEq(mob_perfcells, b_perfcells, reservoir_state, state, well_state);
        }

        // the perforation flux here are different
        // it is related to the segment location
        V aliveWells;
        std::vector<ADB> cq_s;
        wellModel().computeWellFlux(state, mob_perfcells, b_perfcells, aliveWells, cq_s);
        wellModel().updatePerfPhaseRatesAndPressures(cq_s, state, well_state);
        wellModel().addWellFluxEq(cq_s, state, residual_);
        asImpl().addWellContributionToMassBalanceEq(cq_s, state, well_state);
        wellModel().addWellControlEq(state, well_state, aliveWells, residual_);
        return report;
    }





    template <class Grid>
    SimulatorReport
    BlackoilMultiSegmentModel<Grid>::solveWellEq(const std::vector<ADB>& mob_perfcells,
                                                 const std::vector<ADB>& b_perfcells,
                                                 const ReservoirState& reservoir_state,
                                                 SolutionState& state,
                                                 WellState& well_state)
    {
        SimulatorReport report = Base::solveWellEq(mob_perfcells, b_perfcells, reservoir_state, state, well_state);

        if (report.converged) {
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

        return report;
    }




    template <class Grid>
    void
    BlackoilMultiSegmentModel<Grid>::
    computeWellConnectionPressures(const SolutionState& state,
                                   const WellState& well_state)
    {
        const int np = numPhases();
        const std::vector<ADB> kr_adb = Base::computeRelPerm(state);
        std::vector<ADB> fluid_density(np, ADB::null());
        // TODO: make sure the order of the density and the order of the kr are the same.
        for (int phaseIdx = 0; phaseIdx < np; ++phaseIdx) {
            const int canonicalPhaseIdx = canph_[phaseIdx];
            fluid_density[phaseIdx] = fluidDensity(canonicalPhaseIdx, sd_.rq[phaseIdx].b, state.rs, state.rv);
         }
         wellModel().computeWellConnectionPressures(state, well_state, kr_adb, fluid_density);
    }

} // namespace Opm

#endif // OPM_BLACKOILMODELBASE_IMPL_HEADER_INCLUDED
