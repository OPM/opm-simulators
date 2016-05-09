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


    template <class Grid>
    BlackoilMultiSegmentModel<Grid>::
    BlackoilMultiSegmentModel(const typename Base::ModelParameters&  param,
                  const Grid&                     grid ,
                  const BlackoilPropsAdInterface& fluid,
                  const DerivedGeology&           geo  ,
                  const RockCompressibility*      rock_comp_props,
                  const MultisegmentWells&        well_model,
                  const NewtonIterationBlackoilInterface&    linsolver,
                  Opm::EclipseStateConstPtr eclState,
                  const bool has_disgas,
                  const bool has_vapoil,
                  const bool terminal_output)
        : Base(param, grid, fluid, geo, rock_comp_props, well_model, linsolver,
               eclState, has_disgas, has_vapoil, terminal_output)
        , ms_wells_(well_model)
    {
        const double gravity = detail::getGravity(geo_.gravity(), UgGridHelpers::dimensions(grid_));
        const V depth = Opm::AutoDiffGrid::cellCentroidsZToEigen(grid_);

        ms_wells_.init(&fluid_, &active_, &phaseCondition_, &vfp_properties_, gravity, depth);
        // TODO: there should be a better way do the following
        ms_wells_.setWellsActive(Base::wellsActive());
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

        const int nw = wellsMultiSegment().size();

        if ( !msWellOps().has_multisegment_wells ) {
            msWells().segVDt() = V::Zero(nw);
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
        msWells().segVDt() = Eigen::Map<V>(segment_volume.data(), nseg_total) / dt;
    }





    template <class Grid>
    int
    BlackoilMultiSegmentModel<Grid>::numWellVars() const
    {
        // For each segment, we have a pressure variable, and one flux per phase.
        const int nseg = msWellOps().p2s.rows();
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
        const int nw = msWells().topWellSegments().size();
        state.qs = ADB::constant(ADB::V::Zero(np*nw));
        for (int phase = 0; phase < np; ++phase) {
            // Extract segment fluxes for this phase (ns consecutive elements).
            ADB segqs_phase = subset(state.segqs, Span(ns, 1, ns*phase));
            // Extract top segment fluxes (= well fluxes)
            ADB wellqs_phase = subset(segqs_phase, msWells().topWellSegments());
            // Expand to full size of qs (which contains all phases) and add.
            state.qs += superset(wellqs_phase, Span(nw, 1, nw*phase), nw*np);
        }
        state.bhp = subset(state.segp, msWells().topWellSegments());
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
        msWells().updateWellControls(terminal_output_, well_state);

        // Create the primary variables.
        SolutionState state = asImpl().variableState(reservoir_state, well_state);

        if (initial_assembly) {
            // Create the (constant, derivativeless) initial state.
            SolutionState state0 = state;
            asImpl().makeConstantState(state0);
            // Compute initial accumulation contributions
            // and well connection pressures.
            asImpl().computeAccum(state0, 0);
            msWells().computeSegmentFluidProperties(state0);
            const int np = numPhases();
            assert(np == int(msWells().segmentCompSurfVolumeInitial().size()));
            for (int phase = 0; phase < np; ++phase) {
                msWells().segmentCompSurfVolumeInitial()[phase] = msWells().segmentCompSurfVolumeCurrent()[phase].value();
            }

            const std::vector<ADB> kr_adb = Base::computeRelPerm(state0);
            std::vector<ADB> fluid_density(numPhases(), ADB::null());
            // TODO: make sure the order of the density and the order of the kr are the same.
            for (int phaseIdx = 0; phaseIdx < fluid_.numPhases(); ++phaseIdx) {
                const int canonicalPhaseIdx = canph_[phaseIdx];
                fluid_density[phaseIdx] = fluidDensity(canonicalPhaseIdx, rq_[phaseIdx].b, state0.rs, state0.rv);
            }
            msWells().computeWellConnectionPressures(state0, well_state, kr_adb, fluid_density);
            // asImpl().computeWellConnectionPressures(state0, well_state);
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

        // asImpl().computeSegmentFluidProperties(state);
        msWells().computeSegmentFluidProperties(state);

        // asImpl().computeSegmentPressuresDelta(state);
        const double gravity = detail::getGravity(geo_.gravity(), UgGridHelpers::dimensions(grid_));
        msWells().computeSegmentPressuresDelta(gravity);

        std::vector<ADB> mob_perfcells;
        std::vector<ADB> b_perfcells;
        msWells().extractWellPerfProperties(state, rq_, mob_perfcells, b_perfcells);
        if (param_.solve_welleq_initially_ && initial_assembly) {
            // solve the well equations as a pre-processing step
            asImpl().solveWellEq(mob_perfcells, b_perfcells, state, well_state);
        }

        // the perforation flux here are different
        // it is related to the segment location
        V aliveWells;
        std::vector<ADB> cq_s;
        msWells().computeWellFlux(state, mob_perfcells, b_perfcells, aliveWells, cq_s);
        msWells().updatePerfPhaseRatesAndPressures(cq_s, state, well_state);
        msWells().addWellFluxEq(cq_s, state, residual_);
        asImpl().addWellContributionToMassBalanceEq(cq_s, state, well_state);
        // asImpl().addWellControlEq(state, well_state, aliveWells);
        msWells().addWellControlEq(state, well_state, aliveWells, residual_);
    }





    template <class Grid>
    bool BlackoilMultiSegmentModel<Grid>::solveWellEq(const std::vector<ADB>& mob_perfcells,
                                                      const std::vector<ADB>& b_perfcells,
                                                      SolutionState& state,
                                                      WellState& well_state)
    {
        const bool converged = baseSolveWellEq(mob_perfcells, b_perfcells, state, well_state);

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
            const std::vector<ADB> kr_adb = Base::computeRelPerm(state);
            std::vector<ADB> fluid_density(np, ADB::null());
            // TODO: make sure the order of the density and the order of the kr are the same.
            for (int phaseIdx = 0; phaseIdx < np; ++phaseIdx) {
                const int canonicalPhaseIdx = canph_[phaseIdx];
                fluid_density[phaseIdx] = fluidDensity(canonicalPhaseIdx, rq_[phaseIdx].b, state.rs, state.rv);
            }
            msWells().computeWellConnectionPressures(state, well_state, kr_adb, fluid_density);
            // asImpl().computeWellConnectionPressures(state, well_state);
        }

        return converged;
    }





    template <class Grid>
    void
    BlackoilMultiSegmentModel<Grid>::updateWellState(const V& dwells,
                                                     WellState& well_state)
    {
        msWells().updateWellState(dwells, dpMaxRel(), well_state);
    }





        /// added to fixing the flow_multisegment running
    template <class Grid>
    bool
    BlackoilMultiSegmentModel<Grid>::baseSolveWellEq(const std::vector<ADB>& mob_perfcells,
                                                     const std::vector<ADB>& b_perfcells,
                                                     SolutionState& state,
                                                     WellState& well_state) {
        V aliveWells;
        const int np = msWells().numPhases();
        std::vector<ADB> cq_s(np, ADB::null());
        std::vector<int> indices = msWells().variableWellStateIndices();
        SolutionState state0 = state;
        WellState well_state0 = well_state;
        makeConstantState(state0);

        std::vector<ADB> mob_perfcells_const(np, ADB::null());
        std::vector<ADB> b_perfcells_const(np, ADB::null());

        if ( Base::localWellsActive() ){
            // If there are non well in the sudomain of the process
            // thene mob_perfcells_const and b_perfcells_const would be empty
            for (int phase = 0; phase < np; ++phase) {
                mob_perfcells_const[phase] = ADB::constant(mob_perfcells[phase].value());
                b_perfcells_const[phase] = ADB::constant(b_perfcells[phase].value());
            }
        }

        int it  = 0;
        bool converged;
        do {
            // bhp and Q for the wells
            std::vector<V> vars0;
            vars0.reserve(2);
            msWells().variableWellStateInitials(well_state, vars0);
            std::vector<ADB> vars = ADB::variables(vars0);

            SolutionState wellSolutionState = state0;
            variableStateExtractWellsVars(indices, vars, wellSolutionState);

            msWells().computeWellFlux(wellSolutionState, mob_perfcells_const, b_perfcells_const, aliveWells, cq_s);

            msWells().updatePerfPhaseRatesAndPressures(cq_s, wellSolutionState, well_state);
            msWells().addWellFluxEq(cq_s, wellSolutionState, residual_);
            // addWellControlEq(wellSolutionState, well_state, aliveWells);
            msWells().addWellControlEq(wellSolutionState, well_state, aliveWells, residual_);
            converged = Base::getWellConvergence(it);

            if (converged) {
                break;
            }

            ++it;
            if( Base::localWellsActive() )
            {
                std::vector<ADB> eqs;
                eqs.reserve(2);
                eqs.push_back(residual_.well_flux_eq);
                eqs.push_back(residual_.well_eq);
                ADB total_residual = vertcatCollapseJacs(eqs);
                const std::vector<M>& Jn = total_residual.derivative();
                typedef Eigen::SparseMatrix<double> Sp;
                Sp Jn0;
                Jn[0].toSparse(Jn0);
                const Eigen::SparseLU< Sp > solver(Jn0);
                ADB::V total_residual_v = total_residual.value();
                const Eigen::VectorXd& dx = solver.solve(total_residual_v.matrix());
                assert(dx.size() == total_residual_v.size());
                asImpl().updateWellState(dx.array(), well_state);
                msWells().updateWellControls(terminal_output_, well_state);
            }
        } while (it < 15);

        if (converged) {
            if ( terminal_output_ ) {
                std::cout << "well converged iter: " << it << std::endl;
            }
            const int nw = msWells().numWells();
            {
                // We will set the bhp primary variable to the new ones,
                // but we do not change the derivatives here.
                ADB::V new_bhp = Eigen::Map<ADB::V>(well_state.bhp().data(), nw);
                // Avoiding the copy below would require a value setter method
                // in AutoDiffBlock.
                std::vector<ADB::M> old_derivs = state.bhp.derivative();
                state.bhp = ADB::function(std::move(new_bhp), std::move(old_derivs));
            }
            {
                // Need to reshuffle well rates, from phase running fastest
                // to wells running fastest.
                // The transpose() below switches the ordering.
                const DataBlock wrates = Eigen::Map<const DataBlock>(well_state.wellRates().data(), nw, np).transpose();
                ADB::V new_qs = Eigen::Map<const V>(wrates.data(), nw*np);
                std::vector<ADB::M> old_derivs = state.qs.derivative();
                state.qs = ADB::function(std::move(new_qs), std::move(old_derivs));
            }

            const std::vector<ADB> kr_adb = Base::computeRelPerm(state);
            std::vector<ADB> fluid_density(np, ADB::null());
            // TODO: make sure the order of the density and the order of the kr are the same.
            for (int phaseIdx = 0; phaseIdx < np; ++phaseIdx) {
                const int canonicalPhaseIdx = canph_[phaseIdx];
                fluid_density[phaseIdx] = fluidDensity(canonicalPhaseIdx, rq_[phaseIdx].b, state.rs, state.rv);
            }
            msWells().computeWellConnectionPressures(state, well_state, kr_adb, fluid_density);
            // computeWellConnectionPressures(state, well_state);
        }

        if (!converged) {
            well_state = well_state0;
        }

        return converged;
    }





    template <class Grid>
    std::vector<V>
    BlackoilMultiSegmentModel<Grid>::
    variableStateInitials(const ReservoirState& x,
                          const WellState&     xw) const
    {
        assert(active_[ Oil ]);

        const int np = x.numPhases();

        std::vector<V> vars0;
        // p, Sw and Rs, Rv or Sg is used as primary depending on solution conditions
        // and bhp and Q for the wells
        vars0.reserve(np + 1);
        variableReservoirStateInitials(x, vars0);
        msWells().variableWellStateInitials(xw, vars0);
        return vars0;
    }

} // namespace Opm

#endif // OPM_BLACKOILMODELBASE_IMPL_HEADER_INCLUDED
