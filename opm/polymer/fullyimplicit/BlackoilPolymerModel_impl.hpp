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

#ifndef OPM_BLACKOILPOLYMERMODEL_IMPL_HEADER_INCLUDED
#define OPM_BLACKOILPOLYMERMODEL_IMPL_HEADER_INCLUDED

#include <opm/polymer/fullyimplicit/BlackoilPolymerModel.hpp>

#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/autodiff/GridHelpers.hpp>
#include <opm/autodiff/BlackoilPropsAdInterface.hpp>
#include <opm/autodiff/GeoProps.hpp>
#include <opm/autodiff/WellDensitySegmented.hpp>

#include <opm/core/grid.h>
#include <opm/core/linalg/LinearSolverInterface.hpp>
#include <opm/core/linalg/ParallelIstlInformation.hpp>
#include <opm/core/props/rock/RockCompressibility.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/well_controls.h>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>

namespace Opm {



    namespace detail {

        template <class PU>
        int polymerPos(const PU& pu)
        {
            const int maxnp = Opm::BlackoilPhases::MaxNumPhases;
            int pos = 0;
            for (int phase = 0; phase < maxnp; ++phase) {
                if (pu.phase_used[phase]) {
                    pos++;
                }
            }

            return pos;
        }

    } // namespace detail



    template <class Grid>
    BlackoilPolymerModel<Grid>::BlackoilPolymerModel(const typename Base::ModelParameters&   param,
                                                     const Grid&                             grid,
                                                     const BlackoilPropsAdInterface&         fluid,
                                                     const DerivedGeology&                   geo,
                                                     const RockCompressibility*              rock_comp_props,
                                                     const PolymerPropsAd&                   polymer_props_ad,
                                                     const Wells*                            wells,
                                                     const NewtonIterationBlackoilInterface& linsolver,
                                                     EclipseStateConstPtr                    eclipse_state,
                                                     const bool                              has_disgas,
                                                     const bool                              has_vapoil,
                                                     const bool                              has_polymer,
                                                     const bool                              has_plyshlog,
                                                     const bool                              has_shrate,
                                                     const std::vector<double>&              wells_rep_radius,
                                                     const std::vector<double>&              wells_perf_length,
                                                     const std::vector<double>&              wells_bore_diameter,
                                                     const bool                              terminal_output)
        : Base(param, grid, fluid, geo, rock_comp_props, wells, linsolver, eclipse_state,
               has_disgas, has_vapoil, terminal_output),
          polymer_props_ad_(polymer_props_ad),
          has_polymer_(has_polymer),
          has_plyshlog_(has_plyshlog),
          has_shrate_(has_shrate),
          poly_pos_(detail::polymerPos(fluid.phaseUsage())),
          wells_rep_radius_(wells_rep_radius),
          wells_perf_length_(wells_perf_length),
          wells_bore_diameter_(wells_bore_diameter)
    {
        if (has_polymer_) {
            if (!active_[Water]) {
                OPM_THROW(std::logic_error, "Polymer must solved in water!\n");
            }
            // If deck has polymer, residual_ should contain polymer equation.
            rq_.resize(fluid_.numPhases() + 1);
            residual_.material_balance_eq.resize(fluid_.numPhases() + 1, ADB::null());
            Base::material_name_.push_back("Polymer");
            assert(poly_pos_ == fluid_.numPhases());
        }
    }



    template <class Grid>
    void
    BlackoilPolymerModel<Grid>::
    prepareStep(const double dt,
                ReservoirState& reservoir_state,
                WellState& well_state)
    {
        Base::prepareStep(dt, reservoir_state, well_state);
        auto& max_concentration = reservoir_state.getCellData( reservoir_state.CMAX );
        // Initial max concentration of this time step from PolymerBlackoilState.

        cmax_ = Eigen::Map<const V>(max_concentration.data(), Opm::AutoDiffGrid::numCells(grid_));
    }




    template <class Grid>
    void
    BlackoilPolymerModel<Grid>::
    afterStep(const double /* dt */,
              ReservoirState& reservoir_state,
              WellState& /* well_state */)
    {
        computeCmax(reservoir_state);
    }





    template <class Grid>
    void
    BlackoilPolymerModel<Grid>::makeConstantState(SolutionState& state) const
    {
        Base::makeConstantState(state);
        state.concentration = ADB::constant(state.concentration.value());
    }





    template <class Grid>
    std::vector<V>
    BlackoilPolymerModel<Grid>::variableStateInitials(const ReservoirState& x,
                                                      const WellState& xw) const
    {
        std::vector<V> vars0 = Base::variableStateInitials(x, xw);
        assert(int(vars0.size()) == fluid_.numPhases() + 2);

        // Initial polymer concentration.
        if (has_polymer_) {
            const auto& concentration = x.getCellData( x.CONCENTRATION );
            assert (not concentration.empty());
            const int nc = concentration.size();
            const V c = Eigen::Map<const V>(concentration.data() , nc);
            // Concentration belongs after other reservoir vars but before well vars.
            auto concentration_pos = vars0.begin() + fluid_.numPhases();
            assert(concentration_pos == vars0.end() - 2);
            vars0.insert(concentration_pos, c);
        }
        return vars0;
    }





    template <class Grid>
    std::vector<int>
    BlackoilPolymerModel<Grid>::variableStateIndices() const
    {
        std::vector<int> ind = Base::variableStateIndices();
        assert(ind.size() == 5);
        if (has_polymer_) {
            ind.resize(6);
            // Concentration belongs after other reservoir vars but before well vars.
            ind[Concentration] = fluid_.numPhases();
            // Concentration is pushing back the well vars.
            ++ind[Qs];
            ++ind[Bhp];
        }
        return ind;
    }




    template <class Grid>
    typename BlackoilPolymerModel<Grid>::SolutionState
    BlackoilPolymerModel<Grid>::variableStateExtractVars(const ReservoirState& x,
                                                         const std::vector<int>& indices,
                                                         std::vector<ADB>& vars) const
    {
        SolutionState state = Base::variableStateExtractVars(x, indices, vars);
        if (has_polymer_) {
            state.concentration = std::move(vars[indices[Concentration]]);
        }
        return state;
    }





    template <class Grid>
    void
    BlackoilPolymerModel<Grid>::computeAccum(const SolutionState& state,
                                             const int            aix  )
    {
        Base::computeAccum(state, aix);

        // Compute accumulation of polymer equation only if needed.
        if (has_polymer_) {
            const ADB& press = state.pressure;
            const std::vector<ADB>& sat = state.saturation;
            const ADB& c = state.concentration;
            const ADB pv_mult = poroMult(press); // also computed in Base::computeAccum, could be optimized.
            const Opm::PhaseUsage& pu = fluid_.phaseUsage();
            // compute polymer properties.
            const ADB cmax = ADB::constant(cmax_, state.concentration.blockPattern());
            const ADB ads  = polymer_props_ad_.adsorption(state.concentration, cmax);
            const double rho_rock = polymer_props_ad_.rockDensity();
            const V phi = Eigen::Map<const V>(&fluid_.porosity()[0], AutoDiffGrid::numCells(grid_));
            const double dead_pore_vol = polymer_props_ad_.deadPoreVol();
            // Compute polymer accumulation term.
            rq_[poly_pos_].accum[aix] = pv_mult * rq_[pu.phase_pos[Water]].b * sat[pu.phase_pos[Water]] * c * (1. - dead_pore_vol) 
                                        + pv_mult * rho_rock * (1. - phi) / phi * ads;
        }
 
    }






    template <class Grid>
    void BlackoilPolymerModel<Grid>::computeCmax(ReservoirState& state)
    {
        auto& max_concentration = state.getCellData( state.CMAX );
        const auto& concentration = state.getCellData( state.CONCENTRATION );
        std::transform( max_concentration.begin() ,
                        max_concentration.end() ,
                        concentration.begin() ,
                        max_concentration.begin() ,
                        [](double c_max , double c) { return std::max( c_max , c ); });

    }





    template <class Grid>
    void
    BlackoilPolymerModel<Grid>::
    assembleMassBalanceEq(const SolutionState& state)
    {
        // Base::assembleMassBalanceEq(state);

        // Compute b_p and the accumulation term b_p*s_p for each phase,
        // except gas. For gas, we compute b_g*s_g + Rs*b_o*s_o.
        // These quantities are stored in rq_[phase].accum[1].
        // The corresponding accumulation terms from the start of
        // the timestep (b^0_p*s^0_p etc.) were already computed
        // on the initial call to assemble() and stored in rq_[phase].accum[0].
        computeAccum(state, 1);

        // Set up the common parts of the mass balance equations
        // for each active phase.
        const V transi = subset(geo_.transmissibility(), ops_.internal_faces);
        const std::vector<ADB> kr = computeRelPerm(state);


        if (has_plyshlog_) {
            std::vector<double> water_vel;
            std::vector<double> visc_mult;

            computeWaterShearVelocityFaces(transi, kr, state.canonical_phase_pressures, state, water_vel, visc_mult);
            if ( !polymer_props_ad_.computeShearMultLog(water_vel, visc_mult, shear_mult_faces_) ) {
                // std::cerr << " failed in calculating the shear-multiplier " << std::endl;
                OPM_THROW(std::runtime_error, " failed in calculating the shear-multiplier. ");
            }
        }

        for (int phaseIdx = 0; phaseIdx < fluid_.numPhases(); ++phaseIdx) {
            const std::vector<PhasePresence>& cond = phaseCondition();
            const ADB mu = fluidViscosity(canph_[phaseIdx], state.canonical_phase_pressures[canph_[phaseIdx]], state.temperature, state.rs, state.rv, cond);
            const ADB rho = fluidDensity(canph_[phaseIdx], rq_[phaseIdx].b, state.rs, state.rv);
            computeMassFlux(phaseIdx, transi, kr[canph_[phaseIdx]], mu, rho, state.canonical_phase_pressures[canph_[phaseIdx]], state);

            residual_.material_balance_eq[ phaseIdx ] =
                pvdt_ * (rq_[phaseIdx].accum[1] - rq_[phaseIdx].accum[0])
                + ops_.div*rq_[phaseIdx].mflux;
        }

        // -------- Extra (optional) rs and rv contributions to the mass balance equations --------

        // Add the extra (flux) terms to the mass balance equations
        // From gas dissolved in the oil phase (rs) and oil vaporized in the gas phase (rv)
        // The extra terms in the accumulation part of the equation are already handled.
        if (active_[ Oil ] && active_[ Gas ]) {
            const int po = fluid_.phaseUsage().phase_pos[ Oil ];
            const int pg = fluid_.phaseUsage().phase_pos[ Gas ];

            const UpwindSelector<double> upwindOil(grid_, ops_,
                                                rq_[po].dh.value());
            const ADB rs_face = upwindOil.select(state.rs);

            const UpwindSelector<double> upwindGas(grid_, ops_,
                                                rq_[pg].dh.value());
            const ADB rv_face = upwindGas.select(state.rv);

            residual_.material_balance_eq[ pg ] += ops_.div * (rs_face * rq_[po].mflux);
            residual_.material_balance_eq[ po ] += ops_.div * (rv_face * rq_[pg].mflux);

            // OPM_AD_DUMP(residual_.material_balance_eq[ Gas ]);

        }

        // Add polymer equation.
        if (has_polymer_) {
            residual_.material_balance_eq[poly_pos_] = pvdt_ * (rq_[poly_pos_].accum[1] - rq_[poly_pos_].accum[0])
                                               + ops_.div*rq_[poly_pos_].mflux;
        }
    }





    template <class Grid>
    void BlackoilPolymerModel<Grid>::addWellContributionToMassBalanceEq(const std::vector<ADB>& cq_s,
                                                                        const SolutionState& state,
                                                                        WellState& xw)

    {
        Base::addWellContributionToMassBalanceEq(cq_s, state, xw);

        // Add well contributions to polymer mass balance equation
        if (has_polymer_) {
            const ADB mc = computeMc(state);
            const int nc = xw.polymerInflow().size();
            const V polyin = Eigen::Map<const V>(xw.polymerInflow().data(), nc);
            const int nperf = wells().well_connpos[wells().number_of_wells];
            const std::vector<int> well_cells(wells().well_cells, wells().well_cells + nperf);
            const V poly_in_perf = subset(polyin, well_cells);
            const V poly_mc_perf = subset(mc.value(), well_cells);
            const ADB& cq_s_water = cq_s[fluid_.phaseUsage().phase_pos[Water]];
            Selector<double> injector_selector(cq_s_water.value());
            const V poly_perf = injector_selector.select(poly_in_perf, poly_mc_perf);
            const ADB cq_s_poly =  cq_s_water * poly_perf;
            residual_.material_balance_eq[poly_pos_] -= superset(cq_s_poly, well_cells, nc);
        }
    }






    template <class Grid>
    void BlackoilPolymerModel<Grid>::updateState(const V& dx,
                                                 ReservoirState& reservoir_state,
                                                 WellState& well_state)
    {
        if (has_polymer_) {
            // Extract concentration change.
            const int np = fluid_.numPhases();
            const int nc = Opm::AutoDiffGrid::numCells(grid_);
            const V zero = V::Zero(nc);
            const int concentration_start = nc * np;
            const V dc = subset(dx, Span(nc, 1, concentration_start));

            // Create new dx with the dc part deleted.
            V modified_dx = V::Zero(dx.size() - nc);
            modified_dx.head(concentration_start) = dx.head(concentration_start);
            const int tail_len = dx.size() - concentration_start - nc;
            modified_dx.tail(tail_len) = dx.tail(tail_len);

            // Call base version.
            Base::updateState(modified_dx, reservoir_state, well_state);

            {
                auto& concentration = reservoir_state.getCellData( reservoir_state.CONCENTRATION );
                // Update concentration.
                const V c_old = Eigen::Map<const V>(concentration.data(), nc, 1);
                const V c = (c_old - dc).max(zero);
                std::copy(&c[0], &c[0] + nc, concentration.begin());
            }
        } else {
            // Just forward call to base version.
            Base::updateState(dx, reservoir_state, well_state);
        }
    }





    template <class Grid>
    void
    BlackoilPolymerModel<Grid>::computeMassFlux(const int               actph ,
                                                const V&                transi,
                                                const ADB&              kr    ,
                                                const ADB&              mu    ,
                                                const ADB&              rho   ,
                                                const ADB&              phasePressure,
                                                const SolutionState&    state)
    {
        Base::computeMassFlux(actph, transi, kr, mu, rho, phasePressure, state);

        // Polymer treatment.
        const int canonicalPhaseIdx = canph_[ actph ];
        if (canonicalPhaseIdx == Water) {
            if (has_polymer_) {
                const ADB tr_mult = transMult(state.pressure);
                const ADB cmax = ADB::constant(cmax_, state.concentration.blockPattern());
                const ADB mc = computeMc(state);
                const ADB krw_eff = polymer_props_ad_.effectiveRelPerm(state.concentration, cmax, kr);
                const ADB inv_wat_eff_visc = polymer_props_ad_.effectiveInvWaterVisc(state.concentration, mu.value().data());
                // Reduce mobility of water phase by relperm reduction and effective viscosity increase.
                rq_[actph].mob = tr_mult * krw_eff * inv_wat_eff_visc;
                // Compute polymer mobility.
                rq_[poly_pos_].mob = tr_mult * mc * krw_eff * inv_wat_eff_visc;
                rq_[poly_pos_].b = rq_[actph].b;
                rq_[poly_pos_].dh = rq_[actph].dh;
                UpwindSelector<double> upwind(grid_, ops_, rq_[poly_pos_].dh.value());
                // Compute polymer flux.
                rq_[poly_pos_].mflux = upwind.select(rq_[poly_pos_].b * rq_[poly_pos_].mob) * (transi * rq_[poly_pos_].dh);
                // Must recompute water flux since we have to use modified mobilities.
                rq_[ actph ].mflux = upwind.select(rq_[actph].b * rq_[actph].mob) * (transi * rq_[actph].dh);

                // applying the shear-thinning factors
                if (has_plyshlog_) {
                    V shear_mult_faces_v = Eigen::Map<V>(shear_mult_faces_.data(), shear_mult_faces_.size());
                    ADB shear_mult_faces_adb = ADB::constant(shear_mult_faces_v);
                    rq_[poly_pos_].mflux = rq_[poly_pos_].mflux / shear_mult_faces_adb;
                    rq_[actph].mflux = rq_[actph].mflux / shear_mult_faces_adb;
                }
            }
        }
    }




    template <class Grid>
    void
    BlackoilPolymerModel<Grid>::assemble(const ReservoirState& reservoir_state,
                                         WellState& well_state,
                                         const bool initial_assembly)
    {
        using namespace Opm::AutoDiffGrid;

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
            computeAccum(state0, 0);
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
        assembleMassBalanceEq(state);

        // -------- Well equations ----------
        if ( ! wellsActive() ) {
            return;
        }

        V aliveWells;

        const int np = wells().number_of_phases;
        std::vector<ADB> cq_s(np, ADB::null());

        const int nw = wells().number_of_wells;
        const int nperf = wells().well_connpos[nw];
        const std::vector<int> well_cells(wells().well_cells, wells().well_cells + nperf);

        std::vector<ADB> mob_perfcells(np, ADB::null());
        std::vector<ADB> b_perfcells(np, ADB::null());
        for (int phase = 0; phase < np; ++phase) {
            mob_perfcells[phase] = subset(rq_[phase].mob, well_cells);
            b_perfcells[phase] = subset(rq_[phase].b, well_cells);
        }
        if (param_.solve_welleq_initially_ && initial_assembly) {
            // solve the well equations as a pre-processing step
            Base::solveWellEq(mob_perfcells, b_perfcells, state, well_state);
        }

        Base::computeWellFlux(state, mob_perfcells, b_perfcells, aliveWells, cq_s);

        if (has_plyshlog_) {
            std::vector<double> water_vel_wells;
            std::vector<double> visc_mult_wells;

            const int water_pos = fluid_.phaseUsage().phase_pos[Water];
            computeWaterShearVelocityWells(state, well_state, cq_s[water_pos], water_vel_wells, visc_mult_wells);

            if ( !polymer_props_ad_.computeShearMultLog(water_vel_wells, visc_mult_wells, shear_mult_wells_) ) {
                OPM_THROW(std::runtime_error, " failed in calculating the shear factors for wells ");
            }

            // applying the shear-thinning to the water phase
            V shear_mult_wells_v = Eigen::Map<V>(shear_mult_wells_.data(), shear_mult_wells_.size());
            ADB shear_mult_wells_adb = ADB::constant(shear_mult_wells_v);
            mob_perfcells[water_pos] = mob_perfcells[water_pos] / shear_mult_wells_adb;
        }

        Base::computeWellFlux(state, mob_perfcells, b_perfcells, aliveWells, cq_s);
        Base::updatePerfPhaseRatesAndPressures(cq_s, state, well_state);
        Base::addWellFluxEq(cq_s, state);
        addWellContributionToMassBalanceEq(cq_s, state, well_state);
        addWellControlEq(state, well_state, aliveWells);
    }




    template <class Grid>
    ADB
    BlackoilPolymerModel<Grid>::computeMc(const SolutionState& state) const
    {
        return polymer_props_ad_.polymerWaterVelocityRatio(state.concentration);
    }




    template<class Grid>
    void
    BlackoilPolymerModel<Grid>::computeWaterShearVelocityFaces(const V& transi, const std::vector<ADB>& kr,
                                                               const std::vector<ADB>& phasePressure, const SolutionState& state,
                                                               std::vector<double>& water_vel, std::vector<double>& visc_mult)
    {

        const int phase = fluid_.phaseUsage().phase_pos[Water]; // water position

        const int canonicalPhaseIdx = canph_[phase];

        const std::vector<PhasePresence> cond = phaseCondition();

        const ADB tr_mult = transMult(state.pressure);
        const ADB mu    = fluidViscosity(canonicalPhaseIdx, phasePressure[canonicalPhaseIdx], state.temperature, state.rs, state.rv, cond);
        rq_[phase].mob = tr_mult * kr[canonicalPhaseIdx] / mu;

        // compute gravity potensial using the face average as in eclipse and MRST
        const ADB rho   = fluidDensity(canonicalPhaseIdx, rq_[phase].b, state.rs, state.rv);
        const ADB rhoavg = ops_.caver * rho;
        rq_[ phase ].dh = ops_.ngrad * phasePressure[ canonicalPhaseIdx ] - geo_.gravity()[2] * (rhoavg * (ops_.ngrad * geo_.z().matrix()));
        if (use_threshold_pressure_) {
            applyThresholdPressures(rq_[ phase ].dh);
        }

        const ADB& b   = rq_[ phase ].b;
        const ADB& mob = rq_[ phase ].mob;
        const ADB& dh  = rq_[ phase ].dh;
        UpwindSelector<double> upwind(grid_, ops_, dh.value());

        const ADB cmax = ADB::constant(cmax_, state.concentration.blockPattern());
        const ADB mc = computeMc(state);
        ADB krw_eff = polymer_props_ad_.effectiveRelPerm(state.concentration,
                                                         cmax,
                                                         kr[canonicalPhaseIdx]);
        ADB inv_wat_eff_visc = polymer_props_ad_.effectiveInvWaterVisc(state.concentration, mu.value().data());
        rq_[ phase ].mob = tr_mult * krw_eff * inv_wat_eff_visc;

        const V& polymer_conc = state.concentration.value();
        V visc_mult_cells = polymer_props_ad_.viscMult(polymer_conc);
        V visc_mult_faces = upwind.select(visc_mult_cells);

        size_t nface = visc_mult_faces.size();
        visc_mult.resize(nface);
        std::copy(visc_mult_faces.data(), visc_mult_faces.data() + nface, visc_mult.begin());

        rq_[ phase ].mflux = (transi * upwind.select(b * mob)) * dh;


        const auto& b_faces_adb = upwind.select(b);
        std::vector<double> b_faces(b_faces_adb.value().data(), b_faces_adb.value().data() + b_faces_adb.size());

        const auto& internal_faces = ops_.internal_faces;

        std::vector<double> internal_face_areas;
        internal_face_areas.resize(internal_faces.size());

        for (int i = 0; i < internal_faces.size(); ++i) {
            internal_face_areas[i] = grid_.face_areas[internal_faces[i]];
        }

        const ADB phi = Opm::AutoDiffBlock<double>::constant(Eigen::Map<const V>(& fluid_.porosity()[0], AutoDiffGrid::numCells(grid_), 1));
        const ADB phiavg_adb  = ops_.caver * phi;

        std::vector<double> phiavg(phiavg_adb.value().data(), phiavg_adb.value().data() + phiavg_adb.size());

        water_vel.resize(nface);
        std::copy(rq_[0].mflux.value().data(), rq_[0].mflux.value().data() + nface, water_vel.begin());

        for (size_t i = 0; i < nface; ++i) {
            water_vel[i] = water_vel[i] / (b_faces[i] * phiavg[i] * internal_face_areas[i]);
        }

        // for SHRATE keyword treatment
        if (has_shrate_) {

            // get the upwind water saturation
            const Opm::PhaseUsage pu = fluid_.phaseUsage();
            const ADB& sw = state.saturation[pu.phase_pos[ Water ]];
            const ADB& sw_upwind_adb = upwind.select(sw);
            std::vector<double> sw_upwind(sw_upwind_adb.value().data(), sw_upwind_adb.value().data() + sw_upwind_adb.size());

            // get the absolute permeability for the faces
            std::vector<double> perm;
            perm.resize(transi.size());

            for (int i = 0; i < transi.size(); ++i) {
                perm[i] = transi[i] / internal_faces[i];
            }

            // get the upwind krw_eff
            const ADB& krw_adb = upwind.select(krw_eff);
            std::vector<double> krw_upwind(krw_adb.value().data(), krw_adb.value().data() + krw_adb.size());

            const double& shrate_const = polymer_props_ad_.shrate();

            const double epsilon = std::numeric_limits<double>::epsilon();
            // std::cout << "espilon is " << epsilon << std::endl;
            // std::cin.ignore();

            for (size_t i = 0; i < water_vel.size(); ++i) {
                // assuming only when upwinding water saturation is not zero
                // there will be non-zero water velocity
                if (std::abs(water_vel[i]) < epsilon) {
                    continue;
                }

                water_vel[i] *= shrate_const * std::sqrt(phiavg[i] / (perm[i] * sw_upwind[i] * krw_upwind[i]));

            }
        }

    }




    template<class Grid>
    void
    BlackoilPolymerModel<Grid>::computeWaterShearVelocityWells(const SolutionState& state, WellState& xw, const ADB& cq_sw,
                                                               std::vector<double>& water_vel_wells, std::vector<double>& visc_mult_wells)
    {
        if( ! wellsActive() ) return ;

        const int nw = wells().number_of_wells;
        const int nperf = wells().well_connpos[nw];
        const std::vector<int> well_cells(wells().well_cells, wells().well_cells + nperf);

        water_vel_wells.resize(cq_sw.size());
        std::copy(cq_sw.value().data(), cq_sw.value().data() + cq_sw.size(), water_vel_wells.begin());

        const V& polymer_conc = state.concentration.value();

        V visc_mult_cells = polymer_props_ad_.viscMult(polymer_conc);
        V visc_mult_wells_v = subset(visc_mult_cells, well_cells);

        visc_mult_wells.resize(visc_mult_wells_v.size());
        std::copy(visc_mult_wells_v.data(), visc_mult_wells_v.data() + visc_mult_wells_v.size(), visc_mult_wells.begin());

        const int water_pos = fluid_.phaseUsage().phase_pos[Water];
        ADB b_perfcells = subset(rq_[water_pos].b, well_cells);

        const ADB& p_perfcells = subset(state.pressure, well_cells);
        const V& cdp = well_perforation_pressure_diffs_;
        const ADB perfpressure = (wops_.w2p * state.bhp) + cdp;
        // Pressure drawdown (also used to determine direction of flow)
        const ADB drawdown =  p_perfcells - perfpressure;

        // selects injection perforations
        V selectInjectingPerforations = V::Zero(nperf);
        for (int c = 0; c < nperf; ++c) {
            if (drawdown.value()[c] < 0) {
                selectInjectingPerforations[c] = 1;
            }
        }

        // for the injection wells
        for (size_t i = 0; i < well_cells.size(); ++i) {
            if (xw.polymerInflow()[well_cells[i]] == 0. && selectInjectingPerforations[i] == 1) { // maybe comparison with epsilon threshold
                visc_mult_wells[i] = 1.;
            }
        }

        const ADB phi = Opm::AutoDiffBlock<double>::constant(Eigen::Map<const V>(& fluid_.porosity()[0], AutoDiffGrid::numCells(grid_), 1));
        const ADB phi_wells_adb = subset(phi, well_cells);

        std::vector<double> phi_wells(phi_wells_adb.value().data(), phi_wells_adb.value().data() + phi_wells_adb.size());

        std::vector<double> b_wells(b_perfcells.value().data(), b_perfcells.value().data() + b_perfcells.size());

        for (size_t i = 0; i < water_vel_wells.size(); ++i) {
            water_vel_wells[i] = b_wells[i] * water_vel_wells[i] / (phi_wells[i] * 2. * M_PI * wells_rep_radius_[i] * wells_perf_length_[i]);
            // TODO: CHECK to make sure this formulation is corectly used. Why muliplied by bW.
            // Although this formulation works perfectly with the tests compared with other formulations
        }

        // for SHRATE treatment
        if (has_shrate_) {
            const double& shrate_const = polymer_props_ad_.shrate();
            for (size_t i = 0; i < water_vel_wells.size(); ++i) {
                water_vel_wells[i] = shrate_const * water_vel_wells[i] / wells_bore_diameter_[i];
            }
        }

        return;
    }

} // namespace Opm

#endif // OPM_BLACKOILPOLYMERMODEL_IMPL_HEADER_INCLUDED
