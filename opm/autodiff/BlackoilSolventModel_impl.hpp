/*
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

#ifndef OPM_BLACKOILSOLVENTMODEL_IMPL_HEADER_INCLUDED
#define OPM_BLACKOILSOLVENTMODEL_IMPL_HEADER_INCLUDED

#include <opm/autodiff/BlackoilSolventModel.hpp>

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
        int solventPos(const PU& pu)
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
    BlackoilSolventModel<Grid>::BlackoilSolventModel(const typename Base::ModelParameters&   param,
                                                     const Grid&                             grid,
                                                     const BlackoilPropsAdInterface&         fluid,
                                                     const DerivedGeology&                   geo,
                                                     const RockCompressibility*              rock_comp_props,
                                                     const SolventPropsAdFromDeck&           solvent_props,
                                                     const Wells*                            wells_arg,
                                                     const NewtonIterationBlackoilInterface& linsolver,
                                                     const EclipseStateConstPtr              eclState,
                                                     const bool                              has_disgas,
                                                     const bool                              has_vapoil,
                                                     const bool                              terminal_output,
                                                     const bool                              has_solvent,
                                                     const bool                              is_miscible)
        : Base(param, grid, fluid, geo, rock_comp_props, wells_arg, linsolver,
               eclState, has_disgas, has_vapoil, terminal_output),
          has_solvent_(has_solvent),
          solvent_pos_(detail::solventPos(fluid.phaseUsage())),
          solvent_props_(solvent_props),
          is_miscible_(is_miscible)

    {
        if (has_solvent_) {

            // If deck has solvent, residual_ should contain solvent equation.
            rq_.resize(fluid_.numPhases() + 1);
            residual_.material_balance_eq.resize(fluid_.numPhases() + 1, ADB::null());
            Base::material_name_.push_back("Solvent");
            assert(solvent_pos_ == fluid_.numPhases());
            if (has_vapoil_) {
                OPM_THROW(std::runtime_error, "Solvent option only works with dead gas\n");
            }

            residual_.matbalscale.resize(fluid_.numPhases() + 1, 0.0031); // use the same as gas
        }
        if (is_miscible_) {
            mu_eff_.resize(fluid_.numPhases() + 1, ADB::null());
            b_eff_.resize(fluid_.numPhases() + 1, ADB::null());
        }
    }





    template <class Grid>
    void
    BlackoilSolventModel<Grid>::makeConstantState(SolutionState& state) const
    {
        Base::makeConstantState(state);
        state.solvent_saturation = ADB::constant(state.solvent_saturation.value());
    }





    template <class Grid>
    std::vector<V>
    BlackoilSolventModel<Grid>::variableStateInitials(const ReservoirState& x,
                                                      const WellState& xw) const
    {
        std::vector<V> vars0 = Base::variableStateInitials(x, xw);
        assert(int(vars0.size()) == fluid_.numPhases() + 2);

        // Initial polymer concentration.
        if (has_solvent_) {
            const auto& solvent_saturation = x.getCellData( BlackoilSolventState::SSOL );
            const int nc = solvent_saturation.size();
            const V ss = Eigen::Map<const V>(solvent_saturation.data() , nc);

            // This is some insanely detailed flickety flackety code;
            // Solvent belongs after other reservoir vars but before well vars.
            auto solvent_pos = vars0.begin() + fluid_.numPhases();
            assert (not solvent_saturation.empty());
            assert(solvent_pos == vars0.end() - 2);
            vars0.insert(solvent_pos, ss);
        }
        return vars0;
    }





    template <class Grid>
    std::vector<int>
    BlackoilSolventModel<Grid>::variableStateIndices() const
    {
        std::vector<int> ind = Base::variableStateIndices();
        assert(ind.size() == 5);
        if (has_solvent_) {
            ind.resize(6);
            // Solvent belongs after other reservoir vars but before well vars.
            ind[Solvent] = fluid_.numPhases();
            // Solvent is pushing back the well vars.
            ++ind[Qs];
            ++ind[Bhp];
        }
        return ind;
    }




    template <class Grid>
    typename BlackoilSolventModel<Grid>::SolutionState
    BlackoilSolventModel<Grid>::variableStateExtractVars(const ReservoirState& x,
                                                         const std::vector<int>& indices,
                                                         std::vector<ADB>& vars) const
    {         
        // This is more or less a copy of the base class. Refactoring is needed in the base class
        // to avoid unnecessary copying.

        //using namespace Opm::AutoDiffGrid;
        const int nc = Opm::AutoDiffGrid::numCells(grid_);
        const Opm::PhaseUsage pu = fluid_.phaseUsage();

        SolutionState state(fluid_.numPhases());

        // Pressure.
        state.pressure = std::move(vars[indices[Pressure]]);

        // Temperature cannot be a variable at this time (only constant).
        const V temp = Eigen::Map<const V>(& x.temperature()[0], x.temperature().size());
        state.temperature = ADB::constant(temp);

        // Saturations
        {
            ADB so = ADB::constant(V::Ones(nc, 1));

            if (active_[ Water ]) {
                state.saturation[pu.phase_pos[ Water ]] = std::move(vars[indices[Sw]]);
                const ADB& sw = state.saturation[pu.phase_pos[ Water ]];
                so -= sw;
            }
            if (has_solvent_) {
                state.solvent_saturation = std::move(vars[indices[Solvent]]);
                so -= state.solvent_saturation;
            }

            if (active_[ Gas ]) {
                // Define Sg Rs and Rv in terms of xvar.
                // Xvar is only defined if gas phase is active
                const ADB& xvar = vars[indices[Xvar]];
                ADB& sg = state.saturation[ pu.phase_pos[ Gas ] ];
                sg = Base::isSg_*xvar + Base::isRv_*so;
                so -= sg;

                if (active_[ Oil ]) {
                    // RS and RV is only defined if both oil and gas phase are active.
                    state.canonical_phase_pressures = computePressures(state.pressure, state.saturation[pu.phase_pos[ Water ]], so, sg, state.solvent_saturation);
                    const ADB rsSat = fluidRsSat(state.canonical_phase_pressures[ Oil ], so , cells_);
                    if (has_disgas_) {
                        state.rs = (1-Base::isRs_)*rsSat + Base::isRs_*xvar;
                    } else {
                        state.rs = rsSat;
                    }
                    const ADB rvSat = fluidRvSat(state.canonical_phase_pressures[ Gas ], so , cells_);
                    if (has_vapoil_) {
                        state.rv = (1-Base::isRv_)*rvSat + Base::isRv_*xvar;
                    } else {
                        state.rv = rvSat;
                    }
                }
            }

            if (active_[ Oil ]) {
                // Note that so is never a primary variable.
                state.saturation[pu.phase_pos[ Oil ]] = std::move(so);
            }
        }
        // wells
        Base::variableStateExtractWellsVars(indices, vars, state);
        return state;

    }





    template <class Grid>
    void
    BlackoilSolventModel<Grid>::computeAccum(const SolutionState& state,
                                             const int            aix  )
    {
        Base::computeAccum(state, aix);

        // Compute accumulation of the solvent
        if (has_solvent_) {
            const ADB& press = state.pressure;
            const ADB& ss = state.solvent_saturation;
            const ADB pv_mult = poroMult(press); // also computed in Base::computeAccum, could be optimized.
            const Opm::PhaseUsage& pu = fluid_.phaseUsage();

            const ADB& pg = state.canonical_phase_pressures[pu.phase_pos[Gas]];
            const std::vector<PhasePresence>& cond = phaseCondition();
            rq_[solvent_pos_].b = fluidReciprocFVF(Solvent, pg, state.temperature, state.rs, state.rv,cond);
            rq_[solvent_pos_].accum[aix] = pv_mult * rq_[solvent_pos_].b * ss;
        }
    }





    template <class Grid>
    void
    BlackoilSolventModel<Grid>::
    assembleMassBalanceEq(const SolutionState& state)
    {

        Base::assembleMassBalanceEq(state);

        if (has_solvent_) {
            residual_.material_balance_eq[ solvent_pos_ ] =
                pvdt_ * (rq_[solvent_pos_].accum[1] - rq_[solvent_pos_].accum[0])
                + ops_.div*rq_[solvent_pos_].mflux;
        }

    }

    template <class Grid>
    void
    BlackoilSolventModel<Grid>::updateEquationsScaling()
    {
        Base::updateEquationsScaling();
        assert(MaxNumPhases + 1 == residual_.matbalscale.size());
        if (has_solvent_) {
            const ADB& temp_b = rq_[solvent_pos_].b;
            ADB::V B = 1. / temp_b.value();
#if HAVE_MPI
            if ( linsolver_.parallelInformation().type() == typeid(ParallelISTLInformation) )
            {
                const ParallelISTLInformation& real_info =
                    boost::any_cast<const ParallelISTLInformation&>(linsolver_.parallelInformation());
                double B_global_sum = 0;
                real_info.computeReduction(B, Reduction::makeGlobalSumFunctor<double>(), B_global_sum);
                residual_.matbalscale[solvent_pos_] = B_global_sum / Base::global_nc_;
            }
            else
#endif
            {
                residual_.matbalscale[solvent_pos_] = B.mean();
            }
        }
    }

    template <class Grid>
    void BlackoilSolventModel<Grid>::addWellContributionToMassBalanceEq(const std::vector<ADB>& cq_s,
                                                                        const SolutionState& state,
                                                                        WellState& xw)

    {

        // Add well contributions to solvent mass balance equation

        Base::addWellContributionToMassBalanceEq(cq_s, state, xw);

        if (has_solvent_) {
            const int nperf = wells().well_connpos[wells().number_of_wells];
            const int nc = Opm::AutoDiffGrid::numCells(grid_);

            const Opm::PhaseUsage& pu = fluid_.phaseUsage();
            const ADB zero = ADB::constant(V::Zero(nc));
            const ADB& ss = state.solvent_saturation;
            const ADB& sg = (active_[ Gas ]
                             ? state.saturation[ pu.phase_pos[ Gas ] ]
                             : zero);

            const std::vector<int> well_cells(wells().well_cells, wells().well_cells + nperf);
            Selector<double> zero_selector(ss.value() + sg.value(), Selector<double>::Zero);
            ADB F_solvent = subset(zero_selector.select(ss, ss / (ss + sg)),well_cells);

            const int nw = wells().number_of_wells;
            V injectedSolventFraction = Eigen::Map<const V>(&xw.solventFraction()[0], nperf);

            V isProducer = V::Zero(nperf);
            V ones = V::Constant(nperf,1.0);
            for (int w = 0; w < nw; ++w) {
                if(wells().type[w] == PRODUCER) {
                    for (int perf = wells().well_connpos[w]; perf < wells().well_connpos[w+1]; ++perf) {
                        isProducer[perf] = 1;
                    }
                }
            }

            const ADB& rs_perfcells = subset(state.rs, well_cells);
            int gas_pos = fluid_.phaseUsage().phase_pos[Gas];
            int oil_pos = fluid_.phaseUsage().phase_pos[Oil];
            // remove contribution from the dissolved gas.
            // TODO compensate for gas in the oil phase
            assert(!has_vapoil_);
            const ADB cq_s_solvent = (isProducer * F_solvent + (ones - isProducer) * injectedSolventFraction) * (cq_s[gas_pos] - rs_perfcells * cq_s[oil_pos]);

            // Solvent contribution to the mass balance equation is given as a fraction
            // of the gas contribution.
            residual_.material_balance_eq[solvent_pos_] -= superset(cq_s_solvent, well_cells, nc);

            // The gas contribution must be reduced accordingly for the total contribution to be
            // the same.
            residual_.material_balance_eq[gas_pos] += superset(cq_s_solvent, well_cells, nc);

        }
    }

    template <class Grid>
    void BlackoilSolventModel<Grid>::computeWellConnectionPressures(const SolutionState& state,
                                                                        const WellState& xw)
    {
        if( ! Base::localWellsActive() ) return ;

        using namespace Opm::AutoDiffGrid;
        // 1. Compute properties required by computeConnectionPressureDelta().
        //    Note that some of the complexity of this part is due to the function
        //    taking std::vector<double> arguments, and not Eigen objects.
        const int nperf = wells().well_connpos[wells().number_of_wells];
        const int nw = wells().number_of_wells;
        const std::vector<int> well_cells(wells().well_cells, wells().well_cells + nperf);

        // Compute the average pressure in each well block
        const V perf_press = Eigen::Map<const V>(xw.perfPress().data(), nperf);
        V avg_press = perf_press*0;
        for (int w = 0; w < nw; ++w) {
            for (int perf = wells().well_connpos[w]; perf < wells().well_connpos[w+1]; ++perf) {
                const double p_above = perf == wells().well_connpos[w] ? state.bhp.value()[w] : perf_press[perf - 1];
                const double p_avg = (perf_press[perf] + p_above)/2;
                avg_press[perf] = p_avg;
            }
        }

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
        const V bw = fluid_.bWat(avg_press_ad, perf_temp, well_cells).value();
        if (pu.phase_used[BlackoilPhases::Aqua]) {
            b.col(pu.phase_pos[BlackoilPhases::Aqua]) = bw;
        }

        assert(active_[Oil]);
        assert(active_[Gas]);
        const ADB perf_rv = subset(state.rv, well_cells);
        const ADB perf_rs = subset(state.rs, well_cells);
        const V perf_so =  subset(state.saturation[pu.phase_pos[Oil]].value(), well_cells);
        if (pu.phase_used[BlackoilPhases::Liquid]) {
            const V bo = fluid_.bOil(avg_press_ad, perf_temp, perf_rs, perf_cond, well_cells).value();
            //const V bo_eff = subset(rq_[pu.phase_pos[Oil] ].b , well_cells).value();
            b.col(pu.phase_pos[BlackoilPhases::Liquid]) = bo;
            const V rssat = fluidRsSat(avg_press, perf_so, well_cells);
            rsmax_perf.assign(rssat.data(), rssat.data() + nperf);
        }
        V surf_dens_copy = superset(fluid_.surfaceDensity(0, well_cells), Span(nperf, pu.num_phases, 0), nperf*pu.num_phases);
        for (int phase = 1; phase < pu.num_phases; ++phase) {
            if ( phase == pu.phase_pos[BlackoilPhases::Vapour]) {
                continue; // the gas surface density is added after the solvent is accounted for.
            }
            surf_dens_copy += superset(fluid_.surfaceDensity(phase, well_cells), Span(nperf, pu.num_phases, phase), nperf*pu.num_phases);
        }

        if (pu.phase_used[BlackoilPhases::Vapour]) {
            // Unclear wether the effective or the pure values should be used for the wells
            // the current usage of unmodified properties values gives best match.
            //V bg_eff = subset(rq_[pu.phase_pos[Gas]].b,well_cells).value();
            V bg = fluid_.bGas(avg_press_ad, perf_temp, perf_rv, perf_cond, well_cells).value();
            V rhog = fluid_.surfaceDensity(pu.phase_pos[BlackoilPhases::Vapour], well_cells);
            if (has_solvent_) {

                const V bs = solvent_props_.bSolvent(avg_press_ad,well_cells).value();
                //const V bs_eff = subset(rq_[solvent_pos_].b,well_cells).value();

                // A weighted sum of the b-factors of gas and solvent are used.
                const int nc = Opm::AutoDiffGrid::numCells(grid_);

                const ADB zero = ADB::constant(V::Zero(nc));
                const ADB& ss = state.solvent_saturation;
                const ADB& sg = (active_[ Gas ]
                                 ? state.saturation[ pu.phase_pos[ Gas ] ]
                                 : zero);

                Selector<double> zero_selector(ss.value() + sg.value(), Selector<double>::Zero);
                V F_solvent = subset(zero_selector.select(ss, ss / (ss + sg)),well_cells).value();

                V injectedSolventFraction = Eigen::Map<const V>(&xw.solventFraction()[0], nperf);

                V isProducer = V::Zero(nperf);
                V ones = V::Constant(nperf,1.0);
                for (int w = 0; w < nw; ++w) {
                    if(wells().type[w] == PRODUCER) {
                        for (int perf = wells().well_connpos[w]; perf < wells().well_connpos[w+1]; ++perf) {
                            isProducer[perf] = 1;
                        }
                    }
                }

                F_solvent = isProducer * F_solvent + (ones - isProducer) * injectedSolventFraction;

                bg = bg * (ones - F_solvent);
                bg = bg + F_solvent * bs;

                const V& rhos = solvent_props_.solventSurfaceDensity(well_cells);
                rhog = ( (ones - F_solvent) * rhog ) + (F_solvent * rhos);
            }
            b.col(pu.phase_pos[BlackoilPhases::Vapour]) = bg;
            surf_dens_copy += superset(rhog, Span(nperf, pu.num_phases, pu.phase_pos[BlackoilPhases::Vapour]), nperf*pu.num_phases);

            const V rvsat = fluidRvSat(avg_press, perf_so, well_cells);
            rvmax_perf.assign(rvsat.data(), rvsat.data() + nperf);
        }

        // b and surf_dens_perf is row major, so can just copy data.
        std::vector<double> b_perf(b.data(), b.data() + nperf * pu.num_phases);        
        std::vector<double> surf_dens_perf(surf_dens_copy.data(), surf_dens_copy.data() + nperf * pu.num_phases);

        // Extract well connection depths.
        const V depth = cellCentroidsZToEigen(grid_);
        const V pdepth = subset(depth, well_cells);
        std::vector<double> perf_depth(pdepth.data(), pdepth.data() + nperf);

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
        Base::well_perforation_densities_ = Eigen::Map<const V>(cd.data(), nperf);
        Base::well_perforation_pressure_diffs_ = Eigen::Map<const V>(cdp.data(), nperf);
    }






    template <class Grid>
    void BlackoilSolventModel<Grid>::updateState(const V& dx,
                                                  ReservoirState& reservoir_state,
                                                  WellState& well_state)
    {

        if (has_solvent_) {
            // Extract solvent change.
            const int np = fluid_.numPhases();
            const int nc = Opm::AutoDiffGrid::numCells(grid_);
            const V zero = V::Zero(nc);
            const int solvent_start = nc * np;
            const V dss = subset(dx, Span(nc, 1, solvent_start));

            // Create new dx with the dss part deleted.
            V modified_dx = V::Zero(dx.size() - nc);
            modified_dx.head(solvent_start) = dx.head(solvent_start);
            const int tail_len = dx.size() - solvent_start - nc;
            modified_dx.tail(tail_len) = dx.tail(tail_len);

            // Call base version.
            Base::updateState(modified_dx, reservoir_state, well_state);

            // Update solvent.
            auto& solvent_saturation = reservoir_state.getCellData( reservoir_state.SSOL );
            const V ss_old = Eigen::Map<const V>(&solvent_saturation[0], nc, 1);
            const V ss = (ss_old - dss).max(zero);
            std::copy(&ss[0], &ss[0] + nc, solvent_saturation.begin());

            // adjust oil saturation
            const Opm::PhaseUsage& pu = fluid_.phaseUsage();
            const int oilpos = pu.phase_pos[ Oil ];
            for (int c = 0; c < nc; ++c) {
                reservoir_state.saturation()[c*np + oilpos] = 1 - ss[c];
                if (pu.phase_used[ Gas ]) {
                    const int gaspos = pu.phase_pos[ Gas ];
                    reservoir_state.saturation()[c*np + oilpos] -= reservoir_state.saturation()[c*np + gaspos];
                }
                if (pu.phase_used[ Water ]) {
                    const int waterpos = pu.phase_pos[ Water ];
                    reservoir_state.saturation()[c*np + oilpos] -= reservoir_state.saturation()[c*np + waterpos];
                }
            }

        } else {
            // Just forward call to base version.
            Base::updateState(dx, reservoir_state, well_state);
        }
    }





    template <class Grid>
    void
    BlackoilSolventModel<Grid>::computeMassFlux(const int               actph ,
                                                const V&                transi,
                                                const ADB&              kr    ,
                                                const ADB&              mu    ,
                                                const ADB&              rho    ,
                                                const ADB&              phasePressure,
                                                const SolutionState&    state)
    {

        const int canonicalPhaseIdx = canph_[ actph ];
        // make a copy to make it possible to modify it
        ADB kr_mod = kr;
        if (canonicalPhaseIdx == Gas) {
            if (has_solvent_) {
                const int  nc   = Opm::UgGridHelpers::numCells(grid_);
                const Opm::PhaseUsage& pu = fluid_.phaseUsage();
                const ADB zero = ADB::constant(V::Zero(nc));
                const V ones = V::Constant(nc, 1.0);

                const ADB& ss = state.solvent_saturation;
                const ADB& sg = (active_[ Gas ]
                                 ? state.saturation[ pu.phase_pos[ Gas ] ]
                                 : zero);
                Selector<double> zero_selector(ss.value() + sg.value(), Selector<double>::Zero);
                const ADB F_solvent = zero_selector.select(zero, ss / (ss + sg));

                // Compute solvent properties
                const std::vector<PhasePresence>& cond = phaseCondition();
                ADB mu_s = fluidViscosity(Solvent, phasePressure,state.temperature, state.rs, state.rv, cond);
                ADB rho_s = fluidDensity(Solvent,rq_[solvent_pos_].b, state.rs, state.rv);

                // Compute solvent relperm and mass flux
                ADB krs = solvent_props_.solventRelPermMultiplier(F_solvent, cells_) * kr_mod;
                Base::computeMassFlux(solvent_pos_, transi, krs, mu_s, rho_s, phasePressure, state);

                // Modify gas relperm
                kr_mod = solvent_props_.gasRelPermMultiplier( (ones - F_solvent) , cells_) * kr_mod;

            }
        }
        // Compute mobility and flux
        Base::computeMassFlux(actph, transi, kr_mod, mu, rho, phasePressure, state);

    }

    template <class Grid>
    ADB
    BlackoilSolventModel<Grid>::fluidViscosity(const int               phase,
                                                            const ADB&              p    ,
                                                            const ADB&              temp ,
                                                            const ADB&              rs   ,
                                                            const ADB&              rv   ,
                                                            const std::vector<PhasePresence>& cond) const
    {
        const Opm::PhaseUsage& pu = fluid_.phaseUsage();
        if (phase == Solvent) {
            if (!is_miscible_) {
                return solvent_props_.muSolvent(p, cells_);
            } else {
                return mu_eff_[solvent_pos_];
            }

        } else {
            if (!is_miscible_) {
                return Base::fluidViscosity(phase, p, temp, rs, rv, cond);
            } else {
                return mu_eff_[pu.phase_pos[ phase ]];
            }
        }
     }

    template <class Grid>
    ADB
    BlackoilSolventModel<Grid>::fluidReciprocFVF(const int               phase,
                                                 const ADB&              p    ,
                                                 const ADB&              temp ,
                                                 const ADB&              rs   ,
                                                 const ADB&              rv   ,
                                                 const std::vector<PhasePresence>& cond) const
    {
        const Opm::PhaseUsage& pu = fluid_.phaseUsage();
        if (phase == Solvent) {
            if (!is_miscible_) {
                return solvent_props_.bSolvent(p, cells_);
            } else {
                return b_eff_[solvent_pos_];
            }

        } else {
            if (!is_miscible_) {
                return Base::fluidReciprocFVF(phase, p, temp, rs, rv, cond);
            } else {
                return b_eff_[pu.phase_pos[ phase ]];
            }
        }
    }

    template <class Grid>
    ADB
    BlackoilSolventModel<Grid>::fluidDensity(const int  phase,
                                                          const ADB& b,
                                                          const ADB& rs,
                                                          const ADB& rv) const
    {
        if (phase == Solvent && has_solvent_) {
            return solvent_props_.solventSurfaceDensity(cells_) * b;
        } else {
            return Base::fluidDensity(phase, b, rs, rv);
        }
    }

    template <class Grid>
    std::vector<ADB>
    BlackoilSolventModel<Grid>::computeRelPerm(const SolutionState& state) const
    {
        using namespace Opm::AutoDiffGrid;
        const int               nc   = numCells(grid_);

        const ADB zero = ADB::constant(V::Zero(nc));

        const Opm::PhaseUsage& pu = fluid_.phaseUsage();
        const ADB& sw = (active_[ Water ]
                         ? state.saturation[ pu.phase_pos[ Water ] ]
                         : zero);

        const ADB& so = (active_[ Oil ]
                         ? state.saturation[ pu.phase_pos[ Oil ] ]
                         : zero);

        const ADB& sg = (active_[ Gas ]
                         ? state.saturation[ pu.phase_pos[ Gas ] ]
                         : zero);

        if (has_solvent_) {
            const ADB& ss = state.solvent_saturation;
            if (is_miscible_) {

                std::vector<ADB> relperm = fluid_.relperm(sw, so, sg+ss, cells_);

                Selector<double> zero_selector(ss.value() + sg.value(), Selector<double>::Zero);
                ADB F_solvent = zero_selector.select(ss, ss / (ss + sg));
                const ADB& po = state.canonical_phase_pressures[ Oil ];
                const ADB misc = solvent_props_.miscibilityFunction(F_solvent, cells_)
                        * solvent_props_.pressureMiscibilityFunction(po, cells_);

                assert(active_[ Oil ]);
                assert(active_[ Gas ]);

                const ADB sn = ss + so + sg;

                // adjust endpoints
                const V sgcr = fluid_.scaledCriticalGasSaturations(cells_);
                const V sogcr = fluid_.scaledCriticalOilinGasSaturations(cells_);
                const ADB sorwmis = solvent_props_.miscibleResidualOilSaturationFunction(sw, cells_);
                const ADB sgcwmis = solvent_props_.miscibleCriticalGasSaturationFunction(sw, cells_);

                const V ones = V::Constant(nc, 1.0);
                ADB sor = misc * sorwmis + (ones - misc) * sogcr;
                ADB sgc = misc * sgcwmis + (ones - misc) * sgcr;

                // avoid negative values
                Selector<double> negative_selector(sgc.value(), Selector<double>::LessEqualZero);
                sgc = negative_selector.select(zero, sgc);

                const ADB ssg = ss + sg - sgc;
                const ADB sn_eff = sn - sor - sgc;

                // avoid negative values
                Selector<double> zeroSn_selector(sn_eff.value(), Selector<double>::Zero);
                const ADB F_totalGas = zeroSn_selector.select( zero, ssg / sn_eff);

                const ADB mkrgt = solvent_props_.miscibleSolventGasRelPermMultiplier(F_totalGas, cells_) * solvent_props_.misicibleHydrocarbonWaterRelPerm(sn, cells_);
                const ADB mkro = solvent_props_.miscibleOilRelPermMultiplier(ones - F_totalGas, cells_) * solvent_props_.misicibleHydrocarbonWaterRelPerm(sn, cells_);

                const V eps = V::Constant(nc, 1e-7);
                Selector<double> noOil_selector(so.value()-eps, Selector<double>::LessEqualZero);
                relperm[Gas] = noOil_selector.select(relperm[Gas], (ones - misc) * relperm[Gas] + misc * mkrgt);
                relperm[Oil] = noOil_selector.select(relperm[Oil], (ones - misc) * relperm[Oil] + misc * mkro);
                return relperm;
            } else {
                return fluid_.relperm(sw, so, sg+ss, cells_);
            }
        } else {
            return fluid_.relperm(sw, so, sg, cells_);
        }

    }

    template <class Grid>
    void
    BlackoilSolventModel<Grid>::computeEffectiveProperties(const SolutionState&    state)
    {
        // Viscosity
        const Opm::PhaseUsage& pu = fluid_.phaseUsage();
        const int np = fluid_.numPhases();
        const int nc   = Opm::UgGridHelpers::numCells(grid_);
        const ADB zero = ADB::constant(V::Zero(nc));

        const ADB& pw = state.canonical_phase_pressures[pu.phase_pos[Water]];
        const ADB& po = state.canonical_phase_pressures[pu.phase_pos[Oil]];
        const ADB& pg = state.canonical_phase_pressures[pu.phase_pos[Gas]];
        const std::vector<PhasePresence>& cond = phaseCondition();

        const ADB mu_w = fluid_.muWat(pw, state.temperature, cells_);
        const ADB mu_o = fluid_.muOil(po, state.temperature, state.rs, cond, cells_);
        const ADB mu_g = fluid_.muGas(pg, state.temperature, state.rv, cond, cells_);
        const ADB mu_s = solvent_props_.muSolvent(pg,cells_);
        std::vector<ADB> viscosity(np + 1, ADB::null());
        viscosity[pu.phase_pos[Oil]] = mu_o;
        viscosity[pu.phase_pos[Gas]] = mu_g;
        viscosity[pu.phase_pos[Water]] = mu_w;
        viscosity[solvent_pos_] = mu_s;

        // Density
        const ADB bw = fluid_.bWat(pw, state.temperature, cells_);
        const ADB bo = fluid_.bOil(po, state.temperature, state.rs, cond, cells_);
        const ADB bg = fluid_.bGas(pg, state.temperature, state.rv, cond, cells_);
        const ADB bs = solvent_props_.bSolvent(pg, cells_);

        std::vector<ADB> density(np + 1, ADB::null());
        density[pu.phase_pos[Oil]] = fluidDensity(Oil, bo, state.rs, state.rv);
        density[pu.phase_pos[Gas]] = fluidDensity(Gas, bg, state.rs, state.rv);
        density[pu.phase_pos[Water]] = fluidDensity(Water, bw, state.rs, state.rv);
        density[solvent_pos_] = fluidDensity(Solvent, bs, state.rs, state.rv);

        // Effective saturations
        const ADB& ss = state.solvent_saturation;
        const ADB& so = state.saturation[ pu.phase_pos[ Oil ] ];
        const ADB& sg = (active_[ Gas ]
                         ? state.saturation[ pu.phase_pos[ Gas ] ]
                         : zero);
        const ADB& sw = (active_[ Water ]
                         ? state.saturation[ pu.phase_pos[ Water ] ]
                         : zero);

        const ADB sorwmis = solvent_props_.miscibleResidualOilSaturationFunction(sw, cells_);
        const ADB sgcwmis = solvent_props_.miscibleCriticalGasSaturationFunction(sw, cells_);

        std::vector<ADB> effective_saturations (np + 1, ADB::null());
        effective_saturations[pu.phase_pos[Oil]] = so - sorwmis;
        effective_saturations[pu.phase_pos[Gas]] = sg - sgcwmis;
        effective_saturations[pu.phase_pos[Water]] = sw;
        effective_saturations[solvent_pos_] = ss - sgcwmis;

        // Compute effective viscosities and densities
        computeToddLongstaffMixing(viscosity, density, effective_saturations, pu);

        // compute the volume factors from the densities
        const ADB b_eff_o = density[pu.phase_pos[ Oil ]] / (fluid_.surfaceDensity(pu.phase_pos[ Oil ],  cells_) + fluid_.surfaceDensity(pu.phase_pos[ Gas ], cells_) * state.rs);
        const ADB b_eff_g = density[pu.phase_pos[ Gas ]] / (fluid_.surfaceDensity(pu.phase_pos[ Gas ],  cells_) + fluid_.surfaceDensity(pu.phase_pos[ Oil ], cells_) * state.rv);
        const ADB b_eff_s = density[solvent_pos_] / solvent_props_.solventSurfaceDensity(cells_);

        // account for pressure effects and store the computed volume factors and viscosities
        const V ones = V::Constant(nc, 1.0);
        const ADB pmisc = solvent_props_.pressureMiscibilityFunction(po, cells_);

        b_eff_[pu.phase_pos[ Oil ]] = pmisc * b_eff_o + (ones - pmisc) * bo;
        b_eff_[pu.phase_pos[ Gas ]] = pmisc * b_eff_g + (ones - pmisc) * bg;
        b_eff_[solvent_pos_] = pmisc * b_eff_s + (ones - pmisc) * bs;

        // keep the mu*b interpolation
        mu_eff_[pu.phase_pos[ Oil ]] = b_eff_[pu.phase_pos[ Oil ]] / (pmisc * b_eff_o / viscosity[pu.phase_pos[ Oil ]] + (ones - pmisc) * bo / mu_o);
        mu_eff_[pu.phase_pos[ Gas ]] = b_eff_[pu.phase_pos[ Gas ]] / (pmisc * b_eff_g / viscosity[pu.phase_pos[ Gas ]] + (ones - pmisc) * bg / mu_g);
        mu_eff_[solvent_pos_] = b_eff_[solvent_pos_] / (pmisc * b_eff_s / viscosity[solvent_pos_] + (ones - pmisc) * bs / mu_s);

        // for water the pure values are used
        mu_eff_[pu.phase_pos[ Water ]] = mu_w;
        b_eff_[pu.phase_pos[ Water ]] = bw;
    }

    template <class Grid>
    void
    BlackoilSolventModel<Grid>::computeToddLongstaffMixing(std::vector<ADB>& viscosity, std::vector<ADB>& density, const std::vector<ADB>& saturations, const Opm::PhaseUsage pu)
    {
        const int  nc = cells_.size();
        const V ones = V::Constant(nc, 1.0);
        const ADB zero = ADB::constant(V::Zero(nc));

        // Calculation of effective saturations
        ADB so_eff = saturations[pu.phase_pos[ Oil ]];
        ADB sg_eff = saturations[pu.phase_pos[ Gas ]];
        ADB ss_eff = saturations[solvent_pos_];

        // Avoid negative values
        Selector<double> negative_selectorSo(so_eff.value(), Selector<double>::LessZero);
        Selector<double> negative_selectorSg(sg_eff.value(), Selector<double>::LessZero);
        Selector<double> negative_selectorSs(ss_eff.value(), Selector<double>::LessZero);
        so_eff = negative_selectorSo.select(zero, so_eff);
        sg_eff = negative_selectorSg.select(zero, sg_eff);
        ss_eff = negative_selectorSs.select(zero, ss_eff);

        // Make copy of the pure viscosities
        const ADB mu_o = viscosity[pu.phase_pos[ Oil ]];
        const ADB mu_g = viscosity[pu.phase_pos[ Gas ]];
        const ADB mu_s = viscosity[solvent_pos_];

        const ADB sn_eff =  so_eff + sg_eff + ss_eff;
        const ADB sos_eff = so_eff + ss_eff;
        const ADB ssg_eff = ss_eff + sg_eff;

        // Avoid division by zero
        Selector<double> zero_selectorSos(sos_eff.value(), Selector<double>::Zero);
        Selector<double> zero_selectorSsg(ssg_eff.value(), Selector<double>::Zero);
        Selector<double> zero_selectorSn(sn_eff.value(), Selector<double>::Zero);

        const ADB mu_s_pow = pow(mu_s, 0.25);
        const ADB mu_o_pow = pow(mu_o, 0.25);
        const ADB mu_g_pow = pow(mu_g, 0.25);

        const ADB mu_mos = zero_selectorSos.select(mu_o + mu_s, mu_o * mu_s / pow( ( (so_eff / sos_eff) * mu_s_pow) + ( (ss_eff / sos_eff) * mu_o_pow) , 4.0));
        const ADB mu_msg = zero_selectorSsg.select(mu_g + mu_s , mu_g * mu_s / pow( ( (sg_eff / ssg_eff) * mu_s_pow) + ( (ss_eff / ssg_eff) * mu_g_pow) , 4.0));
        const ADB mu_m = zero_selectorSn.select(mu_s + mu_o + mu_g, mu_o * mu_s * mu_g / pow( ( (so_eff / sn_eff) * mu_s_pow *  mu_g_pow)
                                                       + ( (ss_eff / sn_eff) * mu_o_pow *  mu_g_pow) + ( (sg_eff / sn_eff) * mu_s_pow * mu_o_pow), 4.0));
        // Mixing parameter for viscosity
        const V mix_param_mu = solvent_props_.mixingParameterViscosity(cells_);

        // Update viscosities
        viscosity[pu.phase_pos[ Oil ]] = pow(mu_o,ones - mix_param_mu) * pow(mu_mos, mix_param_mu);
        viscosity[pu.phase_pos[ Gas ]] = pow(mu_g,ones - mix_param_mu) * pow(mu_msg, mix_param_mu);
        viscosity[solvent_pos_] = pow(mu_s,ones - mix_param_mu) * pow(mu_m, mix_param_mu);

        // Density
        ADB& rho_o = density[pu.phase_pos[ Oil ]];
        ADB& rho_g = density[pu.phase_pos[ Gas ]];
        ADB& rho_s = density[solvent_pos_];

        // mixing parameter for density
        const V mix_param_rho = solvent_props_.mixingParameterDensity(cells_);

        // compute effective viscosities for density calculations. These have to
        // be recomputed as a different mixing parameter may be used.
        const ADB mu_o_eff = pow(mu_o,ones - mix_param_rho) * pow(mu_mos, mix_param_rho);
        const ADB mu_g_eff = pow(mu_g,ones - mix_param_rho) * pow(mu_msg, mix_param_rho);
        const ADB mu_s_eff = pow(mu_s,ones - mix_param_rho) * pow(mu_m, mix_param_rho);

        const ADB sog_eff = so_eff + sg_eff;
        const ADB sof = so_eff / sog_eff;
        const ADB sgf = sg_eff / sog_eff;

        // Effective densities
        const ADB mu_sog_pow = mu_s_pow * ( (sgf * mu_o_pow) + (sof * mu_g_pow) );
        const ADB mu_o_eff_pow = pow(mu_o_eff, 0.25);
        const ADB mu_g_eff_pow = pow(mu_g_eff, 0.25);
        const ADB mu_s_eff_pow = pow(mu_s_eff, 0.25);
        const ADB sfraction_oe = (mu_o_pow * (mu_o_eff_pow - mu_s_pow)) / (mu_o_eff_pow * (mu_o_pow - mu_s_pow));
        const ADB sfraction_ge = (mu_g_pow * (mu_s_pow - mu_g_eff_pow)) / (mu_g_eff_pow * (mu_s_pow - mu_g_pow));
        const ADB sfraction_se = (mu_sog_pow - ( mu_o_pow * mu_g_pow * mu_s_pow / mu_s_eff_pow) ) / ( mu_sog_pow - (mu_o_pow * mu_g_pow));
        const ADB rho_o_eff = (rho_o * sfraction_oe) + (rho_s * (ones - sfraction_oe));
        const ADB rho_g_eff = (rho_g * sfraction_ge) + (rho_s * (ones - sfraction_ge));
        const ADB rho_s_eff = (rho_s * sfraction_se) + (rho_g * sgf * (ones - sfraction_se)) + (rho_o * sof * (ones - sfraction_se));

        // Avoid division by zero for equal mobilities. For equal mobilities the effecitive density is calculated
        // based on the saturation fraction directly.
        Selector<double> unitGasSolventMobilityRatio_selector(mu_s.value() - mu_g.value(), Selector<double>::Zero);
        Selector<double> unitOilSolventMobilityRatio_selector(mu_s.value() - mu_o.value(), Selector<double>::Zero);

        // Effective densities when the mobilities are equal
        const ADB rho_m = zero_selectorSn.select(zero, (rho_o * so_eff / sn_eff) + (rho_g * sg_eff / sn_eff) + (rho_s * ss_eff / sn_eff));
        const ADB rho_o_eff_simple = ((ones - mix_param_rho) * rho_o) + (mix_param_rho * rho_m);
        const ADB rho_g_eff_simple = ((ones - mix_param_rho) * rho_g) + (mix_param_rho * rho_m);
        //const ADB rho_s_eff_simple = ((ones - mix_param_rho) * rho_s) + (mix_param_rho * rho_m);

        // Update densities
        rho_o = unitOilSolventMobilityRatio_selector.select(rho_o_eff_simple, rho_o_eff);
        rho_g = unitGasSolventMobilityRatio_selector.select(rho_g_eff_simple, rho_g_eff);
        rho_s = rho_s_eff;

    }


    template <class Grid>
    std::vector<ADB>
    BlackoilSolventModel<Grid>::
    computePressures(const ADB& po,
                     const ADB& sw,
                     const ADB& so,
                     const ADB& sg,
                     const ADB& ss) const
    {
        std::vector<ADB> pressures = Base::computePressures(po, sw, so, sg);
        // The imiscible capillary pressure is evaluated using the total gas saturation (sg + ss)
        std::vector<ADB> pressures_imisc = Base::computePressures(po, sw, so, sg + ss);

        // Pressure effects on capillary pressure miscibility
        const ADB pmisc = solvent_props_.pressureMiscibilityFunction(po, cells_);
        // Only the pcog is effected by the miscibility. Since pg = po + pcog, changing pg is eqvivalent
        // to changing the gas pressure directly.
        const int  nc = cells_.size();
        const V ones = V::Constant(nc, 1.0);
        pressures[Gas] = ( pmisc * pressures[Gas] + ((ones - pmisc) * pressures_imisc[Gas]));
        return pressures;
    }


    template <class Grid>
    void
    BlackoilSolventModel<Grid>::assemble(const ReservoirState& reservoir_state,
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
            if (is_miscible_) {
                computeEffectiveProperties(state0);
            }

            computeAccum(state0, 0);
            computeWellConnectionPressures(state0, well_state);
        }
        if (is_miscible_) {
            computeEffectiveProperties(state);
        }

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

        if (has_solvent_) {
            int gas_pos = fluid_.phaseUsage().phase_pos[Gas];
            // Gas and solvent is combinded and solved together
            // The input in the well equation is then the
            // total gas phase = hydro carbon gas + solvent gas

            // The total mobility is the sum of the solvent and gas mobiliy
            mob_perfcells[gas_pos] += subset(rq_[solvent_pos_].mob, well_cells);

            // A weighted sum of the b-factors of gas and solvent are used.
            const int nc = Opm::AutoDiffGrid::numCells(grid_);

            const Opm::PhaseUsage& pu = fluid_.phaseUsage();
            const ADB zero = ADB::constant(V::Zero(nc));
            const ADB& ss = state.solvent_saturation;
            const ADB& sg = (active_[ Gas ]
                             ? state.saturation[ pu.phase_pos[ Gas ] ]
                             : zero);

            Selector<double> zero_selector(ss.value() + sg.value(), Selector<double>::Zero);
            ADB F_solvent = subset(zero_selector.select(ss, ss / (ss + sg)),well_cells);
            V ones = V::Constant(nperf,1.0);

            b_perfcells[gas_pos] = (ones - F_solvent) * b_perfcells[gas_pos];
            b_perfcells[gas_pos] += (F_solvent * subset(rq_[solvent_pos_].b, well_cells));

        }
        if (param_.solve_welleq_initially_ && initial_assembly) {
            // solve the well equations as a pre-processing step
            Base::solveWellEq(mob_perfcells, b_perfcells, state, well_state);
        }
        Base::computeWellFlux(state, mob_perfcells, b_perfcells, aliveWells, cq_s);
        Base::updatePerfPhaseRatesAndPressures(cq_s, state, well_state);
        Base::addWellFluxEq(cq_s, state);
        addWellContributionToMassBalanceEq(cq_s, state, well_state);
        Base::addWellControlEq(state, well_state, aliveWells);

    }


}


#endif // OPM_BLACKOILSOLVENT_IMPL_HEADER_INCLUDED
