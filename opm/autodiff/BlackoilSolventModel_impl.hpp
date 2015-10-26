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
                                                     const bool                              has_solvent)
        : Base(param, grid, fluid, geo, rock_comp_props, wells_arg, linsolver,
               eclState, has_disgas, has_vapoil, terminal_output),
          has_solvent_(has_solvent),
          solvent_pos_(detail::solventPos(fluid.phaseUsage())),
          solvent_props_(solvent_props)
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
            assert (not x.solvent_saturation().empty());
            const int nc = x.solvent_saturation().size();
            const V ss = Eigen::Map<const V>(&x.solvent_saturation()[0], nc);
            // Solvent belongs after other reservoir vars but before well vars.
            auto solvent_pos = vars0.begin() + fluid_.numPhases();
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
        SolutionState state = Base::variableStateExtractVars(x, indices, vars);
        if (has_solvent_) {
            state.solvent_saturation = std::move(vars[indices[Solvent]]);
            if (active_[ Oil ]) {
                // Note that so is never a primary variable.
                const Opm::PhaseUsage pu = fluid_.phaseUsage();
                state.saturation[pu.phase_pos[ Oil ]] -= state.solvent_saturation;
            }
        }
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
            rq_[solvent_pos_].b = solvent_props_.bSolvent(pg,cells_);
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

        // Surface density.
        std::vector<V> surf_dens = fluid_.surfaceDensity(well_cells);

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
            V bg = fluid_.bGas(avg_press_ad, perf_temp, perf_rv, perf_cond, well_cells).value();

            if (has_solvent_) {
                const V bs = solvent_props_.bSolvent(avg_press_ad,well_cells).value();
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

                const V& rhog = surf_dens[pu.phase_pos[BlackoilPhases::Vapour]];
                const V& rhos = solvent_props_.solventSurfaceDensity(well_cells);
                surf_dens[pu.phase_pos[BlackoilPhases::Vapour]] = ( (ones - F_solvent) * rhog ) + (F_solvent * rhos);
            }
            b.col(pu.phase_pos[BlackoilPhases::Vapour]) = bg;

            const V rvsat = fluidRvSat(avg_press, perf_so, well_cells);
            rvmax_perf.assign(rvsat.data(), rvsat.data() + nperf);
        }

        // b and surf_dens_perf is row major, so can just copy data.
        V surf_dens_copy = superset(surf_dens[0], Span(nperf, pu.num_phases, 0), nperf*pu.num_phases);
        for (int phase = 1; phase < pu.num_phases; ++phase) {
            surf_dens_copy += superset(surf_dens[phase], Span(nperf, pu.num_phases, phase), nperf*pu.num_phases);
        }

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
            const V ss_old = Eigen::Map<const V>(&reservoir_state.solvent_saturation()[0], nc, 1);
            const V ss = (ss_old - dss).max(zero);
            std::copy(&ss[0], &ss[0] + nc, reservoir_state.solvent_saturation().begin());

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
                                                const ADB&              phasePressure,
                                                const SolutionState&    state)
    {
        Base::computeMassFlux(actph, transi, kr, phasePressure, state);

        const int canonicalPhaseIdx = canph_[ actph ];
        if (canonicalPhaseIdx == Gas) {
            if (has_solvent_) {
                const int  nc   = Opm::UgGridHelpers::numCells(grid_);

                const Opm::PhaseUsage& pu = fluid_.phaseUsage();
                const ADB zero = ADB::constant(V::Zero(nc));
                const ADB& ss = state.solvent_saturation;
                const ADB& sg = (active_[ Gas ]
                                 ? state.saturation[ pu.phase_pos[ Gas ] ]
                                 : zero);

                Selector<double> zero_selector(ss.value() + sg.value(), Selector<double>::Zero);
                ADB F_solvent = zero_selector.select(ss, ss / (ss + sg));
                V ones = V::Constant(nc, 1.0);

                const ADB tr_mult = transMult(state.pressure);
                const ADB mu = solvent_props_.muSolvent(phasePressure,cells_);               

                rq_[solvent_pos_].mob = solvent_props_.solventRelPermMultiplier(F_solvent, cells_) * tr_mult * kr / mu;

                const ADB rho_solvent = solvent_props_.solventSurfaceDensity(cells_) * rq_[solvent_pos_].b;
                const ADB rhoavg_solvent = ops_.caver * rho_solvent;
                rq_[ solvent_pos_ ].dh = ops_.ngrad * phasePressure - geo_.gravity()[2] * (rhoavg_solvent * (ops_.ngrad * geo_.z().matrix()));

                UpwindSelector<double> upwind_solvent(grid_, ops_, rq_[solvent_pos_].dh.value());
                // Compute solvent flux.
                rq_[solvent_pos_].mflux = upwind_solvent.select(rq_[solvent_pos_].b * rq_[solvent_pos_].mob) * (transi * rq_[solvent_pos_].dh);

                // Update gas mobility and flux
                rq_[actph].mob = solvent_props_.gasRelPermMultiplier( (ones - F_solvent) , cells_) * rq_[actph].mob;

                const ADB& b   = rq_[ actph ].b;
                const ADB& mob = rq_[ actph ].mob;
                const ADB& dh  = rq_[ actph ].dh;
                UpwindSelector<double> upwind_gas(grid_, ops_, dh.value());
                rq_[ actph ].mflux = upwind_gas.select(b * mob) * (transi * dh);
            }
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
            return fluid_.relperm(sw, so, sg+ss, cells_);
        } else {
            return fluid_.relperm(sw, so, sg, cells_);
        }

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
            computeAccum(state0, 0);
            computeWellConnectionPressures(state0, well_state);
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
