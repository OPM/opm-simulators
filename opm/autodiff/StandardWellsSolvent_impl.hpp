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


#include <opm/autodiff/StandardWellsSolvent.hpp>




namespace Opm
{






    StandardWellsSolvent::StandardWellsSolvent(const Wells* wells_arg, WellCollection* well_collection)
        : Base(wells_arg, well_collection)
        , solvent_props_(nullptr)
        , solvent_pos_(-1)
        , has_solvent_(false)
    {
    }





    void
    StandardWellsSolvent::initSolvent(const SolventPropsAdFromDeck* solvent_props,
                                      const int solvent_pos,
                                      const bool has_solvent)
    {
        solvent_props_ = solvent_props;
        solvent_pos_ = solvent_pos;
        has_solvent_ = has_solvent;
    }





    template<class SolutionState, class WellState>
    void
    StandardWellsSolvent::
    computePropertiesForWellConnectionPressures(const SolutionState& state,
                                                const WellState& xw,
                                                std::vector<double>& b_perf,
                                                std::vector<double>& rsmax_perf,
                                                std::vector<double>& rvmax_perf,
                                                std::vector<double>& surf_dens_perf)
    {
        // 1. Compute properties required by computeConnectionPressureDelta().
        //    Note that some of the complexity of this part is due to the function
        //    taking std::vector<double> arguments, and not Eigen objects.
        const int nperf = wells().well_connpos[wells().number_of_wells];
        const int nw = wells().number_of_wells;

        // Compute the average pressure in each well block
        const Vector perf_press = Eigen::Map<const V>(xw.perfPress().data(), nperf);
        Vector avg_press = perf_press*0;
        for (int w = 0; w < nw; ++w) {
            for (int perf = wells().well_connpos[w]; perf < wells().well_connpos[w+1]; ++perf) {
                const double p_above = perf == wells().well_connpos[w] ? state.bhp.value()[w] : perf_press[perf - 1];
                const double p_avg = (perf_press[perf] + p_above)/2;
                avg_press[perf] = p_avg;
            }
        }

        const std::vector<int>& well_cells = wellOps().well_cells;

        // Use cell values for the temperature as the wells don't knows its temperature yet.
        const ADB perf_temp = subset(state.temperature, well_cells);

        // Compute b, rsmax, rvmax values for perforations.
        // Evaluate the properties using average well block pressures
        // and cell values for rs, rv, phase condition and temperature.
        const ADB avg_press_ad = ADB::constant(avg_press);
        std::vector<PhasePresence> perf_cond(nperf);
        for (int perf = 0; perf < nperf; ++perf) {
            perf_cond[perf] = (*phase_condition_)[well_cells[perf]];
        }

        const PhaseUsage& pu = fluid_->phaseUsage();
        DataBlock b(nperf, pu.num_phases);

        const Vector bw = fluid_->bWat(avg_press_ad, perf_temp, well_cells).value();
        if (pu.phase_used[BlackoilPhases::Aqua]) {
            b.col(pu.phase_pos[BlackoilPhases::Aqua]) = bw;
        }

        assert((*active_)[Oil]);
        assert((*active_)[Gas]);
        const ADB perf_rv = subset(state.rv, well_cells);
        const ADB perf_rs = subset(state.rs, well_cells);
        const Vector perf_so =  subset(state.saturation[pu.phase_pos[Oil]].value(), well_cells);
        if (pu.phase_used[BlackoilPhases::Liquid]) {
            const Vector bo = fluid_->bOil(avg_press_ad, perf_temp, perf_rs, perf_cond, well_cells).value();
            //const V bo_eff = subset(rq_[pu.phase_pos[Oil] ].b , well_cells).value();
            b.col(pu.phase_pos[BlackoilPhases::Liquid]) = bo;
            // const Vector rssat = fluidRsSat(avg_press, perf_so, well_cells);
            const Vector rssat = fluid_->rsSat(ADB::constant(avg_press), ADB::constant(perf_so), well_cells).value();
            rsmax_perf.assign(rssat.data(), rssat.data() + nperf);
        } else {
            rsmax_perf.assign(0.0, nperf);
        }
        V surf_dens_copy = superset(fluid_->surfaceDensity(0, well_cells), Span(nperf, pu.num_phases, 0), nperf*pu.num_phases);
        for (int phase = 1; phase < pu.num_phases; ++phase) {
            if ( phase == pu.phase_pos[BlackoilPhases::Vapour]) {
                continue; // the gas surface density is added after the solvent is accounted for.
            }
            surf_dens_copy += superset(fluid_->surfaceDensity(phase, well_cells), Span(nperf, pu.num_phases, phase), nperf*pu.num_phases);
        }

        if (pu.phase_used[BlackoilPhases::Vapour]) {
            // Unclear wether the effective or the pure values should be used for the wells
            // the current usage of unmodified properties values gives best match.
            //V bg_eff = subset(rq_[pu.phase_pos[Gas]].b,well_cells).value();
            Vector bg = fluid_->bGas(avg_press_ad, perf_temp, perf_rv, perf_cond, well_cells).value();
            Vector rhog = fluid_->surfaceDensity(pu.phase_pos[BlackoilPhases::Vapour], well_cells);
            // to handle solvent related
            if (has_solvent_) {

                const Vector bs = solvent_props_->bSolvent(avg_press_ad,well_cells).value();
                //const V bs_eff = subset(rq_[solvent_pos_].b,well_cells).value();

                // number of cells
                const int nc = state.pressure.size();

                const ADB zero = ADB::constant(Vector::Zero(nc));
                const ADB& ss = state.solvent_saturation;
                const ADB& sg = ((*active_)[ Gas ]
                                 ? state.saturation[ pu.phase_pos[ Gas ] ]
                                 : zero);

                Selector<double> zero_selector(ss.value() + sg.value(), Selector<double>::Zero);
                Vector F_solvent = subset(zero_selector.select(ss, ss / (ss + sg)),well_cells).value();

                Vector injectedSolventFraction = Eigen::Map<const Vector>(&xw.solventFraction()[0], nperf);

                Vector isProducer = Vector::Zero(nperf);
                Vector ones = Vector::Constant(nperf,1.0);
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

                const Vector& rhos = solvent_props_->solventSurfaceDensity(well_cells);
                rhog = ( (ones - F_solvent) * rhog ) + (F_solvent * rhos);
            }
            b.col(pu.phase_pos[BlackoilPhases::Vapour]) = bg;
            surf_dens_copy += superset(rhog, Span(nperf, pu.num_phases, pu.phase_pos[BlackoilPhases::Vapour]), nperf*pu.num_phases);

            // const Vector rvsat = fluidRvSat(avg_press, perf_so, well_cells);
            const Vector rvsat = fluid_->rvSat(ADB::constant(avg_press), ADB::constant(perf_so), well_cells).value();
            rvmax_perf.assign(rvsat.data(), rvsat.data() + nperf);
        } else {
            rvmax_perf.assign(0.0, nperf);
        }

        // b and surf_dens_perf is row major, so can just copy data.
        b_perf.assign(b.data(), b.data() + nperf * pu.num_phases);
        surf_dens_perf.assign(surf_dens_copy.data(), surf_dens_copy.data() + nperf * pu.num_phases);
    }

    template <class SolutionState>
    void
    StandardWellsSolvent::
    computeWellFlux(const SolutionState& state,
                    const std::vector<ADB>& mob_perfcells,
                    const std::vector<ADB>& b_perfcells,
                    Vector& aliveWells,
                    std::vector<ADB>& cq_s) const
    {
        if( ! localWellsActive() ) return ;

        const int np = wells().number_of_phases;
        const int nw = wells().number_of_wells;
        const int nperf = wells().well_connpos[nw];
        Vector Tw = Eigen::Map<const Vector>(wells().WI, nperf);
        const std::vector<int>& well_cells = wellOps().well_cells;

        // pressure diffs computed already (once per step, not changing per iteration)
        const Vector& cdp = wellPerforationPressureDiffs();
        // Extract needed quantities for the perforation cells
        const ADB& p_perfcells = subset(state.pressure, well_cells);

        // Perforation pressure
        const ADB perfpressure = (wellOps().w2p * state.bhp) + cdp;

        // Pressure drawdown (also used to determine direction of flow)
        const ADB drawdown =  p_perfcells - perfpressure;

        // Compute vectors with zero and ones that
        // selects the wanted quantities.

        // selects injection perforations
        Vector selectInjectingPerforations = Vector::Zero(nperf);
        // selects producing perforations
        Vector selectProducingPerforations = Vector::Zero(nperf);
        for (int c = 0; c < nperf; ++c){
            if (drawdown.value()[c] < 0)
                selectInjectingPerforations[c] = 1;
            else
                selectProducingPerforations[c] = 1;
        }

        // Handle cross flow
        const Vector numInjectingPerforations = (wellOps().p2w * ADB::constant(selectInjectingPerforations)).value();
        const Vector numProducingPerforations = (wellOps().p2w * ADB::constant(selectProducingPerforations)).value();
        for (int w = 0; w < nw; ++w) {
            if (!wells().allow_cf[w]) {
                for (int perf = wells().well_connpos[w] ; perf < wells().well_connpos[w+1]; ++perf) {
                    // Crossflow is not allowed; reverse flow is prevented.
                    // At least one of the perforation must be open in order to have a meeningful
                    // equation to solve. For the special case where all perforations have reverse flow,
                    // and the target rate is non-zero all of the perforations are keept open.
                    if (wells().type[w] == INJECTOR && numInjectingPerforations[w] > 0) {
                        selectProducingPerforations[perf] = 0.0;
                    } else if (wells().type[w] == PRODUCER && numProducingPerforations[w] > 0 ){
                        selectInjectingPerforations[perf] = 0.0;
                    }
                }
            }
        }

        // HANDLE FLOW INTO WELLBORE
        // compute phase volumetric rates at standard conditions
        std::vector<ADB> cq_p(np, ADB::null());
        std::vector<ADB> cq_ps(np, ADB::null());
        for (int phase = 0; phase < np; ++phase) {
            cq_p[phase] = -(selectProducingPerforations * Tw) * (mob_perfcells[phase] * drawdown);
            cq_ps[phase] = b_perfcells[phase] * cq_p[phase];
        }

        Vector ones = Vector::Constant(nperf,1.0);
        ADB F_gas = ADB::constant(ones);

        const Opm::PhaseUsage& pu = fluid_->phaseUsage();
        if ((*active_)[Oil] && (*active_)[Gas]) {
            const int oilpos = pu.phase_pos[Oil];
            const int gaspos = pu.phase_pos[Gas];
            const ADB cq_psOil = cq_ps[oilpos];
            ADB cq_psGas = cq_ps[gaspos];
            const ADB& rv_perfcells = subset(state.rv, well_cells);
            const ADB& rs_perfcells = subset(state.rs, well_cells);
            cq_ps[gaspos] += rs_perfcells * cq_psOil;

            if(has_solvent_) {
                // The solvent gas need to be removed from the gas
                // before multiplied with rv.
                const ADB& ss = state.solvent_saturation;
                const ADB& sg = state.saturation[ pu.phase_pos[ Gas ] ];

                Selector<double> zero_selector(ss.value() + sg.value(), Selector<double>::Zero);
                F_gas -= subset(zero_selector.select(ss, ss / (ss + sg)),well_cells);
                cq_psGas = cq_psGas * F_gas;
            }
            cq_ps[oilpos] += rv_perfcells * cq_psGas;
        }

        // HANDLE FLOW OUT FROM WELLBORE
        // Using total mobilities
        ADB total_mob = mob_perfcells[0];
        for (int phase = 1; phase < np; ++phase) {
            total_mob += mob_perfcells[phase];
        }
        // injection perforations total volume rates
        const ADB cqt_i = -(selectInjectingPerforations * Tw) * (total_mob * drawdown);

        // Store well perforation total fluxes (reservor volumes) if requested.
        if (store_well_perforation_fluxes_) {
            // Ugly const-cast, but unappealing alternatives.
            Vector& wf = const_cast<Vector&>(well_perforation_fluxes_);
            wf = cqt_i.value();
            for (int phase = 0; phase < np; ++phase) {
                wf += cq_p[phase].value();
            }
        }

        // compute wellbore mixture for injecting perforations
        // The wellbore mixture depends on the inflow from the reservoar
        // and the well injection rates.

        // compute avg. and total wellbore phase volumetric rates at standard conds
        const DataBlock compi = Eigen::Map<const DataBlock>(wells().comp_frac, nw, np);
        std::vector<ADB> wbq(np, ADB::null());
        ADB wbqt = ADB::constant(Vector::Zero(nw));
        for (int phase = 0; phase < np; ++phase) {
            const ADB& q_ps = wellOps().p2w * cq_ps[phase];
            const ADB& q_s = subset(state.qs, Span(nw, 1, phase*nw));
            Selector<double> injectingPhase_selector(q_s.value(), Selector<double>::GreaterZero);
            const int pos = pu.phase_pos[phase];
            wbq[phase] = (compi.col(pos) * injectingPhase_selector.select(q_s,ADB::constant(Vector::Zero(nw))))  - q_ps;
            wbqt += wbq[phase];
        }
        // compute wellbore mixture at standard conditions.
        Selector<double> notDeadWells_selector(wbqt.value(), Selector<double>::Zero);
        std::vector<ADB> cmix_s(np, ADB::null());
        for (int phase = 0; phase < np; ++phase) {
            const int pos = pu.phase_pos[phase];
            cmix_s[phase] = wellOps().w2p * notDeadWells_selector.select(ADB::constant(compi.col(pos)), wbq[phase]/wbqt);
        }

        // compute volume ratio between connection at standard conditions
        ADB volumeRatio = ADB::constant(Vector::Zero(nperf));

        if ((*active_)[Water]) {
            const int watpos = pu.phase_pos[Water];
            volumeRatio += cmix_s[watpos] / b_perfcells[watpos];
        }

        if ((*active_)[Oil] && (*active_)[Gas]) {
            // Incorporate RS/RV factors if both oil and gas active
            const ADB& rv_perfcells = subset(state.rv, well_cells);
            const ADB& rs_perfcells = subset(state.rs, well_cells);
            const ADB d = Vector::Constant(nperf,1.0) - rv_perfcells * rs_perfcells;

            const int oilpos = pu.phase_pos[Oil];
            const int gaspos = pu.phase_pos[Gas];

            const ADB tmp_oil = (cmix_s[oilpos] - rv_perfcells * F_gas * cmix_s[gaspos]) / d;
            volumeRatio += tmp_oil / b_perfcells[oilpos];

            const ADB tmp_gas = (cmix_s[gaspos] - rs_perfcells * cmix_s[oilpos]) / d;
            volumeRatio += tmp_gas / b_perfcells[gaspos];
        }
        else {
            if ((*active_)[Oil]) {
                const int oilpos = pu.phase_pos[Oil];
                volumeRatio += cmix_s[oilpos] / b_perfcells[oilpos];
            }
            if ((*active_)[Gas]) {
                const int gaspos = pu.phase_pos[Gas];
                volumeRatio += cmix_s[gaspos] / b_perfcells[gaspos];
            }
        }


        // injecting connections total volumerates at standard conditions
        ADB cqt_is = cqt_i/volumeRatio;

        // connection phase volumerates at standard conditions
        cq_s.resize(np, ADB::null());
        for (int phase = 0; phase < np; ++phase) {
            cq_s[phase] = cq_ps[phase] + cmix_s[phase]*cqt_is;
        }

        // check for dead wells (used in the well controll equations)
        aliveWells = Vector::Constant(nw, 1.0);
        for (int w = 0; w < nw; ++w) {
            if (wbqt.value()[w] == 0) {
                aliveWells[w] = 0.0;
            }
        }
    }






    template <class SolutionState, class WellState>
    void
    StandardWellsSolvent::
    computeWellConnectionPressures(const SolutionState& state,
                                   const WellState& xw)
    {
        if( ! localWellsActive() ) return ;
        // 1. Compute properties required by computeConnectionPressureDelta().
        //    Note that some of the complexity of this part is due to the function
        //    taking std::vector<double> arguments, and not Eigen objects.
        std::vector<double> b_perf;
        std::vector<double> rsmax_perf;
        std::vector<double> rvmax_perf;
        std::vector<double> surf_dens_perf;
        computePropertiesForWellConnectionPressures(state, xw, b_perf, rsmax_perf, rvmax_perf, surf_dens_perf);

        const Vector pdepth = perf_cell_depth_;
        const int nperf = wells().well_connpos[wells().number_of_wells];
        const std::vector<double> depth_perf(pdepth.data(), pdepth.data() + nperf);

        computeWellConnectionDensitesPressures(xw, b_perf, rsmax_perf, rvmax_perf, surf_dens_perf, depth_perf, gravity_);

    }





    template <class ReservoirResidualQuant, class SolutionState>
    void
    StandardWellsSolvent::
    extractWellPerfProperties(const SolutionState& state,
                              const std::vector<ReservoirResidualQuant>& rq,
                              std::vector<ADB>& mob_perfcells,
                              std::vector<ADB>& b_perfcells) const
    {
        Base::extractWellPerfProperties(state, rq, mob_perfcells, b_perfcells);
        // handle the solvent related
        if (has_solvent_) {
            const Opm::PhaseUsage& pu = fluid_->phaseUsage();
            int gas_pos = pu.phase_pos[Gas];
            const std::vector<int>& well_cells = wellOps().well_cells;
            const int nperf = well_cells.size();
            // Gas and solvent is combinded and solved together
            // The input in the well equation is then the
            // total gas phase = hydro carbon gas + solvent gas

            // The total mobility is the sum of the solvent and gas mobiliy
            mob_perfcells[gas_pos] += subset(rq[solvent_pos_].mob, well_cells);

            // A weighted sum of the b-factors of gas and solvent are used.
            const int nc = rq[solvent_pos_].mob.size();

            const ADB zero = ADB::constant(Vector::Zero(nc));
            const ADB& ss = state.solvent_saturation;
            const ADB& sg = ((*active_)[ Gas ]
                             ? state.saturation[ pu.phase_pos[ Gas ] ]
                             : zero);

            Selector<double> zero_selector(ss.value() + sg.value(), Selector<double>::Zero);
            ADB F_solvent = subset(zero_selector.select(ss, ss / (ss + sg)),well_cells);
            Vector ones = Vector::Constant(nperf,1.0);

            b_perfcells[gas_pos] = (ones - F_solvent) * b_perfcells[gas_pos];
            b_perfcells[gas_pos] += (F_solvent * subset(rq[solvent_pos_].b, well_cells));
        }
    }


}
