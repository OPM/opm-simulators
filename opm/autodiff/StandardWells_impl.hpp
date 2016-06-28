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


#include <opm/autodiff/StandardWells.hpp>
#include <opm/autodiff/WellDensitySegmented.hpp>

#include <opm/autodiff/VFPInjProperties.hpp>
#include <opm/autodiff/VFPProdProperties.hpp>
#include <opm/autodiff/WellHelpers.hpp>




namespace Opm
{


    StandardWells::
    WellOps::WellOps(const Wells* wells)
      : w2p(),
        p2w(),
        well_cells()
    {
        if( wells )
        {
            w2p = Eigen::SparseMatrix<double>(wells->well_connpos[ wells->number_of_wells ], wells->number_of_wells);
            p2w = Eigen::SparseMatrix<double>(wells->number_of_wells, wells->well_connpos[ wells->number_of_wells ]);

            const int        nw   = wells->number_of_wells;
            const int* const wpos = wells->well_connpos;

            typedef Eigen::Triplet<double> Tri;

            std::vector<Tri> scatter, gather;
            scatter.reserve(wpos[nw]);
            gather .reserve(wpos[nw]);

            for (int w = 0, i = 0; w < nw; ++w) {
                for (; i < wpos[ w + 1 ]; ++i) {
                    scatter.push_back(Tri(i, w, 1.0));
                    gather .push_back(Tri(w, i, 1.0));
                }
            }

            w2p.setFromTriplets(scatter.begin(), scatter.end());
            p2w.setFromTriplets(gather .begin(), gather .end());

            well_cells.assign(wells->well_cells, wells->well_cells + wells->well_connpos[wells->number_of_wells]);
        }
    }





    StandardWells::StandardWells(const Wells* wells_arg)
      : wells_active_(wells_arg!=nullptr)
      , wells_(wells_arg)
      , wops_(wells_arg)
      , fluid_(nullptr)
      , active_(nullptr)
      , phase_condition_(nullptr)
      , vfp_properties_(nullptr)
      , well_perforation_densities_(Vector())
      , well_perforation_pressure_diffs_(Vector())
      , store_well_perforation_fluxes_(false)
    {
    }





    void
    StandardWells::init(const BlackoilPropsAdInterface* fluid_arg,
                        const std::vector<bool>* active_arg,
                        const std::vector<PhasePresence>* pc_arg,
                        const VFPProperties*  vfp_properties_arg,
                        const double gravity_arg,
                        const Vector& depth_arg)
    {
        fluid_ = fluid_arg;
        active_ = active_arg;
        phase_condition_ = pc_arg;
        vfp_properties_ = vfp_properties_arg;
        gravity_ = gravity_arg;
        perf_cell_depth_ = subset(depth_arg, wellOps().well_cells);;
    }





    const Wells& StandardWells::wells() const
    {
        assert(wells_ != 0);
        return *(wells_);
    }


    const Wells* StandardWells::wellsPointer() const
    {
        return wells_;
    }



    bool StandardWells::wellsActive() const
    {
        return wells_active_;
    }





    void StandardWells::setWellsActive(const bool wells_active)
    {
        wells_active_ = wells_active;
    }





    bool StandardWells::localWellsActive() const
    {
        return wells_ ? (wells_->number_of_wells > 0 ) : false;
    }





    int
    StandardWells::numWellVars() const
    {
        if ( !localWellsActive() )
        {
            return 0;
        }

        // For each well, we have a bhp variable, and one flux per phase.
        const int nw = wells().number_of_wells;
        return (numPhases() + 1) * nw;
    }





    const StandardWells::WellOps&
    StandardWells::wellOps() const
    {
        return wops_;
    }





    StandardWells::Vector& StandardWells::wellPerforationDensities()
    {
        return well_perforation_densities_;
    }





    const StandardWells::Vector&
    StandardWells::wellPerforationDensities() const
    {
        return well_perforation_densities_;
    }





    StandardWells::Vector&
    StandardWells::wellPerforationPressureDiffs()
    {
        return well_perforation_pressure_diffs_;
    }





    const StandardWells::Vector&
    StandardWells::wellPerforationPressureDiffs() const
    {
        return well_perforation_pressure_diffs_;
    }




    template<class SolutionState, class WellState>
    void
    StandardWells::
    computePropertiesForWellConnectionPressures(const SolutionState& state,
                                                const WellState& xw,
                                                std::vector<double>& b_perf,
                                                std::vector<double>& rsmax_perf,
                                                std::vector<double>& rvmax_perf,
                                                std::vector<double>& surf_dens_perf)
    {
        const int nperf = wells().well_connpos[wells().number_of_wells];
        const int nw = wells().number_of_wells;

        // Compute the average pressure in each well block
        const Vector perf_press = Eigen::Map<const Vector>(xw.perfPress().data(), nperf);
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
        // const std::vector<PhasePresence>& pc = phaseCondition();
        for (int perf = 0; perf < nperf; ++perf) {
            perf_cond[perf] = (*phase_condition_)[well_cells[perf]];
        }
        const PhaseUsage& pu = fluid_->phaseUsage();
        DataBlock b(nperf, pu.num_phases);
        if (pu.phase_used[BlackoilPhases::Aqua]) {
            const Vector bw = fluid_->bWat(avg_press_ad, perf_temp, well_cells).value();
            b.col(pu.phase_pos[BlackoilPhases::Aqua]) = bw;
        }
        assert((*active_)[Oil]);
        const Vector perf_so =  subset(state.saturation[pu.phase_pos[Oil]].value(), well_cells);
        if (pu.phase_used[BlackoilPhases::Liquid]) {
            const ADB perf_rs = (state.rs.size() > 0) ? subset(state.rs, well_cells) : ADB::null();
            const Vector bo = fluid_->bOil(avg_press_ad, perf_temp, perf_rs, perf_cond, well_cells).value();
            b.col(pu.phase_pos[BlackoilPhases::Liquid]) = bo;
        }
        if (pu.phase_used[BlackoilPhases::Vapour]) {
            const ADB perf_rv = (state.rv.size() > 0) ? subset(state.rv, well_cells) : ADB::null();
            const Vector bg = fluid_->bGas(avg_press_ad, perf_temp, perf_rv, perf_cond, well_cells).value();
            b.col(pu.phase_pos[BlackoilPhases::Vapour]) = bg;
        }
        if (pu.phase_used[BlackoilPhases::Liquid] && pu.phase_used[BlackoilPhases::Vapour]) {
            const Vector rssat = fluid_->rsSat(ADB::constant(avg_press), ADB::constant(perf_so), well_cells).value();
            rsmax_perf.assign(rssat.data(), rssat.data() + nperf);

            const Vector rvsat = fluid_->rvSat(ADB::constant(avg_press), ADB::constant(perf_so), well_cells).value();
            rvmax_perf.assign(rvsat.data(), rvsat.data() + nperf);
        }

        // b is row major, so can just copy data.
        b_perf.assign(b.data(), b.data() + nperf * pu.num_phases);

        // Surface density.
        // The compute density segment wants the surface densities as
        // an np * number of wells cells array
        Vector rho = superset(fluid_->surfaceDensity(0 , well_cells), Span(nperf, pu.num_phases, 0), nperf*pu.num_phases);
        for (int phase = 1; phase < pu.num_phases; ++phase) {
            rho += superset(fluid_->surfaceDensity(phase , well_cells), Span(nperf, pu.num_phases, phase), nperf*pu.num_phases);
        }
        surf_dens_perf.assign(rho.data(), rho.data() + nperf * pu.num_phases);

    }





    template <class WellState>
    void
    StandardWells::
    computeWellConnectionDensitesPressures(const WellState& xw,
                                           const std::vector<double>& b_perf,
                                           const std::vector<double>& rsmax_perf,
                                           const std::vector<double>& rvmax_perf,
                                           const std::vector<double>& surf_dens_perf,
                                           const std::vector<double>& depth_perf,
                                           const double grav)
    {
        // Compute densities
        std::vector<double> cd =
                WellDensitySegmented::computeConnectionDensities(
                        wells(), xw, fluid_->phaseUsage(),
                        b_perf, rsmax_perf, rvmax_perf, surf_dens_perf);

        const int nperf = wells().well_connpos[wells().number_of_wells];

        // Compute pressure deltas
        std::vector<double> cdp =
                WellDensitySegmented::computeConnectionPressureDelta(
                        wells(), depth_perf, cd, grav);

        // Store the results
        well_perforation_densities_ = Eigen::Map<const Vector>(cd.data(), nperf);
        well_perforation_pressure_diffs_ = Eigen::Map<const Vector>(cdp.data(), nperf);
    }





    template <class SolutionState, class WellState>
    void
    StandardWells::
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

        const Vector& pdepth = perf_cell_depth_;
        const int nperf = wells().well_connpos[wells().number_of_wells];
        const std::vector<double> depth_perf(pdepth.data(), pdepth.data() + nperf);

        computeWellConnectionDensitesPressures(xw, b_perf, rsmax_perf, rvmax_perf, surf_dens_perf, depth_perf, gravity_);

    }





    template <class ReservoirResidualQuant, class SolutionState>
    void
    StandardWells::
    extractWellPerfProperties(const SolutionState& /* state */,
                              const std::vector<ReservoirResidualQuant>& rq,
                              std::vector<ADB>& mob_perfcells,
                              std::vector<ADB>& b_perfcells) const
    {
        // If we have wells, extract the mobilities and b-factors for
        // the well-perforated cells.
        if ( !localWellsActive() ) {
            mob_perfcells.clear();
            b_perfcells.clear();
            return;
        } else {
            const std::vector<int>& well_cells = wellOps().well_cells;
            const int num_phases = wells().number_of_phases;
            mob_perfcells.resize(num_phases, ADB::null());
            b_perfcells.resize(num_phases, ADB::null());
            for (int phase = 0; phase < num_phases; ++phase) {
                mob_perfcells[phase] = subset(rq[phase].mob, well_cells);
                b_perfcells[phase] = subset(rq[phase].b, well_cells);
            }
        }
    }






    template <class SolutionState>
    void
    StandardWells::
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
        const Opm::PhaseUsage& pu = fluid_->phaseUsage();
        if ((*active_)[Oil] && (*active_)[Gas]) {
            const int oilpos = pu.phase_pos[Oil];
            const int gaspos = pu.phase_pos[Gas];
            const ADB cq_psOil = cq_ps[oilpos];
            const ADB cq_psGas = cq_ps[gaspos];
            const ADB& rv_perfcells = subset(state.rv, well_cells);
            const ADB& rs_perfcells = subset(state.rs, well_cells);
            cq_ps[gaspos] += rs_perfcells * cq_psOil;
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

            const ADB tmp_oil = (cmix_s[oilpos] - rv_perfcells * cmix_s[gaspos]) / d;
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
    StandardWells::
    updatePerfPhaseRatesAndPressures(const std::vector<ADB>& cq_s,
                                     const SolutionState& state,
                                     WellState& xw) const
    {
        if ( !localWellsActive() )
        {
            // If there are no wells in the subdomain of the proces then
            // cq_s has zero size and will cause a segmentation fault below.
            return;
        }

        // Update the perforation phase rates (used to calculate the pressure drop in the wellbore).
        const int np = wells().number_of_phases;
        const int nw = wells().number_of_wells;
        const int nperf = wells().well_connpos[nw];
        Vector cq = superset(cq_s[0].value(), Span(nperf, np, 0), nperf*np);
        for (int phase = 1; phase < np; ++phase) {
            cq += superset(cq_s[phase].value(), Span(nperf, np, phase), nperf*np);
        }
        xw.perfPhaseRates().assign(cq.data(), cq.data() + nperf*np);

        // Update the perforation pressures.
        const Vector& cdp = wellPerforationPressureDiffs();
        const Vector perfpressure = (wellOps().w2p * state.bhp.value().matrix()).array() + cdp;
        xw.perfPress().assign(perfpressure.data(), perfpressure.data() + nperf);
    }





    template <class WellState>
    void
    StandardWells::
    updateWellState(const Vector& dwells,
                    const double dpmaxrel,
                    WellState& well_state)
    {
        if( localWellsActive() )
        {
            const int np = wells().number_of_phases;
            const int nw = wells().number_of_wells;

            // Extract parts of dwells corresponding to each part.
            int varstart = 0;
            const Vector dqs = subset(dwells, Span(np*nw, 1, varstart));
            varstart += dqs.size();
            const Vector dbhp = subset(dwells, Span(nw, 1, varstart));
            varstart += dbhp.size();
            assert(varstart == dwells.size());


            // Qs update.
            // Since we need to update the wellrates, that are ordered by wells,
            // from dqs which are ordered by phase, the simplest is to compute
            // dwr, which is the data from dqs but ordered by wells.
            const DataBlock wwr = Eigen::Map<const DataBlock>(dqs.data(), np, nw).transpose();
            const Vector dwr = Eigen::Map<const Vector>(wwr.data(), nw*np);
            const Vector wr_old = Eigen::Map<const Vector>(&well_state.wellRates()[0], nw*np);
            const Vector wr = wr_old - dwr;
            std::copy(&wr[0], &wr[0] + wr.size(), well_state.wellRates().begin());

            // Bhp update.
            const Vector bhp_old = Eigen::Map<const Vector>(&well_state.bhp()[0], nw, 1);
            const Vector dbhp_limited = sign(dbhp) * dbhp.abs().min(bhp_old.abs()*dpmaxrel);
            const Vector bhp = bhp_old - dbhp_limited;
            std::copy(&bhp[0], &bhp[0] + bhp.size(), well_state.bhp().begin());


            const Opm::PhaseUsage& pu = fluid_->phaseUsage();
            //Loop over all wells
#pragma omp parallel for schedule(static)
            for (int w = 0; w < nw; ++w) {
                const WellControls* wc = wells().ctrls[w];
                const int nwc = well_controls_get_num(wc);
                //Loop over all controls until we find a THP control
                //that specifies what we need...
                //Will only update THP for wells with THP control
                for (int ctrl_index=0; ctrl_index < nwc; ++ctrl_index) {
                    if (well_controls_iget_type(wc, ctrl_index) == THP) {
                        double aqua = 0.0;
                        double liquid = 0.0;
                        double vapour = 0.0;

                        if ((*active_)[ Water ]) {
                            aqua = wr[w*np + pu.phase_pos[ Water ] ];
                        }
                        if ((*active_)[ Oil ]) {
                            liquid = wr[w*np + pu.phase_pos[ Oil ] ];
                        }
                        if ((*active_)[ Gas ]) {
                            vapour = wr[w*np + pu.phase_pos[ Gas ] ];
                        }

                        double alq = well_controls_iget_alq(wc, ctrl_index);
                        int table_id = well_controls_iget_vfp(wc, ctrl_index);

                        const WellType& well_type = wells().type[w];
                        if (well_type == INJECTOR) {
                            double dp = wellhelpers::computeHydrostaticCorrection(
                                    wells(), w, vfp_properties_->getInj()->getTable(table_id)->getDatumDepth(),
                                    wellPerforationDensities(), gravity_);

                            well_state.thp()[w] = vfp_properties_->getInj()->thp(table_id, aqua, liquid, vapour, bhp[w] + dp);
                        }
                        else if (well_type == PRODUCER) {
                            double dp = wellhelpers::computeHydrostaticCorrection(
                                    wells(), w, vfp_properties_->getProd()->getTable(table_id)->getDatumDepth(),
                                    wellPerforationDensities(), gravity_);

                            well_state.thp()[w] = vfp_properties_->getProd()->thp(table_id, aqua, liquid, vapour, bhp[w] + dp, alq);
                        }
                        else {
                            OPM_THROW(std::logic_error, "Expected INJECTOR or PRODUCER well");
                        }

                        //Assume only one THP control specified for each well
                        break;
                    }
                }
            }
        }
    }





    template <class WellState>
    void
    StandardWells::
    updateWellControls(WellState& xw) const
    {
        if( !localWellsActive() ) return ;

        std::string modestring[4] = { "BHP", "THP", "RESERVOIR_RATE", "SURFACE_RATE" };
        // Find, for each well, if any constraints are broken. If so,
        // switch control to first broken constraint.
        const int np = wells().number_of_phases;
        const int nw = wells().number_of_wells;
#pragma omp parallel for schedule(dynamic)
        for (int w = 0; w < nw; ++w) {
            const WellControls* wc = wells().ctrls[w];
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
                if (wellhelpers::constraintBroken(
                        xw.bhp(), xw.thp(), xw.wellRates(),
                        w, np, wells().type[w], wc, ctrl_index)) {
                    // ctrl_index will be the index of the broken constraint after the loop.
                    break;
                }
            }
            if (ctrl_index != nwc) {
                // Constraint number ctrl_index was broken, switch to it.
                // We disregard terminal_ouput here as with it only messages
                // for wells on one process will be printed.
                std::ostringstream ss;
                ss << "Switching control mode for well " << wells().name[w]
                   << " from " << modestring[well_controls_iget_type(wc, current)]
                   << " to " << modestring[well_controls_iget_type(wc, ctrl_index)] << std::endl;
                OpmLog::info(ss.str());
                xw.currentControls()[w] = ctrl_index;
                current = xw.currentControls()[w];
            }

            // Updating well state and primary variables.
            // Target values are used as initial conditions for BHP, THP, and SURFACE_RATE
            const double target = well_controls_iget_target(wc, current);
            const double* distr = well_controls_iget_distr(wc, current);
            switch (well_controls_iget_type(wc, current)) {
            case BHP:
                xw.bhp()[w] = target;
                break;

            case THP: {
                double aqua = 0.0;
                double liquid = 0.0;
                double vapour = 0.0;

                const Opm::PhaseUsage& pu = fluid_->phaseUsage();

                if ((*active_)[ Water ]) {
                    aqua = xw.wellRates()[w*np + pu.phase_pos[ Water ] ];
                }
                if ((*active_)[ Oil ]) {
                    liquid = xw.wellRates()[w*np + pu.phase_pos[ Oil ] ];
                }
                if ((*active_)[ Gas ]) {
                    vapour = xw.wellRates()[w*np + pu.phase_pos[ Gas ] ];
                }

                const int vfp        = well_controls_iget_vfp(wc, current);
                const double& thp    = well_controls_iget_target(wc, current);
                const double& alq    = well_controls_iget_alq(wc, current);

                //Set *BHP* target by calculating bhp from THP
                const WellType& well_type = wells().type[w];

                if (well_type == INJECTOR) {
                    double dp = wellhelpers::computeHydrostaticCorrection(
                            wells(), w, vfp_properties_->getInj()->getTable(vfp)->getDatumDepth(),
                            wellPerforationDensities(), gravity_);

                    xw.bhp()[w] = vfp_properties_->getInj()->bhp(vfp, aqua, liquid, vapour, thp) - dp;
                }
                else if (well_type == PRODUCER) {
                    double dp = wellhelpers::computeHydrostaticCorrection(
                            wells(), w, vfp_properties_->getProd()->getTable(vfp)->getDatumDepth(),
                            wellPerforationDensities(), gravity_);

                    xw.bhp()[w] = vfp_properties_->getProd()->bhp(vfp, aqua, liquid, vapour, thp, alq) - dp;
                }
                else {
                    OPM_THROW(std::logic_error, "Expected PRODUCER or INJECTOR type of well");
                }
                break;
            }

            case RESERVOIR_RATE:
                // No direct change to any observable quantity at
                // surface condition.  In this case, use existing
                // flow rates as initial conditions as reservoir
                // rate acts only in aggregate.
                break;

            case SURFACE_RATE:
                // assign target value as initial guess for injectors and
                // single phase producers (orat, grat, wrat)
                const WellType& well_type = wells().type[w];
                if (well_type == INJECTOR) {
                    for (int phase = 0; phase < np; ++phase) {
                        const double& compi = wells().comp_frac[np * w + phase];
                        if (compi > 0.0) {
                            xw.wellRates()[np*w + phase] = target * compi;
                        }
                    }
                } else if (well_type == PRODUCER) {

                    // only set target as initial rates for single phase
                    // producers. (orat, grat and wrat, and not lrat)
                    // lrat will result in numPhasesWithTargetsUnderThisControl == 2
                    int numPhasesWithTargetsUnderThisControl = 0;
                    for (int phase = 0; phase < np; ++phase) {
                        if (distr[phase] > 0.0) {
                            numPhasesWithTargetsUnderThisControl += 1;
                        }
                    }
                    for (int phase = 0; phase < np; ++phase) {
                        if (distr[phase] > 0.0 && numPhasesWithTargetsUnderThisControl < 2 ) {
                            xw.wellRates()[np*w + phase] = target * distr[phase];
                        }
                    }
                } else {
                    OPM_THROW(std::logic_error, "Expected PRODUCER or INJECTOR type of well");
                }


                break;
            }
        }

    }





    template <class SolutionState>
    void
    StandardWells::
    addWellFluxEq(const std::vector<ADB>& cq_s,
                  const SolutionState& state,
                  LinearisedBlackoilResidual& residual)
    {
        if( !localWellsActive() )
        {
            // If there are no wells in the subdomain of the proces then
            // cq_s has zero size and will cause a segmentation fault below.
            return;
        }

        const int np = wells().number_of_phases;
        const int nw = wells().number_of_wells;
        ADB qs = state.qs;
        for (int phase = 0; phase < np; ++phase) {
            qs -= superset(wellOps().p2w * cq_s[phase], Span(nw, 1, phase*nw), nw*np);

        }

        residual.well_flux_eq = qs;
    }





    template <class SolutionState, class WellState>
    void
    StandardWells::addWellControlEq(const SolutionState& state,
                                    const WellState& xw,
                                    const Vector& aliveWells,
                                    LinearisedBlackoilResidual& residual)
    {
        if( ! localWellsActive() ) return;

        const int np = wells().number_of_phases;
        const int nw = wells().number_of_wells;

        ADB aqua   = ADB::constant(Vector::Zero(nw));
        ADB liquid = ADB::constant(Vector::Zero(nw));
        ADB vapour = ADB::constant(Vector::Zero(nw));

        if ((*active_)[Water]) {
            aqua += subset(state.qs, Span(nw, 1, BlackoilPhases::Aqua*nw));
        }
        if ((*active_)[Oil]) {
            liquid += subset(state.qs, Span(nw, 1, BlackoilPhases::Liquid*nw));
        }
        if ((*active_)[Gas]) {
            vapour += subset(state.qs, Span(nw, 1, BlackoilPhases::Vapour*nw));
        }

        //THP calculation variables
        std::vector<int> inj_table_id(nw, -1);
        std::vector<int> prod_table_id(nw, -1);
        Vector thp_inj_target_v = Vector::Zero(nw);
        Vector thp_prod_target_v = Vector::Zero(nw);
        Vector alq_v = Vector::Zero(nw);

        //Hydrostatic correction variables
        Vector rho_v = Vector::Zero(nw);
        Vector vfp_ref_depth_v = Vector::Zero(nw);

        //Target vars
        Vector bhp_targets  = Vector::Zero(nw);
        Vector rate_targets = Vector::Zero(nw);
        Eigen::SparseMatrix<double> rate_distr(nw, np*nw);

        //Selection variables
        std::vector<int> bhp_elems;
        std::vector<int> thp_inj_elems;
        std::vector<int> thp_prod_elems;
        std::vector<int> rate_elems;

        //Run through all wells to calculate BHP/RATE targets
        //and gather info about current control
        for (int w = 0; w < nw; ++w) {
            auto wc = wells().ctrls[w];

            // The current control in the well state overrides
            // the current control set in the Wells struct, which
            // is instead treated as a default.
            const int current = xw.currentControls()[w];

            switch (well_controls_iget_type(wc, current)) {
            case BHP:
            {
                bhp_elems.push_back(w);
                bhp_targets(w)  = well_controls_iget_target(wc, current);
                rate_targets(w) = -1e100;
            }
            break;

            case THP:
            {
                const int perf = wells().well_connpos[w];
                rho_v[w] = wellPerforationDensities()[perf];

                const int table_id = well_controls_iget_vfp(wc, current);
                const double target = well_controls_iget_target(wc, current);

                const WellType& well_type = wells().type[w];
                if (well_type == INJECTOR) {
                    inj_table_id[w]  = table_id;
                    thp_inj_target_v[w] = target;
                    alq_v[w]     = -1e100;

                    vfp_ref_depth_v[w] = vfp_properties_->getInj()->getTable(table_id)->getDatumDepth();

                    thp_inj_elems.push_back(w);
                }
                else if (well_type == PRODUCER) {
                    prod_table_id[w]  = table_id;
                    thp_prod_target_v[w] = target;
                    alq_v[w]      = well_controls_iget_alq(wc, current);

                    vfp_ref_depth_v[w] =  vfp_properties_->getProd()->getTable(table_id)->getDatumDepth();

                    thp_prod_elems.push_back(w);
                }
                else {
                    OPM_THROW(std::logic_error, "Expected INJECTOR or PRODUCER type well");
                }
                bhp_targets(w)  = -1e100;
                rate_targets(w) = -1e100;
            }
            break;

            case RESERVOIR_RATE: // Intentional fall-through
            case SURFACE_RATE:
            {
                rate_elems.push_back(w);
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
        }

        //Calculate BHP target from THP
        const ADB thp_inj_target = ADB::constant(thp_inj_target_v);
        const ADB thp_prod_target = ADB::constant(thp_prod_target_v);
        const ADB alq = ADB::constant(alq_v);
        const ADB bhp_from_thp_inj = vfp_properties_->getInj()->bhp(inj_table_id, aqua, liquid, vapour, thp_inj_target);
        const ADB bhp_from_thp_prod = vfp_properties_->getProd()->bhp(prod_table_id, aqua, liquid, vapour, thp_prod_target, alq);

        //Perform hydrostatic correction to computed targets
        const Vector dp_v = wellhelpers::computeHydrostaticCorrection(wells(), vfp_ref_depth_v, wellPerforationDensities(), gravity_);
        const ADB dp = ADB::constant(dp_v);
        const ADB dp_inj = superset(subset(dp, thp_inj_elems), thp_inj_elems, nw);
        const ADB dp_prod = superset(subset(dp, thp_prod_elems), thp_prod_elems, nw);

        //Calculate residuals
        const ADB thp_inj_residual = state.bhp - bhp_from_thp_inj + dp_inj;
        const ADB thp_prod_residual = state.bhp - bhp_from_thp_prod + dp_prod;
        const ADB bhp_residual = state.bhp - bhp_targets;
        const ADB rate_residual = rate_distr * state.qs - rate_targets;

        //Select the right residual for each well
        residual.well_eq = superset(subset(bhp_residual, bhp_elems), bhp_elems, nw) +
                superset(subset(thp_inj_residual, thp_inj_elems), thp_inj_elems, nw) +
                superset(subset(thp_prod_residual, thp_prod_elems), thp_prod_elems, nw) +
                superset(subset(rate_residual, rate_elems), rate_elems, nw);

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
        residual.well_eq = alive_selector.select(residual.well_eq, rate_summer * state.qs);
        // OPM_AD_DUMP(residual_.well_eq);
    }





    template <class SolutionState, class WellState>
    void
    StandardWells::computeWellPotentials(const std::vector<ADB>& mob_perfcells,
                                         const std::vector<ADB>& b_perfcells,
                                         SolutionState& state0,
                                         WellState& well_state)
    {
        const int nw = wells().number_of_wells;
        const int np = wells().number_of_phases;
        const Opm::PhaseUsage& pu = fluid_->phaseUsage();

        Vector bhps = Vector::Zero(nw);
        for (int w = 0; w < nw; ++w) {
            const WellControls* ctrl = wells().ctrls[w];
            const int nwc = well_controls_get_num(ctrl);
            //Loop over all controls until we find a BHP control
            //or a THP control that specifies what we need.
            //Pick the value that gives the most restrictive flow
            for (int ctrl_index=0; ctrl_index < nwc; ++ctrl_index) {

                if (well_controls_iget_type(ctrl, ctrl_index) == BHP) {
                    bhps[w] = well_controls_iget_target(ctrl, ctrl_index);
                }

                if (well_controls_iget_type(ctrl, ctrl_index) == THP) {
                    double aqua = 0.0;
                    double liquid = 0.0;
                    double vapour = 0.0;

                    if ((*active_)[ Water ]) {
                        aqua = well_state.wellRates()[w*np + pu.phase_pos[ Water ] ];
                    }
                    if ((*active_)[ Oil ]) {
                        liquid = well_state.wellRates()[w*np + pu.phase_pos[ Oil ] ];
                    }
                    if ((*active_)[ Gas ]) {
                        vapour = well_state.wellRates()[w*np + pu.phase_pos[ Gas ] ];
                    }

                    const int vfp        = well_controls_iget_vfp(ctrl, ctrl_index);
                    const double& thp    = well_controls_iget_target(ctrl, ctrl_index);
                    const double& alq    = well_controls_iget_alq(ctrl, ctrl_index);

                    //Set *BHP* target by calculating bhp from THP
                    const WellType& well_type = wells().type[w];

                    if (well_type == INJECTOR) {
                        double dp = wellhelpers::computeHydrostaticCorrection(
                                    wells(), w, vfp_properties_->getInj()->getTable(vfp)->getDatumDepth(),
                                    wellPerforationDensities(), gravity_);
                        const double bhp = vfp_properties_->getInj()->bhp(vfp, aqua, liquid, vapour, thp) - dp;
                        // apply the strictest of the bhp controlls i.e. smallest bhp for injectors
                        if ( bhp < bhps[w]) {
                            bhps[w] = bhp;
                        }
                    }
                    else if (well_type == PRODUCER) {
                        double dp = wellhelpers::computeHydrostaticCorrection(
                                    wells(), w, vfp_properties_->getProd()->getTable(vfp)->getDatumDepth(),
                                    wellPerforationDensities(), gravity_);

                        const double bhp = vfp_properties_->getProd()->bhp(vfp, aqua, liquid, vapour, thp, alq) - dp;
                        // apply the strictest of the bhp controlls i.e. largest bhp for producers
                        if ( bhp > bhps[w]) {
                            bhps[w] = bhp;
                        }
                    }
                    else {
                        OPM_THROW(std::logic_error, "Expected PRODUCER or INJECTOR type of well");
                    }
                }
            }

        }

        // use bhp limit from control
        state0.bhp = ADB::constant(bhps);

        // compute well potentials
        Vector aliveWells;
        std::vector<ADB> well_potentials;
        computeWellFlux(state0, mob_perfcells,  b_perfcells, aliveWells, well_potentials);

        // store well potentials in the well state
        // transform to a single vector instead of separate vectors pr phase
        const int nperf = wells().well_connpos[nw];
        Vector cq = superset(well_potentials[0].value(), Span(nperf, np, 0), nperf*np);
        for (int phase = 1; phase < np; ++phase) {
            cq += superset(well_potentials[phase].value(), Span(nperf, np, phase), nperf*np);
        }
        well_state.wellPotentials().assign(cq.data(), cq.data() + nperf*np);
    }





    void
    StandardWells::variableStateWellIndices(std::vector<int>& indices,
                                            int& next) const
    {
        indices[Qs] = next++;
        indices[Bhp] = next++;
    }





    template <class SolutionState>
    void
    StandardWells::
    variableStateExtractWellsVars(const std::vector<int>& indices,
                                  std::vector<ADB>& vars,
                                  SolutionState& state) const
    {
        // Qs.
        state.qs = std::move(vars[indices[Qs]]);

        // Bhp.
        state.bhp = std::move(vars[indices[Bhp]]);
    }





    std::vector<int>
    StandardWells::variableWellStateIndices() const
    {
        // Black oil model standard is 5 equation.
        // For the pure well solve, only the well equations are picked.
        std::vector<int> indices(5, -1);
        int next = 0;

        variableStateWellIndices(indices, next);

        assert(next == 2);
        return indices;
    }





    template <class WellState>
    void
    StandardWells::variableWellStateInitials(const WellState& xw,
                                             std::vector<Vector>& vars0) const
    {
        // Initial well rates.
        if ( localWellsActive() )
        {
            // Need to reshuffle well rates, from phase running fastest
            // to wells running fastest.
            const int nw = wells().number_of_wells;
            const int np = wells().number_of_phases;

            // The transpose() below switches the ordering.
            const DataBlock wrates = Eigen::Map<const DataBlock>(& xw.wellRates()[0], nw, np).transpose();
            const Vector qs = Eigen::Map<const V>(wrates.data(), nw*np);
            vars0.push_back(qs);

            // Initial well bottom-hole pressure.
            assert (not xw.bhp().empty());
            const Vector bhp = Eigen::Map<const V>(& xw.bhp()[0], xw.bhp().size());
            vars0.push_back(bhp);
        }
        else
        {
            // push null states for qs and bhp
            vars0.push_back(V());
            vars0.push_back(V());
        }
    }





    void
    StandardWells::setStoreWellPerforationFluxesFlag(const bool store_fluxes)
    {
        store_well_perforation_fluxes_ = store_fluxes;
    }





    const StandardWells::Vector&
    StandardWells::getStoredWellPerforationFluxes() const
    {
        assert(store_well_perforation_fluxes_);
        return well_perforation_fluxes_;
    }






   template<class WellState>
   void
   StandardWells::
   updateListEconLimited(ScheduleConstPtr schedule,
                         const int current_step,
                         const Wells* wells_struct,
                         const WellState& well_state,
                         DynamicListEconLimited& list_econ_limited) const
   {
       const int nw = wells_struct->number_of_wells;
       const int np = wells_struct->number_of_phases;

       for (int w = 0; w < nw; ++w) {
           // flag to check if the mim oil/gas rate limit is violated
           bool rate_limit_violated = false;
           const std::string& well_name = wells_struct->name[w];
           const Well* well_ecl = schedule->getWell(well_name);
           const WellEconProductionLimits& econ_production_limits = well_ecl->getEconProductionLimits(current_step);

           // if no limit is effective here, then continue to the next well
           if ( !econ_production_limits.onAnyEffectiveLimit() ) {
               continue;
           }
           // for the moment, we only handle rate limits, not handling potential limits
           // the potential limits should not be difficult to add
           const WellEcon::QuantityLimitEnum& quantity_limit = econ_production_limits.quantityLimit();
           if (quantity_limit == WellEcon::POTN) {
               OPM_THROW(std::logic_error, "Only RATE limit is supported for the moment");
           }

           const WellMapType& well_map = well_state.wellMap();
           typename WellMapType::const_iterator i = well_map.find(well_name);

           assert(i != well_map.end()); // should always be found?

           const int well_number = (i->second)[0];

           const Opm::PhaseUsage& pu = fluid_->phaseUsage();

           // oil rate limit
           // TODO: SHOULD we give a message when closing a well due to economic limit reason?
           if (econ_production_limits.onMinOilRate()) {
               assert((*active_)[Oil]);
               const double oil_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Oil ] ];
               const double min_oil_rate = econ_production_limits.minOilRate();
               if (std::abs(oil_rate) < min_oil_rate) {
                   rate_limit_violated = true;
               }
           }
           // gas rate limit
           if (econ_production_limits.onMinGasRate()) {
               assert((*active_)[Gas]);
               const double gas_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Gas ] ];
               const double min_gas_rate = econ_production_limits.minGasRate();
               if (std::abs(gas_rate) < min_gas_rate) {
                   rate_limit_violated = true;
               }
           }

           if (rate_limit_violated) {
               if (econ_production_limits.endRun()) {
                   std::cerr << "WARNING: ending run after well closed due to economic limits is not supported yet"
                             << std::endl
                             << "the program will keep running after " << well_name << " is closed " << std::endl;
               }

               if (econ_production_limits.validFollowonWell()) {
                   std::cerr << "WARNING: opening following on well after well closed is not supported yet" << std::endl;
               }

               list_econ_limited.addShuttedWell(well_name);
               // the well is closed, not need to check other limits
               continue;
           }

           // checking for other limits, mostly all kinds of ratio.
           // TODO: not sure when more than one ratio limit is violated, what is the definition of the
           // worst-offending connection
           // Should not be more than one connection be closed each time
           // Maybe we can compare the extent of the different violations to give a worst one.
           bool water_cut_limit_violated = false;
           if (econ_production_limits.onMaxWaterCut()) {
               assert((*active_)[Oil]);
               assert((*active_)[Water]);
               const double oil_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Oil ] ];
               const double water_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Water ] ];
               const double liquid_rate = oil_rate + water_rate;
               if (std::abs(liquid_rate) == 0.) {
                   continue;
               }
               const double water_cut = water_rate / liquid_rate;
               const double max_water_cut = econ_production_limits.maxWaterCut();
               if (water_cut > max_water_cut) {
                   water_cut_limit_violated = true;
               }
           }

           if (econ_production_limits.workover() != WellEcon::CON) {
               std::cerr << "WARNING: only CON workover over exceeding limit is supported for the moment." << std::endl;
           }

           const int perf_start = (i->second)[1];
           const int perf_number = (i->second)[2];

           std::vector<double> oil_perf_rate(perf_number);
           std::vector<double> gas_perf_rate(perf_number);
           std::vector<double> water_perf_rate(perf_number);

           for (int perf = 0; perf < perf_number; ++perf) {
               const int i_perf = perf_start + perf;
               oil_perf_rate[perf] = well_state.perfPhaseRates()[i_perf * np + pu.phase_pos[ Oil ] ];
               water_perf_rate[perf] = well_state.perfPhaseRates()[i_perf * np + pu.phase_pos[ Water ] ];
               gas_perf_rate[perf] = well_state.perfPhaseRates()[i_perf * np + pu.phase_pos[ Gas ] ];
           }

           std::vector<double> water_cut_perf(perf_number);
           double max_water_cut_perf = 0.0;
           int worst_offending_perf = -1;
           for (int perf = 0; perf < perf_number; ++perf) {
               const double liquid_rate = water_perf_rate[perf] + oil_perf_rate[perf];
               if (std::abs(liquid_rate) == 0.0) {
                   water_cut_perf[perf] = 0.0;
               } else {
                   water_cut_perf[perf] = water_perf_rate[perf] / liquid_rate;
               }

               if (water_cut_perf[perf] > max_water_cut_perf) {
                   worst_offending_perf = perf;
                   max_water_cut_perf = water_cut_perf[perf];
               }
           }

           assert((worst_offending_perf >= 0) && (worst_offending_perf <= perf_number));

           const int cell_worst_offending_perf = wells_struct->well_cells[perf_start + worst_offending_perf];

           list_econ_limited.addClosedConnectionsForWell(well_name, cell_worst_offending_perf);

       }
   }





    template <class WellState>
    bool
    StandardWells::
    checkRateEconLimits(const WellEconProductionLimits& econ_production_limits,
                        const WellState& well_state,
                        const int well_number) const
    {
        const Opm::PhaseUsage& pu = fluid_->phaseUsage();
        const int np = well_state.numPhases();

        if (econ_production_limits.onMinOilRate()) {
            assert((*active_)[Oil]);
            const double oil_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Oil ] ];
            const double min_oil_rate = econ_production_limits.minOilRate();
            if (std::abs(oil_rate) < min_oil_rate) {
                return true;
            }
        }

        if (econ_production_limits.onMinGasRate() ) {
            assert((*active_)[Gas]);
            const double gas_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Gas ] ];
            const double min_gas_rate = econ_production_limits.minGasRate();
            if (std::abs(gas_rate) < min_gas_rate) {
                return true;
            }
        }

        if (econ_production_limits.onMinLiquidRate() ) {
            assert((*active_)[Oil]);
            assert((*active_)[Water]);
            const double oil_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Oil ] ];
            const double water_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Water ] ];
            const double liquid_rate = oil_rate + water_rate;
            const double min_liquid_rate = econ_production_limits.minLiquidRate();
            if (std::abs(liquid_rate) < min_liquid_rate) {
                return true;
            }
        }

        if (econ_production_limits.onMinReservoirFluidRate()) {
            std::cerr << "WARNING: Minimum reservoir fluid production rate limit is not supported yet" << std::endl;
        }

        return false;
    }







    template <class WellState>
    bool
    StandardWells::
    checkMaxWaterCutLimit(const WellEconProductionLimits& econ_production_limits,
                          const WellState& well_state,
                          const typename WellMapType::const_iterator& i_well,
                          int& worst_offending_connection,
                          double& violation_extent,
                          bool& last_connection) const
    {
        bool water_cut_limit_violated = false;
        worst_offending_connection = -1;
        violation_extent = -1.0;
        last_connection = false;

        const int np = well_state.numPhases();
        const Opm::PhaseUsage& pu = fluid_->phaseUsage();
        const int well_number = (i_well->second)[0];

        assert((*active_)[Oil]);
        assert((*active_)[Water]);

        const double oil_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Oil ] ];
        const double water_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Water ] ];
        const double liquid_rate = oil_rate + water_rate;
        double water_cut;
        if (std::abs(liquid_rate) != 0.) {
            water_cut = water_rate / liquid_rate;
        } else {
            water_cut = 0.0;
        }

        const double max_water_cut_limit = econ_production_limits.maxWaterCut();
        if (water_cut > max_water_cut_limit) {
            water_cut_limit_violated = true;
        }

        if (water_cut_limit_violated) {
            // need to handle the worst_offending_connection
            const int perf_start = (i_well->second)[1];
            const int perf_number = (i_well->second)[2];

            std::vector<double> water_cut_perf(perf_number);
            for (int perf = 0; perf < perf_number; ++perf) {
                const int i_perf = perf_start + perf;
                const double oil_perf_rate = well_state.perfPhaseRates()[i_perf * np + pu.phase_pos[ Oil ] ];
                const double water_perf_rate = well_state.perfPhaseRates()[i_perf * np + pu.phase_pos[ Water ] ];
                const double liquid_perf_rate = oil_perf_rate + water_perf_rate;
                if (std::abs(liquid_perf_rate) != 0.) {
                    water_cut_perf[perf] = water_perf_rate / liquid_perf_rate;
                } else {
                    water_cut_perf[perf] = 0.;
                }
            }

            if (perf_number == 1) {
                last_connection = true;
                worst_offending_connection = 0;
                violation_extent = water_cut_perf[0] / max_water_cut_limit;
                return water_cut_limit_violated;
            }

            double max_water_cut_perf = 0.;
            for (int perf = 0; perf < perf_number; ++perf) {
                if (water_cut_perf[perf] > max_water_cut_perf) {
                    worst_offending_connection = perf;
                    max_water_cut_perf = water_cut_perf[perf];
                }
            }

            assert(max_water_cut_perf != 0.);
            assert(worst_offending_connection >= 0);

            violation_extent = max_water_cut_perf / max_water_cut_limit;
        }

        return water_cut_limit_violated;
    }

} // namespace Opm
