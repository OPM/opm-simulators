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

#ifndef OPM_MULTISEGMENTWELLS_IMPL_HEADER_INCLUDED
#define OPM_MULTISEGMENTWELLS_IMPL_HEADER_INCLUDED


namespace Opm
{



    namespace wellhelpers {

        using ADB = MultisegmentWells::ADB;
        using Vector = MultisegmentWells::Vector;

        inline
        ADB onlyWellDerivs(const ADB& x)
        {
            Vector val = x.value();
            const int nb = x.numBlocks();
            if (nb < 2) {
                OPM_THROW(std::logic_error, "Called onlyWellDerivs() with argument that has " << nb << " blocks.");
            }
            std::vector<ADB::M> derivs = { x.derivative()[nb - 2], x.derivative()[nb - 1] };
            return ADB::function(std::move(val), std::move(derivs));
        }
    }



    template <class WellState>
    void
    MultisegmentWells::
    updateWellState(const Vector& dwells,
                    const double dpmaxrel,
                    WellState& well_state) const
    {
        if (!wells().empty())
        {
            const int nw = wells().size();
            const int nseg_total = nseg_total_;
            const int np = numPhases();

            // Extract parts of dwells corresponding to each part.
            int varstart = 0;
            const Vector dsegqs = subset(dwells, Span(np * nseg_total, 1, varstart));
            varstart += dsegqs.size();
            const Vector dsegp = subset(dwells, Span(nseg_total, 1, varstart));
            varstart += dsegp.size();
            assert(varstart == dwells.size());


            // segment phase rates update
            // in dwells, the phase rates are ordered by phase.
            // while in WellStateMultiSegment, the phase rates are ordered by segments
            const DataBlock wsr = Eigen::Map<const DataBlock>(dsegqs.data(), np, nseg_total).transpose();
            const Vector dwsr = Eigen::Map<const Vector>(wsr.data(), nseg_total * np);
            const Vector wsr_old = Eigen::Map<const Vector>(&well_state.segPhaseRates()[0], nseg_total * np);
            const Vector sr = wsr_old - dwsr;
            std::copy(&sr[0], &sr[0] + sr.size(), well_state.segPhaseRates().begin());


            // segment pressure updates
            const Vector segp_old = Eigen::Map<const Vector>(&well_state.segPress()[0], nseg_total, 1);
            // TODO: applying the pressure change limiter to all the segments, not sure if it is the correct thing to do
            const Vector dsegp_limited = sign(dsegp) * dsegp.abs().min(segp_old.abs() * dpmaxrel);
            const Vector segp = segp_old - dsegp_limited;
            std::copy(&segp[0], &segp[0] + segp.size(), well_state.segPress().begin());

            // update the well rates and bhps, which are not anymore primary vabriables.
            // they are updated directly from the updated segment phase rates and segment pressures.

            // Bhp update.
            Vector bhp = Vector::Zero(nw);
            Vector wr = Vector::Zero(nw * np);
            // it is better to use subset

            int start_segment = 0;
            for (int w = 0; w < nw; ++w) {
                bhp[w] = well_state.segPress()[start_segment];
                // insert can be faster
                for (int p = 0; p < np; ++p) {
                    wr[p + np * w] = well_state.segPhaseRates()[p + np * start_segment];
                }

                const int nseg = wells()[w]->numberOfSegments();
                start_segment += nseg;
            }

            assert(start_segment == nseg_total);
            std::copy(&bhp[0], &bhp[0] + bhp.size(), well_state.bhp().begin());
            std::copy(&wr[0], &wr[0] + wr.size(), well_state.wellRates().begin());

            // TODO: handling the THP control related.
        }
    }





    template <class SolutionState>
    void
    MultisegmentWells::
    computeWellFlux(const SolutionState& state,
                    const Vector& well_perforation_pressure_diffs,
                    const DataBlock& compi,
                    const std::vector<ADB>& mob_perfcells,
                    const std::vector<ADB>& b_perfcells,
                    Vector& aliveWells,
                    std::vector<ADB>& cq_s) const
    {
        if (wells().size() == 0) return;

        const int np = numPhases();
        const int nw = wells().size();

        aliveWells = Vector::Constant(nw, 1.0);

        const int nseg = nseg_total_;
        const int nperf = nperf_total_;

        const Opm::PhaseUsage& pu = fluid_->phaseUsage();

        cq_s.resize(np, ADB::null());

        {
            const Vector& Tw = wellOps().conn_trans_factors;
            const std::vector<int>& well_cells = wellOps().well_cells;

            // determining in-flow (towards well-bore) or out-flow (towards reservoir)
            // for mutli-segmented wells and non-segmented wells, the calculation of the drawdown are different.
            const ADB& p_perfcells = subset(state.pressure, well_cells);
            const ADB& rs_perfcells = subset(state.rs, well_cells);
            const ADB& rv_perfcells = subset(state.rv, well_cells);

            const ADB& seg_pressures = state.segp;

            const ADB seg_pressures_perf = wellOps().s2p * seg_pressures;

            // Create selector for perforations of multi-segment vs. regular wells.
            Vector is_multisegment_well(nw);
            for (int w = 0; w < nw; ++w) {
                is_multisegment_well[w] = double(wells()[w]->isMultiSegmented());
            }
            // Take one flag per well and expand to one flag per perforation.
            Vector is_multisegment_perf = wellOps().w2p * is_multisegment_well.matrix();
            Selector<double> msperf_selector(is_multisegment_perf, Selector<double>::NotEqualZero);

            // Compute drawdown.
            ADB h_nc = msperf_selector.select(well_segment_perforation_pressure_diffs_,
                                              ADB::constant(well_perforation_pressure_diffs));
            const Vector h_cj = msperf_selector.select(well_perforation_cell_pressure_diffs_, Vector::Zero(nperf));

            // Special handling for when we are called from solveWellEq().
            // TODO: restructure to eliminate need for special treatmemt.
            if ((h_nc.numBlocks() != 0) && (h_nc.numBlocks() != seg_pressures_perf.numBlocks())) {
                assert(seg_pressures_perf.numBlocks() == 2);
                assert(h_nc.numBlocks() > 2);
                h_nc = wellhelpers::onlyWellDerivs(h_nc);
                assert(h_nc.numBlocks() == 2);
            }

            ADB drawdown = (p_perfcells + h_cj - seg_pressures_perf - h_nc);

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

            // handling flow into wellbore
            // maybe there are something to do there make the procedure easier.
            std::vector<ADB> cq_ps(np, ADB::null());
            for (int phase = 0; phase < np; ++phase) {
                const ADB cq_p = -(selectProducingPerforations * Tw) * (mob_perfcells[phase] * drawdown);
                cq_ps[phase] = b_perfcells[phase] * cq_p;
            }

            if ((*active_)[Oil] && (*active_)[Gas]) {
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
            std::vector<ADB> wbq(np, ADB::null());
            ADB wbqt = ADB::constant(Vector::Zero(nseg));

            for (int phase = 0; phase < np; ++phase) {
                const ADB& q_ps = wellOps().p2s * cq_ps[phase];
                const ADB& q_s = subset(state.segqs, Span(nseg, 1, phase * nseg));
                Selector<double> injectingPhase_selector(q_s.value(), Selector<double>::GreaterZero);

                const int pos = pu.phase_pos[phase];

                // this is per segment
                wbq[phase] = (wellOps().w2s * ADB::constant(compi.col(pos)) * injectingPhase_selector.select(q_s, ADB::constant(Vector::Zero(nseg)))) - q_ps;

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
                    topseg += wells()[w]->numberOfSegments();
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
                const ADB phase_fraction = wellOps().topseg2w * (wbq[phase] / wbqt);
                cmix_s[phase] = wellOps().w2p * aliveWells_selector.select(phase_fraction, ADB::constant(compi.col(pos)));
            }

            // compute volume ration between connection at standard conditions
            ADB volumeRatio = ADB::constant(Vector::Zero(nperf));
            const ADB d = Vector::Constant(nperf,1.0) -  rv_perfcells * rs_perfcells;

            for (int phase = 0; phase < np; ++phase) {
                ADB tmp = cmix_s[phase];
                if (phase == Oil && (*active_)[Gas]) {
                    const int gaspos = pu.phase_pos[Gas];
                    tmp = tmp - rv_perfcells * cmix_s[gaspos] / d;
                }
                if (phase == Gas && (*active_)[Oil]) {
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





    template <class SolutionState>
    void
    MultisegmentWells::
    computeSegmentFluidProperties(const SolutionState& state)
    {
        const int np = numPhases();
        const int nw = wells().size();
        const int nseg_total = nseg_total_;

        if ( !wellOps().has_multisegment_wells ){
            // not sure if this is needed actually
            // TODO: to check later if this is really necessary.
            well_segment_densities_ = ADB::constant(Vector::Zero(nseg_total));
            segment_mass_flow_rates_ = ADB::constant(Vector::Zero(nseg_total));
            segment_viscosities_ = ADB::constant(Vector::Zero(nseg_total));
            for (int phase = 0; phase < np; ++phase) {
                segment_comp_surf_volume_current_[phase] = ADB::constant(Vector::Zero(nseg_total));
                segmentCompSurfVolumeInitial()[phase] = Vector::Zero(nseg_total);
            }
            return;
        }

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
            const std::vector<int>& segment_cells_well = wells()[w]->segmentCells();
            segment_cells.insert(segment_cells.end(), segment_cells_well.begin(), segment_cells_well.end());
        }
        assert(int(segment_cells.size()) == nseg_total);

        const ADB segment_temp = subset(state.temperature, segment_cells);
        // using the segment pressure or the average pressure
        // using the segment pressure first
        const ADB& segment_press = state.segp;

        // Compute PVT properties for segments.
        std::vector<PhasePresence> segment_cond(nseg_total);
        for (int s = 0; s < nseg_total; ++s) {
            segment_cond[s] = (*phase_condition_)[segment_cells[s]];
        }
        std::vector<ADB> b_seg(np, ADB::null());
        // Viscosities for different phases
        std::vector<ADB> mu_seg(np, ADB::null());
        ADB rsmax_seg = ADB::null();
        ADB rvmax_seg = ADB::null();
        const PhaseUsage& pu = fluid_->phaseUsage();
        if (pu.phase_used[Water]) {
            b_seg[pu.phase_pos[Water]] = fluid_->bWat(segment_press, segment_temp, segment_cells);
            mu_seg[pu.phase_pos[Water]] = fluid_->muWat(segment_press, segment_temp, segment_cells);
        }
        assert((*active_)[Oil]);
        const ADB segment_so = subset(state.saturation[pu.phase_pos[Oil]], segment_cells);
        if (pu.phase_used[Oil]) {
            const ADB segment_rs = subset(state.rs, segment_cells);
            b_seg[pu.phase_pos[Oil]] = fluid_->bOil(segment_press, segment_temp, segment_rs,
                                                   segment_cond, segment_cells);
            // rsmax_seg = fluidRsSat(segment_press, segment_so, segment_cells);
            rsmax_seg = fluid_->rsSat(segment_press, segment_so, segment_cells);
            mu_seg[pu.phase_pos[Oil]] = fluid_->muOil(segment_press, segment_temp, segment_rs,
                                                     segment_cond, segment_cells);
        }
        assert((*active_)[Gas]);
        if (pu.phase_used[Gas]) {
            const ADB segment_rv = subset(state.rv, segment_cells);
            b_seg[pu.phase_pos[Gas]] = fluid_->bGas(segment_press, segment_temp, segment_rv,
                                                   segment_cond, segment_cells);
            // rvmax_seg = fluidRvSat(segment_press, segment_so, segment_cells);
            rvmax_seg = fluid_->rvSat(segment_press, segment_so, segment_cells);
            mu_seg[pu.phase_pos[Gas]] = fluid_->muGas(segment_press, segment_temp, segment_rv,
                                                   segment_cond, segment_cells);
        }

        // Extract segment flow by phase (segqs) and compute total surface rate.
        ADB tot_surface_rate = ADB::constant(Vector::Zero(nseg_total));
        std::vector<ADB> segqs(np, ADB::null());
        for (int phase = 0; phase < np; ++phase) {
            segqs[phase] = subset(state.segqs, Span(nseg_total, 1, phase * nseg_total));
            tot_surface_rate += segqs[phase];
        }

        // TODO: later this will be implmented as a global mapping
        std::vector<std::vector<double>> comp_frac(np, std::vector<double>(nseg_total, 0.0));
        int start_segment = 0;
        for (int w = 0; w < nw; ++w) {
            WellMultiSegmentConstPtr well = wells()[w];
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
            mix[phase] = ADB::constant(Eigen::Map<Vector>(comp_frac[phase].data(), nseg_total));
        }
        // There should be a better way to do this.
        Selector<double> non_zero_tot_rate(tot_surface_rate.value(), Selector<double>::NotEqualZero);
        for (int phase = 0; phase < np; ++phase) {
            mix[phase] = non_zero_tot_rate.select(segqs[phase] / tot_surface_rate, mix[phase]);
        }

        // Calculate rs and rv.
        ADB rs = ADB::constant(Vector::Zero(nseg_total));
        ADB rv = rs;
        const int gaspos = pu.phase_pos[Gas];
        const int oilpos = pu.phase_pos[Oil];
        Selector<double> non_zero_mix_oilpos(mix[oilpos].value(), Selector<double>::GreaterZero);
        Selector<double> non_zero_mix_gaspos(mix[gaspos].value(), Selector<double>::GreaterZero);
        // What is the better way to do this?
        // big values should not be necessary
        ADB big_values = ADB::constant(Vector::Constant(nseg_total, 1.e100));
        ADB mix_gas_oil = non_zero_mix_oilpos.select(mix[gaspos] / mix[oilpos], big_values);
        ADB mix_oil_gas = non_zero_mix_gaspos.select(mix[oilpos] / mix[gaspos], big_values);
        if ((*active_)[Oil]) {
            Vector selectorUnderRsmax = Vector::Zero(nseg_total);
            Vector selectorAboveRsmax = Vector::Zero(nseg_total);
            for (int s = 0; s < nseg_total; ++s) {
                if (mix_gas_oil.value()[s] > rsmax_seg.value()[s]) {
                    selectorAboveRsmax[s] = 1.0;
                } else {
                    selectorUnderRsmax[s] = 1.0;
                }
            }
            rs = non_zero_mix_oilpos.select(selectorAboveRsmax * rsmax_seg + selectorUnderRsmax * mix_gas_oil, rs);
        }
        if ((*active_)[Gas]) {
            Vector selectorUnderRvmax = Vector::Zero(nseg_total);
            Vector selectorAboveRvmax = Vector::Zero(nseg_total);
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
        if ((*active_)[Gas] && (*active_)[Oil]) {
            x[gaspos] = (mix[gaspos] - mix[oilpos] * rs) / (Vector::Ones(nseg_total) - rs * rv);
            x[oilpos] = (mix[oilpos] - mix[gaspos] * rv) / (Vector::Ones(nseg_total) - rs * rv);
        }

        // Compute total reservoir volume to surface volume ratio.
        ADB volrat = ADB::constant(Vector::Zero(nseg_total));
        for (int phase = 0; phase < np; ++phase) {
            volrat += x[phase] / b_seg[phase];
        }

        // Compute segment densities.
        ADB dens = ADB::constant(Vector::Zero(nseg_total));
        for (int phase = 0; phase < np; ++phase) {
            const Vector surface_density = fluid_->surfaceDensity(phase, segment_cells);
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
        segment_mass_flow_rates_ = ADB::constant(Vector::Zero(nseg_total));
        for (int phase = 0; phase < np; ++phase) {
            // TODO: how to remove one repeated surfaceDensity()
            const Vector surface_density = fluid_->surfaceDensity(phase, segment_cells);
            segment_mass_flow_rates_ += surface_density * segqs[phase];
        }

        // Viscosity of the fluid mixture in the segments
        segment_viscosities_ = ADB::constant(Vector::Zero(nseg_total));
        for (int phase = 0; phase < np; ++phase) {
            segment_viscosities_ += x[phase] * mu_seg[phase];
        }
    }





    template <class SolutionState>
    void
    MultisegmentWells::
    addWellFluxEq(const std::vector<ADB>& cq_s,
                  const SolutionState& state,
                  LinearisedBlackoilResidual& residual)
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
        const int nseg_total = nseg_total_;

        ADB segqs = state.segqs;

        std::vector<ADB> segment_volume_change_dt(np, ADB::null());
        for (int phase = 0; phase < np; ++phase) {
            if ( wellOps().has_multisegment_wells ) {
                // Gain of the surface volume of each component in the segment by dt
                segment_volume_change_dt[phase] = segment_comp_surf_volume_current_[phase] -
                                                  segmentCompSurfVolumeInitial()[phase];

                // Special handling for when we are called from solveWellEq().
                // TODO: restructure to eliminate need for special treatmemt.
                if (segment_volume_change_dt[phase].numBlocks() != segqs.numBlocks()) {
                    assert(segment_volume_change_dt[phase].numBlocks() > 2);
                    assert(segqs.numBlocks() == 2);
                    segment_volume_change_dt[phase] = wellhelpers::onlyWellDerivs(segment_volume_change_dt[phase]);
                    assert(segment_volume_change_dt[phase].numBlocks() == 2);
                }

                const ADB cq_s_seg = wellOps().p2s * cq_s[phase];
                const ADB segqs_phase = subset(segqs, Span(nseg_total, 1, phase * nseg_total));
                segqs -= superset(cq_s_seg + wellOps().s2s_inlets * segqs_phase + segment_volume_change_dt[phase],
                                  Span(nseg_total, 1, phase * nseg_total), np * nseg_total);
            } else {
                segqs -= superset(wellOps().p2s * cq_s[phase], Span(nseg_total, 1, phase * nseg_total), np * nseg_total);
            }
        }

        residual.well_flux_eq = segqs;
    }





    template <class SolutionState, class WellState>
    void
    MultisegmentWells::
    addWellControlEq(const SolutionState& state,
                     const WellState& xw,
                     const Vector& aliveWells,
                     LinearisedBlackoilResidual& residual)
    {
        // the name of the function is a a little misleading.
        // Basically it is the function for the pressure equation.
        // And also, it work as the control equation when it is the segment
        if( wells().empty() ) return;

        const int np = numPhases();
        const int nw = wells().size();
        const int nseg_total = nseg_total_;

        ADB aqua   = ADB::constant(Vector::Zero(nseg_total));
        ADB liquid = ADB::constant(Vector::Zero(nseg_total));
        ADB vapour = ADB::constant(Vector::Zero(nseg_total));

        if ((*active_)[Water]) {
            aqua += subset(state.segqs, Span(nseg_total, 1, BlackoilPhases::Aqua * nseg_total));
        }
        if ((*active_)[Oil]) {
            liquid += subset(state.segqs, Span(nseg_total, 1, BlackoilPhases::Liquid * nseg_total));
        }
        if ((*active_)[Gas]) {
            vapour += subset(state.segqs, Span(nseg_total, 1, BlackoilPhases::Vapour * nseg_total));
        }

        // THP control is not implemented for the moment.

        // Hydrostatic correction variables
        Vector rho_v = Vector::Zero(nw);
        Vector vfp_ref_depth_v = Vector::Zero(nw);

        // Target vars
        Vector bhp_targets  = Vector::Zero(nw);
        Vector rate_targets = Vector::Zero(nw);
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
            const struct WellControls* wc = wells()[w]->wellControls();

            // The current control in the well state overrides
            // the current control set in the Wells struct, which
            // is instead treated as a default.
            const int current = xw.currentControls()[w];

            const int nseg = wells()[w]->numberOfSegments();

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

        ADB others_residual = ADB::constant(Vector::Zero(nseg_total));

        if ( wellOps().has_multisegment_wells ) {
            // Special handling for when we are called from solveWellEq().
            // TODO: restructure to eliminate need for special treatmemt.
            ADB wspd = (state.segp.numBlocks() == 2)
                ? wellhelpers::onlyWellDerivs(well_segment_pressures_delta_)
                : well_segment_pressures_delta_;

            others_residual = wellOps().eliminate_topseg * (state.segp - wellOps().s2s_outlet * state.segp + wspd);
        } else {
            others_residual = wellOps().eliminate_topseg * (state.segp - wellOps().s2s_outlet * state.segp);
        }

        //       all the control equations
        // TODO: can be optimized better
        ADB well_eq_topsegment = subset(superset(bhp_residual, bhp_top_elems, nseg_total) +
                                        superset(rate_residual, rate_top_elems, nseg_total), top_well_segments_);

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
        residual.well_eq = superset(well_eq_topsegment, top_well_segments_, nseg_total) +
                            others_residual;
    }





    template <class WellState>
    void
    MultisegmentWells::
    updateWellControls(const bool terminal_output,
                       WellState& xw) const
    {
        if( wells().empty() ) return ;

        std::string modestring[4] = { "BHP", "THP", "RESERVOIR_RATE", "SURFACE_RATE" };
        // Find, for each well, if any constraints are broken. If so,
        // switch control to first broken constraint.
        const int np = numPhases();
        const int nw = wells().size();
        for (int w = 0; w < nw; ++w) {
            const WellControls* wc = wells()[w]->wellControls();
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
                        w, np, wells()[w]->wellType(), wc, ctrl_index)) {
                    // ctrl_index will be the index of the broken constraint after the loop.
                    break;
                }
            }

            if (ctrl_index != nwc) {
                // Constraint number ctrl_index was broken, switch to it.
                if (terminal_output)
                {
                    std::cout << "Switching control mode for well " << wells()[w]->name()
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
                xw.segPress()[top_well_segments_[w]] = target;
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


}
#endif // OPM_MULTISEGMENTWELLS_IMPL_HEADER_INCLUDED
