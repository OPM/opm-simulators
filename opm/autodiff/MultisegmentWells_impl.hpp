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
                    const int np,
                    const double dpmaxrel,
                    WellState& well_state) const
    {
        if (!wells().empty())
        {
            const int nw = wells().size();
            const int nseg_total = nseg_total_;

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
                    const Opm::PhaseUsage& pu,
                    const std::vector<bool>& active,
                    const Vector& well_perforation_pressure_diffs,
                    const DataBlock& compi,
                    const std::vector<ADB>& mob_perfcells,
                    const std::vector<ADB>& b_perfcells,
                    const int np,
                    Vector& aliveWells,
                    std::vector<ADB>& cq_s) const
    {
        if (wells().size() == 0) return;

        const int nw = wells().size();

        aliveWells = Vector::Constant(nw, 1.0);

        const int nseg = nseg_total_;
        const int nperf = nperf_total_;

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
            ADB h_nc = msperf_selector.select(wellSegmentPerforationPressureDiffs(),
                                              ADB::constant(well_perforation_pressure_diffs));
            const Vector h_cj = msperf_selector.select(wellPerforationCellPressureDiffs(), Vector::Zero(nperf));

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

            if (active[Oil] && active[Gas]) {
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
                if (phase == Oil && active[Gas]) {
                    const int gaspos = pu.phase_pos[Gas];
                    tmp = tmp - rv_perfcells * cmix_s[gaspos] / d;
                }
                if (phase == Gas && active[Oil]) {
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

}
#endif // OPM_MULTISEGMENTWELLS_IMPL_HEADER_INCLUDED
