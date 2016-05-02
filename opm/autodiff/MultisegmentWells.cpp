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

#include <config.h>

#include <opm/autodiff/MultisegmentWells.hpp>

namespace Opm {


    MultisegmentWells::
    MultisegmentWellOps::MultisegmentWellOps(const std::vector<WellMultiSegmentConstPtr>& wells_ms)
    {
        // no multi-segment wells are involved by default.
        has_multisegment_wells = false;

        if (wells_ms.empty()) {
            return;
        }

        // Count the total number of perforations and segments.
        const int nw = wells_ms.size();
        int total_nperf = 0;
        int total_nseg = 0;
        for (int w = 0; w < nw; ++w) {
            total_nperf += wells_ms[w]->numberOfPerforations();
            total_nseg += wells_ms[w]->numberOfSegments();
            if (wells_ms[w]->isMultiSegmented()) {
                has_multisegment_wells = true;
            }
        }

        // Create well_cells and conn_trans_factors.
        well_cells.reserve(total_nperf);
        conn_trans_factors.resize(total_nperf);
        int well_perf_start = 0;
        for (int w = 0; w < nw; ++w) {
            WellMultiSegmentConstPtr well = wells_ms[w];
            well_cells.insert(well_cells.end(), well->wellCells().begin(), well->wellCells().end());
            const std::vector<double>& perf_trans = well->wellIndex();
            std::copy(perf_trans.begin(), perf_trans.end(), conn_trans_factors.data() + well_perf_start);
            well_perf_start += well->numberOfPerforations();
        }
        assert(well_perf_start == total_nperf);
        assert(int(well_cells.size()) == total_nperf);

        // Create all the operator matrices,
        // using the setFromTriplets() method.
        s2s_inlets = Eigen::SparseMatrix<double>(total_nseg, total_nseg);
        s2s_outlet = Eigen::SparseMatrix<double>(total_nseg, total_nseg);
        s2w        = Eigen::SparseMatrix<double>(nw, total_nseg);
        w2s        = Eigen::SparseMatrix<double>(total_nseg, nw);
        topseg2w   = Eigen::SparseMatrix<double>(nw, total_nseg);
        s2p = Eigen::SparseMatrix<double>(total_nperf, total_nseg);
        p2s = Eigen::SparseMatrix<double>(total_nseg, total_nperf);
        typedef Eigen::Triplet<double> Tri;
        std::vector<Tri> s2s_inlets_vector;
        std::vector<Tri> s2s_outlet_vector;
        std::vector<Tri> s2w_vector;
        std::vector<Tri> w2s_vector;
        std::vector<Tri> topseg2w_vector;
        std::vector<Tri> s2p_vector;
        std::vector<Tri> p2s_vector;
        Vector topseg_zero = Vector::Ones(total_nseg);
        s2s_inlets_vector.reserve(total_nseg);
        s2s_outlet_vector.reserve(total_nseg);
        s2w_vector.reserve(total_nseg);
        w2s_vector.reserve(total_nseg);
        topseg2w_vector.reserve(nw);
        s2p_vector.reserve(total_nperf);
        p2s_vector.reserve(total_nperf);
        int seg_start = 0;
        int perf_start = 0;
        for (int w = 0; w < nw; ++w) {
            const int ns = wells_ms[w]->numberOfSegments();
            const int np = wells_ms[w]->numberOfPerforations();
            for (int seg = 0; seg < ns; ++seg) {
                const int seg_ind = seg_start + seg;
                w2s_vector.push_back(Tri(seg_ind, w, 1.0));
                s2w_vector.push_back(Tri(w, seg_ind, 1.0));
                if (seg == 0) {
                    topseg2w_vector.push_back(Tri(w, seg_ind, 1.0));
                    topseg_zero(seg_ind) = 0.0;
                }
                int seg_outlet = wells_ms[w]->outletSegment()[seg];
                if (seg_outlet >= 0) {
                    const int outlet_ind = seg_start + seg_outlet;
                    s2s_inlets_vector.push_back(Tri(outlet_ind, seg_ind, 1.0));
                    s2s_outlet_vector.push_back(Tri(seg_ind, outlet_ind, 1.0));
                }

                const auto& seg_perf = wells_ms[w]->segmentPerforations()[seg];
                // the number of perforations related to this segment
                const int npseg = seg_perf.size();
                for (int perf = 0; perf < npseg; ++perf) {
                    const int perf_ind = perf_start + seg_perf[perf];
                    s2p_vector.push_back(Tri(perf_ind, seg_ind, 1.0));
                    p2s_vector.push_back(Tri(seg_ind, perf_ind, 1.0));
                }
            }
            seg_start += ns;
            perf_start += np;
        }

        s2s_inlets.setFromTriplets(s2s_inlets_vector.begin(), s2s_inlets_vector.end());
        s2s_outlet.setFromTriplets(s2s_outlet_vector.begin(), s2s_outlet_vector.end());
        w2s.setFromTriplets(w2s_vector.begin(), w2s_vector.end());
        s2w.setFromTriplets(s2w_vector.begin(), s2w_vector.end());
        topseg2w.setFromTriplets(topseg2w_vector.begin(), topseg2w_vector.end());
        s2p.setFromTriplets(s2p_vector.begin(), s2p_vector.end());
        p2s.setFromTriplets(p2s_vector.begin(), p2s_vector.end());

        w2p = Eigen::SparseMatrix<double>(total_nperf, nw);
        p2w = Eigen::SparseMatrix<double>(nw, total_nperf);
        w2p = s2p * w2s;
        p2w = s2w * p2s;

        eliminate_topseg = AutoDiffMatrix(topseg_zero.matrix().asDiagonal());
    }





    MultisegmentWells::
    MultisegmentWells(const std::vector<WellMultiSegmentConstPtr>& wells_ms)
      : wells_multisegment_(wells_ms)
      , wops_ms_(wells_ms)
      , num_phases_(wells_ms.empty()? 0 : wells_ms[0]->numberOfPhases())
      , well_segment_perforation_pressure_diffs_(ADB::null())
      , well_segment_densities_(ADB::null())
      , well_segment_pressures_delta_(ADB::null())
      , segment_comp_surf_volume_initial_(num_phases_)
      , segment_comp_surf_volume_current_(num_phases_, ADB::null())
      , segment_mass_flow_rates_(ADB::null())
      , segment_viscosities_(ADB::null())
    {
        // TODO: repeated with the wellOps's initialization, delete one soon.
        // Count the total number of perforations and segments.
        const int nw = wells_ms.size();
        top_well_segments_.resize(nw);
        int nperf_total = 0;
        int nseg_total = 0;
        for (int w = 0; w < nw; ++w) {
            top_well_segments_[w] = nseg_total;
            nperf_total += wells_ms[w]->numberOfPerforations();
            nseg_total += wells_ms[w]->numberOfSegments();
        }

        nperf_total_ = nperf_total;
        nseg_total_ = nseg_total;
    }





    void
    MultisegmentWells::init(const BlackoilPropsAdInterface* fluid_arg,
                            const std::vector<bool>* active_arg,
                            const std::vector<PhasePresence>* pc_arg)
    {
        fluid_ = fluid_arg;
        active_ = active_arg;
        phase_condition_ = pc_arg;
    }





    const std::vector<WellMultiSegmentConstPtr>&
    MultisegmentWells::wells() const
    {
        return wells_multisegment_;
    }





    const MultisegmentWells::MultisegmentWellOps&
    MultisegmentWells::wellOps() const
    {
        return wops_ms_;
    }





    void
    MultisegmentWells::
    computeSegmentPressuresDelta(const double grav)
    {
        const int nw = wells().size();
        const int nseg_total = nseg_total_;

        if ( !wellOps().has_multisegment_wells ) {
            well_segment_pressures_delta_ = ADB::constant(Vector::Zero(nseg_total));
            well_segment_perforation_pressure_diffs_ = wellOps().s2p * well_segment_pressures_delta_;
            return;
        }

        // calculate the depth difference of the segments
        // TODO: we need to store the following values somewhere to avoid recomputation.
        Vector segment_depth_delta = Vector::Zero(nseg_total);
        int start_segment = 0;
        for (int w = 0; w < nw; ++w) {
            WellMultiSegmentConstPtr well = wells()[w];
            const int nseg = well->numberOfSegments();
            for (int s = 1; s < nseg; ++s) {
                const int s_outlet = well->outletSegment()[s];
                assert(s_outlet >= 0 && s_outlet < nseg);
                segment_depth_delta[s + start_segment] = well->segmentDepth()[s_outlet] - well->segmentDepth()[s];
            }
            start_segment += nseg;
        }
        assert(start_segment == nseg_total);

        const ADB grav_adb = ADB::constant(Vector::Constant(nseg_total, grav));
        well_segment_pressures_delta_ = segment_depth_delta * grav_adb * well_segment_densities_;

        ADB well_segment_perforation_densities = wellOps().s2p * well_segment_densities_;
        well_segment_perforation_pressure_diffs_ = grav * well_segment_perforation_depth_diffs_ * well_segment_perforation_densities;
    }





    void
    MultisegmentWells::
    variableStateWellIndices(std::vector<int>& indices,
                             int& next) const
    {
        indices[Qs] = next++;
        indices[Bhp] = next++;
    }





    std::vector<int>
    MultisegmentWells::
    variableWellStateIndices() const
    {
        // Black oil model standard is 5 equation.
        // For the pure well solve, only the well equations are picked.
        std::vector<int> indices(5, -1);
        int next = 0;

        variableStateWellIndices(indices, next);

        assert(next == 2);
        return indices;
    }

} // end of namespace Opm



