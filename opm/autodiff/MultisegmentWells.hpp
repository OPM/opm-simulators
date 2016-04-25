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


#ifndef OPM_MULTISEGMENTWELLS_HEADER_INCLUDED
#define OPM_MULTISEGMENTWELLS_HEADER_INCLUDED

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <cassert>

#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/WellMultiSegment.hpp>


namespace Opm {


    /// Class for handling the multi-segment well model
    class MultisegmentWells {
    public:
            // ---------      Types      ---------
            using ADB = AutoDiffBlock<double>;
            using Vector = ADB::V;

            // Well operations and data needed.
            struct MultisegmentWellOps {
                explicit MultisegmentWellOps(const std::vector<WellMultiSegmentConstPtr>& wells_ms);
                Eigen::SparseMatrix<double> w2p;              // well -> perf (scatter)
                Eigen::SparseMatrix<double> p2w;              // perf -> well (gather)
                Eigen::SparseMatrix<double> w2s;              // well -> segment (scatter)
                Eigen::SparseMatrix<double> s2w;              // segment -> well (gather)
                Eigen::SparseMatrix<double> s2p;              // segment -> perf (scatter)
                Eigen::SparseMatrix<double> p2s;              // perf -> segment (gather)
                Eigen::SparseMatrix<double> s2s_inlets;       // segment -> its inlet segments
                Eigen::SparseMatrix<double> s2s_outlet;       // segment -> its outlet segment
                Eigen::SparseMatrix<double> topseg2w;         // top segment -> well
                AutoDiffMatrix eliminate_topseg;              // change the top segment related to be zero
                std::vector<int> well_cells;                  // the set of perforated cells
                Vector conn_trans_factors;                         // connection transmissibility factors
                bool has_multisegment_wells;                  // flag indicating whether there is any muli-segment well
            };

            // copied from BlackoilModelBase
            // should put to somewhere better
            using DataBlock =  Eigen::Array<double,
                                            Eigen::Dynamic,
                                            Eigen::Dynamic,
                                            Eigen::RowMajor>;

            // ---------  Public methods  ---------
            // TODO: using a vector of WellMultiSegmentConstPtr for now
            // TODO: it should use const Wells or something else later.
            explicit MultisegmentWells(const std::vector<WellMultiSegmentConstPtr>& wells_multisegment);

            const std::vector<WellMultiSegmentConstPtr>& wells() const;
            const MultisegmentWellOps& wellOps() const;

            int numSegment() const { return nseg_total_; };
            int numPerf() const { return nperf_total_; };

    protected:
        const std::vector<WellMultiSegmentConstPtr> wells_multisegment_;
        const MultisegmentWellOps wops_ms_;
        int nseg_total_;
        int nperf_total_;
    };

} // namespace Opm

#endif // OPM_MULTISEGMENTWELLS_HEADER_INCLUDED
