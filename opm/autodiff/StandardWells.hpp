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


#ifndef OPM_STANDARDWELLS_HEADER_INCLUDED
#define OPM_STANDARDWELLS_HEADER_INCLUDED

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <cassert>

#include <opm/core/wells.h>
#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/autodiff/BlackoilPropsAdInterface.hpp>


namespace Opm {

        // ---------      Types      ---------
        typedef AutoDiffBlock<double> ADB;
        typedef ADB::V Vector;

        // copied from BlackoilModelBase
        // should put to somewhere better
        typedef Eigen::Array<double,
                     Eigen::Dynamic,
                     Eigen::Dynamic,
                     Eigen::RowMajor> DataBlock;

        /// Class for handling the standard well model.
        class StandardWells {
        protected:
            struct WellOps {
                explicit WellOps(const Wells* wells);
                Eigen::SparseMatrix<double> w2p;              // well -> perf (scatter)
                Eigen::SparseMatrix<double> p2w;              // perf -> well (gather)
                std::vector<int> well_cells;                  // the set of perforated cells
            };

        public:
            // ---------  Public methods  ---------
            explicit StandardWells(const Wells* wells);

            const Wells& wells() const;

            /// return true if wells are available in the reservoir
            bool wellsActive() const;
            void setWellsActive(const bool wells_active);
            /// return true if wells are available on this process
            bool localWellsActive() const;

            const WellOps& wellOps() const;

            /// Density of each well perforation
            Vector& wellPerforationDensities();
            const Vector& wellPerforationDensities() const;

            /// Diff to bhp for each well perforation.
            Vector& wellPerforationPressureDiffs();
            const Vector& wellPerforationPressureDiffs() const;

            template <class SolutionState, class WellState>
            void computePropertiesForWellConnectionPressures(const SolutionState& state,
                                                             const WellState& xw,
                                                             const BlackoilPropsAdInterface& fluid,
                                                             const std::vector<bool>& active,
                                                             const std::vector<PhasePresence>& pc,
                                                             std::vector<double>& b_perf,
                                                             std::vector<double>& rsmax_perf,
                                                             std::vector<double>& rvmax_perf,
                                                             std::vector<double>& surf_dens_perf);

            template <class WellState>
            void computeWellConnectionDensitesPressures(const WellState& xw,
                                                        const BlackoilPropsAdInterface& fluid,
                                                        const std::vector<double>& b_perf,
                                                        const std::vector<double>& rsmax_perf,
                                                        const std::vector<double>& rvmax_perf,
                                                        const std::vector<double>& surf_dens_perf,
                                                        const std::vector<double>& depth_perf,
                                                        const double grav);

            template <class ReservoirResidualQuant>
            void extractWellPerfProperties(const std::vector<ReservoirResidualQuant>& rq,
                                           const int np,
                                           std::vector<ADB>& mob_perfcells,
                                           std::vector<ADB>& b_perfcells) const;

            template <class SolutionState>
            void computeWellFlux(const SolutionState& state,
                                 const Opm::PhaseUsage& phase_usage,
                                 const std::vector<bool>& active,
                                 const std::vector<ADB>& mob_perfcells,
                                 const std::vector<ADB>& b_perfcells,
                                 Vector& aliveWells,
                                 std::vector<ADB>& cq_s) const;



        protected:
            bool wells_active_;
            const Wells*   wells_;
            const WellOps  wops_;
            Vector well_perforation_densities_;
            Vector well_perforation_pressure_diffs_;
        };


} // namespace Opm

#include "StandardWells_impl.hpp"

#endif
