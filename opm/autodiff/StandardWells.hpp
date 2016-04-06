/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.

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


namespace Opm {

        // ---------      Types      ---------
        typedef AutoDiffBlock<double> ADB;
        typedef ADB::V Vector;
        // typedef ADB::M Matrix;

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

        protected:
            bool wells_active_;
            const Wells*   wells_;
            const WellOps  wops_;
            Vector well_perforation_densities_;
            Vector well_perforation_pressure_diffs_;
        };


} // namespace Opm
#endif
