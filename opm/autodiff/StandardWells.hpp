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
            // ---------      Types      ---------
            using ADB = AutoDiffBlock<double>;
            using Vector = ADB::V;

            // copied from BlackoilModelBase
            // should put to somewhere better
            using DataBlock =  Eigen::Array<double,
                                            Eigen::Dynamic,
                                            Eigen::Dynamic,
                                            Eigen::RowMajor>;
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
            Vector& wellPerforationDensities(); // mutable version kept for BlackoilMultisegmentModel
            const Vector& wellPerforationDensities() const;

            /// Diff to bhp for each well perforation.
            Vector& wellPerforationPressureDiffs(); // mutable version kept for BlackoilMultisegmentModel
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

            template <class ReservoirResidualQuant, class SolutionState>
            void extractWellPerfProperties(const SolutionState& state,
                                           const std::vector<ReservoirResidualQuant>& rq,
                                           const int np,
                                           const BlackoilPropsAdInterface& fluid,
                                           const std::vector<bool>& active,
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

            template <class SolutionState, class WellState>
            void updatePerfPhaseRatesAndPressures(const std::vector<ADB>& cq_s,
                                                  const SolutionState& state,
                                                  WellState& xw) const;

            template <class WellState>
            void updateWellState(const Vector& dwells,
                                 const double gravity,
                                 const double dpmaxrel,
                                 const Opm::PhaseUsage& pu,
                                 const std::vector<bool>& active,
                                 const VFPProperties& vfp_properties,
                                 WellState& well_state);

           template <class WellState>
           void updateWellControls(const Opm::PhaseUsage& pu,
                                   const double gravity,
                                   const VFPProperties& vfp_properties,
                                   const bool terminal_output,
                                   const std::vector<bool>& active,
                                   WellState& xw) const;

           // TODO: should LinearisedBlackoilResidual also be a template class?
           template <class SolutionState>
           void addWellFluxEq(const std::vector<ADB>& cq_s,
                              const SolutionState& state,
                              LinearisedBlackoilResidual& residual);

           // TODO: some parameters, like gravity, maybe it is better to put in the member list
           template <class SolutionState, class WellState>
           void addWellControlEq(const SolutionState& state,
                                 const WellState& xw,
                                 const Vector& aliveWells,
                                 const std::vector<bool> active,
                                 const VFPProperties& vfp_properties,
                                 const double gravity,
                                 LinearisedBlackoilResidual& residual);



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
