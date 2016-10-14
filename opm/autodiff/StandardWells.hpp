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

#include <dune/common/parallel/mpihelper.hh>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <cassert>
#include <tuple>

#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>

#include <opm/core/wells.h>
#include <opm/core/wells/DynamicListEconLimited.hpp>
#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/autodiff/BlackoilPropsAdInterface.hpp>
#include <opm/simulators/WellSwitchingLogger.hpp>

namespace Opm {


        /// Class for handling the standard well model.
        class StandardWells {
        public:
            struct WellOps {
                explicit WellOps(const Wells* wells);
                Eigen::SparseMatrix<double> w2p;              // well -> perf (scatter)
                Eigen::SparseMatrix<double> p2w;              // perf -> well (gather)
                std::vector<int> well_cells;                  // the set of perforated cells
            };

            // ---------      Types      ---------
            using ADB = AutoDiffBlock<double>;
            using Vector = ADB::V;
            using Communication =
                Dune::CollectiveCommunication<typename Dune::MPIHelper
                                              ::MPICommunicator>;

            // copied from BlackoilModelBase
            // should put to somewhere better
            using DataBlock =  Eigen::Array<double,
                                            Eigen::Dynamic,
                                            Eigen::Dynamic,
                                            Eigen::RowMajor>;
            // ---------  Public methods  ---------
            explicit StandardWells(const Wells* wells_arg);

            void init(const BlackoilPropsAdInterface* fluid_arg,
                      const std::vector<bool>* active_arg,
                      const std::vector<PhasePresence>* pc_arg,
                      const VFPProperties*  vfp_properties_arg,
                      const double gravity_arg,
                      const Vector& depth_arg);

            const WellOps& wellOps() const;

            int numPhases() const { return wells().number_of_phases; };

            const Wells& wells() const;

            const Wells* wellsPointer() const;

            /// return true if wells are available in the reservoir
            bool wellsActive() const;
            void setWellsActive(const bool wells_active);
            /// return true if wells are available on this process
            bool localWellsActive() const;

            int numWellVars() const;

            /// Density of each well perforation
            Vector& wellPerforationDensities(); // mutable version kept for BlackoilMultisegmentModel
            const Vector& wellPerforationDensities() const;

            /// Diff to bhp for each well perforation.
            Vector& wellPerforationPressureDiffs(); // mutable version kept for BlackoilMultisegmentModel
            const Vector& wellPerforationPressureDiffs() const;

            template <class ReservoirResidualQuant, class SolutionState>
            void extractWellPerfProperties(const SolutionState& state,
                                           const std::vector<ReservoirResidualQuant>& rq,
                                           std::vector<ADB>& mob_perfcells,
                                           std::vector<ADB>& b_perfcells) const;

            template <class SolutionState>
            void computeWellFlux(const SolutionState& state,
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
                                 const double dpmaxrel,
                                 WellState& well_state);

            template <class WellState>
            void updateWellControls(WellState& xw) const;

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
                                  LinearisedBlackoilResidual& residual);

            template <class SolutionState, class WellState>
            void computeWellConnectionPressures(const SolutionState& state,
                                                const WellState& xw);

            // state0 is non-constant, while it will not be used outside of the function
            template <class SolutionState, class WellState>
            void
            computeWellPotentials(const std::vector<ADB>& mob_perfcells,
                                  const std::vector<ADB>& b_perfcells,
                                  SolutionState& state0,
                                  WellState& well_state);

            template <class SolutionState>
            void
            variableStateExtractWellsVars(const std::vector<int>& indices,
                                          std::vector<ADB>& vars,
                                          SolutionState& state) const;

            void
            variableStateWellIndices(std::vector<int>& indices,
                                     int& next) const;

            std::vector<int>
            variableWellStateIndices() const;

            template <class WellState>
            void
            variableWellStateInitials(const WellState& xw,
                                      std::vector<Vector>& vars0) const;

            /// If set, computeWellFlux() will additionally store the
            /// total reservoir volume perforation fluxes.
            void setStoreWellPerforationFluxesFlag(const bool store_fluxes);

            /// Retrieves the stored fluxes. It is an error to call this
            /// unless setStoreWellPerforationFluxesFlag(true) has been
            /// called.
            const Vector& getStoredWellPerforationFluxes() const;

            /// upate the dynamic lists related to economic limits
            template<class WellState>
            void
            updateListEconLimited(const Schedule& schedule,
                                  const int current_step,
                                  const Wells* wells,
                                  const WellState& well_state,
                                  DynamicListEconLimited& list_econ_limited) const;

        protected:
            bool wells_active_;
            const Wells*   wells_;
            const WellOps  wops_;

            const BlackoilPropsAdInterface* fluid_;
            const std::vector<bool>*  active_;
            const std::vector<PhasePresence>*  phase_condition_;
            const VFPProperties* vfp_properties_;
            double gravity_;
            // the depth of the all the cell centers
            // for standard Wells, it the same with the perforation depth
            Vector perf_cell_depth_;

            Vector well_perforation_densities_;
            Vector well_perforation_pressure_diffs_;

            bool store_well_perforation_fluxes_;
            Vector well_perforation_fluxes_;

            // protected methods
            template <class SolutionState, class WellState>
            void computePropertiesForWellConnectionPressures(const SolutionState& state,
                                                             const WellState& xw,
                                                             std::vector<double>& b_perf,
                                                             std::vector<double>& rsmax_perf,
                                                             std::vector<double>& rvmax_perf,
                                                             std::vector<double>& surf_dens_perf);

            template <class WellState>
            void computeWellConnectionDensitesPressures(const WellState& xw,
                                                        const std::vector<double>& b_perf,
                                                        const std::vector<double>& rsmax_perf,
                                                        const std::vector<double>& rvmax_perf,
                                                        const std::vector<double>& surf_dens_perf,
                                                        const std::vector<double>& depth_perf,
                                                        const double grav);


            template <class WellState>
            bool checkRateEconLimits(const WellEconProductionLimits& econ_production_limits,
                                     const WellState& well_state,
                                     const int well_number) const;

            using WellMapType = typename WellState::WellMapType;
            using WellMapEntryType = typename WellState::mapentry_t;

            // a tuple type for ratio limit check.
            // first value indicates whether ratio limit is violated, when the ratio limit is not violated, the following three
            // values should not be used.
            // second value indicates whehter there is only one connection left.
            // third value indicates the indx of the worst-offending connection.
            // the last value indicates the extent of the violation for the worst-offending connection, which is defined by
            // the ratio of the actual value to the value of the violated limit.
            using RatioCheckTuple = std::tuple<bool, bool, int, double>;

            enum ConnectionIndex {
                INVALIDCONNECTION = -10000
            };


            template <class WellState>
            RatioCheckTuple checkRatioEconLimits(const WellEconProductionLimits& econ_production_limits,
                                                 const WellState& well_state,
                                                 const WellMapEntryType& map_entry) const;

            template <class WellState>
            RatioCheckTuple checkMaxWaterCutLimit(const WellEconProductionLimits& econ_production_limits,
                                                  const WellState& well_state,
                                                  const WellMapEntryType& map_entry) const;

        };


} // namespace Opm

#include "StandardWells_impl.hpp"

#endif
