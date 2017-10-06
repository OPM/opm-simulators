/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 - 2017 Statoil ASA.
  Copyright 2017 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2016 - 2017 IRIS AS

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


#ifndef OPM_STANDARDWELLSDENSE_HEADER_INCLUDED
#define OPM_STANDARDWELLSDENSE_HEADER_INCLUDED

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <cassert>
#include <tuple>

#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>

#include <opm/core/wells.h>
#include <opm/core/wells/DynamicListEconLimited.hpp>
#include <opm/core/wells/WellCollection.hpp>
#include <opm/core/simulator/SimulatorReport.hpp>
#include <opm/autodiff/VFPProperties.hpp>
#include <opm/autodiff/WellHelpers.hpp>
#include <opm/autodiff/BlackoilModelEnums.hpp>
#include <opm/autodiff/WellDensitySegmented.hpp>
#include <opm/autodiff/BlackoilPropsAdFromDeck.hpp>
#include <opm/autodiff/BlackoilDetails.hpp>
#include <opm/autodiff/BlackoilModelParameters.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>
#include <opm/autodiff/RateConverter.hpp>
#include <opm/autodiff/WellInterface.hpp>
#include <opm/autodiff/StandardWell.hpp>
#include<dune/common/fmatrix.hh>
#include<dune/istl/bcrsmatrix.hh>
#include<dune/istl/matrixmatrix.hh>

#include <opm/material/densead/Math.hpp>

#include <opm/simulators/WellSwitchingLogger.hpp>


namespace Opm {

        /// Class for handling the standard well model.
        template<typename TypeTag>
        class StandardWellsDense {
        public:
            // ---------      Types      ---------
            typedef WellStateFullyImplicitBlackoil WellState;
            typedef BlackoilModelParameters ModelParameters;

            typedef typename GET_PROP_TYPE(TypeTag, Grid)                Grid;
            typedef typename GET_PROP_TYPE(TypeTag, FluidSystem)         FluidSystem;
            typedef typename GET_PROP_TYPE(TypeTag, ElementContext)      ElementContext;
            typedef typename GET_PROP_TYPE(TypeTag, Indices)             BlackoilIndices;
            typedef typename GET_PROP_TYPE(TypeTag, Simulator)           Simulator;
            typedef typename GET_PROP_TYPE(TypeTag, Scalar)              Scalar;

            static const int numEq = BlackoilIndices::numEq;
            static const int solventSaturationIdx = BlackoilIndices::solventSaturationIdx;

            // TODO: where we should put these types, WellInterface or Well Model?
            // or there is some other strategy, like TypeTag
            typedef Dune::FieldVector<Scalar, numEq    > VectorBlockType;
            typedef Dune::BlockVector<VectorBlockType> BVector;

            typedef Ewoms::BlackOilPolymerModule<TypeTag> PolymerModule;

            // For the conversion between the surface volume rate and resrevoir voidage rate
            using RateConverterType = RateConverter::
                SurfaceToReservoirVoidage<BlackoilPropsAdFromDeck::FluidSystem, std::vector<int> >;

            // ---------  Public methods  ---------
            StandardWellsDense(const Wells* wells_arg,
                               WellCollection* well_collection,
                               const std::vector< const Well* >& wells_ecl,
                               const ModelParameters& param,
                               const RateConverterType& rate_converter,
                               const bool terminal_output,
                               const int current_index,
                               std::vector<int>& pvt_region_idx);

            void init(const PhaseUsage phase_usage_arg,
                      const std::vector<bool>& active_arg,
                      const double gravity_arg,
                      const std::vector<double>& depth_arg,
                      long int global_nc,
                      const Grid& grid);

            void setVFPProperties(const VFPProperties*  vfp_properties_arg);


            SimulatorReport assemble(Simulator& ebosSimulator,
                                     const int iterationIdx,
                                     const double dt,
                                     WellState& well_state);

            // substract Binv(D)rw from r;
            void apply( BVector& r) const;

            // subtract B*inv(D)*C * x from A*x
            void apply(const BVector& x, BVector& Ax) const;

            // apply well model with scaling of alpha
            void applyScaleAdd(const Scalar alpha, const BVector& x, BVector& Ax) const;

            // using the solution x to recover the solution xw for wells and applying
            // xw to update Well State
            void recoverWellSolutionAndUpdateWellState(const BVector& x, WellState& well_state) const;

            int numWells() const;

            /// return true if wells are available in the reservoir
            bool wellsActive() const;

            void setWellsActive(const bool wells_active);

            /// return true if wells are available on this process
            bool localWellsActive() const;

            bool getWellConvergence(Simulator& ebosSimulator,
                                    const std::vector<Scalar>& B_avg) const;

            /// upate the dynamic lists related to economic limits
            void updateListEconLimited(const Schedule& schedule,
                                       const int current_step,
                                       const Wells* wells_struct,
                                       const WellState& well_state,
                                       DynamicListEconLimited& list_econ_limited) const;

            WellCollection* wellCollection() const;


        protected:
            bool wells_active_;
            const Wells*   wells_;
            const std::vector< const Well* > wells_ecl_;

            // the number of wells in this process
            // trying not to use things from Wells struct
            // TODO: maybe a better name to emphasize it is local?
            const int number_of_wells_;

            const int number_of_phases_;

            using WellInterfacePtr = std::unique_ptr<WellInterface<TypeTag> >;
            // a vector of all the wells.
            // eventually, the wells_ above should be gone.
            // the name is just temporary
            // later, might make share_ptr const later.
            std::vector<WellInterfacePtr > well_container_;

            using ConvergenceReport = typename WellInterface<TypeTag>::ConvergenceReport;

            // create the well container
            static std::vector<WellInterfacePtr > createWellContainer(const Wells* wells,
                                                                      const std::vector<const Well*>& wells_ecl,
                                                                      const int time_step);

            // Well collection is used to enforce the group control
            WellCollection* well_collection_;

            ModelParameters param_;
            bool terminal_output_;
            bool has_solvent_;
            bool has_polymer_;
            int current_timeIdx_;

            PhaseUsage phase_usage_;
            std::vector<bool>  active_;
            const RateConverterType& rate_converter_;
            std::vector<int> pvt_region_idx_;

            // the number of the cells in the local grid
            int number_of_cells_;

            long int global_nc_;

            // used to better efficiency of calcuation
            mutable BVector scaleAddRes_;

            void updateWellControls(WellState& xw) const;

            void updateGroupControls(WellState& well_state) const;

            // setting the well_solutions_ based on well_state.
            void updatePrimaryVariables(const WellState& well_state) const;

            void setupCompressedToCartesian(const int* global_cell, int number_of_cells, std::map<int,int>& cartesian_to_compressed ) const;

            void computeRepRadiusPerfLength(const Grid& grid);


            void computeAverageFormationFactor(Simulator& ebosSimulator,
                                               std::vector<double>& B_avg) const;

            void applyVREPGroupControl(WellState& well_state) const;

            void computeWellVoidageRates(const WellState& well_state,
                                         std::vector<double>& well_voidage_rates,
                                         std::vector<double>& voidage_conversion_coeffs) const;

            // Calculating well potentials for each well
            // TODO: getBhp() will be refactored to reduce the duplication of the code calculating the bhp from THP.
            void computeWellPotentials(const Simulator& ebosSimulator,
                                       const WellState& well_state,
                                       std::vector<double>& well_potentials) const;

            const std::vector<double>& wellPerfEfficiencyFactors() const;

            void calculateEfficiencyFactors();

            void computeWellConnectionPressures(const Simulator& ebosSimulator,
                                                const WellState& xw) const;

            SimulatorReport solveWellEq(Simulator& ebosSimulator,
                                        const double dt,
                                        WellState& well_state) const;

            void computeAccumWells() const;

            void initPrimaryVariablesEvaluation() const;

            // The number of components in the model.
            int numComponents() const
            {
                if (numPhases() == 2) {
                    return 2;
                }
                int numComp = FluidSystem::numComponents;
                if (has_solvent_) {
                    numComp ++;
                }

                return numComp;
            }

            int numPhases() const;

            int flowPhaseToEbosPhaseIdx( const int phaseIdx ) const;

            void resetWellControlFromState(const WellState& xw) const;

            void assembleWellEq(Simulator& ebosSimulator,
                                const double dt,
                                WellState& well_state,
                                bool only_wells) const;

            // some preparation work, mostly related to group control and RESV,
            // at the beginning of each time step (Not report step)
            void prepareTimeStep(const Simulator& ebos_simulator,
                                 WellState& well_state);
        };


} // namespace Opm

#include "StandardWellsDense_impl.hpp"
#endif
