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


#ifndef OPM_BLACKOILWELLMODEL_HEADER_INCLUDED
#define OPM_BLACKOILWELLMODEL_HEADER_INCLUDED

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
#include <opm/autodiff/WellDensitySegmented.hpp>
#include <opm/autodiff/BlackoilPropsAdFromDeck.hpp>
#include <opm/autodiff/BlackoilDetails.hpp>
#include <opm/autodiff/BlackoilModelParameters.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>
#include <opm/autodiff/RateConverter.hpp>
#include <opm/autodiff/WellInterface.hpp>
#include <opm/autodiff/StandardWell.hpp>
#include <opm/autodiff/MultisegmentWell.hpp>
#include<opm/autodiff/SimFIBODetails.hpp>
#include<dune/common/fmatrix.hh>
#include<dune/istl/bcrsmatrix.hh>
#include<dune/istl/matrixmatrix.hh>

#include <opm/material/densead/Math.hpp>

#include <opm/simulators/WellSwitchingLogger.hpp>


namespace Opm {

        /// Class for handling the blackoil well model.
        template<typename TypeTag>
        class BlackoilWellModel {
        public:
            // ---------      Types      ---------
            typedef WellStateFullyImplicitBlackoil WellState;
            typedef BlackoilModelParameters ModelParameters;

            typedef typename GET_PROP_TYPE(TypeTag, Grid)                Grid;
            typedef typename GET_PROP_TYPE(TypeTag, FluidSystem)         FluidSystem;
            typedef typename GET_PROP_TYPE(TypeTag, ElementContext)      ElementContext;
            typedef typename GET_PROP_TYPE(TypeTag, Indices)             Indices;
            typedef typename GET_PROP_TYPE(TypeTag, Simulator)           Simulator;
            typedef typename GET_PROP_TYPE(TypeTag, Scalar)              Scalar;

            static const int numEq = Indices::numEq;
            static const int solventSaturationIdx = Indices::solventSaturationIdx;

            // TODO: where we should put these types, WellInterface or Well Model?
            // or there is some other strategy, like TypeTag
            typedef Dune::FieldVector<Scalar, numEq    > VectorBlockType;
            typedef Dune::BlockVector<VectorBlockType> BVector;

            typedef Ewoms::BlackOilPolymerModule<TypeTag> PolymerModule;

            // For the conversion between the surface volume rate and resrevoir voidage rate
            using RateConverterType = RateConverter::
                SurfaceToReservoirVoidage<FluidSystem, std::vector<int> >;

            BlackoilWellModel(Simulator& ebosSimulator,
                              const ModelParameters& param,
                              const bool terminal_output);

            // compute the well fluxes and assemble them in to the reservoir equations as source terms
            // and in the well equations.
            void assemble(const int iterationIdx,
                                     const double dt);

            // substract Binv(D)rw from r;
            void apply( BVector& r) const;

            // subtract B*inv(D)*C * x from A*x
            void apply(const BVector& x, BVector& Ax) const;

            // apply well model with scaling of alpha
            void applyScaleAdd(const Scalar alpha, const BVector& x, BVector& Ax) const;

            // using the solution x to recover the solution xw for wells and applying
            // xw to update Well State
            void recoverWellSolutionAndUpdateWellState(const BVector& x);

            // Check if well equations is converged.
            bool getWellConvergence(const std::vector<Scalar>& B_avg) const;

            // return all the wells.
            const WellCollection& wellCollection() const;
            // return non const reference to all the wells.
            WellCollection& wellCollection();

            // return the internal well state, ignore the passed one.
            // Used by the legacy code to make it compatible with the legacy well models.
            const WellState& wellState(const WellState& well_state OPM_UNUSED) const;

            // return the internal well state
            const WellState& wellState() const;

            // only use this for restart.
            void setRestartWellState(const WellState& well_state);

            // called at the beginning of a time step
            void beginTimeStep();
            // called at the end of a time step
            void timeStepSucceeded();

            // called at the beginning of a report step
            void beginReportStep(const int time_step);

            // called at the end of a report step
            void endReportStep();

            const SimulatorReport& lastReport() const;

        protected:

            Simulator& ebosSimulator_;
            std::unique_ptr<WellsManager> wells_manager_;
            std::vector< const Well* > wells_ecl_;

            bool wells_active_;

            using WellInterfacePtr = std::unique_ptr<WellInterface<TypeTag> >;
            // a vector of all the wells.
            // eventually, the wells_ above should be gone.
            // the name is just temporary
            // later, might make share_ptr const later.
            std::vector<WellInterfacePtr > well_container_;

            using ConvergenceReport = typename WellInterface<TypeTag>::ConvergenceReport;

            // create the well container
            std::vector<WellInterfacePtr > createWellContainer(const int time_step,
                                   const std::map<std::string, std::vector<int> >& perforation_mapping) const;

            WellState well_state_;
            WellState previous_well_state_;

            const ModelParameters param_;
            bool terminal_output_;
            bool has_solvent_;
            bool has_polymer_;
            std::vector<int> pvt_region_idx_;
            PhaseUsage phase_usage_;
            size_t global_nc_;
            // the number of the cells in the local grid
            size_t number_of_cells_;
            double gravity_;
            std::vector<double> depth_;

            DynamicListEconLimited dynamic_list_econ_limited_;
            std::unique_ptr<RateConverterType> rateConverter_;
            std::unique_ptr<VFPProperties> vfp_properties_;

            SimulatorReport last_report_;

            // used to better efficiency of calcuation
            mutable BVector scaleAddRes_;

            const Wells* wells() const { return wells_manager_->c_wells(); }

            const Schedule& schedule() const
            { return ebosSimulator_.gridManager().schedule(); }

            void updateWellControls();

            void updateGroupControls();

            // setting the well_solutions_ based on well_state.
            void updatePrimaryVariables();

            void setupCompressedToCartesian(const int* global_cell, int number_of_cells, std::map<int,int>& cartesian_to_compressed ) const;

            void computeRepRadiusPerfLength(const Grid& grid);


            void computeAverageFormationFactor(std::vector<double>& B_avg) const;

            void applyVREPGroupControl();

            void computeWellVoidageRates(std::vector<double>& well_voidage_rates,
                                         std::vector<double>& voidage_conversion_coeffs) const;

            // Calculating well potentials for each well
            void computeWellPotentials(std::vector<double>& well_potentials);

            const std::vector<double>& wellPerfEfficiencyFactors() const;

            void calculateEfficiencyFactors();

            // it should be able to go to prepareTimeStep(), however, the updateWellControls() and initPrimaryVariablesEvaluation()
            // makes it a little more difficult. unless we introduce if (iterationIdx != 0) to avoid doing the above functions
            // twice at the beginning of the time step
            /// Calculating the explict quantities used in the well calculation. By explicit, we mean they are cacluated
            /// at the beginning of the time step and no derivatives are included in these quantities
            void calculateExplicitQuantities() const;

            SimulatorReport solveWellEq(const double dt);

            void initPrimaryVariablesEvaluation() const;

            // The number of components in the model.
            int numComponents() const;

            int numWells() const;

            int numPhases() const;

            void resetWellControlFromState() const;

            void assembleWellEq(const double dt,
                                bool only_wells);

            // some preparation work, mostly related to group control and RESV,
            // at the beginning of each time step (Not report step)
            void prepareTimeStep();

            void prepareGroupControl();

            void computeRESV(const std::size_t step);

            void extractLegacyCellPvtRegionIndex_();

            void extractLegacyDepth_();

            /// return true if wells are available in the reservoir
            bool wellsActive() const;

            void setWellsActive(const bool wells_active);

            /// return true if wells are available on this process
            bool localWellsActive() const;

            /// upate the dynamic lists related to economic limits
            void updateListEconLimited(DynamicListEconLimited& list_econ_limited) const;

            void updatePerforationIntensiveQuantities();

        };


} // namespace Opm

#include "BlackoilWellModel_impl.hpp"
#endif
