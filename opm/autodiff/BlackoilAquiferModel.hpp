/*
<<<<<<< HEAD
  File adapted from BlackoilWellModel.hpp

  Copyright 2017 TNO - Heat Transfer & Fluid Dynamics, Modelling & Optimization of the Subsurface
  Copyright 2017 Statoil ASA.
=======
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 - 2017 Statoil ASA.
  Copyright 2017 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2016 - 2017 IRIS AS
>>>>>>> 9ccee28... First addition of the class BlackoilAquiferModel.

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


#ifndef OPM_BLACKOILAQUIFERMODEL_HEADER_INCLUDED
#define OPM_BLACKOILAQUIFERMODEL_HEADER_INCLUDED

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <cassert>
#include <tuple>

#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>

#include <opm/parser/eclipse/EclipseState/AquiferCT.hpp>

#include <opm/core/simulator/SimulatorReport.hpp>

#include <opm/simulators/timestepping/SimulatorTimer.hpp>

#include <opm/autodiff/BlackoilPropsAdFromDeck.hpp>
#include <opm/autodiff/BlackoilDetails.hpp>
#include <opm/autodiff/BlackoilModelParameters.hpp>
#include <opm/autodiff/RateConverter.hpp>

#include <opm/autodiff/AquiferCarterTracy.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Deck/DeckRecord.hpp>
#include <opm/parser/eclipse/Deck/DeckKeyword.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixmatrix.hh>


#include <opm/material/densead/Math.hpp>


namespace Opm {

        /// Class for handling the blackoil well model.
        template<typename TypeTag>
        class BlackoilAquiferModel {
        
        public:


            // ---------      Types      ---------
            typedef BlackoilModelParameters ModelParameters;


            typedef typename GET_PROP_TYPE(TypeTag, Grid)                Grid;
            typedef typename GET_PROP_TYPE(TypeTag, GridView)            GridView;
            typedef typename GET_PROP_TYPE(TypeTag, FluidSystem)         FluidSystem;
            typedef typename GET_PROP_TYPE(TypeTag, ElementContext)      ElementContext;
            typedef typename GET_PROP_TYPE(TypeTag, Indices)             BlackoilIndices;
            typedef typename GET_PROP_TYPE(TypeTag, Simulator)           Simulator;
            typedef typename GET_PROP_TYPE(TypeTag, Scalar)              Scalar;

            static const int numEq = BlackoilIndices::numEq;
            static const int solventSaturationIdx = BlackoilIndices::solventSaturationIdx;

            typedef Ewoms::BlackOilPolymerModule<TypeTag> PolymerModule;

            typedef AquiferCarterTracy<TypeTag> Aquifer_object;

            BlackoilAquiferModel(Simulator& ebosSimulator,
                              const ModelParameters& param,
                              const bool terminal_output);

            // compute the well fluxes and assemble them in to the reservoir equations as source terms
            // and in the well equations.
            void assemble( const SimulatorTimerInterface& timer,
                           const int iterationIdx                );

            // called at the beginning of a time step
            void beginTimeStep();
            // called at the end of a time step
            void timeStepSucceeded();

            // called at the beginning of a report step
            void beginReportStep(const int time_step);

            // called at the end of a report step
            void endReportStep();

            const SimulatorReport& lastReport() const;

            inline const Simulator& simulator() const
            {
                return ebosSimulator_;
            }

            /// Hack function to get what I need from parser
            void init(const Simulator& ebosSimulator, std::vector<Aquifer_object>& aquifers);

        protected:

            Simulator& ebosSimulator_;

            const ModelParameters param_;
            bool terminal_output_;
            bool has_solvent_;
            bool has_polymer_;
            std::vector<int> pvt_region_idx_;
            PhaseUsage phase_usage_;
            std::vector<bool>  active_;
            size_t global_nc_;
            // the number of the cells in the local grid
            size_t number_of_cells_;
            double gravity_;
            std::vector<double> depth_;
            std::vector<Aquifer_object> aquifers_;


            SimulatorReport last_report_;

            const Schedule& schedule() const
            { return ebosSimulator_.gridManager().schedule(); }

            void updatePrimaryVariables();

            void initPrimaryVariablesEvaluation() const;

            void updateConnectionIntensiveQuantities() const;

            void calculateExplicitQuantities();

            // The number of components in the model.
            int numComponents() const;

            int numAquifers() const;

            int numPhases() const;

            int flowPhaseToEbosPhaseIdx( const int phaseIdx ) const;

            void assembleAquiferEq(const SimulatorTimerInterface& timer);

            SimulatorReport solveAquiferEq(const SimulatorTimerInterface& timer);

            // some preparation work, mostly related to group control and RESV,
            // at the beginning of each time step (Not report step)
            void prepareTimeStep(const SimulatorTimerInterface& timer);

            const std::vector<Aquifer_object>& aquifers();

        };


} // namespace Opm

#include "BlackoilAquiferModel_impl.hpp"
#endif
