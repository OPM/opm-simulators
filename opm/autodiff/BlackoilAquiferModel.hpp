/*
  File adapted from BlackoilWellModel.hpp

  Copyright 2017 TNO - Heat Transfer & Fluid Dynamics, Modelling & Optimization of the Subsurface
  Copyright 2017 Statoil ASA.

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
#include <opm/parser/eclipse/EclipseState/Aquancon.hpp>

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

            typedef AquiferCarterTracy<TypeTag> Aquifer_object;

            BlackoilAquiferModel(Simulator& ebosSimulator,
                              const ModelParameters& param,
                              const bool terminal_output);

            // compute the well fluxes and assemble them in to the reservoir equations as source terms
            // and in the well equations.
            void assemble( const SimulatorTimerInterface& timer,
                           const int iterationIdx                );

            // called at the end of a time step
            void timeStepSucceeded(const SimulatorTimerInterface& timer);

            inline const Simulator& simulator() const
            {
                return ebosSimulator_;
            }

            // This initialization function is used to connect the parser objects with the ones needed by AquiferCarterTracy
            void init(const Simulator& ebosSimulator, std::vector<Aquifer_object>& aquifers);

        protected:

            Simulator& ebosSimulator_;

            const ModelParameters param_;
            bool terminal_output_;

            double gravity_;
            std::vector<Aquifer_object> aquifers_;


            void updateConnectionIntensiveQuantities() const;

            int numAquifers() const;

            void assembleAquiferEq(const SimulatorTimerInterface& timer);

            // at the beginning of each time step (Not report step)
            void prepareTimeStep(const SimulatorTimerInterface& timer);

            const std::vector<Aquifer_object>& aquifers();

        };


} // namespace Opm

#include "BlackoilAquiferModel_impl.hpp"
#endif