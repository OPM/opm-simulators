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

#include <opm/parser/eclipse/EclipseState/AquiferCT.hpp>
#include <opm/parser/eclipse/EclipseState/Aquancon.hpp>
#include <opm/simulators/timestepping/SimulatorTimer.hpp>
#include <opm/autodiff/AquiferCarterTracy.hpp>
#include <opm/material/densead/Math.hpp>

namespace Opm {

        /// Class for handling the blackoil well model.
        template<typename TypeTag>
        class BlackoilAquiferModel {
        
        public:

            // ---------      Types      ---------
            typedef typename GET_PROP_TYPE(TypeTag, ElementContext)      ElementContext;
            typedef typename GET_PROP_TYPE(TypeTag, Simulator)           Simulator;
            typedef typename GET_PROP_TYPE(TypeTag, Scalar)              Scalar;

            typedef AquiferCarterTracy<TypeTag> Aquifer_object;

            explicit BlackoilAquiferModel(Simulator& ebosSimulator);

            // compute the well fluxes and assemble them in to the reservoir equations as source terms
            // and in the well equations.
            void assemble( const SimulatorTimerInterface& timer,
                           const int iterationIdx                );

            // called at the end of a time step
            void timeStepSucceeded(const SimulatorTimerInterface& timer);

        protected:

            Simulator& ebosSimulator_;

            std::vector<Aquifer_object> aquifers_;

            // This initialization function is used to connect the parser objects with the ones needed by AquiferCarterTracy
            void init();

            void updateConnectionIntensiveQuantities() const;

            void assembleAquiferEq(const SimulatorTimerInterface& timer);

            // at the beginning of each time step (Not report step)
            void prepareTimeStep(const SimulatorTimerInterface& timer);

            bool aquiferActive() const;

        };


} // namespace Opm

#include "BlackoilAquiferModel_impl.hpp"
#endif