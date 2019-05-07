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

#include <ebos/eclbaseaquifermodel.hh>

#include <opm/parser/eclipse/EclipseState/AquiferCT.hpp>
#include <opm/parser/eclipse/EclipseState/Aquifetp.hpp>
#include <opm/parser/eclipse/EclipseState/Aquancon.hpp>
#include <opm/simulators/timestepping/SimulatorTimer.hpp>
#include <opm/simulators/aquifers/AquiferCarterTracy.hpp>
#include <opm/simulators/aquifers/AquiferFetkovich.hpp>
#include <opm/material/densead/Math.hpp>

namespace Opm {

        /// Class for handling the blackoil well model.
        template<typename TypeTag>
        class BlackoilAquiferModel
        {
            typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
            typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;

        public:
            explicit BlackoilAquiferModel(Simulator& simulator);

            void initialSolutionApplied();
            void beginEpisode();
            void beginTimeStep();
            void beginIteration();
            // add the water rate due to aquifers to the source term.
            template <class Context>
            void addToSource(RateVector& rates,
                             const Context& context,
                             unsigned spaceIdx,
                             unsigned timeIdx) const;
            void endIteration();
            void endTimeStep();
            void endEpisode();

            template <class Restarter>
            void serialize(Restarter& res);

            template <class Restarter>
            void deserialize(Restarter& res);

        protected:
            // ---------      Types      ---------
            typedef typename GET_PROP_TYPE(TypeTag, ElementContext)      ElementContext;
            typedef typename GET_PROP_TYPE(TypeTag, Scalar)              Scalar;

            typedef AquiferCarterTracy<TypeTag> AquiferCarterTracy_object;
            typedef AquiferFetkovich<TypeTag> AquiferFetkovich_object;

            Simulator& simulator_;

            std::unordered_map<int, int> cartesian_to_compressed_;
            mutable  std::vector<AquiferCarterTracy_object> aquifers_CarterTracy;
            mutable  std::vector<AquiferFetkovich_object> aquifers_Fetkovich;

            // This initialization function is used to connect the parser objects with the ones needed by AquiferCarterTracy
            void init();

            bool aquiferActive() const;
            bool aquiferCarterTracyActive() const;
            bool aquiferFetkovichActive() const;

        };


} // namespace Opm

#include "BlackoilAquiferModel_impl.hpp"
#endif
