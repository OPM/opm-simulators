/*
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

#include "config.h"

#include <opm/material/common/ResetLocale.hpp>
#include <opm/simulators/flow/SimulatorFullyImplicitBlackoilEbos.hpp>
#include <opm/simulators/flow/Main.hpp>
// modifications from standard
#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/blackoil/blackoilintensivequantitiessimple.hh>
#include <opm/models/discretization/common/smallelementcontext.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>

namespace Opm {
    namespace Properties {
        namespace TTag {
            struct EclFlowProblemTPFA {
            using InheritsFrom = std::tuple<EclFlowProblem>;
            };
        }
   }
}

namespace Opm {
    namespace Properties {

        template<class TypeTag>
        struct Linearizer<TypeTag, TTag::EclFlowProblemTPFA> { using type = TpfaLinearizer<TypeTag>; };

        template<class TypeTag>
        struct LocalResidual<TypeTag, TTag::EclFlowProblemTPFA> { using type = BlackOilLocalResidualTPFA<TypeTag>; };

        template<class TypeTag>
        //struct ElementContext<TypeTag, TTag::EclFlowProblemTPFA> { using type = SmallElementContext<TypeTag>; };
        struct ElementContext<TypeTag, TTag::EclFlowProblemTPFA> { using type = FvBaseElementContext<TypeTag>; };

        template<class TypeTag>
        struct IntensiveQuantities<TypeTag, TTag::EclFlowProblemTPFA> {
            //using type = BlackOilIntensiveQuantitiesSimple<TypeTag>;
            using type = BlackOilIntensiveQuantities<TypeTag>;
        };

    }
}

int main(int argc, char** argv)
{
    using TypeTag = Opm::Properties::TTag::EclFlowProblemTPFA;
    auto mainObject = Opm::Main(argc, argv);
    return mainObject.runStatic<TypeTag>();
}
