/*
  Copyright 2020, NORCE AS

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
#include <opm/simulators/flow/Main.hpp>
#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>
#include "tracy/Tracy.hpp"
#include "tracy/TracyC.h"
#include <opm/models/blackoil/blackoilintensivequantitiessimple.hh>
//sudo apt install libcapstone-dev
#define OPM_TIME_BLOCK(blockname) ZoneNamedN(blockname, #blockname, true);
namespace Opm {
namespace Properties {
namespace TTag {
struct EclFlowProblemTPFAFast {
    using InheritsFrom = std::tuple<EclFlowProblem>;
};
}
    template<class TypeTag>
    struct Linearizer<TypeTag, TTag::EclFlowProblemTPFAFast> { using type = TpfaLinearizer<TypeTag>; };
    
    template<class TypeTag>
    struct LocalResidual<TypeTag, TTag::EclFlowProblemTPFAFast> { using type = BlackOilLocalResidualTPFA<TypeTag>; };
    
    template<class TypeTag>
    struct EnableDiffusion<TypeTag, TTag::EclFlowProblemTPFAFast> { static constexpr bool value = false; };

    template<class TypeTag>
    struct IntensiveQuantities<TypeTag, TTag::EclFlowProblemTPFAFast> {
        using type = BlackOilIntensiveQuantitiesSimple<TypeTag>;
    };

}
}

int main(int argc, char** argv)
{
    OPM_TIME_BLOCK(FullSimulator);
    using TypeTag = Opm::Properties::TTag::EclFlowProblemTPFAFast;
    auto mainObject = Opm::Main(argc, argv);
    return mainObject.runStatic<TypeTag>();
}
