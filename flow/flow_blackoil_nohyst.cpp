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

#include <flow/flow_blackoil.hpp>

#include <opm/material/common/ResetLocale.hpp>
#include <opm/grid/CpGrid.hpp>
#include <opm/simulators/flow/SimulatorFullyImplicitBlackoil.hpp>
#include <opm/simulators/flow/Main.hpp>

#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>


namespace Opm {
namespace Properties {
    namespace TTag {
        struct FlowProblemNoHystNoEps {
            using InheritsFrom = std::tuple<FlowProblem>;
        };
    }
    template<class TypeTag>
    struct Linearizer<TypeTag, TTag::FlowProblemNoHystNoEps> { using type = TpfaLinearizer<TypeTag>; };

    template<class TypeTag>
    struct LocalResidual<TypeTag, TTag::FlowProblemNoHystNoEps> { using type = BlackOilLocalResidualTPFA<TypeTag>; };

    template<class TypeTag>
    struct EnableDiffusion<TypeTag, TTag::FlowProblemNoHystNoEps> { static constexpr bool value = false; };

    template<class TypeTag>
    struct EnableHysteresis<TypeTag, TTag::FlowProblemNoHystNoEps> { static constexpr bool value = false; };

    template<class TypeTag>
    struct EnableEndpointScaling<TypeTag, TTag::FlowProblemNoHystNoEps> { static constexpr bool value = true; };

    template<class TypeTag>
    struct AvoidElementContext<TypeTag, TTag::FlowProblemNoHystNoEps> { static constexpr bool value = true; };

}
} // namespace Opm


int main(int argc, char** argv)
{
    using TypeTag = Opm::Properties::TTag::FlowProblemNoHystNoEps;
    auto mainObject = std::make_unique<Opm::Main>(argc, argv);
    auto ret = mainObject->runStatic<TypeTag>();
    // Destruct mainObject as the destructor calls MPI_Finalize!
    mainObject.reset();
    return ret;
}
