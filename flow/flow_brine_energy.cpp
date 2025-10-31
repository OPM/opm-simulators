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

#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>
#include <opm/simulators/flow/Main.hpp>
#include <opm/material/thermal/EnergyModuleType.hpp>

namespace Opm {
namespace Properties {
namespace TTag {
struct FlowBrineEnergyProblem {
    using InheritsFrom = std::tuple<FlowProblem>;
};
}

template<class TypeTag>
struct EnableBrine<TypeTag, TTag::FlowBrineEnergyProblem> {
    static constexpr bool value = true;
};

template<class TypeTag>
struct Linearizer<TypeTag, TTag::FlowBrineEnergyProblem> {
    using type = TpfaLinearizer<TypeTag>;
};

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::FlowBrineEnergyProblem> {
    using type = BlackOilLocalResidualTPFA<TypeTag>;
};
template<class TypeTag>
struct EnergyModuleType<TypeTag, TTag::FlowBrineEnergyProblem>
{ static constexpr EnergyModules value = EnergyModules::FullyImplicitThermal; };
}

int flowBrineEnergyMain(int argc, char** argv)
{
    using TypeTag = Opm::Properties::TTag::FlowBrineEnergyProblem;
    auto mainObject = std::make_unique<Opm::Main>(argc, argv);
    auto ret = mainObject->runStatic<TypeTag>();
    // Destruct mainObject as the destructor calls MPI_Finalize!
    mainObject.reset();
    return ret;
}

}
