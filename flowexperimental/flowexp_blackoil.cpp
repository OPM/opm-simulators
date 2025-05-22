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
#include "FlowExpNewtonMethod.hpp"
#include "flowexp.hpp"

#include <opm/models/discretization/common/tpfalinearizer.hh>
#include <opm/models/utils/parametersystem.hpp>
#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/flow/FlowProblemBlackoil.hpp>

#include "BlackOilIntensiveQuantitiesGlobalIndex.hpp"
#include "FIBlackOilModelNoCache.hpp"

namespace Opm::Properties {
namespace TTag {
struct FlowExpProblemBlackOil
{
    using InheritsFrom = std::tuple<FlowExpTypeTag>;
};
}

template<class TypeTag>
struct Model<TypeTag, TTag::FlowExpProblemBlackOil>
{
    using type = FIBlackOilModelNoCache<TypeTag>;
};

template<class TypeTag>
struct IntensiveQuantities<TypeTag, TTag::FlowExpProblemBlackOil>
{
    using type = BlackOilIntensiveQuantitiesGlobalIndex<TypeTag>;
};

// Set the problem class
template<class TypeTag>
struct Problem<TypeTag, TTag::FlowExpProblemBlackOil>
{
    using type = FlowExpProblem<TypeTag>;
};

template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::FlowExpProblemBlackOil>
{
    static constexpr bool value = false;
};

template<class TypeTag>
struct EnableDisgasInWater<TypeTag, TTag::FlowExpProblemBlackOil>
{
    static constexpr bool value = false;
};

template<class TypeTag>
struct Simulator<TypeTag, TTag::FlowExpProblemBlackOil>
{
    using type = Opm::Simulator<TypeTag>;
};

} // namespace Opm::Properties

int main(int argc, char** argv)
{
    using TypeTag = Opm::Properties::TTag::FlowExpProblemBlackOil;
    Opm::registerEclTimeSteppingParameters<double>();
    return Opm::start<TypeTag>(argc, argv, true);
}
