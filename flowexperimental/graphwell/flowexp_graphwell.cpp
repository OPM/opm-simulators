/*
  Copyright 2026 SINTEF Digital, Mathematics and Cybernetics.

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

#include <flowexperimental/FlowExpNewtonMethod.hpp>
#include <flowexperimental/flowexp.hpp>

#include <opm/models/discretization/common/tpfalinearizer.hh>
#include <opm/models/utils/parametersystem.hpp>
#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/flow/FlowProblemBlackoil.hpp>

#include <flowexperimental/BlackOilIntensiveQuantitiesGlobalIndex.hpp>
#include <flowexperimental/FIBlackOilModelNoCache.hpp>
#include <flowexperimental/graphwell/FlowExpGraphWellModel.hpp>

namespace Opm::Properties {
namespace TTag {
struct FlowExpGraphWellProblem
{
    using InheritsFrom = std::tuple<FlowExpTypeTag>;
};
}

template<class TypeTag>
struct Model<TypeTag, TTag::FlowExpGraphWellProblem>
{
    using type = FIBlackOilModelNoCache<TypeTag>;
};

template<class TypeTag>
struct IntensiveQuantities<TypeTag, TTag::FlowExpGraphWellProblem>
{
    using type = BlackOilIntensiveQuantitiesGlobalIndex<TypeTag>;
};

template<class TypeTag>
struct Problem<TypeTag, TTag::FlowExpGraphWellProblem>
{
    using type = FlowExpProblem<TypeTag>;
};

// The GraphWell-based well model.
template<class TypeTag>
struct WellModel<TypeTag, TTag::FlowExpGraphWellProblem>
{
    using type = FlowExpGraphWellModel<TypeTag>;
};

template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::FlowExpGraphWellProblem>
{
    static constexpr bool value = false;
};

template<class TypeTag>
struct EnableDisgasInWater<TypeTag, TTag::FlowExpGraphWellProblem>
{
    static constexpr bool value = false;
};

template<class TypeTag>
struct Simulator<TypeTag, TTag::FlowExpGraphWellProblem>
{
    using type = Opm::Simulator<TypeTag>;
};

} // namespace Opm::Properties

int main(int argc, char** argv)
{
    using TypeTag = Opm::Properties::TTag::FlowExpGraphWellProblem;
    Opm::registerEclTimeSteppingParameters<double>();
    return Opm::start<TypeTag>(argc, argv, true);
}
