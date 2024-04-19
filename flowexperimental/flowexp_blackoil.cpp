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
#include <opm/simulators/flow/FlowProblem.hpp>
#include "eclnewtonmethod.hh"
#include "flowexp.hpp"
#include <opm/simulators/flow/Main.hpp>

#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>
#include "blackoilintensivequantitiessimple.hh"
#include "BlackOilModelFvNoCache.hpp"
// the current code use eclnewtonmethod adding other conditions to proceed_ should do the trick for KA
// adding linearshe sould be chaning the update_ function in the same class with condition that the error is reduced.
// the trick is to be able to recalculate the residual from here.
// unsure where the timestepping is done from suggestedtime??
// suggestTimeStep is taken from newton solver in problem.limitTimestep
/* namespace Opm{
    template<typename TypeTag>
    class OutputAuxModule : public BaseAuxiliaryModule<TypeTag>
    {

    };

} */


namespace Opm::Properties {
namespace TTag {
struct EclFlowProblemEbos {
    using InheritsFrom = std::tuple<EbosTypeTag>;
};
}

template<class TypeTag>
struct Model<TypeTag, TTag::EclFlowProblemEbos> {
    using type = BlackOilModelFvNoCache<TypeTag>;
};
template<class TypeTag>
struct IntensiveQuantities<TypeTag, TTag::EclFlowProblemEbos> {
     using type = BlackOilIntensiveQuantitiesSimple<TypeTag>;
};
// Set the problem class
template<class TypeTag>
struct Problem<TypeTag, TTag::EclFlowProblemEbos> {
    using type = EbosProblem<TypeTag>;
};


template<class TypeTag>
struct ThreadsPerProcess<TypeTag, TTag::EclFlowProblemEbos> {
    static constexpr int value = 1;
};

template<class TypeTag>
struct ContinueOnConvergenceError<TypeTag, TTag::EclFlowProblemEbos> {
    static constexpr bool value = false;
};

template<class TypeTag>
struct EclNewtonSumTolerance<TypeTag, TTag::EclFlowProblemEbos> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-5;
};

// the default for the allowed volumetric error for oil per second
template<class TypeTag>
struct NewtonTolerance<TypeTag, TTag::EclFlowProblemEbos> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-2;
};

// set fraction of the pore volume where the volumetric residual may be violated during
// strict Newton iterations
template<class TypeTag>
struct EclNewtonRelaxedVolumeFraction<TypeTag, TTag::EclFlowProblemEbos> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.0;
};

template<class TypeTag>
struct EclNewtonRelaxedTolerance<TypeTag, TTag::EclFlowProblemEbos> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 10*getPropValue<TypeTag, Properties::NewtonTolerance>();
};

//template<class TypeTag>
//struct Linearizer<TypeTag, TTag::EclFlowProblemEbos> { using type = TpfaLinearizer<TypeTag>; };

// template<class TypeTag>
// struct LocalResidual<TypeTag, TTag::EclFlowProblemEbos> { using type = BlackOilLocalResidualTPFA<TypeTag>; };

template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::EclFlowProblemEbos> { static constexpr bool value = false; };

template<class TypeTag>
struct EnableDisgasInWater<TypeTag, TTag::EclFlowProblemEbos> { static constexpr bool value = false; };

//static constexpr bool has_disgas_in_water = getPropValue<TypeTag, Properties::EnableDisgasInWater>();

template<class TypeTag>
struct Simulator<TypeTag, TTag::EclFlowProblemEbos> { using type = Opm::Simulator<TypeTag>; };

// template<class TypeTag>
// struct LinearSolverBackend<TypeTag, TTag::EclFlowProblemEbos> {
//     using type = ISTLSolver<TypeTag>;
// };

// // Set the problem class
// template<class TypeTag>
// struct Problem<TypeTag, TTag::EbosTypeTag> {
//     using type = EbosProblem<TypeTag>;
// };

// template<class TypeTag>
// struct EclEnableAquifers<TypeTag, TTag::EclFlowProblemTest> {
//     static constexpr bool value = false;
// };
}



int main(int argc, char** argv)
{
    using TypeTag = Opm::Properties::TTag::EclFlowProblemEbos;
    Opm::registerEclTimeSteppingParameters<TypeTag>();
    return Opm::start<TypeTag>(argc, argv);
}
