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

// Define making clear that the simulator supports AMG
#define FLOW_SUPPORT_AMG 1

#include <flow/flow_ebos_gaswater_dissolution_diffuse.hpp>

#include <opm/material/common/ResetLocale.hpp>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/grid/CpGrid.hpp>
#include <opm/simulators/flow/SimulatorFullyImplicitBlackoilEbos.hpp>
#include <opm/simulators/flow/Main.hpp>

#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>

namespace Opm {
namespace Properties {
namespace TTag {
struct EclFlowGasWaterDissolutionDiffuseProblem {
    using InheritsFrom = std::tuple<EclFlowProblem>;
};
}

template<class TypeTag>
struct Linearizer<TypeTag, TTag::EclFlowGasWaterDissolutionDiffuseProblem> { using type = TpfaLinearizer<TypeTag>; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::EclFlowGasWaterDissolutionDiffuseProblem> { using type = BlackOilLocalResidualTPFA<TypeTag>; };

template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::EclFlowGasWaterDissolutionDiffuseProblem> { static constexpr bool value = true; };

template<class TypeTag>
struct EnableDisgasInWater<TypeTag, TTag::EclFlowGasWaterDissolutionDiffuseProblem> {
    static constexpr bool value = true;
};

template<class TypeTag>
struct EnableVapwat<TypeTag, TTag::EclFlowGasWaterDissolutionDiffuseProblem> {
    static constexpr bool value = true;
};

//! The indices required by the model
template<class TypeTag>
struct Indices<TypeTag, TTag::EclFlowGasWaterDissolutionDiffuseProblem>
{
private:
    // it is unfortunately not possible to simply use 'TypeTag' here because this leads
    // to cyclic definitions of some properties. if this happens the compiler error
    // messages unfortunately are *really* confusing and not really helpful.
    using BaseTypeTag = TTag::EclFlowProblem;
    using FluidSystem = GetPropType<BaseTypeTag, Properties::FluidSystem>;

public:
    typedef BlackOilTwoPhaseIndices<getPropValue<TypeTag, Properties::EnableSolvent>(),
                                    getPropValue<TypeTag, Properties::EnableExtbo>(),
                                    getPropValue<TypeTag, Properties::EnablePolymer>(),
                                    getPropValue<TypeTag, Properties::EnableEnergy>(),
                                    getPropValue<TypeTag, Properties::EnableFoam>(),
                                    getPropValue<TypeTag, Properties::EnableBrine>(),
                                    /*PVOffset=*/0,
                                    /*disabledCompIdx=*/FluidSystem::oilCompIdx,
                                    getPropValue<TypeTag, Properties::EnableMICP>()> type;
};
}}

namespace Opm {


// ----------------- Main program -----------------
int flowEbosGasWaterDissolutionDiffuseMain(int argc, char** argv, bool outputCout, bool outputFiles)
{
    // we always want to use the default locale, and thus spare us the trouble
    // with incorrect locale settings.
    resetLocale();

    FlowMainEbos<Properties::TTag::EclFlowGasWaterDissolutionDiffuseProblem>
        mainfunc {argc, argv, outputCout, outputFiles} ;
    return mainfunc.execute();
}

int flowEbosGasWaterDissolutionDiffuseMainStandalone(int argc, char** argv)
{
    using TypeTag = Properties::TTag::EclFlowGasWaterDissolutionDiffuseProblem;
    auto mainObject = std::make_unique<Opm::Main>(argc, argv);
    auto ret = mainObject->runStatic<TypeTag>();
    // Destruct mainObject as the destructor calls MPI_Finalize!
    mainObject.reset();
    return ret;
}

}
