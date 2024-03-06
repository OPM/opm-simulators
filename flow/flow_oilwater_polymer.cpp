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

#include <flow/flow_oilwater_polymer.hpp>

#include <opm/material/common/ResetLocale.hpp>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/grid/CpGrid.hpp>
#include <opm/simulators/flow/SimulatorFullyImplicitBlackoil.hpp>
#include <opm/simulators/flow/Main.hpp>

namespace Opm {
namespace Properties {
namespace TTag {
struct FlowOilWaterPolymerProblem {
    using InheritsFrom = std::tuple<FlowProblem>;
};
}
template<class TypeTag>
struct EnablePolymer<TypeTag, TTag::FlowOilWaterPolymerProblem> {
    static constexpr bool value = true;
};
//! The indices required by the model
template<class TypeTag>
struct Indices<TypeTag, TTag::FlowOilWaterPolymerProblem>
{
private:
    // it is unfortunately not possible to simply use 'TypeTag' here because this leads
    // to cyclic definitions of some properties. if this happens the compiler error
    // messages unfortunately are *really* confusing and not really helpful.
    using BaseTypeTag = TTag::FlowProblem;
    using FluidSystem = GetPropType<BaseTypeTag, Properties::FluidSystem>;

public:
    using type = BlackOilTwoPhaseIndices<getPropValue<TypeTag, Properties::EnableSolvent>(),
                                         getPropValue<TypeTag, Properties::EnableExtbo>(),
                                         getPropValue<TypeTag, Properties::EnablePolymer>(),
                                         getPropValue<TypeTag, Properties::EnableEnergy>(),
                                         getPropValue<TypeTag, Properties::EnableFoam>(),
                                         getPropValue<TypeTag, Properties::EnableBrine>(),
                                         /*PVOffset=*/0,
                                         /*disabledCompIdx=*/FluidSystem::gasCompIdx,
                                         getPropValue<TypeTag, Properties::EnableMICP>()>;
};
}}

namespace Opm {

// ----------------- Main program -----------------
int flowOilWaterPolymerMain(int argc, char** argv, bool outputCout, bool outputFiles)
{
    // we always want to use the default locale, and thus spare us the trouble
    // with incorrect locale settings.
    resetLocale();

    FlowMain<Properties::TTag::FlowOilWaterPolymerProblem>
        mainfunc {argc, argv, outputCout, outputFiles};
    return mainfunc.execute();
}

int flowOilWaterPolymerMainStandalone(int argc, char** argv)
{
    using TypeTag = Properties::TTag::FlowOilWaterPolymerProblem;
    auto mainObject = std::make_unique<Opm::Main>(argc, argv);
    auto ret = mainObject->runStatic<TypeTag>();
    // Destruct mainObject as the destructor calls MPI_Finalize!
    mainObject.reset();
    return ret;
}

}
