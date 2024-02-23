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

#include <flow/flow_ebos_oilwater_polymer_injectivity.hpp>

#include <opm/material/common/ResetLocale.hpp>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/grid/CpGrid.hpp>
#include <opm/simulators/flow/SimulatorFullyImplicitBlackoil.hpp>
#include <opm/simulators/flow/Main.hpp>

namespace Opm {
namespace Properties {
namespace TTag {
struct FlowOilWaterPolymerInjectivityProblem {
    using InheritsFrom = std::tuple<FlowProblem>;
};
}
template<class TypeTag>
struct EnablePolymer<TypeTag, TTag::FlowOilWaterPolymerInjectivityProblem> {
    static constexpr bool value = true;
};
template<class TypeTag>
struct EnablePolymerMW<TypeTag, TTag::FlowOilWaterPolymerInjectivityProblem> {
    static constexpr bool value = true;
};

//! The indices required by the model
// For this case, there will be two primary variables introduced for the polymer
// polymer concentration and polymer molecular weight
template<class TypeTag>
struct Indices<TypeTag, TTag::FlowOilWaterPolymerInjectivityProblem>
{
private:
    // it is unfortunately not possible to simply use 'TypeTag' here because this leads
    // to cyclic definitions of some properties. if this happens the compiler error
    // messages unfortunately are *really* confusing and not really helpful.
    using BaseTypeTag = TTag::FlowProblem;
    using FluidSystem = GetPropType<BaseTypeTag, Properties::FluidSystem>;

public:
    using type = BlackOilTwoPhaseIndices<0,
                                         0,
                                         2,
                                         0,
                                         getPropValue<TypeTag, Properties::EnableFoam>(),
                                         getPropValue<TypeTag, Properties::EnableBrine>(),
                                         /*PVOffset=*/0,
                                         /*disabledCompIdx=*/FluidSystem::gasCompIdx,
                                         getPropValue<TypeTag, Properties::EnableMICP>()>;
};
}}

namespace Opm {
/* void flowEbosOilWaterPolymerInjectivitySetDeck(Deck& deck, EclipseState& eclState)
{
    using TypeTag = Properties::TTag::EclFlowOilWaterPolymerInjectivityProblem;
    using Vanguard = GetPropType<TypeTag, Properties::Vanguard>;

    Vanguard::setExternalDeck(std::move(deck, &eclState));
} */

// ----------------- Main program -----------------
int flowEbosOilWaterPolymerInjectivityMain(int argc, char** argv, bool outputCout, bool outputFiles)
{
    // we always want to use the default locale, and thus spare us the trouble
    // with incorrect locale settings.
    resetLocale();

    FlowMain<Properties::TTag::FlowOilWaterPolymerInjectivityProblem>
        mainfunc {argc, argv, outputCout, outputFiles};
    return mainfunc.execute();
}

int flowEbosOilWaterPolymerInjectivityMainStandalone(int argc, char** argv)
{
    using TypeTag = Properties::TTag::FlowOilWaterPolymerInjectivityProblem;
    auto mainObject = std::make_unique<Opm::Main>(argc, argv);
    auto ret = mainObject->runStatic<TypeTag>();
    // Destruct mainObject as the destructor calls MPI_Finalize!
    mainObject.reset();
    return ret;
}

}
