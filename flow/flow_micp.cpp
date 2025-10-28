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

#include <flow/flow_micp.hpp>

#include <opm/grid/CpGrid.hpp>

#include <opm/material/common/ResetLocale.hpp>
#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>

#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/flow/SimulatorFullyImplicitBlackoil.hpp>

namespace Opm {
/*!
 * \brief Single-phase (water) flow problem including microbially induced calcite precipitation (MICP)
 * effects (three transported quantities: suspended microbes, oxygen, and urea, and two solid phases:
 * biofilm and calcite).
 */
namespace Properties {
namespace TTag {
struct FlowMICPProblem {
    using InheritsFrom = std::tuple<FlowProblem>;
};
}
template<class TypeTag>
struct EnableBioeffects<TypeTag, TTag::FlowMICPProblem> {
    static constexpr bool value = true;
};
template<class TypeTag>
struct Linearizer<TypeTag, TTag::FlowMICPProblem> { using type = TpfaLinearizer<TypeTag>; };
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::FlowMICPProblem> { using type = BlackOilLocalResidualTPFA<TypeTag>; };
//! The indices required by the model
template<class TypeTag>
struct Indices<TypeTag, TTag::FlowMICPProblem>
{
private:
    // it is unfortunately not possible to simply use 'TypeTag' here because this leads
    // to cyclic definitions of some properties. if this happens the compiler error
    // messages unfortunately are *really* confusing and not really helpful.
    using BaseTypeTag = TTag::FlowProblem;
    using FluidSystem = GetPropType<BaseTypeTag, Properties::FluidSystem>;

public:
    using type = BlackOilOnePhaseIndices<getPropValue<TypeTag, Properties::EnableSolvent>(),
                                         getPropValue<TypeTag, Properties::EnableExtbo>(),
                                         getPropValue<TypeTag, Properties::EnablePolymer>(),
                                         getPropValue<TypeTag, Properties::EnergyModuleType>() == EnergyModules::FullyImplicitThermal,
                                         getPropValue<TypeTag, Properties::EnableFoam>(),
                                         getPropValue<TypeTag, Properties::EnableBrine>(),
                                         /*PVOffset=*/0,
                                         /*enabledCompIdx=*/FluidSystem::waterCompIdx,
                                         5>; //Five MICP components
};

template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::FlowMICPProblem> { static constexpr bool value = true; };

template<class TypeTag>
struct EnableDispersion<TypeTag, TTag::FlowMICPProblem> { static constexpr bool value = true; };

}}

namespace Opm {

// ----------------- Main program -----------------
int flowMICPMain(int argc, char** argv, bool outputCout, bool outputFiles)
{
    // we always want to use the default locale, and thus spare us the trouble
    // with incorrect locale settings.
    resetLocale();

    FlowMain<Properties::TTag::FlowMICPProblem>
        mainfunc {argc, argv, outputCout, outputFiles};
    return mainfunc.execute();
}

int flowMICPMainStandalone(int argc, char** argv)
{
    using TypeTag = Properties::TTag::FlowMICPProblem;
    auto mainObject = std::make_unique<Opm::Main>(argc, argv);
    auto ret = mainObject->runStatic<TypeTag>();
    // Destruct mainObject as the destructor calls MPI_Finalize!
    mainObject.reset();
    return ret;
}

}
