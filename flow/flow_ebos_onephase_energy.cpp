/*
  Copyright 2013, 2014, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015, 2017 IRIS AS

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

#include <flow/flow_ebos_onephase_energy.hpp>

#include <opm/simulators/flow/Main.hpp>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/material/common/ResetLocale.hpp>

namespace Opm {
namespace Properties {

namespace TTag {
struct FlowWaterOnlyEnergyProblem {
    using InheritsFrom = std::tuple<FlowProblem>;
};
}
template<class TypeTag>
struct EnableEnergy<TypeTag, TTag::FlowWaterOnlyEnergyProblem> {
    static constexpr bool value = true;
};
//! The indices required by the model
template<class TypeTag>
struct Indices<TypeTag, TTag::FlowWaterOnlyEnergyProblem>
{
private:
    // it is unfortunately not possible to simply use 'TypeTag' here because this leads
    // to cyclic definitions of some properties. if this happens the compiler error
    // messages unfortunately are *really* confusing and not really helpful.
    using BaseTypeTag = TTag::FlowProblem;
    using FluidSystem = GetPropType<BaseTypeTag, Properties::FluidSystem>;

public:
    using type = Opm::BlackOilOnePhaseIndices<getPropValue<TypeTag, Properties::EnableSolvent>(),
                                              getPropValue<TypeTag, Properties::EnableExtbo>(),
                                              getPropValue<TypeTag, Properties::EnablePolymer>(),
                                              getPropValue<TypeTag, Properties::EnableEnergy>(),
                                              getPropValue<TypeTag, Properties::EnableFoam>(),
                                              getPropValue<TypeTag, Properties::EnableBrine>(),
                                              /*PVOffset=*/0,
                                              /*enebledCompIdx=*/FluidSystem::waterCompIdx,
                                              getPropValue<TypeTag, Properties::EnableMICP>()>;
};

} // namespace Opm::Properties

// ----------------- Main program -----------------
int flowEbosWaterOnlyEnergyMain(int argc, char** argv, bool outputCout, bool outputFiles)
{
    // we always want to use the default locale, and thus spare us the trouble
    // with incorrect locale settings.
    resetLocale();

    FlowMain<Properties::TTag::FlowWaterOnlyEnergyProblem>
        mainfunc {argc, argv, outputCout, outputFiles};
    return mainfunc.execute();
}

int flowEbosWaterOnlyEnergyMainStandalone(int argc, char** argv)
{
    using TypeTag = Opm::Properties::TTag::FlowWaterOnlyEnergyProblem;
    auto mainObject = std::make_unique<Opm::Main>(argc, argv);
    auto ret = mainObject->runStatic<TypeTag>();
    // Destruct mainObject as the destructor calls MPI_Finalize!
    mainObject.reset();
    return ret;
}

}
