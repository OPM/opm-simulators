// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2025, NORCE AS

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
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>

#include <opm/simulators/flow/BlackoilModelProperties.hpp>
#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/flow/TTagFlowProblemTPSA.hpp>

#include <tuple>


namespace Opm::Properties {

namespace TTag {

struct FlowWaterOnlyProblemTPSA
{
    using InheritsFrom = std::tuple<FlowProblem, FlowProblemTpsa>;
};

}  // namespace Opm::Properties::TTag

// ///
// Flow related properties
// ///
template<class TypeTag>
struct Linearizer<TypeTag, TTag::FlowWaterOnlyProblemTPSA>
{ using type = TpfaLinearizer<TypeTag>; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::FlowWaterOnlyProblemTPSA>
{ using type = BlackOilLocalResidualTPFA<TypeTag>; };

template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::FlowWaterOnlyProblemTPSA>
{ static constexpr bool value = false; };

template<class TypeTag>
struct Indices<TypeTag, TTag::FlowWaterOnlyProblemTPSA>
{
private:
    // it is unfortunately not possible to simply use 'TypeTag' here because this leads
    // to cyclic definitions of some properties. if this happens the compiler error
    // messages unfortunately are *really* confusing and not really helpful.
    using BaseTypeTag = TTag::FlowProblem;
    using FluidSystem = GetPropType<BaseTypeTag, Properties::FluidSystem>;
    static constexpr EnergyModules energyModuleType = getPropValue<TypeTag, Properties::EnergyModuleType>();
    static constexpr int numEnergyVars = energyModuleType == EnergyModules::FullyImplicitThermal;
    static constexpr bool enableSeqImpEnergy = energyModuleType == EnergyModules::SequentialImplicitThermal;

public:
    using type = BlackOilOnePhaseIndices<getPropValue<TypeTag, Properties::EnableSolvent>(),
                                         getPropValue<TypeTag, Properties::EnableExtbo>(),
                                         getPropValue<TypeTag, Properties::EnablePolymer>(),
                                         numEnergyVars,
                                         enableSeqImpEnergy,
                                         getPropValue<TypeTag, Properties::EnableFoam>(),
                                         getPropValue<TypeTag, Properties::EnableBrine>(),
                                         /*PVOffset=*/0,
                                         /*enabledCompIdx=*/FluidSystem::waterCompIdx,
                                         getPropValue<TypeTag, Properties::EnableBioeffects>()>;
};

// ///
// TPSA related properties
// ///
template <class TypeTag>
struct EnableMech<TypeTag, TTag::FlowWaterOnlyProblemTPSA>
{ static constexpr bool value = true; };

template <class TypeTag>
struct Problem<TypeTag, TTag::FlowWaterOnlyProblemTPSA>
{ using type = FlowProblemTPSA<TypeTag>; };

template <class TypeTag>
struct NonlinearSystem<TypeTag, TTag::FlowWaterOnlyProblemTPSA>
{ using type = BlackoilModelTPSA<TypeTag>; };

}  // namespace Opm::Properties

namespace Opm {

// ----------------- Main program -----------------
int flowWaterOnlyTpsaMain(int argc, char** argv, bool outputCout, bool outputFiles)
{
    // we always want to use the default locale, and thus spare us the trouble
    // with incorrect locale settings.
    resetLocale();

    FlowMain<Properties::TTag::FlowWaterOnlyProblemTPSA>
        mainfunc {argc, argv, outputCout, outputFiles};
    return mainfunc.execute();
}

int flowWaterOnlyTpsaMainStandalone(int argc, char** argv)
{
    using TypeTag = Opm::Properties::TTag::FlowWaterOnlyProblemTPSA;
    auto mainObject = std::make_unique<Opm::Main>(argc, argv);
    auto ret = mainObject->runStatic<TypeTag>();
    // Destruct mainObject as the destructor calls MPI_Finalize!
    mainObject.reset();
    return ret;
}

}