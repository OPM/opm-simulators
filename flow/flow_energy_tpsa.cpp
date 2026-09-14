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

#include <flow/flow_energy_tpsa.hpp>

#include <opm/models/blackoil/blackoilconvectivemixingmodule.hh>
#include <opm/models/blackoil/blackoildiffusionmodule.hh>
#include <opm/models/blackoil/blackoilenergymodules.hh>
#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>

#include <opm/simulators/flow/BlackoilModelProperties.hpp>
#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/flow/TTagFlowProblemTPSA.hpp>

#include <tuple>


namespace Opm::Properties {

namespace TTag {

struct FlowEnergyProblemTPSA
{
    using InheritsFrom = std::tuple<FlowProblem, FlowProblemTpsa>;
};

}  // namespace Opm::Properties::TTag

// ///
// Flow related properties
// ///
template<class TypeTag>
struct Linearizer<TypeTag, TTag::FlowEnergyProblemTPSA>
{ using type = TpfaLinearizer<TypeTag>; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::FlowEnergyProblemTPSA>
{ using type = BlackOilLocalResidualTPFA<TypeTag>; };

template<class TypeTag>
struct EnergyModuleType<TypeTag, TTag::FlowEnergyProblemTPSA>
{ static constexpr EnergyModules value = EnergyModules::FullyImplicitThermal; };

// ///
// TPSA related properties
// ///
template <class TypeTag>
struct EnableMech<TypeTag, TTag::FlowEnergyProblemTPSA>
{ static constexpr bool value = true; };

template <class TypeTag>
struct Problem<TypeTag, TTag::FlowEnergyProblemTPSA>
{ using type = FlowProblemTPSA<TypeTag>; };

template <class TypeTag>
struct NonlinearSystem<TypeTag, TTag::FlowEnergyProblemTPSA>
{ using type = NonlinearSystemBlackOilReservoirTPSA<TypeTag>; };

}  // namespace Opm::Properties

namespace Opm {

// ----------------- Main program -----------------
int flowEnergyTpsaMain(int argc, char** argv, bool outputCout, bool outputFiles)
{
    // we always want to use the default locale, and thus spare us the trouble
    // with incorrect locale settings.
    resetLocale();

    FlowMain<Properties::TTag::FlowEnergyProblemTPSA>
        mainfunc {argc, argv, outputCout, outputFiles};
    return mainfunc.execute();
}

int flowEnergyTpsaMainStandalone(int argc, char** argv)
{
    using TypeTag = Properties::TTag::FlowEnergyProblemTPSA;
    auto mainObject = std::make_unique<Opm::Main>(argc, argv);
    auto ret = mainObject->runStatic<TypeTag>();
    // Destruct mainObject as the destructor calls MPI_Finalize!
    mainObject.reset();
    return ret;
}

}  // namespace Opm
