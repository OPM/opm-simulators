/*
  Copyright 2025 Equinor ASA.

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

#include <flow/flow_blackoil.hpp>

#include <opm/material/common/ResetLocale.hpp>
#include <opm/grid/CpGrid.hpp>
#include <opm/simulators/flow/SimulatorFullyImplicitBlackoil.hpp>
#include <opm/simulators/flow/Main.hpp>

#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>


namespace Opm {
namespace Properties {
    namespace TTag {
        struct FlowProblemNohyst {
            using InheritsFrom = std::tuple<FlowProblem>;
        };
    }
    template<class TypeTag>
    struct Linearizer<TypeTag, TTag::FlowProblemNohyst> { using type = TpfaLinearizer<TypeTag>; };

    template<class TypeTag>
    struct LocalResidual<TypeTag, TTag::FlowProblemNohyst> { using type = BlackOilLocalResidualTPFA<TypeTag>; };

    template<class TypeTag>
    struct EnableDiffusion<TypeTag, TTag::FlowProblemNohyst> { static constexpr bool value = false; };

    template<class TypeTag>
    struct EnableHysteresis<TypeTag, TTag::FlowProblemNohyst> { static constexpr bool value = false; };

    template<class TypeTag>
    struct EnableEndpointScaling<TypeTag, TTag::FlowProblemNohyst> { static constexpr bool value = true; };

    template<class TypeTag>
    struct AvoidElementContext<TypeTag, TTag::FlowProblemNohyst> { static constexpr bool value = true; };

}

std::unique_ptr<FlowMain<Properties::TTag::FlowProblemNohyst>>
flowBlackoilTpfaNohystMainInit(int argc, char** argv, bool outputCout, bool outputFiles)
{
    // we always want to use the default locale, and thus spare us the trouble
    // with incorrect locale settings.
    resetLocale();

    return std::make_unique<FlowMain<Properties::TTag::FlowProblemNohyst>>(
        argc, argv, outputCout, outputFiles);
}

// ----------------- Main program -----------------
int flowBlackoilTpfaNohystMain(int argc, char** argv, bool outputCout, bool outputFiles)
{
    // we always want to use the default locale, and thus spare us the trouble
    // with incorrect locale settings.
    resetLocale();

    FlowMain<Properties::TTag::FlowProblemNohyst>
        mainfunc {argc, argv, outputCout, outputFiles};
    return mainfunc.execute();
}

int flowBlackoilTpfaNohystMainStandalone(int argc, char** argv)
{
    using TypeTag = Opm::Properties::TTag::FlowProblemNohyst;
    auto mainObject = std::make_unique<Opm::Main>(argc, argv);
    auto ret = mainObject->runStatic<TypeTag>();
    // Destruct mainObject as the destructor calls MPI_Finalize!
    mainObject.reset();
    return ret;
}

} // namespace Opm
