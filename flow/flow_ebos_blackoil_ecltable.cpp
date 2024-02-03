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

#include <flow/flow_ebos_blackoil.hpp>

#include <opm/material/common/ResetLocale.hpp>
#include <opm/grid/CpGrid.hpp>
#include <opm/simulators/flow/SimulatorFullyImplicitBlackoilEbos.hpp>
#include <opm/simulators/flow/Main.hpp>

#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>
#include <ebos/BlackoilModelFv.hpp>
#include <ebos/EclProblemStandard.hpp>
#include <ebos/eclblackoilintensivequantitiestable.hh>
namespace Opm {
    namespace Properties {

        template<class TypeTag>
        struct Linearizer<TypeTag, TTag::EclFlowProblemTABLE> { using type = TpfaLinearizer<TypeTag>; };

        template<class TypeTag>
        struct LocalResidual<TypeTag, TTag::EclFlowProblemTABLE> { using type = BlackOilLocalResidualTPFA<TypeTag>; };
        template<class TypeTag>
        struct EnableDiffusion<TypeTag, TTag::EclFlowProblemTABLE> { static constexpr bool value = false; };

        template<class TypeTag>
        struct Model<TypeTag, TTag::EclFlowProblemTABLE> {
            using type = BlackOilModelFv<TypeTag>;
        };

        template<class TypeTag>
        struct Problem<TypeTag, TTag::EclFlowProblemTABLE> {
            using type = EclProblemStandard<TypeTag>;
        };

        template<class TypeTag>
        struct IntensiveQuantities<TypeTag, TTag::EclFlowProblemTABLE> {
            using type = EclBlackOilIntensiveQuantitiesTABLE<TypeTag>;
        };
    }


}


namespace Opm
{
std::unique_ptr<FlowMainEbos<Properties::TTag::EclFlowProblemTABLE>>
flowEbosBlackoilECLTABLEMainInit(int argc, char** argv, bool outputCout, bool outputFiles)
{
    // we always want to use the default locale, and thus spare us the trouble
    // with incorrect locale settings.
    resetLocale();

    return std::make_unique<FlowMainEbos<Properties::TTag::EclFlowProblemTABLE>>(
        argc, argv, outputCout, outputFiles);
}

// ----------------- Main program -----------------
int flowEbosBlackoilECLTABLEMain(int argc, char** argv, bool outputCout, bool outputFiles)
{
    // we always want to use the default locale, and thus spare us the trouble
    // with incorrect locale settings.
    resetLocale();

    FlowMainEbos<Properties::TTag::EclFlowProblemTABLE>
        mainfunc {argc, argv, outputCout, outputFiles};
    return mainfunc.execute();
}

int flowEbosBlackoilECLTABLEMainStandalone(int argc, char** argv)
{
    using TypeTag = Properties::TTag::EclFlowProblemTABLE;
    auto mainObject = std::make_unique<Opm::Main>(argc, argv);
    auto ret = mainObject->runStatic<TypeTag>();
    // Destruct mainObject as the destructor calls MPI_Finalize!
    mainObject.reset();
    return ret;
}

} // namespace Opm
