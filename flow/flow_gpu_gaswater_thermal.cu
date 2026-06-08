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

#undef HAVE_DUNE_ALUGRID // Disable ALUGRID if it is available as nvcc struggles with it
#undef HAVE_DUNE_FEM // Disable DUNE-FEM if it is available as nvcc struggles with it

// For now flow_gpu is developed to support SPE11 simulations, that is 2-phase gas-water flow with
// thermal effects The goal is to support both property-evaluation and matrix assembly on the GPU,
// in addition to the linear solver which already works on the GPU. The CPU would still have to
// manage well contributions, although there are no wells in SPE11, only source-terms.

#include <flow/flow_gpu.hpp>

#include <opm/common/utility/gpuDecorators.hpp>

#include <opm/material/common/ResetLocale.hpp>
#include <opm/material/thermal/EnergyModuleType.hpp>

#include <opm/models/blackoil/blackoilenergymodules.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/grid/CpGrid.hpp>
#include <opm/simulators/flow/FlowGasWaterEnergyTypeTag.hpp>
#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/flow/SimpleFIBlackOilModel.hpp>
#include <opm/simulators/flow/SimulatorFullyImplicit.hpp>

#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>

namespace Opm
{

int
flowGasWaterEnergyMainGPU(int argc, char** argv, bool outputCout, bool outputFiles)
{
    // we always want to use the default locale, and thus spare us the trouble
    // with incorrect locale settings.
    resetLocale();

    FlowMain<Properties::TTag::FlowGasWaterEnergyProblemGPU> mainfunc {
        argc, argv, outputCout, outputFiles};
    return mainfunc.execute();
}

int
flowGasWaterEnergyMainGPUStandalone(int argc, char** argv)
{
    using TypeTag = Properties::TTag::FlowGasWaterEnergyProblemGPU;
    auto mainObject = std::make_unique<::Opm::Main>(argc, argv);
    auto ret = mainObject->runStatic<TypeTag>();
    // Destruct mainObject as the destructor calls MPI_Finalize!
    mainObject.reset();
    return ret;
}

} // namespace Opm
