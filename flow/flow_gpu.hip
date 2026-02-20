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


// For now flow_gpu is developed to support SPE11 simulations, that is 2-phase gas-water flow with thermal effects
// The goal is to support both property-evaluation and matrix assembly on the GPU, in addition to the linear
// solver which already works on the GPU. The CPU would still have to manage well contributions, although
// there are no wells in SPE11, only source-terms.

// Define making clear that the simulator supports AMG
#define FLOW_SUPPORT_AMG 1

#include <flow/flow_gpu.hpp>

#include <opm/common/utility/gpuDecorators.hpp>

#include <opm/material/common/ResetLocale.hpp>
#include <opm/material/thermal/EnergyModuleType.hpp>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/grid/CpGrid.hpp>
#include <opm/simulators/flow/SimulatorFullyImplicitBlackoil.hpp>
#include <opm/simulators/flow/SimplifiedGpuBlackOilModel.hpp>
#include <opm/simulators/flow/Main.hpp>

#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>

namespace Opm {

namespace Properties {
namespace TTag {
struct FlowGasWaterEnergyProblemGPU {
    using InheritsFrom = std::tuple<FlowProblem>;
};

template<template<class> class Storage>
struct FlowGasWaterEnergyProblemGPUTrue {
    using InheritsFrom = std::tuple<FlowGasWaterEnergyProblemGPU>;
};

template<template<class> class Storage>
struct to_gpu_type<FlowGasWaterEnergyProblemGPU, Storage> {
    using type = FlowGasWaterEnergyProblemGPUTrue<Storage>;
};
}

template<class TypeTag>
struct RunAssemblyOnGpu<TypeTag, TTag::FlowGasWaterEnergyProblemGPU> { static constexpr bool value = true; };

template<class TypeTag>
struct Scalar<TypeTag, TTag::FlowGasWaterEnergyProblemGPU> { using type = double; };

template<class TypeTag>
struct Linearizer<TypeTag, TTag::FlowGasWaterEnergyProblemGPU> { using type = TpfaLinearizer<TypeTag>; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::FlowGasWaterEnergyProblemGPU> { using type = BlackOilLocalResidualTPFA<TypeTag>; };

template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::FlowGasWaterEnergyProblemGPU> { static constexpr bool value = true; };

template<class TypeTag>
struct EnableDispersion<TypeTag, TTag::FlowGasWaterEnergyProblemGPU> { static constexpr bool value = false; };

template<class TypeTag>
struct EnableEnergy<TypeTag, TTag::FlowGasWaterEnergyProblemGPU> {
    static constexpr bool value = true;
};

// After a rebase the lines above were not sufficient to run a simulation with thermal effects
template<class TypeTag>
struct EnergyModuleType<TypeTag, TTag::FlowGasWaterEnergyProblemGPU>
{ static constexpr EnergyModules value = EnergyModules::FullyImplicitThermal; };

template<class TypeTag>
struct EnableDisgasInWater<TypeTag, TTag::FlowGasWaterEnergyProblemGPU> {
    static constexpr bool value = true;
};

template<class TypeTag>
struct EnableVapwat<TypeTag, TTag::FlowGasWaterEnergyProblemGPU> {
    static constexpr bool value = true;
};

template<class TypeTag>
struct GpuFIBlackOilModel<TypeTag, TTag::FlowGasWaterEnergyProblemGPU> { using type = SimplifiedGpuFIBlackOilModel<TypeTag>; };

template<class TypeTag, template<class> class Storage>
struct FluidSystem<TypeTag, TTag::FlowGasWaterEnergyProblemGPUTrue<Storage>>
{
    using type = Opm::BlackOilFluidSystemNonStatic<
        double,
        Opm::BlackOilDefaultFluidSystemIndices,
        Storage
    >;
};

//! The indices required by the model
template<class TypeTag>
struct Indices<TypeTag, TTag::FlowGasWaterEnergyProblemGPU>
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
                                         /*disabledCompIdx=*/FluidSystem::oilCompIdx,
                                         getPropValue<TypeTag, Properties::EnableBioeffects>()>;
};
}

// ----------------- Main program -----------------
int flowGasWaterEnergyMainGPU(int argc, char** argv, bool outputCout, bool outputFiles)
{
    // we always want to use the default locale, and thus spare us the trouble
    // with incorrect locale settings.
    resetLocale();

    FlowMain<Properties::TTag::FlowGasWaterEnergyProblemGPU>
        mainfunc {argc, argv, outputCout, outputFiles} ;
    return mainfunc.execute();
}

int flowGasWaterEnergyMainGPUStandalone(int argc, char** argv)
{
    using TypeTag = Properties::TTag::FlowGasWaterEnergyProblemGPU;
    auto mainObject = std::make_unique<::Opm::Main>(argc, argv);
    auto ret = mainObject->runStatic<TypeTag>();
    // Destruct mainObject as the destructor calls MPI_Finalize!
    mainObject.reset();
    return ret;
}

}
