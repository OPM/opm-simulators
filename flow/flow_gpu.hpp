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
#ifndef FLOW_GPU_HIP_HPP
#define FLOW_GPU_HIP_HPP

#include <opm/simulators/flow/FlowGasWaterEnergyTypeTag.hpp>

namespace Opm
{

namespace Properties
{
    namespace TTag
    {
        // FlowGasWaterEnergyProblemGPU is declared in FlowGasWaterEnergyTypeTag.hpp
        // (InheritsFrom = FlowGasWaterEnergyProblem).  The template below maps it
        // to FlowGasWaterEnergyProblemGPUTrue<Storage> for GPU-storage variants.

        template <template <class> class Storage>
        struct FlowGasWaterEnergyProblemGPUTrue {
            using InheritsFrom = std::tuple<FlowGasWaterEnergyProblemGPU>;
        };

        template <template <class> class Storage>
        struct to_gpu_type<FlowGasWaterEnergyProblemGPU, Storage> {
            using type = FlowGasWaterEnergyProblemGPUTrue<Storage>;
        };
    } // namespace TTag

    template <class TypeTag>
    struct RunAssemblyOnGpu<TypeTag, TTag::FlowGasWaterEnergyProblemGPU> {
        static constexpr bool value = true;
    };

    template <class TypeTag>
    struct Scalar<TypeTag, TTag::FlowGasWaterEnergyProblemGPU> {
        using type = double;
    };

    template <class TypeTag>
    struct GpuFIBlackOilModel<TypeTag, TTag::FlowGasWaterEnergyProblemGPU> {
        using type = SimplifiedGpuFIBlackOilModel<TypeTag>;
    };

    template <class TypeTag, template <class> class Storage>
    struct FluidSystem<TypeTag, TTag::FlowGasWaterEnergyProblemGPUTrue<Storage>> {
        using type = Opm::
            BlackOilFluidSystemNonStatic<double, Opm::BlackOilDefaultFluidSystemIndices, Storage>;
    };
} // namespace Properties

//! \brief Main function used in flow binary.
int flowGasWaterEnergyMainGPU(int argc, char** argv, bool outputCout, bool outputFiles);

//! \brief Main function used in flow_gaswater binary.
int flowGasWaterEnergyMainGPUStandalone(int argc, char** argv);

} // namespace Opm

#endif // FLOW_GPU_HIP_HPP
