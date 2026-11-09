/*
  Copyright 2026 Equinor ASA

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
/*!
 * \file
 *
 * GPU intensive-quantities dispatcher for the CO2STORE configuration
 * (water + gas + thermal, \c FlowGasWaterEnergyProblem).
 */
#include <config.h>

// Tricks for older versions of ROCm that do not support alugrid or dune-fem
// (mirrors test_blackoilintensivequantities_gpu.cu).
#ifdef HAVE_DUNE_ALUGRID
#undef HAVE_DUNE_ALUGRID
#endif
#define HAVE_DUNE_ALUGRID 0
#ifdef HAVE_DUNE_FEM
#undef HAVE_DUNE_FEM
#endif
#define HAVE_DUNE_FEM 0

#include <opm/simulators/linalg/gpuistl/GpuBlackoilIntensiveQuantitiesDispatcher_impl.hpp>

#include <opm/material/fluidmatrixinteractions/EclTwoPhaseMaterialParams.hpp>
#include <opm/material/thermal/EclSolidEnergyLawMultiplexerParams.hpp>
#include <opm/material/thermal/EclThermalConductionLawMultiplexerParams.hpp>

#include <opm/simulators/flow/FlowGasWaterEnergyTypeTag.hpp>
#include <opm/simulators/linalg/gpuistl/GpuFlowGasWaterEnergyTypeTags.hpp>

namespace Opm::gpuistl {

namespace dispatcher_detail {

struct GasWaterEnergyDispatcherTraits
{
    using CpuTag = Properties::TTag::FlowGasWaterEnergyCpuKernelBase;
    using GpuTag = Properties::TTag::FlowGasWaterEnergyDummyProblemGPU;

    using FluidSystemView = GetPropType<GpuTag, Properties::FluidSystem>;
    using GpuPrimaryVariables = GetPropType<GpuTag, Properties::PrimaryVariables>;
    using GpuIntensiveQuantities = BlackOilIntensiveQuantities<GpuTag>;
    using GpuFlowProblemBuf = dispatcher_detail::GpuFlowProblemBufferFor<GpuTag>;

    static constexpr const char* name = "CO2STORE: gas-water with energy";

    template <class Problem>
    static void validate(const Problem& problem)
    {
        dispatcher_detail::assertNoDiffusionOrDispersion<GpuTag>();
        dispatcher_detail::rejectGeoModifier(problem);

        const auto& materialLawManager = problem.materialLawManager();
        if (materialLawManager.twoPhaseApproach() != EclTwoPhaseApproach::GasWater) {
            OPM_THROW(std::logic_error,
                      "GPU property evaluation requires the gas-water two-phase approach");
        }
        if (materialLawManager.hasOil() || !materialLawManager.hasGas()
            || !materialLawManager.hasWater()) {
            OPM_THROW(std::logic_error,
                      "GPU property evaluation requires water and gas without an active oil phase");
        }
        dispatcher_detail::rejectUnsupportedSaturationFunctions(materialLawManager);

        const auto thermalLawManager = problem.thermalLawManager();
        if (thermalLawManager == nullptr) {
            OPM_THROW(std::invalid_argument,
                      "GPU property evaluation requires a thermal-law manager");
        }
        if (thermalLawManager->solidEnergyApproach() != EclSolidEnergyApproach::Specrock
            || thermalLawManager->thermalConductionApproach()
                != EclThermalConductionApproach::Thconr) {
            OPM_THROW(std::logic_error,
                      "GPU property evaluation requires SPECROCK solid energy and THCONR "
                      "thermal conduction");
        }
    }
};

} // namespace dispatcher_detail

template <>
struct GpuDispatcherTraits<Properties::TTag::FlowGasWaterEnergyProblem>
    : dispatcher_detail::GasWaterEnergyDispatcherTraits
{
};

// The GPU-assembly TypeTag of flow_gpu inherits all physics from
// FlowGasWaterEnergyProblem, so it has the same CPU-side Problem,
// PrimaryVariables and IntensiveQuantities, and the same GPU types.
template <>
struct GpuDispatcherTraits<Properties::TTag::FlowGasWaterEnergyProblemGPU>
    : dispatcher_detail::GasWaterEnergyDispatcherTraits
{
};

// Used by the Flow simulator binary when a CO2STORE deck with
// WATER+GAS+THERMAL is loaded.
template class GpuBlackoilIntensiveQuantitiesDispatcher<
    Properties::TTag::FlowGasWaterEnergyProblem>;

// Used by the flow_gpu binary, which runs both the assembly and the
// intensive quantities on the GPU.
template class GpuBlackoilIntensiveQuantitiesDispatcher<
    Properties::TTag::FlowGasWaterEnergyProblemGPU>;

} // namespace Opm::gpuistl
