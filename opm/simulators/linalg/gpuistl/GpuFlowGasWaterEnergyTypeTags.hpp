// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
 * Shared GPU TypeTag derivations on top of the simulation-side
 * \c FlowGasWaterEnergyProblem TypeTag (CO2STORE configuration).
 *
 * The intent is that GPU code (the dispatcher and the
 * \c test_blackoilintensivequantities_gpu test) reuses the *exact* same
 * TypeTag as the simulator and only overrides the small set of properties
 * that have to change for the GPU build (storage container of the primary
 * variables, the material law manager, the fluid system, the problem and
 * the element context).
 */
#ifndef OPM_GPU_FLOW_GASWATER_ENERGY_TYPETAGS_HPP
#define OPM_GPU_FLOW_GASWATER_ENERGY_TYPETAGS_HPP

#if HAVE_CUDA

#include <opm/material/fluidmatrixinteractions/EclMaterialLawTwoPhaseTypes.hpp>
#include <opm/material/fluidmatrixinteractions/EclTwoPhaseMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/PiecewiseLinearTwoPhaseMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/PiecewiseLinearTwoPhaseMaterialParams.hpp>

#include <opm/material/fluidsystems/BlackOilFluidSystemNonStatic.hpp>

#include <opm/models/blackoil/blackoilprimaryvariables.hh>
#include <opm/models/discretization/common/fvbaseelementcontextgpu.hh>

#include <opm/material/thermal/EclSpecrockLaw.hpp>
#include <opm/material/thermal/EclSpecrockLawParams.hpp>
#include <opm/material/thermal/EclThconrLaw.hpp>
#include <opm/material/thermal/EclThconrLawParams.hpp>

#include <opm/simulators/flow/FlowGasWaterEnergyTypeTag.hpp>
#include <opm/simulators/flow/GpuEclMaterialLawManager.hpp>
#include <opm/simulators/flow/GpuEclThermalLawManager.hpp>
#include <opm/simulators/flow/GpuFlowProblem.hpp>
#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>
#include <opm/simulators/linalg/gpuistl/MiniVector.hpp>
#include <opm/simulators/linalg/gpuistl/gpu_smart_pointer.hpp>

#include <tuple>

namespace Opm::Properties {

namespace TTag {

/// The GPU sibling of \c FlowGasWaterEnergyProblem. Inherits all property
/// bindings from the simulation TypeTag and only overrides the storage of
/// the primary variables, the material law and the fluid system to their
/// GPU equivalents.
struct FlowGasWaterEnergyProblemGPU {
    using InheritsFrom = std::tuple<FlowGasWaterEnergyProblem>;
};

/// A second derived TypeTag used for the per-cell intensive-quantities
/// kernel: it additionally swaps in \c GpuFlowProblem and
/// \c FvBaseElementContextGpu so that the GPU device side has trivially
/// copyable \c Problem / \c ElementContext types.
struct FlowGasWaterEnergyDummyProblemGPU {
    using InheritsFrom = std::tuple<FlowGasWaterEnergyProblemGPU>;
};

} // namespace TTag

// -----------------------------------------------------------------------
// Diffusion / dispersion: not yet implemented in the GPU IntensiveQuantities
// update path; force them off on the GPU TypeTag.
// -----------------------------------------------------------------------
template <class TypeTag>
struct EnableDiffusion<TypeTag, TTag::FlowGasWaterEnergyProblemGPU>
{
    static constexpr bool value = false;
};

template <class TypeTag>
struct EnableDispersion<TypeTag, TTag::FlowGasWaterEnergyProblemGPU>
{
    static constexpr bool value = false;
};

// -----------------------------------------------------------------------
// MaterialLaw: GPU-friendly EclTwoPhaseMaterial backed by GpuView storage.
// -----------------------------------------------------------------------
template <class TypeTag>
struct MaterialLaw<TypeTag, TTag::FlowGasWaterEnergyProblemGPU>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TTag::FlowProblem, Properties::FluidSystem>;
    using Traits = ThreePhaseMaterialTraits<Scalar,
                                            /*wettingPhaseIdx=*/FluidSystem::waterPhaseIdx,
                                            /*nonWettingPhaseIdx=*/FluidSystem::oilPhaseIdx,
                                            /*gasPhaseIdx=*/FluidSystem::gasPhaseIdx,
                                            getPropValue<TypeTag, Properties::EnableHysteresis>(),
                                            getPropValue<TypeTag, Properties::EnableEndpointScaling>()>;
    using TwoPhaseTraits = TwoPhaseMaterialTraits<Scalar,
                                                  /*wettingPhaseIdx=*/Traits::wettingPhaseIdx,
                                                  /*nonWettingPhaseIdx=*/Traits::nonWettingPhaseIdx>;
    using TwoPhaseParams = ::Opm::PiecewiseLinearTwoPhaseMaterialParams<
        TwoPhaseTraits, ::Opm::gpuistl::GpuView<const Scalar>>;
    using TwoPhaseLaw = ::Opm::PiecewiseLinearTwoPhaseMaterial<TwoPhaseTraits, TwoPhaseParams>;
    using GpuMaterialLawParams = ::Opm::EclTwoPhaseMaterialParams<
        Traits, TwoPhaseParams, TwoPhaseParams, TwoPhaseParams, ::Opm::gpuistl::ValueAsPointer>;
    using GpuMaterialLaw = ::Opm::EclTwoPhaseMaterial<
        Traits, TwoPhaseLaw, TwoPhaseLaw, TwoPhaseLaw, GpuMaterialLawParams>;

public:
    using EclMaterialLawManager = ::Opm::EclMaterialLaw::GpuManager<
        Traits, TwoPhaseLaw, TwoPhaseLaw, ::Opm::VectorWithDefaultAllocator, GpuMaterialLaw>;
    using type = typename EclMaterialLawManager::MaterialLaw;
};

// -----------------------------------------------------------------------
// PrimaryVariables: GPU build uses MiniVector storage so the type is
// trivially copyable to the device.
// -----------------------------------------------------------------------
template <class TypeTag>
struct PrimaryVariables<TypeTag, TTag::FlowGasWaterEnergyProblemGPU>
{
    using type = ::Opm::BlackOilPrimaryVariables<TypeTag, ::Opm::gpuistl::MiniVector>;
};

// -----------------------------------------------------------------------
// FluidSystem: non-static, GpuView-backed instantiation.
// -----------------------------------------------------------------------
template <class TypeTag>
struct FluidSystem<TypeTag, TTag::FlowGasWaterEnergyProblemGPU>
{
    using type = ::Opm::BlackOilFluidSystemNonStatic<
        GetPropType<TTag::FlowProblem, Properties::Scalar>,
        ::Opm::BlackOilDefaultFluidSystemIndices,
        ::Opm::gpuistl::GpuView>;
};

template <class TypeTag>
struct FluidSystem<TypeTag, TTag::FlowGasWaterEnergyDummyProblemGPU>
{
    using type = ::Opm::BlackOilFluidSystemNonStatic<
        GetPropType<TTag::FlowProblem, Properties::Scalar>,
        ::Opm::BlackOilDefaultFluidSystemIndices,
        ::Opm::gpuistl::GpuView>;
};

// -----------------------------------------------------------------------
// SolidEnergyLaw / ThermalConductionLaw: bind to the GPU-portable laws
// (EclSpecrockLaw / EclThconrLaw). Both expose an
// \c EclThermalLawManager nested typedef matching the convention used by
// FlowBaseProblemProperties so existing GetProp<...>::EclThermalLawManager
// chains keep working.
// -----------------------------------------------------------------------
template <class TypeTag>
struct SolidEnergyLaw<TypeTag, TTag::FlowGasWaterEnergyProblemGPU>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

public:
    using EclThermalLawManager = ::Opm::EclThermalLaw::GpuManager<
        Scalar, FluidSystem, ::Opm::gpuistl::GpuView, ::Opm::gpuistl::GpuView>;

    using type = ::Opm::EclSpecrockLaw<
        Scalar,
        ::Opm::EclSpecrockLawParams<Scalar, ::Opm::gpuistl::GpuView>>;
};

template <class TypeTag>
struct ThermalConductionLaw<TypeTag, TTag::FlowGasWaterEnergyProblemGPU>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

public:
    using EclThermalLawManager = ::Opm::EclThermalLaw::GpuManager<
        Scalar, FluidSystem, ::Opm::gpuistl::GpuView, ::Opm::gpuistl::GpuView>;

    using type = ::Opm::EclThconrLaw<Scalar, FluidSystem>;
};

// -----------------------------------------------------------------------
// Problem: GpuFlowProblem (GpuView-backed) instead of FlowProblemBlackoil.
// -----------------------------------------------------------------------
template <class TypeTag>
struct Problem<TypeTag, TTag::FlowGasWaterEnergyDummyProblemGPU>
{
private:
    using ScalarT = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using CpuMaterialLawManager =
        typename Opm::GetProp<TypeTag, Opm::Properties::MaterialLaw>::EclMaterialLawManager;
    using GpuViewMaterialLawManager =
        ::Opm::EclMaterialLaw::GpuManager<typename CpuMaterialLawManager::Traits,
                                          typename CpuMaterialLawManager::GasOilLaw,
                                          typename CpuMaterialLawManager::OilWaterLaw,
                                          ::Opm::gpuistl::GpuView,
                                          typename CpuMaterialLawManager::MaterialLaw>;
    using GpuViewThermalLawManager = ::Opm::EclThermalLaw::GpuManager<
        ScalarT, FluidSystem, ::Opm::gpuistl::GpuView, ::Opm::gpuistl::GpuView>;

public:
    using type = ::Opm::GpuFlowProblem<ScalarT,
                                       GpuViewMaterialLawManager,
                                       ::Opm::gpuistl::GpuView,
                                       GpuViewThermalLawManager>;
};

// -----------------------------------------------------------------------
// ElementContext: trivially copyable GPU element context.
// -----------------------------------------------------------------------
template <class TypeTag>
struct ElementContext<TypeTag, TTag::FlowGasWaterEnergyDummyProblemGPU>
{
    using type = ::Opm::FvBaseElementContextGpu<TypeTag>;
};

} // namespace Opm::Properties

#endif // HAVE_CUDA

#endif // OPM_GPU_FLOW_GASWATER_ENERGY_TYPETAGS_HPP
