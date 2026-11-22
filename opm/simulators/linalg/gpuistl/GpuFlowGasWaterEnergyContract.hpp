// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2026 Equinor ASA

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
*/
/*!
 * \file
 *
 * \brief Dependency-light GPU TypeTag contract for the CO2STORE flow path.
 *
 * The property-evaluation and TPFA-assembly paths must instantiate exactly the
 * same device intensive-quantity type.  This header deliberately does not
 * include the executable TypeTag header or FIBlackOilModel.hpp: it is included
 * by both paths while those headers are still being assembled.
 */
#ifndef OPM_GPU_FLOW_GASWATER_ENERGY_CONTRACT_HPP
#define OPM_GPU_FLOW_GASWATER_ENERGY_CONTRACT_HPP

#if HAVE_CUDA

#include <opm/material/fluidmatrixinteractions/EclMaterialLawTwoPhaseTypes.hpp>
#include <opm/material/fluidmatrixinteractions/EclTwoPhaseMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/GpuEclMaterialLawManager.hpp>
#include <opm/material/fluidmatrixinteractions/PiecewiseLinearTwoPhaseMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/PiecewiseLinearTwoPhaseMaterialParams.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystemNonStatic.hpp>
#include <opm/material/thermal/EclSpecrockLaw.hpp>
#include <opm/material/thermal/EclSpecrockLawParams.hpp>
#include <opm/material/thermal/EclThconrLaw.hpp>
#include <opm/material/thermal/EclThconrLawParams.hpp>
#include <opm/material/thermal/GpuEclThermalLawManager.hpp>

#include <opm/models/blackoil/blackoilprimaryvariables.hh>
#include <opm/models/blackoil/blackoilproperties.hh>
#include <opm/models/common/multiphasebaseproperties.hh>
#include <opm/models/utils/propertysystem.hh>

#include <opm/simulators/flow/GpuFlowProblem.hpp>
#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>
#include <opm/simulators/linalg/gpuistl/MiniVector.hpp>
#include <opm/simulators/linalg/gpuistl/gpu_smart_pointer.hpp>

#include <tuple>

namespace Opm::Properties {

// These Flow properties are normally declared by
// FlowBaseProblemProperties.hpp. Forward-declaring them here avoids pulling
// that header (and its tpfalinearizer dependency) into the contract cycle.
template <class TypeTag, class MyTypeTag>
struct EnableHysteresis;
template <class TypeTag, class MyTypeTag>
struct EnableEndpointScaling;

namespace TTag {

// The definitions live in FlowGasWaterEnergyTypeTag.hpp.  Forward declarations
// keep this contract usable from tpfalinearizer.hh without introducing its
// dependency cycle through the Flow executable TypeTags.
struct FlowProblem;
struct FlowGasWaterEnergyProblem;
struct FlowGasWaterEnergyProblemGPU;

/*!
 * \brief The one GPU TypeTag family shared by property evaluation and TPFA.
 *
 * The Storage parameter is intentionally part of the type identity.  Device
 * kernels use the GpuView instantiation; retaining the parameter makes the
 * mapping through to_gpu_type explicit and prevents accidental reuse of a
 * similarly laid-out, but differently bound, device IQ type.
 */
template <template <class> class Storage>
struct FlowGasWaterEnergyDeviceTypeTag {
    using InheritsFrom = std::tuple<FlowGasWaterEnergyProblem>;
};

/*!
 * \brief CPU-side source TypeTag for constructing device IQ prototypes.
 *
 * Its feature set is identical to FlowGasWaterEnergyDeviceTypeTag so the
 * BlackOilIntensiveQuantities converting constructor remains well defined.
 */
template <class ParentTypeTag>
struct FlowGasWaterEnergyCpuDeviceContract {
    using InheritsFrom = std::tuple<FlowGasWaterEnergyProblem>;
};

// Compatibility names for the pre-contract GPU IQ test and external users.
using FlowGasWaterEnergyKernelBaseGPU = FlowGasWaterEnergyDeviceTypeTag<gpuistl::GpuView>;
using FlowGasWaterEnergyDummyProblemGPU = FlowGasWaterEnergyDeviceTypeTag<gpuistl::GpuView>;
using FlowGasWaterEnergyCpuKernelBase = FlowGasWaterEnergyCpuDeviceContract<FlowGasWaterEnergyProblem>;

template <template <class> class Storage>
struct to_gpu_type<FlowGasWaterEnergyProblemGPU, Storage> {
    using type = FlowGasWaterEnergyDeviceTypeTag<Storage>;
};

template <template <class> class Storage>
struct to_gpu_type<FlowGasWaterEnergyProblem, Storage> {
    using type = FlowGasWaterEnergyDeviceTypeTag<Storage>;
};

} // namespace TTag

// Unsupported BlackOil modules must not change the device IQ layout.
template <class TypeTag, template <class> class Storage>
struct EnableDiffusion<TypeTag, TTag::FlowGasWaterEnergyDeviceTypeTag<Storage>> {
    static constexpr bool value = false;
};

template <class TypeTag, template <class> class Storage>
struct EnableDispersion<TypeTag, TTag::FlowGasWaterEnergyDeviceTypeTag<Storage>> {
    static constexpr bool value = false;
};

template <class TypeTag, template <class> class Storage>
struct EnableSolvent<TypeTag, TTag::FlowGasWaterEnergyDeviceTypeTag<Storage>> {
    static constexpr bool value = false;
};

template <class TypeTag, template <class> class Storage>
struct EnableExtbo<TypeTag, TTag::FlowGasWaterEnergyDeviceTypeTag<Storage>> {
    static constexpr bool value = false;
};

template <class TypeTag, template <class> class Storage>
struct EnablePolymer<TypeTag, TTag::FlowGasWaterEnergyDeviceTypeTag<Storage>> {
    static constexpr bool value = false;
};

template <class TypeTag, template <class> class Storage>
struct EnablePolymerMW<TypeTag, TTag::FlowGasWaterEnergyDeviceTypeTag<Storage>> {
    static constexpr bool value = false;
};

template <class TypeTag, template <class> class Storage>
struct EnableFoam<TypeTag, TTag::FlowGasWaterEnergyDeviceTypeTag<Storage>> {
    static constexpr bool value = false;
};

template <class TypeTag, template <class> class Storage>
struct EnableBrine<TypeTag, TTag::FlowGasWaterEnergyDeviceTypeTag<Storage>> {
    static constexpr bool value = false;
};

template <class TypeTag, template <class> class Storage>
struct EnableSaltPrecipitation<TypeTag, TTag::FlowGasWaterEnergyDeviceTypeTag<Storage>> {
    static constexpr bool value = false;
};

template <class TypeTag, template <class> class Storage>
struct EnableBioeffects<TypeTag, TTag::FlowGasWaterEnergyDeviceTypeTag<Storage>> {
    static constexpr bool value = false;
};

template <class TypeTag, template <class> class Storage>
struct EnableMech<TypeTag, TTag::FlowGasWaterEnergyDeviceTypeTag<Storage>> {
    static constexpr bool value = false;
};

template <class TypeTag, template <class> class Storage>
struct EnableGeochemistry<TypeTag, TTag::FlowGasWaterEnergyDeviceTypeTag<Storage>> {
    static constexpr bool value = false;
};

template <class TypeTag, template <class> class Storage>
struct EnableConvectiveMixing<TypeTag, TTag::FlowGasWaterEnergyDeviceTypeTag<Storage>> {
    static constexpr bool value = false;
};

template <class TypeTag, template <class> class Storage>
struct EnableEnergy<TypeTag, TTag::FlowGasWaterEnergyDeviceTypeTag<Storage>> {
    static constexpr bool value = true;
};

// Match the device contract when making an IQ prototype on the host.
template <class TypeTag, class ParentTypeTag>
struct EnableDiffusion<TypeTag, TTag::FlowGasWaterEnergyCpuDeviceContract<ParentTypeTag>> {
    static constexpr bool value = false;
};

template <class TypeTag, class ParentTypeTag>
struct EnableDispersion<TypeTag, TTag::FlowGasWaterEnergyCpuDeviceContract<ParentTypeTag>> {
    static constexpr bool value = false;
};

// Material law: GPU-friendly EclTwoPhaseMaterial backed by device views.
template <class TypeTag, template <class> class Storage>
struct MaterialLaw<TypeTag, TTag::FlowGasWaterEnergyDeviceTypeTag<Storage>>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Traits = ThreePhaseMaterialTraits<Scalar,
                                             FluidSystem::waterPhaseIdx,
                                             FluidSystem::oilPhaseIdx,
                                             FluidSystem::gasPhaseIdx,
                                             getPropValue<TypeTag, Properties::EnableHysteresis>(),
                                             getPropValue<TypeTag, Properties::EnableEndpointScaling>()>;
    using TwoPhaseTraits = TwoPhaseMaterialTraits<Scalar,
                                                   Traits::wettingPhaseIdx,
                                                   Traits::nonWettingPhaseIdx>;
    using TwoPhaseParams =
        ::Opm::PiecewiseLinearTwoPhaseMaterialParams<TwoPhaseTraits, ::Opm::gpuistl::GpuView<const Scalar>>;
    using TwoPhaseLaw = ::Opm::PiecewiseLinearTwoPhaseMaterial<TwoPhaseTraits, TwoPhaseParams>;
    using GpuMaterialLawParams =
        ::Opm::EclTwoPhaseMaterialParams<Traits,
                                         TwoPhaseParams,
                                         TwoPhaseParams,
                                         TwoPhaseParams,
                                         ::Opm::gpuistl::ValueAsPointer>;
    using GpuMaterialLaw =
        ::Opm::EclTwoPhaseMaterial<Traits, TwoPhaseLaw, TwoPhaseLaw, TwoPhaseLaw, GpuMaterialLawParams>;

public:
    using EclMaterialLawManager =
        ::Opm::EclMaterialLaw::GpuManager<Traits,
                                           TwoPhaseLaw,
                                           TwoPhaseLaw,
                                           ::Opm::VectorWithDefaultAllocator,
                                           GpuMaterialLaw>;
    using type = typename EclMaterialLawManager::MaterialLaw;
};

template <class TypeTag, template <class> class Storage>
struct PrimaryVariables<TypeTag, TTag::FlowGasWaterEnergyDeviceTypeTag<Storage>> {
    using type = ::Opm::BlackOilPrimaryVariables<TypeTag, ::Opm::gpuistl::MiniVector>;
};

template <class TypeTag, template <class> class Storage>
struct FluidSystem<TypeTag, TTag::FlowGasWaterEnergyDeviceTypeTag<Storage>> {
    using type = ::Opm::BlackOilFluidSystemNonStatic<
        GetPropType<TypeTag, Properties::Scalar>,
        ::Opm::BlackOilDefaultFluidSystemIndices,
        Storage>;
};

template <class TypeTag, template <class> class Storage>
struct SolidEnergyLaw<TypeTag, TTag::FlowGasWaterEnergyDeviceTypeTag<Storage>>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

public:
    using EclThermalLawManager =
        ::Opm::EclThermalLaw::GpuManager<Scalar, FluidSystem, ::Opm::gpuistl::GpuView, ::Opm::gpuistl::GpuView>;
    using type = ::Opm::EclSpecrockLaw<Scalar, ::Opm::EclSpecrockLawParams<Scalar, ::Opm::gpuistl::GpuView>>;
};

template <class TypeTag, template <class> class Storage>
struct ThermalConductionLaw<TypeTag, TTag::FlowGasWaterEnergyDeviceTypeTag<Storage>>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

public:
    using EclThermalLawManager =
        ::Opm::EclThermalLaw::GpuManager<Scalar, FluidSystem, ::Opm::gpuistl::GpuView, ::Opm::gpuistl::GpuView>;
    using type = ::Opm::EclThconrLaw<Scalar, FluidSystem>;
};

template <class TypeTag, template <class> class Storage>
struct Problem<TypeTag, TTag::FlowGasWaterEnergyDeviceTypeTag<Storage>>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using MaterialLawManager = typename GetProp<TypeTag, Properties::MaterialLaw>::EclMaterialLawManager;
    using GpuViewMaterialLawManager =
        ::Opm::EclMaterialLaw::GpuManager<typename MaterialLawManager::Traits,
                                           typename MaterialLawManager::GasOilLaw,
                                           typename MaterialLawManager::OilWaterLaw,
                                           ::Opm::gpuistl::GpuView,
                                           typename MaterialLawManager::MaterialLaw>;
    using GpuViewThermalLawManager =
        ::Opm::EclThermalLaw::GpuManager<Scalar, FluidSystem, ::Opm::gpuistl::GpuView, ::Opm::gpuistl::GpuView>;

public:
    using type = ::Opm::GpuFlowProblem<
        Scalar, GpuViewMaterialLawManager, ::Opm::gpuistl::GpuView, GpuViewThermalLawManager>;
};

template <class TypeTag, template <class> class Storage>
struct ElementContext<TypeTag, TTag::FlowGasWaterEnergyDeviceTypeTag<Storage>> {
    using type = std::nullptr_t;
};

} // namespace Opm::Properties

#endif // HAVE_CUDA

#endif // OPM_GPU_FLOW_GASWATER_ENERGY_CONTRACT_HPP
