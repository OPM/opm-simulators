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
 * Shared TypeTag declaration for the gas+water+energy (CO2STORE-style)
 * blackoil flow problem. Originally inlined into
 * \c flow/flow_gaswater_energy.cpp; extracted into a header so that other
 * translation units (e.g. the GPU intensive-quantities dispatcher) can
 * obtain matching property bindings via \c GetPropType.
 */
#ifndef OPM_FLOW_GASWATER_ENERGY_TYPETAG_HPP
#define OPM_FLOW_GASWATER_ENERGY_TYPETAG_HPP

#include <opm/material/thermal/EnergyModuleType.hpp>

#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>

#include <opm/simulators/flow/FlowProblemBlackoilProperties.hpp>

namespace Opm::Properties {
namespace TTag {
struct FlowGasWaterEnergyProblem {
    using InheritsFrom = std::tuple<FlowProblem>;
};

/// GPU-simulation variant: all physics properties are inherited from
/// \c FlowGasWaterEnergyProblem; the GPU-specific property overrides
/// (RunAssemblyOnGpu, GpuFIBlackOilModel, …) are added in \c flow_gpu.cu.
struct FlowGasWaterEnergyProblemGPU {
    using InheritsFrom = std::tuple<FlowGasWaterEnergyProblem>;
};

} // namespace TTag

// -----------------------------------------------------------------------
// FlowGasWaterEnergyProblemGPU property overrides that affect the layout
// of \c BlackOilIntensiveQuantities (and therefore must be visible to
// every translation unit that instantiates that class with this TypeTag,
// otherwise we get an ODR violation: e.g. the dispatcher TU and
// \c flow_gpu.cu would compile differently sized / ordered IQ structs).
//
// Dispersion is not yet supported on the GPU assembly path; disable it
// even though the CPU parent TypeTag enables it.
// -----------------------------------------------------------------------
template <class TypeTag>
struct EnableDispersion<TypeTag, TTag::FlowGasWaterEnergyProblemGPU> {
    static constexpr bool value = false;
};

template <class TypeTag>
struct EnableDiffusion<TypeTag, TTag::FlowGasWaterEnergyProblemGPU> {
    static constexpr bool value = false;
};

template <class TypeTag>
struct EnableEnergy<TypeTag, TTag::FlowGasWaterEnergyProblemGPU> {
    static constexpr bool value = true;
};

template<class TypeTag>
struct Linearizer<TypeTag, TTag::FlowGasWaterEnergyProblem> { using type = TpfaLinearizer<TypeTag>; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::FlowGasWaterEnergyProblem> { using type = BlackOilLocalResidualTPFA<TypeTag>; };

template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::FlowGasWaterEnergyProblem> { static constexpr bool value = true; };

template<class TypeTag>
struct EnableDispersion<TypeTag, TTag::FlowGasWaterEnergyProblem> { static constexpr bool value = true; };

template<class TypeTag>
struct EnergyModuleType<TypeTag, TTag::FlowGasWaterEnergyProblem>
{ static constexpr EnergyModules value = EnergyModules::FullyImplicitThermal; };

template<class TypeTag>
struct EnableDisgasInWater<TypeTag, TTag::FlowGasWaterEnergyProblem> {
    static constexpr bool value = true;
};

template<class TypeTag>
struct EnableVapwat<TypeTag, TTag::FlowGasWaterEnergyProblem> {
    static constexpr bool value = true;
};

//! The indices required by the model
template<class TypeTag>
struct Indices<TypeTag, TTag::FlowGasWaterEnergyProblem>
{
private:
    // it is unfortunately not possible to simply use 'TypeTag' here because
    // this leads to cyclic definitions of some properties.
    using BaseTypeTag = TTag::FlowProblem;
    using FluidSystem = GetPropType<BaseTypeTag, Properties::FluidSystem>;
    static constexpr EnergyModules energyModuleType = getPropValue<TypeTag, Properties::EnergyModuleType>();
    static constexpr int numEnergyVars = energyModuleType == EnergyModules::FullyImplicitThermal;
public:
    using type = BlackOilTwoPhaseIndices<getPropValue<TypeTag, Properties::EnableSolvent>(),
                                         getPropValue<TypeTag, Properties::EnableExtbo>(),
                                         getPropValue<TypeTag, Properties::EnablePolymer>(),
                                         numEnergyVars,
                                         getPropValue<TypeTag, Properties::EnableFoam>(),
                                         getPropValue<TypeTag, Properties::EnableBrine>(),
                                         /*PVOffset=*/0,
                                         /*disabledCompIdx=*/FluidSystem::oilCompIdx,
                                         getPropValue<TypeTag, Properties::EnableBioeffects>()>;
};

} // namespace Opm::Properties

#endif // OPM_FLOW_GASWATER_ENERGY_TYPETAG_HPP
