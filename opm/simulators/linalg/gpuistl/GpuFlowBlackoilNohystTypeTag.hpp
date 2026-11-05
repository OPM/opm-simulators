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
 * GPU TypeTag for the three-phase black-oil configuration with live oil and
 * dry gas (e.g. SPE1), derived from the \c FlowProblemNohyst TypeTag that
 * \c flow uses for such decks.
 *
 * Only the properties that must differ on the GPU are overridden: the storage
 * of the primary variables, the fluid system, the material law, the problem
 * and the element context. Everything else, in particular the module switches
 * that determine the layout of \c BlackOilIntensiveQuantities, is inherited
 * from \c FlowProblemNohyst. The CPU TypeTag can therefore be the source of
 * the cross-TypeTag converting constructors directly.
 */
#ifndef OPM_GPU_FLOW_BLACKOIL_NOHYST_TYPETAG_HPP
#define OPM_GPU_FLOW_BLACKOIL_NOHYST_TYPETAG_HPP

#if HAVE_CUDA

#include <opm/common/utility/VectorWithDefaultAllocator.hpp>

#include <opm/material/fluidmatrixinteractions/EclMaterialLawTwoPhaseTypes.hpp>
#include <opm/material/fluidmatrixinteractions/GpuEclMaterialLawManager.hpp>
#include <opm/material/fluidmatrixinteractions/PiecewiseLinearTwoPhaseMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/PiecewiseLinearTwoPhaseMaterialParams.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystemNonStatic.hpp>
#include <opm/material/fluidsystems/blackoilpvt/BlackOilPvtApproaches.hpp>

#include <opm/models/blackoil/blackoilintensivequantities.hh>
#include <opm/models/blackoil/blackoilprimaryvariables.hh>

#include <opm/simulators/flow/FlowProblemNohystTypeTag.hpp>
#include <opm/simulators/flow/GpuFlowProblem.hpp>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>
#include <opm/simulators/linalg/gpuistl/MiniVector.hpp>

#include <cstddef>
#include <tuple>
#include <type_traits>

namespace Opm::Properties {

namespace TTag {

/// The TypeTag of the GPU intensive quantities: \c FlowProblemNohyst with
/// GPU-compatible primary variables, fluid system, material law and problem.
struct FlowBlackoilNohystGPU {
    using InheritsFrom = std::tuple<FlowProblemNohyst>;
};

} // namespace TTag

// -----------------------------------------------------------------------
// MaterialLaw: EclDefaultMaterial with GpuView-backed piecewise-linear
// two-phase laws. The traits are those of the CPU law (no hysteresis,
// endpoint scaling enabled), so the GpuManager accepts the CPU manager.
// The GPU manager does not apply endpoint scaling; decks that need it
// must be rejected before reaching the GPU (see OPM/opm-simulators#7644).
// -----------------------------------------------------------------------
template <class TypeTag>
struct MaterialLaw<TypeTag, TTag::FlowBlackoilNohystGPU>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Traits = typename GetPropType<TTag::FlowProblemNohyst, Properties::MaterialLaw>::Traits;
    using TwoPhaseTypes = ::Opm::EclMaterialLaw::TwoPhaseTypes<Traits>;
    using GasOilTraits = typename TwoPhaseTypes::GasOilTraits;
    using OilWaterTraits = typename TwoPhaseTypes::OilWaterTraits;
    using GasOilParams = ::Opm::PiecewiseLinearTwoPhaseMaterialParams<
        GasOilTraits, ::Opm::gpuistl::GpuView<const Scalar>>;
    using OilWaterParams = ::Opm::PiecewiseLinearTwoPhaseMaterialParams<
        OilWaterTraits, ::Opm::gpuistl::GpuView<const Scalar>>;
    using GasOilLaw = ::Opm::PiecewiseLinearTwoPhaseMaterial<GasOilTraits, GasOilParams>;
    using OilWaterLaw = ::Opm::PiecewiseLinearTwoPhaseMaterial<OilWaterTraits, OilWaterParams>;

public:
    using EclMaterialLawManager = ::Opm::EclMaterialLaw::GpuManager<
        Traits, GasOilLaw, OilWaterLaw, ::Opm::VectorWithDefaultAllocator>;
    using type = typename EclMaterialLawManager::MaterialLaw;
};

template <class TypeTag>
struct PrimaryVariables<TypeTag, TTag::FlowBlackoilNohystGPU>
{
    using type = ::Opm::BlackOilPrimaryVariables<TypeTag, ::Opm::gpuistl::MiniVector>;
};

template <class TypeTag>
struct FluidSystem<TypeTag, TTag::FlowBlackoilNohystGPU>
{
    using type = ::Opm::BlackOilFluidSystemNonStatic<
        GetPropType<TypeTag, Properties::Scalar>,
        ::Opm::BlackOilDefaultFluidSystemIndices,
        ::Opm::gpuistl::GpuView,
        ::Opm::LiveOilDryGasPvtApproaches>;
};

// -----------------------------------------------------------------------
// Problem: GpuFlowProblem (GpuView-backed) without thermal laws.
// -----------------------------------------------------------------------
template <class TypeTag>
struct Problem<TypeTag, TTag::FlowBlackoilNohystGPU>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using CpuStorageMaterialLawManager =
        typename GetProp<TypeTag, Properties::MaterialLaw>::EclMaterialLawManager;
    using GpuViewMaterialLawManager =
        ::Opm::EclMaterialLaw::GpuManager<typename CpuStorageMaterialLawManager::Traits,
                                          typename CpuStorageMaterialLawManager::GasOilLaw,
                                          typename CpuStorageMaterialLawManager::OilWaterLaw,
                                          ::Opm::gpuistl::GpuView,
                                          typename CpuStorageMaterialLawManager::MaterialLaw>;

public:
    using type = ::Opm::GpuFlowProblem<Scalar, GpuViewMaterialLawManager, ::Opm::gpuistl::GpuView>;
};

template <class TypeTag>
struct ElementContext<TypeTag, TTag::FlowBlackoilNohystGPU>
{
    using type = std::nullptr_t;
};

} // namespace Opm::Properties

namespace Opm::detail {

// Consistency checks between the CPU TypeTag and the GPU TypeTag.
struct CheckFlowBlackoilNohystGpuTypeTags
{
    using CpuTag = Properties::TTag::FlowProblemNohyst;
    using GpuTag = Properties::TTag::FlowBlackoilNohystGPU;

    using GpuMaterialLaw = GetPropType<GpuTag, Properties::MaterialLaw>;
    using CpuMaterialLaw = GetPropType<CpuTag, Properties::MaterialLaw>;

    // The intensive quantities use the MaterialLaw property for the capillary
    // pressures; GpuFlowProblem uses the law of its manager for the relative
    // permeabilities and for the parameters returned by materialLawParams().
    static_assert(std::is_same_v<GpuMaterialLaw,
                                 typename GetPropType<GpuTag, Properties::Problem>::
                                     EclMaterialLawManager::MaterialLaw>);
    static_assert(std::is_same_v<typename GpuMaterialLaw::Traits,
                                 typename CpuMaterialLaw::Traits>);

    static_assert(std::is_same_v<GetPropType<GpuTag, Properties::Indices>,
                                 GetPropType<CpuTag, Properties::Indices>>);
    static_assert(std::is_same_v<GetPropType<GpuTag, Properties::Evaluation>,
                                 GetPropType<CpuTag, Properties::Evaluation>>);
    static_assert(getPropValue<GpuTag, Properties::NumPhases>()
                  == getPropValue<CpuTag, Properties::NumPhases>());
    static_assert(getPropValue<GpuTag, Properties::EnergyModuleType>()
                  == getPropValue<CpuTag, Properties::EnergyModuleType>());

    static_assert(std::is_trivially_copyable_v<
                  typename BlackOilIntensiveQuantities<GpuTag>::FluidState>,
                  "the fluid state is copied byte-wise to and from the device");
};

} // namespace Opm::detail

#endif // HAVE_CUDA

#endif // OPM_GPU_FLOW_BLACKOIL_NOHYST_TYPETAG_HPP
