// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \copydoc Opm::FlowBaseProblemComp
 */
#ifndef OPM_FLOW_PROBLEM_COMP_PROPERTIES_HPP
#define OPM_FLOW_PROBLEM_COMP_PROPERTIES_HPP


#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>
#include <opm/material/thermal/EclThermalLawManager.hpp>

#include <opm/models/discretization/ecfv/ecfvdiscretization.hh>
#include <opm/models/utils/propertysystem.hh>

#include <opm/simulators/flow/BaseAquiferModel.hpp>
#include <opm/simulators/flow/CpGridVanguard.hpp>
#include <opm/simulators/flow/DummyGradientCalculator.hpp>
#include <opm/simulators/flow/EclWriter.hpp>
#include <opm/simulators/flow/FlowProblemParameters.hpp>
#include <opm/simulators/flow/FIBlackoilModel.hpp>
#include <opm/simulators/flow/NewTranFluxModule.hpp>
#include <opm/simulators/flow/OutputBlackoilModule.hpp>
#include <opm/simulators/flow/VtkTracerModule.hpp>
#include <opm/simulators/flow/TracerModel.hpp> // we have an empty TracerModel later
#if HAVE_DAMARIS
#include <opm/simulators/flow/DamarisWriter.hpp>
#endif

#include <tuple>

namespace Opm {
template <class TypeTag>
class FlowProblemComp;
}

namespace Opm::Properties {

namespace TTag {

struct FlowBaseProblemComp {
    using InheritsFrom = std::tuple<CpGridVanguard>;
};

}

// The class which deals with ECL aquifers
template<class TypeTag, class MyTypeTag>
struct AquiferModel { using type = UndefinedProperty; };

// Specify whether API tracking should be enabled (replaces PVT regions).
// TODO: This is not yet implemented
template<class TypeTag, class MyTypeTag>
struct EnableApiTracking { using type = UndefinedProperty; };

// Enable the additional checks even if compiled in debug mode (i.e., with the NDEBUG
// macro undefined). Next to a slightly better performance, this also eliminates some
// print statements in debug mode.
template<class TypeTag, class MyTypeTag>
struct EnableDebuggingChecks { using type = Properties::UndefinedProperty; };
// if thermal flux boundaries are enabled an effort is made to preserve the initial
// thermal gradient specified via the TEMPVD keyword
template<class TypeTag, class MyTypeTag>
struct EnableThermalFluxBoundaries { using type = UndefinedProperty; };

// The class which deals with wells
template<class TypeTag, class MyTypeTag>
struct WellModel { using type = UndefinedProperty; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::FlowBaseProblemComp>
{ using type = FlowProblemComp<TypeTag>; };

/* template<class TypeTag>
struct Model<TypeTag, TTag::FlowBaseProblem>
{ using type = FIBlackOilModel<TypeTag>; }; */

// TODO: for the purpose of the compositiona model
template<class TypeTag, class MyTypeTag>
struct TracerModelDef {
    using type = UndefinedProperty;
};

template<class TypeTag>
struct TracerModelDef<TypeTag, TTag::FlowBaseProblemComp> {
    using type = ::Opm::TracerModel<TypeTag>;
};

// Select the element centered finite volume method as spatial discretization
template<class TypeTag>
struct SpatialDiscretizationSplice<TypeTag, TTag::FlowBaseProblemComp>
{ using type = TTag::EcfvDiscretization; };

// use automatic differentiation to linearize the system of PDEs
template<class TypeTag>
struct LocalLinearizerSplice<TypeTag, TTag::FlowBaseProblemComp>
{ using type = TTag::AutoDiffLocalLinearizer; };

template<class TypeTag>
struct BaseDiscretizationType<TypeTag, TTag::FlowBaseProblemComp>
{ using type = FvBaseDiscretizationNoAdapt<TypeTag>; };

template<class TypeTag>
struct DiscreteFunction<TypeTag, TTag::FlowBaseProblemComp>
{
    using BaseDiscretization = FvBaseDiscretization<TypeTag>;
    using type = typename BaseDiscretization::BlockVectorWrapper;
};

template<class TypeTag>
struct GridView<TypeTag, TTag::FlowBaseProblemComp>
{ using type = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView; };

// Set the material law for fluid fluxes
template<class TypeTag>
struct MaterialLaw<TypeTag, TTag::FlowBaseProblemComp>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    // using Traits = ThreePhaseMaterialTraits<Scalar,
    //                                         /*wettingPhaseIdx=*/FluidSystem::waterPhaseIdx,
    //                                         /*nonWettingPhaseIdx=*/FluidSystem::oilPhaseIdx,
    //                                         /*gasPhaseIdx=*/FluidSystem::gasPhaseIdx>;

    // TODO: We should be able to use FluidSystem here and using Indices to handle the active phases    
        using Traits = ThreePhaseMaterialTraits<Scalar,
                                            /*wettingPhaseIdx=*/ 0,
                                            /*nonWettingPhaseIdx=*/ 1,
                                            /*gasPhaseIdx=*/ 2>;

public:
    using EclMaterialLawManager = ::Opm::EclMaterialLawManager<Traits>;

    using type = typename EclMaterialLawManager::MaterialLaw;
};

// Set the material law for energy storage in rock
template<class TypeTag>
struct SolidEnergyLaw<TypeTag, TTag::FlowBaseProblemComp>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

public:
    using EclThermalLawManager = ::Opm::EclThermalLawManager<Scalar, FluidSystem>;

    using type = typename EclThermalLawManager::SolidEnergyLaw;
};

// Set the material law for thermal conduction
template<class TypeTag>
struct ThermalConductionLaw<TypeTag, TTag::FlowBaseProblemComp>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

public:
    using EclThermalLawManager = ::Opm::EclThermalLawManager<Scalar, FluidSystem>;

    using type = typename EclThermalLawManager::ThermalConductionLaw;
};

// use a slightly faster stencil class because it does not need the normals and
// the integration points of intersections
template<class TypeTag>
struct Stencil<TypeTag, TTag::FlowBaseProblemComp>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;

public:
    using type = EcfvStencil<Scalar,
                             GridView,
                             /*needIntegrationPos=*/false,
                             /*needNormal=*/false>;
};

// by default use the dummy aquifer "model"
template<class TypeTag>
struct AquiferModel<TypeTag, TTag::FlowBaseProblemComp> {
    using type = BaseAquiferModel<TypeTag>;
};


// Enable diffusion // TODO:
template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::FlowBaseProblemComp>
{ static constexpr bool value = false; };

// Disable dispersion
template<class TypeTag>
struct EnableDispersion<TypeTag, TTag::FlowBaseProblemComp>
{ static constexpr bool value = false; };

// Enable Convective Mixing
template<class TypeTag>
struct EnableConvectiveMixing<TypeTag, TTag::FlowBaseProblemComp>
{ static constexpr bool value = true; };

// disable API tracking
template<class TypeTag>
struct EnableApiTracking<TypeTag, TTag::FlowBaseProblemComp>
{ static constexpr bool value = false; };

// Use the "velocity module" which uses the Eclipse "NEWTRAN" transmissibilities
/* template<class TypeTag>
struct FluxModule<TypeTag, TTag::FlowBaseProblemComp>
{ using type = NewTranFluxModule<TypeTag>; }; */ // TODO: check whether it works

// Use the dummy gradient calculator in order not to do unnecessary work.
// TODO: remove this for the comp modeling
// template<class TypeTag>
// struct GradientCalculator<TypeTag, TTag::FlowBaseProblemComp>
// { using type = DummyGradientCalculator<TypeTag>; };

// store temperature (but do not conserve energy, as long as EnableEnergy is false)
template<class TypeTag>
struct EnableTemperature<TypeTag, TTag::FlowBaseProblemComp>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableMech<TypeTag, TTag::FlowBaseProblemComp>
{ static constexpr bool value = false; };

// disable all extensions supported by black oil model. this should not really be
// necessary but it makes things a bit more explicit
template<class TypeTag>
struct EnablePolymer<TypeTag, TTag::FlowBaseProblemComp>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableSolvent<TypeTag, TTag::FlowBaseProblemComp>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableEnergy<TypeTag, TTag::FlowBaseProblemComp>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableFoam<TypeTag, TTag::FlowBaseProblemComp>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableExtbo<TypeTag, TTag::FlowBaseProblemComp>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableMICP<TypeTag, TTag::FlowBaseProblemComp>
{ static constexpr bool value = false; };

// disable thermal flux boundaries by default
template<class TypeTag>
struct EnableThermalFluxBoundaries<TypeTag, TTag::FlowBaseProblemComp>
{ static constexpr bool value = false; };

// By default, simulators derived from the FlowBaseProblemComp are production simulators,
// i.e., experimental features must be explicitly enabled at compile time
template<class TypeTag>
struct EnableExperiments<TypeTag, TTag::FlowBaseProblemComp>
{ static constexpr bool value = false; };
// By default, we enable the debugging checks if we're compiled in debug mode
template<class TypeTag>
struct EnableDebuggingChecks<TypeTag, TTag::FlowBaseProblemComp>
{ static constexpr bool value = true; };

} // namespace Opm::Properties


#endif // OPM_FLOW_PROBLEM_COMP_PROPERTIES_HPP
