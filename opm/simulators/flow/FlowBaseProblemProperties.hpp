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
 * \copydoc Opm::FlowBaseProblem
 */
#ifndef OPM_FLOW_BASE_PROBLEM_PROPERTIES_HPP
#define OPM_FLOW_BASE_PROBLEM_PROPERTIES_HPP


#include <opm/material/thermal/EclThermalLawManager.hpp>

#include <opm/models/discretization/ecfv/ecfvdiscretization.hh>
#include <opm/models/utils/propertysystem.hh>

#include <opm/simulators/flow/BaseAquiferModel.hpp>
#include <opm/simulators/flow/CpGridVanguard.hpp>
#include <opm/simulators/flow/DummyGradientCalculator.hpp>
#include <opm/simulators/flow/EclWriter.hpp>
#include <opm/simulators/flow/FlowProblemParameters.hpp>
#include <opm/simulators/flow/TracerModel.hpp>

#if HAVE_DAMARIS
#include <opm/simulators/flow/DamarisWriter.hpp>
#endif

#include <tuple>

namespace Opm::Properties {

namespace TTag {

struct FlowBaseProblem {
    using InheritsFrom = std::tuple<CpGridVanguard>;
};

}

template<class TypeTag, class MyTypeTag>
struct NonlinearSystem { using type = UndefinedProperty; };



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

// Avoid using ElementContext-based code if possible.
template<class TypeTag, class MyTypeTag>
struct AvoidElementContext { using type = Properties::UndefinedProperty; };

// if thermal flux boundaries are enabled an effort is made to preserve the initial
// thermal gradient specified via the TEMPVD keyword
template<class TypeTag, class MyTypeTag>
struct EnableThermalFluxBoundaries { using type = UndefinedProperty; };

// The class which deals with wells
template<class TypeTag, class MyTypeTag>
struct WellModel { using type = UndefinedProperty; };

// Tracer might be moved to the blackoil side
// The class that deals with the tracer
template<class TypeTag, class MyTypeTag>
struct TracerModel {  using type = UndefinedProperty; };

template <class TypeTag>
struct TracerModel<TypeTag, TTag::FlowBaseProblem>
{ using type =  ::Opm::TracerModel<TypeTag>; };

// Select the element centered finite volume method as spatial discretization
template<class TypeTag>
struct SpatialDiscretizationSplice<TypeTag, TTag::FlowBaseProblem>
{ using type = TTag::EcfvDiscretization; };

// use automatic differentiation to linearize the system of PDEs
template<class TypeTag>
struct LocalLinearizerSplice<TypeTag, TTag::FlowBaseProblem>
{ using type = TTag::AutoDiffLocalLinearizer; };

template<class TypeTag>
struct BaseDiscretizationType<TypeTag, TTag::FlowBaseProblem>
{ using type = FvBaseDiscretizationNoAdapt<TypeTag>; };

template<class TypeTag>
struct DiscreteFunction<TypeTag, TTag::FlowBaseProblem>
{
    using BaseDiscretization = FvBaseDiscretization<TypeTag>;
    using type = typename BaseDiscretization::BlockVectorWrapper;
};

template<class TypeTag>
struct GridView<TypeTag, TTag::FlowBaseProblem>
{ using type = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView; };

// Set the material law for energy storage in rock
template<class TypeTag>
struct SolidEnergyLaw<TypeTag, TTag::FlowBaseProblem>
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
struct ThermalConductionLaw<TypeTag, TTag::FlowBaseProblem>
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
struct Stencil<TypeTag, TTag::FlowBaseProblem>
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
struct AquiferModel<TypeTag, TTag::FlowBaseProblem> {
    using type = BaseAquiferModel<TypeTag>;
};

// Enable diffusion
template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::FlowBaseProblem>
{ static constexpr bool value = true; };

// Disable dispersion
template<class TypeTag>
struct EnableDispersion<TypeTag, TTag::FlowBaseProblem>
{ static constexpr bool value = false; };

// Enable Convective Mixing
template<class TypeTag>
struct EnableConvectiveMixing<TypeTag, TTag::FlowBaseProblem>
{ static constexpr bool value = true; };

// disable API tracking
template<class TypeTag>
struct EnableApiTracking<TypeTag, TTag::FlowBaseProblem>
{ static constexpr bool value = false; };

// store temperature (but do not conserve energy, as long as EnableEnergy is false)
template<class TypeTag>
struct EnableTemperature<TypeTag, TTag::FlowBaseProblem>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableMech<TypeTag, TTag::FlowBaseProblem>
{ static constexpr bool value = false; };

// disable all extensions supported by black oil model. this should not really be
// necessary but it makes things a bit more explicit
template<class TypeTag>
struct EnablePolymer<TypeTag, TTag::FlowBaseProblem>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableSolvent<TypeTag, TTag::FlowBaseProblem>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableEnergy<TypeTag, TTag::FlowBaseProblem>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableFoam<TypeTag, TTag::FlowBaseProblem>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableExtbo<TypeTag, TTag::FlowBaseProblem>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableMICP<TypeTag, TTag::FlowBaseProblem>
{ static constexpr bool value = false; };

// disable thermal flux boundaries by default
template<class TypeTag>
struct EnableThermalFluxBoundaries<TypeTag, TTag::FlowBaseProblem>
{ static constexpr bool value = false; };

// By default, simulators derived from the FlowBaseProblem are production simulators,
// i.e., experimental features must be explicitly enabled at compile time
template<class TypeTag>
struct EnableExperiments<TypeTag, TTag::FlowBaseProblem>
{ static constexpr bool value = false; };

// By default, we enable the debugging checks if we're compiled in debug mode
template<class TypeTag>
struct EnableDebuggingChecks<TypeTag, TTag::FlowBaseProblem>
{ static constexpr bool value = true; };

// Most modules are implemented only in terms of element contexts,
// so this must default to false.
template<class TypeTag>
struct AvoidElementContext<TypeTag, TTag::FlowBaseProblem>
{ static constexpr bool value = false; };

} // namespace Opm::Properties

#endif // OPM_BASE_FLOW_PROBLEM_PROPERTIES_HPP
