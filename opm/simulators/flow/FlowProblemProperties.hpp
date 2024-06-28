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
 * \copydoc Opm::FlowProblem
 */
#ifndef OPM_FLOW_PROBLEM_PROPERTIES_HPP
#define OPM_FLOW_PROBLEM_PROPERTIES_HPP

#include <opm/input/eclipse/Parser/ParserKeywords/E.hpp>

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

#if HAVE_DAMARIS
#include <opm/simulators/flow/DamarisWriter.hpp>
#endif

#include <tuple>

namespace Opm {
template <class TypeTag>
class FlowProblem;
}

namespace Opm::Properties {

namespace TTag {

struct FlowBaseProblem {
    using InheritsFrom = std::tuple<VtkTracer, OutputBlackOil, CpGridVanguard>;
};

}

// The class which deals with ECL aquifers
template<class TypeTag, class MyTypeTag>
struct AquiferModel { using type = UndefinedProperty; };

// Specify whether API tracking should be enabled (replaces PVT regions).
// TODO: This is not yet implemented
template<class TypeTag, class MyTypeTag>
struct EnableApiTracking { using type = UndefinedProperty; };

// if thermal flux boundaries are enabled an effort is made to preserve the initial
// thermal gradient specified via the TEMPVD keyword
template<class TypeTag, class MyTypeTag>
struct EnableThermalFluxBoundaries { using type = UndefinedProperty; };

// The class which deals with wells
template<class TypeTag, class MyTypeTag>
struct WellModel { using type = UndefinedProperty; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::FlowBaseProblem>
{ using type = FlowProblem<TypeTag>; };

template<class TypeTag>
struct Model<TypeTag, TTag::FlowBaseProblem>
{ using type = FIBlackOilModel<TypeTag>; };

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

// Set the material law for fluid fluxes
template<class TypeTag>
struct MaterialLaw<TypeTag, TTag::FlowBaseProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    using Traits = ThreePhaseMaterialTraits<Scalar,
                                            /*wettingPhaseIdx=*/FluidSystem::waterPhaseIdx,
                                            /*nonWettingPhaseIdx=*/FluidSystem::oilPhaseIdx,
                                            /*gasPhaseIdx=*/FluidSystem::gasPhaseIdx>;

public:
    using EclMaterialLawManager = ::Opm::EclMaterialLawManager<Traits>;

    using type = typename EclMaterialLawManager::MaterialLaw;
};

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

// Enable dispersion
template<class TypeTag>
struct EnableDispersion<TypeTag, TTag::FlowBaseProblem>
{ static constexpr bool value = false; };

// disable API tracking
template<class TypeTag>
struct EnableApiTracking<TypeTag, TTag::FlowBaseProblem>
{ static constexpr bool value = false; };

// Use the "velocity module" which uses the Eclipse "NEWTRAN" transmissibilities
template<class TypeTag>
struct FluxModule<TypeTag, TTag::FlowBaseProblem>
{ using type = NewTranFluxModule<TypeTag>; };

// Use the dummy gradient calculator in order not to do unnecessary work.
template<class TypeTag>
struct GradientCalculator<TypeTag, TTag::FlowBaseProblem>
{ using type = DummyGradientCalculator<TypeTag>; };

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

} // namespace Opm::Properties

namespace Opm::Parameters {

#ifdef HAVE_DAMARIS
//! Disable the Damaris HDF5 output by default
template<class TypeTag>
struct EnableDamarisOutput<TypeTag, Properties::TTag::FlowBaseProblem>
{ static constexpr bool value = false; };

// If Damaris is available, write specific variable output in parallel
template<class TypeTag>
struct DamarisOutputHdfCollective<TypeTag, Properties::TTag::FlowBaseProblem>
{ static constexpr bool value = true; };

// Save the reservoir model mesh data to the HDF5 file
// (even if field data HDF5 output is disabled)
template<class TypeTag>
struct DamarisSaveMeshToHdf<TypeTag, Properties::TTag::FlowBaseProblem>
{ static constexpr bool value = false; };

// Save the simulation fields (currently only PRESSURE) variables to HDF5 file
template<class TypeTag>
struct DamarisSaveToHdf<TypeTag, Properties::TTag::FlowBaseProblem>
{ static constexpr bool value = true; };

// Specify path and filename of a Python script to run on each end of iteration output
template<class TypeTag>
struct DamarisPythonScript<TypeTag, Properties::TTag::FlowBaseProblem>
{ static constexpr auto value = ""; };

// Specifiy a Paraview Catalyst in situ visualisation script
// (if Paraview is enabled in Damaris)
template<class TypeTag>
struct DamarisPythonParaviewScript<TypeTag, Properties::TTag::FlowBaseProblem>
{ static constexpr auto value = ""; };

// Specify a unique name for the Damaris simulation (used as prefix to HDF5 filenames)
template<class TypeTag>
struct DamarisSimName<TypeTag, Properties::TTag::FlowBaseProblem>
{ static constexpr auto value = ""; };

// Specify the number of Damaris cores (dc) to create (per-node).
// Must divide into the remaining ranks
// equally, e.g. mpirun -np 16 ... -> (if running on one node)
// The following are allowed:
// 1 dc + 15 sim ranks
// or 2 dc + 14 sim
// or 4 dc + 12 sim
// *not* 3 dc + 13 sim ranks
template<class TypeTag>
struct DamarisDedicatedCores<TypeTag, Properties::TTag::FlowBaseProblem>
{ static constexpr int value = 1; };

// Specify the number of Damaris nodes to create
template<class TypeTag>
struct DamarisDedicatedNodes<TypeTag, Properties::TTag::FlowBaseProblem>
{ static constexpr int value = 0; };

// Specify a name for the Damaris shared memory file
// (a unique name will be created by default)
template<class TypeTag>
struct DamarisSharedMemoryName<TypeTag, Properties::TTag::FlowBaseProblem>
{ static constexpr auto value = "" ; };

// Specify the shared memory file size
template<class TypeTag>
struct DamarisSharedMemorySizeBytes<TypeTag, Properties::TTag::FlowBaseProblem>
{ static constexpr long value = 536870912; }; // 512 MB

// Specify the Damaris log level - if set to debug then log is flushed regularly
template<class TypeTag>
struct DamarisLogLevel<TypeTag, Properties::TTag::FlowBaseProblem>
{ static constexpr auto value = "info"; };

// Specify the dask file jason file that specifies the Dask scheduler etc.
template<class TypeTag>
struct DamarisDaskFile<TypeTag, Properties::TTag::FlowBaseProblem>
{ static constexpr auto value = ""; };

// Specify the the exact variables to be passed through
// to Damaris (must exist in the XML file / intiDamarisXmlFile.cpp)
template<class TypeTag>
struct DamarisLimitVariables<TypeTag, Properties::TTag::FlowBaseProblem>
{ static constexpr auto value = ""; };
#endif

// By default, use single precision for the ECL formated results
template<class TypeTag>
struct EclOutputDoublePrecision<TypeTag, Properties::TTag::FlowBaseProblem>
{ static constexpr bool value = false; };

// If available, write the ECL output in a non-blocking manner
template<class TypeTag>
struct EnableAsyncEclOutput<TypeTag, Properties::TTag::FlowBaseProblem>
{ static constexpr bool value = true; };

// By default, we enable the debugging checks if we're compiled in debug mode
template<class TypeTag>
struct EnableDebuggingChecks<TypeTag, Properties::TTag::FlowBaseProblem>
{ static constexpr bool value = true; };

// Drift compensation is an experimental feature, i.e., systematic errors in the
// conservation quantities are only compensated for
// as default if experimental mode is enabled.
template<class TypeTag>
struct EnableDriftCompensation<TypeTag, Properties::TTag::FlowBaseProblem>
{ static constexpr bool value = true; };

// enable the ECL output by default
template<class TypeTag>
struct EnableEclOutput<TypeTag,Properties::TTag::FlowBaseProblem>
{ static constexpr bool value = true; };

// Write ESMRY file for fast loading of summary data
template<class TypeTag>
struct EnableEsmry<TypeTag, Properties::TTag::FlowBaseProblem>
{ static constexpr bool value = false; };

// Enable gravity
template<class TypeTag>
struct EnableGravity<TypeTag, Properties::TTag::FlowBaseProblem>
{ static constexpr bool value = true; };

// the cache for intensive quantities can be used for ECL problems and also yields a
// decent speedup...
template<class TypeTag>
struct EnableIntensiveQuantityCache<TypeTag, Properties::TTag::FlowBaseProblem>
{ static constexpr bool value = true; };

// the cache for the storage term can also be used and also yields a decent speedup
template<class TypeTag>
struct EnableStorageCache<TypeTag, Properties::TTag::FlowBaseProblem>
{ static constexpr bool value = true; };

// Disable the VTK output by default for this problem ...
template<class TypeTag>
struct EnableVtkOutput<TypeTag, Properties::TTag::FlowBaseProblem>
{ static constexpr bool value = false; };

// only write the solutions for the report steps to disk
template<class TypeTag>
struct EnableWriteAllSolutions<TypeTag, Properties::TTag::FlowBaseProblem>
{ static constexpr bool value = false; };

// The default for the end time of the simulation [s]
//
// By default, stop it after the universe will probably have stopped
// to exist. (the ECL problem will finish the simulation explicitly
// after it simulated the last episode specified in the deck.)
template<class TypeTag>
struct EndTime<TypeTag, Properties::TTag::FlowBaseProblem>
{
    using type = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr type value = 1e100;
};

// By default, use implicit pressure in rock compaction
template<class TypeTag>
struct ExplicitRockCompaction<TypeTag, Properties::TTag::FlowBaseProblem>
{ static constexpr bool value = false; };

// The default for the initial time step size of the simulation [s].
//
// The chosen value means that the size of the first time step is the
// one of the initial episode (if the length of the initial episode is
// not millions of trillions of years, that is...)
template<class TypeTag>
struct InitialTimeStepSize<TypeTag, Properties::TTag::FlowBaseProblem>
{
    using type = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr type value = 3600*24;
};

// the default for the allowed volumetric error for oil per second
template<class TypeTag>
struct NewtonTolerance<TypeTag, Properties::TTag::FlowBaseProblem>
{
    using type = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr type value = 1e-2;
};

// Parameterize equilibration accuracy
template<class TypeTag>
struct NumPressurePointsEquil<TypeTag, Properties::TTag::FlowBaseProblem>
{ static constexpr int value = ParserKeywords::EQLDIMS::DEPTH_NODES_P::defaultValue; };

// The default location for the ECL output files
template<class TypeTag>
struct OutputDir<TypeTag, Properties::TTag::FlowBaseProblem>
{ static constexpr auto value = "."; };

template<class TypeTag>
struct OutputMode<TypeTag, Properties::TTag::FlowBaseProblem>
{ static constexpr auto value = "all"; };

// The frequency of writing restart (*.ers) files. This is the number of time steps
// between writing restart files
template<class TypeTag>
struct RestartWritingInterval<TypeTag, Properties::TTag::FlowBaseProblem>
{ static constexpr int value = 0xffffff; }; // disable

} // namespace Opm::Parameters

#endif // OPM_FLOW_PROBLEM_PROPERTIES_HPP
