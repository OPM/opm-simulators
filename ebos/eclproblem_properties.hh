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
 * \copydoc Opm::EclProblem
 */
#ifndef ECL_PROBLEM_PROPERTIES_HH
#define ECL_PROBLEM_PROPERTIES_HH

#include <ebos/eclbaseaquifermodel.hh>
#include <ebos/eclcpgridvanguard.hh>
#include <ebos/ecldummygradientcalculator.hh>
#include <ebos/eclfluxmodule.hh>
#include <ebos/eclnewtonmethod.hh>
#include <ebos/ecloutputblackoilmodule.hh>
#include <ebos/eclwriter.hh>
#if HAVE_DAMARIS
#include <ebos/damariswriter.hh>
#endif
#include <ebos/FIBlackOilModel.hpp>
#include <ebos/vtkecltracermodule.hh>

#include <opm/input/eclipse/Parser/ParserKeywords/E.hpp>

#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>
#include <opm/material/thermal/EclThermalLawManager.hpp>

#include <opm/models/discretization/ecfv/ecfvdiscretization.hh>
#include <opm/models/utils/propertysystem.hh>

#include <tuple>

namespace Opm {
template <class TypeTag>
class EclProblem;
}

namespace Opm::Properties {

namespace TTag {

struct EclBaseProblem {
  using InheritsFrom = std::tuple<VtkEclTracer, EclOutputBlackOil, EclCpGridVanguard>;
};
}

// The class which deals with ECL wells
template<class TypeTag, class MyTypeTag>
struct EclWellModel {
    using type = UndefinedProperty;
};

// Write all solutions for visualization, not just the ones for the
// report steps...
template<class TypeTag, class MyTypeTag>
struct EnableWriteAllSolutions {
    using type = UndefinedProperty;
};

// The number of time steps skipped between writing two consequtive restart files
template<class TypeTag, class MyTypeTag>
struct RestartWritingInterval {
    using type = UndefinedProperty;
};

// Enable partial compensation of systematic mass losses via the source term of the next time
// step
template<class TypeTag, class MyTypeTag>
struct EclEnableDriftCompensation {
    using type = UndefinedProperty;
};

// Enable the additional checks even if compiled in debug mode (i.e., with the NDEBUG
// macro undefined). Next to a slightly better performance, this also eliminates some
// print statements in debug mode.
template<class TypeTag, class MyTypeTag>
struct EnableDebuggingChecks {
    using type = UndefinedProperty;
};

// if thermal flux boundaries are enabled an effort is made to preserve the initial
// thermal gradient specified via the TEMPVD keyword
template<class TypeTag, class MyTypeTag>
struct EnableThermalFluxBoundaries {
    using type = UndefinedProperty;
};

// Specify whether API tracking should be enabled (replaces PVT regions).
// TODO: This is not yet implemented
template<class TypeTag, class MyTypeTag>
struct EnableApiTracking {
    using type = UndefinedProperty;
};

// The class which deals with ECL aquifers
template<class TypeTag, class MyTypeTag>
struct EclAquiferModel {
    using type = UndefinedProperty;
};

// In experimental mode, decides if the aquifer model should be enabled or not
template<class TypeTag, class MyTypeTag>
struct EclEnableAquifers {
    using type = UndefinedProperty;
};

// time stepping parameters
template<class TypeTag, class MyTypeTag>
struct EclEnableTuning {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct OutputMode {
    using type = UndefinedProperty;
};
// Parameterize equilibration accuracy
template<class TypeTag, class MyTypeTag>
struct NumPressurePointsEquil {
    using type = UndefinedProperty;
};



// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::EclBaseProblem> {
    using type = EclProblem<TypeTag>;
};

template<class TypeTag>
struct Model<TypeTag, TTag::EclBaseProblem> {
    using type = FIBlackOilModel<TypeTag>;
};

// Select the element centered finite volume method as spatial discretization
template<class TypeTag>
struct SpatialDiscretizationSplice<TypeTag, TTag::EclBaseProblem> {
    using type = TTag::EcfvDiscretization;
};

//! for ebos, use automatic differentiation to linearize the system of PDEs
template<class TypeTag>
struct LocalLinearizerSplice<TypeTag, TTag::EclBaseProblem> {
    using type = TTag::AutoDiffLocalLinearizer;
};

template<class TypeTag>
struct BaseDiscretizationType<TypeTag, TTag::EclBaseProblem> {
    using type = FvBaseDiscretizationNoAdapt<TypeTag>;
};

template<class TypeTag>
struct DiscreteFunction<TypeTag, TTag::EclBaseProblem> {
    using BaseDiscretization = FvBaseDiscretization<TypeTag>;
    using type = typename BaseDiscretization::BlockVectorWrapper;
};

template<class TypeTag>
struct GridView<TypeTag, TTag::EclBaseProblem>
{
    using type = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
};

// Set the material law for fluid fluxes
template<class TypeTag>
struct MaterialLaw<TypeTag, TTag::EclBaseProblem>
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
struct SolidEnergyLaw<TypeTag, TTag::EclBaseProblem>
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
struct ThermalConductionLaw<TypeTag, TTag::EclBaseProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

public:
    using EclThermalLawManager = ::Opm::EclThermalLawManager<Scalar, FluidSystem>;

    using type = typename EclThermalLawManager::ThermalConductionLaw;
};

// ebos can use a slightly faster stencil class because it does not need the normals and
// the integration points of intersections
template<class TypeTag>
struct Stencil<TypeTag, TTag::EclBaseProblem>
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
struct EclAquiferModel<TypeTag, TTag::EclBaseProblem> {
    using type = EclBaseAquiferModel<TypeTag>;
};

// Enable aquifers by default in experimental mode
template<class TypeTag>
struct EclEnableAquifers<TypeTag, TTag::EclBaseProblem> {
    static constexpr bool value = true;
};

// Enable gravity
template<class TypeTag>
struct EnableGravity<TypeTag, TTag::EclBaseProblem> {
    static constexpr bool value = true;
};

// Enable diffusion
template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::EclBaseProblem> {
    static constexpr bool value = true;
};

// Enable dispersion
template<class TypeTag>
struct EnableDispersion<TypeTag, TTag::EclBaseProblem> {
    static constexpr bool value = false;
};

// only write the solutions for the report steps to disk
template<class TypeTag>
struct EnableWriteAllSolutions<TypeTag, TTag::EclBaseProblem> {
    static constexpr bool value = false;
};

// disable API tracking
template<class TypeTag>
struct EnableApiTracking<TypeTag, TTag::EclBaseProblem> {
    static constexpr bool value = false;
};

// The default for the end time of the simulation [s]
//
// By default, stop it after the universe will probably have stopped
// to exist. (the ECL problem will finish the simulation explicitly
// after it simulated the last episode specified in the deck.)
template<class TypeTag>
struct EndTime<TypeTag, TTag::EclBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e100;
};

// The default for the initial time step size of the simulation [s].
//
// The chosen value means that the size of the first time step is the
// one of the initial episode (if the length of the initial episode is
// not millions of trillions of years, that is...)
template<class TypeTag>
struct InitialTimeStepSize<TypeTag, TTag::EclBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 3600*24;
};

// the default for the allowed volumetric error for oil per second
template<class TypeTag>
struct NewtonTolerance<TypeTag, TTag::EclBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-2;
};

// the tolerated amount of "incorrect" amount of oil per time step for the complete
// reservoir. this is scaled by the pore volume of the reservoir, i.e., larger reservoirs
// will tolerate larger residuals.
template<class TypeTag>
struct EclNewtonSumTolerance<TypeTag, TTag::EclBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-4;
};

// set the exponent for the volume scaling of the sum tolerance: larger reservoirs can
// tolerate a higher amount of mass lost per time step than smaller ones! since this is
// not linear, we use the cube root of the overall pore volume by default, i.e., the
// value specified by the NewtonSumTolerance parameter is the "incorrect" mass per
// timestep for an reservoir that exhibits 1 m^3 of pore volume. A reservoir with a total
// pore volume of 10^3 m^3 will tolerate 10 times as much.
template<class TypeTag>
struct EclNewtonSumToleranceExponent<TypeTag, TTag::EclBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.0/3.0;
};

// set number of Newton iterations where the volumetric residual is considered for
// convergence
template<class TypeTag>
struct EclNewtonStrictIterations<TypeTag, TTag::EclBaseProblem> {
    static constexpr int value = 8;
};

// set fraction of the pore volume where the volumetric residual may be violated during
// strict Newton iterations
template<class TypeTag>
struct EclNewtonRelaxedVolumeFraction<TypeTag, TTag::EclBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.03;
};

// the maximum volumetric error of a cell in the relaxed region
template<class TypeTag>
struct EclNewtonRelaxedTolerance<TypeTag, TTag::EclBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e9;
};

// Ignore the maximum error mass for early termination of the newton method.
template<class TypeTag>
struct NewtonMaxError<TypeTag, TTag::EclBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 10e9;
};

// set the maximum number of Newton iterations to 14 because the likelyhood that a time
// step succeeds at more than 14 Newton iteration is rather small
template<class TypeTag>
struct NewtonMaxIterations<TypeTag, TTag::EclBaseProblem> {
    static constexpr int value = 14;
};

// also, reduce the target for the "optimum" number of Newton iterations to 6. Note that
// this is only relevant if the time step is reduced from the report step size for some
// reason. (because ebos first tries to do a report step using a single time step.)
template<class TypeTag>
struct NewtonTargetIterations<TypeTag, TTag::EclBaseProblem> {
    static constexpr int value = 6;
};

// Disable the VTK output by default for this problem ...
template<class TypeTag>
struct EnableVtkOutput<TypeTag, TTag::EclBaseProblem> {
    static constexpr bool value = false;
};

// ... but enable the ECL output by default
template<class TypeTag>
struct EnableEclOutput<TypeTag,TTag::EclBaseProblem> {
    static constexpr bool value = true;
};
#ifdef HAVE_DAMARIS
//! Enable the Damaris output by default
template<class TypeTag>
struct EnableDamarisOutput<TypeTag, TTag::EclBaseProblem> {
    static constexpr bool value = false;
};

// If Damaris is available, write specific variable output in parallel
template<class TypeTag>
struct EnableDamarisOutputCollective<TypeTag, TTag::EclBaseProblem> {
    static constexpr bool value = true;
};
#endif
// If available, write the ECL output in a non-blocking manner
template<class TypeTag>
struct EnableAsyncEclOutput<TypeTag, TTag::EclBaseProblem> {
    static constexpr bool value = true;
};

// Write ESMRY file for fast loading of summary data
template<class TypeTag>
struct EnableEsmry<TypeTag, TTag::EclBaseProblem> {
    static constexpr bool value = false;
};

// By default, use single precision for the ECL formated results
template<class TypeTag>
struct EclOutputDoublePrecision<TypeTag, TTag::EclBaseProblem> {
    static constexpr bool value = false;
};

// The default location for the ECL output files
template<class TypeTag>
struct OutputDir<TypeTag, TTag::EclBaseProblem> {
    static constexpr auto value = ".";
};

// the cache for intensive quantities can be used for ECL problems and also yields a
// decent speedup...
template<class TypeTag>
struct EnableIntensiveQuantityCache<TypeTag, TTag::EclBaseProblem> {
    static constexpr bool value = true;
};

// the cache for the storage term can also be used and also yields a decent speedup
template<class TypeTag>
struct EnableStorageCache<TypeTag, TTag::EclBaseProblem> {
    static constexpr bool value = true;
};

// Use the "velocity module" which uses the Eclipse "NEWTRAN" transmissibilities
template<class TypeTag>
struct FluxModule<TypeTag, TTag::EclBaseProblem> {
    using type = EclTransFluxModule<TypeTag>;
};

// Use the dummy gradient calculator in order not to do unnecessary work.
template<class TypeTag>
struct GradientCalculator<TypeTag, TTag::EclBaseProblem> {
    using type = EclDummyGradientCalculator<TypeTag>;
};

// Use a custom Newton-Raphson method class for ebos in order to attain more
// sophisticated update and error computation mechanisms
template<class TypeTag>
struct NewtonMethod<TypeTag, TTag::EclBaseProblem> {
    using type = EclNewtonMethod<TypeTag>;
};

// The frequency of writing restart (*.ers) files. This is the number of time steps
// between writing restart files
template<class TypeTag>
struct RestartWritingInterval<TypeTag, TTag::EclBaseProblem> {
    static constexpr int value = 0xffffff; // disable
};

// Drift compensation is an experimental feature, i.e., systematic errors in the
// conservation quantities are only compensated for
// as default if experimental mode is enabled.
template<class TypeTag>
struct EclEnableDriftCompensation<TypeTag, TTag::EclBaseProblem> {
    static constexpr bool value = true;

};

// By default, we enable the debugging checks if we're compiled in debug mode
template<class TypeTag>
struct EnableDebuggingChecks<TypeTag, TTag::EclBaseProblem> {
    static constexpr bool value = true;
};

// store temperature (but do not conserve energy, as long as EnableEnergy is false)
template<class TypeTag>
struct EnableTemperature<TypeTag, TTag::EclBaseProblem> {
    static constexpr bool value = true;
};

template<class TypeTag>
struct EnableMech<TypeTag, TTag::EclBaseProblem> {
    static constexpr bool value = false;
};
// disable all extensions supported by black oil model. this should not really be
// necessary but it makes things a bit more explicit
template<class TypeTag>
struct EnablePolymer<TypeTag, TTag::EclBaseProblem> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EnableSolvent<TypeTag, TTag::EclBaseProblem> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EnableEnergy<TypeTag, TTag::EclBaseProblem> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EnableFoam<TypeTag, TTag::EclBaseProblem> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EnableExtbo<TypeTag, TTag::EclBaseProblem> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EnableMICP<TypeTag, TTag::EclBaseProblem> {
    static constexpr bool value = false;
};

// disable thermal flux boundaries by default
template<class TypeTag>
struct EnableThermalFluxBoundaries<TypeTag, TTag::EclBaseProblem> {
    static constexpr bool value = false;
};

// By default, simulators derived from the EclBaseProblem are production simulators,
// i.e., experimental features must be explicitly enabled at compile time
template<class TypeTag>
struct EnableExperiments<TypeTag, TTag::EclBaseProblem> {
    static constexpr bool value = false;
};

// set defaults for the time stepping parameters
template<class TypeTag>
struct EclEnableTuning<TypeTag, TTag::EclBaseProblem> {
    static constexpr bool value = false;
};

template<class TypeTag>
struct OutputMode<TypeTag, TTag::EclBaseProblem> {
    static constexpr auto value = "all";
};
// Parameterize equilibration accuracy
template<class TypeTag>
struct NumPressurePointsEquil<TypeTag, TTag::EclBaseProblem> {
    static constexpr int value = ParserKeywords::EQLDIMS::DEPTH_NODES_P::defaultValue;
};

} // namespace Opm::Properties

#endif // ECL_PROBLEM_PROPERTIES_HH
