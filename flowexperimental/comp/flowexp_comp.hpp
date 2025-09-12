/*
  Copyright 2024, SINTEF Digital

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
#ifndef FLOWEXP_COMP_HPP
#define FLOWEXP_COMP_HPP

#include <opm/material/constraintsolvers/PTFlash.hpp>
#include <opm/material/fluidsystems/GenericOilGasWaterFluidSystem.hpp>

#include <opm/models/discretization/common/baseauxiliarymodule.hh>
#include <opm/models/ptflash/flashmodel.hh>

#include <opm/simulators/flow/FlowProblemComp.hpp>
#include <opm/simulators/flow/FlowProblemCompProperties.hpp>

#include <opm/simulators/linalg/parallelbicgstabbackend.hh>

#include <flowexperimental/comp/EmptyModel.hpp>
#include <flowexperimental/comp/wells/CompWellModel.hpp>

// // the current code use eclnewtonmethod adding other conditions to proceed_ should do the trick for KA
// // adding linearshe sould be chaning the update_ function in the same class with condition that the error is reduced.
// the trick is to be able to recalculate the residual from here.
// unsure where the timestepping is done from suggestedtime??
// suggestTimeStep is taken from newton solver in problem.limitTimestep
namespace Opm {

template<int numComp, bool EnableWater>
int dispatchFlowExpComp(int argc, char** argv);

}

namespace Opm::Properties {
namespace TTag {

template<int NumComp, bool EnableWater>
struct FlowExpCompProblem {
   using InheritsFrom = std::tuple<FlowBaseProblemComp, FlashModel>;
};

}

template<class TypeTag, int NumComp, bool EnableWater>
struct SparseMatrixAdapter<TypeTag, TTag::FlowExpCompProblem<NumComp, EnableWater>>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    using Block = MatrixBlock<Scalar, numEq, numEq>;

public:
    using type = typename Linear::IstlSparseMatrixAdapter<Block>;
};

#if 0
template<class TypeTag>
struct SolidEnergyLaw<TypeTag, TTag::FlowExpCompProblem>
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
struct ThermalConductionLaw<TypeTag, TTag::FlowExpCompProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

public:
    using EclThermalLawManager = ::Opm::EclThermalLawManager<Scalar, FluidSystem>;

    using type = typename EclThermalLawManager::ThermalConductionLaw;
};


template <class TypeTag>
struct SpatialDiscretizationSplice<TypeTag, TTag::FlowExpCompProblem>
{
    using type = TTag::EcfvDiscretization;
};

template <class TypeTag>
struct LocalLinearizerSplice<TypeTag, TTag::FlowExpCompProblem>
{
    using type = TTag::AutoDiffLocalLinearizer;
};
#endif

// Set the problem property
template <class TypeTag, int NumComp, bool EnableWater>
struct Problem<TypeTag, TTag::FlowExpCompProblem<NumComp, EnableWater>>
{
    using type = FlowProblemComp<TypeTag>;
};

template<class TypeTag, int NumComp, bool EnableWater>
struct AquiferModel<TypeTag, TTag::FlowExpCompProblem<NumComp, EnableWater>> {
    using type = EmptyModel<TypeTag>;
};

template<class TypeTag, int NumComp, bool EnableWater>
struct WellModel<TypeTag, TTag::FlowExpCompProblem<NumComp, EnableWater>> {
    using type = CompWellModel<TypeTag>;
};

template<class TypeTag, int NumComp, bool EnableWater>
struct TracerModel<TypeTag, TTag::FlowExpCompProblem<NumComp, EnableWater>> {
    using type = EmptyModel<TypeTag>;
};


template <class TypeTag, int NumComp, bool EnableWater>
struct FlashSolver<TypeTag, TTag::FlowExpCompProblem<NumComp, EnableWater>> {
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;

public:
    using type = Opm::PTFlash<Scalar, FluidSystem>;
};


template <class TypeTag, class MyTypeTag>
struct NumComp { using type = UndefinedProperty; };

// TODO: this is unfortunate, have to check why we need to hard-code it
template <class TypeTag, int NumComp_, bool EnableWater_>
struct NumComp<TypeTag, TTag::FlowExpCompProblem<NumComp_, EnableWater_>> {
    static constexpr int value = NumComp_;
};

template <class TypeTag, class MyTypeTag>
struct EnableDummyWater { using type = UndefinedProperty; };

template <class TypeTag, int NumComp_, bool EnableWater_>
struct EnableDummyWater<TypeTag, TTag::FlowExpCompProblem<NumComp_, EnableWater_>> {
    static constexpr bool value = EnableWater_;
};
#if 0
struct Temperature { using type = UndefinedProperty; };

 template <class TypeTag>
 struct Temperature<TypeTag, TTag::FlowExpCompProblem> {
     using type = GetPropType<TypeTag, Scalar>;
     static constexpr type value = 423.25;
 };
#endif

template <class TypeTag, int NumComp_, bool EnableWater_>
struct FluidSystem<TypeTag, TTag::FlowExpCompProblem<NumComp_, EnableWater_>>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr int num_comp = getPropValue<TypeTag, Properties::NumComp>();
    static constexpr bool enable_water = getPropValue<TypeTag, Properties::EnableDummyWater>();

public:
    using type = Opm::GenericOilGasWaterFluidSystem<Scalar, num_comp, enable_water>;
};
template<class TypeTag, int NumComp, bool EnableWater>
struct EnableMech<TypeTag, TTag::FlowExpCompProblem<NumComp, EnableWater>> {
    static constexpr bool value = false;
};

template<class TypeTag, int NumComp, bool EnableWater>
struct EnableDisgasInWater<TypeTag, TTag::FlowExpCompProblem<NumComp, EnableWater>> { 
    static constexpr bool value = false; 
};

template<class TypeTag, int NumComp, bool EnableWater>
struct Stencil<TypeTag, TTag::FlowExpCompProblem<NumComp, EnableWater>>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;

public:
    using type = EcfvStencil<Scalar, GridView>;
};

template<class TypeTag, int NumComp, bool EnableWater>
struct EnableApiTracking<TypeTag, TTag::FlowExpCompProblem<NumComp, EnableWater>> {
    static constexpr bool value = false;
};

template<class TypeTag, int NumComp, bool EnableWater>
struct EnableTemperature<TypeTag, TTag::FlowExpCompProblem<NumComp, EnableWater>> {
    static constexpr bool value = false;
};

template<class TypeTag, int NumComp, bool EnableWater>
struct EnableSaltPrecipitation<TypeTag, TTag::FlowExpCompProblem<NumComp, EnableWater>> {
    static constexpr bool value = false;
};
template<class TypeTag, int NumComp, bool EnableWater>
struct EnablePolymerMW<TypeTag, TTag::FlowExpCompProblem<NumComp, EnableWater>> {
    static constexpr bool value = false;
};

template<class TypeTag, int NumComp, bool EnableWater>
struct EnablePolymer<TypeTag, TTag::FlowExpCompProblem<NumComp, EnableWater>> {
    static constexpr bool value = false;
};

template<class TypeTag, int NumComp, bool EnableWater>
struct EnableDispersion<TypeTag, TTag::FlowExpCompProblem<NumComp, EnableWater>> {
    static constexpr bool value = false;
};

template<class TypeTag, int NumComp, bool EnableWater>
struct EnableBrine<TypeTag, TTag::FlowExpCompProblem<NumComp, EnableWater>> {
    static constexpr bool value = false;
};
template<class TypeTag, int NumComp, bool EnableWater>
struct EnableVapwat<TypeTag, TTag::FlowExpCompProblem<NumComp, EnableWater>> {
    static constexpr bool value = false;
};

template<class TypeTag, int NumComp, bool EnableWater>
struct EnableSolvent<TypeTag, TTag::FlowExpCompProblem<NumComp, EnableWater>> {
    static constexpr bool value = false;
};
template<class TypeTag, int NumComp, bool EnableWater>
struct EnableEnergy<TypeTag, TTag::FlowExpCompProblem<NumComp, EnableWater>> {
    static constexpr bool value = false;
};
template<class TypeTag, int NumComp, bool EnableWater>
struct EnableFoam<TypeTag, TTag::FlowExpCompProblem<NumComp, EnableWater>> {
    static constexpr bool value = false;
};
template<class TypeTag, int NumComp, bool EnableWater>
struct EnableExtbo<TypeTag, TTag::FlowExpCompProblem<NumComp, EnableWater>> {
    static constexpr bool value = false;
};
template<class TypeTag, int NumComp, bool EnableWater>
struct EnableBioeffects<TypeTag, TTag::FlowExpCompProblem<NumComp, EnableWater>> {
    static constexpr bool value = false;
};

// disable thermal flux boundaries by default
#if 0
template<class TypeTag>
struct EnableThermalFluxBoundaries<TypeTag, TTag::FlowExpCompProblem> {
    static constexpr bool value = false;
};
#endif

} // namespace Opm::Properties

#endif
