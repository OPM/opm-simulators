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
#include "config.h"
#include <opm/models/utils/start.hh>
#include <opm/simulators/flow/FlowProblemComp.hpp>
#include <opm/material/constraintsolvers/PTFlash.hpp>
#include "../FlowExpNewtonMethod.hpp"
#include <opm/simulators/flow/FlowProblemComp.hpp>
#include <opm/models/ptflash/flashmodel.hh>
#include <opm/material/fluidsystems/GenericOilGasFluidSystem.hpp>
// #include <opm/simulators/flow/Main.hpp>
// #include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>
// #include <flowexperimental/blackoilintensivequantitiessimple.hh>
// #include "BlackOilModelFvNoCache.hpp"
// #include "co2ptflowproblem.hh"
// #include <tests/problems/co2ptflashproblem.hh>
#include <opm/models/discretization/common/baseauxiliarymodule.hh>

#include <opm/simulators/flow/FlowGenericProblem.hpp>
#include <opm/simulators/flow/FlowGenericProblem_impl.hpp>

#include <opm/simulators/flow/FlowProblemCompProperties.hpp>

#include <opm/simulators/linalg/parallelbicgstabbackend.hh>

// // the current code use eclnewtonmethod adding other conditions to proceed_ should do the trick for KA
// // adding linearshe sould be chaning the update_ function in the same class with condition that the error is reduced.
// the trick is to be able to recalculate the residual from here.
// unsure where the timestepping is done from suggestedtime??
// suggestTimeStep is taken from newton solver in problem.limitTimestep
namespace Opm{
    template<typename TypeTag>
    class EmptyModel : public BaseAuxiliaryModule<TypeTag>
    {
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        using GridView = GetPropType<TypeTag, Properties::GridView>;
        using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
        using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    public:
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        EmptyModel(Simulator& /*simulator*/){
        };
        void init(){}
        template<class Something>
        void init(Something /*A*/){}
        void prepareTracerBatches(){};
        using NeighborSet = std::set<unsigned>;
        void linearize(SparseMatrixAdapter& /*matrix*/, GlobalEqVector& /*residual*/){};
        unsigned numDofs() const{return 0;};
        void addNeighbors(std::vector<NeighborSet>& /*neighbors*/) const{};
        //void applyInitial(){};
        void initialSolutionApplied(){};
        //void initFromRestart(const data::Aquifers& aquiferSoln);
        template <class Restarter>
        void serialize(Restarter& /*res*/){};

        template <class Restarter>
        void deserialize(Restarter& /*res*/){};

        void beginEpisode(){};
        void beginTimeStep(){};
        void beginIteration(){};
        // add the water rate due to aquifers to the source term.
        template<class RateVector, class Context>
        void addToSource(RateVector& rates, const Context& context, unsigned spaceIdx, unsigned timeIdx) const{};
        template<class RateVector>
        void addToSource(RateVector& rates, unsigned globalSpaceIdx, unsigned timeIdx) const{};
        void endIteration()const{};
        void endTimeStep(){};
        void endEpisode(){};
        void applyInitial(){};
        template<class RateType>
        void computeTotalRatesForDof(RateType& /*rate*/, unsigned /*globalIdx*/) const{};
    };

}


namespace Opm::Properties {

    template<class TypeTag, class MyTypeTag>
    struct EnableTerminalOutput {
        using type = UndefinedProperty;
    };

   namespace TTag {
   struct FlowExpCompProblem {
       using InheritsFrom = std::tuple<FlowBaseProblemComp, FlashModel>;
   };
   }
#if 0
    template<class TypeTag, class MyTypeTag>
    struct ExpliciteRockCompaction{
        using type = UndefinedProperty;
    };
#endif

#if 0
    template<class TypeTag>
    struct MaterialLaw<TypeTag, TTag::FlowExpCompProblem>
    {
    private:
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        using Indices = GetPropType<TypeTag, Properties::Indices>;

        using Traits = ThreePhaseMaterialTraits<Scalar,
                                                /*wettingPhaseIdx=*/0,
                                                /*nonWettingPhaseIdx=*/1,
                                                /*gasPhaseIdx=*/2>;

    public:
        using EclMaterialLawManager = ::Opm::EclMaterialLawManager<Traits>;
        //using EclMaterialLawManager = ::Opm::EclMaterialLawManager<Traits>;

        using type = typename EclMaterialLawManager::MaterialLaw;
    };
#endif
    template<class TypeTag>
    struct SparseMatrixAdapter<TypeTag, TTag::FlowExpCompProblem>
    {
    private:
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
        using Block = MatrixBlock<Scalar, numEq, numEq>;

    public:
        using type = typename Linear::IstlSparseMatrixAdapter<Block>;
    };

#if 1
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
#endif

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

// Set the problem property
template <class TypeTag>
struct Problem<TypeTag, TTag::FlowExpCompProblem>
{
    // using type = EbosProblemComp<TypeTag>;
    using type = FlowProblemComp<TypeTag>;
};

template<class TypeTag>
struct AquiferModel<TypeTag, TTag::FlowExpCompProblem> {
    using type = EmptyModel<TypeTag>;
};

template<class TypeTag>
struct WellModel<TypeTag, TTag::FlowExpCompProblem> {
    using type = EmptyModel<TypeTag>;
};

template<class TypeTag>
struct TracerModelDef<TypeTag, TTag::FlowExpCompProblem> {
    using type = EmptyModel<TypeTag>;
};


template <class TypeTag>
struct FlashSolver<TypeTag, TTag::FlowExpCompProblem> {
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
template <class TypeTag>
struct NumComp<TypeTag, TTag::FlowExpCompProblem> {
    static constexpr int value = 3;
};

// set the defaults for the problem specific properties
// TODO: should it be here?
template<class TypeTag, class MyTypeTag>
struct Temperature { using type = UndefinedProperty; };

 template <class TypeTag>
 struct Temperature<TypeTag, TTag::FlowExpCompProblem> {
     using type = GetPropType<TypeTag, Scalar>;
     static constexpr type value = 423.25;//TODO
 };

/* template <class TypeTag>
struct SimulationName<TypeTag, TTag::FlowExpCompProblem> {
    static constexpr auto value = "co2_ptflash";
}; */


template <class TypeTag>
struct FluidSystem<TypeTag, TTag::FlowExpCompProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr int num_comp = getPropValue<TypeTag, Properties::NumComp>();

public:
    using type = Opm::GenericOilGasFluidSystem<Scalar, num_comp>;
};
template<class TypeTag>
struct EnableMech<TypeTag, TTag::FlowExpCompProblem> {
    static constexpr bool value = false;
};

template<class TypeTag>
struct EnableDisgasInWater<TypeTag, TTag::FlowExpCompProblem> { static constexpr bool value = false; };

template<class TypeTag>
struct Stencil<TypeTag, TTag::FlowExpCompProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;

public:
    using type = EcfvStencil<Scalar, GridView>;
};


template<class TypeTag>
struct EnableApiTracking<TypeTag, TTag::FlowExpCompProblem> {
    static constexpr bool value = false;
};

template<class TypeTag>
struct EnableTemperature<TypeTag, TTag::FlowExpCompProblem> {
    static constexpr bool value = false;
};

template<class TypeTag>
struct EnableSaltPrecipitation<TypeTag, TTag::FlowExpCompProblem> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EnablePolymerMW<TypeTag, TTag::FlowExpCompProblem> {
    static constexpr bool value = false;
};



template<class TypeTag>
struct EnablePolymer<TypeTag, TTag::FlowExpCompProblem> {
    static constexpr bool value = false;
};


template<class TypeTag>
struct EnableDispersion<TypeTag, TTag::FlowExpCompProblem> {
    static constexpr bool value = false;
};

template<class TypeTag>
struct EnableBrine<TypeTag, TTag::FlowExpCompProblem> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EnableVapwat<TypeTag, TTag::FlowExpCompProblem> {
    static constexpr bool value = false;
};

template<class TypeTag>
struct EnableSolvent<TypeTag, TTag::FlowExpCompProblem> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EnableEnergy<TypeTag, TTag::FlowExpCompProblem> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EnableFoam<TypeTag, TTag::FlowExpCompProblem> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EnableExtbo<TypeTag, TTag::FlowExpCompProblem> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EnableMICP<TypeTag, TTag::FlowExpCompProblem> {
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


int main(int argc, char** argv)
{
    //using TypeTag = Opm::Properties::TTag::EclFlowProblemEbos;
    using TypeTag = Opm::Properties::TTag::FlowExpCompProblem;
    Opm::registerEclTimeSteppingParameters<double>();
    return Opm::start<TypeTag>(argc, argv);
}
