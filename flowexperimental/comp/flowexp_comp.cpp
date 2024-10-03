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

#include <dune/common/parallel/mpihelper.hh>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <opm/models/utils/start.hh>
#include <opm/material/constraintsolvers/PTFlash.hpp>
#include "../FlowExpNewtonMethod.hpp"
#include <opm/models/ptflash/flashmodel.hh>
#include <opm/material/fluidsystems/GenericOilGasFluidSystem.hpp>
#include <opm/models/discretization/common/baseauxiliarymodule.hh>

#include <opm/simulators/flow/FlowProblemComp.hpp>
#include <opm/simulators/flow/FlowProblemCompProperties.hpp>
// TODO: not understanding why we need FlowGenericProblem here
#include <opm/simulators/flow/FlowGenericProblem.hpp>
#include <opm/simulators/flow/FlowGenericProblem_impl.hpp>
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
        void addToSource(RateVector& /*rates*/, const Context& /*context*/,
                         unsigned /*spaceIdx*/, unsigned /*timeIdx*/) const {}
        template<class RateVector>
        void addToSource(RateVector& /*rates*/, unsigned /*globalSpaceIdx*/,
                         unsigned /*timeIdx*/) const {}
        void endIteration()const{};
        void endTimeStep(){};
        void endEpisode(){};
        void applyInitial(){};
        template<class RateType>
        void computeTotalRatesForDof(RateType& /*rate*/, unsigned /*globalIdx*/) const{};
    };

}


namespace Opm::Properties {

   namespace TTag {
   template<int NumComp>
   struct FlowExpCompProblem {
       using InheritsFrom = std::tuple<FlowBaseProblemComp, FlashModel>;
   };
   }

    template<class TypeTag, int NumComp>
    struct SparseMatrixAdapter<TypeTag, TTag::FlowExpCompProblem<NumComp>>
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
template <class TypeTag, int NumComp>
struct Problem<TypeTag, TTag::FlowExpCompProblem<NumComp>>
{
    using type = FlowProblemComp<TypeTag>;
};

template<class TypeTag, int NumComp>
struct AquiferModel<TypeTag, TTag::FlowExpCompProblem<NumComp>> {
    using type = EmptyModel<TypeTag>;
};

template<class TypeTag, int NumComp>
struct WellModel<TypeTag, TTag::FlowExpCompProblem<NumComp>> {
    using type = EmptyModel<TypeTag>;
};

template<class TypeTag, int NumComp>
struct TracerModel<TypeTag, TTag::FlowExpCompProblem<NumComp>> {
    using type = EmptyModel<TypeTag>;
};


template <class TypeTag, int NumComp>
struct FlashSolver<TypeTag, TTag::FlowExpCompProblem<NumComp>> {
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
template <class TypeTag, int NumComp_>
struct NumComp<TypeTag, TTag::FlowExpCompProblem<NumComp_>> {
    static constexpr int value = NumComp_;
};
#if 0
struct Temperature { using type = UndefinedProperty; };

 template <class TypeTag>
 struct Temperature<TypeTag, TTag::FlowExpCompProblem> {
     using type = GetPropType<TypeTag, Scalar>;
     static constexpr type value = 423.25;
 };
#endif

template <class TypeTag, int NumComp_>
struct FluidSystem<TypeTag, TTag::FlowExpCompProblem<NumComp_>>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr int num_comp = getPropValue<TypeTag, Properties::NumComp>();

public:
    using type = Opm::GenericOilGasFluidSystem<Scalar, num_comp>;
};
template<class TypeTag, int NumComp>
struct EnableMech<TypeTag, TTag::FlowExpCompProblem<NumComp>> {
    static constexpr bool value = false;
};

template<class TypeTag, int NumComp>
struct EnableDisgasInWater<TypeTag, TTag::FlowExpCompProblem<NumComp>> { static constexpr bool value = false; };

template<class TypeTag, int NumComp>
struct Stencil<TypeTag, TTag::FlowExpCompProblem<NumComp>>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;

public:
    using type = EcfvStencil<Scalar, GridView>;
};

template<class TypeTag, int NumComp>
struct EnableApiTracking<TypeTag, TTag::FlowExpCompProblem<NumComp>> {
    static constexpr bool value = false;
};

template<class TypeTag, int NumComp>
struct EnableTemperature<TypeTag, TTag::FlowExpCompProblem<NumComp>> {
    static constexpr bool value = false;
};

template<class TypeTag, int NumComp>
struct EnableSaltPrecipitation<TypeTag, TTag::FlowExpCompProblem<NumComp>> {
    static constexpr bool value = false;
};
template<class TypeTag, int NumComp>
struct EnablePolymerMW<TypeTag, TTag::FlowExpCompProblem<NumComp>> {
    static constexpr bool value = false;
};

template<class TypeTag, int NumComp>
struct EnablePolymer<TypeTag, TTag::FlowExpCompProblem<NumComp>> {
    static constexpr bool value = false;
};

template<class TypeTag, int NumComp>
struct EnableDispersion<TypeTag, TTag::FlowExpCompProblem<NumComp>> {
    static constexpr bool value = false;
};

template<class TypeTag, int NumComp>
struct EnableBrine<TypeTag, TTag::FlowExpCompProblem<NumComp>> {
    static constexpr bool value = false;
};
template<class TypeTag, int NumComp>
struct EnableVapwat<TypeTag, TTag::FlowExpCompProblem<NumComp>> {
    static constexpr bool value = false;
};

template<class TypeTag, int NumComp>
struct EnableSolvent<TypeTag, TTag::FlowExpCompProblem<NumComp>> {
    static constexpr bool value = false;
};
template<class TypeTag, int NumComp>
struct EnableEnergy<TypeTag, TTag::FlowExpCompProblem<NumComp>> {
    static constexpr bool value = false;
};
template<class TypeTag, int NumComp>
struct EnableFoam<TypeTag, TTag::FlowExpCompProblem<NumComp>> {
    static constexpr bool value = false;
};
template<class TypeTag, int NumComp>
struct EnableExtbo<TypeTag, TTag::FlowExpCompProblem<NumComp>> {
    static constexpr bool value = false;
};
template<class TypeTag, int NumComp>
struct EnableMICP<TypeTag, TTag::FlowExpCompProblem<NumComp>> {
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

//! @brief Runs the simulator with the correct number of components
//!
//! This function will compile-time recurse from numComponentsCompileTime
//! down to (inclusive) minimumNumberOfComponents and check if
//!
//!     numComponentsCompileTime == numComponentsRuntime
//!
//! and if this is the case, it will star the simulator with numComponentsRuntime.
//!
//! @tparam minimumNumberOfComponents the minimum number of components supported
//! @tparam numCompnentsCompileTime the maximum number of components supported
//! @param numComponentsRuntime the number of components found the in the deck file
//! @param argc argc argument from CLI
//! @param argv argv argument from CLI
//!
//! @return error code
template <int minimumNumberOfComponents, int numComponentsCompileTime>
int
startSimulationComponents(int numComponentsRuntime, int argc, char** argv)
{
    OPM_ERROR_IF(numComponentsCompileTime < numComponentsRuntime || numComponentsRuntime < minimumNumberOfComponents,
                 fmt::format("Deck has {} components, not supported. We support a maximum of {} components, "
                             "and a minimum of {}.",
                             numComponentsRuntime,
                             numComponentsCompileTime,
                             minimumNumberOfComponents));

    if (numComponentsCompileTime == numComponentsRuntime) {
        return Opm::start<Opm::Properties::TTag::FlowExpCompProblem<numComponentsCompileTime>>(argc, argv, false);
    }
    if constexpr (numComponentsCompileTime > minimumNumberOfComponents) {
        return startSimulationComponents<minimumNumberOfComponents, numComponentsCompileTime - 1>(numComponentsRuntime, argc, argv);
    }
    // It will never actually reach this, but the compiler does not seem to realize, so keeping
    // this to avoid warnings.
    return EXIT_FAILURE;
}
int
main(int argc, char** argv)
{
    using TypeTag = Opm::Properties::TTag::FlowExpCompProblem<0>;
    Opm::registerEclTimeSteppingParameters<double>();

    // This is a bit cumbersome, but we need to read the input file
    // first to figure out the number of components (I think) in order
    // to select the correct type tag.
    //
    // TODO: Do a more dynamic dispatch approach similar to the normal
    //       flow application
    auto comm = Dune::MPIHelper::instance(argc, argv).getCommunication();
    auto commPtr = std::make_unique<decltype(comm)>(comm);
    Opm::setupParameters_<TypeTag>(argc, const_cast<const char**>(argv), true);

    auto inputFilename
        = Opm::FlowGenericVanguard::canonicalDeckPath(Opm::Parameters::Get<Opm::Parameters::EclDeckFileName>());
    Opm::FlowGenericVanguard::setCommunication(std::move(commPtr));
    Opm::FlowGenericVanguard::readDeck(inputFilename);
    Opm::FlowGenericVanguard vanguard;
    const auto numComps = vanguard.eclState().compositionalConfig().numComps();
    return startSimulationComponents<2, 7>(numComps, argc, argv);
}

