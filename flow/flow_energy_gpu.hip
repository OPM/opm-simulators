/*
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

// Define making clear that the simulator supports AMG
#define FLOW_SUPPORT_AMG 1

#include <flow/flow_gaswater_energy.hpp>

#include <opm/material/common/ResetLocale.hpp>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/grid/CpGrid.hpp>
#include <opm/simulators/flow/SimulatorFullyImplicitBlackoil.hpp>
#include <opm/simulators/flow/Main.hpp>

#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>

namespace Opm {

// Simple dummy fluid state
class DummyFluidState {
public:
    __host__ __device__ double saturation(unsigned /*phaseIdx*/) const { return 0.5; }
    __host__ __device__ double invB(unsigned /*phaseIdx*/) const { return 1.0; }
    __host__ __device__ double pressure(unsigned /*phaseIdx*/) const { return 1.0e5; }
    __host__ __device__ double Rs() const { return 0.0; }
    __host__ __device__ double Rsw() const { return 0.0; }
    __host__ __device__ double Rv() const { return 0.0; }
    __host__ __device__ double Rvw() const { return 0.0; }
    __host__ __device__ double temperature() const { return 300.0; }
};

// Dummy intensive quantities class
template <class TypeTag>
class DummyIntensiveQuantities {
public:
    using FluidState = DummyFluidState;
    __host__ __device__ const FluidState& fluidState() const { return fluidState_; }
    __host__ __device__ double porosity() const { return 0.2; }
    __host__ __device__ unsigned pvtRegionIndex() const { return 0; }
    __host__ __device__ double mobility(unsigned /*phaseIdx*/, FaceDir::DirEnum /*faceDir*/) const { return 1.0; }
    __host__ __device__ double rockCompTransMultiplier() const { return 1.0; }
    __host__ __device__ double permFactor() const { return 1.0; }
private:
    FluidState fluidState_;
};

template <typename TypeTag>
class DummyFIBlackOilModel
{
public:
    __host__ __device__ explicit DummyFIBlackOilModel(int /* simulator */) {}
    __host__ __device__ void invalidateAndUpdateIntensiveQuantities(unsigned /* timeIdx */) const {}
    __host__ __device__ void invalidateAndUpdateIntensiveQuantitiesOverlap(unsigned /* timeIdx */) const {}
    template <class GridSubDomain>
    __host__ __device__ void invalidateAndUpdateIntensiveQuantities(unsigned /* timeIdx */, const GridSubDomain& /* gridSubDomain */) const {}
    __host__ __device__ void updateFailed() {}
    __host__ __device__ int intensiveQuantities(unsigned /* globalIdx */, unsigned /* timeIdx */) const { return 42; }
protected:
    __host__ __device__ void updateCachedIntQuants(const unsigned /* timeIdx */) const {}
    template <class EMDArg>
    __host__ __device__ void updateCachedIntQuants1(const unsigned /* timeIdx */) const {}
    template <class ...Args>
    __host__ __device__ void updateCachedIntQuantsLoop(const unsigned /* timeIdx */) const {}
    template <class ...Args>
    __host__ __device__ void updateSingleCachedIntQuantUnchecked(const unsigned /* globalIdx */, const unsigned /* timeIdx */) const {}
};


template <class TypeTag>
class DummyBlackOilLocalResidualTPFA
{
public:
    using RateVector = std::array<double, 3>;
    using IntensiveQuantities = int;
    using ElementContext = int;
    using Problem = int;
    struct ResidualNBInfo
    {
        double trans = 0.0;
        double faceArea = 0.0;
        double thpres = 0.0;
        double dZg = 0.0;
        int faceDir = 0;
        double Vin = 0.0;
        double Vex = 0.0;
        double inAlpha = 0.0;
        double outAlpha = 0.0;
        double diffusivity = 0.0;
        double dispersivity = 0.0;
    };

    struct ModuleParams { int dummy = 0; };

    // Dummy storage computation
    template <class LhsEval>
    __host__ __device__ void computeStorage(std::array<LhsEval, 3>& storage,
                        const ElementContext& /* elemCtx */,
                        unsigned /* dofIdx */,
                        unsigned /* timeIdx */) const
    { storage = {0.0, 0.0, 0.0}; }

    template <class LhsEval>
    __host__ __device__ static void computeStorage(std::array<LhsEval, 3>& storage,
                               const IntensiveQuantities& /* intQuants */)
    { storage = {0.0, 0.0, 0.0}; }

    // Dummy flux computation
    __host__ __device__ static void computeFlux(RateVector& flux,
                            RateVector& darcy,
                            const unsigned /* globalIndexIn */,
                            const unsigned /* globalIndexEx */,
                            const IntensiveQuantities& /* intQuantsIn */,
                            const IntensiveQuantities& /* intQuantsEx */,
                            const ResidualNBInfo& /* nbInfo */,
                            const ModuleParams& /* moduleParams */)
    {
        flux = {0.0, 0.0, 0.0};
        darcy = {0.0, 0.0, 0.0};
    }

    __host__ __device__ static void computeFlux(RateVector& flux,
                            const ElementContext& /* elemCtx */,
                            unsigned /* scvfIdx */,
                            unsigned /* timeIdx */)
    { flux = {0.0, 0.0, 0.0}; }

    // Dummy boundary flux computation
    template <class BoundaryConditionData>
    __host__ __device__ static void computeBoundaryFlux(RateVector& bdyFlux,
                                    const Problem& /* problem */,
                                    const BoundaryConditionData& /* bdyInfo */,
                                    const IntensiveQuantities& /* insideIntQuants */,
                                    unsigned /* globalSpaceIdx */)
    { bdyFlux = {0.0, 0.0, 0.0}; }

    template <class BoundaryConditionData>
    __host__ __device__ static void computeBoundaryFluxRate(RateVector& bdyFlux,
                                        const BoundaryConditionData& /* bdyInfo */)
    { bdyFlux = {0.0, 0.0, 0.0}; }

    template <class BoundaryConditionData>
    __host__ __device__ static void computeBoundaryFluxFree(const Problem& /* problem */,
                                        RateVector& bdyFlux,
                                        const BoundaryConditionData& /* bdyInfo */,
                                        const IntensiveQuantities& /* insideIntQuants */,
                                        unsigned /* globalSpaceIdx */)
    { bdyFlux = {0.0, 0.0, 0.0}; }

    template <class BoundaryConditionData>
    __host__ __device__ static void computeBoundaryThermal(const Problem& /* problem */,
                                       RateVector& bdyFlux,
                                       const BoundaryConditionData& /* bdyInfo */,
                                       const IntensiveQuantities& /* insideIntQuants */,
                                       unsigned /* globalSpaceIdx */)
    { bdyFlux = {0.0, 0.0, 0.0}; }

    // Dummy source computation
    __host__ __device__ static void computeSource(RateVector& source,
                              const Problem& /* problem */,
                              const IntensiveQuantities& /* insideIntQuants */,
                              unsigned /* globalSpaceIdex */,
                              unsigned /* timeIdx */)
    { source = {0.0, 0.0, 0.0};}

    __host__ __device__ static void computeSourceDense(RateVector& source,
                                   const Problem& /* problem */,
                                   const IntensiveQuantities& /* insideIntQuants */,
                                   unsigned /* globalSpaceIdex */,
                                   unsigned /* timeIdx */)
    { source = {0.0, 0.0, 0.0}; }

    __host__ __device__ void computeSource(RateVector& source,
                       const ElementContext& /* elemCtx */,
                       unsigned /* dofIdx */,
                       unsigned /* timeIdx */) const
    { assert(false && "DO NOT USE ELEMENT CONTEXT ON GPU"); } // Should be caught in the CPU version of the assembly
};

namespace Properties {
namespace TTag {
struct FlowGasWaterEnergyProblem {
    using InheritsFrom = std::tuple<FlowProblem>;
};
}

template<class TypeTag>
struct Linearizer<TypeTag, TTag::FlowGasWaterEnergyProblem> { using type = TpfaLinearizer<TypeTag>; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::FlowGasWaterEnergyProblem> { using type = BlackOilLocalResidualTPFA<TypeTag>; };

template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::FlowGasWaterEnergyProblem> { static constexpr bool value = true; };

template<class TypeTag>
struct EnableDispersion<TypeTag, TTag::FlowGasWaterEnergyProblem> { static constexpr bool value = true; };

template<class TypeTag>
struct EnableEnergy<TypeTag, TTag::FlowGasWaterEnergyProblem> {
    static constexpr bool value = true;
};

template<class TypeTag>
struct EnableDisgasInWater<TypeTag, TTag::FlowGasWaterEnergyProblem> {
    static constexpr bool value = true;
};

template<class TypeTag>
struct EnableVapwat<TypeTag, TTag::FlowGasWaterEnergyProblem> {
    static constexpr bool value = true;
};

template<class TypeTag>
struct GPUIntensiveQuantities<TypeTag, TTag::FlowGasWaterEnergyProblem> { using type = DummyIntensiveQuantities<TypeTag>; };

template<class TypeTag>
struct GPUFIBlackOilModel<TypeTag, TTag::FlowGasWaterEnergyProblem> { using type = DummyFIBlackOilModel<TypeTag>; };

template<class TypeTag>
struct GPUBlackOilLocalResidualTPFA<TypeTag, TTag::FlowGasWaterEnergyProblem> { using type = DummyBlackOilLocalResidualTPFA<TypeTag>; };

//! The indices required by the model
template<class TypeTag>
struct Indices<TypeTag, TTag::FlowGasWaterEnergyProblem>
{
private:
    // it is unfortunately not possible to simply use 'TypeTag' here because this leads
    // to cyclic definitions of some properties. if this happens the compiler error
    // messages unfortunately are *really* confusing and not really helpful.
    using BaseTypeTag = TTag::FlowProblem;
    using FluidSystem = GetPropType<BaseTypeTag, Properties::FluidSystem>;

public:
    using type = BlackOilTwoPhaseIndices<getPropValue<TypeTag, Properties::EnableSolvent>(),
                                         getPropValue<TypeTag, Properties::EnableExtbo>(),
                                         getPropValue<TypeTag, Properties::EnablePolymer>(),
                                         getPropValue<TypeTag, Properties::EnableEnergy>(),
                                         getPropValue<TypeTag, Properties::EnableFoam>(),
                                         getPropValue<TypeTag, Properties::EnableBrine>(),
                                         /*PVOffset=*/0,
                                         /*disabledCompIdx=*/FluidSystem::oilCompIdx,
                                         getPropValue<TypeTag, Properties::EnableBioeffects>()>;
};
}}

namespace Opm {


// ----------------- Main program -----------------
int flowGasWaterEnergyMain(int argc, char** argv, bool outputCout, bool outputFiles)
{
    // we always want to use the default locale, and thus spare us the trouble
    // with incorrect locale settings.
    resetLocale();

    FlowMain<Properties::TTag::FlowGasWaterEnergyProblem>
        mainfunc {argc, argv, outputCout, outputFiles} ;
    return mainfunc.execute();
}

int flowGasWaterEnergyMainStandalone(int argc, char** argv)
{
    using TypeTag = Properties::TTag::FlowGasWaterEnergyProblem;
    auto mainObject = std::make_unique<::Opm::Main>(argc, argv);
    auto ret = mainObject->runStatic<TypeTag>();
    // Destruct mainObject as the destructor calls MPI_Finalize!
    mainObject.reset();
    return ret;
}

}
