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

#include <opm/common/utility/gpuDecorators.hpp>

#include <opm/material/common/ResetLocale.hpp>
#include <opm/material/thermal/EnergyModuleType.hpp>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/grid/CpGrid.hpp>
#include <opm/simulators/flow/SimulatorFullyImplicitBlackoil.hpp>
#include <opm/simulators/flow/Main.hpp>

#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>

namespace Opm {

// Simple dummy fluid state
class FakeFluidState {
public:
    OPM_HOST_DEVICE FakeFluidState() {}
    OPM_HOST_DEVICE double saturation(unsigned /*phaseIdx*/) const { return 0.5; }
    OPM_HOST_DEVICE double invB(unsigned /*phaseIdx*/) const { return 1.0; }
    OPM_HOST_DEVICE double pressure(unsigned /*phaseIdx*/) const { return 1.0e5; }
    OPM_HOST_DEVICE double Rs() const { return 0.0; }
    OPM_HOST_DEVICE double Rsw() const { return 0.0; }
    OPM_HOST_DEVICE double Rv() const { return 0.0; }
    OPM_HOST_DEVICE double Rvw() const { return 0.0; }
    OPM_HOST_DEVICE double temperature() const { return 300.0; }
};

// Fake intensive quantities class
template <class TypeTag>
class FakeIntensiveQuantities {
public:
    using FluidState = FakeFluidState;
    OPM_HOST_DEVICE FakeIntensiveQuantities() : fluidState_() {}
    OPM_HOST_DEVICE const FluidState& fluidState() const { return fluidState_; }
    OPM_HOST_DEVICE double porosity() const { return 0.2; }
    OPM_HOST_DEVICE unsigned pvtRegionIndex() const { return 0; }
    OPM_HOST_DEVICE double mobility(unsigned /*phaseIdx*/, FaceDir::DirEnum /*faceDir*/) const { return 1.0; }
    OPM_HOST_DEVICE double rockCompTransMultiplier() const { return 1.0; }
    OPM_HOST_DEVICE double permFactor() const { return 1.0; }
private:
    FluidState fluidState_;
};

template <typename TypeTag>
class FakeFIBlackOilModel
{
public:
    OPM_HOST_DEVICE FakeFIBlackOilModel() {}
    OPM_HOST_DEVICE explicit FakeFIBlackOilModel(int /* simulator */) {}
    OPM_HOST_DEVICE void invalidateAndUpdateIntensiveQuantities(unsigned /* timeIdx */) const {}
    OPM_HOST_DEVICE void invalidateAndUpdateIntensiveQuantitiesOverlap(unsigned /* timeIdx */) const {}
    template <class GridSubDomain>
    OPM_HOST_DEVICE void invalidateAndUpdateIntensiveQuantities(unsigned /* timeIdx */, const GridSubDomain& /* gridSubDomain */) const {}
    OPM_HOST_DEVICE void updateFailed() {}
    OPM_HOST_DEVICE FakeIntensiveQuantities<TypeTag> intensiveQuantities(unsigned /* globalIdx */, unsigned /* timeIdx */) const { return FakeIntensiveQuantities<TypeTag>{}; }
protected:
    OPM_HOST_DEVICE void updateCachedIntQuants(const unsigned /* timeIdx */) const {}
    template <class EMDArg>
    OPM_HOST_DEVICE void updateCachedIntQuants1(const unsigned /* timeIdx */) const {}
    template <class ...Args>
    OPM_HOST_DEVICE void updateCachedIntQuantsLoop(const unsigned /* timeIdx */) const {}
    template <class ...Args>
    OPM_HOST_DEVICE void updateSingleCachedIntQuantUnchecked(const unsigned /* globalIdx */, const unsigned /* timeIdx */) const {}
};


template <class TypeTag>
class FakeBlackOilLocalResidualTPFA
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

    // Fake storage computation
    template <class LhsEval>
    OPM_HOST_DEVICE void computeStorage(std::array<LhsEval, 3>& storage,
                        const ElementContext& /* elemCtx */,
                        unsigned /* dofIdx */,
                        unsigned /* timeIdx */) const
    { storage = {0.0, 0.0, 0.0}; }

    template <class LhsEval>
    OPM_HOST_DEVICE static void computeStorage(std::array<LhsEval, 3>& storage,
                               const IntensiveQuantities& /* intQuants */)
    { storage = {0.0, 0.0, 0.0}; }

    // Fake flux computation
    OPM_HOST_DEVICE static void computeFlux(RateVector& flux,
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

    OPM_HOST_DEVICE static void computeFlux(RateVector& flux,
                            const ElementContext& /* elemCtx */,
                            unsigned /* scvfIdx */,
                            unsigned /* timeIdx */)
    { flux = {0.0, 0.0, 0.0}; }

    // Fake boundary flux computation
    template <class BoundaryConditionData>
    OPM_HOST_DEVICE static void computeBoundaryFlux(RateVector& bdyFlux,
                                    const Problem& /* problem */,
                                    const BoundaryConditionData& /* bdyInfo */,
                                    const IntensiveQuantities& /* insideIntQuants */,
                                    unsigned /* globalSpaceIdx */)
    { bdyFlux = {0.0, 0.0, 0.0}; }

    template <class BoundaryConditionData>
    OPM_HOST_DEVICE static void computeBoundaryFluxRate(RateVector& bdyFlux,
                                        const BoundaryConditionData& /* bdyInfo */)
    { bdyFlux = {0.0, 0.0, 0.0}; }

    template <class BoundaryConditionData>
    OPM_HOST_DEVICE static void computeBoundaryFluxFree(const Problem& /* problem */,
                                        RateVector& bdyFlux,
                                        const BoundaryConditionData& /* bdyInfo */,
                                        const IntensiveQuantities& /* insideIntQuants */,
                                        unsigned /* globalSpaceIdx */)
    { bdyFlux = {0.0, 0.0, 0.0}; }

    template <class BoundaryConditionData>
    OPM_HOST_DEVICE static void computeBoundaryThermal(const Problem& /* problem */,
                                       RateVector& bdyFlux,
                                       const BoundaryConditionData& /* bdyInfo */,
                                       const IntensiveQuantities& /* insideIntQuants */,
                                       unsigned /* globalSpaceIdx */)
    { bdyFlux = {0.0, 0.0, 0.0}; }

    // Fake source computation
    OPM_HOST_DEVICE static void computeSource(RateVector& source,
                              const Problem& /* problem */,
                              const IntensiveQuantities& /* insideIntQuants */,
                              unsigned /* globalSpaceIdex */,
                              unsigned /* timeIdx */)
    { source = {0.0, 0.0, 0.0};}

    OPM_HOST_DEVICE static void computeSourceDense(RateVector& source,
                                   const Problem& /* problem */,
                                   const IntensiveQuantities& /* insideIntQuants */,
                                   unsigned /* globalSpaceIdex */,
                                   unsigned /* timeIdx */)
    { source = {0.0, 0.0, 0.0}; }

    OPM_HOST_DEVICE void computeSource(RateVector& source,
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

// After a rebase the lines above were not sufficient to run a simulation with thermal effects
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

template<class TypeTag>
struct TrivialIntensiveQuantities<TypeTag, TTag::FlowGasWaterEnergyProblem> { using type = FakeIntensiveQuantities<TypeTag>; };

template<class TypeTag>
struct TrivialFIBlackOilModel<TypeTag, TTag::FlowGasWaterEnergyProblem> { using type = FakeFIBlackOilModel<TypeTag>; };

template<class TypeTag>
struct TrivialBlackOilLocalResidualTPFA<TypeTag, TTag::FlowGasWaterEnergyProblem> { using type = FakeBlackOilLocalResidualTPFA<TypeTag>; };

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
