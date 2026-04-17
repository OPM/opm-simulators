/*
  Copyright 2025 Equinor ASA

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
#include <config.h>

// Tricks for older versions of ROCm that do not support alugrid
// or dune fem. Note that these changes are only valid for this test file.
#ifdef HAVE_DUNE_ALUGRID
#undef HAVE_DUNE_ALUGRID
#endif
#define HAVE_DUNE_ALUGRID 0
#ifdef HAVE_DUNE_FEM
#undef HAVE_DUNE_FEM
#endif
#define HAVE_DUNE_FEM 0

#include <opm/common/utility/gpuDecorators.hpp>
#include <opm/material/common/ResetLocale.hpp>
#include <opm/material/fluidmatrixinteractions/EclDefaultMaterial.hpp>
#define HAVE_ECL_INPUT 1


#define BOOST_TEST_MODULE TestBlackOilIntensiveQuantitiesGPU

#include <boost/test/unit_test.hpp>
#include <cuda.h>
#include <cuda_runtime.h>
#include <opm/common/ErrorMacros.hpp>
#include <opm/models/blackoil/blackoilmodel.hh>
#include <opm/models/blackoil/blackoilprimaryvariables.hh>
#include <opm/models/discretization/common/fvbaseprimaryvariables.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>
#include <opm/models/io/dgfvanguard.hh>
#include <opm/models/utils/start.hh>
#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>
#include <opm/simulators/linalg/gpuistl/gpu_smart_pointer.hpp>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Python/Python.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystemNonStatic.hpp>
#include <opm/models/blackoil/blackoilintensivequantities.hh>
#include <opm/models/nonlinear/newtonmethodparams.hpp>
#include <opm/simulators/flow/equil/EquilibrationHelpers.hpp>
#include <opm/simulators/flow/equil/InitStateEquil.hpp>
#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/models/utils/simulator.hh>
#include <opm/simulators/flow/FlowProblemBlackoil.hpp>
#include <opm/simulators/flow/FlowProblemBlackoilProperties.hpp>
#include <opm/simulators/utils/moduleVersion.hpp>
#include <opm/simulators/flow/BlackoilModelParameters.hpp>
#include <opm/simulators/flow/FlowGenericVanguard.hpp>
#include <opm/simulators/linalg/parallelbicgstabbackend.hh>
#include <opm/simulators/wells/BlackoilWellModel.hpp>
#include <opm/models/discretization/common/fvbaseelementcontextgpu.hh>

#include <opm/material/fluidmatrixinteractions/EclMaterialLawTwoPhaseTypes.hpp>
#include <opm/material/fluidmatrixinteractions/EclTwoPhaseMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/PiecewiseLinearTwoPhaseMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/PiecewiseLinearTwoPhaseMaterialParams.hpp>
#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>

#include <opm/simulators/flow/GpuEclMaterialLawManager.hpp>
#include <opm/simulators/flow/GpuEclThermalLawManager.hpp>
#include <opm/simulators/flow/GpuFlowProblem.hpp>

#include <opm/simulators/flow/FlowGasWaterEnergyTypeTag.hpp>
#include <opm/simulators/linalg/gpuistl/GpuFlowGasWaterEnergyTypeTags.hpp>

#include <chrono>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <sstream>


/// Build a CO2STORE WATER+GAS deck with the requested DIMENS. The grid is
/// a uniform unit-cube grid; PERMX/PORO are constant. EQUIL gives the
/// pressure profile so different cells (different depths) get different
/// fluid states, exercising the per-cell update on both CPU and GPU.
static std::string makeDeckString(int nx, int ny, int nz)
{
    const long total = static_cast<long>(nx) * ny * nz;
    const long areal = static_cast<long>(nx) * ny;
    std::ostringstream o;
    o << "RUNSPEC\n"
      << "DIMENS\n" << nx << ' ' << ny << ' ' << nz << " /\n"
      << "EQLDIMS\n/\n"
      << "TABDIMS\n/\n"
      << "WATER\nGAS\nCO2STORE\nTHERMAL\nMETRIC\n"
      << "GRID\n"
      << "GRIDFILE\n0 0 /\n"
      << "DX\n" << total << "*1 /\n"
      << "DY\n" << total << "*1 /\n"
      << "DZ\n" << total << "*1 /\n"
      << "TOPS\n" << areal << "*0 /\n"
      << "PERMX\n" << total << "*1013.25 /\n"
      << "PORO\n" << total << "*0.25 /\n"
      << "COPY\nPERMX PERMY /\nPERMX PERMZ /\n/\n"
      << "THCONR\n" << total << "*1.5 /\n"
      << "PROPS\n"
      << "SGWFN\n"
      << "0.000000E+00 0.000000E+00 1.000000E+00 3.060000E-02\n"
      << "1.000000E+00 1.000000E+00 0.000000E+00 3.060000E-01 /\n"
      << "SPECROCK\n20 2125.0  100 2125.0 /\n"
      << "SOLUTION\n"
      << "RPTRST\n'BASIC=0' /\n"
      << "EQUIL\n0 300 100 0 0 0 1 1 0 /\n"
      << "RTEMPVD\n0 36.12\n1000 70.0 /\n"
      << "SCHEDULE\n"
      << "RPTRST\n'BASIC=0' /\n"
      << "TSTEP\n1 /";
    return o.str();
}


static constexpr const char* deckString1 = "-- =============== RUNSPEC\n"
                                           "RUNSPEC\n"
                                           "DIMENS\n"
                                           "3 3 3 /\n"
                                           "EQLDIMS\n"
                                           "/\n"
                                           "TABDIMS\n"
                                           "/\n"
                                           "WATER\n"
                                           "GAS\n"
                                           "CO2STORE\n"
                                           "THERMAL\n"
                                           "METRIC\n"
                                           "-- =============== GRID\n"
                                           "GRID\n"
                                           "GRIDFILE\n"
                                           "0 0 /\n"
                                           "DX\n"
                                           "27*1 /\n"
                                           "DY\n"
                                           "27*1 /\n"
                                           "DZ\n"
                                           "27*1 /\n"
                                           "TOPS\n"
                                           "9*0 /\n"
                                           "PERMX\n"
                                           "27*1013.25 /\n"
                                           "PORO\n"
                                           "27*0.25 /\n"
                                           "COPY\n"
                                           "PERMX PERMY /\n"
                                           "PERMX PERMZ /\n"
                                           "/\n"
                                           "THCONR\n"
                                           "27*1.5 /\n"
                                           "-- =============== PROPS\n"
                                           "PROPS\n"
                                           "SGWFN\n"
                                           "0.000000E+00 0.000000E+00 1.000000E+00 3.060000E-02\n"
                                           "1.000000E+00 1.000000E+00 0.000000E+00 3.060000E-01 /\n"
                                           "SPECROCK\n"
                                           "20 2125.0  100 2125.0 /\n"
                                           "-- =============== SOLUTION\n"
                                           "SOLUTION\n"
                                           "RPTRST\n"
                                           "'BASIC=0' /\n"
                                           "EQUIL\n"
                                           "0 300 100 0 0 0 1 1 0 /\n"
                                           "RTEMPVD\n"
                                           "0 36.12\n"
                                           "1000 70.0 /\n"
                                           "-- =============== SCHEDULE\n"
                                           "SCHEDULE\n"
                                           "RPTRST\n"
                                           "'BASIC=0' /\n"
                                           "TSTEP\n"
                                           "1 /";

using ScalarToUse = Opm::GetPropType<Opm::Properties::TTag::FlowProblem, Opm::Properties::Scalar>;

using GpuB = Opm::gpuistl::GpuBuffer<double>;
using GpuV = Opm::gpuistl::GpuView<double>;
using FluidSystem = Opm::BlackOilFluidSystem<ScalarToUse>;
using Evaluation = Opm::DenseAd::Evaluation<double, 2>;
using BlackOilFluidSystemView
    = Opm::BlackOilFluidSystemNonStatic<ScalarToUse,
                                        Opm::BlackOilDefaultFluidSystemIndices,
                                        Opm::gpuistl::GpuView>;

template <class Value>
OPM_HOST_DEVICE double asDouble(const Value& value)
{
    if constexpr (requires { value.value(); }) {
        return static_cast<double>(value.value());
    }
    else {
        return static_cast<double>(value);
    }
}

/// Compare two scalar-or-Evaluation quantities. Always compares the value
/// part. When \p checkDerivatives is true and both operands are
/// `DenseAd::Evaluation`s, also checks every partial derivative using
/// \c BOOST_CHECK_CLOSE with the supplied tolerance and increments
/// \p derivativeComparisons by the number of partial derivatives compared.
/// If \p checkDerivatives is true but the operands are not Evaluations,
/// the test fails with \c BOOST_FAIL so a misconfigured derivative test
/// cannot silently pass on plain scalars.
template <class CpuValue, class GpuValue>
void checkValueAndDerivatives(const CpuValue& cpuValue,
                              const GpuValue& gpuValue,
                              double tol,
                              bool checkDerivatives,
                              const char* label,
                              std::size_t& derivativeComparisons)
{
    BOOST_CHECK_CLOSE(asDouble(cpuValue), asDouble(gpuValue), tol);

    if (!checkDerivatives) {
        return;
    }

    if constexpr (requires {
                      CpuValue::numVars;
                      cpuValue.derivative(0);
                      GpuValue::numVars;
                      gpuValue.derivative(0);
                  })
    {
        static_assert(static_cast<int>(CpuValue::numVars)
                          == static_cast<int>(GpuValue::numVars),
                      "CPU and GPU Evaluation types must have the same numVars");
        static_assert(static_cast<int>(CpuValue::numVars) > 0,
                      "Derivative comparison requires at least one partial derivative");
        for (int derivIdx = 0; derivIdx < CpuValue::numVars; ++derivIdx) {
            const double cpuDerivative = static_cast<double>(cpuValue.derivative(derivIdx));
            const double gpuDerivative = static_cast<double>(gpuValue.derivative(derivIdx));
            BOOST_CHECK_MESSAGE(std::abs(cpuDerivative - gpuDerivative)
                                    <= tol * (1.0 + std::abs(cpuDerivative)),
                                label << " derivative[" << derivIdx
                                      << "] cpu=" << cpuDerivative
                                      << " gpu=" << gpuDerivative);
            ++derivativeComparisons;
        }
    }
    else {
        BOOST_FAIL(std::string("checkDerivatives=true but ") + label
                   + " is not a DenseAd::Evaluation");
    }
}

template <class TypeTag>
struct DummyProblem {
 

    OPM_HOST_DEVICE DummyProblem() {
        // empty
    }

    OPM_HOST_DEVICE ~DummyProblem() {
        // empty
    }
     
    using EclMaterialLawManager =
        typename Opm::GetProp<TypeTag, Opm::Properties::MaterialLaw>::EclMaterialLawManager;
    using EclThermalLawManager =
        typename Opm::GetProp<TypeTag, Opm::Properties::SolidEnergyLaw>::EclThermalLawManager;
    struct MaterialLawParams {};

    struct {
        struct {
            OPM_HOST_DEVICE Opm::LinearizationType getLinearizationType() const
            {
                return Opm::LinearizationType();
            }
        } lin_;

        OPM_HOST_DEVICE auto linearizer() const
        {
            return lin_;
        }
    } model_;

    OPM_HOST_DEVICE auto model() const
    {
        return model_;
    }

    OPM_HOST_DEVICE int satnumRegionIndex(std::size_t) const
    {
        return 0;
    }
    OPM_HOST_DEVICE const MaterialLawParams& materialLawParams(std::size_t) const
    {
        return materialLawParams_;
    }
    MaterialLawParams materialLawParams_;
    OPM_HOST_DEVICE double rockCompressibility(std::size_t) const
    {
        return 0.0;
    }
    OPM_HOST_DEVICE double rockReferencePressure(std::size_t) const
    {
        return 0.0;
    }
    OPM_HOST_DEVICE double porosity(std::size_t, unsigned int) const
    {
        return 0.0;
    }
    OPM_HOST_DEVICE double maxOilVaporizationFactor(unsigned int, std::size_t) const
    {
        return 0.0;
    }
    OPM_HOST_DEVICE double maxGasDissolutionFactor(unsigned int, std::size_t) const
    {
        return 0.0;
    }
    OPM_HOST_DEVICE double maxOilSaturation(std::size_t) const
    {
        return 0.0;
    }

    template <class Evaluation>
    OPM_HOST_DEVICE Evaluation rockCompPoroMultiplier(const auto&, std::size_t) const
    {
        return Evaluation(0.0);
    }

    template <class FluidState, class ...Args>
    OPM_HOST_DEVICE void updateRelperms(auto& mobility, auto& dirMob, FluidState& fluidState, unsigned globalSpaceIdx) const
    {
    }

    template <class Evaluation>
    OPM_HOST_DEVICE Evaluation rockCompTransMultiplier(const auto&, std::size_t) const
    {
        return Evaluation(0.0);
    }
};

namespace Opm
{
namespace Properties
{
namespace TTag
{

/// CPU side mirror of the simulation TypeTag, with diffusion/dispersion
/// disabled because the simpler \c BlackOilIntensiveQuantities::update()
/// overload used in this test does not support them. All other property
/// bindings (Indices, MaterialLaw, FluidSystem, etc.) are inherited as-is
/// from \c FlowGasWaterEnergyProblem.
struct FlowGasWaterEnergyProblemTest {
    using InheritsFrom = std::tuple<FlowGasWaterEnergyProblem>;
};

/// CPU TypeTag used by the real-deck test: identical to the production
/// \c FlowGasWaterEnergyProblem (does not force-disable diffusion or
/// dispersion) so that decks enabling those keywords can construct the
/// simulator successfully.
struct FlowGasWaterEnergyProblemTestRealDeck {
    using InheritsFrom = std::tuple<FlowGasWaterEnergyProblem>;
};

} // namespace TTag

template <class TypeTag>
struct EnableDiffusion<TypeTag, TTag::FlowGasWaterEnergyProblemTest>
{ static constexpr bool value = false; };

template <class TypeTag>
struct EnableDispersion<TypeTag, TTag::FlowGasWaterEnergyProblemTest>
{ static constexpr bool value = false; };

} // namespace Properties

} // namespace Opm

using TypeTag = Opm::Properties::TTag::FlowGasWaterEnergyProblemTest;
using TypeNacht = Opm::Properties::TTag::FlowGasWaterEnergyDummyProblemGPU;
// FlowGasWaterEnergyKernelBaseGPU is the GPU-kernel-side TypeTag with GpuView
// FluidSystem; previously named FlowGasWaterEnergyProblemGPU in the dispatcher
// internals.  Use it wherever a GpuView-backed ScalarFluidState is needed.
using TypeTagGPU = Opm::Properties::TTag::FlowGasWaterEnergyKernelBaseGPU;

template<class IndexTraits, class GpuProblem>
__global__ void
testUsingOnGPU(Opm::BlackOilFluidSystemNonStatic<double, IndexTraits, Opm::gpuistl::GpuView> fluidSystem,
               Opm::BlackOilIntensiveQuantities<TypeNacht> intensiveQuantities,
               Opm::BlackOilPrimaryVariables<TypeNacht, Opm::gpuistl::MiniVector> primaryVariables,
               GpuProblem problem)
{
    printf("fluidSystem.phaseIsActive(0): %d\n", fluidSystem.phaseIsActive(0));
    printf("fluidSystem.phaseIsActive(1): %d\n", fluidSystem.phaseIsActive(1));
    printf("fluidSystem.phaseIsActive(2): %d\n", fluidSystem.phaseIsActive(2));
    using ScalarFluidState = typename Opm::BlackOilIntensiveQuantities<TypeTagGPU>::ScalarFluidState;
    ScalarFluidState state(fluidSystem);
    printf("BlackOilState density before update: %f\n", asDouble(state.density(0)));

    ScalarFluidState state2(fluidSystem);
    primaryVariables.assignNaive(state2);
    printf("BlackOilState density before update: %f\n", asDouble(state.density(0)));
    intensiveQuantities.updateSaturations(primaryVariables, 0, Opm::LinearizationType{});
    intensiveQuantities.update(problem, primaryVariables, 0, 0);

    printf("BlackOilState density after update: %f\n", asDouble(state.density(0)));
}

__global__ void dummykernel()
{
    printf("Hello from dummy kernel!\n");
}

template <class GpuProblem, class PrimaryVariables, class IntensiveQuantities>
__global__ void
updateAllCellsKernel(GpuProblem problem,
                     Opm::gpuistl::GpuView<const PrimaryVariables> primaryVariables,
                     Opm::gpuistl::GpuView<IntensiveQuantities> outIntensiveQuantities,
                     std::size_t numCells)
{
    const std::size_t i = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= numCells) {
        return;
    }
    IntensiveQuantities intensiveQuantities = outIntensiveQuantities[i];
    intensiveQuantities.updateSaturations(primaryVariables[i], 0, Opm::LinearizationType{});
    intensiveQuantities.update(problem, primaryVariables[i], static_cast<unsigned>(i), 0);
    intensiveQuantities.updateEnergyQuantities_(problem, static_cast<unsigned>(i), 0u);
    outIntensiveQuantities[i] = intensiveQuantities;
}

/// One-shot global initialization of MPI, communicator, and parameter
/// registration. Safe to call from multiple BOOST test cases; only the
/// first call performs the work.
static void initSimulatorOnce()
{
    static bool done = false;
    if (done) {
        return;
    }
    using namespace Opm;
    int argc1 = boost::unit_test::framework::master_test_suite().argc;
    char** argv1 = boost::unit_test::framework::master_test_suite().argv;
#if HAVE_DUNE_FEM
    Dune::Fem::MPIManager::initialize(argc1, argv1);
#else
    Dune::MPIHelper::instance(argc1, argv1);
#endif
    FlowGenericVanguard::setCommunication(std::make_unique<Opm::Parallel::Communication>());
    Opm::ThreadManager::registerParameters();
    Opm::NewtonMethodParams<double>::registerParameters();
    BlackoilModelParameters<ScalarToUse>::registerParameters();
    AdaptiveTimeStepping<TypeTag>::registerParameters();
    Parameters::Register<Parameters::EnableTerminalOutput>(
        "Dummy added for the well model to compile.");
    registerAllParameters_<TypeTag>(/*finalizeRegistration=*/false);
    // Also register parameters for the RealDeck TypeTag here so that all
    // tests in the same process share a single (still-open) registration
    // pass. Once endRegistration() (called via setupParameters_) closes
    // registration, no further parameters can be added globally.
    using RealDeckTypeTag = Opm::Properties::TTag::FlowGasWaterEnergyProblemTestRealDeck;
    AdaptiveTimeStepping<RealDeckTypeTag>::registerParameters();
    registerAllParameters_<RealDeckTypeTag>(/*finalizeRegistration=*/false);
    Parameters::endRegistration();
    done = true;
}

/// Run the per-cell intensive quantities update on both GPU and CPU for the
/// deck stored at \p deckPath, optionally measuring wall-clock time on both
/// sides.
///
/// Correctness is checked by comparing CPU vs GPU intensive quantities
/// cell-by-cell at indices \p i = 0, \p sampleStride, 2*\p sampleStride, ...
///
/// \param deckPath        path to a CO2STORE WATER+GAS deck written to disk
/// \param expectedNumCells number of grid cells the deck describes
/// \param sampleStride    stride for sampled correctness; 1 means full check
/// \param measureTiming   if true, prints CPU and GPU wall-clock to stdout
static void runIntensiveQuantitiesTestForDeck(const std::string& deckPath,
                                              std::size_t expectedNumCells,
                                              std::size_t sampleStride,
                                              bool measureTiming,
                                              bool checkDerivatives = false)
{
    Opm::resetLocale();
    initSimulatorOnce();

    using namespace Opm;
    using Simulator = Opm::GetPropType<TypeTag, Opm::Properties::Simulator>;

    const auto filenameArg = std::string{"--ecl-deck-file-name="} + deckPath;
    const char* argv2[] = {
        "test_gpuflowproblem",
        filenameArg.c_str(),
        "--check-satfunc-consistency=false",
    };
    Opm::setupParameters_<TypeTag>(/*argc=*/sizeof(argv2) / sizeof(argv2[0]),
                                   argv2,
                                   /*registerParams=*/false,
                                   /*allowUnused=*/true,
                                   /*handleHelp=*/false,
                                   /*myRank=*/0);

    Opm::FlowGenericVanguard::readDeck(deckPath);
    auto sim = std::make_unique<Simulator>();

    auto& cpuProblem = sim->problem();
    auto& dynamicFluidSystem = FluidSystem::getNonStaticInstance();

    auto dynamicGpuFluidSystemBuffer  = ::Opm::gpuistl::copy_to_gpu(dynamicFluidSystem);
    auto dynamicGpuFluidSystemView = ::Opm::gpuistl::make_view(dynamicGpuFluidSystemBuffer);

    // Place the FluidSystemView in unified memory so that the device-side
    // BlackOilFluidState pointer dereference is valid both on host and on
    // device.
    using FluidSystemViewType = std::decay_t<decltype(dynamicGpuFluidSystemView)>;
    FluidSystemViewType* managedFluidSystemView = nullptr;
    OPM_GPU_SAFE_CALL(cudaMallocManaged(&managedFluidSystemView, sizeof(FluidSystemViewType)));
    new (managedFluidSystemView) FluidSystemViewType(dynamicGpuFluidSystemView);

    using CpuMaterialLawManager = typename Opm::GetProp<TypeTag, Opm::Properties::MaterialLaw>::EclMaterialLawManager;
    using Traits = typename CpuMaterialLawManager::MaterialLaw::Traits;
    using TwoPhaseTraits = Opm::TwoPhaseMaterialTraits<ScalarToUse,
                                                       Traits::wettingPhaseIdx,
                                                       Traits::nonWettingPhaseIdx>;
    using GpuPiecewiseLinearParams = Opm::PiecewiseLinearTwoPhaseMaterialParams<
        TwoPhaseTraits, Opm::gpuistl::GpuView<const ScalarToUse>>;
    using GpuPiecewiseLinearLaw = Opm::PiecewiseLinearTwoPhaseMaterial<TwoPhaseTraits, GpuPiecewiseLinearParams>;
    using GpuMaterialLawParams = Opm::EclTwoPhaseMaterialParams<
        Traits, GpuPiecewiseLinearParams, GpuPiecewiseLinearParams, GpuPiecewiseLinearParams,
        Opm::gpuistl::ValueAsPointer>;
    using GpuMaterialLaw = Opm::EclTwoPhaseMaterial<
        Traits, GpuPiecewiseLinearLaw, GpuPiecewiseLinearLaw, GpuPiecewiseLinearLaw,
        GpuMaterialLawParams>;
    using GpuManagerBuffer = Opm::EclMaterialLaw::GpuManager<
        Traits, GpuPiecewiseLinearLaw, GpuPiecewiseLinearLaw,
        Opm::gpuistl::GpuBuffer, GpuMaterialLaw>;
    using GpuThermalManagerBuffer = Opm::EclThermalLaw::GpuManager<
        ScalarToUse, FluidSystemViewType,
        Opm::gpuistl::GpuBuffer, Opm::gpuistl::GpuView>;
    using GpuProblemBuffer = Opm::GpuFlowProblem<ScalarToUse, GpuManagerBuffer, Opm::gpuistl::GpuBuffer, GpuThermalManagerBuffer>;

    GpuProblemBuffer gpuProblemBuffer(cpuProblem);
    auto gpuProblemView = Opm::gpuistl::make_view(gpuProblemBuffer);

    const std::size_t numCells = cpuProblem.model().numGridDof();
    BOOST_REQUIRE_EQUAL(numCells, expectedNumCells);
    BOOST_REQUIRE_GT(numCells, 0u);

    using PrimaryVariablesCpu = Opm::GetPropType<TypeTag, Opm::Properties::PrimaryVariables>;
    using PrimaryVariablesGpu = Opm::BlackOilPrimaryVariables<TypeNacht, Opm::gpuistl::MiniVector>;
    using IntensiveQuantitiesCpu = Opm::BlackOilIntensiveQuantities<TypeTag>;
    using IntensiveQuantitiesGpu = Opm::BlackOilIntensiveQuantities<TypeNacht>;

    PrimaryVariablesCpu primaryVariablesCpu;
    primaryVariablesCpu.setPrimaryVarsMeaningPressure(Opm::BlackOil::PressureMeaning::Pg);
    PrimaryVariablesGpu primaryVariablesGpu(primaryVariablesCpu);
    std::vector<PrimaryVariablesGpu> hostPrimaryVariablesGpu(numCells, primaryVariablesGpu);
    Opm::gpuistl::GpuBuffer<PrimaryVariablesGpu> primaryVariablesBuffer(hostPrimaryVariablesGpu);

    IntensiveQuantitiesCpu cpuIntensiveQuantitiesPrototype;
    IntensiveQuantitiesGpu gpuIntensiveQuantitiesPrototype = cpuIntensiveQuantitiesPrototype.template withOtherFluidSystem<TypeNacht>(*managedFluidSystemView);
    if (measureTiming) {
        std::cout << std::format("[timing] sizeof(IntensiveQuantitiesGpu)={}  sizeof(IntensiveQuantitiesCpu)={}  "
                                 "sizeof(PrimaryVariablesGpu)={}  cells={}  "
                                 "approx GPU mem for IQ buffer={:.2f} MiB\n",
                                 sizeof(IntensiveQuantitiesGpu), sizeof(IntensiveQuantitiesCpu),
                                 sizeof(PrimaryVariablesGpu), numCells,
                                 (sizeof(IntensiveQuantitiesGpu) * numCells) / (1024.0 * 1024.0));
    }
    std::vector<IntensiveQuantitiesGpu> hostIntensiveQuantities(numCells, gpuIntensiveQuantitiesPrototype);
    Opm::gpuistl::GpuBuffer<IntensiveQuantitiesGpu> intensiveQuantitiesBuffer(hostIntensiveQuantities);

    const unsigned blockSize = 64u;
    const unsigned gridSize = static_cast<unsigned>((numCells + blockSize - 1u) / blockSize);

    // Time the GPU kernel using CUDA events (excludes the host-to-device
    // setup and device-to-host copy already done via GpuBuffer ctors above).
    cudaEvent_t eventStart{}, eventStop{};
    OPM_GPU_SAFE_CALL(cudaEventCreate(&eventStart));
    OPM_GPU_SAFE_CALL(cudaEventCreate(&eventStop));
    OPM_GPU_SAFE_CALL(cudaEventRecord(eventStart));
    updateAllCellsKernel<<<gridSize, blockSize>>>(
        gpuProblemView,
        Opm::gpuistl::GpuView<const PrimaryVariablesGpu>(primaryVariablesBuffer.data(), primaryVariablesBuffer.size()),
        Opm::gpuistl::GpuView<IntensiveQuantitiesGpu>(intensiveQuantitiesBuffer.data(), intensiveQuantitiesBuffer.size()),
        numCells);
    OPM_GPU_SAFE_CALL(cudaEventRecord(eventStop));
    OPM_GPU_SAFE_CALL(cudaEventSynchronize(eventStop));
    OPM_GPU_SAFE_CALL(cudaGetLastError());
    float gpuMilliseconds = 0.0f;
    OPM_GPU_SAFE_CALL(cudaEventElapsedTime(&gpuMilliseconds, eventStart, eventStop));
    OPM_GPU_SAFE_CALL(cudaEventDestroy(eventStart));
    OPM_GPU_SAFE_CALL(cudaEventDestroy(eventStop));

    // CPU reference: same per-cell update on the host side using the full
    // CPU problem; time it with std::chrono.
    std::vector<IntensiveQuantitiesCpu> cpuIntensiveQuantities(numCells);
    const auto cpuT0 = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < numCells; ++i) {
        cpuIntensiveQuantities[i].updateSaturations(primaryVariablesCpu, 0, Opm::LinearizationType{});
        cpuIntensiveQuantities[i].update(cpuProblem, primaryVariablesCpu, static_cast<unsigned>(i), 0);
    }
    const auto cpuT1 = std::chrono::steady_clock::now();
    const double cpuMilliseconds =
        std::chrono::duration<double, std::milli>(cpuT1 - cpuT0).count();

    if (measureTiming) {
        std::cout << std::format("[timing] cells={}  CPU={:.3f} ms  GPU(kernel)={:.3f} ms  "
                                 "speedup={:.2f}x\n",
                                 numCells, cpuMilliseconds, gpuMilliseconds,
                                 cpuMilliseconds / gpuMilliseconds);
    }

    // Bring GPU intensive quantities back and compare cells (sampled by sampleStride).
    std::vector<IntensiveQuantitiesGpu> gpuIntensiveQuantitiesHost(numCells, gpuIntensiveQuantitiesPrototype);
    OPM_GPU_SAFE_CALL(cudaMemcpy(gpuIntensiveQuantitiesHost.data(),
                                 intensiveQuantitiesBuffer.data(),
                                 numCells * sizeof(IntensiveQuantitiesGpu),
                                 cudaMemcpyDeviceToHost));
    BOOST_REQUIRE_EQUAL(gpuIntensiveQuantitiesHost.size(), numCells);

    constexpr double tol = 1e-6;
    // Active phase indices for CO2STORE WATER+GAS (oil disabled).
    constexpr unsigned waterPhaseIdx = FluidSystem::waterPhaseIdx;
    constexpr unsigned gasPhaseIdx   = FluidSystem::gasPhaseIdx;
    const unsigned activePhases[] = {waterPhaseIdx, gasPhaseIdx};

    std::size_t checked = 0;
    std::size_t derivativeComparisons = 0;
    for (std::size_t i = 0; i < numCells; i += sampleStride) {
        const auto& cpuIQ = cpuIntensiveQuantities[i];
        const auto& gpuIQ = gpuIntensiveQuantitiesHost[i];
        const auto& cpuFluidState = cpuIQ.fluidState();
        const auto& gpuFluidState = gpuIQ.fluidState();

        // Per-phase fluid state quantities for the active phases. Restricted
        // to fields that are stored directly on the FluidState (no call into
        // the underlying FluidSystem), since the GPU-side FluidState's
        // FluidSystem holds device pointers (GpuView) that are not
        // dereferenceable on the host.
        for (unsigned p : activePhases) {
            checkValueAndDerivatives(cpuFluidState.saturation(p),
                                     gpuFluidState.saturation(p),
                                     tol, checkDerivatives, "saturation",
                                     derivativeComparisons);
            checkValueAndDerivatives(cpuFluidState.pressure(p),
                                     gpuFluidState.pressure(p),
                                     tol, checkDerivatives, "pressure",
                                     derivativeComparisons);
            checkValueAndDerivatives(cpuFluidState.density(p),
                                     gpuFluidState.density(p),
                                     tol, checkDerivatives, "density",
                                     derivativeComparisons);
            checkValueAndDerivatives(cpuFluidState.invB(p),
                                     gpuFluidState.invB(p),
                                     tol, checkDerivatives, "invB",
                                     derivativeComparisons);
            // Mobility now goes through the opm-common EclTwoPhaseMaterial
            // (relativePermeabilities) on both CPU and GPU, so the values
            // must agree to the same tolerance as the rest of the IQ.
            checkValueAndDerivatives(cpuIQ.mobility(p),
                                     gpuIQ.mobility(p),
                                     tol, checkDerivatives, "mobility",
                                     derivativeComparisons);
        }

        // Scalar / per-cell quantities. Both CPU and GPU FluidState configs
        // have enableVapwat=true and enableDisgasInWater=true so Rsw/Rvw are
        // stored on both sides.
        checkValueAndDerivatives(cpuFluidState.Rsw(), gpuFluidState.Rsw(),
                                 tol, checkDerivatives, "Rsw",
                                 derivativeComparisons);
        checkValueAndDerivatives(cpuFluidState.Rvw(), gpuFluidState.Rvw(),
                                 tol, checkDerivatives, "Rvw",
                                 derivativeComparisons);
        BOOST_CHECK_EQUAL(static_cast<unsigned>(cpuFluidState.pvtRegionIndex()),
                          static_cast<unsigned>(gpuFluidState.pvtRegionIndex()));

        checkValueAndDerivatives(cpuIQ.porosity(), gpuIQ.porosity(),
                                 tol, checkDerivatives, "porosity",
                                 derivativeComparisons);
        BOOST_CHECK_CLOSE(static_cast<double>(cpuIQ.referencePorosity()),
                          static_cast<double>(gpuIQ.referencePorosity()), tol);
        checkValueAndDerivatives(cpuIQ.rockCompTransMultiplier(),
                                 gpuIQ.rockCompTransMultiplier(),
                                 tol, checkDerivatives, "rockCompTransMultiplier",
                                 derivativeComparisons);

        ++checked;
    }
    BOOST_TEST_MESSAGE("Per-cell GPU vs CPU IQ comparison: checked "
                       << checked << " / " << numCells << " cells");
    if (checkDerivatives) {
        BOOST_TEST_MESSAGE("Derivative comparisons performed: " << derivativeComparisons);
        // Guard against silent regressions: if checkDerivatives is requested
        // but the helper never actually compared a derivative (e.g. because
        // the IQ types degraded to plain scalars), fail loudly.
        BOOST_REQUIRE_GT(derivativeComparisons, 0u);
    }
}

/// RAII wrapper around a temporary file. The file is deleted when the
/// object goes out of scope.
struct TemporaryFile {
    std::filesystem::path path;
    explicit TemporaryFile(std::string_view filename)
        : path(std::filesystem::temp_directory_path() / filename) {}
    ~TemporaryFile() { std::filesystem::remove(path); }
    TemporaryFile(const TemporaryFile&) = delete;
    TemporaryFile& operator=(const TemporaryFile&) = delete;
};

/// Variant of \c runIntensiveQuantitiesTestForDeck that takes the per-cell
/// primary variables from the simulator's initial (EQUIL-driven) solution
/// vector instead of using a single uniform prototype. This exercises the
/// exact code path used by the dispatcher in production:
/// \c BlackOilPrimaryVariables<TypeNacht>(cpuPrimaryVariables[i]) -> GPU
/// kernel update -> compare CPU vs GPU IQ field-by-field. \p sampleStride
/// controls how often we check (1 means every cell).
static void runIntensiveQuantitiesTestFromSimulatorSolution(const std::string& deckPath,
                                                            std::size_t sampleStride)
{
    Opm::resetLocale();

    using namespace Opm;
    using LocalTypeTag    = Opm::Properties::TTag::FlowGasWaterEnergyProblemTestRealDeck;
    using LocalGpuTypeTag = Opm::Properties::TTag::FlowGasWaterEnergyDummyProblemGPU;
    using Simulator = Opm::GetPropType<LocalTypeTag, Opm::Properties::Simulator>;

    // Reuse the global init: it registers parameters for both TypeTags so
    // this test works whether it runs first or after another test in the
    // same process (parameter registration is closed after the first
    // setupParameters_ call).
    initSimulatorOnce();

    const auto filenameArg = std::string{"--ecl-deck-file-name="} + deckPath;
    const char* argv2[] = {
        "test_gpuflowproblem",
        filenameArg.c_str(),
        "--check-satfunc-consistency=false",
        "--enable-tuning=false",
    };
    Opm::setupParameters_<LocalTypeTag>(/*argc=*/sizeof(argv2) / sizeof(argv2[0]),
                                        argv2,
                                        /*registerParams=*/false,
                                        /*allowUnused=*/true,
                                        /*handleHelp=*/false,
                                        /*myRank=*/0);

    Opm::FlowGenericVanguard::readDeck(deckPath);
    auto sim = std::make_unique<Simulator>();

    // The Simulator constructor only calls finishInit; the solution vector
    // is not yet EQUIL-populated. Trigger initial-solution application so
    // model().solution(0)[i] holds the real per-cell primary variables.
    sim->model().applyInitialSolution();

    auto& cpuProblem = sim->problem();
    auto& dynamicFluidSystem = FluidSystem::getNonStaticInstance();

    BOOST_TEST_MESSAGE(std::format(
        "FluidSystem phase activation: static[O={}, W={}, G={}]  dynamic[O={}, W={}, G={}]",
        FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx),
        FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx),
        FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx),
        dynamicFluidSystem.phaseIsActive(FluidSystem::oilPhaseIdx),
        dynamicFluidSystem.phaseIsActive(FluidSystem::waterPhaseIdx),
        dynamicFluidSystem.phaseIsActive(FluidSystem::gasPhaseIdx)));

    auto dynamicGpuFluidSystemBuffer  = ::Opm::gpuistl::copy_to_gpu(dynamicFluidSystem);
    auto dynamicGpuFluidSystemView = ::Opm::gpuistl::make_view(dynamicGpuFluidSystemBuffer);

    using FluidSystemViewType = std::decay_t<decltype(dynamicGpuFluidSystemView)>;
    FluidSystemViewType* managedFluidSystemView = nullptr;
    OPM_GPU_SAFE_CALL(cudaMallocManaged(&managedFluidSystemView, sizeof(FluidSystemViewType)));
    new (managedFluidSystemView) FluidSystemViewType(dynamicGpuFluidSystemView);

    BOOST_TEST_MESSAGE(std::format(
        "GPU fluid-system view phase activation: O={} W={} G={}",
        managedFluidSystemView->phaseIsActive(FluidSystem::oilPhaseIdx),
        managedFluidSystemView->phaseIsActive(FluidSystem::waterPhaseIdx),
        managedFluidSystemView->phaseIsActive(FluidSystem::gasPhaseIdx)));

    using CpuMaterialLawManager = typename Opm::GetProp<LocalTypeTag, Opm::Properties::MaterialLaw>::EclMaterialLawManager;
    using Traits = typename CpuMaterialLawManager::MaterialLaw::Traits;
    using TwoPhaseTraits = Opm::TwoPhaseMaterialTraits<ScalarToUse,
                                                       Traits::wettingPhaseIdx,
                                                       Traits::nonWettingPhaseIdx>;
    using GpuPiecewiseLinearParams = Opm::PiecewiseLinearTwoPhaseMaterialParams<
        TwoPhaseTraits, Opm::gpuistl::GpuView<const ScalarToUse>>;
    using GpuPiecewiseLinearLaw = Opm::PiecewiseLinearTwoPhaseMaterial<TwoPhaseTraits, GpuPiecewiseLinearParams>;
    using GpuMaterialLawParams = Opm::EclTwoPhaseMaterialParams<
        Traits, GpuPiecewiseLinearParams, GpuPiecewiseLinearParams, GpuPiecewiseLinearParams,
        Opm::gpuistl::ValueAsPointer>;
    using GpuMaterialLaw = Opm::EclTwoPhaseMaterial<
        Traits, GpuPiecewiseLinearLaw, GpuPiecewiseLinearLaw, GpuPiecewiseLinearLaw,
        GpuMaterialLawParams>;
    using GpuManagerBuffer = Opm::EclMaterialLaw::GpuManager<
        Traits, GpuPiecewiseLinearLaw, GpuPiecewiseLinearLaw,
        Opm::gpuistl::GpuBuffer, GpuMaterialLaw>;
    using GpuThermalManagerBuffer = Opm::EclThermalLaw::GpuManager<
        ScalarToUse, FluidSystemViewType,
        Opm::gpuistl::GpuBuffer, Opm::gpuistl::GpuView>;
    using GpuProblemBuffer = Opm::GpuFlowProblem<ScalarToUse, GpuManagerBuffer, Opm::gpuistl::GpuBuffer, GpuThermalManagerBuffer>;

    GpuProblemBuffer gpuProblemBuffer(cpuProblem);
    auto gpuProblemView = Opm::gpuistl::make_view(gpuProblemBuffer);

    const std::size_t numCells = cpuProblem.model().numGridDof();
    BOOST_REQUIRE_GT(numCells, 0u);
    BOOST_TEST_MESSAGE("Real-deck test: numCells=" << numCells);

    using PrimaryVariablesCpu = Opm::GetPropType<LocalTypeTag, Opm::Properties::PrimaryVariables>;
    using PrimaryVariablesGpu = Opm::BlackOilPrimaryVariables<LocalGpuTypeTag, Opm::gpuistl::MiniVector>;
    using IntensiveQuantitiesCpu = Opm::BlackOilIntensiveQuantities<LocalTypeTag>;
    using IntensiveQuantitiesGpu = Opm::BlackOilIntensiveQuantities<LocalGpuTypeTag>;

    // Pull the EQUIL-initialized per-cell primary variables straight from
    // the simulator's solution vector and convert them to the GPU
    // PrimaryVariables type cell-by-cell.
    const auto& cpuSolution = cpuProblem.model().solution(/*timeIdx=*/0);
    BOOST_REQUIRE_EQUAL(cpuSolution.size(), numCells);
    std::vector<PrimaryVariablesCpu> cpuPrimaryVariables(numCells);
    std::vector<PrimaryVariablesGpu> hostPrimaryVariablesGpu;
    hostPrimaryVariablesGpu.reserve(numCells);
    for (std::size_t i = 0; i < numCells; ++i) {
        cpuPrimaryVariables[i] = cpuSolution[i];
        hostPrimaryVariablesGpu.emplace_back(cpuPrimaryVariables[i]);
    }
    Opm::gpuistl::GpuBuffer<PrimaryVariablesGpu> primaryVariablesBuffer(hostPrimaryVariablesGpu);

    // Diagnostic: distribution of pressure-meaning values across cells.
    {
        std::array<std::size_t, 3> meaningCounts{0, 0, 0};
        for (std::size_t i = 0; i < numCells; ++i) {
            const auto m = cpuPrimaryVariables[i].primaryVarsMeaningPressure();
            switch (m) {
                case Opm::BlackOil::PressureMeaning::Pg: ++meaningCounts[0]; break;
                case Opm::BlackOil::PressureMeaning::Po: ++meaningCounts[1]; break;
                case Opm::BlackOil::PressureMeaning::Pw: ++meaningCounts[2]; break;
                default: break;
            }
        }
        BOOST_TEST_MESSAGE(std::format(
            "PressureMeaning distribution (cells={}): Pg={} Po={} Pw={}",
            numCells, meaningCounts[0], meaningCounts[1], meaningCounts[2]));
    }
    if (numCells > 0) {
        const auto& pv0 = cpuPrimaryVariables[0];
        BOOST_TEST_MESSAGE(std::format(
            "Cell 0 PV: meaning(P={},W={},G={}) pvtRegion={}",
            static_cast<int>(pv0.primaryVarsMeaningPressure()),
            static_cast<int>(pv0.primaryVarsMeaningWater()),
            static_cast<int>(pv0.primaryVarsMeaningGas()),
            pv0.pvtRegionIndex()));
        for (unsigned k = 0; k < pv0.size(); ++k) {
            BOOST_TEST_MESSAGE(std::format("  pv[{}] = {}", k, double(pv0[k])));
        }
    }


    IntensiveQuantitiesCpu cpuIntensiveQuantitiesPrototype;
    IntensiveQuantitiesGpu gpuIntensiveQuantitiesPrototype =
        cpuIntensiveQuantitiesPrototype.template withOtherFluidSystem<TypeNacht>(*managedFluidSystemView);
    std::vector<IntensiveQuantitiesGpu> hostIntensiveQuantities(numCells, gpuIntensiveQuantitiesPrototype);
    Opm::gpuistl::GpuBuffer<IntensiveQuantitiesGpu> intensiveQuantitiesBuffer(hostIntensiveQuantities);

    const unsigned blockSize = 64u;
    const unsigned gridSize = static_cast<unsigned>((numCells + blockSize - 1u) / blockSize);

    updateAllCellsKernel<<<gridSize, blockSize>>>(
        gpuProblemView,
        Opm::gpuistl::GpuView<const PrimaryVariablesGpu>(
            primaryVariablesBuffer.data(), primaryVariablesBuffer.size()),
        Opm::gpuistl::GpuView<IntensiveQuantitiesGpu>(
            intensiveQuantitiesBuffer.data(), intensiveQuantitiesBuffer.size()),
        numCells);
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
    OPM_GPU_SAFE_CALL(cudaGetLastError());

    // CPU reference: use the IntensiveQuantities the simulator itself
    // computed during invalidateAndUpdateIntensiveQuantities(0). This is the
    // exact CPU value that the production dispatcher path would overlay
    // GPU-computed BlackOil fields onto, so it is a meaningful reference for
    // the GPU comparison without requiring the (diffusion-incompatible)
    // simpler IQ::update overload to compile here.
    sim->model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);
    std::vector<IntensiveQuantitiesCpu> cpuIntensiveQuantities(numCells);
    for (std::size_t i = 0; i < numCells; ++i) {
        const auto* cached = sim->model().cachedIntensiveQuantities(static_cast<unsigned>(i), 0);
        BOOST_REQUIRE(cached != nullptr);
        cpuIntensiveQuantities[i] = *cached;
    }

    std::vector<IntensiveQuantitiesGpu> gpuIntensiveQuantitiesHost(numCells, gpuIntensiveQuantitiesPrototype);
    OPM_GPU_SAFE_CALL(cudaMemcpy(gpuIntensiveQuantitiesHost.data(),
                                 intensiveQuantitiesBuffer.data(),
                                 numCells * sizeof(IntensiveQuantitiesGpu),
                                 cudaMemcpyDeviceToHost));

    constexpr double tol = 1e-6;
    constexpr unsigned waterPhaseIdx = FluidSystem::waterPhaseIdx;
    constexpr unsigned gasPhaseIdx   = FluidSystem::gasPhaseIdx;
    const unsigned activePhases[] = {waterPhaseIdx, gasPhaseIdx};

    std::size_t checked = 0;
    std::size_t derivativeComparisons = 0;
    std::size_t firstBadCell = static_cast<std::size_t>(-1);
    for (std::size_t i = 0; i < numCells; i += sampleStride) {
        const auto& cpuIQ = cpuIntensiveQuantities[i];
        const auto& gpuIQ = gpuIntensiveQuantitiesHost[i];
        const auto& cpuFs = cpuIQ.fluidState();
        const auto& gpuFs = gpuIQ.fluidState();

        // Probe non-finite values explicitly; report the first offender so
        // we can debug the production simulator failure mode.
        for (unsigned p : activePhases) {
            const double cpuRho = asDouble(cpuFs.density(p));
            const double gpuRho = asDouble(gpuFs.density(p));
            const double cpuP   = asDouble(cpuFs.pressure(p));
            const double gpuP   = asDouble(gpuFs.pressure(p));
            const double cpuS   = asDouble(cpuFs.saturation(p));
            const double gpuS   = asDouble(gpuFs.saturation(p));
            if (!std::isfinite(gpuRho) || !std::isfinite(gpuP) || !std::isfinite(gpuS)) {
                if (firstBadCell == static_cast<std::size_t>(-1)) {
                    firstBadCell = i;
                    BOOST_TEST_MESSAGE(std::format(
                        "First non-finite GPU IQ at cell {} (phase {}): "
                        "S(cpu={}, gpu={}) P(cpu={}, gpu={}) rho(cpu={}, gpu={})",
                        i, p, cpuS, gpuS, cpuP, gpuP, cpuRho, gpuRho));
                }
            }
        }

        for (unsigned p : activePhases) {
            checkValueAndDerivatives(cpuFs.saturation(p), gpuFs.saturation(p),
                                     tol, /*checkDerivatives=*/false, "saturation",
                                     derivativeComparisons);
            checkValueAndDerivatives(cpuFs.pressure(p), gpuFs.pressure(p),
                                     tol, /*checkDerivatives=*/false, "pressure",
                                     derivativeComparisons);
            checkValueAndDerivatives(cpuFs.density(p), gpuFs.density(p),
                                     tol, /*checkDerivatives=*/false, "density",
                                     derivativeComparisons);
            checkValueAndDerivatives(cpuFs.invB(p), gpuFs.invB(p),
                                     tol, /*checkDerivatives=*/false, "invB",
                                     derivativeComparisons);
            // mobility is intentionally NOT compared: the GPU dispatcher
            // does not overlay it (left to the CPU-computed value in
            // production).
        }
        checkValueAndDerivatives(cpuFs.Rsw(), gpuFs.Rsw(),
                                 tol, /*checkDerivatives=*/false, "Rsw",
                                 derivativeComparisons);
        checkValueAndDerivatives(cpuFs.Rvw(), gpuFs.Rvw(),
                                 tol, /*checkDerivatives=*/false, "Rvw",
                                 derivativeComparisons);
        BOOST_CHECK_EQUAL(static_cast<unsigned>(cpuFs.pvtRegionIndex()),
                          static_cast<unsigned>(gpuFs.pvtRegionIndex()));
        checkValueAndDerivatives(cpuIQ.porosity(), gpuIQ.porosity(),
                                 tol, /*checkDerivatives=*/false, "porosity",
                                 derivativeComparisons);
        BOOST_CHECK_CLOSE(static_cast<double>(cpuIQ.referencePorosity()),
                          static_cast<double>(gpuIQ.referencePorosity()), tol);

        // Thermal IQ fields: temperature in the fluid state, and the three
        // BlackOilEnergyIntensiveQuantities scalars (rockInternalEnergy,
        // totalThermalConductivity, rockFraction). For FullyImplicitThermal
        // the GPU dispatcher writes all of these via overlayBlackOilFieldsFrom,
        // so they should match the CPU value field-by-field.
        checkValueAndDerivatives(cpuFs.temperature(/*phaseIdx=*/0),
                                 gpuFs.temperature(/*phaseIdx=*/0),
                                 tol, /*checkDerivatives=*/false, "temperature",
                                 derivativeComparisons);
        checkValueAndDerivatives(cpuIQ.rockInternalEnergy(),
                                 gpuIQ.rockInternalEnergy(),
                                 tol, /*checkDerivatives=*/false, "rockInternalEnergy",
                                 derivativeComparisons);
        checkValueAndDerivatives(cpuIQ.totalThermalConductivity(),
                                 gpuIQ.totalThermalConductivity(),
                                 tol, /*checkDerivatives=*/false, "totalThermalConductivity",
                                 derivativeComparisons);
        BOOST_CHECK_CLOSE(static_cast<double>(cpuIQ.rockFraction()),
                          static_cast<double>(gpuIQ.rockFraction()), tol);
        for (unsigned p : activePhases) {
            checkValueAndDerivatives(cpuFs.enthalpy(p), gpuFs.enthalpy(p),
                                     tol, /*checkDerivatives=*/false, "enthalpy",
                                     derivativeComparisons);
        }
        ++checked;
    }
    BOOST_TEST_MESSAGE("Real-deck per-cell GPU vs CPU IQ comparison: checked "
                       << checked << " / " << numCells << " cells");
}

BOOST_AUTO_TEST_CASE(TestRealDeckGpuVsCpuFromSimulatorSolution)
{
    // Drive the per-cell IQ update with the EQUIL-initialized primary
    // variables of the production CO2STORE deck the user reports diverges
    // when the dispatcher is enabled. This reproduces the exact code path
    // taken by the dispatcher in production (per-cell PrimaryVariables
    // upload + GPU update + readback) and compares every IQ field against
    // a CPU reference, so we can identify which field becomes non-finite.
    //
    // NOTE: this test must run BEFORE the synthetic-deck tests because the
    // BlackOilFluidSystem singleton is global state that gets re-set by
    // every readDeck call; running this after the synthetic tests in the
    // same process leaves the static FluidSystem in a state that is
    // inconsistent with the per-cell GPU FluidSystem copy taken here.
    const std::string deckPath = "/workspaces/opm/thecaseiwant/deck/THECASEIWANT.DATA";
    if (!std::filesystem::exists(deckPath)) {
        BOOST_TEST_MESSAGE("Skipping: deck not found at " << deckPath);
        return;
    }
    runIntensiveQuantitiesTestFromSimulatorSolution(deckPath, /*sampleStride=*/1u);
}

BOOST_AUTO_TEST_CASE(TestInstantiateGpuFlowProblem)
{
    // Persist the in-source two-phase WATER+GAS+CO2STORE deck to a file so
    // the standard Vanguard read path can ingest it.
    const TemporaryFile tempFile("test_blackoilintensivequantities_gpu.DATA");
    {
        std::ofstream f(tempFile.path);
        f << deckString1;
    }
    runIntensiveQuantitiesTestForDeck(tempFile.path.string(),
                                      /*expectedNumCells=*/27u,
                                      /*sampleStride=*/1u,
                                      /*measureTiming=*/false);
}

BOOST_AUTO_TEST_CASE(TestPerformanceGpuVsCpuLargeGrid)
{
    // 100^3 = 1'000'000 cells. Build a CO2STORE WATER+GAS deck with a uniform
    // unit-cube grid and EQUIL initialization, then time the per-cell IQ
    // update on both GPU and CPU. Correctness is checked at every 1000th
    // cell (1000 sample points) so the test stays reasonably fast.
    // 64^3 = 262'144 cells. Build a CO2STORE WATER+GAS deck with a uniform
    // unit-cube grid and EQUIL initialization, then time the per-cell IQ
    // update on both GPU and CPU. Correctness is checked at every 256th cell
    // (~1024 sample points) so the test stays reasonably fast.
    //
    // NOTE: the upper grid size here is bounded by the per-cell device
    // allocation pattern in GpuEclMaterialLawManager (one cudaMalloc per
    // piecewise-linear sample array per cell, i.e. ~12 cudaMalloc per cell).
    // With ROCm's typical 4 KiB allocation granularity, ~16 GiB of VRAM
    // limits us to a few hundred thousand cells; a 1M-cell deck would
    // require consolidating per-cell sample buffers (left as a future
    // optimisation, since for a homogeneous deck all cells share the same
    // SGWFN table).
    constexpr int dim = 64;
    constexpr std::size_t expected =
        static_cast<std::size_t>(dim) * dim * dim;
    const TemporaryFile tempFile("test_blackoilintensivequantities_gpu_1M.DATA");
    {
        std::ofstream f(tempFile.path);
        f << makeDeckString(dim, dim, dim);
    }
    runIntensiveQuantitiesTestForDeck(tempFile.path.string(),
                                      /*expectedNumCells=*/expected,
                                      /*sampleStride=*/256u,
                                      /*measureTiming=*/true);
}

BOOST_AUTO_TEST_CASE(TestGpuVsCpuIntensiveQuantitiesEvaluationDerivatives)
{
    // Re-runs the small (3x3x3) WATER+GAS+CO2STORE deck used by
    // TestInstantiateGpuFlowProblem, but in addition to the value of every
    // intensive-quantity field this test also compares every partial
    // derivative of the underlying DenseAd::Evaluation between CPU and GPU.
    // This guarantees that the GPU material-law / fluid-system evaluation
    // chain is differentiated identically to the CPU one.
    const TemporaryFile tempFile("test_blackoilintensivequantities_gpu_deriv.DATA");
    {
        std::ofstream f(tempFile.path);
        f << deckString1;
    }
    runIntensiveQuantitiesTestForDeck(tempFile.path.string(),
                                      /*expectedNumCells=*/27u,
                                      /*sampleStride=*/1u,
                                      /*measureTiming=*/false,
                                      /*checkDerivatives=*/true);
}

