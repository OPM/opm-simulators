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
#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>


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
                                           "-- =============== PROPS\n"
                                           "PROPS\n"
                                           "SGWFN\n"
                                           "0.000000E+00 0.000000E+00 1.000000E+00 3.060000E-02\n"
                                           "1.000000E+00 1.000000E+00 0.000000E+00 3.060000E-01 /\n"
                                           "-- =============== SOLUTION\n"
                                           "SOLUTION\n"
                                           "RPTRST\n"
                                           "'BASIC=0' /\n"
                                           "EQUIL\n"
                                           "0 300 100 0 0 0 1 1 0 /\n"
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
        struct FlowSimpleProblem {
            using InheritsFrom = std::tuple<FlowProblem>;
        };

        struct FlowSimpleProblemGPU {
            using InheritsFrom = std::tuple<FlowSimpleProblem>;
        };

        struct FlowSimpleDummyProblemGPU {
            using InheritsFrom = std::tuple<FlowSimpleProblemGPU>;
        };

    } // namespace TTag

    // Indices for two-phase gas-water.
    template <class TypeTag>
    struct Indices<TypeTag, TTag::FlowSimpleProblem> {
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
                                             0 /*numBioCompV*/>;
    };

    // SPE11C requires thermal/energy
    template <class TypeTag>
    struct EnableEnergy<TypeTag, TTag::FlowSimpleProblem> {
        static constexpr bool value = true;
    };

    // SPE11C requires dispersion
    template <class TypeTag>
    struct EnableDispersion<TypeTag, TTag::FlowSimpleProblem> {
        static constexpr bool value = false;
    };

    // Use the simple material law.
    template <class TypeTag>
    struct MaterialLaw<TypeTag, TTag::FlowSimpleProblem> {
    private:
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

        using Traits = ThreePhaseMaterialTraits<Scalar,
                                                /*wettingPhaseIdx=*/FluidSystem::waterPhaseIdx,
                                                /*nonWettingPhaseIdx=*/FluidSystem::oilPhaseIdx,
                                                /*gasPhaseIdx=*/FluidSystem::gasPhaseIdx,
                                                /*hysteresis=*/false,
                                                /*enableEndPointScaling=*/true>;

    public:
        using EclMaterialLawManager = ::Opm::EclMaterialLaw::Manager<Traits>;
        using type = typename EclMaterialLawManager::MaterialLaw;
    };

    // Use the TPFA linearizer.
    template <class TypeTag>
    struct Linearizer<TypeTag, TTag::FlowSimpleProblem> {
        using type = TpfaLinearizer<TypeTag>;
    };

    template <class TypeTag>
    struct LocalResidual<TypeTag, TTag::FlowSimpleProblem> {
        using type = BlackOilLocalResidualTPFA<TypeTag>;
    };

    // Diffusion.
    template <class TypeTag>
    struct EnableDiffusion<TypeTag, TTag::FlowSimpleProblem> {
        static constexpr bool value = false;
    };

    template <class TypeTag>
    struct EnableDisgasInWater<TypeTag, TTag::FlowSimpleProblem> {
        static constexpr bool value = true;
    };

    template <class TypeTag>
    struct EnableVapwat<TypeTag, TTag::FlowSimpleProblem> {
        static constexpr bool value = true;
    };

    template <class TypeTag>
    struct PrimaryVariables<TypeTag, TTag::FlowSimpleProblem> {
        using type = BlackOilPrimaryVariables<TypeTag>;
    };

    template <class TypeTag>
    struct PrimaryVariables<TypeTag, TTag::FlowSimpleProblemGPU> {
        using type = BlackOilPrimaryVariables<TypeTag, Opm::gpuistl::MiniVector>;
    };

    template <class TypeTag>
    struct IntensiveQuantities<TypeTag, TTag::FlowSimpleProblem> {
        using type = BlackOilIntensiveQuantities<TypeTag>;
    };

    template <class TypeTag>
    struct Problem<TypeTag, TTag::FlowSimpleDummyProblemGPU> {
        using type = DummyProblem<TypeTag>;
    };

    template <class TypeTag>
    struct FluidSystem<TypeTag, TTag::FlowSimpleDummyProblemGPU> {
        using type = BlackOilFluidSystemView;
    };

    template <class TypeTag>
    struct FluidSystem<TypeTag, TTag::FlowSimpleProblemGPU> {
        using type = BlackOilFluidSystemView;
    };

    template<class TypeTag>
    struct ElementContext<TypeTag, TTag::FlowSimpleDummyProblemGPU>
    { using type = FvBaseElementContextGpu<TypeTag>; };

}; // namespace Properties

} // namespace Opm

using TypeTag = Opm::Properties::TTag::FlowSimpleProblem;
using TypeNacht = Opm::Properties::TTag::FlowSimpleDummyProblemGPU;
using TypeTagGPU = Opm::Properties::TTag::FlowSimpleProblemGPU;

template<class IndexTraits>
__global__ void
testUsingOnGPU(Opm::BlackOilFluidSystemNonStatic<double, IndexTraits, Opm::gpuistl::GpuView> fs,  Opm::BlackOilIntensiveQuantities<TypeNacht> intensiveQuantities, Opm::BlackOilPrimaryVariables<TypeNacht, Opm::gpuistl::MiniVector> primaryVariables)
{
    printf("fs.phaseIsActive(0): %d\n", fs.phaseIsActive(0));
    printf("fs.phaseIsActive(1): %d\n", fs.phaseIsActive(1));
    printf("fs.phaseIsActive(2): %d\n", fs.phaseIsActive(2));
    using ScalarFluidState = typename Opm::BlackOilIntensiveQuantities<TypeTagGPU>::ScalarFluidState;
    ScalarFluidState state(fs);
    printf("BlackOilState density before update: %f\n", asDouble(state.density(0)));

    ScalarFluidState state2(fs);
    primaryVariables.assignNaive(state2);
    printf("BlackOilState density before update: %f\n", asDouble(state.density(0)));
    printf("BlackOilState density after update: %f\n", asDouble(state.density(0)));
}

__global__ void dummykernel()
{
    printf("Hello from dummy kernel!\n");
}

BOOST_AUTO_TEST_CASE(TestPrimaryVariablesCreationGPU)
{
    Opm::Parser parser;

    auto deck = parser.parseString(deckString1);
    auto python = std::make_shared<Opm::Python>();
    Opm::EclipseState eclState(deck);
    Opm::Schedule schedule(deck, eclState, python);

    FluidSystem::initFromState(eclState, schedule);

    auto& dynamicFluidSystem = FluidSystem::getNonStaticInstance();

    auto dynamicGpuFluidSystemBuffer = ::Opm::gpuistl::copy_to_gpu(dynamicFluidSystem);
    auto dynamicGpuFluidSystemView = ::Opm::gpuistl::make_view(dynamicGpuFluidSystemBuffer);
    std::cout << "phaseIsActive: " << dynamicFluidSystem.phaseIsActive(0) << ", "
              << dynamicFluidSystem.phaseIsActive(1) << ", " << dynamicFluidSystem.phaseIsActive(2)
              << std::endl;
    std::cout << "blackoil is active" << dynamicFluidSystem.phaseIsActive(FluidSystem::oilPhaseIdx)
              << std::endl;
    Opm::BlackOilIntensiveQuantities<TypeTag> intensiveQuantities;

    auto& state = intensiveQuantities.fluidState();


    printf("(CPU) BlackOilState density before update: %f\n", asDouble(state.density(0)));
    intensiveQuantities.updatePhaseDensities();
    printf("(CPU) BlackOilState density after update: %f\n", asDouble(state.density(0)));

    using PrimaryVariables = Opm::GetPropType<TypeTag, Opm::Properties::PrimaryVariables>;
    PrimaryVariables primaryVariablesCPU;
    Opm::BlackOilPrimaryVariables<TypeNacht, Opm::gpuistl::MiniVector> primaryVariablesGPU(primaryVariablesCPU);

    Opm::BlackOilIntensiveQuantities<TypeNacht> intensiveQuantitiesGPU = intensiveQuantities.withOtherFluidSystem<TypeNacht>(dynamicGpuFluidSystemView);
    testUsingOnGPU<<<1, 1>>>(dynamicGpuFluidSystemView, intensiveQuantitiesGPU, primaryVariablesGPU);
    dummykernel<<<1,1>>>();
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
    OPM_GPU_SAFE_CALL(cudaGetLastError());
}


BOOST_AUTO_TEST_CASE(TestInstantiateGpuFlowProblem)
{
    Opm::resetLocale();
    using TypeTag = Opm::Properties::TTag::FlowSimpleProblem;
    int argc1 = boost::unit_test::framework::master_test_suite().argc;
    char** argv1 = boost::unit_test::framework::master_test_suite().argv;
    std::cout << "Initializing MPI..." << std::endl;
#if HAVE_DUNE_FEM
    Dune::Fem::MPIManager::initialize(argc1, argv1);
#else
    Dune::MPIHelper::instance(argc1, argv1);
#endif
    std::cout << "Initialized MPI..." << std::endl;

    using namespace Opm;
    FlowGenericVanguard::setCommunication(std::make_unique<Opm::Parallel::Communication>());
    Opm::ThreadManager::registerParameters();
    Opm::NewtonMethodParams<double>::registerParameters();
    std::cout << "Registered ThreadManager parameters." << std::endl;
    BlackoilModelParameters<ScalarToUse>::registerParameters();
    AdaptiveTimeStepping<TypeTag>::registerParameters();
    Parameters::Register<Parameters::EnableTerminalOutput>(
        "Dummy added for the well model to compile.");

    using Simulator = Opm::GetPropType<TypeTag, Opm::Properties::Simulator>;

    const std::string filename = "very_simple_deck.DATA";
    const auto filenameArg = std::string {"--ecl-deck-file-name="} + filename;

    const char* argv2[] = {
        "test_gpuflowproblem",
        filenameArg.c_str(),
        "--check-satfunc-consistency=false",
    };
    registerAllParameters_<TypeTag>(true);
    Opm::setupParameters_<TypeTag>(/*argc=*/sizeof(argv2) / sizeof(argv2[0]),
                                   argv2,
                                   /*registerParams=*/false,
                                   /*allowUnused=*/true,
                                   /*handleHelp=*/false,
                                   /*myRank=*/0);

    std::cout << "Registered all parameters." << std::endl;
    std::cout << "Reading deck" << std::endl;

    if (!std::filesystem::exists(filename)) {
        throw std::runtime_error(std::format("Missing file {}.", filename));
    }


    Opm::FlowGenericVanguard::readDeck(filename);
    std::cout << "Read deck." << std::endl;
    std::cout << "Initializing simulator..." << std::endl;
    auto sim = std::make_unique<Simulator>();
    std::cout << "Created simulator." << std::endl;

    auto& dynamicFluidSystem = FluidSystem::getNonStaticInstance();
    std::cout << "Initialized dynamic fluid system." << std::endl;

    
    std::cout << "phaseIsActive: " << dynamicFluidSystem.phaseIsActive(0) << ", " << dynamicFluidSystem.phaseIsActive(1)
              << ", " << dynamicFluidSystem.phaseIsActive(2) << std::endl;
    auto dynamicGpuFluidSystemBuffer
        = ::Opm::gpuistl::copy_to_gpu(dynamicFluidSystem);
    std::cout << "Copied dynamic fluid system to GPU." << std::endl;
    auto dynamicGpuFluidSystemView = ::Opm::gpuistl::make_view(
        dynamicGpuFluidSystemBuffer);
    std::cout << "Created GPU view of dynamic fluid system." << std::endl;

}
BOOST_AUTO_TEST_CASE(TestIntensiveQuantitiesCreationGPU)
{
    // This test is currently just to check that the code compiles and can be launched on the GPU
    // without errors. It does not actually check that the intensive quantities are correct. Adding
    // such checks would require implementing a CPU version of the intensive quantities, which is
    // non-trivial and currently not available.
}
