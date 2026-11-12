/*
  Copyright 2026 Equinor ASA

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
/*!
 * \file
 *
 * \brief Compares BlackOilIntensiveQuantities::update() on the GPU with the CPU
 *        for the TypeTag that flow uses for three-phase black-oil decks with live
 *        oil and dry gas (e.g. SPE1) or wet gas (e.g. Norne).
 *
 * The CPU simulator is set up with \c FlowProblemNohyst from an SPE1-like deck
 * whose gas-oil and oil-water contacts both lie inside the grid, and from a
 * variant of it with vaporized oil. The GPU side uses the TypeTags of
 * GpuFlowBlackoilNohystTypeTag.hpp, and every cell is compared, values and
 * derivatives.
 *
 * Also checks that the production dispatcher for this TypeTag rejects a deck
 * with endpoint scaling, which the GPU material law manager does not apply.
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

#define BOOST_TEST_MODULE TestBlackOilIntensiveQuantitiesGpuThreePhase

#include <boost/test/unit_test.hpp>

#include <cuda_runtime.h>

#include <opm/common/utility/gpuDecorators.hpp>

#include <opm/material/common/ResetLocale.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystemNonStatic.hpp>
#include <opm/material/fluidsystems/blackoilpvt/BlackOilPvtApproaches.hpp>

#include <opm/models/blackoil/blackoilintensivequantities.hh>
#include <opm/models/nonlinear/newtonmethodparams.hpp>
#include <opm/models/utils/start.hh>

#include <opm/simulators/flow/BlackoilModelParameters.hpp>
#include <opm/simulators/flow/FlowGenericVanguard.hpp>
#include <opm/simulators/flow/FlowProblemNohystTypeTag.hpp>
#include <opm/simulators/flow/GpuFlowProblem.hpp>
#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/linalg/gpuistl/GpuBlackoilIntensiveQuantitiesDispatcher.hpp>
#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuFlowBlackoilNohystTypeTag.hpp>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>
#include <opm/simulators/linalg/gpuistl/MiniVector.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>
#include <opm/simulators/linalg/gpuistl/gpu_smart_pointer.hpp>

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

namespace Opm {
    // Definition to resolve linking error in test
    std::string moduleVersionName() { return "test-version"; }
}

namespace {

using CpuTypeTag = Opm::Properties::TTag::FlowProblemNohyst;
using Scalar = Opm::GetPropType<CpuTypeTag, Opm::Properties::Scalar>;
using FluidSystem = Opm::GetPropType<CpuTypeTag, Opm::Properties::FluidSystem>;

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
/// `DenseAd::Evaluation`s, also checks every partial derivative with the
/// supplied tolerance and increments \p derivativeComparisons by the number of
/// partial derivatives compared. If \p checkDerivatives is true but the
/// operands are not Evaluations, the test fails with \c BOOST_FAIL so a
/// misconfigured derivative test cannot silently pass on plain scalars.
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

/// Update the intensive quantities of every cell on the GPU, in the same way
/// as the GPU intensive-quantities dispatcher. The TypeTag is not thermal, so
/// there is no energy update.
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
    intensiveQuantities.update(problem, primaryVariables[i], static_cast<unsigned>(i), 0);
    outIntensiveQuantities[i] = intensiveQuantities;
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

// SPE1-like deck on a 3x3x6 grid. The gas-oil contact (2015 m) and the
// water-oil contact (2045 m) both lie inside the grid (2000-2060 m), so EQUIL
// gives cells in the gas cap (Sg meaning, saturated Rs), in the oil zone (Rs
// meaning, undersaturated) and in the water zone. Both SATNUM regions have
// nonzero capillary pressures, and the PVTO table has saturated-only and
// undersaturated rows.
constexpr const char* deckString = R"(
RUNSPEC
DIMENS
 3 3 6 /
TABDIMS
 2 1 /
EQLDIMS
/
OIL
WATER
GAS
DISGAS
METRIC
GRID
DX
 54*100 /
DY
 54*100 /
DZ
 54*10 /
TOPS
 9*2000 /
PORO
 54*0.25 /
PERMX
 54*100 /
COPY
 PERMX PERMY /
 PERMX PERMZ /
/
PROPS
SWOF
-- Sw   krw   krow  Pcow
 0.15  0.0   1.0   0.40
 0.30  0.05  0.60  0.20
 0.50  0.20  0.25  0.10
 0.70  0.45  0.05  0.04
 0.85  0.70  0.0   0.01
 1.00  1.00  0.0   0.0 /
 0.20  0.0   1.0   0.30
 0.40  0.08  0.45  0.15
 0.60  0.25  0.12  0.06
 0.80  0.55  0.0   0.02
 1.00  1.00  0.0   0.0 /
SGOF
-- Sg   krg   krog  Pcog
 0.00  0.0   1.0   0.0
 0.10  0.02  0.65  0.05
 0.30  0.15  0.25  0.12
 0.50  0.40  0.05  0.20
 0.70  0.70  0.0   0.30
 0.85  0.90  0.0   0.35 /
 0.00  0.0   1.0   0.0
 0.15  0.05  0.55  0.08
 0.40  0.30  0.12  0.18
 0.60  0.60  0.0   0.25
 0.80  0.85  0.0   0.30 /
DENSITY
 850.0 1030.0 0.90 /
PVTW
 250.0 1.03 4.5E-5 0.40 0.0 /
ROCK
 200.0 4.0E-5 /
PVDG
  20.0 0.0600 0.012
 100.0 0.0120 0.016
 300.0 0.0045 0.025 /
PVTO
-- RS     P      BO     MU
   20.0   50.0  1.10   1.20 /
   50.0  100.0  1.20   1.00 /
   90.0  150.0  1.30   0.80
         200.0  1.28   0.85
         250.0  1.26   0.90 /
  130.0  200.0  1.40   0.70 /
  170.0  250.0  1.50   0.60
         300.0  1.48   0.62
         350.0  1.46   0.64 /
/
REGIONS
SATNUM
 1 1 2 1 1 2 1 1 2
 1 1 2 1 1 2 1 1 2
 1 1 2 1 1 2 1 1 2
 1 1 2 1 1 2 1 1 2
 1 1 2 1 1 2 1 1 2
 1 1 2 1 1 2 1 1 2 /
SOLUTION
EQUIL
 2030 200 2045 0 2015 0 1 0 0 /
RSVD
 2000 120.0
 2060 120.0 /
SCHEDULE
DRSDT
 0 /
TSTEP
 1 /
)";

void initSimulatorOnce()
{
    static bool done = false;
    if (done) {
        return;
    }
    int argc = boost::unit_test::framework::master_test_suite().argc;
    char** argv = boost::unit_test::framework::master_test_suite().argv;
    Dune::MPIHelper::instance(argc, argv);
    Opm::FlowGenericVanguard::setCommunication(std::make_unique<Opm::Parallel::Communication>());
    Opm::ThreadManager::registerParameters();
    Opm::NewtonMethodParams<double>::registerParameters();
    Opm::BlackoilModelParameters<Scalar>::registerParameters();
    Opm::Parameters::Register<Opm::Parameters::EnableTerminalOutput>(
        "Dummy added for the well model to compile.");
    Opm::AdaptiveTimeStepping<CpuTypeTag>::registerParameters();
    Opm::registerAllParameters_<CpuTypeTag>(/*finalizeRegistration=*/false);
    Opm::Parameters::endRegistration();
    done = true;
}

using Simulator = Opm::GetPropType<CpuTypeTag, Opm::Properties::Simulator>;

/// Write \p deck to a temporary file and set up a simulator with the initial
/// solution applied.
std::unique_ptr<Simulator> makeSimulator(const std::string& deck)
{
    const TemporaryFile deckFile("test_blackoilintensivequantities_gpu_threephase.DATA");
    {
        std::ofstream f(deckFile.path);
        f << deck;
    }
    const std::string deckPath = deckFile.path.string();

    Opm::resetLocale();
    initSimulatorOnce();

    const auto filenameArg = std::string{"--ecl-deck-file-name="} + deckPath;
    const char* argv[] = {
        "test_blackoilintensivequantities_gpu_threephase",
        filenameArg.c_str(),
    };
    Opm::setupParameters_<CpuTypeTag>(/*argc=*/sizeof(argv) / sizeof(argv[0]),
                                      argv,
                                      /*registerParams=*/false,
                                      /*allowUnused=*/true,
                                      /*handleHelp=*/false,
                                      /*myRank=*/0);

    Opm::FlowGenericVanguard::readDeck(deckPath, /*ignoredKeywords=*/"");
    auto sim = std::make_unique<Simulator>();
    sim->model().applyInitialSolution();
    return sim;
}

/// Run the production GPU intensive-quantities dispatcher once on the initial
/// solution of \p sim, as flow does for the first linearization.
void runDispatcher(const Simulator& sim)
{
    using Dispatcher = Opm::gpuistl::GpuBlackoilIntensiveQuantitiesDispatcher<CpuTypeTag>;
    using PrimaryVariables = typename Dispatcher::PrimaryVariables;
    using IntensiveQuantities = typename Dispatcher::IntensiveQuantities;

    const auto& problem = sim.problem();
    const auto& solution = problem.model().solution(/*timeIdx=*/0);
    const std::size_t numCells = problem.model().numGridDof();

    std::vector<IntensiveQuantities> intensiveQuantities(numCells);
    std::vector<const PrimaryVariables*> primaryVariablesPtrs(numCells);
    std::vector<IntensiveQuantities*> intensiveQuantitiesPtrs(numCells);
    for (std::size_t i = 0; i < numCells; ++i) {
        primaryVariablesPtrs[i] = &solution[i];
        intensiveQuantitiesPtrs[i] = &intensiveQuantities[i];
    }

    Dispatcher dispatcher;
    dispatcher.update(problem,
                      primaryVariablesPtrs.data(),
                      intensiveQuantitiesPtrs.data(),
                      numCells,
                      /*timeIdx=*/0u);
}

/// Replace the text from \p begin up to, but not including, \p end.
void replaceSection(std::string& deck, std::string_view begin, std::string_view end,
                    std::string_view text)
{
    const auto first = deck.find(begin);
    BOOST_REQUIRE(first != std::string::npos);
    const auto last = deck.find(end, first);
    BOOST_REQUIRE(last != std::string::npos);
    deck.replace(first, last - first, text);
}

void replaceOnce(std::string& deck, std::string_view from, std::string_view to)
{
    const auto pos = deck.find(from);
    BOOST_REQUIRE(pos != std::string::npos);
    BOOST_REQUIRE(deck.find(from, pos + 1) == std::string::npos);
    deck.replace(pos, from.size(), to);
}

/// Add VAPOIL, with RVVD for the initial Rv (EQUIL item 8). Without RVVD,
/// EQUIL requires the datum depth at the gas-oil contact.
void addVaporizedOil(std::string& deck)
{
    replaceOnce(deck, "DISGAS\n", "DISGAS\nVAPOIL\n");
    replaceOnce(deck, " 2030 200 2045 0 2015 0 1 0 0 /", " 2030 200 2045 0 2015 0 1 1 0 /");
    replaceOnce(deck, " 2060 120.0 /\n", " 2060 120.0 /\nRVVD\n 2000 0.0003\n 2060 0.0003 /\n");
}

/// The test deck with wet gas: PVDG replaced by a PVTG table with the same
/// saturated Bg, which covers the pressures of the deck. With \p vaporizedOil,
/// VAPOIL is added, with an Rv below the saturated value (5e-4 at 200 bar),
/// and the gas-oil contact moves 3 m down. The top layer then holds no oil
/// (Rv meaning), and the cell centres of the second layer lie in the
/// capillary transition zone (oil and gas, Sg meaning).
std::string deckWithWetGas(bool vaporizedOil)
{
    std::string deck = deckString;
    replaceSection(deck, "PVDG\n", "PVTO\n",
                   "PVTG\n"
                   "-- PG     RV        BG       MUG\n"
                   "   20.0   0.00004   0.0600   0.0118\n"
                   "          0.00000   0.0598   0.0115 /\n"
                   "  100.0   0.00025   0.0120   0.0160\n"
                   "          0.00000   0.0118   0.0155 /\n"
                   "  200.0   0.00050   0.0060   0.0210\n"
                   "          0.00025   0.0059   0.0205\n"
                   "          0.00000   0.0058   0.0200 /\n"
                   "  300.0   0.00070   0.0045   0.0250\n"
                   "          0.00000   0.0044   0.0240 /\n"
                   "/\n");
    if (vaporizedOil) {
        addVaporizedOil(deck);
        replaceOnce(deck, " 2030 200 2045 0 2015 0 1 1 0 /", " 2030 200 2045 0 2018 0 1 1 0 /");
    }
    return deck;
}

/// The test deck with endpoint scaling: ENDSCALE in RUNSPEC and a SWL array,
/// but no ENKRVD. SWL equals the connate water saturation of each cell's
/// SATNUM table (0.15 and 0.20), so the scaled end-points equal the unscaled
/// ones and pass the saturation function consistency checks.
std::string deckWithEndpointScaling()
{
    std::string deck = deckString;
    const auto insertAfter = [&deck](std::string_view anchor, std::string_view text) {
        const auto pos = deck.find(anchor);
        BOOST_REQUIRE(pos != std::string::npos);
        deck.insert(pos + anchor.size(), text);
    };
    insertAfter("METRIC\n", "ENDSCALE\n/\n");
    // Same layout as SATNUM: 1 1 2 repeated.
    insertAfter("PROPS\n",
                "SWL\n"
                " 2*0.15 0.20 2*0.15 0.20 2*0.15 0.20\n"
                " 2*0.15 0.20 2*0.15 0.20 2*0.15 0.20\n"
                " 2*0.15 0.20 2*0.15 0.20 2*0.15 0.20\n"
                " 2*0.15 0.20 2*0.15 0.20 2*0.15 0.20\n"
                " 2*0.15 0.20 2*0.15 0.20 2*0.15 0.20\n"
                " 2*0.15 0.20 2*0.15 0.20 2*0.15 0.20 /\n");
    return deck;
}

} // namespace

namespace {

/// Update the intensive quantities of every cell of \p deck on the GPU with
/// \p GpuTypeTag and on the CPU, and compare them cell by cell.
template <class GpuTypeTag, class PvtApproaches>
void compareGpuWithCpu(const std::string& deck)
{
    const auto sim = makeSimulator(deck);
    const bool vaporizedOil = FluidSystem::enableVaporizedOil();
    BOOST_REQUIRE_EQUAL(vaporizedOil, PvtApproaches::gas == Opm::GasPvtApproach::WetGas);
    const auto& cpuProblem = sim->problem();

    // GPU fluid system. The fluid state holds a pointer to it and is read back
    // on the host, so the view lives in managed memory.
    auto gpuFluidSystemBuffer = Opm::gpuistl::copy_to_gpu<PvtApproaches>(
        FluidSystem::getNonStaticInstance());
    auto gpuFluidSystemView = Opm::gpuistl::make_view(gpuFluidSystemBuffer);
    using GpuFluidSystemView = std::decay_t<decltype(gpuFluidSystemView)>;
    static_assert(std::is_same_v<GpuFluidSystemView,
                                 Opm::GetPropType<GpuTypeTag, Opm::Properties::FluidSystem>>);
    auto managedFluidSystemView
        = Opm::gpuistl::make_gpu_managed_unique_ptr<GpuFluidSystemView>(gpuFluidSystemView);

    // GPU problem, with the GpuBuffer-storage types derived from the GPU TypeTag.
    using GpuProblemView = Opm::GetPropType<GpuTypeTag, Opm::Properties::Problem>;
    using GpuViewMaterialLawManager = typename GpuProblemView::EclMaterialLawManager;
    using GpuBufferMaterialLawManager = Opm::EclMaterialLaw::GpuManager<
        typename GpuViewMaterialLawManager::Traits,
        typename GpuViewMaterialLawManager::GasOilLaw,
        typename GpuViewMaterialLawManager::OilWaterLaw,
        Opm::gpuistl::GpuBuffer,
        typename GpuViewMaterialLawManager::MaterialLaw>;
    using GpuProblemBuffer = Opm::GpuFlowProblem<Scalar, GpuBufferMaterialLawManager, Opm::gpuistl::GpuBuffer>;
    GpuProblemBuffer gpuProblemBuffer(cpuProblem);
    auto gpuProblemView = Opm::gpuistl::make_view(gpuProblemBuffer);
    static_assert(std::is_same_v<decltype(gpuProblemView), GpuProblemView>);

    const std::size_t numCells = cpuProblem.model().numGridDof();
    BOOST_REQUIRE_EQUAL(numCells, 54u);

    // Primary variables
    using PrimaryVariablesCpu = Opm::GetPropType<CpuTypeTag, Opm::Properties::PrimaryVariables>;
    using PrimaryVariablesGpu = Opm::GetPropType<GpuTypeTag, Opm::Properties::PrimaryVariables>;
    const auto& cpuSolution = cpuProblem.model().solution(/*timeIdx=*/0);
    BOOST_REQUIRE_EQUAL(cpuSolution.size(), numCells);

    // Check that the deck gives both the saturated and the undersaturated
    // branch of the oil, and with vaporized oil also gas-only cells.
    using GasMeaning = typename PrimaryVariablesCpu::GasMeaning;
    using Indices = Opm::GetPropType<CpuTypeTag, Opm::Properties::Indices>;
    std::size_t numSgCells = 0;
    std::size_t numRsCells = 0;
    std::size_t numRvCells = 0;
    std::vector<PrimaryVariablesGpu> hostPrimaryVariablesGpu;
    hostPrimaryVariablesGpu.reserve(numCells);
    for (std::size_t i = 0; i < numCells; ++i) {
        const PrimaryVariablesCpu& priVars = cpuSolution[i];
        // The water zone also gets the Sg meaning, but with Sg = 0.
        if (priVars.primaryVarsMeaningGas() == GasMeaning::Sg
            && priVars[Indices::compositionSwitchIdx] > 0.0) {
            ++numSgCells;
        }
        else if (priVars.primaryVarsMeaningGas() == GasMeaning::Rs) {
            ++numRsCells;
        }
        else if (priVars.primaryVarsMeaningGas() == GasMeaning::Rv) {
            ++numRvCells;
        }
        hostPrimaryVariablesGpu.emplace_back(priVars);
    }
    BOOST_TEST_MESSAGE("Gas meaning: Sg > 0 in " << numSgCells << " cells, Rs in "
                       << numRsCells << " cells, Rv in " << numRvCells << " cells");
    BOOST_REQUIRE_GT(numSgCells, 0u);
    BOOST_REQUIRE_GT(numRsCells, 0u);
    if (vaporizedOil) {
        BOOST_REQUIRE_GT(numRvCells, 0u);
    }
    else {
        BOOST_REQUIRE_EQUAL(numRvCells, 0u);
    }
    Opm::gpuistl::GpuBuffer<PrimaryVariablesGpu> primaryVariablesBuffer(hostPrimaryVariablesGpu);

    // Intensive quantities
    using IntensiveQuantitiesCpu = Opm::BlackOilIntensiveQuantities<CpuTypeTag>;
    using IntensiveQuantitiesGpu = Opm::BlackOilIntensiveQuantities<GpuTypeTag>;
    const IntensiveQuantitiesCpu cpuIntensiveQuantitiesPrototype;
    const IntensiveQuantitiesGpu gpuIntensiveQuantitiesPrototype
        = cpuIntensiveQuantitiesPrototype.template withOtherFluidSystem<GpuTypeTag>(*managedFluidSystemView);
    Opm::gpuistl::GpuBuffer<IntensiveQuantitiesGpu> intensiveQuantitiesBuffer(
        std::vector<IntensiveQuantitiesGpu>(numCells, gpuIntensiveQuantitiesPrototype));

    const unsigned blockSize = 64u;
    const unsigned gridSize = static_cast<unsigned>((numCells + blockSize - 1u) / blockSize);
    updateAllCellsKernel<<<gridSize, blockSize>>>(
        gpuProblemView,
        Opm::gpuistl::GpuView<const PrimaryVariablesGpu>(primaryVariablesBuffer.data(),
                                                         primaryVariablesBuffer.size()),
        Opm::gpuistl::GpuView<IntensiveQuantitiesGpu>(intensiveQuantitiesBuffer.data(),
                                                      intensiveQuantitiesBuffer.size()),
        numCells);
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
    OPM_GPU_SAFE_CALL(cudaGetLastError());

    std::vector<IntensiveQuantitiesCpu> cpuIntensiveQuantities(numCells);
    for (std::size_t i = 0; i < numCells; ++i) {
        cpuIntensiveQuantities[i].update(cpuProblem, cpuSolution[i], static_cast<unsigned>(i), 0);
    }

    std::vector<IntensiveQuantitiesGpu> gpuIntensiveQuantities(numCells, gpuIntensiveQuantitiesPrototype);
    OPM_GPU_SAFE_CALL(cudaMemcpy(gpuIntensiveQuantities.data(),
                                 intensiveQuantitiesBuffer.data(),
                                 numCells * sizeof(IntensiveQuantitiesGpu),
                                 cudaMemcpyDeviceToHost));

    // Only fields stored in the fluid state are compared: the GPU fluid state's
    // fluid system holds device pointers.
    constexpr double tol = 1e-6;
    constexpr bool checkDerivatives = true;
    constexpr unsigned phases[] = {FluidSystem::waterPhaseIdx,
                                   FluidSystem::oilPhaseIdx,
                                   FluidSystem::gasPhaseIdx};
    std::size_t derivativeComparisons = 0;
    for (std::size_t i = 0; i < numCells; ++i) {
        BOOST_TEST_CONTEXT("cell " << i)
        {
            const auto& cpuIQ = cpuIntensiveQuantities[i];
            const auto& gpuIQ = gpuIntensiveQuantities[i];
            const auto& cpuFs = cpuIQ.fluidState();
            const auto& gpuFs = gpuIQ.fluidState();

            for (const unsigned p : phases) {
                BOOST_TEST_CONTEXT("phase " << p)
                {
                    checkValueAndDerivatives(cpuFs.saturation(p), gpuFs.saturation(p),
                                             tol, checkDerivatives, "saturation",
                                             derivativeComparisons);
                    checkValueAndDerivatives(cpuFs.pressure(p), gpuFs.pressure(p),
                                             tol, checkDerivatives, "pressure",
                                             derivativeComparisons);
                    checkValueAndDerivatives(cpuFs.density(p), gpuFs.density(p),
                                             tol, checkDerivatives, "density",
                                             derivativeComparisons);
                    checkValueAndDerivatives(cpuFs.invB(p), gpuFs.invB(p),
                                             tol, checkDerivatives, "invB",
                                             derivativeComparisons);
                    checkValueAndDerivatives(cpuIQ.mobility(p), gpuIQ.mobility(p),
                                             tol, checkDerivatives, "mobility",
                                             derivativeComparisons);
                }
            }

            checkValueAndDerivatives(cpuFs.Rs(), gpuFs.Rs(),
                                     tol, checkDerivatives, "Rs", derivativeComparisons);
            if (vaporizedOil) {
                // Rv is about 1e-4, so the absolute part of the derivative
                // tolerance would hide errors. Compare Rv relative to its value.
                const double rvScale = asDouble(cpuFs.Rv()) > 0.0 ? asDouble(cpuFs.Rv()) : 1.0;
                checkValueAndDerivatives(cpuFs.Rv() / rvScale, gpuFs.Rv() / rvScale,
                                         tol, checkDerivatives, "Rv", derivativeComparisons);
            }
            else {
                BOOST_CHECK_EQUAL(asDouble(cpuFs.Rv()), 0.0);
                BOOST_CHECK_EQUAL(asDouble(gpuFs.Rv()), 0.0);
            }
            BOOST_CHECK_EQUAL(static_cast<unsigned>(cpuFs.pvtRegionIndex()),
                              static_cast<unsigned>(gpuFs.pvtRegionIndex()));

            checkValueAndDerivatives(cpuIQ.porosity(), gpuIQ.porosity(),
                                     tol, checkDerivatives, "porosity", derivativeComparisons);
            BOOST_CHECK_CLOSE(static_cast<double>(cpuIQ.referencePorosity()),
                              static_cast<double>(gpuIQ.referencePorosity()), tol);
            checkValueAndDerivatives(cpuIQ.rockCompTransMultiplier(),
                                     gpuIQ.rockCompTransMultiplier(),
                                     tol, checkDerivatives, "rockCompTransMultiplier",
                                     derivativeComparisons);
        }
    }
    BOOST_REQUIRE_GT(derivativeComparisons, 0u);
}

} // namespace

BOOST_AUTO_TEST_CASE(DispatcherRejectsEndpointScaling)
{
    // Control: the deck without endpoint scaling is accepted.
    BOOST_CHECK_NO_THROW(runDispatcher(*makeSimulator(deckString)));

    // The GPU material law manager ignores endpoint scaling, so a deck that
    // uses it must be rejected instead of being evaluated with unscaled tables.
    const auto sim = makeSimulator(deckWithEndpointScaling());
    BOOST_REQUIRE(sim->problem().materialLawManager().oilWaterConfig().enableSatScaling());
    BOOST_CHECK_EXCEPTION(runDispatcher(*sim),
                          std::logic_error,
                          [](const std::logic_error& e) {
                              return std::string_view(e.what()).find("endpoint scaling")
                                  != std::string_view::npos;
                          });
}

BOOST_AUTO_TEST_CASE(ThreePhaseLiveOilDryGasGpuVsCpu)
{
    compareGpuWithCpu<Opm::Properties::TTag::FlowBlackoilNohystGPU,
                      Opm::LiveOilDryGasPvtApproaches>(deckString);
}

BOOST_AUTO_TEST_CASE(ThreePhaseLiveOilWetGasGpuVsCpu)
{
    compareGpuWithCpu<Opm::Properties::TTag::FlowBlackoilNohystWetGasGPU,
                      Opm::LiveOilWetGasPvtApproaches>(deckWithWetGas(true));
}
