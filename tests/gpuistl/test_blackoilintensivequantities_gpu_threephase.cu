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
 * whose gas-oil and oil-water contacts both lie inside the grid, from a
 * variant of it with vaporized oil, and from variants with two- and
 * three-point endpoint scaling and SWATINIT. The GPU side uses the TypeTags of
 * GpuFlowBlackoilNohystTypeTag.hpp, and every cell is compared, values and
 * derivatives.
 *
 * Also checks which configuration the production dispatcher for this TypeTag
 * chooses, and that it rejects decks outside them.
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

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
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
/// part: relative to \p tol percent of the smaller magnitude, as
/// \c BOOST_CHECK_CLOSE, with an absolute floor for values that are zero on
/// one side and round-off on the other. When \p checkDerivatives is true and both operands are
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
    constexpr double absoluteFloor = 1e-12;
    const double cpu = asDouble(cpuValue);
    const double gpu = asDouble(gpuValue);
    BOOST_CHECK_MESSAGE(std::abs(cpu - gpu)
                            <= 1e-2 * tol * std::min(std::abs(cpu), std::abs(gpu)) + absoluteFloor,
                        label << " cpu=" << cpu << " gpu=" << gpu);

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
/// solution of \p sim, as flow does for the first linearization, and return
/// the name of the configuration it chose.
std::string runDispatcher(const Simulator& sim)
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
    BOOST_REQUIRE(dispatcher.configurationName() != nullptr);
    return dispatcher.configurationName();
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

/// The test deck with VAPOIL but dry gas (PVDG).
std::string deckWithVaporizedOilAndDryGas()
{
    std::string deck = deckString;
    addVaporizedOil(deck);
    return deck;
}

/// The test deck with dead oil (PVDO) instead of live oil.
std::string deckWithDeadOil()
{
    std::string deck = deckString;
    replaceOnce(deck, "DISGAS\n", "");
    replaceSection(deck, "PVTO\n", "REGIONS\n",
                   "PVDO\n"
                   "  50.0 1.10 1.20\n"
                   " 150.0 1.08 1.25\n"
                   " 350.0 1.05 1.35 /\n");
    replaceSection(deck, "RSVD\n", "SCHEDULE\n", "");
    replaceOnce(deck, " 2030 200 2045 0 2015 0 1 0 0 /", " 2030 200 2045 0 2015 0 0 0 0 /");
    replaceOnce(deck, "DRSDT\n 0 /\n", "");
    return deck;
}

bool contains(std::string_view text, std::string_view part)
{
    return text.find(part) != std::string_view::npos;
}

/// \p keyword with the nine values of \p layer in each of the six layers.
std::string sameInEachLayer(std::string_view keyword, std::string_view layer)
{
    std::string text = std::string(keyword) + "\n";
    for (int k = 0; k < 6; ++k) {
        text += " ";
        text += layer;
        text += "\n";
    }
    return text + "/\n";
}

void insertAfter(std::string& deck, std::string_view anchor, std::string_view text)
{
    const auto pos = deck.find(anchor);
    BOOST_REQUIRE(pos != std::string::npos);
    deck.insert(pos + anchor.size(), text);
}

/// Add ENDSCALE and per-cell saturation endpoints to \p deck. With
/// \p threePointAndSwatinit, also add SCALECRS YES and SWATINIT.
///
/// The endpoints are the same in every layer. SATNUM is 1 1 2 in each row of
/// a layer, so the two region 1 cells of a row, and the region 2 cells of
/// different rows, have different endpoints. All cells pass the saturation
/// function consistency checks, including those of three-point scaling. One
/// cell has a nonzero SGL, and its SWU is below 1 - SGL. SGU is 1 - SWL, as
/// in the tables, so the gas cap of the wet gas deck keeps cells without oil.
///
/// SWATINIT is set per layer, below every SWL in the top layer. Above the
/// water-oil contact it gives each cell a new maximum Pcow in applySwatinit();
/// in the water zone Pcow is zero and the maximum is unchanged.
void addEndpointScaling(std::string& deck, bool threePointAndSwatinit)
{
    insertAfter(deck, "METRIC\n", "ENDSCALE\n/\n");
    std::string props
        = sameInEachLayer("SWL", "0.12 0.17 0.18  0.14 0.16 0.22  0.13 0.18 0.20")
        + sameInEachLayer("SWCR", "0.16 0.20 0.24  0.18 0.19 0.28  0.15 0.23 0.25")
        + sameInEachLayer("SWU", "1.0 0.98 1.0  0.97 0.97 1.0  1.0 0.99 0.98")
        + sameInEachLayer("SGL", "0.0 0.0 0.0  0.0 0.01 0.0  0.0 0.0 0.0")
        + sameInEachLayer("SGCR", "0.03 0.05 0.04  0.02 0.06 0.07  0.04 0.03 0.05")
        + sameInEachLayer("SGU", "0.88 0.83 0.82  0.86 0.84 0.78  0.87 0.82 0.80")
        + sameInEachLayer("SOWCR", "0.14 0.18 0.22  0.16 0.12 0.19  0.17 0.15 0.24")
        + sameInEachLayer("SOGCR", "0.10 0.12 0.18  0.13 0.09 0.21  0.14 0.11 0.16");
    if (threePointAndSwatinit) {
        props += "SCALECRS\n YES /\n"
                 "SWATINIT\n"
                 " 9*0.10 9*0.22 9*0.30 9*0.45 9*0.75 9*1.0 /\n";
    }
    insertAfter(deck, "PROPS\n", props);
}

/// The test deck with two-point endpoint scaling.
std::string deckWithTwoPointScaling()
{
    std::string deck = deckString;
    addEndpointScaling(deck, /*threePointAndSwatinit=*/false);
    return deck;
}

/// The test deck with three-point endpoint scaling and SWATINIT, the
/// combination that Norne uses.
std::string deckWithThreePointScalingAndSwatinit()
{
    std::string deck = deckString;
    addEndpointScaling(deck, /*threePointAndSwatinit=*/true);
    return deck;
}

/// The wet gas deck with vaporized oil, three-point endpoint scaling and
/// SWATINIT.
std::string deckWithWetGasAndThreePointScalingAndSwatinit()
{
    std::string deck = deckWithWetGas(/*vaporizedOil=*/true);
    addEndpointScaling(deck, /*threePointAndSwatinit=*/true);
    return deck;
}

/// Whether \p sw is strictly inside the scaled mobile water range of a cell
/// with endpoints \p info, and \p sg inside its scaled mobile gas range.
template <class EpsInfo>
bool insideScaledWaterRange(const EpsInfo& info, double sw)
{
    return info.Swl < sw && sw < 1.0 - info.Sowcr;
}

template <class EpsInfo>
bool insideScaledGasRange(const EpsInfo& info, double sg)
{
    return info.Sgcr < sg && sg < info.Sgu;
}

/// The water and gas saturations of \p priVars, read from the primary
/// variables.
template <class PrimaryVariables>
std::pair<double, double> waterAndGasSaturation(const PrimaryVariables& priVars)
{
    using GasMeaning = typename PrimaryVariables::GasMeaning;
    using WaterMeaning = typename PrimaryVariables::WaterMeaning;
    using Indices = Opm::GetPropType<CpuTypeTag, Opm::Properties::Indices>;
    const double sw = priVars.primaryVarsMeaningWater() == WaterMeaning::Sw
        ? static_cast<double>(priVars[Indices::waterSwitchIdx])
        : 0.0;
    switch (priVars.primaryVarsMeaningGas()) {
    case GasMeaning::Sg:
        return {sw, static_cast<double>(priVars[Indices::compositionSwitchIdx])};
    case GasMeaning::Rv:
        return {sw, 1.0 - sw};
    default:
        return {sw, 0.0};
    }
}

/// Make sure that some cell has its water saturation, and some cell its gas
/// saturation, strictly inside the scaled mobile range, so that the
/// comparison covers the scaled segments of the curves. Equilibration is
/// used as it is if it already gives such cells. Otherwise one cell of each
/// SATNUM region with the Sg meaning and gas is moved into both ranges.
template <class PrimaryVariables, class MaterialLawManager>
void ensureCellsInsideScaledRanges(std::vector<PrimaryVariables>& solution,
                                   const MaterialLawManager& materialLawManager)
{
    bool waterCovered = false;
    bool gasCovered = false;
    for (std::size_t i = 0; i < solution.size(); ++i) {
        const auto& info = materialLawManager.oilWaterScaledEpsInfoDrainage(i);
        const auto [sw, sg] = waterAndGasSaturation(solution[i]);
        waterCovered = waterCovered || insideScaledWaterRange(info, sw);
        gasCovered = gasCovered || insideScaledGasRange(info, sg);
    }
    BOOST_TEST_MESSAGE("Equilibration gives a cell inside the scaled water range: "
                       << waterCovered << ", gas range: " << gasCovered);
    if (waterCovered && gasCovered) {
        return;
    }

    using GasMeaning = typename PrimaryVariables::GasMeaning;
    using WaterMeaning = typename PrimaryVariables::WaterMeaning;
    using Indices = Opm::GetPropType<CpuTypeTag, Opm::Properties::Indices>;
    constexpr int numRegions = 2;
    std::vector<bool> moved(numRegions, false);
    for (std::size_t i = 0; i < solution.size(); ++i) {
        auto& priVars = solution[i];
        const int region = materialLawManager.satnumRegionIdx(static_cast<unsigned>(i));
        if (moved[region]
            || priVars.primaryVarsMeaningGas() != GasMeaning::Sg
            || priVars.primaryVarsMeaningWater() != WaterMeaning::Sw
            || !(priVars[Indices::compositionSwitchIdx] > 0.0)) {
            continue;
        }
        // Both saturations are set, even if only one range is missing: a
        // gas-cap cell has Sg near Sgu, which leaves no room for a larger Sw.
        // Sw + Sg stays below 1 for the endpoints of the test decks.
        const auto& info = materialLawManager.oilWaterScaledEpsInfoDrainage(i);
        priVars[Indices::waterSwitchIdx] = info.Swl + 0.25 * (1.0 - info.Sowcr - info.Swl);
        priVars[Indices::compositionSwitchIdx] = 0.5 * (info.Sgcr + info.Sgu);
        moved[region] = true;
    }
    for (int region = 0; region < numRegions; ++region) {
        BOOST_TEST_INFO("SATNUM region " << region + 1);
        BOOST_REQUIRE(moved[region]);
    }
}

/// Whether the water saturation of a cell lies on its scaled maximum SWU,
/// where the scaled curves end. Equilibration puts the cells of the water
/// zone there. Mapping SWU back to the end of the table is a ratio that the
/// CPU and the GPU may round to different sides of the end, which gives the
/// same values but different one-sided derivatives. The lower endpoints map
/// exactly and need no exception.
template <class EpsInfo>
bool onScaledMaximumWaterSaturation(const EpsInfo& info, double sw)
{
    return std::abs(sw - info.Swu) < 1e-10;
}

/// Check that \p deck exercises endpoint scaling, so that the GPU-CPU
/// comparison cannot pass on unscaled curves: the scaled endpoints vary
/// within a SATNUM region, some cells have water and gas saturations strictly
/// inside their scaled ranges, and with Pc scaling, SWATINIT has changed the
/// maximum Pcow of some cell.
template <class MaterialLawManager, class IntensiveQuantities>
void checkEndpointScalingCoverage(const MaterialLawManager& materialLawManager,
                                  const std::vector<IntensiveQuantities>& intensiveQuantities)
{
    const std::size_t numCells = intensiveQuantities.size();
    bool swlVariesWithinRegion = false;
    bool maxPcowScaled = false;
    std::size_t numSwInside = 0;
    std::size_t numSgInside = 0;
    for (std::size_t i = 0; i < numCells; ++i) {
        const auto& info = materialLawManager.oilWaterScaledEpsInfoDrainage(i);
        const int region = materialLawManager.satnumRegionIdx(static_cast<unsigned>(i));
        for (std::size_t j = 0; j < i; ++j) {
            if (materialLawManager.satnumRegionIdx(static_cast<unsigned>(j)) == region
                && materialLawManager.oilWaterScaledEpsInfoDrainage(j).Swl != info.Swl) {
                swlVariesWithinRegion = true;
            }
        }
        if (info.maxPcow != materialLawManager.unscaledEpsInfo(region).maxPcow) {
            maxPcowScaled = true;
        }

        const auto& fs = intensiveQuantities[i].fluidState();
        const double sw = asDouble(fs.saturation(FluidSystem::waterPhaseIdx));
        const double sg = asDouble(fs.saturation(FluidSystem::gasPhaseIdx));
        if (insideScaledWaterRange(info, sw)) {
            ++numSwInside;
        }
        if (insideScaledGasRange(info, sg)) {
            ++numSgInside;
        }
    }
    BOOST_TEST_MESSAGE("Endpoint scaling: Sw inside its scaled range in " << numSwInside
                       << " cells, Sg in " << numSgInside << " cells");
    BOOST_CHECK(swlVariesWithinRegion);
    BOOST_CHECK_GT(numSwInside, 0u);
    BOOST_CHECK_GT(numSgInside, 0u);
    if (materialLawManager.oilWaterConfig().enablePcScaling()) {
        BOOST_CHECK(maxPcowScaled);
    }
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
    const auto& cpuMaterialLawManager = cpuProblem.materialLawManager();
    const bool endpointScaling = cpuMaterialLawManager.oilWaterConfig().enableSatScaling();

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
    BOOST_CHECK_EQUAL(gpuProblemBuffer.materialLawManager().hasPerCellEndpoints(),
                      endpointScaling);

    const std::size_t numCells = cpuProblem.model().numGridDof();
    BOOST_REQUIRE_EQUAL(numCells, 54u);

    // Primary variables
    using PrimaryVariablesCpu = Opm::GetPropType<CpuTypeTag, Opm::Properties::PrimaryVariables>;
    using PrimaryVariablesGpu = Opm::GetPropType<GpuTypeTag, Opm::Properties::PrimaryVariables>;
    const auto& initialSolution = cpuProblem.model().solution(/*timeIdx=*/0);
    BOOST_REQUIRE_EQUAL(initialSolution.size(), numCells);
    std::vector<PrimaryVariablesCpu> cpuSolution;
    cpuSolution.reserve(numCells);
    for (std::size_t i = 0; i < numCells; ++i) {
        cpuSolution.push_back(initialSolution[i]);
    }
    if (endpointScaling) {
        ensureCellsInsideScaledRanges(cpuSolution, cpuMaterialLawManager);
    }

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
    std::size_t numCellsOnScaledEndpoint = 0;
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

            // At the end of the scaled curves, only the values of the
            // saturation-dependent quantities are compared.
            const bool onKink = endpointScaling
                && onScaledMaximumWaterSaturation(
                    cpuMaterialLawManager.oilWaterScaledEpsInfoDrainage(i),
                    asDouble(cpuFs.saturation(FluidSystem::waterPhaseIdx)));
            if (onKink) {
                ++numCellsOnScaledEndpoint;
            }
            const bool checkPhaseDerivatives = checkDerivatives && !onKink;

            for (const unsigned p : phases) {
                BOOST_TEST_CONTEXT("phase " << p)
                {
                    checkValueAndDerivatives(cpuFs.saturation(p), gpuFs.saturation(p),
                                             tol, checkPhaseDerivatives, "saturation",
                                             derivativeComparisons);
                    checkValueAndDerivatives(cpuFs.pressure(p), gpuFs.pressure(p),
                                             tol, checkPhaseDerivatives, "pressure",
                                             derivativeComparisons);
                    checkValueAndDerivatives(cpuFs.density(p), gpuFs.density(p),
                                             tol, checkPhaseDerivatives, "density",
                                             derivativeComparisons);
                    checkValueAndDerivatives(cpuFs.invB(p), gpuFs.invB(p),
                                             tol, checkPhaseDerivatives, "invB",
                                             derivativeComparisons);
                    checkValueAndDerivatives(cpuIQ.mobility(p), gpuIQ.mobility(p),
                                             tol, checkPhaseDerivatives, "mobility",
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

    if (endpointScaling) {
        BOOST_TEST_MESSAGE("Cells at the scaled SWU, without phase derivative checks: "
                           << numCellsOnScaledEndpoint << " of " << numCells);
        checkEndpointScalingCoverage(cpuMaterialLawManager, cpuIntensiveQuantities);
    }
}

} // namespace

BOOST_AUTO_TEST_CASE(DispatcherAcceptsEndpointScaling)
{
    BOOST_CHECK(contains(runDispatcher(*makeSimulator(deckWithThreePointScalingAndSwatinit())),
                         "dry gas"));
    BOOST_CHECK(contains(
        runDispatcher(*makeSimulator(deckWithWetGasAndThreePointScalingAndSwatinit())),
        "wet gas"));
}

BOOST_AUTO_TEST_CASE(DispatcherRejectsUnsupportedEndpointScaling)
{
    // The GPU material law manager rejects the options it does not implement.
    // Each has its own case in test_gpu_ecl_material_law_manager; this one
    // shows that the message reaches the user through the dispatcher.
    std::string deck = deckWithTwoPointScaling();
    insertAfter(deck, "PROPS\n", "KRW\n 54*0.9 /\n");
    const auto sim = makeSimulator(deck);
    BOOST_CHECK_EXCEPTION(runDispatcher(*sim),
                          std::logic_error,
                          [](const std::logic_error& e) {
                              BOOST_TEST_MESSAGE(e.what());
                              return contains(e.what(), "relative permeability value scaling");
                          });
}

BOOST_AUTO_TEST_CASE(DispatcherChoosesConfigurationFromPvt)
{
    // Each case starts from its own deck: the CPU fluid system is static.
    BOOST_CHECK(contains(runDispatcher(*makeSimulator(deckString)), "dry gas"));
    BOOST_CHECK(contains(runDispatcher(*makeSimulator(deckWithWetGas(true))), "wet gas"));
    BOOST_CHECK(contains(runDispatcher(*makeSimulator(deckString)), "dry gas"));
}

BOOST_AUTO_TEST_CASE(DispatcherRejectsUnsupportedPvt)
{
    const auto throwsWith = [](const std::string& deck, std::string_view part) {
        const auto sim = makeSimulator(deck);
        BOOST_CHECK_EXCEPTION(runDispatcher(*sim),
                              std::logic_error,
                              [part](const std::logic_error& e) {
                                  BOOST_TEST_MESSAGE(e.what());
                                  return contains(e.what(), part);
                              });
    };
    // Dry gas configuration, which rejects VAPOIL.
    throwsWith(deckWithVaporizedOilAndDryGas(), "VAPOIL");
    // Wet gas configuration, which requires VAPOIL.
    throwsWith(deckWithWetGas(false), "without VAPOIL");
    // No configuration has dead oil.
    throwsWith(deckWithDeadOil(), "PVDO");
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

BOOST_AUTO_TEST_CASE(TwoPointScalingDryGasGpuVsCpu)
{
    compareGpuWithCpu<Opm::Properties::TTag::FlowBlackoilNohystGPU,
                      Opm::LiveOilDryGasPvtApproaches>(deckWithTwoPointScaling());
}

// The simulator applies SWATINIT during equilibration, before the GPU problem
// is built from it, as in flow. This is the in-process check that the GPU
// material law manager is a snapshot taken after applySwatinit().
BOOST_AUTO_TEST_CASE(ThreePointScalingAndSwatinitDryGasGpuVsCpu)
{
    compareGpuWithCpu<Opm::Properties::TTag::FlowBlackoilNohystGPU,
                      Opm::LiveOilDryGasPvtApproaches>(deckWithThreePointScalingAndSwatinit());
}

BOOST_AUTO_TEST_CASE(ThreePointScalingAndSwatinitWetGasGpuVsCpu)
{
    compareGpuWithCpu<Opm::Properties::TTag::FlowBlackoilNohystWetGasGPU,
                      Opm::LiveOilWetGasPvtApproaches>(
        deckWithWetGasAndThreePointScalingAndSwatinit());
}
