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
 * \brief Compares the three-phase EclDefaultMaterial law evaluated on the GPU
 *        through EclMaterialLaw::GpuManager against the CPU
 *        EclMaterialLaw::Manager, for an SPE1-like deck with two SATNUM
 *        regions.
 */

#include <config.h>

#define BOOST_TEST_MODULE TestGpuEclMaterialLawManager

#include <cuda_runtime.h>

#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FieldPropsManager.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/fluidmatrixinteractions/EclDefaultMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>
#include <opm/material/fluidmatrixinteractions/EclMaterialLawTwoPhaseTypes.hpp>
#include <opm/material/fluidmatrixinteractions/EclMultiplexerMaterialParams.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/fluidmatrixinteractions/PiecewiseLinearTwoPhaseMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/PiecewiseLinearTwoPhaseMaterialParams.hpp>

#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>

#include <opm/material/fluidmatrixinteractions/GpuEclMaterialLawManager.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

namespace {

// ---------------------------------------------------------------------------
// Deck
// ---------------------------------------------------------------------------

// Four cells in a row with SATNUM 1 2 2 1. Table 1 holds SPE1's SWOF/SGOF
// (from opm-data/spe1/SPE1CASE1.DATA). Both SPE1 capillary pressures are zero,
// so table 2 uses a different Swl, nonzero Pcow and Pcog, and a different
// number of rows to exercise region indexing and the per-region sample buffers.
constexpr std::array<double, 15> table1SwofSw
    = {0.12, 0.18, 0.24, 0.3, 0.36, 0.42, 0.48, 0.54, 0.6, 0.66, 0.72, 0.78, 0.84, 0.91, 1.0};
constexpr std::array<double, 15> table1SgofSg
    = {0.0, 0.001, 0.02, 0.05, 0.12, 0.2, 0.25, 0.3, 0.4, 0.45, 0.5, 0.6, 0.7, 0.85, 0.88};
constexpr std::array<double, 6> table2SwofSw = {0.2, 0.3, 0.45, 0.6, 0.8, 1.0};
constexpr std::array<double, 7> table2SgofSg = {0.0, 0.05, 0.15, 0.3, 0.5, 0.65, 0.8};

const std::string tablesString =
    "SWOF\n"
    "0.12  0                      1       0\n"
    "0.18  4.64876033057851E-008  1       0\n"
    "0.24  0.000000186            0.997   0\n"
    "0.3   4.18388429752066E-007  0.98    0\n"
    "0.36  7.43801652892562E-007  0.7     0\n"
    "0.42  1.16219008264463E-006  0.35    0\n"
    "0.48  1.67355371900826E-006  0.2     0\n"
    "0.54  2.27789256198347E-006  0.09    0\n"
    "0.6   2.97520661157025E-006  0.021   0\n"
    "0.66  3.7654958677686E-006   0.01    0\n"
    "0.72  4.64876033057851E-006  0.001   0\n"
    "0.78  0.000005625            0.0001  0\n"
    "0.84  6.69421487603306E-006  0       0\n"
    "0.91  8.05914256198347E-006  0       0\n"
    "1     0.00001                0       0 /\n"
    "0.2   0.0   1.0   8.0\n"
    "0.3   0.02  0.75  4.0\n"
    "0.45  0.08  0.40  2.0\n"
    "0.6   0.20  0.15  1.0\n"
    "0.8   0.50  0.02  0.3\n"
    "1.0   1.00  0.0   0.0 /\n"
    "\n"
    "SGOF\n"
    "0      0      1       0\n"
    "0.001  0      1       0\n"
    "0.02   0      0.997   0\n"
    "0.05   0.005  0.980   0\n"
    "0.12   0.025  0.700   0\n"
    "0.2    0.075  0.350   0\n"
    "0.25   0.125  0.200   0\n"
    "0.3    0.190  0.090   0\n"
    "0.4    0.410  0.021   0\n"
    "0.45   0.60   0.010   0\n"
    "0.5    0.72   0.001   0\n"
    "0.6    0.87   0.0001  0\n"
    "0.7    0.94   0.000   0\n"
    "0.85   0.98   0.000   0\n"
    "0.88   0.984  0.000   0 /\n"
    "0.0   0.0   1.0   0.0\n"
    "0.05  0.0   0.85  0.1\n"
    "0.15  0.04  0.55  0.4\n"
    "0.3   0.15  0.25  0.9\n"
    "0.5   0.40  0.06  1.6\n"
    "0.65  0.62  0.01  2.2\n"
    "0.8   0.85  0.0   3.0 /\n"
    "\n";

const std::string letTablesString =
    "SWOFLET\n"
    " 0.1  0.2 1.5 7.0 1.5 0.5    0.05 0.15 3.5 3.0 1.3 1.0    0.7 17.0 0.95 3.8 0.04 /\n"
    " 0.1  0.2 1.5 7.0 1.5 0.5    0.05 0.15 3.5 3.0 1.3 1.0    0.7 17.0 0.95 3.8 0.04 /\n"
    "\n"
    "SGOFLET\n"
    " 0.0  0.03 1.8 1.9 1.0 0.95   0.0  0.01 3.5 4.0 1.1 1.0   1.0 1.0 1.0 0.2 0.01 /\n"
    " 0.0  0.03 1.8 1.9 1.0 0.95   0.0  0.01 3.5 4.0 1.1 1.0   1.0 1.0 1.0 0.2 0.01 /\n"
    "\n";

/*!
 * \brief Build the test deck.
 *
 * \param runspecExtra Extra RUNSPEC keywords (e.g. ENDSCALE).
 * \param propsTables  The saturation-function keywords (two tables each).
 * \param propsExtra   Extra PROPS keywords (e.g. STONE1, SWL).
 */
std::string makeDeckString(const std::string& runspecExtra,
                           const std::string& propsTables,
                           const std::string& propsExtra)
{
    return "RUNSPEC\n"
           "DIMENS\n"
           "  4 1 1 /\n"
           "TABDIMS\n"
           "  2 /\n"
           "OIL\n"
           "GAS\n"
           "WATER\n"
           "DISGAS\n"
           "FIELD\n"
        + runspecExtra
        + "GRID\n"
          "DX\n"
          "  4*1000 /\n"
          "DY\n"
          "  4*1000 /\n"
          "DZ\n"
          "  4*20 /\n"
          "TOPS\n"
          "  4*8325 /\n"
          "PORO\n"
          "  4*0.3 /\n"
          "PROPS\n"
        + propsTables + propsExtra
        + "REGIONS\n"
          "SATNUM\n"
          "  1 2 2 1 /\n";
}

const std::string defaultDeckString = makeDeckString("", tablesString, "");

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

constexpr int waterPhaseIdx = 0;
constexpr int oilPhaseIdx = 1;
constexpr int gasPhaseIdx = 2;
constexpr int numPhases = 3;

template <class Scalar>
struct TestTypes {
    // Same traits as flow's black-oil problem (hysteresis and endpoint
    // scaling enabled at compile time), so the CPU manager has the same
    // parameter layering that GpuManager has to unwrap.
    using Traits = Opm::ThreePhaseMaterialTraits<Scalar,
                                                 waterPhaseIdx,
                                                 oilPhaseIdx,
                                                 gasPhaseIdx,
                                                 /*enableHysteresis=*/true,
                                                 /*enableEndpointScaling=*/true>;
    using CpuManager = Opm::EclMaterialLaw::Manager<Traits>;
    using CpuMaterialLaw = typename CpuManager::MaterialLaw;

    using TwoPhaseTypes = Opm::EclMaterialLaw::TwoPhaseTypes<Traits>;
    using GasOilTraits = typename TwoPhaseTypes::GasOilTraits;
    using OilWaterTraits = typename TwoPhaseTypes::OilWaterTraits;

    using GasOilParams = Opm::PiecewiseLinearTwoPhaseMaterialParams<
        GasOilTraits, Opm::gpuistl::GpuView<const Scalar>>;
    using OilWaterParams = Opm::PiecewiseLinearTwoPhaseMaterialParams<
        OilWaterTraits, Opm::gpuistl::GpuView<const Scalar>>;
    using GasOilLaw = Opm::PiecewiseLinearTwoPhaseMaterial<GasOilTraits, GasOilParams>;
    using OilWaterLaw = Opm::PiecewiseLinearTwoPhaseMaterial<OilWaterTraits, OilWaterParams>;

    using GpuManagerBuf = Opm::EclMaterialLaw::
        GpuManager<Traits, GasOilLaw, OilWaterLaw, Opm::gpuistl::GpuBuffer>;
    using GpuManagerView = Opm::EclMaterialLaw::
        GpuManager<Traits, GasOilLaw, OilWaterLaw, Opm::gpuistl::GpuView>;
    using GpuMaterialLaw = typename GpuManagerView::MaterialLaw;

    static_assert(std::is_same_v<GpuMaterialLaw,
                                 Opm::EclDefaultMaterial<Traits, GasOilLaw, OilWaterLaw>>,
                  "GpuManager is expected to default to EclDefaultMaterial");

    // Derivatives w.r.t. Sw (0) and Sg (1), with So = 1 - Sw - Sg, mimicking
    // how the intensive quantities build saturations from primary variables.
    using Evaluation = Opm::DenseAd::Evaluation<Scalar, 2>;
    using PhaseValues = std::array<Evaluation, numPhases>;
};

/*!
 * \brief Minimal fluid state: the three-phase material law only reads
 *        saturations.
 */
template <class Evaluation>
struct SaturationFluidState {
    using ValueType = Evaluation;

    Evaluation saturation_[numPhases];

    OPM_HOST_DEVICE const Evaluation& saturation(int phaseIdx) const
    {
        return saturation_[phaseIdx];
    }
};

template <class Scalar>
struct SaturationPoint {
    unsigned cell;
    Scalar sw;
    Scalar sg;
};

template <class Evaluation, class Scalar>
OPM_HOST_DEVICE SaturationFluidState<Evaluation> makeFluidState(const SaturationPoint<Scalar>& p)
{
    SaturationFluidState<Evaluation> fs;
    fs.saturation_[waterPhaseIdx] = Evaluation::createVariable(p.sw, 0);
    fs.saturation_[gasPhaseIdx] = Evaluation::createVariable(p.sg, 1);
    fs.saturation_[oilPhaseIdx]
        = Evaluation(Scalar(1)) - fs.saturation_[waterPhaseIdx] - fs.saturation_[gasPhaseIdx];
    return fs;
}

// The functions used by the CPU manager's initParamsForElements() to map
// field properties to cells. No LGRs here, so these are the trivial versions
// also used by opm-common's test_eclmateriallawmanager.cpp.
std::vector<int> lookupFieldPropInt(const Opm::FieldPropsManager& fieldProps,
                                    const std::string& keyword,
                                    bool needsTranslation)
{
    const auto& raw = fieldProps.get_int(keyword);
    std::vector<int> dest(raw.size());
    for (std::size_t i = 0; i < raw.size(); ++i) {
        dest[i] = raw[i] - static_cast<int>(needsTranslation);
    }
    return dest;
}

unsigned identityIndex(unsigned elemIdx)
{
    return elemIdx;
}

template <class CpuManager>
std::size_t initCpuManager(CpuManager& manager, const std::string& deckString)
{
    Opm::Parser parser;
    const auto deck = parser.parseString(deckString);
    const Opm::EclipseState eclState(deck);
    const std::size_t numCells = eclState.getInputGrid().getCartesianSize();

    manager.initFromState(eclState);
    manager.initParamsForElements(eclState, numCells, lookupFieldPropInt, identityIndex);
    return numCells;
}

template <class Types>
typename Types::GpuManagerBuf buildGpuManager(const typename Types::CpuManager& cpuManager,
                                              std::size_t numCells)
{
    return typename Types::GpuManagerBuf(cpuManager, numCells);
}

// ---------------------------------------------------------------------------
// Kernels
// ---------------------------------------------------------------------------

template <class Types>
__global__ void evaluateMaterialLawKernel(typename Types::GpuManagerView manager,
                                          const SaturationPoint<typename Types::Traits::Scalar>* points,
                                          std::size_t numPoints,
                                          typename Types::PhaseValues* kr,
                                          typename Types::PhaseValues* pc)
{
    using Evaluation = typename Types::Evaluation;
    using Law = typename Types::GpuMaterialLaw;

    const std::size_t i = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= numPoints) {
        return;
    }
    const auto fs = makeFluidState<Evaluation>(points[i]);
    const auto& params = manager.materialLawParams(points[i].cell);
    Law::relativePermeabilities(kr[i], params, fs);
    Law::capillaryPressures(pc[i], params, fs);
}

template <class ManagerView>
__global__ void satnumRegionKernel(ManagerView manager, int* regions, std::size_t numCells)
{
    const std::size_t i = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < numCells) {
        regions[i] = manager.satnumRegionIdx(static_cast<unsigned>(i));
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

template <class Scalar>
constexpr Scalar tolerance()
{
    // The CPU and the GPU run the same code, so the only expected differences
    // come from FMA contraction on the device.
    if constexpr (std::is_same_v<Scalar, float>) {
        return Scalar(1e-4);
    } else {
        return Scalar(1e-11);
    }
}

template <class Scalar>
bool isClose(Scalar a, Scalar b, double toleranceFactor = 1.0)
{
    const Scalar scale = std::max({Scalar(1), std::abs(a), std::abs(b)});
    return std::abs(a - b) <= tolerance<Scalar>() * toleranceFactor * scale;
}

/*!
 * \brief Tolerance factor for the oil relperm of EclDefaultMaterial.
 *
 * Close to connate water, kro is a weighted average divided by
 * Sw + Sg - Swl. That denominator carries a rounding error of about one
 * ulp of Swl, so the relative error of kro, and even more of its derivatives,
 * grows like 1/(Sw + Sg - Swl). FMA contraction on the device rounds
 * differently from the host, which makes this visible. Below eps/2 = 5e-6 the
 * law no longer divides, so the factor is capped there.
 */
double oilKrToleranceFactor(double sw, double sg, double swl)
{
    const double denominator = std::max(sw, swl) + sg - swl;
    return std::max(1.0, 1e-3 / std::max(denominator, 5e-6));
}

/*!
 * \brief Whether any kind of endpoint scaling is active.
 *
 * Manager::enableEndPointScaling() only reports ENKRVD tables. ENDSCALE is
 * recorded in the per-system EPS configs.
 */
template <class CpuManager>
bool usesEndpointScaling(const CpuManager& manager)
{
    return manager.enableEndPointScaling()
        || manager.oilWaterConfig().enableSatScaling()
        || manager.gasOilConfig().enableSatScaling()
        || manager.gasWaterConfig().enableSatScaling();
}

template <class Scalar>
void addPoint(std::vector<SaturationPoint<Scalar>>& points, unsigned cell, double sw, double sg)
{
    points.push_back({cell, static_cast<Scalar>(sw), static_cast<Scalar>(sg)});
}

/*!
 * \brief Saturation points for one cell, built from that cell's table.
 */
template <class Scalar, class SwBreaks, class SgBreaks>
void addPointsForCell(std::vector<SaturationPoint<Scalar>>& points,
                      unsigned cell,
                      double swl,
                      const SwBreaks& swBreaks,
                      const SgBreaks& sgBreaks)
{
    // Interior grid, including Sw < Swl where the law clamps Sw to Swl.
    for (int i = 0; i <= 20; ++i) {
        for (int j = 0; i + j <= 20; ++j) {
            addPoint<Scalar>(points, cell, 0.05 * i, 0.05 * j);
        }
    }

    // Exact table breakpoints. At a kink, the CPU and the GPU must pick the
    // same one-sided derivative.
    for (const double sw : swBreaks) {
        addPoint<Scalar>(points, cell, sw, 0.0);
        if (sw + 0.1 <= 1.0) {
            addPoint<Scalar>(points, cell, sw, 0.1);
        }
    }
    for (const double sg : sgBreaks) {
        addPoint<Scalar>(points, cell, swl, sg);
    }

    // Connate water and its neighbourhood.
    addPoint<Scalar>(points, cell, swl, 0.0);
    addPoint<Scalar>(points, cell, swl - 0.05, 0.0);
    addPoint<Scalar>(points, cell, 0.0, 0.0);
    addPoint<Scalar>(points, cell, 0.0, 0.3);

    // So = 0.
    addPoint<Scalar>(points, cell, swl, 1.0 - swl);
    addPoint<Scalar>(points, cell, 1.0, 0.0);

    // The krn regularization band of EclDefaultMaterial: Sw + Sg - Swl in
    // [0, eps/2), (eps/2, eps) and just above eps, with eps = 1e-5.
    for (const double sg : {2e-6, 4e-6, 7e-6, 9e-6, 1.2e-5, 3e-5}) {
        addPoint<Scalar>(points, cell, swl, sg);
    }
    addPoint<Scalar>(points, cell, swl + 3e-6, 4e-6);
    addPoint<Scalar>(points, cell, swl - 0.01, 3e-6); // sw is clamped to Swl

    // Slightly outside the physical range.
    addPoint<Scalar>(points, cell, 1.01, 0.0);
    addPoint<Scalar>(points, cell, 0.5, -0.01);
    addPoint<Scalar>(points, cell, -0.01, 0.2);
    addPoint<Scalar>(points, cell, swl, 1.0 - swl + 0.01);
}

/*!
 * \brief Predicate for BOOST_CHECK_EXCEPTION: the exception message contains
 *        \p expected, so the test pins which check in GpuManager rejected the
 *        deck.
 */
auto messageContains(std::string_view expected)
{
    return [expected](const std::logic_error& e) {
        return std::string_view(e.what()).find(expected) != std::string_view::npos;
    };
}

template <class Evaluation>
void checkPhaseValues(const std::array<Evaluation, numPhases>& cpu,
                      const std::array<Evaluation, numPhases>& gpu,
                      const char* quantity,
                      unsigned cell,
                      double sw,
                      double sg,
                      double oilToleranceFactor = 1.0)
{
    for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
        const auto& c = cpu[phaseIdx];
        const auto& g = gpu[phaseIdx];
        const double factor = phaseIdx == oilPhaseIdx ? oilToleranceFactor : 1.0;
        bool ok = isClose(c.value(), g.value(), factor);
        for (int d = 0; d < Evaluation::numVars; ++d) {
            ok = ok && isClose(c.derivative(d), g.derivative(d), factor);
        }
        if (!ok) {
            std::ostringstream msg;
            msg.precision(17);
            msg << quantity << " mismatch for phase " << phaseIdx << " in cell " << cell
                << " at Sw = " << sw << ", Sg = " << sg << ": CPU " << c << ", GPU " << g;
            BOOST_ERROR(msg.str());
        }
    }
}

} // namespace

using ScalarTypes = boost::mpl::list<float, double>;

BOOST_AUTO_TEST_CASE_TEMPLATE(SatnumRegionLookupMatchesCpu, Scalar, ScalarTypes)
{
    using Types = TestTypes<Scalar>;

    typename Types::CpuManager cpuManager;
    const std::size_t numCells = initCpuManager(cpuManager, defaultDeckString);
    BOOST_REQUIRE_EQUAL(numCells, 4u);

    typename Types::GpuManagerBuf gpuManagerBuf(cpuManager, numCells);
    auto gpuManagerView = Opm::gpuistl::make_view(gpuManagerBuf);
    BOOST_CHECK_EQUAL(gpuManagerBuf.materialLawParamsRegionCount(), 2u);

    Opm::gpuistl::GpuBuffer<int> regions(numCells);
    satnumRegionKernel<<<1, 32>>>(gpuManagerView, regions.data(), numCells);
    OPM_GPU_SAFE_CALL(cudaGetLastError());
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());

    const auto hostRegions = regions.asStdVector();
    for (std::size_t cell = 0; cell < numCells; ++cell) {
        BOOST_CHECK_EQUAL(hostRegions[cell], cpuManager.satnumRegionIdx(static_cast<unsigned>(cell)));
    }
    // SATNUM 1 2 2 1, zero-based.
    BOOST_CHECK_EQUAL(hostRegions[0], 0);
    BOOST_CHECK_EQUAL(hostRegions[1], 1);
    BOOST_CHECK_EQUAL(hostRegions[2], 1);
    BOOST_CHECK_EQUAL(hostRegions[3], 0);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(ThreePhaseRelpermAndPcMatchCpu, Scalar, ScalarTypes)
{
    using Types = TestTypes<Scalar>;
    using Evaluation = typename Types::Evaluation;
    using PhaseValues = typename Types::PhaseValues;

    typename Types::CpuManager cpuManager;
    const std::size_t numCells = initCpuManager(cpuManager, defaultDeckString);

    // The GPU manager only supports this configuration. Check that the deck
    // really produces it, so the comparison below is meaningful.
    BOOST_REQUIRE(!usesEndpointScaling(cpuManager));
    BOOST_REQUIRE(!cpuManager.enableHysteresis());
    BOOST_REQUIRE(cpuManager.satCurveIsAllPiecewiseLinear());

    typename Types::GpuManagerBuf gpuManagerBuf(cpuManager, numCells);
    auto gpuManagerView = Opm::gpuistl::make_view(gpuManagerBuf);

    std::vector<SaturationPoint<Scalar>> points;
    std::vector<double> cellSwl(numCells);
    for (unsigned cell = 0; cell < numCells; ++cell) {
        const auto& cpuParams = cpuManager.materialLawParams(cell);
        BOOST_REQUIRE(cpuParams.approach() == Opm::EclMultiplexerApproach::Default);
        const double swl
            = cpuParams.template getRealParams<Opm::EclMultiplexerApproach::Default>().Swl();
        cellSwl[cell] = swl;

        if (cpuManager.satnumRegionIdx(cell) == 0) {
            BOOST_CHECK_CLOSE(swl, table1SwofSw.front(), 1e-4);
            addPointsForCell(points, cell, swl, table1SwofSw, table1SgofSg);
        } else {
            BOOST_CHECK_CLOSE(swl, table2SwofSw.front(), 1e-4);
            addPointsForCell(points, cell, swl, table2SwofSw, table2SgofSg);
        }
    }

    // Evaluate on the GPU.
    const std::size_t numPoints = points.size();
    const Opm::gpuistl::GpuBuffer<SaturationPoint<Scalar>> gpuPoints(points);
    Opm::gpuistl::GpuBuffer<PhaseValues> gpuKr(numPoints);
    Opm::gpuistl::GpuBuffer<PhaseValues> gpuPc(numPoints);

    constexpr unsigned threadsPerBlock = 128;
    const unsigned numBlocks = static_cast<unsigned>((numPoints + threadsPerBlock - 1) / threadsPerBlock);
    evaluateMaterialLawKernel<Types><<<numBlocks, threadsPerBlock>>>(
        gpuManagerView, gpuPoints.data(), numPoints, gpuKr.data(), gpuPc.data());
    OPM_GPU_SAFE_CALL(cudaGetLastError());
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());

    const auto hostKr = gpuKr.asStdVector();
    const auto hostPc = gpuPc.asStdVector();

    // Compare with the CPU multiplexer, i.e. the path flow uses.
    bool anyNonZeroPcow = false;
    bool anyNonZeroPcgo = false;
    for (std::size_t i = 0; i < numPoints; ++i) {
        const auto& p = points[i];
        const auto fs = makeFluidState<Evaluation>(p);
        const auto& cpuParams = cpuManager.materialLawParams(p.cell);

        PhaseValues cpuKr;
        PhaseValues cpuPc;
        Types::CpuMaterialLaw::relativePermeabilities(cpuKr, cpuParams, fs);
        Types::CpuMaterialLaw::capillaryPressures(cpuPc, cpuParams, fs);

        checkPhaseValues(cpuKr, hostKr[i], "kr", p.cell, p.sw, p.sg,
                         oilKrToleranceFactor(p.sw, p.sg, cellSwl[p.cell]));
        checkPhaseValues(cpuPc, hostPc[i], "pc", p.cell, p.sw, p.sg);

        anyNonZeroPcow = anyNonZeroPcow || cpuPc[waterPhaseIdx].value() != 0;
        anyNonZeroPcgo = anyNonZeroPcgo || cpuPc[gasPhaseIdx].value() != 0;
    }

    // Guard against a deck change that would leave the Pc paths untested.
    BOOST_CHECK(anyNonZeroPcow);
    BOOST_CHECK(anyNonZeroPcgo);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(RejectsStone1, Scalar, ScalarTypes)
{
    using Types = TestTypes<Scalar>;

    typename Types::CpuManager cpuManager;
    const std::size_t numCells
        = initCpuManager(cpuManager, makeDeckString("", tablesString, "STONE1\n"));
    BOOST_REQUIRE(cpuManager.materialLawParams(0).approach()
                  == Opm::EclMultiplexerApproach::Stone1);

    BOOST_CHECK_EXCEPTION(buildGpuManager<Types>(cpuManager, numCells),
                          std::logic_error,
                          messageContains("only supports the Default three-phase approach"));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(RejectsNonPiecewiseLinearCurves, Scalar, ScalarTypes)
{
    using Types = TestTypes<Scalar>;

    typename Types::CpuManager cpuManager;
    const std::size_t numCells = initCpuManager(cpuManager, makeDeckString("", letTablesString, ""));
    BOOST_REQUIRE(!cpuManager.satCurveIsAllPiecewiseLinear());

    BOOST_CHECK_EXCEPTION(buildGpuManager<Types>(cpuManager, numCells),
                          std::logic_error,
                          messageContains("only supports piecewise-linear saturation curves"));
}

// With endpoint scaling, the Default params carry the cell's scaled Swl. When
// it varies within one SATNUM region, the per-region GPU storage cannot
// represent it.
//
// Note: GpuManager does not apply endpoint scaling or hysteresis at all; it
// takes the unscaled piecewise-linear tables. A deck with ENDSCALE but the
// same scaled Swl in all cells of a region is accepted, and gives wrong
// values whenever other scaled endpoints (SWCR, SOWCR, KRW, ...) differ from
// the unscaled tables, so the caller must reject such decks.
BOOST_AUTO_TEST_CASE_TEMPLATE(RejectsSwlVaryingWithinRegion, Scalar, ScalarTypes)
{
    using Types = TestTypes<Scalar>;

    // Cells 0 and 3 are both in region 1 but get different SWL.
    typename Types::CpuManager cpuManager;
    const std::size_t numCells = initCpuManager(
        cpuManager,
        makeDeckString("ENDSCALE\n/\n", tablesString, "SWL\n  0.12 0.2 0.2 0.15 /\n"));
    BOOST_REQUIRE(usesEndpointScaling(cpuManager));

    BOOST_CHECK_EXCEPTION(buildGpuManager<Types>(cpuManager, numCells),
                          std::logic_error,
                          messageContains("vary within one SATNUM region"));
}
