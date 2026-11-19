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
 *        regions, with and without endpoint scaling.
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
#include <opm/material/fluidmatrixinteractions/EclEpsScalingPoints.hpp>
#include <opm/material/fluidmatrixinteractions/EclEpsTwoPhaseLaw.hpp>
#include <opm/material/fluidmatrixinteractions/EclEpsTwoPhaseLawParams.hpp>
#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>
#include <opm/material/fluidmatrixinteractions/EclMaterialLawTwoPhaseTypes.hpp>
#include <opm/material/fluidmatrixinteractions/EclMultiplexerMaterialParams.hpp>
#include <opm/material/fluidmatrixinteractions/EclTwoPhaseMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/EclTwoPhaseMaterialParams.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/fluidmatrixinteractions/PiecewiseLinearTwoPhaseMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/PiecewiseLinearTwoPhaseMaterialParams.hpp>

#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>
#include <opm/simulators/linalg/gpuistl/gpu_smart_pointer.hpp>

#include <opm/material/fluidmatrixinteractions/GpuEclMaterialLawManager.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
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
 * \param gridExtra    Extra GRID keywords (e.g. JFUNC).
 * \param regionsExtra Extra REGIONS keywords (e.g. KRNUMX).
 */
std::string makeDeckString(const std::string& runspecExtra,
                           const std::string& propsTables,
                           const std::string& propsExtra,
                           const std::string& gridExtra = "",
                           const std::string& regionsExtra = "")
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
        + gridExtra
        + "PROPS\n"
        + propsTables + propsExtra
        + "REGIONS\n"
          "SATNUM\n"
          "  1 2 2 1 /\n"
        + regionsExtra;
}

const std::string defaultDeckString = makeDeckString("", tablesString, "");

// Per-cell endpoints for the endpoint-scaling decks, chosen so that both cells
// of each region differ. Cell 0 has SWCR = SWL, which makes the first
// segment of the three-point oil-water krn scaling degenerate. Cell 2 has a
// nonzero SGL. The values stay off the 0.05 grid of addPointsForCell, so the
// grid does not hit the kinks at the scaled endpoints (see
// markPointsOnKinks).
const std::string endpointArrays =
    "SWL\n   0.11  0.23  0.17  0.14 /\n"
    "SWCR\n  0.11  0.31  0.26  0.21 /\n"
    "SWU\n   0.99  0.98  0.96  0.965 /\n"
    "SGL\n   0.0   0.0   0.015 0.0 /\n"
    "SGCR\n  0.025 0.055 0.085 0.035 /\n"
    "SGU\n   0.87  0.76  0.74  0.83 /\n"
    "SOWCR\n 0.155 0.055 0.105 0.205 /\n"
    "SOGCR\n 0.105 0.055 0.125 0.115 /\n";

const std::string endscale = "ENDSCALE\n/\n";
const std::string twoPointDeckString = makeDeckString(endscale, tablesString, endpointArrays);
const std::string threePointDeckString
    = makeDeckString(endscale, tablesString, endpointArrays + "SCALECRS\n  YES /\n");
const std::string swatinitDeckString = makeDeckString(
    endscale,
    tablesString,
    endpointArrays + "SCALECRS\n  YES /\nSWATINIT\n  0.3 0.4 0.5 0.5 /\n");
// PCW and PCG enable capillary-pressure scaling without SWATINIT. Table 1
// has zero Pc, which can only be scaled to zero.
const std::string pcwPcgDeckString = makeDeckString(
    endscale,
    tablesString,
    endpointArrays + "PCW\n  0.0 12.0 6.0 0.0 /\nPCG\n  0.0 4.0 2.0 0.0 /\n");

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
 * \brief As TestTypes, but the GPU two-phase laws have the endpoint-scaling
 *        layer, with its default parameter type.
 */
template <class Scalar>
struct EpsTestTypes : TestTypes<Scalar> {
    using Base = TestTypes<Scalar>;
    using Traits = typename Base::Traits;

    using GasOilLaw = Opm::EclEpsTwoPhaseLaw<typename Base::GasOilLaw>;
    using OilWaterLaw = Opm::EclEpsTwoPhaseLaw<typename Base::OilWaterLaw>;

    using GpuManagerBuf = Opm::EclMaterialLaw::
        GpuManager<Traits, GasOilLaw, OilWaterLaw, Opm::gpuistl::GpuBuffer>;
    using GpuManagerView = Opm::EclMaterialLaw::
        GpuManager<Traits, GasOilLaw, OilWaterLaw, Opm::gpuistl::GpuView>;
    using GpuMaterialLaw = typename GpuManagerView::MaterialLaw;

    static_assert(std::is_same_v<GpuMaterialLaw,
                                 Opm::EclDefaultMaterial<Traits, GasOilLaw, OilWaterLaw>>,
                  "GpuManager is expected to default to EclDefaultMaterial");
};

/*!
 * \brief The two-phase (gas-water) GPU manager, as flow's
 *        GpuFlowGasWaterEnergyTypeTags builds it: EclTwoPhaseMaterial over
 *        piecewise-linear laws, without endpoint-scaling layer.
 */
template <class Scalar>
struct GasWaterTestTypes {
    using Traits = typename TestTypes<Scalar>::Traits;
    using CpuManager = Opm::EclMaterialLaw::Manager<Traits>;
    using TwoPhaseTypes = Opm::EclMaterialLaw::TwoPhaseTypes<Traits>;

    template <class TwoPhaseTraits>
    using GpuLaw = Opm::PiecewiseLinearTwoPhaseMaterial<
        TwoPhaseTraits,
        Opm::PiecewiseLinearTwoPhaseMaterialParams<TwoPhaseTraits,
                                                   Opm::gpuistl::GpuView<const Scalar>>>;
    using GasOilLaw = GpuLaw<typename TwoPhaseTypes::GasOilTraits>;
    using OilWaterLaw = GpuLaw<typename TwoPhaseTypes::OilWaterTraits>;
    using GasWaterLaw = GpuLaw<typename TwoPhaseTypes::GasWaterTraits>;

    using GpuMaterialLawParams = Opm::EclTwoPhaseMaterialParams<Traits,
                                                                typename GasOilLaw::Params,
                                                                typename OilWaterLaw::Params,
                                                                typename GasWaterLaw::Params,
                                                                Opm::gpuistl::ValueAsPointer>;
    using GpuMaterialLaw = Opm::EclTwoPhaseMaterial<Traits,
                                                    GasOilLaw,
                                                    OilWaterLaw,
                                                    GasWaterLaw,
                                                    GpuMaterialLawParams>;
    using GpuManagerBuf = Opm::EclMaterialLaw::
        GpuManager<Traits, GasOilLaw, OilWaterLaw, Opm::gpuistl::GpuBuffer, GpuMaterialLaw>;
};

/*!
 * \brief What materialLawParams() of one cell holds, as seen on the device.
 */
template <class Scalar>
struct AssembledParams {
    Opm::EclEpsScalingPoints<Scalar> gasOilScaled;
    Opm::EclEpsScalingPoints<Scalar> oilWaterScaled;
    Opm::EclEpsScalingPoints<Scalar> gasOilUnscaled;
    Opm::EclEpsScalingPoints<Scalar> oilWaterUnscaled;
    Scalar swl;
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
    //! At a kink whose side depends on rounding: compare values only.
    bool valueOnly = false;
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

template <class ManagerView, class Scalar>
__global__ void assembledParamsKernel(ManagerView manager,
                                      AssembledParams<Scalar>* out,
                                      std::size_t numCells)
{
    const std::size_t i = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < numCells) {
        const auto params = manager.materialLawParams(static_cast<unsigned>(i));
        out[i].gasOilScaled = params.gasOilParams().scaledPoints();
        out[i].oilWaterScaled = params.oilWaterParams().scaledPoints();
        out[i].gasOilUnscaled = params.gasOilParams().unscaledPoints();
        out[i].oilWaterUnscaled = params.oilWaterParams().unscaledPoints();
        out[i].swl = params.Swl();
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
bool isClose(Scalar a, Scalar b, double toleranceFactor = 1.0, Scalar absoluteSlack = 0)
{
    const Scalar scale = std::max({Scalar(1), std::abs(a), std::abs(b)});
    return std::abs(a - b) <= tolerance<Scalar>() * toleranceFactor * scale + absoluteSlack;
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

template <class Evaluation, class Scalar = typename Evaluation::ValueType>
void checkPhaseValues(const std::array<Evaluation, numPhases>& cpu,
                      const std::array<Evaluation, numPhases>& gpu,
                      const char* quantity,
                      unsigned cell,
                      double sw,
                      double sg,
                      double oilToleranceFactor = 1.0,
                      bool compareDerivatives = true,
                      bool saturationIsMapped = false)
{
    for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
        const auto& c = cpu[phaseIdx];
        const auto& g = gpu[phaseIdx];
        const double factor = phaseIdx == oilPhaseIdx ? oilToleranceFactor : 1.0;
        // With endpoint scaling, the saturation goes through a linear map
        // before the table lookup, which the device rounds differently (FMA).
        // That error of a few ulps of the saturation is amplified by the
        // slope of the curve. At a kink the CPU and the GPU may report
        // different one-sided slopes, so the steeper one is used. Without
        // the map there is no such error, and no slack is given.
        Scalar slack = 0;
        if (saturationIsMapped) {
            constexpr int saturationUlps = 64;
            Scalar slope = 0;
            for (int d = 0; d < Evaluation::numVars; ++d) {
                slope = std::max({slope, std::abs(c.derivative(d)), std::abs(g.derivative(d))});
            }
            slack = slope * saturationUlps * std::numeric_limits<Scalar>::epsilon();
        }
        bool ok = isClose(c.value(), g.value(), factor, slack);
        for (int d = 0; compareDerivatives && d < Evaluation::numVars; ++d) {
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

template <class CpuManager>
const auto& cpuDefaultParams(const CpuManager& cpuManager, unsigned cell)
{
    const auto& cpuParams = cpuManager.materialLawParams(cell);
    BOOST_REQUIRE(cpuParams.approach() == Opm::EclMultiplexerApproach::Default);
    return cpuParams.template getRealParams<Opm::EclMultiplexerApproach::Default>();
}

// The CPU traits enable hysteresis, so the endpoint-scaling parameters sit
// below the hysteresis layer.
template <class CpuManager>
const auto& cpuGasOilEpsParams(const CpuManager& cpuManager, unsigned cell)
{
    return cpuDefaultParams(cpuManager, cell).gasOilParams().drainageParams();
}

template <class CpuManager>
const auto& cpuOilWaterEpsParams(const CpuManager& cpuManager, unsigned cell)
{
    return cpuDefaultParams(cpuManager, cell).oilWaterParams().drainageParams();
}

/*!
 * \brief The scaled endpoints of one cell, sorted: the kinks of its scaled
 *        curves, in the saturation each two-phase law sees.
 */
struct ScaledKinks {
    std::vector<double> oilWater; //!< in Sw
    std::vector<double> gasOil;   //!< in 1 - Swl - Sg
};

template <class CpuManager>
ScaledKinks scaledKinks(const CpuManager& cpuManager, unsigned cell)
{
    auto collect = [](const auto& scaledPoints) {
        std::vector<double> values;
        for (int i = 0; i < 3; ++i) {
            values.push_back(scaledPoints.saturationPcPoints()[i]);
            values.push_back(scaledPoints.saturationKrwPoints()[i]);
            values.push_back(scaledPoints.saturationKrnPoints()[i]);
        }
        std::sort(values.begin(), values.end());
        values.erase(std::unique(values.begin(), values.end()), values.end());
        return values;
    };
    return {collect(cpuOilWaterEpsParams(cpuManager, cell).scaledPoints()),
            collect(cpuGasOilEpsParams(cpuManager, cell).scaledPoints())};
}

/*!
 * \brief Add points at the given kinks, on both sides of them, and inside
 *        each interval between neighbours.
 *
 * The interior point is off-centre: the middle of a scaled interval maps to
 * the middle of the unscaled one, which is often a table breakpoint.
 */
template <class Scalar>
void addPointsAroundKinks(std::vector<SaturationPoint<Scalar>>& points,
                          unsigned cell,
                          const std::vector<double>& kinks,
                          double otherSaturation,
                          bool saturationIsWater,
                          double swl)
{
    constexpr double offset = 1e-4;
    auto add = [&](double value) {
        if (saturationIsWater) {
            addPoint<Scalar>(points, cell, value, otherSaturation);
        } else {
            addPoint<Scalar>(points, cell, swl, 1.0 - swl - value);
        }
    };

    for (std::size_t i = 0; i < kinks.size(); ++i) {
        add(kinks[i]);
        add(kinks[i] - offset);
        add(kinks[i] + offset);
        if (i + 1 < kinks.size()) {
            add(kinks[i] + 0.37 * (kinks[i + 1] - kinks[i]));
        }
    }
}

/*!
 * \brief Saturation points at and around the scaled endpoints of one cell.
 *
 * Three-point scaling has kinks at the scaled critical saturations, which a
 * sweep over the unscaled table breakpoints can miss.
 */
template <class Scalar>
void addScaledEndpointPoints(std::vector<SaturationPoint<Scalar>>& points,
                             const ScaledKinks& kinks,
                             unsigned cell,
                             double swl)
{
    for (const double sg : {0.0, 0.07}) {
        addPointsAroundKinks(points, cell, kinks.oilWater, sg, /*saturationIsWater=*/true, swl);
    }
    addPointsAroundKinks(points, cell, kinks.gasOil, 0.0, /*saturationIsWater=*/false, swl);
}

/*!
 * \brief Mark the points that lie on a kink of the scaled curves.
 *
 * The scaled-to-unscaled map rounds differently on the device (FMA
 * contraction), so on a kink the CPU and the GPU may take different
 * one-sided derivatives. Only the values are compared there. Sw below Swl
 * is clamped to Swl, which is a kink too.
 */
template <class Scalar>
void markPointsOnKinks(std::vector<SaturationPoint<Scalar>>& points,
                       const std::vector<ScaledKinks>& kinks,
                       const std::vector<double>& cellSwl)
{
    constexpr double distance = 1e-6;
    auto near = [](double value, const std::vector<double>& candidates) {
        return std::any_of(candidates.begin(), candidates.end(), [value](double k) {
            return std::abs(value - k) < distance;
        });
    };
    for (auto& p : points) {
        const double swl = cellSwl[p.cell];
        const double sw = std::max(static_cast<double>(p.sw), swl);
        p.valueOnly = near(sw, kinks[p.cell].oilWater)
            || near(1.0 - swl - static_cast<double>(p.sg), kinks[p.cell].gasOil);
    }
}

/*!
 * \brief Saturation points from the cell's unscaled tables (see
 *        addPointsForCell), and the cells' Swl.
 */
template <class Scalar, class CpuManager>
void addTablePoints(std::vector<SaturationPoint<Scalar>>& points,
                    std::vector<double>& cellSwl,
                    const CpuManager& cpuManager,
                    std::size_t numCells)
{
    cellSwl.resize(numCells);
    for (unsigned cell = 0; cell < numCells; ++cell) {
        const double swl = cpuDefaultParams(cpuManager, cell).Swl();
        cellSwl[cell] = swl;
        if (cpuManager.satnumRegionIdx(cell) == 0) {
            addPointsForCell(points, cell, swl, table1SwofSw, table1SgofSg);
        } else {
            addPointsForCell(points, cell, swl, table2SwofSw, table2SgofSg);
        }
    }
}

struct PcCoverage {
    bool anyNonZeroPcow = false;
    bool anyNonZeroPcgo = false;
};

/*!
 * \brief Evaluate the material law at \p points on the GPU and compare with
 *        the CPU manager, values and derivatives.
 *
 * \p saturationIsMapped: whether the law maps scaled to unscaled saturations
 * (endpoint scaling active), which allows a slope-dependent value slack (see
 * checkPhaseValues).
 */
template <class Types>
PcCoverage compareWithCpu(const typename Types::CpuManager& cpuManager,
                          const typename Types::GpuManagerView& gpuManagerView,
                          const std::vector<SaturationPoint<typename Types::Traits::Scalar>>& points,
                          const std::vector<double>& cellSwl,
                          bool saturationIsMapped)
{
    using Scalar = typename Types::Traits::Scalar;
    using Evaluation = typename Types::Evaluation;
    using PhaseValues = typename Types::PhaseValues;

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
    PcCoverage coverage;
    for (std::size_t i = 0; i < numPoints; ++i) {
        const auto& p = points[i];
        const auto fs = makeFluidState<Evaluation>(p);
        const auto& cpuParams = cpuManager.materialLawParams(p.cell);

        PhaseValues cpuKr;
        PhaseValues cpuPc;
        Types::CpuMaterialLaw::relativePermeabilities(cpuKr, cpuParams, fs);
        Types::CpuMaterialLaw::capillaryPressures(cpuPc, cpuParams, fs);

        checkPhaseValues(cpuKr, hostKr[i], "kr", p.cell, p.sw, p.sg,
                         oilKrToleranceFactor(p.sw, p.sg, cellSwl[p.cell]), !p.valueOnly,
                         saturationIsMapped);
        checkPhaseValues(cpuPc, hostPc[i], "pc", p.cell, p.sw, p.sg, 1.0, !p.valueOnly,
                         saturationIsMapped);

        coverage.anyNonZeroPcow = coverage.anyNonZeroPcow || cpuPc[waterPhaseIdx].value() != 0;
        coverage.anyNonZeroPcgo = coverage.anyNonZeroPcgo || cpuPc[gasPhaseIdx].value() != 0;
    }
    return coverage;
}

template <class Scalar>
void checkSaturationPoints(const std::array<Scalar, 3>& cpu,
                           const std::array<Scalar, 3>& gpu,
                           const char* what,
                           unsigned cell)
{
    for (int i = 0; i < 3; ++i) {
        BOOST_TEST_INFO("cell " << cell << ", " << what << "[" << i << "]");
        BOOST_CHECK_EQUAL(cpu[i], gpu[i]);
    }
}

template <class Scalar>
void checkScalingPoints(const Opm::EclEpsScalingPoints<Scalar>& cpu,
                        const Opm::EclEpsScalingPoints<Scalar>& gpu,
                        const Opm::EclEpsConfig& config,
                        const char* system,
                        unsigned cell)
{
    if (!config.enableSatScaling()) {
        return;
    }
    BOOST_TEST_CONTEXT(system)
    {
        checkSaturationPoints(cpu.saturationPcPoints(), gpu.saturationPcPoints(), "pc", cell);
        checkSaturationPoints(cpu.saturationKrwPoints(), gpu.saturationKrwPoints(), "krw", cell);
        checkSaturationPoints(cpu.saturationKrnPoints(), gpu.saturationKrnPoints(), "krn", cell);
        if (config.enablePcScaling()) {
            BOOST_TEST_INFO("cell " << cell << ", maxPcnw");
            BOOST_CHECK_EQUAL(cpu.maxPcnw(), gpu.maxPcnw());
        }
    }
}

/*!
 * \brief Check the GPU manager built from \p cpuManager against it: the
 *        layout, the assembled parameters of every cell, and the material
 *        law at the table breakpoints and the scaled endpoints.
 *
 * The buffer manager is moved before the view is made, so pointers into a
 * moved-from manager would show.
 */
template <class Types>
PcCoverage checkEpsManagerMatchesCpu(const typename Types::CpuManager& cpuManager,
                                     std::size_t numCells,
                                     bool expectPerCellEndpoints)
{
    using Scalar = typename Types::Traits::Scalar;

    BOOST_REQUIRE(!cpuManager.enableHysteresis());
    BOOST_REQUIRE(cpuManager.satCurveIsAllPiecewiseLinear());

    typename Types::GpuManagerBuf builtManager(cpuManager, numCells);
    typename Types::GpuManagerBuf gpuManagerBuf(std::move(builtManager));
    auto gpuManagerView = Opm::gpuistl::make_view(gpuManagerBuf);

    BOOST_CHECK_EQUAL(gpuManagerBuf.hasPerCellEndpoints(), expectPerCellEndpoints);
    BOOST_CHECK_EQUAL(gpuManagerBuf.numMaterialLawRegions(), 2u);

    // SATNUM lookup.
    Opm::gpuistl::GpuBuffer<int> regions(numCells);
    satnumRegionKernel<<<1, 32>>>(gpuManagerView, regions.data(), numCells);
    OPM_GPU_SAFE_CALL(cudaGetLastError());
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
    const auto hostRegions = regions.asStdVector();
    for (unsigned cell = 0; cell < numCells; ++cell) {
        BOOST_CHECK_EQUAL(hostRegions[cell], cpuManager.satnumRegionIdx(cell));
    }

    // Assembled parameters, field by field. They are copies of the CPU
    // values, so they must match exactly.
    Opm::gpuistl::GpuBuffer<AssembledParams<Scalar>> assembled(numCells);
    assembledParamsKernel<<<1, 32>>>(gpuManagerView, assembled.data(), numCells);
    OPM_GPU_SAFE_CALL(cudaGetLastError());
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
    const auto hostAssembled = assembled.asStdVector();
    for (unsigned cell = 0; cell < numCells; ++cell) {
        const auto& gpu = hostAssembled[cell];
        const auto& cpuGasOil = cpuGasOilEpsParams(cpuManager, cell);
        const auto& cpuOilWater = cpuOilWaterEpsParams(cpuManager, cell);
        BOOST_TEST_INFO("cell " << cell);
        BOOST_CHECK_EQUAL(gpu.swl, cpuDefaultParams(cpuManager, cell).Swl());
        checkScalingPoints(cpuGasOil.scaledPoints(), gpu.gasOilScaled,
                           cpuGasOil.config(), "gas-oil scaled", cell);
        checkScalingPoints(cpuOilWater.scaledPoints(), gpu.oilWaterScaled,
                           cpuOilWater.config(), "oil-water scaled", cell);
        checkScalingPoints(cpuGasOil.unscaledPoints(), gpu.gasOilUnscaled,
                           cpuGasOil.config(), "gas-oil unscaled", cell);
        checkScalingPoints(cpuOilWater.unscaledPoints(), gpu.oilWaterUnscaled,
                           cpuOilWater.config(), "oil-water unscaled", cell);
    }

    // The material law.
    std::vector<SaturationPoint<Scalar>> points;
    std::vector<double> cellSwl;
    addTablePoints(points, cellSwl, cpuManager, numCells);
    if (expectPerCellEndpoints) {
        std::vector<ScaledKinks> kinks;
        for (unsigned cell = 0; cell < numCells; ++cell) {
            kinks.push_back(scaledKinks(cpuManager, cell));
            addScaledEndpointPoints(points, kinks.back(), cell, cellSwl[cell]);
        }
        markPointsOnKinks(points, kinks, cellSwl);
    }
    return compareWithCpu<Types>(
        cpuManager, gpuManagerView, points, cellSwl, /*saturationIsMapped=*/expectPerCellEndpoints);
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
    BOOST_CHECK_EQUAL(gpuManagerBuf.numMaterialLawRegions(), 2u);

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
    std::vector<double> cellSwl;
    addTablePoints(points, cellSwl, cpuManager, numCells);
    for (unsigned cell = 0; cell < numCells; ++cell) {
        const double tableSwl = cpuManager.satnumRegionIdx(cell) == 0 ? table1SwofSw.front()
                                                                       : table2SwofSw.front();
        BOOST_CHECK_CLOSE(cellSwl[cell], tableSwl, 1e-4);
    }

    const auto coverage = compareWithCpu<Types>(
        cpuManager, gpuManagerView, points, cellSwl, /*saturationIsMapped=*/false);

    // Guard against a deck change that would leave the Pc paths untested.
    BOOST_CHECK(coverage.anyNonZeroPcow);
    BOOST_CHECK(coverage.anyNonZeroPcgo);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(RejectsStone1, Scalar, ScalarTypes)
{
    using Types = TestTypes<Scalar>;

    typename Types::CpuManager cpuManager;
    const std::size_t numCells
        = initCpuManager(cpuManager, makeDeckString("", tablesString, "STONE1\n"));
    BOOST_REQUIRE(cpuManager.materialLawParams(0).approach()
                  == Opm::EclMultiplexerApproach::Stone1);

    BOOST_CHECK_THROW(buildGpuManager<Types>(cpuManager, numCells),
                      std::logic_error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(RejectsNonPiecewiseLinearCurves, Scalar, ScalarTypes)
{
    using Types = TestTypes<Scalar>;

    typename Types::CpuManager cpuManager;
    const std::size_t numCells = initCpuManager(cpuManager, makeDeckString("", letTablesString, ""));
    BOOST_REQUIRE(!cpuManager.satCurveIsAllPiecewiseLinear());

    BOOST_CHECK_THROW(buildGpuManager<Types>(cpuManager, numCells),
                      std::logic_error);
}

// The material law of TestTypes has no endpoint-scaling layer, so any
// endpoint scaling must be rejected rather than silently ignored.
BOOST_AUTO_TEST_CASE_TEMPLATE(RejectsEndpointScalingWithoutEpsLayer, Scalar, ScalarTypes)
{
    using Types = TestTypes<Scalar>;

    typename Types::CpuManager cpuManager;
    const std::size_t numCells = initCpuManager(
        cpuManager,
        makeDeckString("ENDSCALE\n/\n", tablesString, "SWL\n  0.12 0.2 0.2 0.15 /\n"));
    BOOST_REQUIRE(usesEndpointScaling(cpuManager));

    BOOST_CHECK_THROW(buildGpuManager<Types>(cpuManager, numCells),
                      std::logic_error);
}

namespace {

struct UnsupportedOption {
    const char* name;
    //! RUNSPEC keywords of the base deck, which the GPU manager accepts.
    std::string baseRunspec;
    //! The keywords that add the option to the base deck.
    std::string runspecExtra;
    std::string propsExtra;
    std::string gridExtra;
    std::string regionsExtra;
};

std::vector<UnsupportedOption> unsupportedOptions()
{
    return {
        {"hysteresis", "", "SATOPTS\n  HYSTER /\n", "EHYSTR\n  0.1 0 /\n", "", ""},
        {"directional relperms", "", "SATOPTS\n  DIRECT /\n", "", "", "KRNUMX\n  1 2 2 1 /\n"},
        {"kr value scaling", endscale, "", "KRW\n  4*0.9 /\n", "", ""},
        {"JFUNC", endscale, "", "",
         "PERMX\n  4*100 /\nPERMY\n  4*100 /\nPERMZ\n  4*10 /\nJFUNC\n  WATER 22.0 /\n", ""},
        {"ENKRVD", endscale, "",
         "ENKRVD\n"
         "  8000 1.0 1.0 1.0 0.9 0.9 0.9 0.9\n"
         "  9000 1.0 1.0 1.0 0.9 0.9 0.9 0.9 /\n",
         "", ""},
        {"PPCWMAX", "", "", "PPCWMAX\n  100.0 NO /\n  100.0 NO /\n", "", ""},
    };
}

} // namespace

// Each option the GPU manager does not implement is rejected. The law has the
// endpoint-scaling layer, so the rejection of endpoint scaling without that
// layer cannot fire, and each case depends on the check for its own option.
// The base deck without the option is accepted, so the rejection is caused by
// the option alone.
BOOST_AUTO_TEST_CASE_TEMPLATE(RejectsUnsupportedOptions, Scalar, ScalarTypes)
{
    using Types = EpsTestTypes<Scalar>;

    for (const auto& option : unsupportedOptions()) {
        BOOST_TEST_CONTEXT("Option: " << option.name)
        {
            {
                typename Types::CpuManager cpuManager;
                const std::size_t numCells
                    = initCpuManager(cpuManager, makeDeckString(option.baseRunspec, tablesString, ""));
                BOOST_CHECK_NO_THROW(buildGpuManager<Types>(cpuManager, numCells));
            }

            typename Types::CpuManager cpuManager;
            const std::size_t numCells
                = initCpuManager(cpuManager,
                                 makeDeckString(option.baseRunspec + option.runspecExtra,
                                                tablesString,
                                                option.propsExtra,
                                                option.gridExtra,
                                                option.regionsExtra));
            BOOST_CHECK_THROW(buildGpuManager<Types>(cpuManager, numCells),
                              std::logic_error);
        }
    }
}

// ---------------------------------------------------------------------------
// Two-phase (gas-water) path
// ---------------------------------------------------------------------------

namespace {

/*!
 * \brief A WATER+GAS deck, which the CPU manager runs with the GasWater
 *        two-phase approach.
 */
std::string makeGasWaterDeckString(const std::string& runspecExtra, const std::string& propsExtra)
{
    return "RUNSPEC\n"
           "DIMENS\n"
           "  4 1 1 /\n"
           "TABDIMS\n"
           "  1 /\n"
           "WATER\n"
           "GAS\n"
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
          "SWFN\n"
          "  0.12  0.0  4.0\n"
          "  0.5   0.3  1.0\n"
          "  1.0   1.0  0.0 /\n"
          "SGFN\n"
          "  0.0   0.0  0.0\n"
          "  0.5   0.3  0.0\n"
          "  0.88  1.0  0.0 /\n"
        + propsExtra;
}

} // namespace

// The gas-water path is the one flow's GPU property evaluation uses. Its law
// has no endpoint-scaling layer, so the manager itself must reject a deck
// with ENDSCALE, whatever checks its caller makes.
BOOST_AUTO_TEST_CASE_TEMPLATE(TwoPhaseRejectsEndpointScaling, Scalar, ScalarTypes)
{
    using Types = GasWaterTestTypes<Scalar>;

    // Without endpoint scaling the deck is accepted, so the rejection below
    // is caused by ENDSCALE alone.
    {
        typename Types::CpuManager cpuManager;
        const std::size_t numCells = initCpuManager(cpuManager, makeGasWaterDeckString("", ""));
        BOOST_REQUIRE(cpuManager.twoPhaseApproach() == Opm::EclTwoPhaseApproach::GasWater);
        BOOST_REQUIRE(!usesEndpointScaling(cpuManager));
        BOOST_CHECK_NO_THROW(buildGpuManager<Types>(cpuManager, numCells));
    }

    typename Types::CpuManager cpuManager;
    const std::size_t numCells = initCpuManager(
        cpuManager,
        makeGasWaterDeckString("ENDSCALE\n/\n", "SWL\n  0.12 0.15 0.18 0.2 /\n"));
    BOOST_REQUIRE(cpuManager.twoPhaseApproach() == Opm::EclTwoPhaseApproach::GasWater);
    BOOST_REQUIRE(cpuManager.gasWaterConfig().enableSatScaling());
    BOOST_REQUIRE(!cpuManager.enableEndPointScaling());

    BOOST_CHECK_THROW(buildGpuManager<Types>(cpuManager, numCells),
                      std::logic_error);
}

// ---------------------------------------------------------------------------
// Endpoint scaling
// ---------------------------------------------------------------------------

// The EPS-capable GPU law matches the CPU on a deck without endpoint scaling.
BOOST_AUTO_TEST_CASE_TEMPLATE(EpsLayerWithoutEndpointScalingMatchesCpu, Scalar, ScalarTypes)
{
    using Types = EpsTestTypes<Scalar>;

    typename Types::CpuManager cpuManager;
    const std::size_t numCells = initCpuManager(cpuManager, defaultDeckString);
    BOOST_REQUIRE(!usesEndpointScaling(cpuManager));

    const auto coverage
        = checkEpsManagerMatchesCpu<Types>(cpuManager, numCells, /*expectPerCellEndpoints=*/false);
    BOOST_CHECK(coverage.anyNonZeroPcow);
    BOOST_CHECK(coverage.anyNonZeroPcgo);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(TwoPointScalingMatchesCpu, Scalar, ScalarTypes)
{
    using Types = EpsTestTypes<Scalar>;

    typename Types::CpuManager cpuManager;
    const std::size_t numCells = initCpuManager(cpuManager, twoPointDeckString);
    BOOST_REQUIRE(cpuManager.oilWaterConfig().enableSatScaling());
    BOOST_REQUIRE(!cpuManager.oilWaterConfig().enableThreePointKrSatScaling());
    BOOST_REQUIRE(!cpuManager.oilWaterConfig().enablePcScaling());

    checkEpsManagerMatchesCpu<Types>(cpuManager, numCells, /*expectPerCellEndpoints=*/true);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(ThreePointScalingMatchesCpu, Scalar, ScalarTypes)
{
    using Types = EpsTestTypes<Scalar>;

    typename Types::CpuManager cpuManager;
    const std::size_t numCells = initCpuManager(cpuManager, threePointDeckString);
    BOOST_REQUIRE(cpuManager.oilWaterConfig().enableThreePointKrSatScaling());
    BOOST_REQUIRE(cpuManager.gasOilConfig().enableThreePointKrSatScaling());

    checkEpsManagerMatchesCpu<Types>(cpuManager, numCells, /*expectPerCellEndpoints=*/true);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(SwatinitScalingMatchesCpu, Scalar, ScalarTypes)
{
    using Types = EpsTestTypes<Scalar>;

    typename Types::CpuManager cpuManager;
    const std::size_t numCells = initCpuManager(cpuManager, swatinitDeckString);
    BOOST_REQUIRE(cpuManager.oilWaterConfig().enablePcScaling());
    BOOST_REQUIRE_EQUAL(numCells, 4u);

    // As the equilibration would, before the GPU manager is built. Both
    // table 2 cells get a new maximum Pcow, each its own. Table 1 has zero
    // Pcow, so cell 0 stays below applySwatinit's threshold, and cell 3 has
    // pcow <= 0; neither changes.
    struct SwatinitCall {
        unsigned cell;
        Scalar pcow;
        Scalar sw;
        bool changesMaxPcow;
    };
    const std::array<SwatinitCall, 4> calls{{
        {0, Scalar(1.0e4), Scalar(0.3), false},
        {1, Scalar(3.0e4), Scalar(0.4), true},
        {2, Scalar(2.0e4), Scalar(0.5), true},
        {3, Scalar(-1.0), Scalar(0.5), false},
    }};
    for (const auto& call : calls) {
        const Scalar maxPcowBefore = cpuManager.oilWaterScaledEpsInfoDrainage(call.cell).maxPcow;
        cpuManager.applySwatinit(call.cell, call.pcow, call.sw);
        const Scalar maxPcowAfter = cpuManager.oilWaterScaledEpsInfoDrainage(call.cell).maxPcow;
        BOOST_TEST_INFO("cell " << call.cell);
        BOOST_CHECK_EQUAL(maxPcowAfter != maxPcowBefore, call.changesMaxPcow);
    }
    // The two cells of region 2 must differ, so the per-cell storage is needed.
    BOOST_CHECK_NE(cpuManager.oilWaterScaledEpsInfoDrainage(1).maxPcow,
                   cpuManager.oilWaterScaledEpsInfoDrainage(2).maxPcow);

    const auto coverage
        = checkEpsManagerMatchesCpu<Types>(cpuManager, numCells, /*expectPerCellEndpoints=*/true);
    BOOST_CHECK(coverage.anyNonZeroPcow);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(PcwPcgScalingMatchesCpu, Scalar, ScalarTypes)
{
    using Types = EpsTestTypes<Scalar>;

    typename Types::CpuManager cpuManager;
    const std::size_t numCells = initCpuManager(cpuManager, pcwPcgDeckString);
    BOOST_REQUIRE(cpuManager.oilWaterConfig().enablePcScaling());
    BOOST_REQUIRE(cpuManager.gasOilConfig().enablePcScaling());

    const auto coverage
        = checkEpsManagerMatchesCpu<Types>(cpuManager, numCells, /*expectPerCellEndpoints=*/true);
    BOOST_CHECK(coverage.anyNonZeroPcow);
    BOOST_CHECK(coverage.anyNonZeroPcgo);
}

// The device constructor of EclEpsTwoPhaseLawParams, used on the host with
// the CPU effective law: it must give the same law as the setter-built
// parameters it copies.
BOOST_AUTO_TEST_CASE_TEMPLATE(EpsParamsDeviceConstructorMatchesSetters, Scalar, ScalarTypes)
{
    using Types = TestTypes<Scalar>;
    using Law = typename Types::TwoPhaseTypes::OilWaterEpsLaw;
    using Params = typename Law::Params;

    typename Types::CpuManager cpuManager;
    const std::size_t numCells = initCpuManager(cpuManager, threePointDeckString);

    for (unsigned cell = 0; cell < numCells; ++cell) {
        const Params& setterBuilt = cpuOilWaterEpsParams(cpuManager, cell);
        const Params constructed(&setterBuilt.effectiveLawParams(),
                                 &setterBuilt.unscaledPoints(),
                                 setterBuilt.config(),
                                 setterBuilt.scaledPoints());
        for (int i = 0; i <= 50; ++i) {
            const Scalar sw = Scalar(0.02 * i);
            BOOST_TEST_CONTEXT("cell " << cell << ", Sw = " << sw)
            {
                BOOST_CHECK_EQUAL(Law::twoPhaseSatPcnw(constructed, sw),
                                  Law::twoPhaseSatPcnw(setterBuilt, sw));
                BOOST_CHECK_EQUAL(Law::twoPhaseSatKrw(constructed, sw),
                                  Law::twoPhaseSatKrw(setterBuilt, sw));
                BOOST_CHECK_EQUAL(Law::twoPhaseSatKrn(constructed, sw),
                                  Law::twoPhaseSatKrn(setterBuilt, sw));
            }
        }
    }
}
