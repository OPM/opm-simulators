// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
 * \brief Unit tests for WetGasPvt on the GPU.
 *
 * A CPU WetGasPvt is initialized from a deck with two PVT regions of different
 * table sizes, copied to the GPU with copy_to_gpu()/make_view(), and evaluated
 * in kernels. All GPU results are compared to the same evaluations on the CPU.
 *
 * Unlike LiveOilPvt, the 2D tables of WetGasPvt have pressure on the X axis and Rv on
 * the Y axis, and they use the RightExtreme interpolation guide, so these tests are
 * also the first to run that branch of UniformXTabulated2DFunction::findPoints() on
 * the GPU.
 */
#include <config.h>

#define BOOST_TEST_MODULE TestGpuWetGasPvt
#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/fluidsystems/blackoilpvt/WetGasPvt.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Python/Python.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

template <class Scalar>
using Evaluation = Opm::DenseAd::Evaluation<Scalar, 3>;

template <class Scalar>
using CpuPvt = Opm::WetGasPvt<Scalar>;
template <class Scalar>
using GpuBufPvt = Opm::WetGasPvt<Scalar, Opm::gpuistl::GpuBuffer, Opm::gpuistl::GpuView>;
template <class Scalar>
using GpuViewPvt = Opm::WetGasPvt<Scalar, Opm::gpuistl::GpuView, Opm::gpuistl::GpuView>;

// WetGasPvt<float> is explicitly instantiated in opm-common independently of
// FLOW_INSTANTIATE_FLOAT, so both types are always tested.
using Scalars = boost::mpl::list<double, float>;

template <class Scalar>
struct CheckCopySemantics {
    // The views in the owning GPU object point into the tables of its device table owner,
    // so a copy would have views into the tables of the original. The owning GPU object
    // must be move-only.
    static_assert(!std::is_copy_constructible_v<GpuBufPvt<Scalar>>);
    static_assert(!std::is_copy_assignable_v<GpuBufPvt<Scalar>>);
    static_assert(std::is_move_constructible_v<GpuBufPvt<Scalar>>);
    static_assert(std::is_move_assignable_v<GpuBufPvt<Scalar>>);
    // The CPU object is still copied by GasPvtMultiplexer, and the view is passed by value.
    static_assert(std::is_copy_constructible_v<CpuPvt<Scalar>>);
    static_assert(std::is_trivially_copyable_v<GpuViewPvt<Scalar>>);
};
template struct CheckCopySemantics<double>;
template struct CheckCopySemantics<float>;

// Two PVT regions with different numbers of pressure rows (3 and 8) and different
// numbers of undersaturated points per row. Rows with only the saturated point are
// extended by WetGasPvt::extendPvtgTable_() from the next row with undersaturated
// points; the Rv steps of those rows are small enough that the extended Rv stay >= 0.
//
// In WetGasPvt the 2D tables have the pressure on the X axis and Rv on the Y axis:
// - Region 1 has eight pressure rows, so the bisection in the pressure direction
//   iterates twice and takes both sides of the pivot comparison.
// - Region 1's row at PG = 120 has eight Rv points, so the bisection in the Rv
//   direction does the same when that row brackets the pressure.
//
// The BG values keep the first two rows of each region from triggering the
// low-pressure padding of PVTG tables in the parser, so the tables are as given here.
constexpr const char* deckString = R"(
RUNSPEC
DIMENS
 1 1 1 /
TABDIMS
 * 2 /
OIL
GAS
WATER
VAPOIL
METRIC
GRID
DX
 1*100 /
DY
 1*100 /
DZ
 1*10 /
TOPS
 1*1000 /
PORO
 1*0.2 /
PROPS
DENSITY
 850.0 1030.0 0.90 /
 870.0 1020.0 0.80 /
PVTG
-- PG     RV        BG       MUG
   50.0   0.00010   0.0200   0.0140
          0.00000   0.0198   0.0136 /
  100.0   0.00025   0.0120   0.0160 /
  150.0   0.00040   0.0085   0.0185
          0.00030   0.0084   0.0182
          0.00020   0.0083   0.0179 /
/
   40.0   0.00002   0.0250   0.0130
          0.00000   0.0249   0.0129 /
   60.0   0.00012   0.0180   0.0138
          0.00000   0.0178   0.0135 /
   90.0   0.00025   0.0125   0.0150
          0.00010   0.0123   0.0147 /
  120.0   0.00038   0.0098   0.0163
          0.00030   0.0097   0.0161
          0.00022   0.0096   0.0159
          0.00016   0.0095   0.0157
          0.00010   0.0094   0.0155
          0.00006   0.0093   0.0153
          0.00003   0.0092   0.0151
          0.00000   0.0091   0.0149 /
  150.0   0.00052   0.0082   0.0178 /
  180.0   0.00065   0.0072   0.0192
          0.00030   0.0070   0.0185 /
  210.0   0.00075   0.0065   0.0205 /
  240.0   0.00082   0.0060   0.0218
          0.00050   0.0058   0.0210
          0.00020   0.0056   0.0202 /
/
)";

constexpr unsigned numRegions = 2;
// [K], not used by WetGasPvt. A function rather than a constexpr variable: the PVT
// methods take the temperature by reference, which would require the variable to
// exist in device memory.
template <class Scalar>
OPM_HOST_DEVICE constexpr Scalar temperature()
{ return 350.0; }

// Relative tolerance for comparing CPU and GPU results. The GPU may contract
// multiplications and additions into FMAs, so results can differ in the last bits.
template <class Scalar>
constexpr Scalar relTol = std::is_same_v<Scalar, float> ? 1e-4 : 1e-10;

// Relative tolerance for derivatives. A derivative can be a small difference of much
// larger terms, e.g. d(mu)/dp with mu = invB / invBMu where mu hardly depends on p
// (extrapolated far above the saturated Rv). In float, the rounding of those terms then
// shows up much larger relative to the result. In double, the tolerance stays tight,
// which still catches a different segment or branch on the CPU and the GPU.
template <class Scalar>
constexpr Scalar derivativeRelTol = std::is_same_v<Scalar, float> ? 1e-3 : 1e-10;

template <class Scalar>
CpuPvt<Scalar> makeCpuPvt()
{
    const auto deck = Opm::Parser{}.parseString(deckString);
    const Opm::EclipseState eclState(deck);
    const Opm::Schedule schedule(deck, eclState, std::make_shared<Opm::Python>());

    CpuPvt<Scalar> pvt;
    pvt.initFromState(eclState, schedule);
    return pvt;
}

//! One evaluation point. Pressure in [Pa], Rv in [sm3/sm3].
template <class Scalar>
struct TestPoint {
    unsigned regionIdx;
    Scalar p;
    Scalar rv;
    Scalar so;
    Scalar soMax;
};

//! The parts of a fluid state used by WetGasPvt::inverseFormationVolumeFactorAndViscosity().
template <class Scalar>
struct TestFluidState {
    using ValueType = Evaluation<Scalar>;
    static constexpr int oilPhaseIdx = 0;
    static constexpr int gasPhaseIdx = 1;

    OPM_HOST_DEVICE const ValueType& pressure(unsigned) const
    { return p; }

    OPM_HOST_DEVICE ValueType Rv() const
    { return rv; }

    OPM_HOST_DEVICE const ValueType& saturation(unsigned phaseIdx) const
    { return phaseIdx == static_cast<unsigned>(oilPhaseIdx) ? so : sg; }

    ValueType p;
    ValueType rv;
    ValueType so;
    ValueType sg;
};

//! Result of inverseFormationVolumeFactorAndViscosity(). A plain struct instead of
//! std::pair, because it is copied between host and device with memcpy.
template <class Scalar>
struct InvBAndMu {
    Evaluation<Scalar> invB;
    Evaluation<Scalar> mu;
};

//! \p value as an \c Eval: the plain value if \c Eval is a scalar type, otherwise an AD
//! variable with index \p varIdx.
template <class Eval, class Scalar>
OPM_HOST_DEVICE Eval makeVariable(Scalar value, int varIdx)
{
    if constexpr (std::is_same_v<Eval, Scalar>) {
        return value;
    } else {
        return Eval::createVariable(value, varIdx);
    }
}

// Functors evaluating one WetGasPvt method. They are templated on the PVT type, so the
// same functor is used for the CPU reference and in the kernel. The functors templated
// on Eval pass the state to the method as Eval, so the method is evaluated with plain
// scalars or with AD variables (pressure: 0, Rv or So: 1, SoMax: 2), as the fluid
// system does.

template <class Eval>
struct InverseFormationVolumeFactor {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE Eval operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    {
        return pvt.inverseFormationVolumeFactor(pt.regionIdx, Eval(temperature<Scalar>()),
                                                makeVariable<Eval>(pt.p, 0),
                                                makeVariable<Eval>(pt.rv, 1),
                                                /*Rvw=*/Eval(0.0));
    }
};

template <class Eval>
struct Viscosity {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE Eval operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    {
        return pvt.viscosity(pt.regionIdx, Eval(temperature<Scalar>()),
                             makeVariable<Eval>(pt.p, 0),
                             makeVariable<Eval>(pt.rv, 1),
                             /*Rvw=*/Eval(0.0));
    }
};

template <class Eval>
struct SaturatedInverseFormationVolumeFactor {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE Eval operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    {
        return pvt.saturatedInverseFormationVolumeFactor(pt.regionIdx, Eval(temperature<Scalar>()),
                                                         makeVariable<Eval>(pt.p, 0));
    }
};

template <class Eval>
struct SaturatedViscosity {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE Eval operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    {
        return pvt.saturatedViscosity(pt.regionIdx, Eval(temperature<Scalar>()),
                                      makeVariable<Eval>(pt.p, 0));
    }
};

template <class Eval>
struct SaturatedOilVaporizationFactor {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE Eval operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    {
        return pvt.saturatedOilVaporizationFactor(pt.regionIdx, Eval(temperature<Scalar>()),
                                                  makeVariable<Eval>(pt.p, 0));
    }
};

template <class Eval>
struct SaturatedOilVaporizationFactorVappars {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE Eval operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    {
        return pvt.saturatedOilVaporizationFactor(pt.regionIdx, Eval(temperature<Scalar>()),
                                                  makeVariable<Eval>(pt.p, 0),
                                                  makeVariable<Eval>(pt.so, 1),
                                                  makeVariable<Eval>(pt.soMax, 2));
    }
};

//! The methods that are constant zero for wet gas (no vaporized water), summed. Each
//! term is zero, so the sum is zero exactly if all of them are.
struct WetGasZeros {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE Evaluation<Scalar> operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    {
        using Eval = Evaluation<Scalar>;
        const Eval T(temperature<Scalar>());
        const Eval p = makeVariable<Eval>(pt.p, 0);
        return pvt.saturatedWaterVaporizationFactor(pt.regionIdx, T, p)
            + pvt.saturatedWaterVaporizationFactor(pt.regionIdx, T, p, Eval(0.0));
    }
};

struct GasReferenceDensity {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE Scalar operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    { return pvt.gasReferenceDensity(pt.regionIdx); }
};

struct OilReferenceDensity {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE Scalar operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    { return pvt.oilReferenceDensity(pt.regionIdx); }
};

struct VapPar1 {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE Scalar operator()(const Pvt& pvt, const TestPoint<Scalar>&) const
    { return pvt.vapPar1(); }
};

struct NumRegions {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE unsigned operator()(const Pvt& pvt, const TestPoint<Scalar>&) const
    { return pvt.numRegions(); }
};

//! The hot path, with pressure, Rv and So as AD variables.
struct InverseFormationVolumeFactorAndViscosity {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE InvBAndMu<Scalar> operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    {
        using Eval = Evaluation<Scalar>;
        TestFluidState<Scalar> fs;
        fs.p = Eval::createVariable(pt.p, 0);
        fs.rv = Eval::createVariable(pt.rv, 1);
        fs.so = Eval::createVariable(pt.so, 2);
        fs.sg = Eval(1 - pt.so);
        const auto [invB, mu] = pvt.inverseFormationVolumeFactorAndViscosity(fs, pt.regionIdx);
        return {invB, mu};
    }
};

// Every method used on the GPU must be callable on a const view, since the fluid
// system holding the PVT object may only be available through a const reference.
template <class Scalar>
concept CallableOnConstView = requires(const GpuViewPvt<Scalar>& pvt,
                                       const TestPoint<Scalar>& pt) {
    InverseFormationVolumeFactor<Evaluation<Scalar>>{}(pvt, pt);
    Viscosity<Evaluation<Scalar>>{}(pvt, pt);
    SaturatedInverseFormationVolumeFactor<Evaluation<Scalar>>{}(pvt, pt);
    SaturatedViscosity<Evaluation<Scalar>>{}(pvt, pt);
    SaturatedOilVaporizationFactor<Evaluation<Scalar>>{}(pvt, pt);
    SaturatedOilVaporizationFactorVappars<Evaluation<Scalar>>{}(pvt, pt);
    WetGasZeros{}(pvt, pt);
    GasReferenceDensity{}(pvt, pt);
    OilReferenceDensity{}(pvt, pt);
    VapPar1{}(pvt, pt);
    NumRegions{}(pvt, pt);
    InverseFormationVolumeFactorAndViscosity{}(pvt, pt);
};
static_assert(CallableOnConstView<double>);
static_assert(CallableOnConstView<float>);

template <class Functor, class Pvt, class Point, class Result>
__global__ void evaluateKernel(Functor f, Pvt pvt, const Point* points, Result* results, std::size_t n)
{
    const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        results[i] = f(pvt, points[i]);
    }
}

//! As evaluateKernel(), but with the PVT object in device memory, accessed through a
//! pointer to const (as when it is a member of a fluid system in device memory).
template <class Functor, class Pvt, class Point, class Result>
__global__ void evaluateThroughConstPointerKernel(Functor f, const Pvt* pvt, const Point* points,
                                                  Result* results, std::size_t n)
{
    const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        results[i] = f(*pvt, points[i]);
    }
}

constexpr unsigned threadsPerBlock = 128;

unsigned numBlocks(std::size_t n)
{ return static_cast<unsigned>((n + threadsPerBlock - 1) / threadsPerBlock); }

template <class Functor, class Pvt, class Point>
auto evaluateOnCpu(Functor f, const Pvt& pvt, const std::vector<Point>& points)
{
    using Result = decltype(f(pvt, points.front()));
    std::vector<Result> results;
    results.reserve(points.size());
    for (const auto& pt : points) {
        results.push_back(f(pvt, pt));
    }
    return results;
}

template <class Functor, class Pvt, class Point>
auto evaluateOnGpu(Functor f, Pvt pvt, const std::vector<Point>& points)
{
    using Result = decltype(f(pvt, points.front()));
    const Opm::gpuistl::GpuBuffer<Point> gpuPoints(points);
    Opm::gpuistl::GpuBuffer<Result> gpuResults(points.size());

    evaluateKernel<<<numBlocks(points.size()), threadsPerBlock>>>(f, pvt, gpuPoints.data(),
                                                                  gpuResults.data(), points.size());
    OPM_GPU_SAFE_CALL(cudaPeekAtLastError());
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());

    return gpuResults.asStdVector();
}

template <class Functor, class Pvt, class Point>
auto evaluateOnGpuThroughConstPointer(Functor f, Pvt pvt, const std::vector<Point>& points)
{
    using Result = decltype(f(pvt, points.front()));
    const Opm::gpuistl::GpuBuffer<Pvt> gpuPvt(std::vector<Pvt>{pvt});
    const Opm::gpuistl::GpuBuffer<Point> gpuPoints(points);
    Opm::gpuistl::GpuBuffer<Result> gpuResults(points.size());

    evaluateThroughConstPointerKernel<<<numBlocks(points.size()), threadsPerBlock>>>(
        f, gpuPvt.data(), gpuPoints.data(), gpuResults.data(), points.size());
    OPM_GPU_SAFE_CALL(cudaPeekAtLastError());
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());

    return gpuResults.asStdVector();
}

template <std::integral T>
bool isClose(T a, T b)
{ return a == b; }

template <std::floating_point Scalar>
bool isClose(Scalar a, Scalar b, Scalar tol = relTol<Scalar>)
{
    if (a == b) {
        return true;
    }
    return std::abs(a - b) <= tol * std::max(std::abs(a), std::abs(b));
}

template <class Scalar>
bool isClose(const Evaluation<Scalar>& a, const Evaluation<Scalar>& b)
{
    bool close = isClose(a.value(), b.value());
    for (int varIdx = 0; varIdx < Evaluation<Scalar>::numVars; ++varIdx) {
        close = close && isClose(a.derivative(varIdx), b.derivative(varIdx), derivativeRelTol<Scalar>);
    }
    return close;
}

template <class Scalar>
bool isClose(const InvBAndMu<Scalar>& a, const InvBAndMu<Scalar>& b)
{
    return isClose(a.invB, b.invB) && isClose(a.mu, b.mu);
}

template <class Scalar>
std::ostream& operator<<(std::ostream& os, const InvBAndMu<Scalar>& result)
{
    return os << "(" << result.invB << ", " << result.mu << ")";
}

//! Checks that each GPU result is close to its CPU result. Result is a scalar, an
//! Evaluation or InvBAndMu.
template <class Result, class Scalar>
void checkResults(const std::vector<Result>& cpu,
                  const std::vector<Result>& gpu,
                  const std::vector<TestPoint<Scalar>>& points,
                  const std::string& name)
{
    BOOST_REQUIRE_EQUAL(cpu.size(), gpu.size());
    for (std::size_t i = 0; i < cpu.size(); ++i) {
        BOOST_CHECK_MESSAGE(isClose(cpu[i], gpu[i]),
                            name << " at point " << i << " (region " << points[i].regionIdx
                                 << ", p = " << points[i].p << ", Rv = " << points[i].rv
                                 << ", So = " << points[i].so << ", SoMax = " << points[i].soMax
                                 << "): CPU = " << cpu[i] << ", GPU = " << gpu[i]);
    }
}

template <class Functor, class Scalar>
void compareCpuAndGpu(Functor f,
                      const CpuPvt<Scalar>& cpuPvt,
                      GpuViewPvt<Scalar> gpuPvt,
                      const std::vector<TestPoint<Scalar>>& points,
                      const std::string& name)
{
    checkResults(evaluateOnCpu(f, cpuPvt, points), evaluateOnGpu(f, gpuPvt, points), points, name);
}

template <class Scalar>
Scalar rvSat(const CpuPvt<Scalar>& cpuPvt, unsigned regionIdx, Scalar p)
{ return cpuPvt.saturatedOilVaporizationFactor(regionIdx, temperature<Scalar>(), p); }

//! Points inside the tables, between the pressure rows and outside the tables
//! (extrapolation in both directions).
template <class Scalar>
std::vector<TestPoint<Scalar>> gridPoints(const CpuPvt<Scalar>& cpuPvt)
{
    constexpr Scalar bar = 1e5;
    // 10 and 20 bar extrapolate below both regions. The saturated Rv extrapolated from
    // the first two rows is negative at 10 bar in both regions and at 20 bar in region 1,
    // so the RightExtreme guide falls back to vertical interpolation there; at 20 bar in
    // region 0 it is still positive. 260 bar extrapolates above both. The others hit
    // every pressure segment of both regions, and none is a pressure row (nodePoints()
    // covers those).
    const std::vector<Scalar> pressures = {10 * bar, 20 * bar, 45 * bar, 55 * bar, 75 * bar,
                                           95 * bar, 110 * bar, 135 * bar, 165 * bar,
                                           195 * bar, 225 * bar, 260 * bar};
    // From Rv = 0 to above the saturated Rv of every row (extrapolation in Rv). None of
    // them is a PVTG Rv.
    const std::vector<Scalar> rvValues = {0.0, 0.00004, 0.00011, 0.00027, 0.00045, 0.0007, 0.001};

    std::vector<TestPoint<Scalar>> points;
    for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
        for (const Scalar p : pressures) {
            const Scalar sat = rvSat(cpuPvt, regionIdx, p);
            for (const Scalar rv : rvValues) {
                // The branch in inverseFormationVolumeFactorAndViscosity() must not be
                // decided by the last bits of RvSat; switchPoints() tests the switch with
                // suitable margins.
                if (std::abs(rv - sat) <= Scalar(1e-3) * std::abs(sat)) {
                    continue;
                }
                for (const Scalar so : {Scalar{0.0}, Scalar{0.3}}) {
                    points.push_back({regionIdx, p, rv, so, 1});
                }
            }
        }
    }
    return points;
}

//! Points at the first, a middle and the last pressure row of each region's 2D tables,
//! with Rv at each Rv node of that row and halfway between them. At a pressure row, the
//! guided interpolation looks up exactly the given Rv in that row, so the Rv nodes are
//! hit exactly. The value is continuous at a node, but the derivative w.r.t. Rv depends
//! on the segment that is chosen, so this checks that the CPU and the GPU choose the
//! same one.
template <class Scalar>
std::vector<TestPoint<Scalar>> nodePoints(const CpuPvt<Scalar>& cpuPvt)
{
    std::vector<TestPoint<Scalar>> points;
    for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
        const auto& table = cpuPvt.inverseGasB()[regionIdx];
        const std::size_t lastIdx = table.numX() - 1;
        for (const std::size_t xIdx : {std::size_t{0}, lastIdx / 2, lastIdx}) {
            const Scalar p = table.xAt(xIdx);
            const Scalar sat = rvSat(cpuPvt, regionIdx, p);
            std::vector<Scalar> rvs;
            for (std::size_t yIdx = 0; yIdx < table.numY(xIdx); ++yIdx) {
                rvs.push_back(table.yAt(xIdx, yIdx));
                if (yIdx + 1 < table.numY(xIdx)) {
                    rvs.push_back((table.yAt(xIdx, yIdx) + table.yAt(xIdx, yIdx + 1)) / 2);
                }
            }
            for (const Scalar rv : rvs) {
                points.push_back({regionIdx, p, rv, 0, 1});
                // The last node is the saturated Rv. With So > 0 it would take the
                // saturated branch, which switchPoints() covers.
                if (std::abs(rv - sat) > Scalar(1e-3) * sat) {
                    points.push_back({regionIdx, p, rv, 0.3, 1});
                }
            }
        }
    }
    return points;
}

//! Points at and just below the saturated curve, where
//! inverseFormationVolumeFactorAndViscosity() switches between its saturated
//! (So > 0 and Rv >= (1 - 1e-10) * RvSat) and undersaturated branch.
template <class Scalar>
std::vector<TestPoint<Scalar>> switchPoints(const CpuPvt<Scalar>& cpuPvt)
{
    // Rv relative to RvSat. In float, (1 - 1e-10) rounds to 1, so the switch is at
    // Rv >= RvSat. A point exactly on the curve could then take different branches on
    // the CPU and the GPU if RvSat differs in the last bit, so for float the points
    // keep a safe distance from the curve.
    const std::vector<Scalar> rvFactors = std::is_same_v<Scalar, float>
        ? std::vector<Scalar>{1.0 + 1e-4, 1.0 - 1e-4, 0.5}
        : std::vector<Scalar>{1.0, 1.0 - 1e-11, 1.0 - 1e-9, 0.5};

    constexpr Scalar bar = 1e5;
    std::vector<TestPoint<Scalar>> points;
    for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
        for (const Scalar p : {55 * bar, 110 * bar, 195 * bar}) {
            const Scalar sat = rvSat(cpuPvt, regionIdx, p);
            for (const Scalar factor : rvFactors) {
                // So > 0: saturated branch if Rv is on the curve. So = 0: always the
                // undersaturated branch, also with Rv on the curve.
                points.push_back({regionIdx, p, factor * sat, 0.3, 1});
                points.push_back({regionIdx, p, factor * sat, 0, 1});
            }
        }
    }
    return points;
}

template <class Scalar>
struct Fixture {
    Fixture()
        : cpuPvt(makeCpuPvt<Scalar>())
        , gpuBufPvt(Opm::gpuistl::copy_to_gpu(cpuPvt))
        , gpuViewPvt(Opm::gpuistl::make_view(gpuBufPvt))
    {}

    CpuPvt<Scalar> cpuPvt;
    GpuBufPvt<Scalar> gpuBufPvt;
    GpuViewPvt<Scalar> gpuViewPvt;
};

} // anonymous namespace

BOOST_AUTO_TEST_CASE_TEMPLATE(SizesAfterCopyToGpu, Scalar, Scalars)
{
    const Fixture<Scalar> fx;

    BOOST_REQUIRE_EQUAL(fx.cpuPvt.numRegions(), numRegions);
    BOOST_CHECK_EQUAL(fx.gpuBufPvt.numRegions(), numRegions);
    BOOST_CHECK_EQUAL(fx.gpuViewPvt.numRegions(), numRegions);

    // Only sizes are safe to read on the host; the elements live in device memory.
    BOOST_CHECK_EQUAL(fx.gpuBufPvt.inverseGasB().size(), numRegions);
    BOOST_CHECK_EQUAL(fx.gpuBufPvt.inverseSaturatedGasB().size(), numRegions);
    BOOST_CHECK_EQUAL(fx.gpuBufPvt.gasMu().size(), numRegions);
    BOOST_CHECK_EQUAL(fx.gpuBufPvt.inverseGasBMu().size(), numRegions);
    BOOST_CHECK_EQUAL(fx.gpuBufPvt.inverseSaturatedGasBMu().size(), numRegions);
    BOOST_CHECK_EQUAL(fx.gpuBufPvt.saturatedOilVaporizationFactorTable().size(), numRegions);
    BOOST_CHECK_EQUAL(fx.gpuBufPvt.saturationPressure().size(), numRegions);

    BOOST_CHECK_EQUAL(fx.gpuViewPvt.inverseGasB().size(), numRegions);
    BOOST_CHECK_EQUAL(fx.gpuViewPvt.saturationPressure().size(), numRegions);

    // The table layout the tests below rely on (see the deck).
    BOOST_CHECK_EQUAL(fx.cpuPvt.inverseGasB()[0].numX(), 3u);
    BOOST_CHECK_EQUAL(fx.cpuPvt.inverseGasB()[1].numX(), 8u);
    BOOST_CHECK_EQUAL(fx.cpuPvt.inverseGasB()[1].numY(3), 8u);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(CopyToGpuRejectsUninitializedPvt, Scalar, Scalars)
{
    // No regions
    CpuPvt<Scalar> uninitialized;
    BOOST_CHECK_THROW(Opm::gpuistl::copy_to_gpu(uninitialized), std::logic_error);

    // Regions, but initEnd() not called, so the tables it fills are empty
    CpuPvt<Scalar> withoutInitEnd;
    withoutInitEnd.setNumRegions(numRegions);
    BOOST_CHECK_THROW(Opm::gpuistl::copy_to_gpu(withoutInitEnd), std::logic_error);
}

//! Compares the methods that depend on the state (pressure, Rv), with the state passed
//! as Eval.
template <class Eval, class Scalar>
void compareStateDependentMethods(const Fixture<Scalar>& fx, const std::vector<TestPoint<Scalar>>& points)
{
    compareCpuAndGpu(InverseFormationVolumeFactor<Eval>{}, fx.cpuPvt, fx.gpuViewPvt, points, "inverseFormationVolumeFactor");
    compareCpuAndGpu(Viscosity<Eval>{}, fx.cpuPvt, fx.gpuViewPvt, points, "viscosity");
    compareCpuAndGpu(SaturatedInverseFormationVolumeFactor<Eval>{}, fx.cpuPvt, fx.gpuViewPvt, points, "saturatedInverseFormationVolumeFactor");
    compareCpuAndGpu(SaturatedViscosity<Eval>{}, fx.cpuPvt, fx.gpuViewPvt, points, "saturatedViscosity");
    compareCpuAndGpu(SaturatedOilVaporizationFactor<Eval>{}, fx.cpuPvt, fx.gpuViewPvt, points, "saturatedOilVaporizationFactor");
}

BOOST_AUTO_TEST_CASE_TEMPLATE(ScalarMethodsMatchCpu, Scalar, Scalars)
{
    const Fixture<Scalar> fx;
    const auto points = gridPoints(fx.cpuPvt);
    compareStateDependentMethods<Scalar>(fx, points);
    compareCpuAndGpu(GasReferenceDensity{}, fx.cpuPvt, fx.gpuViewPvt, points, "gasReferenceDensity");
    compareCpuAndGpu(OilReferenceDensity{}, fx.cpuPvt, fx.gpuViewPvt, points, "oilReferenceDensity");
    compareCpuAndGpu(VapPar1{}, fx.cpuPvt, fx.gpuViewPvt, points, "vapPar1");
    compareCpuAndGpu(NumRegions{}, fx.cpuPvt, fx.gpuViewPvt, points, "numRegions");
}

BOOST_AUTO_TEST_CASE_TEMPLATE(EvaluationMethodsMatchCpu, Scalar, Scalars)
{
    // The fluid system calls these methods with AD values, which go through other
    // table functions than inverseFormationVolumeFactorAndViscosity(). Values and
    // derivatives are compared.
    const Fixture<Scalar> fx;
    compareStateDependentMethods<Evaluation<Scalar>>(fx, gridPoints(fx.cpuPvt));
    compareStateDependentMethods<Evaluation<Scalar>>(fx, nodePoints(fx.cpuPvt));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(RightExtremeInterpolationMatchesCpu, Scalar, Scalars)
{
    // With the RightExtreme guide, the shift between the two pressure rows scales with
    // Rv / RvSat(p), so it depends on Rv (and on p), and so do the derivatives. This
    // checks the parts of UniformXTabulated2DFunction::findPoints() that only this guide
    // uses, which LiveOilPvt (LeftExtreme) does not reach.
    const Fixture<Scalar> fx;
    for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
        BOOST_REQUIRE(fx.cpuPvt.inverseGasB()[regionIdx].interpolationGuide()
                      == Opm::UniformXTabulated2DFunctionInterpolationPolicy::RightExtreme);
        BOOST_REQUIRE(fx.cpuPvt.inverseGasBMu()[regionIdx].interpolationGuide()
                      == Opm::UniformXTabulated2DFunctionInterpolationPolicy::RightExtreme);
    }

    constexpr Scalar bar = 1e5;
    std::vector<TestPoint<Scalar>> points;
    for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
        for (const Scalar p : {10 * bar, 20 * bar, 75 * bar, 135 * bar, 195 * bar, 260 * bar}) {
            const Scalar sat = rvSat(fx.cpuPvt, regionIdx, p);
            // Rv = 0 (no shift), Rv between 0 and RvSat, on the saturated curve, and
            // above it (extrapolation in Rv).
            for (const Scalar rv : {Scalar{0.0}, Scalar(0.3) * sat, Scalar(0.8) * sat, sat,
                                    Scalar(1.5) * sat, Scalar{0.001}}) {
                if (rv >= 0) {
                    points.push_back({regionIdx, p, rv, 0, 1});
                }
            }
        }
    }

    // The saturated Rv extrapolated below the first rows is negative at 10 bar in both
    // regions, so the yEnd <= 0 branch of findPoints() is taken there.
    for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
        BOOST_REQUIRE_LT(rvSat(fx.cpuPvt, regionIdx, Scalar(10 * bar)), 0);
    }

    // Inside the tables, the guide follows the saturated curve: at Rv = RvSat(p), the 2D
    // table equals the saturated 1D table, also between the pressure rows.
    for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
        for (const Scalar p : {75 * bar, 135 * bar}) {
            const Scalar sat = rvSat(fx.cpuPvt, regionIdx, p);
            const TestPoint<Scalar> pt{regionIdx, p, sat, 0, 1};
            const Scalar invB = InverseFormationVolumeFactor<Scalar>{}(fx.cpuPvt, pt);
            const Scalar satInvB = SaturatedInverseFormationVolumeFactor<Scalar>{}(fx.cpuPvt, pt);
            BOOST_CHECK_CLOSE(invB, satInvB, 100 * relTol<Scalar>);
        }
    }

    // Between two pressure rows, the derivative w.r.t. Rv is non-zero for Rv > 0.
    const auto dInvB = InverseFormationVolumeFactor<Evaluation<Scalar>>{}(
        fx.cpuPvt, TestPoint<Scalar>{1, 135 * bar, 0.0002, 0, 1});
    BOOST_CHECK_NE(dInvB.derivative(1), 0);

    compareStateDependentMethods<Scalar>(fx, points);
    compareStateDependentMethods<Evaluation<Scalar>>(fx, points);
    compareCpuAndGpu(InverseFormationVolumeFactorAndViscosity{}, fx.cpuPvt, fx.gpuViewPvt,
                     points, "inverseFormationVolumeFactorAndViscosity (RightExtreme)");
}

BOOST_AUTO_TEST_CASE_TEMPLATE(InverseFormationVolumeFactorAndViscosityMatchesCpu, Scalar, Scalars)
{
    const Fixture<Scalar> fx;
    compareCpuAndGpu(InverseFormationVolumeFactorAndViscosity{}, fx.cpuPvt, fx.gpuViewPvt,
                     gridPoints(fx.cpuPvt), "inverseFormationVolumeFactorAndViscosity (grid)");
    compareCpuAndGpu(InverseFormationVolumeFactorAndViscosity{}, fx.cpuPvt, fx.gpuViewPvt,
                     nodePoints(fx.cpuPvt), "inverseFormationVolumeFactorAndViscosity (nodes)");

    const auto points = switchPoints(fx.cpuPvt);

    // Make sure that both branches are exercised: In the saturated branch, the result
    // equals the saturated tables. In the undersaturated branch, it comes from the 2D
    // tables and differs from the saturated tables away from the saturated curve.
    const auto cpuResults = evaluateOnCpu(InverseFormationVolumeFactorAndViscosity{}, fx.cpuPvt, points);
    for (std::size_t i = 0; i < points.size(); ++i) {
        const auto& pt = points[i];
        const Scalar satInvB = fx.cpuPvt.saturatedInverseFormationVolumeFactor(pt.regionIdx, temperature<Scalar>(), pt.p);
        const Scalar sat = rvSat(fx.cpuPvt, pt.regionIdx, pt.p);
        if (pt.so > 0 && pt.rv >= Scalar(1.0 - 1e-10) * sat) {
            BOOST_CHECK_EQUAL(cpuResults[i].invB.value(), satInvB);
            // The saturated tables do not depend on Rv.
            BOOST_CHECK_EQUAL(cpuResults[i].invB.derivative(1), 0);
        } else if (pt.rv < Scalar(0.9) * sat) {
            BOOST_CHECK_NE(cpuResults[i].invB.value(), satInvB);
        }
        if (pt.so == 0) {
            // With So = 0 the undersaturated branch is taken even on the saturated curve,
            // so the result depends on Rv.
            BOOST_CHECK_NE(cpuResults[i].invB.derivative(1), 0);
        }
    }

    compareCpuAndGpu(InverseFormationVolumeFactorAndViscosity{}, fx.cpuPvt, fx.gpuViewPvt,
                     points, "inverseFormationVolumeFactorAndViscosity (branch switch)");
}

BOOST_AUTO_TEST_CASE_TEMPLATE(SaturatedOilVaporizationFactorWithVappars, Scalar, Scalars)
{
    auto cpuPvt = makeCpuPvt<Scalar>();
    cpuPvt.setVapPars(2.0, 0.1);
    auto gpuBufPvt = Opm::gpuistl::copy_to_gpu(cpuPvt);
    const auto gpuViewPvt = Opm::gpuistl::make_view(gpuBufPvt);

    constexpr Scalar p = 120e5;
    std::vector<TestPoint<Scalar>> points;
    for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
        // So < SoMax: the VAPPARS correction pow(So/SoMax, vapPar1) is applied
        points.push_back({regionIdx, p, 0, 0.3, 0.8});
        points.push_back({regionIdx, p, 0, 0.0005, 0.8}); // So below the cut-off eps
        points.push_back({regionIdx, p, 0, 0.5, 1.2}); // SoMax > 1 is clamped to 1
        // So >= SoMax, and SoMax <= 0.01: no correction
        points.push_back({regionIdx, p, 0, 0.8, 0.8});
        points.push_back({regionIdx, p, 0, 0.005, 0.009});
    }

    // The correction is active for the first point of each region, and with AD values it
    // depends on So, so the derivatives of pow() and max() are compared as well.
    const Scalar uncorrected = rvSat(cpuPvt, 0, p);
    BOOST_CHECK_LT(SaturatedOilVaporizationFactorVappars<Scalar>{}(cpuPvt, points[0]), uncorrected);
    BOOST_CHECK_NE(SaturatedOilVaporizationFactorVappars<Evaluation<Scalar>>{}(cpuPvt, points[0]).derivative(1), 0);

    compareCpuAndGpu(SaturatedOilVaporizationFactorVappars<Scalar>{}, cpuPvt, gpuViewPvt, points,
                     "saturatedOilVaporizationFactor (VAPPARS)");
    compareCpuAndGpu(SaturatedOilVaporizationFactorVappars<Evaluation<Scalar>>{}, cpuPvt, gpuViewPvt, points,
                     "saturatedOilVaporizationFactor (VAPPARS, AD)");
    compareCpuAndGpu(VapPar1{}, cpuPvt, gpuViewPvt, points, "vapPar1 (VAPPARS)");
}

BOOST_AUTO_TEST_CASE_TEMPLATE(WetGasMethodsReturnZeroOnGpu, Scalar, Scalars)
{
    // No vaporized water, also on the GPU.
    const Fixture<Scalar> fx;
    for (const auto& r : evaluateOnGpu(WetGasZeros{}, fx.gpuViewPvt, gridPoints(fx.cpuPvt))) {
        BOOST_CHECK(r == Evaluation<Scalar>(0.0));
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(RegionsAreIndexedCorrectly, Scalar, Scalars)
{
    const Fixture<Scalar> fx;

    // The same (p, Rv) in both regions must give different values (the regions have
    // different tables), and each region must match its CPU counterpart.
    constexpr Scalar p = 130e5;
    constexpr Scalar rv = 0.0002;
    const std::vector<TestPoint<Scalar>> points = {{0, p, rv, 0, 1},
                                                   {1, p, rv, 0, 1}};

    const auto cpuInvB = evaluateOnCpu(InverseFormationVolumeFactor<Scalar>{}, fx.cpuPvt, points);
    const auto cpuMu = evaluateOnCpu(Viscosity<Scalar>{}, fx.cpuPvt, points);
    BOOST_CHECK_NE(cpuInvB[0], cpuInvB[1]);
    BOOST_CHECK_NE(cpuMu[0], cpuMu[1]);

    checkResults(cpuInvB, evaluateOnGpu(InverseFormationVolumeFactor<Scalar>{}, fx.gpuViewPvt, points), points,
                 "inverseFormationVolumeFactor (regions)");
    checkResults(cpuMu, evaluateOnGpu(Viscosity<Scalar>{}, fx.gpuViewPvt, points), points,
                 "viscosity (regions)");
    compareCpuAndGpu(GasReferenceDensity{}, fx.cpuPvt, fx.gpuViewPvt, points, "gasReferenceDensity (regions)");
    compareCpuAndGpu(OilReferenceDensity{}, fx.cpuPvt, fx.gpuViewPvt, points, "oilReferenceDensity (regions)");
}

BOOST_AUTO_TEST_CASE_TEMPLATE(CallableThroughConstView, Scalar, Scalars)
{
    // That every method compiles on a const view is checked by CallableOnConstView
    // above. Here the view is placed in device memory and used through a const pointer.
    const Fixture<Scalar> fx;
    const auto points = gridPoints(fx.cpuPvt);
    checkResults(evaluateOnCpu(InverseFormationVolumeFactorAndViscosity{}, fx.cpuPvt, points),
                 evaluateOnGpuThroughConstPointer(InverseFormationVolumeFactorAndViscosity{},
                                                  fx.gpuViewPvt, points),
                 points, "inverseFormationVolumeFactorAndViscosity (const pointer)");
}

BOOST_AUTO_TEST_CASE_TEMPLATE(ViewsStayValidAfterMovingOwner, Scalar, Scalars)
{
    Fixture<Scalar> fx;

    // Moving the owning GPU object moves the pointer to its device table owner and the
    // outer GpuBuffers of views; the tables and their device memory stay in place. Both a
    // view created before the move and one created after it must stay valid.
    const auto viewBeforeMove = fx.gpuViewPvt;
    GpuBufPvt<Scalar> movedPvt = std::move(fx.gpuBufPvt);
    const auto viewAfterMove = Opm::gpuistl::make_view(movedPvt);

    // Move-assigning into an existing owning GPU object frees the tables it owned, but
    // must leave the moved-in tables, and thereby the views above, untouched.
    GpuBufPvt<Scalar> assignedPvt = Opm::gpuistl::copy_to_gpu(fx.cpuPvt);
    assignedPvt = std::move(movedPvt);
    const auto viewAfterAssignment = Opm::gpuistl::make_view(assignedPvt);

    const auto points = gridPoints(fx.cpuPvt);
    compareCpuAndGpu(InverseFormationVolumeFactorAndViscosity{}, fx.cpuPvt, viewBeforeMove, points,
                     "inverseFormationVolumeFactorAndViscosity (view before move)");
    compareCpuAndGpu(InverseFormationVolumeFactorAndViscosity{}, fx.cpuPvt, viewAfterMove, points,
                     "inverseFormationVolumeFactorAndViscosity (view after move)");
    compareCpuAndGpu(InverseFormationVolumeFactorAndViscosity{}, fx.cpuPvt, viewAfterAssignment, points,
                     "inverseFormationVolumeFactorAndViscosity (view after move assignment)");
}
