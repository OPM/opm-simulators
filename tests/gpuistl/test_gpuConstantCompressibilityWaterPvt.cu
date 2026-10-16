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
 * \brief Unit tests for ConstantCompressibilityWaterPvt on the GPU.
 *
 * A CPU ConstantCompressibilityWaterPvt is initialized from a deck with two PVT
 * regions, copied to the GPU with copy_to_gpu()/make_view(), and evaluated in
 * kernels. All GPU results are compared to the same evaluations on the CPU.
 */
#include <config.h>

#define BOOST_TEST_MODULE TestGpuConstantCompressibilityWaterPvt
#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/fluidsystems/blackoilpvt/ConstantCompressibilityWaterPvt.hpp>

#include <opm/common/utility/gpuDecorators.hpp>
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

// One derivative, with respect to pressure.
template <class Scalar>
using Evaluation = Opm::DenseAd::Evaluation<Scalar, 1>;

template <class Scalar>
using CpuPvt = Opm::ConstantCompressibilityWaterPvt<Scalar>;
template <class Scalar>
using GpuBufPvt = Opm::ConstantCompressibilityWaterPvt<Scalar, Opm::gpuistl::GpuBuffer>;
template <class Scalar>
using GpuViewPvt = Opm::ConstantCompressibilityWaterPvt<Scalar, Opm::gpuistl::GpuView>;

// ConstantCompressibilityWaterPvt<float> is explicitly instantiated in opm-common
// independently of FLOW_INSTANTIATE_FLOAT, so both types are always tested.
using Scalars = boost::mpl::list<double, float>;

// The view is passed by value to kernels.
static_assert(std::is_trivially_copyable_v<GpuViewPvt<double>>);
static_assert(std::is_trivially_copyable_v<GpuViewPvt<float>>);

// Two PVT regions with different parameters. Region 0 has a non-zero viscosibility,
// so that the viscosity has its own pressure dependence. Region 1 has zero
// viscosibility.
constexpr const char* deckString = R"(
RUNSPEC
DIMENS
 1 1 1 /
TABDIMS
 * 2 /
OIL
WATER
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
PVTW
-- PREF   BW    CW       MUW   CV
   200.0  1.02  4.5e-5   0.5   1.0e-5 /
   250.0  1.04  6.0e-5   0.8   0.0    /
)";

constexpr unsigned numRegions = 2;
// [K], not used by ConstantCompressibilityWaterPvt.
template <class Scalar>
OPM_HOST_DEVICE constexpr Scalar temperature()
{ return 350.0; }

// Relative tolerance for comparing CPU and GPU results. The GPU may contract
// multiplications and additions into FMAs, so results can differ in the last bits.
template <class Scalar>
constexpr Scalar relTol = std::is_same_v<Scalar, float> ? 1e-4 : 1e-10;

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

//! One evaluation point. Pressure in [Pa].
template <class Scalar>
struct TestPoint {
    unsigned regionIdx;
    Scalar p;
};

//! The parts of a fluid state used by
//! ConstantCompressibilityWaterPvt::inverseFormationVolumeFactorAndViscosity().
template <class Scalar>
struct TestFluidState {
    using ValueType = Evaluation<Scalar>;
    static constexpr int waterPhaseIdx = 0;

    OPM_HOST_DEVICE const ValueType& pressure(unsigned) const
    { return p; }

    ValueType p;
};

//! Result of the methods returning both 1/B and mu. A plain struct instead of
//! std::pair, because it is copied between host and device with memcpy.
template <class Scalar>
struct InvBAndMu {
    Evaluation<Scalar> invB;
    Evaluation<Scalar> mu;
};

// Functors evaluating one ConstantCompressibilityWaterPvt method. They are templated on
// the PVT type, so the same functor is used for the CPU reference and in the kernel.
// All arguments except the pressure are constants; the pressure is the AD variable.

template <class Scalar>
OPM_HOST_DEVICE Evaluation<Scalar> pressureVariable(const TestPoint<Scalar>& pt)
{ return Evaluation<Scalar>::createVariable(pt.p, 0); }

struct InverseFormationVolumeFactor {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE Evaluation<Scalar> operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    {
        using Eval = Evaluation<Scalar>;
        return pvt.inverseFormationVolumeFactor(pt.regionIdx, Eval(temperature<Scalar>()),
                                                pressureVariable(pt), Eval(0.0), Eval(0.0));
    }
};

struct Viscosity {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE Evaluation<Scalar> operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    {
        using Eval = Evaluation<Scalar>;
        return pvt.viscosity(pt.regionIdx, Eval(temperature<Scalar>()),
                             pressureVariable(pt), Eval(0.0), Eval(0.0));
    }
};

struct SaturatedInverseFormationVolumeFactor {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE Evaluation<Scalar> operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    {
        using Eval = Evaluation<Scalar>;
        return pvt.saturatedInverseFormationVolumeFactor(pt.regionIdx, Eval(temperature<Scalar>()),
                                                         pressureVariable(pt), Eval(0.0));
    }
};

struct SaturatedViscosity {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE Evaluation<Scalar> operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    {
        using Eval = Evaluation<Scalar>;
        return pvt.saturatedViscosity(pt.regionIdx, Eval(temperature<Scalar>()),
                                      pressureVariable(pt), Eval(0.0));
    }
};

struct SaturationPressure {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE Evaluation<Scalar> operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    {
        using Eval = Evaluation<Scalar>;
        return pvt.saturationPressure(pt.regionIdx, Eval(temperature<Scalar>()), Eval(0.0), Eval(0.0));
    }
};

struct SaturatedGasDissolutionFactor {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE Evaluation<Scalar> operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    {
        using Eval = Evaluation<Scalar>;
        return pvt.saturatedGasDissolutionFactor(pt.regionIdx, Eval(temperature<Scalar>()),
                                                 pressureVariable(pt), Eval(0.0));
    }
};

struct WaterReferenceDensity {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE Scalar operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    { return pvt.waterReferenceDensity(pt.regionIdx); }
};

struct NumRegions {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE unsigned operator()(const Pvt& pvt, const TestPoint<Scalar>&) const
    { return pvt.numRegions(); }
};

//! inverseBAndMu() with all arguments.
struct InverseBAndMu {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE InvBAndMu<Scalar> operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    {
        using Eval = Evaluation<Scalar>;
        Eval invB;
        Eval mu;
        pvt.inverseBAndMu(invB, mu, pt.regionIdx, Eval(temperature<Scalar>()),
                          pressureVariable(pt), Eval(0.0), Eval(0.0));
        return {invB, mu};
    }
};

//! inverseBAndMu() with the pressure only.
struct InverseBAndMuPressureOnly {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE InvBAndMu<Scalar> operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    {
        Evaluation<Scalar> invB;
        Evaluation<Scalar> mu;
        pvt.inverseBAndMu(invB, mu, pt.regionIdx, pressureVariable(pt));
        return {invB, mu};
    }
};

//! The hot path, called by the fluid system.
struct InverseFormationVolumeFactorAndViscosity {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE InvBAndMu<Scalar> operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    {
        TestFluidState<Scalar> fs;
        fs.p = pressureVariable(pt);
        const auto [invB, mu] = pvt.inverseFormationVolumeFactorAndViscosity(fs, pt.regionIdx);
        return {invB, mu};
    }
};

// Every method used on the GPU must be callable on a const view, since the fluid
// system holding the PVT object may only be available through a const reference.
template <class Scalar>
concept CallableOnConstView = requires(const GpuViewPvt<Scalar>& pvt,
                                       const TestFluidState<Scalar>& fs,
                                       const TestPoint<Scalar>& pt) {
    InverseFormationVolumeFactor{}(pvt, pt);
    Viscosity{}(pvt, pt);
    SaturatedInverseFormationVolumeFactor{}(pvt, pt);
    SaturatedViscosity{}(pvt, pt);
    SaturationPressure{}(pvt, pt);
    SaturatedGasDissolutionFactor{}(pvt, pt);
    WaterReferenceDensity{}(pvt, pt);
    NumRegions{}(pvt, pt);
    InverseBAndMu{}(pvt, pt);
    InverseBAndMuPressureOnly{}(pvt, pt);
    pvt.inverseFormationVolumeFactorAndViscosity(fs, pt.regionIdx);
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
bool isClose(Scalar a, Scalar b)
{
    if (a == b) {
        return true;
    }
    return std::abs(a - b) <= relTol<Scalar> * std::max(std::abs(a), std::abs(b));
}

// [1/Pa] Upper bound of the compressibilities and viscosibilities in the deck, so the
// pressure derivatives are at most of order value * pressureDerivativeScale.
template <class Scalar>
constexpr Scalar pressureDerivativeScale = 1e-9;

template <class Scalar>
bool isClose(const Evaluation<Scalar>& a, const Evaluation<Scalar>& b)
{
    bool close = isClose(a.value(), b.value());
    // A derivative may be zero up to rounding: With zero viscosibility (region 1) the
    // viscosity does not depend on pressure, and its derivative is the cancellation of
    // terms of order value * compressibility. A purely relative tolerance would then
    // compare rounding noise, so the tolerance has an absolute floor relative to that
    // scale.
    const Scalar derivativeFloor = std::abs(a.value()) * pressureDerivativeScale<Scalar>;
    for (int varIdx = 0; varIdx < Evaluation<Scalar>::numVars; ++varIdx) {
        const Scalar da = a.derivative(varIdx);
        const Scalar db = b.derivative(varIdx);
        close = close && std::abs(da - db) <= relTol<Scalar> * std::max({std::abs(da), std::abs(db), derivativeFloor});
    }
    return close;
}

template <class Scalar>
bool isClose(const InvBAndMu<Scalar>& a, const InvBAndMu<Scalar>& b)
{ return isClose(a.invB, b.invB) && isClose(a.mu, b.mu); }

template <class Scalar>
std::ostream& operator<<(std::ostream& os, const InvBAndMu<Scalar>& r)
{ return os << "(" << r.invB << ", " << r.mu << ")"; }

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
                                 << ", p = " << points[i].p << "): CPU = " << cpu[i]
                                 << ", GPU = " << gpu[i]);
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

//! Pressures below, at and above the reference pressures of both regions (200 and
//! 250 bar), so that (p - pRef) takes both signs and zero, plus extreme pressures.
template <class Scalar>
std::vector<TestPoint<Scalar>> pressurePoints()
{
    constexpr Scalar bar = 1e5;
    const std::vector<Scalar> pressures = {1 * bar, 50 * bar, 150 * bar, 200 * bar,
                                           250 * bar, 400 * bar, 1000 * bar};

    std::vector<TestPoint<Scalar>> points;
    for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
        for (const Scalar p : pressures) {
            points.push_back({regionIdx, p});
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

BOOST_AUTO_TEST_CASE_TEMPLATE(SizesAndValuesAfterCopyToGpu, Scalar, Scalars)
{
    const Fixture<Scalar> fx;

    BOOST_REQUIRE_EQUAL(fx.cpuPvt.numRegions(), numRegions);
    BOOST_CHECK_EQUAL(fx.gpuBufPvt.numRegions(), numRegions);
    BOOST_CHECK_EQUAL(fx.gpuViewPvt.numRegions(), numRegions);

    // The parameters are copied, not computed, so they must match exactly.
    auto checkEqual = [](const auto& cpu, const auto& gpu, const std::string& name) {
        const auto gpuOnHost = gpu.asStdVector();
        BOOST_CHECK_MESSAGE(std::ranges::equal(cpu, gpuOnHost), name << " differs after copy_to_gpu()");
    };
    checkEqual(fx.cpuPvt.waterReferencePressure(), fx.gpuBufPvt.waterReferencePressure(),
               "waterReferencePressure");
    checkEqual(fx.cpuPvt.waterReferenceFormationVolumeFactor(),
               fx.gpuBufPvt.waterReferenceFormationVolumeFactor(),
               "waterReferenceFormationVolumeFactor");
    checkEqual(fx.cpuPvt.waterCompressibility(), fx.gpuBufPvt.waterCompressibility(),
               "waterCompressibility");
    checkEqual(fx.cpuPvt.waterViscosity(), fx.gpuBufPvt.waterViscosity(), "waterViscosity");
    checkEqual(fx.cpuPvt.waterViscosibility(), fx.gpuBufPvt.waterViscosibility(),
               "waterViscosibility");

    // There is no accessor for the reference density vector, so it is read in a kernel.
    const auto points = pressurePoints<Scalar>();
    const auto cpuDensity = evaluateOnCpu(WaterReferenceDensity{}, fx.cpuPvt, points);
    const auto gpuDensity = evaluateOnGpu(WaterReferenceDensity{}, fx.gpuViewPvt, points);
    BOOST_CHECK(cpuDensity == gpuDensity);
    compareCpuAndGpu(NumRegions{}, fx.cpuPvt, fx.gpuViewPvt, points, "numRegions");
}

BOOST_AUTO_TEST_CASE_TEMPLATE(CopyToGpuRejectsUninitializedPvt, Scalar, Scalars)
{
    const CpuPvt<Scalar> uninitialized;
    BOOST_CHECK_THROW(Opm::gpuistl::copy_to_gpu(uninitialized), std::logic_error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(MethodsMatchCpu, Scalar, Scalars)
{
    const Fixture<Scalar> fx;
    const auto points = pressurePoints<Scalar>();
    compareCpuAndGpu(InverseFormationVolumeFactor{}, fx.cpuPvt, fx.gpuViewPvt, points, "inverseFormationVolumeFactor");
    compareCpuAndGpu(Viscosity{}, fx.cpuPvt, fx.gpuViewPvt, points, "viscosity");
    compareCpuAndGpu(SaturatedInverseFormationVolumeFactor{}, fx.cpuPvt, fx.gpuViewPvt, points, "saturatedInverseFormationVolumeFactor");
    compareCpuAndGpu(SaturatedViscosity{}, fx.cpuPvt, fx.gpuViewPvt, points, "saturatedViscosity");
    compareCpuAndGpu(InverseBAndMu{}, fx.cpuPvt, fx.gpuViewPvt, points, "inverseBAndMu");
    compareCpuAndGpu(InverseBAndMuPressureOnly{}, fx.cpuPvt, fx.gpuViewPvt, points, "inverseBAndMu (pressure only)");
    compareCpuAndGpu(InverseFormationVolumeFactorAndViscosity{}, fx.cpuPvt, fx.gpuViewPvt, points, "inverseFormationVolumeFactorAndViscosity");
    compareCpuAndGpu(WaterReferenceDensity{}, fx.cpuPvt, fx.gpuViewPvt, points, "waterReferenceDensity");

    // Dead water: no dissolved gas and no saturation pressure, also on the GPU.
    for (const auto& r : evaluateOnGpu(SaturationPressure{}, fx.gpuViewPvt, points)) {
        BOOST_CHECK(r == Evaluation<Scalar>(0.0));
    }
    for (const auto& r : evaluateOnGpu(SaturatedGasDissolutionFactor{}, fx.gpuViewPvt, points)) {
        BOOST_CHECK(r == Evaluation<Scalar>(0.0));
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(RegionsAreIndexedCorrectly, Scalar, Scalars)
{
    const Fixture<Scalar> fx;

    // The same pressure in both regions must give different values (the regions have
    // different parameters), and each region must match its CPU counterpart.
    constexpr Scalar p = 300e5;
    const std::vector<TestPoint<Scalar>> points = {{0, p}, {1, p}};

    const auto cpuInvB = evaluateOnCpu(InverseFormationVolumeFactor{}, fx.cpuPvt, points);
    const auto cpuMu = evaluateOnCpu(Viscosity{}, fx.cpuPvt, points);
    const auto cpuDensity = evaluateOnCpu(WaterReferenceDensity{}, fx.cpuPvt, points);
    BOOST_CHECK_NE(cpuInvB[0].value(), cpuInvB[1].value());
    BOOST_CHECK_NE(cpuInvB[0].derivative(0), cpuInvB[1].derivative(0));
    BOOST_CHECK_NE(cpuMu[0].value(), cpuMu[1].value());
    BOOST_CHECK_NE(cpuMu[0].derivative(0), cpuMu[1].derivative(0));
    BOOST_CHECK_NE(cpuDensity[0], cpuDensity[1]);

    checkResults(cpuInvB, evaluateOnGpu(InverseFormationVolumeFactor{}, fx.gpuViewPvt, points), points,
                 "inverseFormationVolumeFactor (regions)");
    checkResults(cpuMu, evaluateOnGpu(Viscosity{}, fx.gpuViewPvt, points), points,
                 "viscosity (regions)");
    checkResults(cpuDensity, evaluateOnGpu(WaterReferenceDensity{}, fx.gpuViewPvt, points), points,
                 "waterReferenceDensity (regions)");
}

BOOST_AUTO_TEST_CASE_TEMPLATE(CallableThroughConstView, Scalar, Scalars)
{
    // That every method compiles on a const view is checked by CallableOnConstView
    // above. Here the view is placed in device memory and used through a const pointer.
    const Fixture<Scalar> fx;
    const auto points = pressurePoints<Scalar>();
    checkResults(evaluateOnCpu(InverseFormationVolumeFactorAndViscosity{}, fx.cpuPvt, points),
                 evaluateOnGpuThroughConstPointer(InverseFormationVolumeFactorAndViscosity{},
                                                  fx.gpuViewPvt, points),
                 points, "inverseFormationVolumeFactorAndViscosity (const pointer)");
}

BOOST_AUTO_TEST_CASE_TEMPLATE(ViewsStayValidAfterMovingOwner, Scalar, Scalars)
{
    Fixture<Scalar> fx;

    // Moving the owning object moves the GpuBuffers, which keeps their device memory.
    // Both a view created before the move and one created after it must stay valid.
    const auto viewBeforeMove = fx.gpuViewPvt;
    GpuBufPvt<Scalar> movedPvt = std::move(fx.gpuBufPvt);
    const auto viewAfterMove = Opm::gpuistl::make_view(movedPvt);

    const auto points = pressurePoints<Scalar>();
    compareCpuAndGpu(InverseFormationVolumeFactorAndViscosity{}, fx.cpuPvt, viewBeforeMove, points,
                     "inverseFormationVolumeFactorAndViscosity (view before move)");
    compareCpuAndGpu(InverseFormationVolumeFactorAndViscosity{}, fx.cpuPvt, viewAfterMove, points,
                     "inverseFormationVolumeFactorAndViscosity (view after move)");
}
