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
 * \brief Unit tests for DryGasPvt on the GPU.
 *
 * A CPU DryGasPvt is initialized from a deck with two PVT regions of different
 * table sizes, copied to the GPU with copy_to_gpu()/make_view(), and evaluated
 * in kernels. All GPU results are compared to the same evaluations on the CPU.
 */
#include <config.h>

#define BOOST_TEST_MODULE TestGpuDryGasPvt
#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/fluidsystems/blackoilpvt/DryGasPvt.hpp>

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

// One derivative, with respect to pressure. DryGasPvt depends on the pressure only.
template <class Scalar>
using Evaluation = Opm::DenseAd::Evaluation<Scalar, 1>;

template <class Scalar>
using CpuPvt = Opm::DryGasPvt<Scalar>;
template <class Scalar>
using GpuBufPvt = Opm::DryGasPvt<Scalar, Opm::gpuistl::GpuBuffer, Opm::gpuistl::GpuView>;
template <class Scalar>
using GpuViewPvt = Opm::DryGasPvt<Scalar, Opm::gpuistl::GpuView, Opm::gpuistl::GpuView>;

// DryGasPvt<float> is explicitly instantiated in opm-common independently of
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

// Two PVT regions with different table sizes (3 and 7 pressure points). With three
// points, Tabulated1DFunction::findSegmentIndex() never bisects; with seven points the
// bisection iterates twice and takes both sides of the pivot comparison.
constexpr const char* deckString = R"(
RUNSPEC
DIMENS
 1 1 1 /
TABDIMS
 * 2 /
OIL
GAS
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
PVDG
-- P      BG       MUG
   20.0   0.0500   0.0130
  100.0   0.0100   0.0150
  250.0   0.0045   0.0200 /
   10.0   0.1000   0.0120
   40.0   0.0260   0.0128
   80.0   0.0130   0.0140
  130.0   0.0082   0.0158
  190.0   0.0058   0.0180
  260.0   0.0045   0.0205
  350.0   0.0036   0.0235 /
)";

constexpr unsigned numRegions = 2;
// [K], not used by DryGasPvt. A function rather than a constexpr variable: the PVT
// methods take the temperature by reference, which would require the variable to
// exist in device memory.
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

//! The parts of a fluid state used by DryGasPvt::inverseFormationVolumeFactorAndViscosity().
template <class Scalar>
struct TestFluidState {
    using ValueType = Evaluation<Scalar>;
    static constexpr int gasPhaseIdx = 0;

    OPM_HOST_DEVICE const ValueType& pressure(unsigned) const
    { return p; }

    ValueType p;
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

// Functors evaluating one DryGasPvt method. They are templated on the PVT type, so the
// same functor is used for the CPU reference and in the kernel. The functors templated
// on Eval pass the pressure to the method as Eval, so the method is evaluated with a
// plain scalar or with an AD variable, as the fluid system does. Rv and Rvw are ignored
// by DryGasPvt and passed as zero.

template <class Eval>
struct InverseFormationVolumeFactor {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE Eval operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    {
        return pvt.inverseFormationVolumeFactor(pt.regionIdx, Eval(temperature<Scalar>()),
                                                makeVariable<Eval>(pt.p, 0), Eval(0.0), Eval(0.0));
    }
};

template <class Eval>
struct Viscosity {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE Eval operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    {
        return pvt.viscosity(pt.regionIdx, Eval(temperature<Scalar>()),
                             makeVariable<Eval>(pt.p, 0), Eval(0.0), Eval(0.0));
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

//! The methods that are constant zero for dry gas, summed. Each term is zero, so the
//! sum is zero exactly if all of them are.
struct DryGasZeros {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE Evaluation<Scalar> operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    {
        using Eval = Evaluation<Scalar>;
        const Eval T(temperature<Scalar>());
        const Eval p = makeVariable<Eval>(pt.p, 0);
        return pvt.saturationPressure(pt.regionIdx, T, Eval(0.0))
            + pvt.saturatedOilVaporizationFactor(pt.regionIdx, T, p)
            + pvt.saturatedOilVaporizationFactor(pt.regionIdx, T, p, Eval(0.7), Eval(1.0))
            + pvt.saturatedWaterVaporizationFactor(pt.regionIdx, T, p)
            + pvt.saturatedWaterVaporizationFactor(pt.regionIdx, T, p, Eval(0.0));
    }
};

struct GasReferenceDensity {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE Scalar operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    { return pvt.gasReferenceDensity(pt.regionIdx); }
};

struct NumRegions {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE unsigned operator()(const Pvt& pvt, const TestPoint<Scalar>&) const
    { return pvt.numRegions(); }
};

//! The hot path, called by the fluid system.
struct InverseFormationVolumeFactorAndViscosity {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE InvBAndMu<Scalar> operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    {
        TestFluidState<Scalar> fs;
        fs.p = Evaluation<Scalar>::createVariable(pt.p, 0);
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
    DryGasZeros{}(pvt, pt);
    GasReferenceDensity{}(pvt, pt);
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
bool isClose(Scalar a, Scalar b)
{
    if (a == b) {
        return true;
    }
    return std::abs(a - b) <= relTol<Scalar> * std::max(std::abs(a), std::abs(b));
}

template <class Scalar>
bool isClose(const Evaluation<Scalar>& a, const Evaluation<Scalar>& b)
{
    bool close = isClose(a.value(), b.value());
    for (int varIdx = 0; varIdx < Evaluation<Scalar>::numVars; ++varIdx) {
        close = close && isClose(a.derivative(varIdx), b.derivative(varIdx));
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

//! Pressures below both tables, between the nodes and above both tables (extrapolation).
template <class Scalar>
std::vector<TestPoint<Scalar>> gridPoints()
{
    constexpr Scalar bar = 1e5;
    // In region 1, 60-220 bar hit every interior segment of the seven-point table, so
    // the bisection takes both sides of the pivot comparison. 5 bar is below and 400 bar
    // above both tables.
    const std::vector<Scalar> pressures = {5 * bar, 15 * bar, 30 * bar, 60 * bar, 110 * bar,
                                           160 * bar, 220 * bar, 300 * bar, 400 * bar};

    std::vector<TestPoint<Scalar>> points;
    for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
        for (const Scalar p : pressures) {
            points.push_back({regionIdx, p});
        }
    }
    return points;
}

//! Points at every pressure node of each region's tables. The value is continuous at a
//! node, but the derivative depends on the segment that is chosen, so this checks that
//! the CPU and the GPU choose the same one.
template <class Scalar>
std::vector<TestPoint<Scalar>> nodePoints(const CpuPvt<Scalar>& cpuPvt)
{
    std::vector<TestPoint<Scalar>> points;
    for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
        const auto& table = cpuPvt.inverseGasB()[regionIdx];
        for (std::size_t i = 0; i < table.numSamples(); ++i) {
            points.push_back({regionIdx, table.xAt(i)});
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
    BOOST_CHECK_EQUAL(fx.gpuBufPvt.gasMu().size(), numRegions);
    BOOST_CHECK_EQUAL(fx.gpuBufPvt.inverseGasBMu().size(), numRegions);

    BOOST_CHECK_EQUAL(fx.gpuViewPvt.inverseGasB().size(), numRegions);
    BOOST_CHECK_EQUAL(fx.gpuViewPvt.gasMu().size(), numRegions);
    BOOST_CHECK_EQUAL(fx.gpuViewPvt.inverseGasBMu().size(), numRegions);

    // The device code must see the same number of regions as the host.
    compareCpuAndGpu(NumRegions{}, fx.cpuPvt, fx.gpuViewPvt, gridPoints<Scalar>(), "numRegions");
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

//! Compares the methods that depend on the pressure, with the pressure passed as Eval.
template <class Eval, class Scalar>
void comparePressureDependentMethods(const Fixture<Scalar>& fx, const std::vector<TestPoint<Scalar>>& points)
{
    compareCpuAndGpu(InverseFormationVolumeFactor<Eval>{}, fx.cpuPvt, fx.gpuViewPvt, points, "inverseFormationVolumeFactor");
    compareCpuAndGpu(Viscosity<Eval>{}, fx.cpuPvt, fx.gpuViewPvt, points, "viscosity");
    compareCpuAndGpu(SaturatedInverseFormationVolumeFactor<Eval>{}, fx.cpuPvt, fx.gpuViewPvt, points, "saturatedInverseFormationVolumeFactor");
    compareCpuAndGpu(SaturatedViscosity<Eval>{}, fx.cpuPvt, fx.gpuViewPvt, points, "saturatedViscosity");
}

BOOST_AUTO_TEST_CASE_TEMPLATE(ScalarMethodsMatchCpu, Scalar, Scalars)
{
    const Fixture<Scalar> fx;
    const auto points = gridPoints<Scalar>();
    comparePressureDependentMethods<Scalar>(fx, points);
    compareCpuAndGpu(GasReferenceDensity{}, fx.cpuPvt, fx.gpuViewPvt, points, "gasReferenceDensity");
}

BOOST_AUTO_TEST_CASE_TEMPLATE(EvaluationMethodsMatchCpu, Scalar, Scalars)
{
    // The fluid system calls these methods with AD values. Values and derivatives are
    // compared.
    const Fixture<Scalar> fx;
    comparePressureDependentMethods<Evaluation<Scalar>>(fx, gridPoints<Scalar>());
    comparePressureDependentMethods<Evaluation<Scalar>>(fx, nodePoints(fx.cpuPvt));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(InverseFormationVolumeFactorAndViscosityMatchesCpu, Scalar, Scalars)
{
    const Fixture<Scalar> fx;
    compareCpuAndGpu(InverseFormationVolumeFactorAndViscosity{}, fx.cpuPvt, fx.gpuViewPvt,
                     gridPoints<Scalar>(), "inverseFormationVolumeFactorAndViscosity (grid)");
    compareCpuAndGpu(InverseFormationVolumeFactorAndViscosity{}, fx.cpuPvt, fx.gpuViewPvt,
                     nodePoints(fx.cpuPvt), "inverseFormationVolumeFactorAndViscosity (nodes)");
}

BOOST_AUTO_TEST_CASE_TEMPLATE(DryGasMethodsReturnZeroOnGpu, Scalar, Scalars)
{
    // No vaporized oil or water and no saturation pressure, also on the GPU.
    const Fixture<Scalar> fx;
    for (const auto& r : evaluateOnGpu(DryGasZeros{}, fx.gpuViewPvt, gridPoints<Scalar>())) {
        BOOST_CHECK(r == Evaluation<Scalar>(0.0));
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(RegionsAreIndexedCorrectly, Scalar, Scalars)
{
    const Fixture<Scalar> fx;

    // The same pressure in both regions must give different values (the regions have
    // different tables), and each region must match its CPU counterpart.
    constexpr Scalar p = 150e5;
    const std::vector<TestPoint<Scalar>> points = {{0, p}, {1, p}};

    const auto cpuInvB = evaluateOnCpu(InverseFormationVolumeFactor<Scalar>{}, fx.cpuPvt, points);
    const auto cpuMu = evaluateOnCpu(Viscosity<Scalar>{}, fx.cpuPvt, points);
    const auto cpuDensity = evaluateOnCpu(GasReferenceDensity{}, fx.cpuPvt, points);
    BOOST_CHECK_NE(cpuInvB[0], cpuInvB[1]);
    BOOST_CHECK_NE(cpuMu[0], cpuMu[1]);
    BOOST_CHECK_NE(cpuDensity[0], cpuDensity[1]);

    checkResults(cpuInvB, evaluateOnGpu(InverseFormationVolumeFactor<Scalar>{}, fx.gpuViewPvt, points), points,
                 "inverseFormationVolumeFactor (regions)");
    checkResults(cpuMu, evaluateOnGpu(Viscosity<Scalar>{}, fx.gpuViewPvt, points), points,
                 "viscosity (regions)");
    checkResults(cpuDensity, evaluateOnGpu(GasReferenceDensity{}, fx.gpuViewPvt, points), points,
                 "gasReferenceDensity (regions)");
}

BOOST_AUTO_TEST_CASE_TEMPLATE(CallableThroughConstView, Scalar, Scalars)
{
    // That every method compiles on a const view is checked by CallableOnConstView
    // above. Here the view is placed in device memory and used through a const pointer.
    const Fixture<Scalar> fx;
    const auto points = gridPoints<Scalar>();
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

    const auto points = gridPoints<Scalar>();
    compareCpuAndGpu(InverseFormationVolumeFactorAndViscosity{}, fx.cpuPvt, viewBeforeMove, points,
                     "inverseFormationVolumeFactorAndViscosity (view before move)");
    compareCpuAndGpu(InverseFormationVolumeFactorAndViscosity{}, fx.cpuPvt, viewAfterMove, points,
                     "inverseFormationVolumeFactorAndViscosity (view after move)");
    compareCpuAndGpu(InverseFormationVolumeFactorAndViscosity{}, fx.cpuPvt, viewAfterAssignment, points,
                     "inverseFormationVolumeFactorAndViscosity (view after move assignment)");
}
