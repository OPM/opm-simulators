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
 * \brief Unit tests for LiveOilPvt on the GPU.
 *
 * A CPU LiveOilPvt is initialized from a deck with two PVT regions of different
 * table sizes, copied to the GPU with copy_to_gpu()/make_view(), and evaluated
 * in kernels. All GPU results are compared to the same evaluations on the CPU.
 */
#include <config.h>

#define BOOST_TEST_MODULE TestGpuLiveOilPvt
#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/fluidsystems/blackoilpvt/LiveOilPvt.hpp>

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
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

template <class Scalar>
using Evaluation = Opm::DenseAd::Evaluation<Scalar, 3>;

template <class Scalar>
using CpuPvt = Opm::LiveOilPvt<Scalar>;
template <class Scalar>
using GpuBufPvt = Opm::LiveOilPvt<Scalar, Opm::gpuistl::GpuBuffer, Opm::gpuistl::GpuView>;
template <class Scalar>
using GpuViewPvt = Opm::LiveOilPvt<Scalar, Opm::gpuistl::GpuView, Opm::gpuistl::GpuView>;

// LiveOilPvt<float> is explicitly instantiated in opm-common independently of
// FLOW_INSTANTIATE_FLOAT, so both types are always tested.
using Scalars = boost::mpl::list<double, float>;

template <class Scalar>
struct CheckCopySemantics {
    // The views in an owning object point into its own device buffers, so a copy would
    // have views into the buffers of the original. The owning object must be move-only.
    static_assert(!std::is_copy_constructible_v<GpuBufPvt<Scalar>>);
    static_assert(!std::is_copy_assignable_v<GpuBufPvt<Scalar>>);
    static_assert(std::is_move_constructible_v<GpuBufPvt<Scalar>>);
    static_assert(std::is_move_assignable_v<GpuBufPvt<Scalar>>);
    // The CPU object is still copied by OilPvtMultiplexer, and the view is passed by value.
    static_assert(std::is_copy_constructible_v<CpuPvt<Scalar>>);
    static_assert(std::is_trivially_copyable_v<GpuViewPvt<Scalar>>);
};
template struct CheckCopySemantics<double>;
template struct CheckCopySemantics<float>;

// Two PVT regions with different numbers of Rs rows (3 and 5) and different numbers
// of undersaturated points per row. Rows with only the saturated point are extended
// by LiveOilPvt::extendPvtoTable_(). Region 1 has a row with four pressure points, so
// that the bisection in the pressure direction is exercised.
constexpr const char* deckString = R"(
RUNSPEC
DIMENS
 1 1 1 /
TABDIMS
 * 2 /
OIL
GAS
WATER
DISGAS
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
PVTO
-- RS     P      BO     MU
   20.0   50.0  1.10   1.20 /
   50.0  100.0  1.20   1.00 /
   90.0  150.0  1.30   0.80
         200.0  1.28   0.85
         250.0  1.26   0.90 /
/
   10.0   40.0  1.05   1.50
          80.0  1.04   1.55 /
   30.0   70.0  1.12   1.30 /
   60.0  110.0  1.21   1.10
         160.0  1.20   1.12
         210.0  1.19   1.15
         260.0  1.18   1.19 /
   90.0  140.0  1.30   0.95 /
  120.0  170.0  1.40   0.85
         230.0  1.38   0.88 /
/
)";

constexpr unsigned numRegions = 2;
// [K], not used by LiveOilPvt. A function rather than a constexpr variable: the PVT
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

//! One evaluation point. Pressure in [Pa], Rs in [sm3/sm3].
template <class Scalar>
struct TestPoint {
    unsigned regionIdx;
    Scalar p;
    Scalar rs;
    Scalar sg;
    Scalar so;
    Scalar soMax;
};

//! The parts of a fluid state used by LiveOilPvt::inverseFormationVolumeFactorAndViscosity().
template <class Scalar>
struct TestFluidState {
    using ValueType = Evaluation<Scalar>;
    static constexpr int oilPhaseIdx = 0;
    static constexpr int gasPhaseIdx = 1;

    OPM_HOST_DEVICE const ValueType& pressure(unsigned) const
    { return p; }

    OPM_HOST_DEVICE ValueType Rs() const
    { return rs; }

    OPM_HOST_DEVICE const ValueType& saturation(unsigned phaseIdx) const
    { return phaseIdx == static_cast<unsigned>(gasPhaseIdx) ? sg : so; }

    ValueType p;
    ValueType rs;
    ValueType sg;
    ValueType so;
};

//! Result of inverseFormationVolumeFactorAndViscosity(). A plain struct instead of
//! std::pair, because it is copied between host and device with memcpy.
template <class Scalar>
struct InvBAndMu {
    Evaluation<Scalar> invB;
    Evaluation<Scalar> mu;
};

// Functors evaluating one LiveOilPvt method. They are templated on the PVT type, so the
// same functor is used for the CPU reference and in the kernel.

struct InverseFormationVolumeFactor {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE Scalar operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    { return pvt.inverseFormationVolumeFactor(pt.regionIdx, temperature<Scalar>(), pt.p, pt.rs); }
};

struct Viscosity {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE Scalar operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    { return pvt.viscosity(pt.regionIdx, temperature<Scalar>(), pt.p, pt.rs); }
};

struct SaturatedInverseFormationVolumeFactor {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE Scalar operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    { return pvt.saturatedInverseFormationVolumeFactor(pt.regionIdx, temperature<Scalar>(), pt.p); }
};

struct SaturatedViscosity {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE Scalar operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    { return pvt.saturatedViscosity(pt.regionIdx, temperature<Scalar>(), pt.p); }
};

struct SaturatedGasDissolutionFactor {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE Scalar operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    { return pvt.saturatedGasDissolutionFactor(pt.regionIdx, temperature<Scalar>(), pt.p); }
};

struct SaturatedGasDissolutionFactorVappars {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE Scalar operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    { return pvt.saturatedGasDissolutionFactor(pt.regionIdx, temperature<Scalar>(), pt.p, pt.so, pt.soMax); }
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

struct VapPar2 {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE Scalar operator()(const Pvt& pvt, const TestPoint<Scalar>&) const
    { return pvt.vapPar2(); }
};

//! The hot path, with pressure, Rs and Sg as AD variables.
struct InverseFormationVolumeFactorAndViscosity {
    template <class Pvt, class Scalar>
    OPM_HOST_DEVICE InvBAndMu<Scalar> operator()(const Pvt& pvt, const TestPoint<Scalar>& pt) const
    {
        using Eval = Evaluation<Scalar>;
        TestFluidState<Scalar> fs;
        fs.p = Eval::createVariable(pt.p, 0);
        fs.rs = Eval::createVariable(pt.rs, 1);
        fs.sg = Eval::createVariable(pt.sg, 2);
        fs.so = Eval(pt.so);
        const auto [invB, mu] = pvt.inverseFormationVolumeFactorAndViscosity(fs, pt.regionIdx);
        return {invB, mu};
    }
};

template <class Functor, class Pvt, class Point, class Result>
__global__ void evaluateKernel(Functor f, Pvt pvt, const Point* points, Result* results, std::size_t n)
{
    const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        results[i] = f(pvt, points[i]);
    }
}

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

    constexpr unsigned threadsPerBlock = 128;
    const auto numBlocks = static_cast<unsigned>((points.size() + threadsPerBlock - 1) / threadsPerBlock);
    evaluateKernel<<<numBlocks, threadsPerBlock>>>(f, pvt, gpuPoints.data(), gpuResults.data(), points.size());
    OPM_GPU_SAFE_CALL(cudaPeekAtLastError());
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());

    return gpuResults.asStdVector();
}

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

template <std::floating_point Scalar>
void checkResults(const std::vector<Scalar>& cpu,
                  const std::vector<Scalar>& gpu,
                  const std::vector<TestPoint<Scalar>>& points,
                  const std::string& name)
{
    BOOST_REQUIRE_EQUAL(cpu.size(), gpu.size());
    for (std::size_t i = 0; i < cpu.size(); ++i) {
        BOOST_CHECK_MESSAGE(isClose(cpu[i], gpu[i]),
                            name << " at point " << i << " (region " << points[i].regionIdx
                                 << ", p = " << points[i].p << ", Rs = " << points[i].rs
                                 << "): CPU = " << cpu[i] << ", GPU = " << gpu[i]);
    }
}

template <class Scalar>
void checkResults(const std::vector<InvBAndMu<Scalar>>& cpu,
                  const std::vector<InvBAndMu<Scalar>>& gpu,
                  const std::vector<TestPoint<Scalar>>& points,
                  const std::string& name)
{
    BOOST_REQUIRE_EQUAL(cpu.size(), gpu.size());
    for (std::size_t i = 0; i < cpu.size(); ++i) {
        BOOST_CHECK_MESSAGE(isClose(cpu[i].invB, gpu[i].invB) && isClose(cpu[i].mu, gpu[i].mu),
                            name << " at point " << i << " (region " << points[i].regionIdx
                                 << ", p = " << points[i].p << ", Rs = " << points[i].rs
                                 << ", Sg = " << points[i].sg << "): CPU = (" << cpu[i].invB
                                 << ", " << cpu[i].mu << "), GPU = (" << gpu[i].invB
                                 << ", " << gpu[i].mu << ")");
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

//! Points inside the tables, between the Rs rows and outside the tables (extrapolation).
template <class Scalar>
std::vector<TestPoint<Scalar>> gridPoints()
{
    constexpr Scalar bar = 1e5;
    const std::vector<Scalar> pressures = {10 * bar, 45 * bar, 75 * bar, 100 * bar, 125 * bar,
                                           180 * bar, 240 * bar, 320 * bar};
    // None of the Rs values is on the saturated curve at any of the pressures, so that
    // the branch in inverseFormationVolumeFactorAndViscosity() is not decided by the
    // last bits of RsSat (switchPoints() tests the switch with suitable margins).
    const std::vector<Scalar> rsValues = {0.0, 15.0, 37.0, 55.0, 80.0, 110.0, 140.0};

    std::vector<TestPoint<Scalar>> points;
    for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
        for (const Scalar p : pressures) {
            for (const Scalar rs : rsValues) {
                for (const Scalar sg : {Scalar{0.0}, Scalar{0.3}}) {
                    points.push_back({regionIdx, p, rs, sg, 1 - sg, 1});
                }
            }
        }
    }
    return points;
}

//! Points at and just below the saturated curve, where
//! inverseFormationVolumeFactorAndViscosity() switches between its saturated
//! (Sg > 0 and Rs >= (1 - 1e-10) * RsSat) and undersaturated branch.
template <class Scalar>
std::vector<TestPoint<Scalar>> switchPoints(const CpuPvt<Scalar>& cpuPvt)
{
    // Rs relative to RsSat. In float, (1 - 1e-10) rounds to 1, so the switch is at
    // Rs >= RsSat. A point exactly on the curve could then take different branches on
    // the CPU and the GPU if RsSat differs in the last bit, so for float the points
    // keep a safe distance from the curve.
    const std::vector<Scalar> rsFactors = std::is_same_v<Scalar, float>
        ? std::vector<Scalar>{1.0 + 1e-4, 1.0 - 1e-4, 0.5}
        : std::vector<Scalar>{1.0, 1.0 - 1e-11, 1.0 - 1e-9, 0.5};

    constexpr Scalar bar = 1e5;
    std::vector<TestPoint<Scalar>> points;
    for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
        for (const Scalar p : {60 * bar, 120 * bar, 190 * bar}) {
            const Scalar rsSat = cpuPvt.saturatedGasDissolutionFactor(regionIdx, temperature<Scalar>(), p);
            for (const Scalar factor : rsFactors) {
                points.push_back({regionIdx, p, factor * rsSat, 0.3, 0.7, 1});
                points.push_back({regionIdx, p, factor * rsSat, 0, 1, 1});
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
    BOOST_CHECK_EQUAL(fx.gpuBufPvt.inverseOilBTable().size(), numRegions);
    BOOST_CHECK_EQUAL(fx.gpuBufPvt.oilMuTable().size(), numRegions);
    BOOST_CHECK_EQUAL(fx.gpuBufPvt.inverseOilBMuTable().size(), numRegions);
    BOOST_CHECK_EQUAL(fx.gpuBufPvt.saturatedOilMuTable().size(), numRegions);
    BOOST_CHECK_EQUAL(fx.gpuBufPvt.inverseSaturatedOilBTable().size(), numRegions);
    BOOST_CHECK_EQUAL(fx.gpuBufPvt.inverseSaturatedOilBMuTable().size(), numRegions);
    BOOST_CHECK_EQUAL(fx.gpuBufPvt.saturatedGasDissolutionFactorTable().size(), numRegions);
    BOOST_CHECK_EQUAL(fx.gpuBufPvt.saturationPressure().size(), numRegions);

    BOOST_CHECK_EQUAL(fx.gpuViewPvt.inverseOilBTable().size(), numRegions);
    BOOST_CHECK_EQUAL(fx.gpuViewPvt.saturationPressure().size(), numRegions);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(CopyToGpuRejectsUninitializedPvt, Scalar, Scalars)
{
    CpuPvt<Scalar> uninitialized;
    BOOST_CHECK_THROW(Opm::gpuistl::copy_to_gpu(uninitialized), std::logic_error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(ScalarMethodsMatchCpu, Scalar, Scalars)
{
    const Fixture<Scalar> fx;
    const auto points = gridPoints<Scalar>();
    compareCpuAndGpu(InverseFormationVolumeFactor{}, fx.cpuPvt, fx.gpuViewPvt, points, "inverseFormationVolumeFactor");
    compareCpuAndGpu(Viscosity{}, fx.cpuPvt, fx.gpuViewPvt, points, "viscosity");
    compareCpuAndGpu(SaturatedInverseFormationVolumeFactor{}, fx.cpuPvt, fx.gpuViewPvt, points, "saturatedInverseFormationVolumeFactor");
    compareCpuAndGpu(SaturatedViscosity{}, fx.cpuPvt, fx.gpuViewPvt, points, "saturatedViscosity");
    compareCpuAndGpu(SaturatedGasDissolutionFactor{}, fx.cpuPvt, fx.gpuViewPvt, points, "saturatedGasDissolutionFactor");
    compareCpuAndGpu(GasReferenceDensity{}, fx.cpuPvt, fx.gpuViewPvt, points, "gasReferenceDensity");
    compareCpuAndGpu(OilReferenceDensity{}, fx.cpuPvt, fx.gpuViewPvt, points, "oilReferenceDensity");
    compareCpuAndGpu(VapPar2{}, fx.cpuPvt, fx.gpuViewPvt, points, "vapPar2");
}

BOOST_AUTO_TEST_CASE_TEMPLATE(InverseFormationVolumeFactorAndViscosityMatchesCpu, Scalar, Scalars)
{
    const Fixture<Scalar> fx;
    compareCpuAndGpu(InverseFormationVolumeFactorAndViscosity{}, fx.cpuPvt, fx.gpuViewPvt,
                     gridPoints<Scalar>(), "inverseFormationVolumeFactorAndViscosity (grid)");

    const auto points = switchPoints(fx.cpuPvt);

    // Make sure that both branches are exercised: In the saturated branch, the result
    // equals the saturated tables. In the undersaturated branch, it comes from the 2D
    // tables and differs from the saturated tables away from the saturated curve.
    const auto cpuResults = evaluateOnCpu(InverseFormationVolumeFactorAndViscosity{}, fx.cpuPvt, points);
    for (std::size_t i = 0; i < points.size(); ++i) {
        const auto& pt = points[i];
        const Scalar satInvB = fx.cpuPvt.saturatedInverseFormationVolumeFactor(pt.regionIdx, temperature<Scalar>(), pt.p);
        const Scalar rsSat = fx.cpuPvt.saturatedGasDissolutionFactor(pt.regionIdx, temperature<Scalar>(), pt.p);
        if (pt.sg > 0 && pt.rs >= Scalar(1.0 - 1e-10) * rsSat) {
            BOOST_CHECK_EQUAL(cpuResults[i].invB.value(), satInvB);
        } else if (pt.rs < Scalar(0.9) * rsSat) {
            BOOST_CHECK_NE(cpuResults[i].invB.value(), satInvB);
        }
    }

    compareCpuAndGpu(InverseFormationVolumeFactorAndViscosity{}, fx.cpuPvt, fx.gpuViewPvt,
                     points, "inverseFormationVolumeFactorAndViscosity (branch switch)");
}

BOOST_AUTO_TEST_CASE_TEMPLATE(SaturatedGasDissolutionFactorWithVappars, Scalar, Scalars)
{
    auto cpuPvt = makeCpuPvt<Scalar>();
    cpuPvt.setVapPars(2.0, 0.1);
    auto gpuBufPvt = Opm::gpuistl::copy_to_gpu(cpuPvt);
    const auto gpuViewPvt = Opm::gpuistl::make_view(gpuBufPvt);

    constexpr Scalar p = 120e5;
    std::vector<TestPoint<Scalar>> points;
    for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
        // So < SoMax: the VAPPARS correction pow(So/SoMax, vapPar2) is applied
        points.push_back({regionIdx, p, 0, 0, 0.3, 0.8});
        points.push_back({regionIdx, p, 0, 0, 0.0005, 0.8}); // So below the cut-off eps
        // So >= SoMax, and SoMax <= 0.01: no correction
        points.push_back({regionIdx, p, 0, 0, 0.8, 0.8});
        points.push_back({regionIdx, p, 0, 0, 0.005, 0.009});
    }

    // The correction is active for the first point of each region
    const Scalar uncorrected = cpuPvt.saturatedGasDissolutionFactor(0, temperature<Scalar>(), p);
    BOOST_CHECK_LT(SaturatedGasDissolutionFactorVappars{}(cpuPvt, points[0]), uncorrected);

    compareCpuAndGpu(SaturatedGasDissolutionFactorVappars{}, cpuPvt, gpuViewPvt, points,
                     "saturatedGasDissolutionFactor (VAPPARS)");
    compareCpuAndGpu(VapPar2{}, cpuPvt, gpuViewPvt, points, "vapPar2 (VAPPARS)");
}

BOOST_AUTO_TEST_CASE_TEMPLATE(RegionsAreIndexedCorrectly, Scalar, Scalars)
{
    const Fixture<Scalar> fx;

    // The same (p, Rs) in both regions must give different values (the regions have
    // different tables), and each region must match its CPU counterpart.
    constexpr Scalar p = 150e5;
    constexpr Scalar rs = 40.0;
    const std::vector<TestPoint<Scalar>> points = {{0, p, rs, 0, 1, 1},
                                                   {1, p, rs, 0, 1, 1}};

    const auto cpuInvB = evaluateOnCpu(InverseFormationVolumeFactor{}, fx.cpuPvt, points);
    const auto cpuMu = evaluateOnCpu(Viscosity{}, fx.cpuPvt, points);
    BOOST_CHECK_NE(cpuInvB[0], cpuInvB[1]);
    BOOST_CHECK_NE(cpuMu[0], cpuMu[1]);

    checkResults(cpuInvB, evaluateOnGpu(InverseFormationVolumeFactor{}, fx.gpuViewPvt, points), points,
                 "inverseFormationVolumeFactor (regions)");
    checkResults(cpuMu, evaluateOnGpu(Viscosity{}, fx.gpuViewPvt, points), points,
                 "viscosity (regions)");
    compareCpuAndGpu(GasReferenceDensity{}, fx.cpuPvt, fx.gpuViewPvt, points, "gasReferenceDensity (regions)");
    compareCpuAndGpu(OilReferenceDensity{}, fx.cpuPvt, fx.gpuViewPvt, points, "oilReferenceDensity (regions)");
}

BOOST_AUTO_TEST_CASE_TEMPLATE(ViewsStayValidAfterMovingOwner, Scalar, Scalars)
{
    Fixture<Scalar> fx;

    // Moving the owning object moves the GpuBuffers, which keeps their device memory.
    // Both a view created before the move and one created after it must stay valid.
    const auto viewBeforeMove = fx.gpuViewPvt;
    GpuBufPvt<Scalar> movedPvt = std::move(fx.gpuBufPvt);
    const auto viewAfterMove = Opm::gpuistl::make_view(movedPvt);

    const auto points = gridPoints<Scalar>();
    compareCpuAndGpu(InverseFormationVolumeFactorAndViscosity{}, fx.cpuPvt, viewBeforeMove, points,
                     "inverseFormationVolumeFactorAndViscosity (view before move)");
    compareCpuAndGpu(InverseFormationVolumeFactorAndViscosity{}, fx.cpuPvt, viewAfterMove, points,
                     "inverseFormationVolumeFactorAndViscosity (view after move)");
}
