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
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using Evaluation = Opm::DenseAd::Evaluation<double, 3>;

using CpuPvt = Opm::LiveOilPvt<double>;
using GpuBufPvt = Opm::LiveOilPvt<double, Opm::gpuistl::GpuBuffer, Opm::gpuistl::GpuView>;
using GpuViewPvt = Opm::LiveOilPvt<double, Opm::gpuistl::GpuView, Opm::gpuistl::GpuView>;

// The views in an owning object point into its own device buffers, so a copy would
// have views into the buffers of the original. The owning object must be move-only.
static_assert(!std::is_copy_constructible_v<GpuBufPvt>);
static_assert(!std::is_copy_assignable_v<GpuBufPvt>);
static_assert(std::is_move_constructible_v<GpuBufPvt>);
static_assert(std::is_move_assignable_v<GpuBufPvt>);
// The CPU object is still copied by OilPvtMultiplexer, and the view is passed by value.
static_assert(std::is_copy_constructible_v<CpuPvt>);
static_assert(std::is_trivially_copyable_v<GpuViewPvt>);

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
OPM_HOST_DEVICE constexpr double temperature()
{ return 350.0; }
constexpr double relTol = 1e-10;

CpuPvt makeCpuPvt()
{
    const auto deck = Opm::Parser{}.parseString(deckString);
    const Opm::EclipseState eclState(deck);
    const Opm::Schedule schedule(deck, eclState, std::make_shared<Opm::Python>());

    CpuPvt pvt;
    pvt.initFromState(eclState, schedule);
    return pvt;
}

//! One evaluation point. Pressure in [Pa], Rs in [sm3/sm3].
struct TestPoint {
    unsigned regionIdx;
    double p;
    double rs;
    double sg;
    double so;
    double soMax;
};

//! The parts of a fluid state used by LiveOilPvt::inverseFormationVolumeFactorAndViscosity().
struct TestFluidState {
    using ValueType = Evaluation;
    static constexpr int oilPhaseIdx = 0;
    static constexpr int gasPhaseIdx = 1;

    OPM_HOST_DEVICE const Evaluation& pressure(unsigned) const
    { return p; }

    OPM_HOST_DEVICE Evaluation Rs() const
    { return rs; }

    OPM_HOST_DEVICE const Evaluation& saturation(unsigned phaseIdx) const
    { return phaseIdx == static_cast<unsigned>(gasPhaseIdx) ? sg : so; }

    Evaluation p;
    Evaluation rs;
    Evaluation sg;
    Evaluation so;
};

//! Result of inverseFormationVolumeFactorAndViscosity(). A plain struct instead of
//! std::pair, because it is copied between host and device with memcpy.
struct InvBAndMu {
    Evaluation invB;
    Evaluation mu;
};

// Functors evaluating one LiveOilPvt method. They are templated on the PVT type, so the
// same functor is used for the CPU reference and in the kernel.

struct InverseFormationVolumeFactor {
    template <class Pvt>
    OPM_HOST_DEVICE double operator()(const Pvt& pvt, const TestPoint& pt) const
    { return pvt.inverseFormationVolumeFactor(pt.regionIdx, temperature(), pt.p, pt.rs); }
};

struct Viscosity {
    template <class Pvt>
    OPM_HOST_DEVICE double operator()(const Pvt& pvt, const TestPoint& pt) const
    { return pvt.viscosity(pt.regionIdx, temperature(), pt.p, pt.rs); }
};

struct SaturatedInverseFormationVolumeFactor {
    template <class Pvt>
    OPM_HOST_DEVICE double operator()(const Pvt& pvt, const TestPoint& pt) const
    { return pvt.saturatedInverseFormationVolumeFactor(pt.regionIdx, temperature(), pt.p); }
};

struct SaturatedViscosity {
    template <class Pvt>
    OPM_HOST_DEVICE double operator()(const Pvt& pvt, const TestPoint& pt) const
    { return pvt.saturatedViscosity(pt.regionIdx, temperature(), pt.p); }
};

struct SaturatedGasDissolutionFactor {
    template <class Pvt>
    OPM_HOST_DEVICE double operator()(const Pvt& pvt, const TestPoint& pt) const
    { return pvt.saturatedGasDissolutionFactor(pt.regionIdx, temperature(), pt.p); }
};

struct SaturatedGasDissolutionFactorVappars {
    template <class Pvt>
    OPM_HOST_DEVICE double operator()(const Pvt& pvt, const TestPoint& pt) const
    { return pvt.saturatedGasDissolutionFactor(pt.regionIdx, temperature(), pt.p, pt.so, pt.soMax); }
};

struct GasReferenceDensity {
    template <class Pvt>
    OPM_HOST_DEVICE double operator()(const Pvt& pvt, const TestPoint& pt) const
    { return pvt.gasReferenceDensity(pt.regionIdx); }
};

struct OilReferenceDensity {
    template <class Pvt>
    OPM_HOST_DEVICE double operator()(const Pvt& pvt, const TestPoint& pt) const
    { return pvt.oilReferenceDensity(pt.regionIdx); }
};

struct VapPar2 {
    template <class Pvt>
    OPM_HOST_DEVICE double operator()(const Pvt& pvt, const TestPoint&) const
    { return pvt.vapPar2(); }
};

//! The hot path, with pressure, Rs and Sg as AD variables.
struct InverseFormationVolumeFactorAndViscosity {
    template <class Pvt>
    OPM_HOST_DEVICE InvBAndMu operator()(const Pvt& pvt, const TestPoint& pt) const
    {
        TestFluidState fs;
        fs.p = Evaluation::createVariable(pt.p, 0);
        fs.rs = Evaluation::createVariable(pt.rs, 1);
        fs.sg = Evaluation::createVariable(pt.sg, 2);
        fs.so = Evaluation(pt.so);
        const auto [invB, mu] = pvt.inverseFormationVolumeFactorAndViscosity(fs, pt.regionIdx);
        return {invB, mu};
    }
};

template <class Functor, class Result>
__global__ void evaluateKernel(Functor f, GpuViewPvt pvt, const TestPoint* points, Result* results, std::size_t n)
{
    const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        results[i] = f(pvt, points[i]);
    }
}

template <class Functor>
auto evaluateOnCpu(Functor f, const CpuPvt& pvt, const std::vector<TestPoint>& points)
{
    using Result = decltype(f(pvt, points.front()));
    std::vector<Result> results;
    results.reserve(points.size());
    for (const auto& pt : points) {
        results.push_back(f(pvt, pt));
    }
    return results;
}

template <class Functor>
auto evaluateOnGpu(Functor f, GpuViewPvt pvt, const std::vector<TestPoint>& points)
{
    using Result = decltype(f(pvt, points.front()));
    const Opm::gpuistl::GpuBuffer<TestPoint> gpuPoints(points);
    Opm::gpuistl::GpuBuffer<Result> gpuResults(points.size());

    constexpr unsigned threadsPerBlock = 128;
    const auto numBlocks = static_cast<unsigned>((points.size() + threadsPerBlock - 1) / threadsPerBlock);
    evaluateKernel<<<numBlocks, threadsPerBlock>>>(f, pvt, gpuPoints.data(), gpuResults.data(), points.size());
    OPM_GPU_SAFE_CALL(cudaPeekAtLastError());
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());

    return gpuResults.asStdVector();
}

bool isClose(double a, double b)
{
    if (a == b) {
        return true;
    }
    return std::abs(a - b) <= relTol * std::max(std::abs(a), std::abs(b));
}

bool isClose(const Evaluation& a, const Evaluation& b)
{
    bool close = isClose(a.value(), b.value());
    for (int varIdx = 0; varIdx < Evaluation::numVars; ++varIdx) {
        close = close && isClose(a.derivative(varIdx), b.derivative(varIdx));
    }
    return close;
}

void checkResults(const std::vector<double>& cpu,
                  const std::vector<double>& gpu,
                  const std::vector<TestPoint>& points,
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

void checkResults(const std::vector<InvBAndMu>& cpu,
                  const std::vector<InvBAndMu>& gpu,
                  const std::vector<TestPoint>& points,
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

template <class Functor>
void compareCpuAndGpu(Functor f,
                      const CpuPvt& cpuPvt,
                      GpuViewPvt gpuPvt,
                      const std::vector<TestPoint>& points,
                      const std::string& name)
{
    checkResults(evaluateOnCpu(f, cpuPvt, points), evaluateOnGpu(f, gpuPvt, points), points, name);
}

//! Points inside the tables, between the Rs rows and outside the tables (extrapolation).
std::vector<TestPoint> gridPoints()
{
    constexpr double bar = 1e5;
    const std::vector<double> pressures = {10 * bar, 45 * bar, 75 * bar, 100 * bar, 125 * bar,
                                           180 * bar, 240 * bar, 320 * bar};
    const std::vector<double> rsValues = {0.0, 15.0, 35.0, 55.0, 80.0, 110.0, 140.0};

    std::vector<TestPoint> points;
    for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
        for (const double p : pressures) {
            for (const double rs : rsValues) {
                for (const double sg : {0.0, 0.3}) {
                    points.push_back({regionIdx, p, rs, sg, 1.0 - sg, 1.0});
                }
            }
        }
    }
    return points;
}

//! Points at and just below the saturated curve, where
//! inverseFormationVolumeFactorAndViscosity() switches between its saturated
//! (Sg > 0 and Rs >= (1 - 1e-10) * RsSat) and undersaturated branch.
std::vector<TestPoint> switchPoints(const CpuPvt& cpuPvt)
{
    constexpr double bar = 1e5;
    std::vector<TestPoint> points;
    for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
        for (const double p : {60 * bar, 120 * bar, 190 * bar}) {
            const double rsSat = cpuPvt.saturatedGasDissolutionFactor(regionIdx, temperature(), p);
            for (const double rs : {rsSat, rsSat * (1.0 - 1e-11), rsSat * (1.0 - 1e-9), 0.5 * rsSat}) {
                points.push_back({regionIdx, p, rs, 0.3, 0.7, 1.0});
                points.push_back({regionIdx, p, rs, 0.0, 1.0, 1.0});
            }
        }
    }
    return points;
}

struct Fixture {
    Fixture()
        : cpuPvt(makeCpuPvt())
        , gpuBufPvt(Opm::gpuistl::copy_to_gpu(cpuPvt))
        , gpuViewPvt(Opm::gpuistl::make_view(gpuBufPvt))
    {}

    CpuPvt cpuPvt;
    GpuBufPvt gpuBufPvt;
    GpuViewPvt gpuViewPvt;
};

} // anonymous namespace

BOOST_FIXTURE_TEST_CASE(SizesAfterCopyToGpu, Fixture)
{
    BOOST_REQUIRE_EQUAL(cpuPvt.numRegions(), numRegions);
    BOOST_CHECK_EQUAL(gpuBufPvt.numRegions(), numRegions);
    BOOST_CHECK_EQUAL(gpuViewPvt.numRegions(), numRegions);

    // Only sizes are safe to read on the host; the elements live in device memory.
    BOOST_CHECK_EQUAL(gpuBufPvt.inverseOilBTable().size(), numRegions);
    BOOST_CHECK_EQUAL(gpuBufPvt.oilMuTable().size(), numRegions);
    BOOST_CHECK_EQUAL(gpuBufPvt.inverseOilBMuTable().size(), numRegions);
    BOOST_CHECK_EQUAL(gpuBufPvt.saturatedOilMuTable().size(), numRegions);
    BOOST_CHECK_EQUAL(gpuBufPvt.inverseSaturatedOilBTable().size(), numRegions);
    BOOST_CHECK_EQUAL(gpuBufPvt.inverseSaturatedOilBMuTable().size(), numRegions);
    BOOST_CHECK_EQUAL(gpuBufPvt.saturatedGasDissolutionFactorTable().size(), numRegions);
    BOOST_CHECK_EQUAL(gpuBufPvt.saturationPressure().size(), numRegions);

    BOOST_CHECK_EQUAL(gpuViewPvt.inverseOilBTable().size(), numRegions);
    BOOST_CHECK_EQUAL(gpuViewPvt.saturationPressure().size(), numRegions);
}

BOOST_AUTO_TEST_CASE(CopyToGpuRejectsUninitializedPvt)
{
    CpuPvt uninitialized;
    BOOST_CHECK_THROW(Opm::gpuistl::copy_to_gpu(uninitialized), std::logic_error);
}

BOOST_FIXTURE_TEST_CASE(ScalarMethodsMatchCpu, Fixture)
{
    const auto points = gridPoints();
    compareCpuAndGpu(InverseFormationVolumeFactor{}, cpuPvt, gpuViewPvt, points, "inverseFormationVolumeFactor");
    compareCpuAndGpu(Viscosity{}, cpuPvt, gpuViewPvt, points, "viscosity");
    compareCpuAndGpu(SaturatedInverseFormationVolumeFactor{}, cpuPvt, gpuViewPvt, points, "saturatedInverseFormationVolumeFactor");
    compareCpuAndGpu(SaturatedViscosity{}, cpuPvt, gpuViewPvt, points, "saturatedViscosity");
    compareCpuAndGpu(SaturatedGasDissolutionFactor{}, cpuPvt, gpuViewPvt, points, "saturatedGasDissolutionFactor");
    compareCpuAndGpu(GasReferenceDensity{}, cpuPvt, gpuViewPvt, points, "gasReferenceDensity");
    compareCpuAndGpu(OilReferenceDensity{}, cpuPvt, gpuViewPvt, points, "oilReferenceDensity");
    compareCpuAndGpu(VapPar2{}, cpuPvt, gpuViewPvt, points, "vapPar2");
}

BOOST_FIXTURE_TEST_CASE(InverseFormationVolumeFactorAndViscosityMatchesCpu, Fixture)
{
    compareCpuAndGpu(InverseFormationVolumeFactorAndViscosity{}, cpuPvt, gpuViewPvt,
                     gridPoints(), "inverseFormationVolumeFactorAndViscosity (grid)");

    const auto points = switchPoints(cpuPvt);

    // Make sure that both branches are exercised: In the saturated branch, the result
    // equals the saturated tables. In the undersaturated branch, it comes from the 2D
    // tables and differs from the saturated tables away from the saturated curve.
    const auto cpuResults = evaluateOnCpu(InverseFormationVolumeFactorAndViscosity{}, cpuPvt, points);
    for (std::size_t i = 0; i < points.size(); ++i) {
        const auto& pt = points[i];
        const double satInvB = cpuPvt.saturatedInverseFormationVolumeFactor(pt.regionIdx, temperature(), pt.p);
        const double rsSat = cpuPvt.saturatedGasDissolutionFactor(pt.regionIdx, temperature(), pt.p);
        if (pt.sg > 0.0 && pt.rs >= rsSat * (1.0 - 1e-10)) {
            BOOST_CHECK_EQUAL(cpuResults[i].invB.value(), satInvB);
        } else if (pt.rs < 0.9 * rsSat) {
            BOOST_CHECK_NE(cpuResults[i].invB.value(), satInvB);
        }
    }

    compareCpuAndGpu(InverseFormationVolumeFactorAndViscosity{}, cpuPvt, gpuViewPvt,
                     points, "inverseFormationVolumeFactorAndViscosity (branch switch)");
}

BOOST_AUTO_TEST_CASE(SaturatedGasDissolutionFactorWithVappars)
{
    auto cpuPvt = makeCpuPvt();
    cpuPvt.setVapPars(2.0, 0.1);
    auto gpuBufPvt = Opm::gpuistl::copy_to_gpu(cpuPvt);
    const auto gpuViewPvt = Opm::gpuistl::make_view(gpuBufPvt);

    constexpr double p = 120e5;
    std::vector<TestPoint> points;
    for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
        // So < SoMax: the VAPPARS correction pow(So/SoMax, vapPar2) is applied
        points.push_back({regionIdx, p, 0.0, 0.0, 0.3, 0.8});
        points.push_back({regionIdx, p, 0.0, 0.0, 0.0005, 0.8}); // So below the cut-off eps
        // So >= SoMax, and SoMax <= 0.01: no correction
        points.push_back({regionIdx, p, 0.0, 0.0, 0.8, 0.8});
        points.push_back({regionIdx, p, 0.0, 0.0, 0.005, 0.009});
    }

    // The correction is active for the first point of each region
    const double uncorrected = cpuPvt.saturatedGasDissolutionFactor(0, temperature(), p);
    BOOST_CHECK_LT(SaturatedGasDissolutionFactorVappars{}(cpuPvt, points[0]), uncorrected);

    compareCpuAndGpu(SaturatedGasDissolutionFactorVappars{}, cpuPvt, gpuViewPvt, points,
                     "saturatedGasDissolutionFactor (VAPPARS)");
    compareCpuAndGpu(VapPar2{}, cpuPvt, gpuViewPvt, points, "vapPar2 (VAPPARS)");
}

BOOST_FIXTURE_TEST_CASE(RegionsAreIndexedCorrectly, Fixture)
{
    // The same (p, Rs) in both regions must give different values (the regions have
    // different tables), and each region must match its CPU counterpart.
    constexpr double p = 150e5;
    constexpr double rs = 40.0;
    const std::vector<TestPoint> points = {{0, p, rs, 0.0, 1.0, 1.0},
                                           {1, p, rs, 0.0, 1.0, 1.0}};

    const auto cpuInvB = evaluateOnCpu(InverseFormationVolumeFactor{}, cpuPvt, points);
    const auto cpuMu = evaluateOnCpu(Viscosity{}, cpuPvt, points);
    BOOST_CHECK_NE(cpuInvB[0], cpuInvB[1]);
    BOOST_CHECK_NE(cpuMu[0], cpuMu[1]);

    checkResults(cpuInvB, evaluateOnGpu(InverseFormationVolumeFactor{}, gpuViewPvt, points), points,
                 "inverseFormationVolumeFactor (regions)");
    checkResults(cpuMu, evaluateOnGpu(Viscosity{}, gpuViewPvt, points), points,
                 "viscosity (regions)");
    compareCpuAndGpu(GasReferenceDensity{}, cpuPvt, gpuViewPvt, points, "gasReferenceDensity (regions)");
    compareCpuAndGpu(OilReferenceDensity{}, cpuPvt, gpuViewPvt, points, "oilReferenceDensity (regions)");
}

BOOST_FIXTURE_TEST_CASE(ViewsStayValidAfterMovingOwner, Fixture)
{
    // Moving the owning object moves the GpuBuffers, which keeps their device memory.
    // Both a view created before the move and one created after it must stay valid.
    const auto viewBeforeMove = gpuViewPvt;
    GpuBufPvt movedPvt = std::move(gpuBufPvt);
    const auto viewAfterMove = Opm::gpuistl::make_view(movedPvt);

    const auto points = gridPoints();
    compareCpuAndGpu(InverseFormationVolumeFactorAndViscosity{}, cpuPvt, viewBeforeMove, points,
                     "inverseFormationVolumeFactorAndViscosity (view before move)");
    compareCpuAndGpu(InverseFormationVolumeFactorAndViscosity{}, cpuPvt, viewAfterMove, points,
                     "inverseFormationVolumeFactorAndViscosity (view after move)");
}
