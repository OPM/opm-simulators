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
 * \brief Checks that GpuFlowProblem::refreshMixingLimits copies new values of
 *        the time-dependent mixing limits into the existing device buffers, so
 *        that a view taken before the refresh reads the new values.
 */
#include <config.h>

#define BOOST_TEST_MODULE TestGpuFlowProblemRefresh

#include <boost/test/unit_test.hpp>

#include <cuda_runtime.h>

#include <opm/common/utility/gpuDecorators.hpp>

#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>

#include <cstddef>
#include <vector>

namespace {

/// GpuFlowProblem needs its material law manager only for the relative
/// permeabilities and the saturation regions, which this test does not use.
/// The same type serves as the buffer and as the view.
struct StubMaterialLawManager
{
    using MaterialLawParams = int;

    StubMaterialLawManager() = default;

    template <class CpuMaterialLawManager>
    StubMaterialLawManager(const CpuMaterialLawManager&, std::size_t)
    {
    }
};

} // namespace

// Declared before GpuFlowProblem.hpp, whose make_view calls it by qualified name.
namespace Opm::gpuistl {
inline StubMaterialLawManager make_view(StubMaterialLawManager&)
{
    return {};
}
} // namespace Opm::gpuistl

#include <opm/simulators/flow/GpuFlowProblem.hpp>

namespace {

using Scalar = double;

/// The part of the CPU FlowProblem interface that the GpuFlowProblem
/// constructor reads, with mutable mixing limits. timeIdx 1 returns different
/// values, to check that the refresh reads the requested time index.
struct StubCpuProblem
{
    struct Model
    {
        std::size_t numDof;
        std::size_t numGridDof() const { return numDof; }
    };

    Model model_;
    std::vector<Scalar> maxRs;
    std::vector<Scalar> maxRv;
    std::vector<Scalar> maxSo;

    explicit StubCpuProblem(std::size_t n)
        : model_{n}
        , maxRs(n)
        , maxRv(n)
        , maxSo(n)
    {
    }

    const Model& model() const { return model_; }
    int materialLawManager() const { return 0; }
    bool usesDefaultRockCompaction() const { return true; }

    Scalar porosity(unsigned i, unsigned) const { return 0.1 + 0.01 * i; }
    Scalar rockCompressibility(unsigned) const { return 0.0; }
    Scalar rockReferencePressure(unsigned) const { return 1e5; }
    // Not read for a problem without thermal laws, but named by the constructor.
    Scalar rockFraction(unsigned, unsigned) const { return 0.0; }
    unsigned pvtRegionIndex(unsigned) const { return 0; }

    Scalar maxOilSaturation(unsigned i) const { return maxSo[i]; }
    Scalar maxGasDissolutionFactor(unsigned timeIdx, unsigned i) const
    {
        return timeIdx == 0 ? maxRs[i] : -maxRs[i];
    }
    Scalar maxOilVaporizationFactor(unsigned timeIdx, unsigned i) const
    {
        return timeIdx == 0 ? maxRv[i] : -maxRv[i];
    }

    void set(Scalar offset)
    {
        for (std::size_t i = 0; i < maxRs.size(); ++i) {
            maxRs[i] = offset + 10.0 * i;
            maxRv[i] = offset + 0.001 * i;
            maxSo[i] = offset / 1000.0 + 0.01 * i;
        }
    }
};

using GpuProblemBuffer
    = Opm::GpuFlowProblem<Scalar, StubMaterialLawManager, Opm::gpuistl::GpuBuffer>;

constexpr std::size_t numFields = 4;

template <class ProblemView>
__global__ void readProblem(ProblemView problem, Opm::gpuistl::GpuView<Scalar> out, std::size_t n)
{
    const std::size_t i = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) {
        return;
    }
    out[numFields * i + 0] = problem.maxGasDissolutionFactor(0u, i);
    out[numFields * i + 1] = problem.maxOilVaporizationFactor(0u, i);
    out[numFields * i + 2] = problem.maxOilSaturation(i);
    out[numFields * i + 3] = problem.porosity(i, 0u);
}

template <class ProblemView>
std::vector<Scalar> readOnDevice(const ProblemView& view, std::size_t n)
{
    Opm::gpuistl::GpuBuffer<Scalar> out(numFields * n);
    readProblem<<<1, static_cast<unsigned>(n)>>>(
        view, Opm::gpuistl::GpuView<Scalar>(out.data(), out.size()), n);
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
    OPM_GPU_SAFE_CALL(cudaGetLastError());
    return out.asStdVector();
}

void checkValues(const std::vector<Scalar>& values,
                 const StubCpuProblem& cpu,
                 unsigned timeIdx)
{
    const std::size_t n = cpu.model().numGridDof();
    BOOST_REQUIRE_EQUAL(values.size(), numFields * n);
    for (std::size_t i = 0; i < n; ++i) {
        const unsigned u = static_cast<unsigned>(i);
        BOOST_TEST_CONTEXT("cell " << i)
        {
            BOOST_CHECK_EQUAL(values[numFields * i + 0], cpu.maxGasDissolutionFactor(timeIdx, u));
            BOOST_CHECK_EQUAL(values[numFields * i + 1], cpu.maxOilVaporizationFactor(timeIdx, u));
            BOOST_CHECK_EQUAL(values[numFields * i + 2], cpu.maxOilSaturation(u));
            BOOST_CHECK_EQUAL(values[numFields * i + 3], cpu.porosity(u, 0u));
        }
    }
}

} // namespace

BOOST_AUTO_TEST_CASE(RefreshUpdatesExistingView)
{
    constexpr std::size_t n = 7;
    StubCpuProblem cpu(n);
    cpu.set(100.0);

    GpuProblemBuffer buffer(cpu);
    const auto view = Opm::gpuistl::make_view(buffer);
    checkValues(readOnDevice(view, n), cpu, 0u);

    // New mixing limits on the CPU side (as after a new time step): the view
    // taken before the refresh must see them.
    cpu.set(200.0);
    buffer.refreshMixingLimits(cpu, 0u);
    checkValues(readOnDevice(view, n), cpu, 0u);

    // A chopped time step goes back to older values.
    cpu.set(150.0);
    buffer.refreshMixingLimits(cpu, 0u);
    checkValues(readOnDevice(view, n), cpu, 0u);
}

BOOST_AUTO_TEST_CASE(RefreshReadsRequestedTimeIndex)
{
    constexpr std::size_t n = 3;
    StubCpuProblem cpu(n);
    cpu.set(50.0);

    GpuProblemBuffer buffer(cpu);
    buffer.refreshMixingLimits(cpu, 1u);
    checkValues(readOnDevice(Opm::gpuistl::make_view(buffer), n), cpu, 1u);
}
