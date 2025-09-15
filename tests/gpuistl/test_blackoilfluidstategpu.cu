/*
  Copyright 2025 Equinor ASA

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
#include <config.h>

#define BOOST_TEST_MODULE TestBlackOilFluidStateGPU

#include <boost/test/unit_test.hpp>
#include <cuda.h>
#include <cuda_runtime.h>
#include <opm/material/common/HasMemberGeneratorMacros.hpp>
#include <opm/material/fluidstates/BlackOilFluidState.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>
#include <opm/simulators/linalg/gpuistl/gpu_smart_pointer.hpp>
namespace
{

template <class ScalarT>
struct DummyFluidSystem {
    static constexpr auto numPhases = 3u;

    static auto reservoirTemperature(int)
    {
        return ScalarT {0.0};
    }
    static auto enthalpyEqualEnergy()
    {
        return true;
    }
    static auto molarMass(int, int)
    {
        return ScalarT {0.0};
    }
    template <class T>
    static auto viscosity(const T&, int, int)
    {
        return ScalarT {0.0};
    }

    static auto convertRsToXoG(ScalarT, int)
    {
        return ScalarT {0.0};
    }
    static auto convertRvToXg0(ScalarT, int)
    {
        return ScalarT {0.0};
    }
    static auto convertXoGToRs(ScalarT, int)
    {
        return ScalarT {0.0};
    }

    template <class T>
    static auto fugacityCoefficient(const T&, int, int, int)
    {
        return ScalarT {0.0};
    }

    static auto activeToCanonicalPhaseIdx(int)
    {
        return 0u;
    }
    static auto canonicalToActivePhaseIdx(int)
    {
        return 0u;
    }
};

template <class ScalarT>
struct DummyFluidSystemDynamic {
    static constexpr auto numPhases = 3u;

    OPM_HOST_DEVICE auto reservoirTemperature(int) const
    {
        return ScalarT {0.0};
    }
    OPM_HOST_DEVICE auto enthalpyEqualEnergy() const
    {
        return true;
    }
    OPM_HOST_DEVICE auto molarMass(int, int) const
    {
        return ScalarT {0.0};
    }
    template <class T>
    OPM_HOST_DEVICE auto viscosity(const T&, int, int) const
    {
        return ScalarT {someVariable};
    }

    OPM_HOST_DEVICE auto convertRsToXoG(ScalarT, int) const
    {
        return ScalarT {0.0};
    }
    OPM_HOST_DEVICE auto convertRvToXg0(ScalarT, int) const
    {
        return ScalarT {0.0};
    }
    OPM_HOST_DEVICE auto convertXoGToRs(ScalarT, int) const
    {
        return ScalarT {0.0};
    }

    template <class T>
    OPM_HOST_DEVICE auto fugacityCoefficient(const T&, int, int, int) const
    {
        return ScalarT {0.0};
    }

    OPM_HOST_DEVICE auto activeToCanonicalPhaseIdx(int) const
    {
        return 0u;
    }
    OPM_HOST_DEVICE auto canonicalToActivePhaseIdx(int) const
    {
        return 0u;
    }


    double someVariable = 43.2;
};


template <class FluidState>
__global__ void
kernelCreatingBlackoilFluidState()
{
    FluidState state;
}

template <class FluidState, class FluidSystem>
__global__ void
kernelCreatingBlackoilFluidStateDynamic()
{
    FluidSystem system;
    FluidState state(system);
}

template <class FluidState>
__global__ void
kernelSetAndGetTotalSaturation(double saturation, double* readSaturation)
{
    FluidState state;
    state.setTotalSaturation(saturation);
    *readSaturation = state.totalSaturation();
}

template <class FluidState>
__global__ void
getPressure(Opm::gpuistl::PointerView<FluidState> input, std::array<double, 3>* output)
{
    for (int i = 0; i < 3; ++i) {
        (*output)[i] = input->pressure(i);
    }
}

template <class FluidState, class FluidSystem>
__global__ void
getViscosity(FluidSystem input, double* output)
{
    FluidState state(input);
    *output = state.viscosity(0);
}

} // namespace

using ScalarT = double;
using FluidState = Opm::BlackOilFluidState<ScalarT, DummyFluidSystem<ScalarT>>;
using FluidStateDynamic = Opm::BlackOilFluidState<ScalarT, DummyFluidSystemDynamic<ScalarT>>;


BOOST_AUTO_TEST_CASE(TestCreation)
{
    kernelCreatingBlackoilFluidState<FluidState><<<1, 1>>>();
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
    OPM_GPU_SAFE_CALL(cudaGetLastError());
}

BOOST_AUTO_TEST_CASE(TestSaturation)
{
    const double saturation = 0.5;
    auto saturationRead = Opm::gpuistl::make_gpu_unique_ptr<double>(0.0);
    kernelSetAndGetTotalSaturation<FluidState><<<1, 1>>>(saturation, saturationRead.get());
    auto saturationFromGPU = Opm::gpuistl::copyFromGPU(saturationRead);
    BOOST_CHECK_EQUAL(saturationFromGPU, saturation);
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
    OPM_GPU_SAFE_CALL(cudaGetLastError());
}

BOOST_AUTO_TEST_CASE(TestPressure)
{
    FluidState state;
    state.setPressure(0, 1.0);
    state.setPressure(1, 2.0);
    state.setPressure(2, 3.0);

    auto stateGPU = Opm::gpuistl::make_gpu_unique_ptr<FluidState>(state);
    auto output = Opm::gpuistl::make_gpu_unique_ptr<std::array<double, 3>>();

    getPressure<<<1, 1>>>(Opm::gpuistl::make_view(stateGPU), output.get());
    auto outputCPU = Opm::gpuistl::copyFromGPU(output);
    BOOST_CHECK_EQUAL(1.0, outputCPU[0]);
    BOOST_CHECK_EQUAL(2.0, outputCPU[1]);
    BOOST_CHECK_EQUAL(3.0, outputCPU[2]);
}

BOOST_AUTO_TEST_CASE(TestDynamicCreation)
{
    kernelCreatingBlackoilFluidStateDynamic<FluidStateDynamic, DummyFluidSystemDynamic<ScalarT>><<<1, 1>>>();
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
    OPM_GPU_SAFE_CALL(cudaGetLastError());
}

BOOST_AUTO_TEST_CASE(TestPassByValueToGPUDynamic)
{
    DummyFluidSystemDynamic<ScalarT> system;

    system.someVariable = 1234;
    auto output = Opm::gpuistl::make_gpu_unique_ptr<double>();
    getViscosity<FluidStateDynamic><<<1, 1>>>(system, output.get());

    auto outputCPU = Opm::gpuistl::copyFromGPU(output);
    BOOST_CHECK_EQUAL(1234, outputCPU);
}
