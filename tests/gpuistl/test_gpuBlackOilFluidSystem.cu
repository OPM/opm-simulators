// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).
  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.
  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/

/*
    This file is based on a copy of the regular blackoilfluidsystem test
    This file contributes extra assertions that the values match on GPU and CPU
*/


/*!
 * \file
 *
 * \brief This is the unit test for the black oil fluid system on the GPU.
 *
 * This test requires the presence of opm-parser.
 */
#include "config.h"
#include <fmt/format.h>
#include <iostream>

#if !HAVE_ECL_INPUT
#error "The test for the black oil fluid system classes requires ecl input support in opm-common"
#endif

#include <boost/mpl/list.hpp>

#define BOOST_TEST_MODULE EclBlackOilFluidSystem
#include <boost/test/unit_test.hpp>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystemNonStatic.hpp>
#include <opm/material/fluidstates/BlackOilFluidState.hpp>
#include <opm/material/fluidsystems/blackoilpvt/GasPvtMultiplexer.hpp>
#include <opm/material/fluidsystems/blackoilpvt/WaterPvtMultiplexer.hpp>
#include <opm/material/components/CO2Tables.hpp>
#include <opm/material/densead/Evaluation.hpp>

#include <opm/material/common/Valgrind.hpp>

#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Python/Python.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <type_traits>
#include <cmath>

#include <opm/simulators/linalg/gpuistl/GpuView.hpp>
#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/gpu_smart_pointer.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>

static constexpr const char* deckString1 =
"-- =============== RUNSPEC\n"
"RUNSPEC\n"
"DIMENS\n"
"3 3 3 /\n"
"EQLDIMS\n"
"/\n"
"TABDIMS\n"
"/\n"
"WATER\n"
"GAS\n"
"CO2STORE\n"
"METRIC\n"
"-- =============== GRID\n"
"GRID\n"
"GRIDFILE\n"
"0 0 /\n"
"DX\n"
"27*1 /\n"
"DY\n"
"27*1 /\n"
"DZ\n"
"27*1 /\n"
"TOPS\n"
"9*0 /\n"
"PERMX\n"
"27*1013.25 /\n"
"PORO\n"
"27*0.25 /\n"
"COPY\n"
"PERMX PERMY /\n"
"PERMX PERMZ /\n"
"/\n"
"-- =============== PROPS\n"
"PROPS\n"
"SGWFN\n"
"0.000000E+00 0.000000E+00 1.000000E+00 3.060000E-02\n"
"1.000000E+00 1.000000E+00 0.000000E+00 3.060000E-01 /\n"
"-- =============== SOLUTION\n"
"SOLUTION\n"
"RPTRST\n"
"'BASIC=0' /\n"
"EQUIL\n"
"0 300 100 0 0 0 1 1 0 /\n"
"-- =============== SCHEDULE\n"
"SCHEDULE\n"
"RPTRST\n"
"'BASIC=0' /\n"
"TSTEP\n"
"1 /";

using GpuBufCo2Tables = Opm::CO2Tables<double, Opm::gpuistl::GpuBuffer>;
using GpuBufBrineCo2Pvt = Opm::BrineCo2Pvt<double, Opm::gpuistl::GpuBuffer>;
using FluidSystem = Opm::BlackOilFluidSystem<double>;
using Evaluation = Opm::DenseAd::Evaluation<double,2>;
using Scalar = typename Opm::MathToolbox<Evaluation>::Scalar;

// checks that we can access value stored as scalar
template <class IndexTraits>
__global__ void getReservoirTemperature(Opm::BlackOilFluidSystemNonStatic<double, IndexTraits, Opm::gpuistl::GpuView> fs, double* res)
{
  *res = fs.reservoirTemperature();
}

// checks that we can access value stored in vectors/buffer/views
template <class IndexTraits>
__global__ void getReferenceDensity(Opm::BlackOilFluidSystemNonStatic<double, IndexTraits, Opm::gpuistl::GpuView> fs, double* res)
{
  *res = fs.referenceDensity(0, 0);
}

// check that we can correctly compute values that require using pvt multiplexers
template <class IndexTraits>
__global__ void getReferenceDensityFromGasPvt(Opm::BlackOilFluidSystemNonStatic<double, IndexTraits, Opm::gpuistl::GpuView> fs, double* res)
{
  *res = fs.gasPvt().gasReferenceDensity(0);
}

BOOST_AUTO_TEST_CASE(BlackOilFluidSystemOnGpu)
{
    Opm::Parser parser;

    auto deck = parser.parseString(deckString1);
    auto python = std::make_shared<Opm::Python>();
    Opm::EclipseState eclState(deck);
    Opm::Schedule schedule(deck, eclState, python);

    FluidSystem::initFromState(eclState, schedule);

    auto& dynamicFluidSystem = FluidSystem::getNonStaticInstance();

    auto dynamicGpuFluidSystemBuffer = ::Opm::gpuistl::copy_to_gpu(dynamicFluidSystem);
    auto dynamicGpuFluidSystemView = ::Opm::gpuistl::make_view(dynamicGpuFluidSystemBuffer);

    // create a parameter cache
    using ParamCache = typename FluidSystem::template ParameterCache<Scalar>;
    ParamCache paramCache(1);
    BOOST_CHECK_EQUAL(paramCache.regionIndex(), 1);
    BOOST_CHECK_EQUAL(FluidSystem::numRegions(), 1);
    BOOST_CHECK_EQUAL(FluidSystem::numActivePhases(), 2);

    double GpuComputedVal = 0.0;
    double* gpuComputedValPtr = nullptr;
    OPM_GPU_SAFE_CALL(cudaMalloc(&gpuComputedValPtr, sizeof(double)));
    getReservoirTemperature<<<1, 1>>>(dynamicGpuFluidSystemView, gpuComputedValPtr);
    OPM_GPU_SAFE_CALL(cudaMemcpy(&GpuComputedVal, gpuComputedValPtr, sizeof(double), cudaMemcpyDeviceToHost));
    BOOST_CHECK_CLOSE(FluidSystem::reservoirTemperature(), GpuComputedVal, 1e-10);

    getReferenceDensityFromGasPvt<<<1, 1>>>(dynamicGpuFluidSystemView, gpuComputedValPtr);
    OPM_GPU_SAFE_CALL(cudaMemcpy(&GpuComputedVal, gpuComputedValPtr, sizeof(double), cudaMemcpyDeviceToHost));
    BOOST_CHECK_CLOSE(FluidSystem::gasPvt().gasReferenceDensity(0), GpuComputedVal, 1e-10);

    getReferenceDensity<<<1, 1>>>(dynamicGpuFluidSystemView, gpuComputedValPtr);
    OPM_GPU_SAFE_CALL(cudaMemcpy(&GpuComputedVal, gpuComputedValPtr, sizeof(double), cudaMemcpyDeviceToHost));
    BOOST_CHECK_CLOSE(FluidSystem::referenceDensity(0, 0), GpuComputedVal, 1e-10);


    OPM_GPU_SAFE_CALL(cudaFree(gpuComputedValPtr));
}
