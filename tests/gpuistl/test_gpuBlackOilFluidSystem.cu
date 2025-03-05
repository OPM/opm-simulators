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
 * \brief This is the unit test for the black oil fluid system
 *
 * This test requires the presence of opm-parser.
 */
#include "config.h"

#if !HAVE_ECL_INPUT
#error "The test for the black oil fluid system classes requires ecl input support in opm-common"
#endif

#include <boost/mpl/list.hpp>

#define BOOST_TEST_MODULE EclBlackOilFluidSystem
#include <boost/test/unit_test.hpp>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystemNonStatic.hpp>
#include <opm/material/fluidstates/BlackOilFluidState.hpp>
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

using Types = boost::mpl::list<double,Opm::DenseAd::Evaluation<double,2>>;
// using GpuB = ::Opm::gpuistl::GpuBuffer;
// using GpuPointer = ::Opm::gpuistl::ViewPointer;

BOOST_AUTO_TEST_CASE_TEMPLATE(BlackOil, Evaluation, Types)
{
    // test the black-oil specific methods of BlackOilFluidSystem. The generic methods
    // for fluid systems are already tested by the generic test for all fluidsystems.

    using Scalar = typename Opm::MathToolbox<Evaluation>::Scalar;
    using FluidSystem = Opm::BlackOilFluidSystem<double>;

    using GpuB = Opm::gpuistl::GpuBuffer<double>;
    using GpuBufCo2Tables = Opm::CO2Tables<double, GpuB>;
    using GpuBufBrineCo2Pvt = Opm::BrineCo2Pvt<double, GpuBufCo2Tables, GpuB>;

    static constexpr int numPhases = FluidSystem::numPhases;

    static constexpr int gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr int oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static constexpr int waterPhaseIdx = FluidSystem::waterPhaseIdx;

    static constexpr int gasCompIdx = FluidSystem::gasCompIdx;
    static constexpr int oilCompIdx = FluidSystem::oilCompIdx;
    static constexpr int waterCompIdx = FluidSystem::waterCompIdx;

    Opm::Parser parser;

    auto deck = parser.parseString(deckString1);
    auto python = std::make_shared<Opm::Python>();
    Opm::EclipseState eclState(deck);
    Opm::Schedule schedule(deck, eclState, python);

    FluidSystem::initFromState(eclState, schedule);


    // TODO: get the type of the fluid system gasvpt multiplexer approach

    auto& dynamicFluidSystem = FluidSystem::getNonStaticInstance();

    using GpuCo2Tables = Opm::CO2Tables<double, ::Opm::gpuistl::GpuBuffer<double>>;
    using CpuCo2Tables = Opm::CO2Tables<double>;
    using GpuGasPvt = Opm::Co2GasPvt<double, GpuCo2Tables, ::Opm::gpuistl::GpuBuffer<double>>;
    using GpuWaterPvt = Opm::BrineCo2Pvt<double, GpuCo2Tables, ::Opm::gpuistl::GpuBuffer<double>>;

    // auto cpuGasPvt = dynamicFluidSystem.gasPvt().realGasPvt();

    auto gpuCo2GasPvt = ::Opm::gpuistl::copy_to_gpu<::Opm::gpuistl::GpuBuffer, double>(dynamicFluidSystem);
    // auto gpuGasPvt = ::Opm::gpuistl::copy_to_gpu<::Opm::gpuistl::GpuBuffer, CO2Tables<double, ::Opm::gpuistl::GpuBuffer<double>>(dynamicFluidSystem->getPvt());

    // auto dynamicGpuFluidSystem = ::Opm::gpuistl::copy_to_gpu<::Opm::gpuistl::GpuBuffer>(dynamicFluidSystem, gpuGasPvt);
    // auto dynamicGpuFluidSystemView = ::Opm::gpuistl::make_view<::Opm::gpuistl::GpuView, ::Opm::gpuistl::ViewPointer>(dynamicGpuFluidSystem);

    // create a parameter cache
    // using ParamCache = typename FluidSystem::template ParameterCache<Scalar>;
    // ParamCache paramCache(/*maxOilSat=*/0.5, /*regionIdx=*/1);
    // BOOST_CHECK_EQUAL(paramCache.regionIndex(), 1);

    // BOOST_CHECK_SMALL(Opm::abs(FluidSystem::reservoirTemperature() - (273.15 + 15.555)), 1e-10);
    // BOOST_CHECK_EQUAL(FluidSystem::numRegions(), 1);
    // BOOST_CHECK_EQUAL(FluidSystem::numActivePhases(), 2);

    // BOOST_CHECK(FluidSystem::phaseIsActive(0));
    // BOOST_CHECK(FluidSystem::phaseIsActive(2));

    // // make sure that the {oil,gas,water}Pvt() methods are available
    // [[maybe_unused]] const auto& gPvt = FluidSystem::gasPvt();
    // [[maybe_unused]] const auto& oPvt = FluidSystem::oilPvt();
    // [[maybe_unused]] const auto& wPvt = FluidSystem::waterPvt();
}
