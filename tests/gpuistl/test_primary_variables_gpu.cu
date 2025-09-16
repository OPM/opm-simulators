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
#define HAVE_ECL_INPUT 1
#include <config.h>


#include <stdexcept>

#define BOOST_TEST_MODULE TestPrimaryVariablesGPU

#include <boost/test/unit_test.hpp>
#include <cuda.h>
#include <cuda_runtime.h>
#include <opm/common/ErrorMacros.hpp>
#include <opm/models/blackoil/blackoilmodel.hh>
#include <opm/models/blackoil/blackoilprimaryvariables.hh>
#include <opm/models/discretization/common/fvbaseprimaryvariables.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>
#include <opm/models/io/dgfvanguard.hh>
#include <opm/models/utils/start.hh>
#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/linalg/gpuistl/MiniVector.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>
#include <opm/simulators/linalg/gpuistl/gpu_smart_pointer.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Python/Python.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <cmath>
#include <type_traits>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystemNonStatic.hpp>
#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>
#include <opm/simulators/linalg/gpuistl/gpu_smart_pointer.hpp>


static constexpr const char* deckString1 = "-- =============== RUNSPEC\n"
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
using FluidSystem = Opm::BlackOilFluidSystem<double>;
using Evaluation = Opm::DenseAd::Evaluation<double, 2>;
using Scalar = typename Opm::MathToolbox<Evaluation>::Scalar;

using ScalarToUse = Opm::GetPropType<Opm::Properties::TTag::FlowProblem, Opm::Properties::Scalar>;

using FluidSystem = Opm::BlackOilFluidSystem<ScalarToUse>;
using Evaluation = Opm::DenseAd::Evaluation<double, 2>;
// using Scalar = typename Opm::MathToolbox<Evaluation>::Scalar;
using BlackOilFluidSystemView
    = Opm::BlackOilFluidSystemNonStatic<ScalarToUse, Opm::BlackOilDefaultFluidSystemIndices, Opm::gpuistl::GpuView>;

namespace Opm
{
namespace Properties
{
    namespace TTag
    {
        struct FlowProblemGPU {
            using InheritsFrom = std::tuple<FlowProblem>;
        };
    } // namespace TTag
    template <class TypeTag>
    struct FluidSystem<TypeTag, TTag::FlowProblemGPU> {
        using type = BlackOilFluidSystemView;
    };

} // namespace Properties

} // namespace Opm
using TypeTagCPU = Opm::Properties::TTag::FlowProblem;
using TypeTagGPU = Opm::Properties::TTag::FlowProblemGPU;
namespace
{
template <class FluidSystemView>
__global__ void
testCreationGPU(FluidSystemView fluidSystemView)
{
    // for now we just test that we can create a primaryvariables object
    Opm::BlackOilPrimaryVariables<TypeTagGPU, Opm::gpuistl::MiniVector> primaryVariables;
}
} // namespace

BOOST_AUTO_TEST_CASE(TestPrimaryVariablesCreationWithFieldVector)
{
    Opm::BlackOilPrimaryVariables<TypeTagCPU> primaryVariables;

    Opm::BlackOilPrimaryVariables<TypeTagCPU, Opm::gpuistl::MiniVector> primaryVariablesFieldVector;
}

BOOST_AUTO_TEST_CASE(TestPrimaryVariablesCrationGPU)
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
    testCreationGPU<<<1, 1>>>(dynamicGpuFluidSystemView);
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
    OPM_GPU_SAFE_CALL(cudaGetLastError());
}
