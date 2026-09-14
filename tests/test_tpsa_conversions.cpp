// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2026 NORCE AS

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

#define BOOST_TEST_MODULE TestTpsaConversions
#define BOOST_TEST_NO_MAIN

#include <opm/models/common/multiphasebaseparameters.hh>
#include <opm/models/discretization/common/fvbaseparameters.hh>
#include <opm/models/nonlinear/newtonmethodparams.hpp>
#include <opm/models/utils/basicparameters.hh>
#include <opm/models/utils/parametersystem.hpp>
#include <opm/simulators/flow/FlowGenericProblem_impl.hpp>
#include <opm/simulators/flow/FlowProblemParameters.hpp>
#include <opm/simulators/timestepping/EclTimeSteppingParams.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/LookUpData.hh>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Python/Python.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Units/UnitSystem.hpp>
#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <boost/test/unit_test.hpp>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/defaultgridview.hh>
#include <dune/grid/common/gridview.hh>

#include <memory>
#include <stdexcept>
#include <string>

#if HAVE_MPI
struct MPIError
{
    MPIError(std::string s, int e) : errorstring(std::move(s)), errorcode(e){}
    std::string errorstring;
    int errorcode;
};

void MPI_err_handler(MPI_Comm*, int* err_code, ...)
{
    std::vector<char> err_string(MPI_MAX_ERROR_STRING);
    int err_length;
    MPI_Error_string(*err_code, err_string.data(), &err_length);
    std::string s(err_string.data(), err_length);
    std::cerr << "An MPI Error ocurred:" << std::endl << s << std::endl;
    throw MPIError(s, *err_code);
}
#endif

bool init_unit_test_func()
{
    return true;
}

using Grid = Dune::CpGrid;
using GridView = Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>;
using FluidSystem = Opm::BlackOilFluidSystem<double, Opm::BlackOilDefaultFluidSystemIndices>;
using Problem = Opm::FlowGenericProblem<GridView, FluidSystem>;

static Opm::Deck createDeck(const std::string& gridProps)
{
    // Single-cell deck; gridProps supplies the mechanical/thermal
    // field-property keywords under test.
    std::string deck_string = R"(
RUNSPEC
DIMENS
1 1 1 /
EQLDIMS
/
TABDIMS
/
WATER
METRIC

GRID
DX
10.0 /
DY
10.0 /
DZ
10.0 /
TOPS
1000.0 /
PORO
0.3 /
PERMX
100.0 /
COPY
PERMX PERMY /
PERMX PERMZ /
/
)" + gridProps + R"(

PROPS
PVTW
1.0 /
ROCK
1.0132  1.0E-6 /
DENSITY
859.5 1033.0 0.854 /

SOLUTION
EQUIL
1000 100 /

SUMMARY

SCHEDULE
)";
    Opm::Parser parser;
    return parser.parseString(deck_string);
}

#define MAKE_PROBLEM(problemVar, gridProps)                                               \
    Opm::Deck deck = createDeck(gridProps);                                               \
    Opm::EclipseState eclState(deck);                                                     \
    auto python = std::make_shared<Opm::Python>();                                        \
    Opm::Schedule schedule(deck, eclState, python);                                       \
    Grid grid;                                                                            \
    grid.processEclipseFormat(&eclState.getInputGrid(), &eclState, false, false, false);  \
    auto gridView = grid.leafGridView();                                                  \
    Problem problemVar(eclState, schedule, gridView)

namespace {
// Unit conversions
const Opm::UnitSystem metricUnits(Opm::UnitSystem::UnitType::UNIT_TYPE_METRIC);
const double yModuleUnitFactor = metricUnits.to_si("Ymodule", 1.0);
const double pressureUnitFactor = metricUnits.to_si("Pressure", 1.0);

// Input parameters
constexpr double lameVal = 3.5;
constexpr double sModVal = 3.5;
constexpr double yModVal = 8.75;   // = 2*sModVal*(1+pRatioVal)
constexpr double pRatioVal = 0.25;

constexpr double biotCoeffVal = 0.9;
constexpr double poelCoefVal = 0.6;  // biotCoeffVal = poelCoefVal*(1-nu)/(1-2*nu)

constexpr double thermexrVal = 2.0;
constexpr double thelcoefVal = 4.0;

const double lameExpected = lameVal * yModuleUnitFactor;

// Expected THERMEXR and THELCOEF
const double biotTempThermexrExpected = thermexrVal * 17.5 * yModuleUnitFactor;
const double biotTempThelcoefExpected = thelcoefVal * pressureUnitFactor * 1.5;

constexpr double tol = 1.0e-8;
}

// ///
// lame() conversions
// ///

BOOST_AUTO_TEST_CASE(Lame_Direct)
{
    MAKE_PROBLEM(problem, "LAME\n" + std::to_string(lameVal) + " /\n");
    BOOST_CHECK_CLOSE(problem.lame(0), lameExpected, tol);
}

BOOST_AUTO_TEST_CASE(Lame_YmoduleSmodulus)
{
    MAKE_PROBLEM(problem,
                 "YMODULE\n" + std::to_string(yModVal) + " /\n" +
                 "SMODULUS\n" + std::to_string(sModVal) + " /\n");
    BOOST_CHECK_CLOSE(problem.lame(0), lameExpected, tol);
}

BOOST_AUTO_TEST_CASE(Lame_YmodulePratio)
{
    MAKE_PROBLEM(problem,
                 "YMODULE\n" + std::to_string(yModVal) + " /\n" +
                 "PRATIO\n" + std::to_string(pRatioVal) + " /\n");
    BOOST_CHECK_CLOSE(problem.lame(0), lameExpected, tol);
}

BOOST_AUTO_TEST_CASE(Lame_SmodulusPratio)
{
    MAKE_PROBLEM(problem,
                 "SMODULUS\n" + std::to_string(sModVal) + " /\n" +
                 "PRATIO\n" + std::to_string(pRatioVal) + " /\n");
    BOOST_CHECK_CLOSE(problem.lame(0), lameExpected, tol);
}

// ///
// biotCoeff() conversions
// ///

BOOST_AUTO_TEST_CASE(BiotCoeff_Direct)
{
    MAKE_PROBLEM(problem, "BIOTCOEF\n" + std::to_string(biotCoeffVal) + " /\n");
    BOOST_CHECK_CLOSE(problem.biotCoeff(0), biotCoeffVal, tol);
}

BOOST_AUTO_TEST_CASE(BiotCoeff_PoelcoefPratio)
{
    MAKE_PROBLEM(problem,
                 "POELCOEF\n" + std::to_string(poelCoefVal) + " /\n" +
                 "PRATIO\n" + std::to_string(pRatioVal) + " /\n");
    BOOST_CHECK_CLOSE(problem.biotCoeff(0), biotCoeffVal, tol);
}

// ///
// biotTemp() conversions: THERMEXR
// ///

BOOST_AUTO_TEST_CASE(BiotTemp_Thermexr_LameSmodulus)
{
    MAKE_PROBLEM(problem,
                 "THERMEXR\n" + std::to_string(thermexrVal) + " /\n" +
                 "LAME\n" + std::to_string(lameVal) + " /\n" +
                 "SMODULUS\n" + std::to_string(sModVal) + " /\n");
    BOOST_CHECK_CLOSE(problem.biotTemp(0), biotTempThermexrExpected, tol);
}

BOOST_AUTO_TEST_CASE(BiotTemp_Thermexr_YmodulePratio)
{
    MAKE_PROBLEM(problem,
                 "THERMEXR\n" + std::to_string(thermexrVal) + " /\n" +
                 "YMODULE\n" + std::to_string(yModVal) + " /\n" +
                 "PRATIO\n" + std::to_string(pRatioVal) + " /\n");
    BOOST_CHECK_CLOSE(problem.biotTemp(0), biotTempThermexrExpected, tol);
}

BOOST_AUTO_TEST_CASE(BiotTemp_Thermexr_LamePratio)
{
    MAKE_PROBLEM(problem,
                 "THERMEXR\n" + std::to_string(thermexrVal) + " /\n" +
                 "LAME\n" + std::to_string(lameVal) + " /\n" +
                 "PRATIO\n" + std::to_string(pRatioVal) + " /\n");
    BOOST_CHECK_CLOSE(problem.biotTemp(0), biotTempThermexrExpected, tol);
}

BOOST_AUTO_TEST_CASE(BiotTemp_Thermexr_SmodulusPratio)
{
    MAKE_PROBLEM(problem,
                 "THERMEXR\n" + std::to_string(thermexrVal) + " /\n" +
                 "SMODULUS\n" + std::to_string(sModVal) + " /\n" +
                 "PRATIO\n" + std::to_string(pRatioVal) + " /\n");
    BOOST_CHECK_CLOSE(problem.biotTemp(0), biotTempThermexrExpected, tol);
}

BOOST_AUTO_TEST_CASE(BiotTemp_Thermexr_YmoduleSmodulus)
{
    MAKE_PROBLEM(problem,
                 "THERMEXR\n" + std::to_string(thermexrVal) + " /\n" +
                 "YMODULE\n" + std::to_string(yModVal) + " /\n" +
                 "SMODULUS\n" + std::to_string(sModVal) + " /\n");
    BOOST_CHECK_CLOSE(problem.biotTemp(0), biotTempThermexrExpected, tol);
}

// ///
// biotTemp() conversions: THELCOEF
// ///

BOOST_AUTO_TEST_CASE(BiotTemp_Thelcoef_LameSmodulus)
{
    MAKE_PROBLEM(problem,
                 "THELCOEF\n" + std::to_string(thelcoefVal) + " /\n" +
                 "LAME\n" + std::to_string(lameVal) + " /\n" +
                 "SMODULUS\n" + std::to_string(sModVal) + " /\n");
    BOOST_CHECK_CLOSE(problem.biotTemp(0), biotTempThelcoefExpected, tol);
}

BOOST_AUTO_TEST_CASE(BiotTemp_Thelcoef_PratioOnly)
{
    MAKE_PROBLEM(problem,
                 "THELCOEF\n" + std::to_string(thelcoefVal) + " /\n" +
                 "PRATIO\n" + std::to_string(pRatioVal) + " /\n");
    BOOST_CHECK_CLOSE(problem.biotTemp(0), biotTempThelcoefExpected, tol);
}

BOOST_AUTO_TEST_CASE(BiotTemp_Thelcoef_YmoduleSmodulus)
{
    MAKE_PROBLEM(problem,
                 "THELCOEF\n" + std::to_string(thelcoefVal) + " /\n" +
                 "YMODULE\n" + std::to_string(yModVal) + " /\n" +
                 "SMODULUS\n" + std::to_string(sModVal) + " /\n");
    BOOST_CHECK_CLOSE(problem.biotTemp(0), biotTempThelcoefExpected, tol);
}

// ///
// biotTemp() error handling
// ///

BOOST_AUTO_TEST_CASE(BiotTemp_Thermexr_InsufficientParams_Throws)
{
    // THERMEXR alone with only PRATIO available does not match any of the
    // supported (LAME,SMODULUS) / (YMODULE,PRATIO) / (LAME,PRATIO) /
    // (SMODULUS,PRATIO) / (YMODULE,SMODULUS) combinations.
    MAKE_PROBLEM(problem,
                 "THERMEXR\n" + std::to_string(thermexrVal) + " /\n" +
                 "PRATIO\n" + std::to_string(pRatioVal) + " /\n");
    BOOST_CHECK_THROW(problem.biotTemp(0), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(BiotTemp_BothKeywordsGiven_Throws)
{
    MAKE_PROBLEM(problem,
                 "THERMEXR\n" + std::to_string(thermexrVal) + " /\n" +
                 "THELCOEF\n" + std::to_string(thelcoefVal) + " /\n" +
                 "LAME\n" + std::to_string(lameVal) + " /\n" +
                 "SMODULUS\n" + std::to_string(sModVal) + " /\n");
    BOOST_CHECK_THROW(problem.biotTemp(0), std::runtime_error);
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
#if HAVE_MPI
    // register a throwing error handler to allow for
    // debugging with "catch throw" in gdb
    MPI_Errhandler handler;
    MPI_Comm_create_errhandler(MPI_err_handler, &handler);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, handler);
#endif

    // FlowProblemParameters.cpp's registerFlowProblemParameters() calls
    // SetDefault() on several parameters that are normally registered by
    // Simulator<TypeTag>::registerParameters() (which we don't instantiate
    // here, to avoid pulling in the whole type-tag/property system) -
    // register them manually instead.
    Opm::Parameters::Register<Opm::Parameters::EndTime<double>>
        ("The simulation time at which the simulation is finished [s]");
    Opm::Parameters::Register<Opm::Parameters::InitialTimeStepSize<double>>
        ("The size of the initial time step [s]");
    Opm::Parameters::Register<Opm::Parameters::EnableVtkOutput>
        ("Write VTK output files");
    Opm::Parameters::Register<Opm::Parameters::EnableIntensiveQuantityCache>
        ("Cache intensive quantities");
    Opm::Parameters::Register<Opm::Parameters::EnableStorageCache>
        ("Cache storage terms");
    Opm::Parameters::Register<Opm::Parameters::NewtonTolerance<double>>
        ("Newton convergence tolerance");
    Opm::Parameters::Register<Opm::Parameters::EnableGravity>
        ("Enable gravity");

    Opm::registerFlowProblemParameters<double>();
    Opm::registerEclTimeSteppingParameters<double>();
    Opm::Parameters::endRegistration();

    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
