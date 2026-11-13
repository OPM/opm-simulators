// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2026 Equinor ASA.

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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/

#include "config.h"

#define BOOST_TEST_MODULE FlowThresholdPressure

#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if (BOOST_VERSION / 100000 == 1) && ((BOOST_VERSION / 100) % 1000 < 71)
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif

#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/start.hh>

#include <opm/simulators/flow/BlackoilModelParameters.hpp>
#include <opm/simulators/flow/FlowGenericVanguard.hpp>
#include <opm/simulators/flow/FlowProblemBlackoil.hpp>
#include <opm/simulators/flow/FlowProblemBlackoilProperties.hpp>
#include <opm/simulators/linalg/parallelbicgstabbackend.hh>
#include <opm/simulators/timestepping/AdaptiveTimeStepping.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

#include <memory>
#include <stdexcept>
#include <string>

namespace Opm::Properties {

namespace TTag {
struct TestThresholdPressureTypeTag {
    using InheritsFrom = std::tuple<FlowBaseProblemBlackoil, BlackOilModel>;
};
} // namespace TTag

template<class TypeTag>
struct WellModel<TypeTag, TTag::TestThresholdPressureTypeTag>
{ using type = BlackoilWellModel<TypeTag>; };

template<class TypeTag>
struct EnableConvectiveMixing<TypeTag, TTag::TestThresholdPressureTypeTag>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::TestThresholdPressureTypeTag>
{ static constexpr bool value = false; };

} // namespace Opm::Properties

namespace {

using TypeTag = Opm::Properties::TTag::TestThresholdPressureTypeTag;
using Simulator = Opm::GetPropType<TypeTag, Opm::Properties::Simulator>;

constexpr const char* deckName = "thpres_default_phases.DATA";

// Gravity corrected potential differences between the two cells of each
// column of the test deck: incompressible oil (800 kg/m3) and water
// (1000 kg/m3), zero capillary pressure, cell centres 10 m apart, and a
// bottom cell pressure exceeding the top cell pressure by 0.834532 bar.
constexpr double gravity = 9.80665;
constexpr double dz = 10.0;
constexpr double dp = 0.834532e5;
constexpr double dPhiOil = dp - 800.0 * gravity * dz;    //   5000 Pa
constexpr double dPhiWater = dp - 1000.0 * gravity * dz; // -14613.3 Pa

std::unique_ptr<Simulator> initSimulator(const std::string& phases)
{
    const auto filenameArg = std::string{"--ecl-deck-file-name="} + deckName;
    const auto phasesArg = "--threshold-pressure-default-phases=" + phases;

    const char* argv[] = {
        "test_flow_thresholdpressure",
        filenameArg.c_str(),
        phasesArg.c_str(),
        "--check-satfunc-consistency=false",
        "--enable-ecl-output=false",
        "--enable-vtk-output=false",
    };

    Opm::setupParameters_<TypeTag>(/*argc=*/sizeof(argv) / sizeof(argv[0]),
                                   argv,
                                   /*registerParams=*/false,
                                   /*allowUnused=*/false,
                                   /*handleHelp=*/true,
                                   /*myRank=*/0);

    Opm::FlowGenericVanguard::readDeck(deckName, /*ignoredKeywords=*/"");

    auto simulator = std::make_unique<Simulator>();

    // Defaulted threshold pressures are computed once the initial solution
    // has been applied.
    simulator->model().applyInitialSolution();

    return simulator;
}

unsigned elementIndex(const Simulator& simulator, const int i, const int j, const int k)
{
    const auto& vanguard = simulator.vanguard();
    const unsigned cart = vanguard.cartesianIndex({i, j, k});
    for (unsigned elemIdx = 0; elemIdx < simulator.model().numGridDof(); ++elemIdx) {
        if (vanguard.cartesianIndex(elemIdx) == cart) {
            return elemIdx;
        }
    }
    throw std::invalid_argument("Cell is not active");
}

// Threshold pressure between the top and bottom cell of column i (zero based).
double columnThreshold(const Simulator& simulator, const int i)
{
    const auto top = elementIndex(simulator, i, 0, 0);
    const auto bottom = elementIndex(simulator, i, 0, 1);
    const auto& problem = simulator.problem();

    const double thpres = problem.thresholdPressure(top, bottom);
    BOOST_CHECK_EQUAL(thpres, problem.thresholdPressure(bottom, top));

    return thpres;
}

struct ThresholdPressureFixture {
    ThresholdPressureFixture()
    {
        int argc = boost::unit_test::framework::master_test_suite().argc;
        char** argv = boost::unit_test::framework::master_test_suite().argv;
#if HAVE_DUNE_FEM
        Dune::Fem::MPIManager::initialize(argc, argv);
#else
        Dune::MPIHelper::instance(argc, argv);
#endif
        using namespace Opm;
        FlowGenericVanguard::setCommunication(std::make_unique<Opm::Parallel::Communication>());
        ThreadManager::registerParameters();
        // Registered in opm/models/nonlinear/newtonmethodparams.cpp, not reached from here.
        Parameters::Register<Parameters::NewtonMaxIterations>("The maximum number of Newton iterations per time step");
        BlackoilModelParameters<double>::registerParameters();
        AdaptiveTimeStepping<TypeTag>::registerParameters();
        Parameters::Register<Parameters::EnableTerminalOutput>("Dummy added for the well model to compile.");
        registerAllParameters_<TypeTag>(true);
    }
};

} // Anonymous namespace

BOOST_GLOBAL_FIXTURE(ThresholdPressureFixture);

// Default: only phases which are mobile in the upstream cell contribute.  Oil
// is mobile everywhere and water is immobile everywhere, so all thresholds
// are given by the oil potential difference.
BOOST_AUTO_TEST_CASE(MobilePhasesOnly)
{
    const auto simulator = initSimulator("mobile");

    for (const int i : {0, 2, 4, 6, 8}) {
        BOOST_TEST_CONTEXT("Column I=" << i + 1) {
            BOOST_CHECK_CLOSE(columnThreshold(*simulator, i), dPhiOil, 1.0e-6);
        }
    }
}

BOOST_AUTO_TEST_CASE(ImmobilePhaseAboveConnate)
{
    const auto simulator = initSimulator("above-connate");

    // Upstream water saturation 0.25 above the connate saturation 0.2: water
    // contributes although immobile.  The water potential difference is
    // dominated by the gravity term.
    BOOST_CHECK_CLOSE(columnThreshold(*simulator, 0), -dPhiWater, 1.0e-6);
}

BOOST_AUTO_TEST_CASE(PhaseAtConnateExcluded)
{
    const auto simulator = initSimulator("above-connate");

    // Upstream water saturation equal to the connate saturation: only oil
    // contributes.
    BOOST_CHECK_CLOSE(columnThreshold(*simulator, 2), dPhiOil, 1.0e-6);
}

BOOST_AUTO_TEST_CASE(ConnateOfUpstreamCellSaturationRegion)
{
    const auto simulator = initSimulator("above-connate");

    // Upstream (top) cell in SATNUM 1 with water saturation 0.3 above its
    // connate saturation 0.2.  The downstream cell (SATNUM 2) is at its
    // connate saturation 0.4.
    BOOST_CHECK_CLOSE(columnThreshold(*simulator, 4), -dPhiWater, 1.0e-6);

    // Upstream (top) cell in SATNUM 2 at its connate saturation 0.4, which is
    // above the connate saturation of SATNUM 1.  The downstream cell
    // (SATNUM 1) is above its connate saturation.
    BOOST_CHECK_CLOSE(columnThreshold(*simulator, 6), dPhiOil, 1.0e-6);
}

BOOST_AUTO_TEST_CASE(UnscaledConnateSaturation)
{
    const auto simulator = initSimulator("above-connate");

    // Upstream water saturation equal to the cell's scaled connate
    // saturation (SWL 0.25), but above the connate saturation of the
    // unscaled table (0.2).
    BOOST_CHECK_CLOSE(columnThreshold(*simulator, 8), -dPhiWater, 1.0e-6);
}
