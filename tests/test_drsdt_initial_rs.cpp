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

#define BOOST_TEST_MODULE DrsdtInitialRs

#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if (BOOST_VERSION / 100000 == 1) && ((BOOST_VERSION / 100) % 1000 < 71)
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/start.hh>

#include <opm/simulators/flow/FlowProblemBlackoilProperties.hpp>
#include <opm/simulators/flow/BlackoilModelParameters.hpp>
#include <opm/simulators/timestepping/AdaptiveTimeStepping.hpp>
#include <opm/models/utils/parametersystem.hpp>
#include <opm/simulators/flow/FlowGenericVanguard.hpp>
#include <opm/simulators/flow/FlowProblemBlackoil.hpp>
#include <opm/simulators/linalg/parallelbicgstabbackend.hh>
#include <opm/simulators/wells/BlackoilWellModel.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

#include <cstddef>
#include <memory>
#include <string>

namespace Opm::Properties {

namespace TTag {
struct TestDrsdtTypeTag {
    using InheritsFrom = std::tuple<FlowBaseProblemBlackoil, BlackOilModel>;
};
} // namespace TTag

template<class TypeTag>
struct WellModel<TypeTag, TTag::TestDrsdtTypeTag>
{ using type = BlackoilWellModel<TypeTag>; };

template<class TypeTag>
struct EnableConvectiveMixing<TypeTag, TTag::TestDrsdtTypeTag>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::TestDrsdtTypeTag>
{ static constexpr bool value = false; };

} // namespace Opm::Properties

namespace {

template <class TypeTag>
std::unique_ptr<Opm::GetPropType<TypeTag, Opm::Properties::Simulator>>
initSimulator(const char* filename)
{
    using Simulator = Opm::GetPropType<TypeTag, Opm::Properties::Simulator>;

    const auto filenameArg = std::string {"--ecl-deck-file-name="} + filename;

    const char* argv[] = {
        "test_drsdt_initial_rs",
        filenameArg.c_str(),
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

    Opm::FlowGenericVanguard::readDeck(filename);

    return std::make_unique<Simulator>();
}

struct DrsdtFixture {
    using TypeTag = Opm::Properties::TTag::TestDrsdtTypeTag;

    DrsdtFixture()
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

BOOST_GLOBAL_FIXTURE(DrsdtFixture);

// The DRSDT limiter is relative to the previous step's Rs. At timeIdx 1 the
// problem reports exactly that stored value, so it must equal the initial Rs
// before any step has run -- otherwise the first step caps Rs at DRSDT * dt
// and strips an undersaturated reservoir of its dissolved gas.
BOOST_AUTO_TEST_CASE(MaxGasDissolutionFactorStartsFromInitialRs)
{
    using TypeTag = Opm::Properties::TTag::TestDrsdtTypeTag;

    const auto simulator = initSimulator<TypeTag>("drsdt_initial_rs.DATA");
    const auto& problem = simulator->problem();

    const std::size_t numDof = simulator->model().numGridDof();
    BOOST_REQUIRE_GT(numDof, 0u);

    for (std::size_t dofIdx = 0; dofIdx < numDof; ++dofIdx) {
        const auto initialRs = problem.initialFluidState(dofIdx).Rs();

        // The deck is undersaturated with a non-trivial RS; a zero here would
        // make the check vacuous.
        BOOST_REQUIRE_GT(initialRs, 0.0);

        BOOST_CHECK_CLOSE(problem.maxGasDissolutionFactor(/*timeIdx=*/1, dofIdx),
                          initialRs, 1.0e-10);
    }
}
