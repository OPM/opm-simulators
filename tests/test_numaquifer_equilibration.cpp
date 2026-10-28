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
*/

#include "config.h"

#define BOOST_TEST_MODULE NumericalAquiferEquilibration

#include <boost/test/unit_test.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/OpmLog/StreamLog.hpp>

#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/start.hh>

#include <opm/simulators/flow/BlackoilModelParameters.hpp>
#include <opm/simulators/flow/FlowGenericVanguard.hpp>
#include <opm/simulators/flow/FlowProblemBlackoil.hpp>
#include <opm/simulators/flow/FlowProblemBlackoilProperties.hpp>
#include <opm/simulators/flow/NumericalAquiferEquilibrationCheck.hpp>
#include <opm/simulators/linalg/parallelbicgstabbackend.hh>
#include <opm/simulators/timestepping/AdaptiveTimeStepping.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

#include <cmath>
#include <cstddef>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

namespace Opm::Properties {

namespace TTag {
struct TestNumAquEquilTypeTag {
    using InheritsFrom = std::tuple<FlowBaseProblemBlackoil, BlackOilModel>;
};
} // namespace TTag

template<class TypeTag>
struct WellModel<TypeTag, TTag::TestNumAquEquilTypeTag>
{ using type = BlackoilWellModel<TypeTag>; };

template<class TypeTag>
struct EnableConvectiveMixing<TypeTag, TTag::TestNumAquEquilTypeTag>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::TestNumAquEquilTypeTag>
{ static constexpr bool value = false; };

} // namespace Opm::Properties

namespace {

using TypeTag = Opm::Properties::TTag::TestNumAquEquilTypeTag;
using Simulator = Opm::GetPropType<TypeTag, Opm::Properties::Simulator>;

std::unique_ptr<Simulator> initSimulator(const std::string& deck, const bool allowMismatch)
{
    const auto filenameArg = "--ecl-deck-file-name=" + deck;
    const auto allowArg = std::string{"--allow-inconsistent-aquifer-equilibration="}
        + (allowMismatch ? "true" : "false");

    const char* argv[] = {
        "test_numaquifer_equilibration",
        filenameArg.c_str(),
        allowArg.c_str(),
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

    Opm::FlowGenericVanguard::readDeck(deck, /*ignoredKeywords=*/"");

    return std::make_unique<Simulator>();
}

// Captures messages logged through OpmLog while in scope.
class LogCapture
{
public:
    LogCapture()
    {
        Opm::OpmLog::addBackend(name_, std::make_shared<Opm::StreamLog>(stream_, Opm::Log::DefaultMessageTypes));
    }

    ~LogCapture()
    {
        Opm::OpmLog::removeBackend(name_);
    }

    LogCapture(const LogCapture&) = delete;
    LogCapture& operator=(const LogCapture&) = delete;

    std::string text() const
    {
        return stream_.str();
    }

private:
    static constexpr const char* name_ = "NumAquEquil";
    std::ostringstream stream_;
};

unsigned elementIndex(const Simulator& simulator, const int i, const int j, const int k)
{
    const auto& vanguard = simulator.vanguard();
    const auto& dims = vanguard.cartesianDimensions();
    const auto cart = i + dims[0] * (j + dims[1] * k);
    for (unsigned elemIdx = 0; elemIdx < simulator.model().numGridDof(); ++elemIdx) {
        if (vanguard.cartesianIndex(elemIdx) == cart) {
            return elemIdx;
        }
    }
    throw std::invalid_argument("Cell is not active");
}

struct NumAquEquilFixture {
    NumAquEquilFixture()
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

BOOST_GLOBAL_FIXTURE(NumAquEquilFixture);

BOOST_AUTO_TEST_CASE(MismatchStopsSimulation)
{
    const auto simulator = initSimulator("numaquifer_equil_mismatch.DATA", /*allowMismatch=*/false);
    const LogCapture log;

    BOOST_CHECK_THROW(simulator->model().applyInitialSolution(), std::domain_error);

    // The 8 active region 1 cells are connected to the aquifer, with a
    // 10 bar difference between the equilibration regions.
    const auto text = log.text();
    BOOST_CHECK(text.find("Aquifer 1 (cells in EQLNUM region 2): 8 connection(s) "
                          "to cells in EQLNUM region 1") != std::string::npos);
    BOOST_CHECK(text.find("up to 10 BARS") != std::string::npos);
    BOOST_CHECK(text.find("AQUCON") != std::string::npos);
    BOOST_CHECK(text.find("'EQLNUM' 1") != std::string::npos);
    BOOST_CHECK(text.find("--allow-inconsistent-aquifer-equilibration=true") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(MismatchAllowedByOption)
{
    const auto simulator = initSimulator("numaquifer_equil_mismatch.DATA", /*allowMismatch=*/true);
    const LogCapture log;

    BOOST_CHECK_NO_THROW(simulator->model().applyInitialSolution());

    const auto text = log.text();
    BOOST_CHECK(text.find("Aquifer 1 (cells in EQLNUM region 2): 8 connection(s)") != std::string::npos);
    BOOST_CHECK(text.find("This option is set, so the simulation continues.") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(DifferenceBelowToleranceAccepted)
{
    const auto simulator = initSimulator("numaquifer_equil_tolerance.DATA", /*allowMismatch=*/false);
    const LogCapture log;

    BOOST_CHECK_NO_THROW(simulator->model().applyInitialSolution());
    BOOST_CHECK(log.text().find("Numerical aquifer(s)") == std::string::npos);
}

BOOST_AUTO_TEST_CASE(ImmobilePhasesAccepted)
{
    const auto simulator = initSimulator("numaquifer_equil_immobile.DATA", /*allowMismatch=*/false);
    const LogCapture log;

    BOOST_CHECK_NO_THROW(simulator->model().applyInitialSolution());
    BOOST_CHECK(log.text().find("Numerical aquifer(s)") == std::string::npos);

    // Not vacuous: the water potentials of the aquifer cell and a connected
    // region 1 cell, where the water is immobile, differ by more than the
    // tolerance.
    const auto& problem = simulator->problem();
    const auto aquIdx = elementIndex(*simulator, 0, 0, 4);
    const auto cellIdx = elementIndex(*simulator, 0, 1, 0);
    const auto& fsAqu = simulator->model().intensiveQuantities(aquIdx, 0).fluidState();
    const auto& fsCell = simulator->model().intensiveQuantities(cellIdx, 0).fluidState();

    constexpr auto water = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>::waterPhaseIdx;
    const double rho = 0.5 * (fsAqu.density(water).value() + fsCell.density(water).value());
    const double dPhi = (fsAqu.pressure(water).value() - fsCell.pressure(water).value())
        - rho * problem.gravity()[2] * (problem.dofCenterDepth(aquIdx) - problem.dofCenterDepth(cellIdx));
    BOOST_CHECK_GT(std::abs(dPhi), Opm::numericalAquiferEquilibrationTolerance);
}
