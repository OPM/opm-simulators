/*
  Copyright 2026 SINTEF AS

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

#undef HAVE_DUNE_ALUGRID
#define HAVE_DUNE_ALUGRID 0
#undef HAVE_DUNE_FEM
#define HAVE_DUNE_FEM 0

#define BOOST_TEST_MODULE GpuAssemblyBridge
#include <boost/test/unit_test.hpp>

#include <flow/flow_gpu.hpp>
#include <opm/common/CriticalError.hpp>
#include <opm/material/common/ResetLocale.hpp>
#include <opm/models/utils/start.hh>
#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/flow/SimulatorFullyImplicit.hpp>
#include <opm/simulators/linalg/gpuistl/ISTLSolverGPUISTL.hpp>
#include <opm/simulators/utils/moduleVersion.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

namespace Opm {
std::string moduleVersionName() { return "test-version"; }

template <class TypeTag>
class BridgeTestProblem : public FlowProblemBlackoil<TypeTag>
{
    using Parent = FlowProblemBlackoil<TypeTag>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
public:
    using Parent::Parent;
    double sourceValue = 0.0;
    double sourceDerivative = 0.0;
    double boundaryTemperatureOffset = 0.0;

    void addToSourceDense(RateVector& rate, unsigned cell, unsigned) const override
    {
        if (cell == 0) {
            Evaluation source(sourceValue);
            source.setDerivative(0, sourceDerivative);
            rate[0] += source;
        }
    }

    auto boundaryFluidState(unsigned cell, int direction) const
    {
        auto state = Parent::boundaryFluidState(cell, direction);
        if (cell == 0 && direction == 0) {
            state.setTemperature(state.temperature(0) + boundaryTemperatureOffset);
        }
        return state;
    }

};

namespace Properties {
namespace TTag {
struct AssemblyBridgeTest { using InheritsFrom = std::tuple<FlowGasWaterEnergyProblemGPU>; };
template <template <class> class Storage>
struct to_gpu_type<AssemblyBridgeTest, Storage> {
    using type = FlowGasWaterEnergyProblemGPUTrue<Storage>;
};
}
template <class TypeTag>
struct Problem<TypeTag, TTag::AssemblyBridgeTest> { using type = BridgeTestProblem<TypeTag>; };
}
}

namespace {
using TypeTag = Opm::Properties::TTag::AssemblyBridgeTest;
using Simulator = Opm::GetPropType<TypeTag, Opm::Properties::Simulator>;
using Matrix = typename Opm::GetPropType<TypeTag, Opm::Properties::SparseMatrixAdapter>::IstlMatrix;
using Vector = Opm::GetPropType<TypeTag, Opm::Properties::GlobalEqVector>;
using RateVector = Opm::GetPropType<TypeTag, Opm::Properties::RateVector>;
using Indices = Opm::GetPropType<TypeTag, Opm::Properties::Indices>;
using LocalResidual = Opm::GetPropType<TypeTag, Opm::Properties::LocalResidual>;
using CpuLinearizer = Opm::TpfaLinearizer<Opm::Properties::TTag::FlowGasWaterEnergyProblem>;
constexpr unsigned numEq = Opm::getPropValue<TypeTag, Opm::Properties::NumEq>();

// Self-contained, deliberately small; no opm-tests checkout or external deck.
const char* deck = R"(
RUNSPEC
DIMENS
3 1 1 /
WATER
GAS
CO2STORE
THERMAL
METRIC
GRID
FAULTS
'BRIDGE' 1 1 1 1 1 1 X /
/
DX
3*1 /
DY
3*1 /
DZ
3*1 /
TOPS
3*0 /
PERMX
3*1013.25 /
PORO
3*0.25 /
COPY
PERMX PERMY /
PERMX PERMZ /
/
THCONR
3*1.5 /
BCCON
1 1 1 1 1 1 1 X- /
/
PROPS
SGWFN
0 0 1 0.0306
1 1 0 0.306 /
SPECROCK
20 2125 100 2125 /
SOLUTION
EQUIL
0 300 100 0 0 0 1 1 0 /
RTEMPVD
0 36.12
1000 70 /
SCHEDULE
BCPROP
1 THERMAL /
/
TSTEP
1 /
END
)";

struct Fixture {
    std::filesystem::path directory;
    std::unique_ptr<Simulator> simulator;

    Fixture() : Fixture("", "dilu") {}

    explicit Fixture(const std::string& laterControls) : Fixture(laterControls, "dilu") {}

    Fixture(const std::string& laterControls, const std::string& solverName)
    {
        Opm::resetLocale();
        int argc = boost::unit_test::framework::master_test_suite().argc;
        char** argv = boost::unit_test::framework::master_test_suite().argv;
        Dune::MPIHelper::instance(argc, argv);
        Opm::FlowGenericVanguard::setCommunication(std::make_unique<Opm::Parallel::Communication>());
        static bool registered = false;
        if (!registered) {
            Opm::ThreadManager::registerParameters();
            Opm::NewtonMethodParams<double>::registerParameters();
            Opm::BlackoilModelParameters<double>::registerParameters();
            Opm::AdaptiveTimeStepping<TypeTag>::registerParameters();
            Opm::Parameters::Register<Opm::Parameters::EnableTerminalOutput>("Enable terminal output.");
            Opm::registerAllParameters_<TypeTag>(false);
            Opm::Parameters::endRegistration();
            registered = true;
        }
        directory = std::filesystem::current_path() / "gpu-assembly-bridge-test";
        std::filesystem::create_directories(directory);
        const auto filename = directory / "BRIDGE.DATA";
        std::string input(deck);
        if (!laterControls.empty()) {
            input.replace(input.find("TSTEP\n1 /\nEND"), std::string::npos,
                "TSTEP\n1 /\n" + laterControls + "TSTEP\n1 /\nEND\n");
        }
        { std::ofstream out(filename); out << input; }
        const auto deckArg = "--ecl-deck-file-name=" + filename.string();
        const auto outputArg = "--output-dir=" + directory.string();
        const auto solverArg = "--linear-solver=" + solverName;
        const char* args[] = {"test_GpuAssemblyBridge", deckArg.c_str(), outputArg.c_str(),
            solverArg.c_str(), "--linear-solver-accelerator=gpu",
            "--matrix-add-well-contributions=true", "--enable-storage-cache=false",
            "--linear-solver-reduction=1e-10", "--threads-per-process=2",
            "--enable-tuning=false", "--check-satfunc-consistency=false"};
        Opm::setupParameters_<TypeTag>(sizeof(args) / sizeof(args[0]), args, false, true, false, 0);
        Opm::FlowGenericVanguard::readDeck(filename.string());
        simulator = std::make_unique<Simulator>();
        simulator->model().applyInitialSolution();
        simulator->setEpisodeIndex(0);
        simulator->setTimeStepSize(100.0);
        simulator->model().invalidateAndUpdateIntensiveQuantities(0);
        simulator->model().invalidateAndUpdateIntensiveQuantities(1);
    }

    ~Fixture()
    {
        simulator.reset();
        std::error_code error;
        std::filesystem::remove_all(directory, error);
    }
};

// Assemble on the CPU into independent storage. This never copies the device
// Jacobian and does not use the persistent GPU assembly state as its reference.
std::pair<Matrix, Vector> assembleCpu(Simulator& sim)
{
    auto& linearizer = sim.model().linearizer();
    Matrix matrix = linearizer.jacobian().istlMatrix();
    matrix = 0.0;
    Vector residual(sim.model().numGridDof());
    residual = 0.0;
    auto neighbors = linearizer.getNeighborInfo();
    std::vector<typename Matrix::block_type*> diagonal(matrix.N());
    Opm::FullDomain<> domain;
    domain.cells.resize(matrix.N());
    std::iota(domain.cells.begin(), domain.cells.end(), 0);
    for (unsigned cell = 0; cell < matrix.N(); ++cell) {
        diagonal[cell] = &matrix[cell][cell];
        for (auto& neighbor : neighbors[cell]) {
            neighbor.matBlockAddress = &matrix[neighbor.neighbor][cell];
        }
    }
    auto velocity = linearizer.getVelocityInfo();
    for (unsigned cell = 0; cell < matrix.N(); ++cell) {
        CpuLinearizer::template linearize_cell<false, LocalResidual>(cell, domain, neighbors,
            diagonal, residual, sim.model(), velocity, sim.timeStepSize(), sim.problem());
        RateVector source(0.0);
        LocalResidual::computeSource(source, sim.problem(), sim.model().intensiveQuantities(cell, 0), cell, 0);
        source *= -sim.model().dofTotalVolume(cell);
        typename Vector::block_type r(0.0);
        typename Matrix::block_type j(0.0);
        CpuLinearizer::setResAndJacobi(r, j, source);
        residual[cell] += r;
        matrix[cell][cell] += j;
    }
    using IQ = Opm::GetPropType<TypeTag, Opm::Properties::IntensiveQuantities>;
    using Boundary = Opm::BoundaryConditionData<typename Vector::block_type, typename IQ::ScalarFluidState>;
    const auto [type, rate] = sim.problem().boundaryCondition(0, 0);
    typename Vector::block_type massRate(0.0);
    for (unsigned eq = 0; eq < numEq; ++eq) {
        massRate[eq] = rate[eq].value();
    }
    const auto& fluid = sim.problem().boundaryFluidState(0, 0);
    Boundary boundary{type, massRate, fluid.pvtRegionIndex(), 0, 1.0,
                      sim.problem().dofCenterDepth(0), fluid};
    RateVector flux(0.0);
    LocalResidual::computeBoundaryFlux(flux, sim.problem(), boundary,
                                       sim.model().intensiveQuantities(0, 0), 0);
    typename Vector::block_type r(0.0);
    typename Matrix::block_type j(0.0);
    CpuLinearizer::setResAndJacobi(r, j, flux);
    residual[0] += r;
    matrix[0][0] += j;
    return {std::move(matrix), std::move(residual)};
}

void checkClose(double actual, double expected)
{
    BOOST_CHECK_SMALL(actual - expected, 1e-10 + 1e-10 * std::abs(expected));
}

void checkAssembly(Simulator& sim, Opm::TpfaLinearizer<TypeTag>& linearizer)
{
    BOOST_REQUIRE(sim.problem().nonTrivialBoundaryConditions());
    linearizer.linearizeDomain();
    BOOST_REQUIRE(linearizer.gpuSystemActive());
    const auto [matrix, residual] = assembleCpu(sim);
    const auto values = linearizer.gpuJacobian().getNonZeroValues().asStdVector();
    std::size_t offset = 0;
    for (auto row = matrix.begin(); row != matrix.end(); ++row) {
        for (auto block = row->begin(); block != row->end(); ++block) {
            for (unsigned eq = 0; eq < numEq; ++eq) {
                for (unsigned pv = 0; pv < numEq; ++pv) {
                    checkClose(values.at(offset++), (*block)[eq][pv]);
                }
            }
        }
    }
    BOOST_CHECK_EQUAL(offset, values.size());
    const auto flattened = linearizer.flattenedGpuResidual().asStdVector();
    for (unsigned cell = 0; cell < residual.size(); ++cell) {
        for (unsigned eq = 0; eq < numEq; ++eq) {
            checkClose(linearizer.residual()[cell][eq], residual[cell][eq]);
            checkClose(flattened[cell * numEq + eq], residual[cell][eq]);
        }
    }
}

void checkAssembly(Simulator& sim)
{
    checkAssembly(sim, sim.model().linearizer());
}

void checkSolution(Simulator& sim, const Vector& x)
{
    const auto [matrix, residual] = assembleCpu(sim);
    Vector error = residual;
    matrix.mmv(x, error);
    BOOST_REQUIRE_GT(residual.two_norm(), 1e-6);
    BOOST_CHECK_SMALL(error.two_norm() / residual.two_norm(), 1e-7);
}
}

BOOST_FIXTURE_TEST_CASE(PersistentAssemblyAndResidentSolve, Fixture)
{
    auto& sim = *simulator;
    auto& linearizer = sim.model().linearizer();
    sim.problem().sourceValue = 0.25;
    sim.problem().sourceDerivative = 0.125;
    checkAssembly(sim);

    // Refresh both IQ time levels, boundary data and a derivative-only source.
    sim.model().solution(0)[1][Indices::pressureSwitchIdx] *= 1.001;
    sim.model().solution(1)[2][Indices::pressureSwitchIdx] *= 0.999;
    sim.model().invalidateAndUpdateIntensiveQuantities(0);
    sim.model().invalidateAndUpdateIntensiveQuantities(1);
    sim.problem().sourceValue = 0.0;
    const auto beforeBoundary = assembleCpu(sim).second;
    sim.problem().boundaryTemperatureOffset = 1000.0;
    linearizer.updateBoundaryConditionData();
    const auto afterBoundary = assembleCpu(sim).second;
    BOOST_REQUIRE_GT(std::abs(afterBoundary[0][Indices::contiEnergyEqIdx]
                             - beforeBoundary[0][Indices::contiEnergyEqIdx]), 1e-8);
    checkAssembly(sim);

    sim.problem().sourceDerivative = 0.0;
    checkAssembly(sim);
    Opm::Parameters::SetDefault<Opm::Parameters::SeparateSparseSourceTerms>(true);
    {
        Opm::TpfaLinearizer<TypeTag> sparse;
        sparse.init(sim);
        checkAssembly(sim, sparse);
    }
    Opm::Parameters::SetDefault<Opm::Parameters::SeparateSparseSourceTerms>(false);
    auto& solver = sim.model().newtonMethod().linearSolver();
    for (unsigned iteration = 0; iteration < 2; ++iteration) {
        BOOST_REQUIRE(solver.prepareGpu(linearizer.gpuJacobian(), linearizer.flattenedGpuResidual()));
        Vector x(sim.model().numGridDof());
        BOOST_REQUIRE(solver.solveGpu(x));
        checkSolution(sim, x);
        sim.problem().sourceValue += 0.5;
        checkAssembly(sim);
    }
    // Exercise the backend's owned/borrowed transition and explicit invalidation.
    auto [matrix, residual] = assembleCpu(sim);
    solver.prepare(matrix, residual);
    Vector x(residual.size());
    x = 0.0;
    BOOST_REQUIRE(solver.solve(x));
    checkSolution(sim, x);
    sim.problem().sourceValue += 0.5;
    checkAssembly(sim);
    BOOST_REQUIRE(solver.prepareGpu(linearizer.gpuJacobian(), linearizer.flattenedGpuResidual()));
    BOOST_REQUIRE(solver.solveGpu(x));
    checkSolution(sim, x);
    sim.problem().sourceValue += 0.5;
    Opm::Parameters::SetDefault<Opm::Parameters::EnableGpuAssemblySolverBridge>(false);
    linearizer.linearizeDomain();
    BOOST_CHECK(!linearizer.gpuSystemActive());
    const auto hostReference = assembleCpu(sim);
    for (unsigned cell = 0; cell < residual.size(); ++cell) {
        for (unsigned eq = 0; eq < numEq; ++eq) {
            checkClose(linearizer.residual()[cell][eq], hostReference.second[cell][eq]);
        }
    }
    solver.prepare(linearizer.jacobian(), linearizer.residual());
    x = 0.0;
    BOOST_REQUIRE(solver.solve(x));
    checkSolution(sim, x);
    Opm::Parameters::SetDefault<Opm::Parameters::EnableGpuAssemblySolverBridge>(true);
    sim.problem().sourceValue += 0.5;
    checkAssembly(sim);
    solver.eraseMatrix();
    BOOST_REQUIRE(solver.prepareGpu(linearizer.gpuJacobian(), linearizer.flattenedGpuResidual()));
    BOOST_REQUIRE(solver.solveGpu(x));
    checkSolution(sim, x);

    Opm::gpuistl::ISTLSolverGPUISTL<TypeTag> backend(sim);
    Opm::gpuistl::GpuVector<double> wrongSize(1);
    BOOST_CHECK_THROW(backend.prepare(linearizer.gpuJacobian(), wrongSize), Opm::CriticalError);

    Opm::Parameters::SetDefault<Opm::Parameters::EnableGpuAssemblySolverBridge>(true);
}

void checkScheduleRejection(const std::string& controls, const std::string& diagnostic)
{
    Fixture fixture(controls);
    auto& sim = *fixture.simulator;
    sim.setEpisodeIndex(1);
    for (const bool enabled : {true, false}) {
        Opm::Parameters::SetDefault<Opm::Parameters::EnableGpuAssemblySolverBridge>(enabled);
        BOOST_CHECK_EXCEPTION(sim.problem().beginEpisode(), std::runtime_error,
            [&diagnostic](const auto& error) {
                return std::string(error.what()).find(diagnostic) != std::string::npos;
            });
    }
    Opm::Parameters::SetDefault<Opm::Parameters::EnableGpuAssemblySolverBridge>(true);
}

BOOST_AUTO_TEST_CASE(RejectsGeometryChanges)
{
    checkScheduleRejection("MULTFLT\n'BRIDGE' 0.5 /\n/\n", "GEO_MODIFIER");
}

BOOST_AUTO_TEST_CASE(RejectsConvectiveMixingChanges)
{
    checkScheduleRejection("DRSDTCON\n0.01 0.3 1e-7 /\n/\n", "DRSDTCON");
}

BOOST_AUTO_TEST_CASE(RejectsInactiveConvectiveMixing)
{
    checkScheduleRejection("DRSDTCON\n-1 0.3 1e-7 /\n/\n", "DRSDTCON");
}

BOOST_FIXTURE_TEST_CASE(RejectsNlddWithGpuAssembly, Fixture)
{
    auto& sim = *simulator;
    Opm::BlackoilModelParameters<double> parameters;
    parameters.nonlinear_solver_ = "nldd";
    for (const bool enabled : {true, false}) {
        Opm::Parameters::SetDefault<Opm::Parameters::EnableGpuAssemblySolverBridge>(enabled);
        BOOST_CHECK_EXCEPTION(
            Opm::NonlinearSystemBlackOilReservoir<TypeTag>(
                sim, parameters, sim.problem().wellModel(), false),
            std::invalid_argument,
            [](const auto& error) {
                return std::string(error.what()).find("GPU assembly does not support NLDD")
                    != std::string::npos;
            });
    }
    Opm::Parameters::SetDefault<Opm::Parameters::EnableGpuAssemblySolverBridge>(true);
}
