/*
  Copyright 2026 EQUINOR ASA

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

#define BOOST_TEST_MODULE GpuPropertyAssemblyBridge
#include <boost/test/unit_test.hpp>

#include <flow/flow_gpu.hpp>
#include <opm/material/common/ResetLocale.hpp>
#include <opm/models/utils/start.hh>
#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/flow/SimulatorFullyImplicit.hpp>
#include <opm/simulators/utils/moduleVersion.hpp>

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

namespace Opm {
std::string moduleVersionName() { return "test-version"; }
}

namespace {
using TypeTag = Opm::Properties::TTag::FlowGasWaterEnergyProblemGPU;
using Simulator = Opm::GetPropType<TypeTag, Opm::Properties::Simulator>;
using IQ = Opm::GetPropType<TypeTag, Opm::Properties::IntensiveQuantities>;
using Matrix = typename Opm::GetPropType<TypeTag, Opm::Properties::SparseMatrixAdapter>::IstlMatrix;
using Vector = Opm::GetPropType<TypeTag, Opm::Properties::GlobalEqVector>;
using RateVector = Opm::GetPropType<TypeTag, Opm::Properties::RateVector>;
using Indices = Opm::GetPropType<TypeTag, Opm::Properties::Indices>;
using LocalResidual = Opm::GetPropType<TypeTag, Opm::Properties::LocalResidual>;
using CpuLinearizer = Opm::TpfaLinearizer<Opm::Properties::TTag::FlowGasWaterEnergyProblem>;
constexpr unsigned numEq = Opm::getPropValue<TypeTag, Opm::Properties::NumEq>();
static_assert(Opm::gpuistl::GpuBlackoilIntensiveQuantitiesDispatcherSupport<TypeTag>::value);

// Self-contained three-cell thermal problem, including a thermal boundary.
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
    std::unique_ptr<Simulator> gpu;
    std::unique_ptr<Simulator> cpu;

    Fixture()
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
        directory = std::filesystem::current_path() / "gpu-property-assembly-bridge-test";
        std::filesystem::create_directories(directory);
        const auto filename = directory / "BRIDGE.DATA";
        { std::ofstream out(filename); out << deck; }
        const auto deckArg = "--ecl-deck-file-name=" + filename.string();
        const auto outputArg = "--output-dir=" + directory.string();
        // Construct the same production physics twice. The dispatcher choice is
        // captured by each model's constructor; the reference never uses GPU IQs.
        for (const bool propertiesOnGpu : {true, false}) {
            const char* args[] = {"test_GpuPropertyAssemblyBridge", deckArg.c_str(), outputArg.c_str(),
                "--linear-solver=dilu", "--linear-solver-accelerator=gpu",
                "--matrix-add-well-contributions=true", "--enable-storage-cache=false",
                "--enable-gpu-assembly-solver-bridge=true", "--threads-per-process=2",
                "--enable-tuning=false", "--check-satfunc-consistency=false",
                propertiesOnGpu ? "--experimental-compute-properties-on-gpu=true"
                                : "--experimental-compute-properties-on-gpu=false"};
            Opm::setupParameters_<TypeTag>(sizeof(args) / sizeof(args[0]), args, false, true, false, 0);
            Opm::FlowGenericVanguard::readDeck(filename.string(), /*ignoredKeywords=*/"");
            auto sim = std::make_unique<Simulator>();
            // Exercises production startup, including the first GPU update while
            // the allocated CPU history slot has not yet been populated.
            sim->model().applyInitialSolution();
            sim->setEpisodeIndex(0);
            sim->setTimeStepSize(100.0);
            (propertiesOnGpu ? gpu : cpu) = std::move(sim);
        }
        BOOST_REQUIRE(gpu->model().hasGpuPropertyAssemblyBridge());
        BOOST_REQUIRE(!cpu->model().hasGpuPropertyAssemblyBridge());
    }

    ~Fixture()
    {
        cpu.reset();
        gpu.reset();
        std::error_code error;
        std::filesystem::remove_all(directory, error);
    }
};

// Independent host IQs and host assembly into separate matrix/residual storage.
// Only topology and geometric neighbor data come from the tested linearizer; no device values
// or GPU-materialized host IQs enter the reference.
std::pair<Matrix, Vector> assembleCpu(Simulator& sim, Opm::TpfaLinearizer<TypeTag>& linearizer)
{
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
    BOOST_CHECK(std::isfinite(actual));
    BOOST_CHECK_SMALL(actual - expected, 1e-9 + 1e-9 * std::abs(expected));
}

Vector checkAssembly(Fixture& fixture, Opm::TpfaLinearizer<TypeTag>& linearizer)
{
    linearizer.linearizeDomain();
    BOOST_REQUIRE(linearizer.gpuSystemActive());
    // Model availability alone would not detect an always-false assembly binding.
    BOOST_REQUIRE(linearizer.usesGpuPropertyAssemblyBridge());
    const auto [matrix, residual] = assembleCpu(*fixture.cpu, linearizer);
    const auto values = linearizer.gpuJacobian().getNonZeroValues().asStdVector();
    std::size_t offset = 0;
    for (auto row = matrix.begin(); row != matrix.end(); ++row) {
        for (auto block = row->begin(); block != row->end(); ++block) {
            for (unsigned eq = 0; eq < numEq; ++eq) {
                for (unsigned pv = 0; pv < numEq; ++pv) {
                    BOOST_TEST_CONTEXT("Jacobian cell " << row.index() << " neighbor " << block.index()
                                       << " equation " << eq << " derivative " << pv) {
                        checkClose(values.at(offset++), (*block)[eq][pv]);
                    }
                }
            }
        }
    }
    BOOST_CHECK_EQUAL(offset, values.size());
    const auto flattened = linearizer.flattenedGpuResidual().asStdVector();
    BOOST_REQUIRE_EQUAL(flattened.size(), residual.size() * numEq);
    for (unsigned cell = 0; cell < residual.size(); ++cell) {
        for (unsigned eq = 0; eq < numEq; ++eq) {
            BOOST_TEST_CONTEXT("Residual cell " << cell << " equation " << eq) {
                checkClose(linearizer.residual()[cell][eq], residual[cell][eq]);
                checkClose(flattened[cell * numEq + eq], residual[cell][eq]);
            }
        }
    }
    return residual;
}

Vector checkAssembly(Fixture& fixture)
{
    return checkAssembly(fixture, fixture.gpu->model().linearizer());
}

void changeCurrentState(Fixture& fixture, double factor)
{
    for (auto* sim : {fixture.gpu.get(), fixture.cpu.get()}) {
        auto& solution = sim->model().solution(0);
        for (unsigned cell = 0; cell < solution.size(); ++cell) {
            solution[cell][Indices::pressureSwitchIdx] *= 1.0 + factor * (cell + 1);
            solution[cell][Indices::temperatureIdx] += 2.0 * (cell + 1);
        }
        sim->model().invalidateAndUpdateIntensiveQuantities(0);
    }
}
}

BOOST_FIXTURE_TEST_CASE(ProductionHandoffAcrossTimeLevels, Fixture)
{
    checkAssembly(*this);
    const auto* current = gpu->model().gpuPropertyAssemblyBridge().intensiveQuantitiesView(0).data();
    const auto* history = gpu->model().gpuPropertyAssemblyBridge().intensiveQuantitiesView(1).data();
    BOOST_REQUIRE(current != history);

    // B/A: nonzero accumulation and flux, with different per-cell states.
    changeCurrentState(*this, 0.01);
    const auto beforeAdvance = checkAssembly(*this);
    BOOST_REQUIRE_GT(beforeAdvance.two_norm(), 1e-6);

    // Advance through the model hook (rather than calling the dispatcher directly).
    gpu->model().advanceTimeLevel();
    cpu->model().advanceTimeLevel();
    const auto afterAdvance = checkAssembly(*this);
    Vector accumulation = beforeAdvance;
    accumulation -= afterAdvance;
    BOOST_REQUIRE_GT(accumulation.two_norm(), 1e-6);

    // C/B detects stale current data, incorrect history shifts and slot aliasing.
    changeCurrentState(*this, 0.02);
    checkAssembly(*this);
    BOOST_CHECK(current == gpu->model().gpuPropertyAssemblyBridge().intensiveQuantitiesView(0).data());
    BOOST_CHECK(history == gpu->model().gpuPropertyAssemblyBridge().intensiveQuantitiesView(1).data());

    // Recreate assembly storage while preserving the model-owned property bridge.
    Opm::TpfaLinearizer<TypeTag> rebuilt;
    rebuilt.init(*gpu);
    checkAssembly(*this, rebuilt);
}

BOOST_FIXTURE_TEST_CASE(CpuPropertiesUseAssemblyFallback, Fixture)
{
    changeCurrentState(*this, 0.01);
    auto& linearizer = cpu->model().linearizer();
    linearizer.linearizeDomain();
    BOOST_REQUIRE(linearizer.gpuSystemActive());
    BOOST_CHECK(!linearizer.usesGpuPropertyAssemblyBridge());
    // The direct path and fallback use independent IQ calculations but the
    // same input states and must produce the same complete system.
    checkAssembly(*this);
    const auto direct = gpu->model().linearizer().gpuJacobian().getNonZeroValues().asStdVector();
    const auto fallback = linearizer.gpuJacobian().getNonZeroValues().asStdVector();
    BOOST_REQUIRE_EQUAL(direct.size(), fallback.size());
    for (std::size_t i = 0; i < direct.size(); ++i) {
        checkClose(fallback[i], direct[i]);
    }
    const auto directResidual = gpu->model().linearizer().flattenedGpuResidual().asStdVector();
    const auto fallbackResidual = linearizer.flattenedGpuResidual().asStdVector();
    BOOST_REQUIRE_EQUAL(directResidual.size(), fallbackResidual.size());
    for (std::size_t i = 0; i < directResidual.size(); ++i) {
        checkClose(fallbackResidual[i], directResidual[i]);
    }
}

BOOST_FIXTURE_TEST_CASE(HistoryInitializationRequiresEveryCpuEntry, Fixture)
{
    using Bridge = Opm::gpuistl::GpuFlowGasWaterEnergyBridge<TypeTag>;
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;
    auto& model = cpu->model();
    const unsigned numCells = model.numGridDof();
    BOOST_REQUIRE_GT(numCells, 1u);

    // Reject both an immediately invalid cache and a partially valid prefix.
    // Allocated entries must not be confused with computed, valid quantities.
    for (const unsigned invalidCell : {0u, numCells - 1}) {
        model.setIntensiveQuantitiesCacheEntryValidity(invalidCell, 1, false);
        Bridge bridge;
        bridge.updatePrimaryVariables(cpu->problem(), model.solution(0), 0);
        BOOST_CHECK(!bridge.hasIntensiveQuantities(1));
        BOOST_CHECK(!bridge.hasModelView());
        model.setIntensiveQuantitiesCacheEntryValidity(invalidCell, 1, true);
    }

    // A fully valid history must be uploaded, preserving values and derivatives.
    Bridge bridge;
    bridge.updatePrimaryVariables(cpu->problem(), model.solution(0), 0);
    BOOST_REQUIRE(bridge.hasIntensiveQuantities(1));
    std::vector<IQ> materialized(numCells);
    std::vector<IQ*> destination(numCells);
    for (unsigned cell = 0; cell < numCells; ++cell) {
        destination[cell] = &materialized[cell];
    }
    bridge.materializeHostIntensiveQuantities(1, destination.data(), numCells);
    for (unsigned cell = 0; cell < numCells; ++cell) {
        const auto& expected = model.intensiveQuantities(cell, 1);
        checkClose(materialized[cell].referencePorosity(), expected.referencePorosity());
        for (const unsigned phase : {FluidSystem::waterPhaseIdx, FluidSystem::gasPhaseIdx}) {
            const auto& density = materialized[cell].fluidState().density(phase);
            const auto& expectedDensity = expected.fluidState().density(phase);
            checkClose(density.value(), expectedDensity.value());
            for (unsigned pv = 0; pv < numEq; ++pv) {
                checkClose(density.derivative(pv), expectedDensity.derivative(pv));
            }
        }
    }
}
