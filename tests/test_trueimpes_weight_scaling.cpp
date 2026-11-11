/*
  Copyright 2026 Equinor ASA

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

#define BOOST_TEST_MODULE TrueImpesWeightScaling

#include "SimulatorFixture.hpp"

#include <opm/grid/utility/ElementChunks.hpp>
#include <opm/simulators/flow/FlowProblemBlackoil.hpp>
#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/TrueImpesWeightsOptions.hpp>
#include <opm/simulators/linalg/getQuasiImpesWeights.hpp>
#include <opm/simulators/linalg/setupPropertyTree.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <vector>

namespace
{

using TypeTag = Opm::Properties::TTag::TestTypeTag;
using Simulator = Opm::GetPropType<TypeTag, Opm::Properties::Simulator>;
using ElementContext = Opm::GetPropType<TypeTag, Opm::Properties::ElementContext>;
using GridView = Opm::GetPropType<TypeTag, Opm::Properties::GridView>;
using Indices = Opm::GetPropType<TypeTag, Opm::Properties::Indices>;
using Vector = Opm::GetPropType<TypeTag, Opm::Properties::GlobalEqVector>;
using Options = Opm::Amg::TrueImpesWeightsOptions;

constexpr int pressureIndex = Indices::pressureSwitchIdx;

std::unique_ptr<Simulator> initialisedSimulator()
{
    auto simulator = Opm::initSimulator<TypeTag>("trueimpes_weight_scaling.DATA",
                                                 "test_trueimpes_weight_scaling");
    simulator->model().applyInitialSolution();
    simulator->setTimeStepSize(86400.0);
    BOOST_REQUIRE_GT(simulator->model().numGridDof(), 0u);
    return simulator;
}

Vector trueImpesWeights(const Simulator& simulator, const Options& options)
{
    Vector weights(simulator.model().numGridDof());
    const ElementContext elemCtx(simulator);
    const Opm::ElementChunks<GridView, Dune::Partitions::All>
        chunks(simulator.gridView(), Dune::Partitions::all, /*num_chunks=*/1);
    Opm::Amg::getTrueImpesWeights(pressureIndex, weights, elemCtx, simulator.model(),
                                  chunks, /*enable_thread_parallel=*/false, options);
    return weights;
}

Vector analyticWeights(const Simulator& simulator)
{
    Vector weights(simulator.model().numGridDof());
    const ElementContext elemCtx(simulator);
    const Opm::ElementChunks<GridView, Dune::Partitions::All>
        chunks(simulator.gridView(), Dune::Partitions::all, /*num_chunks=*/1);
    Opm::Amg::getTrueImpesWeightsAnalytic(pressureIndex, weights, elemCtx, simulator.model(),
                                          chunks, /*enable_thread_parallel=*/false);
    return weights;
}

Vector::block_type maxAbsNormalised(Vector::block_type block)
{
    double maxAbs = 0.0;
    for (const double w : block) {
        maxAbs = std::max(maxAbs, std::abs(w));
    }
    block /= maxAbs;
    return block;
}

// Weighted storage of the mass conservation equations, and porosity, per cell.
std::vector<std::array<double, 2>> weightedStorageAndPorosity(const Simulator& simulator,
                                                              const Vector& weights)
{
    std::vector<std::array<double, 2>> result(simulator.model().numGridDof());
    ElementContext elemCtx(simulator);
    using Evaluation = Opm::GetPropType<TypeTag, Opm::Properties::Evaluation>;
    Dune::FieldVector<Evaluation, Indices::numEq> storage;
    for (const auto& elem : elements(simulator.gridView())) {
        elemCtx.updatePrimaryStencil(elem);
        elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
        const auto cell = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
        simulator.model().localLinearizer(0).localResidual()
            .computeStorage(storage, elemCtx, /*spaceIdx=*/0, /*timeIdx=*/0);
        double weighted = 0.0;
        for (int eq = Indices::conti0EqIdx; eq < Indices::conti0EqIdx + Indices::numPhases; ++eq) {
            weighted += weights[cell][eq] * storage[eq].value();
        }
        result[cell] = {weighted, Opm::getValue(elemCtx.intensiveQuantities(0, 0).porosity())};
    }
    return result;
}

} // namespace

using SimulatorFixture = Opm::SimulatorFixture;
BOOST_GLOBAL_FIXTURE(SimulatorFixture);

BOOST_AUTO_TEST_CASE(PoreVolumeScalingGivesReservoirVolumeUnits)
{
    const auto simulator = initialisedSimulator();
    const auto numDof = simulator->model().numGridDof();

    // The default options: pore volume scaling.
    const Options defaults{};
    BOOST_REQUIRE(defaults.scaling == Options::Scaling::PoreVolume);

    const auto poreVolume = trueImpesWeights(*simulator, defaults);
    const auto maxAbs = trueImpesWeights(*simulator, Options{Options::Scaling::MaxAbs});
    const auto analytic = analyticWeights(*simulator);
    const auto storage = weightedStorageAndPorosity(*simulator, poreVolume);
    const auto analyticStorage = weightedStorageAndPorosity(*simulator, analytic);

    for (std::size_t cell = 0; cell < numDof; ++cell) {
        BOOST_TEST_CONTEXT("cell " << cell)
        {
            // Reservoir volume units: the weighted storage is the porosity, as
            // it is for the analytic weights.
            BOOST_CHECK_CLOSE(storage[cell][0], storage[cell][1], 1.0e-10);
            BOOST_CHECK_CLOSE(analyticStorage[cell][0], analyticStorage[cell][1], 1.0e-8);

            // Only the scale differs from the max-abs scaled weights.
            const auto normalised = maxAbsNormalised(poreVolume[cell]);
            for (int eq = 0; eq < Indices::numEq; ++eq) {
                BOOST_CHECK_CLOSE(normalised[eq], maxAbs[cell][eq], 1.0e-10);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(OptionsFromPropertyTree)
{
    const auto defaults = Options::fromPropertyTree(Opm::PropertyTree{});
    BOOST_CHECK(defaults.scaling == Options::Scaling::PoreVolume);

    Opm::PropertyTree prm;
    prm.put("weight_scaling", std::string{"max_abs"});
    const auto options = Options::fromPropertyTree(prm);
    BOOST_CHECK(options.scaling == Options::Scaling::MaxAbs);

    // Round trip through put().
    Opm::PropertyTree stored;
    options.put(stored, "preconditioner.");
    const auto restored = Options::fromPropertyTree(stored.get_child("preconditioner"));
    BOOST_CHECK(restored.scaling == options.scaling);

    Opm::PropertyTree bad;
    bad.put("weight_scaling", std::string{"unknown"});
    BOOST_CHECK_THROW(Options::fromPropertyTree(bad), std::invalid_argument);

    // The generated true-IMPES configurations state the defaults explicitly.
    Opm::FlowLinearSolverParameters parameters;
    for (const auto& tree : {Opm::setupCPRW("cprw", parameters),
                             Opm::setupCPR("cpr_trueimpes", parameters)}) {
        BOOST_CHECK_EQUAL(tree.get<std::string>("preconditioner.weight_scaling"), "pore_volume");
    }
    const auto system = Opm::setupSystemCPR("system_cprw", parameters);
    BOOST_CHECK_EQUAL(system.get<std::string>(
                          "preconditioner.reservoir_solver.preconditioner.weight_scaling"),
                      "pore_volume");
}
