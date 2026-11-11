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
#include <opm/models/blackoil/blackoilsolventmodules.hh>
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
#include <tuple>
#include <vector>

namespace Opm::Properties {
namespace TTag {
struct TestSolventTypeTag
{
    using InheritsFrom = std::tuple<TestTypeTag>;
};
} // namespace TTag

template <class TypeTag>
struct EnableSolvent<TypeTag, TTag::TestSolventTypeTag>
{
    static constexpr bool value = true;
};
} // namespace Opm::Properties

namespace
{

using TypeTag = Opm::Properties::TTag::TestTypeTag;
using SolventTypeTag = Opm::Properties::TTag::TestSolventTypeTag;
using Options = Opm::Amg::TrueImpesWeightsOptions;

template <class TT> using SimulatorT = Opm::GetPropType<TT, Opm::Properties::Simulator>;
template <class TT> using ElementContextT = Opm::GetPropType<TT, Opm::Properties::ElementContext>;
template <class TT> using GridViewT = Opm::GetPropType<TT, Opm::Properties::GridView>;
template <class TT> using IndicesT = Opm::GetPropType<TT, Opm::Properties::Indices>;
template <class TT> using VectorT = Opm::GetPropType<TT, Opm::Properties::GlobalEqVector>;

using Indices = IndicesT<TypeTag>;
using Vector = VectorT<TypeTag>;

// Pore volume scaling applies to black-oil models only.  Make sure the test
// models qualify, so that they do not silently fall back to max-abs scaling.
static_assert(Opm::Details::isBlackOilTrueImpesModel<
              Opm::GetPropType<TypeTag, Opm::Properties::Model>, ElementContextT<TypeTag>>);
static_assert(Opm::Details::isBlackOilTrueImpesModel<
              Opm::GetPropType<SolventTypeTag, Opm::Properties::Model>, ElementContextT<SolventTypeTag>>);

template <class TT>
std::unique_ptr<SimulatorT<TT>> initialisedSimulator(const char* deck)
{
    auto simulator = Opm::initSimulator<TT>(deck, "test_trueimpes_weight_scaling");
    simulator->model().applyInitialSolution();
    simulator->setTimeStepSize(86400.0);
    BOOST_REQUIRE_GT(simulator->model().numGridDof(), 0u);
    return simulator;
}

template <class TT>
VectorT<TT> trueImpesWeights(const SimulatorT<TT>& simulator, const Options& options)
{
    VectorT<TT> weights(simulator.model().numGridDof());
    const ElementContextT<TT> elemCtx(simulator);
    const Opm::ElementChunks<GridViewT<TT>, Dune::Partitions::All>
        chunks(simulator.gridView(), Dune::Partitions::all, /*num_chunks=*/1);
    Opm::Amg::getTrueImpesWeights(IndicesT<TT>::pressureSwitchIdx, weights, elemCtx,
                                  simulator.model(), chunks,
                                  /*enable_thread_parallel=*/false, options);
    return weights;
}

Vector analyticWeights(const SimulatorT<TypeTag>& simulator)
{
    Vector weights(simulator.model().numGridDof());
    const ElementContextT<TypeTag> elemCtx(simulator);
    const Opm::ElementChunks<GridViewT<TypeTag>, Dune::Partitions::All>
        chunks(simulator.gridView(), Dune::Partitions::all, /*num_chunks=*/1);
    Opm::Amg::getTrueImpesWeightsAnalytic(Indices::pressureSwitchIdx, weights, elemCtx,
                                          simulator.model(), chunks,
                                          /*enable_thread_parallel=*/false);
    return weights;
}

template <class Block>
Block maxAbsNormalised(Block block)
{
    double maxAbs = 0.0;
    for (const double w : block) {
        maxAbs = std::max(maxAbs, std::abs(w));
    }
    block /= maxAbs;
    return block;
}

struct CellStorage
{
    double phases = 0.0;   // weighted storage of the phase component equations
    double solvent = 0.0;  // weighted storage of the solvent equation
    double porosity = 0.0;
    double solventSaturation = 0.0;
};

template <class TT>
std::vector<CellStorage> weightedStorage(const SimulatorT<TT>& simulator,
                                         const VectorT<TT>& weights)
{
    using Idx = IndicesT<TT>;
    std::vector<CellStorage> result(simulator.model().numGridDof());
    ElementContextT<TT> elemCtx(simulator);
    using Evaluation = Opm::GetPropType<TT, Opm::Properties::Evaluation>;
    Dune::FieldVector<Evaluation, Idx::numEq> storage;
    for (const auto& elem : elements(simulator.gridView())) {
        elemCtx.updatePrimaryStencil(elem);
        elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
        const auto cell = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
        simulator.model().localLinearizer(0).localResidual()
            .computeStorage(storage, elemCtx, /*spaceIdx=*/0, /*timeIdx=*/0);
        auto& r = result[cell];
        for (int eq = Idx::conti0EqIdx; eq < Idx::conti0EqIdx + Idx::numPhases; ++eq) {
            r.phases += weights[cell][eq] * storage[eq].value();
        }
        const auto& intQuants = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
        if constexpr (Idx::enableSolvent) {
            r.solvent = weights[cell][Idx::contiSolventEqIdx]
                * storage[Idx::contiSolventEqIdx].value();
            r.solventSaturation = Opm::getValue(intQuants.solventSaturation());
        }
        r.porosity = Opm::getValue(intQuants.porosity());
    }
    return result;
}

} // namespace

using SimulatorFixture = Opm::SimulatorFixture;
BOOST_GLOBAL_FIXTURE(SimulatorFixture);

BOOST_AUTO_TEST_CASE(PoreVolumeScalingGivesReservoirVolumeUnits)
{
    const auto simulator = initialisedSimulator<TypeTag>("trueimpes_weight_scaling.DATA");
    const auto numDof = simulator->model().numGridDof();

    // The default options: pore volume scaling.
    const Options defaults{};
    BOOST_REQUIRE(defaults.scaling == Options::Scaling::PoreVolume);

    const auto poreVolume = trueImpesWeights<TypeTag>(*simulator, defaults);
    const auto maxAbs = trueImpesWeights<TypeTag>(*simulator, Options{Options::Scaling::MaxAbs});
    const auto analytic = analyticWeights(*simulator);
    const auto storage = weightedStorage<TypeTag>(*simulator, poreVolume);
    const auto analyticStorage = weightedStorage<TypeTag>(*simulator, analytic);

    for (std::size_t cell = 0; cell < numDof; ++cell) {
        BOOST_TEST_CONTEXT("cell " << cell)
        {
            // Reservoir volume units: the weighted storage is the porosity, as
            // it is for the analytic weights.
            BOOST_CHECK_CLOSE(storage[cell].phases, storage[cell].porosity, 1.0e-10);
            BOOST_CHECK_CLOSE(analyticStorage[cell].phases, analyticStorage[cell].porosity, 1.0e-8);

            // Only the scale differs from the max-abs scaled weights.
            const auto normalised = maxAbsNormalised(poreVolume[cell]);
            for (int eq = 0; eq < Indices::numEq; ++eq) {
                BOOST_CHECK_CLOSE(normalised[eq], maxAbs[cell][eq], 1.0e-10);
            }
        }
    }
}

// The solvent occupies pore volume with its own saturation, so its
// conservation equation must be part of the pore volume scaling.  Leaving it
// out would scale the rows by 1/(1 - S_solvent), i.e. differently on either
// side of a solvent front.
BOOST_AUTO_TEST_CASE(PoreVolumeScalingIncludesSolvent)
{
    using SolventIndices = IndicesT<SolventTypeTag>;
    static_assert(SolventIndices::enableSolvent);

    const auto simulator =
        initialisedSimulator<SolventTypeTag>("trueimpes_weight_scaling_solvent.DATA");
    const auto numDof = simulator->model().numGridDof();

    const auto poreVolume = trueImpesWeights<SolventTypeTag>(*simulator, Options{});
    const auto maxAbs =
        trueImpesWeights<SolventTypeTag>(*simulator, Options{Options::Scaling::MaxAbs});
    const auto storage = weightedStorage<SolventTypeTag>(*simulator, poreVolume);

    int numSolventCells = 0;
    for (std::size_t cell = 0; cell < numDof; ++cell) {
        const auto& s = storage[cell];
        BOOST_TEST_CONTEXT("cell " << cell << ", solvent saturation " << s.solventSaturation)
        {
            // All fluids, solvent included, fill the pore volume.
            BOOST_CHECK_CLOSE(s.phases + s.solvent, s.porosity, 1.0e-10);

            // The phases alone fill the remaining pore volume, so the solvent
            // weight is its reservoir volume per surface volume.
            BOOST_CHECK_CLOSE(s.phases, s.porosity * (1.0 - s.solventSaturation), 1.0e-8);

            // Only the scale differs from the max-abs scaled weights.
            const auto normalised = maxAbsNormalised(poreVolume[cell]);
            for (int eq = 0; eq < SolventIndices::numEq; ++eq) {
                BOOST_CHECK_CLOSE(normalised[eq], maxAbs[cell][eq], 1.0e-10);
            }
            numSolventCells += s.solventSaturation > 0.01 ? 1 : 0;
        }
    }

    // The deck has solvent in several cells, otherwise the test is vacuous.
    BOOST_CHECK_GE(numSolventCells, 4);
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
