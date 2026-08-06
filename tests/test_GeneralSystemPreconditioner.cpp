/*
  Copyright Equinor ASA 2026

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
#define BOOST_TEST_MODULE OPM_test_GeneralSystemPreconditioner
#include <boost/test/unit_test.hpp>

#include <opm/simulators/linalg/system/GeneralSystemPreconditioner.hpp>
#include <opm/simulators/linalg/system/SystemPreconditioner.hpp>
#include <opm/simulators/linalg/system/SystemTypes.hpp>

#include <opm/simulators/linalg/PropertyTree.hpp>

#include <cstddef>
#include <functional>
#include <vector>

namespace {

using Scalar = double;
using RRMatrix = Opm::RRMatrix<Scalar>;
using RWMatrix = Opm::RWMatrix<Scalar>;
using WRMatrix = Opm::WRMatrix<Scalar>;
using WWMatrix = Opm::WWMatrix<Scalar>;
using SystemVector = Opm::SystemVector<Scalar>;

using Fixed = Opm::SystemPreconditioner<Scalar, Opm::SeqResOperator<Scalar>>;
using General = Opm::GeneralSystemPreconditioner<Scalar, Opm::SeqResOperator<Scalar>>;

constexpr int numRes = Opm::numResDofs;
constexpr int numWell = Opm::numWellDofs;
constexpr int pressureIndex = 0;

// 4 reservoir cells and 2 wells: well 0 a standard well (1 block row,
// perforating cells 0 and 1), well 1 a multisegment well (3 block rows), so
// the per-well aggregation is exercised.
constexpr std::size_t numCells = 4;
constexpr std::size_t numWellBlocks = 4;

struct BlockSpec
{
    std::size_t column;
    Scalar base;
};

using MatrixPattern = std::vector<std::vector<BlockSpec>>;

// Diagonally dominant so that every sub-solve is well posed, and
// non-symmetric so a transposed application cannot pass by accident.
template <class Block>
Block makeBlock(const Scalar base, const bool diagonal)
{
    Block block;
    for (int row = 0; row < Block::rows; ++row) {
        for (int col = 0; col < Block::cols; ++col) {
            block[row][col] = 0.05 * (base + 3.0 * row + 7.0 * col);
            if (diagonal && row == col) {
                block[row][col] += 40.0 + base;
            }
        }
    }
    return block;
}

template <class Matrix>
Matrix buildMatrix(const std::size_t rows, const std::size_t cols,
                   const MatrixPattern& pattern, const bool diagonal)
{
    std::size_t nonzeroes = 0;
    for (const auto& row : pattern) {
        nonzeroes += row.size();
    }
    Matrix matrix(rows, cols, nonzeroes, Matrix::row_wise);
    for (auto row = matrix.createbegin(); row != matrix.createend(); ++row) {
        for (const auto& e : pattern[row.index()]) {
            row.insert(e.column);
        }
    }
    for (std::size_t row = 0; row < rows; ++row) {
        for (const auto& e : pattern[row]) {
            matrix[row][e.column] =
                makeBlock<typename Matrix::block_type>(e.base, diagonal && e.column == row);
        }
    }
    return matrix;
}

struct Fixture
{
    RRMatrix A;
    RWMatrix C;
    WRMatrix B;
    WWMatrix D;
    Opm::WellDofLayout layout;
    Opm::SystemMatrix<Scalar> S;
    SystemVector weights;

    Fixture()
    {
        A = buildMatrix<RRMatrix>(numCells, numCells,
                                  {{{0, 1.0}, {1, 2.0}},
                                   {{0, 3.0}, {1, 4.0}, {2, 5.0}},
                                   {{1, 6.0}, {2, 7.0}, {3, 8.0}},
                                   {{2, 9.0}, {3, 10.0}}},
                                  /*diagonal=*/true);

        C = buildMatrix<RWMatrix>(numCells, numWellBlocks,
                                  {{{0, 11.0}},
                                   {{0, 12.0}},
                                   {{1, 13.0}, {2, 14.0}},
                                   {{1, 15.0}, {3, 16.0}}},
                                  false);

        B = buildMatrix<WRMatrix>(numWellBlocks, numCells,
                                  {{{0, 17.0}, {1, 18.0}},
                                   {{2, 19.0}},
                                   {{2, 20.0}},
                                   {{3, 21.0}}},
                                  false);

        D = buildMatrix<WWMatrix>(numWellBlocks, numWellBlocks,
                                  {{{0, 22.0}},
                                   {{1, 23.0}, {2, 24.0}},
                                   {{1, 25.0}, {2, 26.0}, {3, 27.0}},
                                   {{2, 28.0}, {3, 29.0}}},
                                  /*diagonal=*/true);

        layout.wellBlockOffsets = {0, 1, 4};
        layout.pressureDofIndex = numWell - 1;

        S.A = &A;
        S.C = &C;
        S.B = &B;
        S.D = &D;
        S.wellLayout = &layout;

        weights[Dune::Indices::_0].resize(numCells);
        for (std::size_t c = 0; c < numCells; ++c) {
            for (int i = 0; i < numRes; ++i) {
                weights[Dune::Indices::_0][c][i] = 0.5 + 0.1 * c + 0.3 * i;
            }
        }
        weights[Dune::Indices::_1].resize(numWellBlocks);
        for (std::size_t wb = 0; wb < numWellBlocks; ++wb) {
            for (int i = 0; i < numWell; ++i) {
                weights[Dune::Indices::_1][wb][i] = 0.2 + 0.7 * wb - 0.15 * i;
            }
        }
    }

    std::function<SystemVector()> weightsCalculator() const
    {
        const auto w = weights;
        return [w]() { return w; };
    }
};

// One "apply once" sub-solve: run the configured preconditioner a single time
// without a Krylov method around it.
void putOnceSolver(Opm::PropertyTree& prm, const std::string& at, const std::string& type)
{
    prm.put(at + ".maxiter", 1);
    prm.put(at + ".tol", 1e-2);
    prm.put(at + ".verbosity", 0);
    prm.put(at + ".solver", std::string{"preconditioner2inverseoperator"});
    prm.put(at + ".preconditioner.type", type);
    prm.put(at + ".preconditioner.relaxation", 1.0);
}

void putReservoirSolver(Opm::PropertyTree& prm, const std::string& at, const bool addWells)
{
    putOnceSolver(prm, at, "cpr");
    prm.put(at + ".preconditioner.use_well_weights", std::string{"false"});
    prm.put(at + ".preconditioner.add_wells", addWells ? std::string{"true"} : std::string{"false"});
    prm.put(at + ".preconditioner.weight_type", std::string{"trueimpes"});
    prm.put(at + ".preconditioner.pre_smooth", 0);
    prm.put(at + ".preconditioner.post_smooth", 0);
    prm.put(at + ".preconditioner.finesmoother.type", std::string{"jac"});
    prm.put(at + ".preconditioner.finesmoother.relaxation", 1.0);
    prm.put(at + ".preconditioner.verbosity", 0);
    // The coarse pressure system is tiny here, so a single ILU0 apply is both
    // adequate and free of any optional dependency.
    prm.put(at + ".preconditioner.coarsesolver.maxiter", 1);
    prm.put(at + ".preconditioner.coarsesolver.tol", 1e-1);
    prm.put(at + ".preconditioner.coarsesolver.solver", std::string{"preconditioner2inverseoperator"});
    prm.put(at + ".preconditioner.coarsesolver.verbosity", 0);
    prm.put(at + ".preconditioner.coarsesolver.preconditioner.type", std::string{"ilu0"});
    prm.put(at + ".preconditioner.coarsesolver.preconditioner.relaxation", 1.0);
}

void putWellOptions(Opm::PropertyTree& prm, const std::string& transfer)
{
    prm.put("well_transfer", transfer);
    prm.put("well_coarse_diagonal", std::string{"contract_d"});
    prm.put("well_identity_on_pressure_control", std::string{"false"});
    prm.put("verbosity", 0);
}

Opm::PropertyTree fixedTree(const bool addWells, const std::string& transfer)
{
    Opm::PropertyTree prm;
    prm.put("type", std::string{"system_cpr"});
    putWellOptions(prm, transfer);
    putReservoirSolver(prm, "reservoir_solver", addWells);
    putOnceSolver(prm, "reservoir_smoother", "ilu0");
    putOnceSolver(prm, "well_solver", "ilu0");
    return prm;
}

Opm::PropertyTree generalTree(const bool addWells, const std::string& transfer,
                              const bool coarseWellSolve = true,
                              const bool smootherWellSolve = true)
{
    Opm::PropertyTree prm;
    prm.put("type", std::string{"general_system_cpr"});
    putWellOptions(prm, transfer);
    putReservoirSolver(prm, "coarse_solver.reservoir_solver", addWells);
    if (coarseWellSolve) {
        putOnceSolver(prm, "coarse_solver.well_solver", "ilu0");
    }
    putOnceSolver(prm, "smoother.reservoir_smoother", "ilu0");
    if (smootherWellSolve) {
        putOnceSolver(prm, "smoother.well_solver", "ilu0");
    }
    return prm;
}

SystemVector makeRhs()
{
    SystemVector d;
    d[Dune::Indices::_0].resize(numCells);
    d[Dune::Indices::_1].resize(numWellBlocks);
    for (std::size_t c = 0; c < numCells; ++c) {
        for (int i = 0; i < numRes; ++i) {
            d[Dune::Indices::_0][c][i] = 1.0 + 0.37 * c - 0.11 * i;
        }
    }
    for (std::size_t wb = 0; wb < numWellBlocks; ++wb) {
        for (int i = 0; i < numWell; ++i) {
            d[Dune::Indices::_1][wb][i] = 2.0 - 0.23 * wb + 0.41 * i;
        }
    }
    return d;
}

SystemVector zeroLike(const SystemVector& d)
{
    SystemVector v;
    v[Dune::Indices::_0].resize(d[Dune::Indices::_0].size());
    v[Dune::Indices::_1].resize(d[Dune::Indices::_1].size());
    v[Dune::Indices::_0] = 0.0;
    v[Dune::Indices::_1] = 0.0;
    return v;
}

// Bit for bit, not merely close: the point of the decomposition is that the
// composed sweep performs the same floating point operations in the same order.
void checkIdentical(const SystemVector& a, const SystemVector& b)
{
    BOOST_REQUIRE_EQUAL(a[Dune::Indices::_0].size(), b[Dune::Indices::_0].size());
    BOOST_REQUIRE_EQUAL(a[Dune::Indices::_1].size(), b[Dune::Indices::_1].size());
    for (std::size_t c = 0; c < a[Dune::Indices::_0].size(); ++c) {
        for (int i = 0; i < numRes; ++i) {
            BOOST_CHECK_EQUAL(a[Dune::Indices::_0][c][i], b[Dune::Indices::_0][c][i]);
        }
    }
    for (std::size_t wb = 0; wb < a[Dune::Indices::_1].size(); ++wb) {
        for (int i = 0; i < numWell; ++i) {
            BOOST_CHECK_EQUAL(a[Dune::Indices::_1][wb][i], b[Dune::Indices::_1][wb][i]);
        }
    }
}

bool differs(const SystemVector& a, const SystemVector& b)
{
    for (std::size_t c = 0; c < a[Dune::Indices::_0].size(); ++c) {
        for (int i = 0; i < numRes; ++i) {
            if (a[Dune::Indices::_0][c][i] != b[Dune::Indices::_0][c][i]) {
                return true;
            }
        }
    }
    for (std::size_t wb = 0; wb < a[Dune::Indices::_1].size(); ++wb) {
        for (int i = 0; i < numWell; ++i) {
            if (a[Dune::Indices::_1][wb][i] != b[Dune::Indices::_1][wb][i]) {
                return true;
            }
        }
    }
    return false;
}

// Both preconditioners applied once to the same right-hand side.
std::pair<SystemVector, SystemVector>
applyBoth(const Fixture& f, const Opm::PropertyTree& fixedPrm, const Opm::PropertyTree& generalPrm)
{
    const auto d = makeRhs();

    Fixed fixed(f.S, f.weightsCalculator(), pressureIndex, fixedPrm);
    General general(f.S, f.weightsCalculator(), pressureIndex, generalPrm);

    auto vFixed = zeroLike(d);
    auto vGeneral = zeroLike(d);
    fixed.apply(vFixed, d);
    general.apply(vGeneral, d);
    return {vFixed, vGeneral};
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(MatchesFixedThreeStageWithCprwPressureStage)
{
    const Fixture f;
    for (const auto* transfer : {"classic", "no_prolongation"}) {
        const auto [vFixed, vGeneral] =
            applyBoth(f, fixedTree(true, transfer), generalTree(true, transfer));
        BOOST_TEST_CONTEXT("well_transfer = " << transfer) {
            checkIdentical(vFixed, vGeneral);
        }
    }
}

// well_transfer = full is the one mode where the coarse correction reaches the
// well unknowns, so the well half of the result is a sum of three terms rather
// than two: the coarse correction and the two well solves.  The fixed sequence
// accumulates them as ((0 + dw) + c1) + c2, while the two-level method adds the
// coarse correction to the update and then adds the smoother's own accumulated
// correction, dw + ((0 + c1) + c2).  Same terms, same order, different
// parenthesisation -- so the two agree to a rounding step and no closer.  The
// reservoir half has only two terms and stays exact, and so do the other
// transfer modes, where dw is exactly zero.
BOOST_AUTO_TEST_CASE(MatchesFixedThreeStageToRoundoffWithProlongedWellPressure)
{
    const Fixture f;
    const auto [vFixed, vGeneral] = applyBoth(f, fixedTree(true, "full"), generalTree(true, "full"));

    for (std::size_t c = 0; c < numCells; ++c) {
        for (int i = 0; i < numRes; ++i) {
            BOOST_CHECK_EQUAL(vFixed[Dune::Indices::_0][c][i], vGeneral[Dune::Indices::_0][c][i]);
        }
    }
    for (std::size_t wb = 0; wb < numWellBlocks; ++wb) {
        for (int i = 0; i < numWell; ++i) {
            BOOST_CHECK_CLOSE(vFixed[Dune::Indices::_1][wb][i],
                              vGeneral[Dune::Indices::_1][wb][i], 1e-12);
        }
    }
}

BOOST_AUTO_TEST_CASE(MatchesFixedThreeStageWithReservoirOnlyPressureStage)
{
    const Fixture f;
    const auto [vFixed, vGeneral] =
        applyBoth(f, fixedTree(false, "classic"), generalTree(false, "classic"));
    checkIdentical(vFixed, vGeneral);
}

BOOST_AUTO_TEST_CASE(DroppingAWellSolveChangesTheResult)
{
    // Guards the equivalence above: if the composed parts were somehow being
    // ignored, leaving one out would not change anything.
    const Fixture f;
    const auto both = applyBoth(f, fixedTree(true, "classic"), generalTree(true, "classic"));
    const auto noCoarseWell =
        applyBoth(f, fixedTree(true, "classic"), generalTree(true, "classic", false, true));
    const auto noSmootherWell =
        applyBoth(f, fixedTree(true, "classic"), generalTree(true, "classic", true, false));

    BOOST_CHECK(differs(both.second, noCoarseWell.second));
    BOOST_CHECK(differs(both.second, noSmootherWell.second));
}

BOOST_AUTO_TEST_CASE(RejectsAnEmptyComposition)
{
    const Fixture f;
    Opm::PropertyTree prm;
    prm.put("type", std::string{"general_system_cpr"});
    putWellOptions(prm, "classic");
    BOOST_CHECK_THROW(General(f.S, f.weightsCalculator(), pressureIndex, prm),
                      std::invalid_argument);
}
