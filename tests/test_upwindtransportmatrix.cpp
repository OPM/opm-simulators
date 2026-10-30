/*
  Copyright 2026 NORCE Norwegian Research Centre AS

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

#define BOOST_TEST_MODULE UpwindTransportMatrixTests

#include <boost/test/unit_test.hpp>

#include <opm/simulators/utils/UpwindTransportSolver.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{

template <int N>
using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, N, N>>;

template <int N>
using Vectors = std::vector<Dune::BlockVector<Dune::FieldVector<double, N>>>;

// A flow of `flux` from cell `from` into cell `to`.
struct Edge
{
    int from;
    int to;
    double flux;
};

// A matrix with a diagonal block in every row, and the off-diagonal blocks of `pattern`. All blocks
// are zero.
template <int N>
Matrix<N>
makeMatrix(const std::vector<std::set<int>>& pattern)
{
    const std::size_t n = pattern.size();
    Matrix<N> matrix(n, n, Matrix<N>::random);
    for (std::size_t r = 0; r < n; ++r) {
        auto row = pattern[r];
        row.insert(static_cast<int>(r));
        matrix.setrowsize(r, row.size());
    }
    matrix.endrowsizes();
    for (std::size_t r = 0; r < n; ++r) {
        matrix.addindex(r, r);
        for (const int c : pattern[r]) {
            matrix.addindex(r, static_cast<std::size_t>(c));
        }
    }
    matrix.endindices();
    matrix = 0.0;
    return matrix;
}

// The pattern of all neighbours in both directions, like a tracer-type matrix.
std::vector<std::set<int>>
symmetricPattern(std::size_t n, const std::vector<Edge>& edges)
{
    std::vector<std::set<int>> pattern(n);
    for (const Edge& edge : edges) {
        pattern[static_cast<std::size_t>(edge.to)].insert(edge.from);
        pattern[static_cast<std::size_t>(edge.from)].insert(edge.to);
    }
    return pattern;
}

// The output of std::mt19937 is fixed by the standard, but that of the distributions is not, so the
// random values are computed from the engine directly to give the same tests on every platform.
int
randomInt(std::mt19937& rng, int low, int high)
{
    return low + static_cast<int>(rng() % static_cast<std::uint32_t>(high - low + 1));
}

// A multiple of 1/1024 in [low, high), where low and high are multiples of 1/1024 too. All
// operations are exact, so the value does not depend on how a compiler rounds or fuses them.
double
randomReal(std::mt19937& rng, double low, double high)
{
    const int steps = static_cast<int>((high - low) * 1024);
    return low + randomInt(rng, 0, steps - 1) / 1024.0;
}

template <int N>
Vectors<N>
randomVectors(std::size_t numVectors, std::size_t n, std::mt19937& rng)
{
    Vectors<N> vectors(numVectors, typename Vectors<N>::value_type(n));
    for (auto& vector : vectors) {
        for (auto& block : vector) {
            for (auto& v : block) {
                v = randomReal(rng, 0.0, 1.0);
            }
        }
    }
    return vectors;
}

// Dense solve of the whole system, by Gaussian elimination with partial pivoting.
template <int N>
Vectors<N>
referenceSolution(const Matrix<N>& matrix, const Vectors<N>& rhs)
{
    const std::size_t size = matrix.N() * N;
    std::vector<double> a(size * size, 0.0);
    for (auto row = matrix.begin(); row != matrix.end(); ++row) {
        for (auto entry = row->begin(); entry != row->end(); ++entry) {
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    a[(row.index() * N + i) * size + entry.index() * N + j] = (*entry)[i][j];
                }
            }
        }
    }
    std::vector<double> b(size * rhs.size());
    for (std::size_t q = 0; q < rhs.size(); ++q) {
        for (std::size_t k = 0; k < size; ++k) {
            b[k * rhs.size() + q] = rhs[q][k / N][k % N];
        }
    }
    Opm::detail::solveDenseSystem(a, b, size, rhs.size());
    Vectors<N> x = rhs;
    for (std::size_t q = 0; q < rhs.size(); ++q) {
        for (std::size_t k = 0; k < size; ++k) {
            x[q][k / N][k % N] = b[k * rhs.size() + q];
        }
    }
    return x;
}

template <int N>
double
maxDifference(const Vectors<N>& x, const Vectors<N>& reference)
{
    double largest = 0.0;
    for (std::size_t q = 0; q < x.size(); ++q) {
        for (std::size_t c = 0; c < x[q].size(); ++c) {
            for (int a = 0; a < N; ++a) {
                largest = std::max(largest, std::abs(x[q][c][a] - reference[q][c][a]));
            }
        }
    }
    return largest;
}

// Random flow between `n` cells, with each cell receiving flow from up to three random cells. With
// `loops` false the flow only goes from lower to higher cell numbers, so there are no loops.
std::vector<Edge>
randomEdges(int n, bool loops, std::mt19937& rng)
{
    std::vector<Edge> edges;
    for (int to = 0; to < n; ++to) {
        const int numUpwind = randomInt(rng, 0, 3);
        for (int k = 0; k < numUpwind; ++k) {
            int from = randomInt(rng, 0, n - 1);
            int target = to;
            if (from == target) {
                continue;
            }
            if (!loops && from > target) {
                std::swap(from, target);
            }
            const double flux = randomReal(rng, 0.125, 2.0);
            edges.push_back({from, target, flux});
        }
    }
    return edges;
}

Opm::UpwindTransportSettings
alwaysDirect()
{
    Opm::UpwindTransportSettings settings;
    settings.maxDirectSize = 100000;
    return settings;
}

Opm::UpwindTransportSettings
alwaysIterate()
{
    Opm::UpwindTransportSettings settings;
    settings.maxDirectSize = 1;
    settings.relTolerance = 1e-14;
    return settings;
}

} // namespace

BOOST_AUTO_TEST_CASE(DependencyGraphIgnoresZeroBlocks)
{
    // Pattern with all neighbours, but only 0 -> 1 -> 2 filled.
    auto matrix = makeMatrix<1>({{1}, {0, 2}, {1}});
    matrix[1][0] = -1.0;
    matrix[2][1] = -1.0;
    const Opm::detail::MatrixCoefficients coefficients(matrix);
    const auto& graph = coefficients.graph();
    BOOST_REQUIRE_EQUAL(graph.size(), 3);
    BOOST_CHECK(graph[0].empty());
    BOOST_REQUIRE_EQUAL(graph[1].size(), 1);
    BOOST_CHECK_EQUAL(graph[1][0], 0);
    BOOST_REQUIRE_EQUAL(graph[2].size(), 1);
    BOOST_CHECK_EQUAL(graph[2][0], 1);
}

BOOST_AUTO_TEST_CASE(DependencyGraphUsesAnyNonzeroEntryOfABlock)
{
    auto matrix = makeMatrix<2>({{1}, {0}});
    matrix[0][1][1][1] = -1.0; // Only the second component depends on cell 1.
    const Opm::detail::MatrixCoefficients coefficients(matrix);
    const auto& graph = coefficients.graph();
    BOOST_REQUIRE_EQUAL(graph[0].size(), 1);
    BOOST_CHECK_EQUAL(graph[0][0], 1);
    BOOST_CHECK(graph[1].empty());

    // The coupling of row 0 is the block A[0][1].
    int numCouplings = 0;
    coefficients.forEachCoupling(0, [&](std::size_t u, const auto& block) {
        ++numCouplings;
        BOOST_CHECK_EQUAL(u, 1u);
        BOOST_CHECK_EQUAL(block[1][1], -1.0);
    });
    BOOST_CHECK_EQUAL(numCouplings, 1);
}

// The 1x1 matrix of an upwind scheme, with and without flow loops.
BOOST_AUTO_TEST_CASE(ScalarMatchesDenseSolve)
{
    std::mt19937 rng(1234);
    const double dt = 0.7;
    for (const bool loops : {false, true}) {
        for (int trial = 0; trial < 20; ++trial) {
            const int n = 60;
            const auto edges = randomEdges(n, loops, rng);

            auto matrix = makeMatrix<1>(symmetricPattern(n, edges));
            for (int c = 0; c < n; ++c) {
                matrix[c][c] = randomReal(rng, 0.5, 2.0);
            }
            for (const Edge& edge : edges) {
                matrix[edge.from][edge.from] += dt * edge.flux;
                matrix[edge.to][edge.from] -= dt * edge.flux;
            }

            const auto rhs = randomVectors<1>(3, n, rng);
            const auto expected = referenceSolution(matrix, rhs);
            for (const auto& settings : {alwaysDirect(), alwaysIterate()}) {
                Vectors<1> x(3, typename Vectors<1>::value_type(n));
                for (auto& v : x) {
                    v = 0.0;
                }
                const auto statistics = Opm::solveReordered(matrix, rhs, x, settings);
                if (!loops) {
                    BOOST_CHECK_EQUAL(statistics.numLoops, 0);
                    BOOST_CHECK_EQUAL(statistics.numComponents, n);
                }
                BOOST_CHECK_SMALL(maxDifference(x, expected), 1e-12);
            }
        }
    }
}

// 2x2 blocks like the tracer: the two components are coupled inside each cell, and flow in
// different directions, so the cells form loops even though each component on its own has none.
BOOST_AUTO_TEST_CASE(BlocksWithOppositeFlowMatchDenseSolve)
{
    std::mt19937 rng(42);
    const double dt = 1.3;
    for (int trial = 0; trial < 20; ++trial) {
        const int n = 40;
        // Component 0 flows from lower to higher cells, component 1 the other way.
        const auto edges0 = randomEdges(n, false, rng);
        auto edges1 = randomEdges(n, false, rng);
        for (auto& edge : edges1) {
            std::swap(edge.from, edge.to);
        }
        auto allEdges = edges0;
        allEdges.insert(allEdges.end(), edges1.begin(), edges1.end());

        auto matrix = makeMatrix<2>(symmetricPattern(n, allEdges));
        for (int c = 0; c < n; ++c) {
            const double s01 = randomReal(rng, 0.0, 0.5);
            const double s10 = randomReal(rng, 0.0, 0.5);
            const double volume0 = randomReal(rng, 0.5, 2.0);
            const double volume1 = randomReal(rng, 0.5, 2.0);
            matrix[c][c][0][0] = volume0 + s10;
            matrix[c][c][1][1] = volume1 + s01;
            matrix[c][c][0][1] = -s01;
            matrix[c][c][1][0] = -s10;
        }
        const std::vector<std::pair<int, const std::vector<Edge>*>> components {{0, &edges0},
                                                                               {1, &edges1}};
        for (const auto& [component, edges] : components) {
            for (const Edge& edge : *edges) {
                matrix[edge.from][edge.from][component][component] += dt * edge.flux;
                matrix[edge.to][edge.from][component][component] -= dt * edge.flux;
            }
        }

        const auto rhs = randomVectors<2>(2, n, rng);
        const auto expected = referenceSolution(matrix, rhs);
        for (const auto& settings : {alwaysDirect(), alwaysIterate()}) {
            Vectors<2> x(2, typename Vectors<2>::value_type(n));
            for (auto& v : x) {
                v = 0.0;
            }
            const auto statistics = Opm::solveReordered(matrix, rhs, x, settings);
            BOOST_CHECK_SMALL(maxDifference(x, expected), 1e-11);
            if (!edges0.empty() && !edges1.empty()) {
                BOOST_CHECK_GT(statistics.numLoops, 0);
            }
        }
    }
}

// A loop of cells with full 2x2 blocks also off the diagonal.
BOOST_AUTO_TEST_CASE(FullBlocksInALoop)
{
    const int n = 5;
    std::vector<std::set<int>> pattern(n);
    for (int c = 0; c < n; ++c) {
        pattern[c].insert((c + n - 1) % n);
    }
    auto matrix = makeMatrix<2>(pattern);
    for (int c = 0; c < n; ++c) {
        matrix[c][c] = Dune::FieldMatrix<double, 2, 2> {{4.0, 0.5}, {-0.3, 3.0}};
        matrix[c][(c + n - 1) % n]
            = Dune::FieldMatrix<double, 2, 2> {{-1.0, -0.2}, {-0.1, -0.8}};
    }
    std::mt19937 rng(7);
    const auto rhs = randomVectors<2>(1, n, rng);
    const auto expected = referenceSolution(matrix, rhs);
    for (const auto& settings : {alwaysDirect(), alwaysIterate()}) {
        Vectors<2> x(1, typename Vectors<2>::value_type(n));
        x[0] = 0.0;
        const auto statistics = Opm::solveReordered(matrix, rhs, x, settings);
        BOOST_CHECK_EQUAL(statistics.numComponents, 1);
        BOOST_CHECK_EQUAL(statistics.numLoops, 1);
        BOOST_CHECK_EQUAL(statistics.largestComponent, n);
        BOOST_CHECK_SMALL(maxDifference(x, expected), 1e-12);
    }
}

BOOST_AUTO_TEST_CASE(NoRightHandSides)
{
    auto matrix = makeMatrix<1>({{}, {0}});
    matrix[0][0] = 1.0;
    matrix[1][1] = 1.0;
    matrix[1][0] = -1.0;
    Vectors<1> rhs;
    Vectors<1> x;
    const auto statistics = Opm::solveReordered(matrix, rhs, x);
    BOOST_CHECK_EQUAL(statistics.numComponents, 2);
}

BOOST_AUTO_TEST_CASE(InvalidInput)
{
    auto matrix = makeMatrix<1>({{}, {0}});
    matrix[0][0] = 1.0;
    matrix[1][1] = 1.0;
    Vectors<1> rhs(1, Dune::BlockVector<Dune::FieldVector<double, 1>>(2));
    rhs[0] = 1.0;

    // Wrong vector lengths
    Vectors<1> x(1, Dune::BlockVector<Dune::FieldVector<double, 1>>(3));
    BOOST_CHECK_THROW(Opm::solveReordered(matrix, rhs, x), std::invalid_argument);
    Vectors<1> tooMany(2, Dune::BlockVector<Dune::FieldVector<double, 1>>(2));
    BOOST_CHECK_THROW(Opm::solveReordered(matrix, rhs, tooMany), std::invalid_argument);

    // Missing diagonal block
    Matrix<1> noDiagonal(2, 2, Matrix<1>::random);
    noDiagonal.setrowsize(0, 1);
    noDiagonal.setrowsize(1, 1);
    noDiagonal.endrowsizes();
    noDiagonal.addindex(0, 0);
    noDiagonal.addindex(1, 0);
    noDiagonal.endindices();
    noDiagonal = 1.0;
    x[0].resize(2);
    BOOST_CHECK_THROW(Opm::solveReordered(noDiagonal, rhs, x), std::invalid_argument);

    // Singular diagonal block
    matrix[1][1] = 0.0;
    BOOST_CHECK_THROW(Opm::solveReordered(matrix, rhs, x), std::runtime_error);

    // Singular loop
    auto loop = makeMatrix<1>({{1}, {0}});
    loop[0][0] = 1.0;
    loop[0][1] = -1.0;
    loop[1][1] = 1.0;
    loop[1][0] = -1.0;
    BOOST_CHECK_THROW(Opm::solveReordered(loop, rhs, x, alwaysDirect()), std::runtime_error);

    // No convergence
    BOOST_CHECK_THROW(Opm::solveReordered(loop, rhs, x, alwaysIterate()), std::runtime_error);
}

// Each component of the blocks converges on its own scale in Gauss-Seidel. Here the second
// component is a million times smaller than the first and converges much more slowly, like a
// solution part next to its free part. With one scale for the whole block, the sweeps would stop
// when the first component has converged, long before the second.
BOOST_AUTO_TEST_CASE(SmallComponentConvergesOnItsOwnScale)
{
    // A chain of cells coupled in both directions, so that all of them are one loop
    const int n = 200;
    std::vector<std::set<int>> pattern(n);
    for (int c = 0; c + 1 < n; ++c) {
        pattern[c].insert(c + 1);
        pattern[c + 1].insert(c);
    }
    auto matrix = makeMatrix<2>(pattern);
    const std::array<double, 2> coupling {0.3, 0.95}; // fraction of the diagonal
    for (int c = 0; c < n; ++c) {
        for (int a = 0; a < 2; ++a) {
            matrix[c][c][a][a] = 1.0;
        }
        for (const int u : pattern[c]) {
            for (int a = 0; a < 2; ++a) {
                matrix[c][u][a][a] = -coupling[a] / 2;
            }
        }
    }
    std::mt19937 rng(3);
    auto rhs = randomVectors<2>(1, n, rng);
    for (auto& block : rhs[0]) {
        block[1] *= 1e-6;
    }
    const auto expected = referenceSolution(matrix, rhs);

    Opm::UpwindTransportSettings settings;
    settings.maxDirectSize = 1;
    settings.relTolerance = 1e-4;
    Vectors<2> x(1, typename Vectors<2>::value_type(n));
    x[0] = 0.0;
    const auto statistics = Opm::solveReordered(matrix, rhs, x, settings);
    BOOST_CHECK_EQUAL(statistics.largestComponent, n);

    for (int a = 0; a < 2; ++a) {
        double error = 0.0;
        double scale = 0.0;
        for (int c = 0; c < n; ++c) {
            error = std::max(error, std::abs(x[0][c][a] - expected[0][c][a]));
            scale = std::max(scale, std::abs(expected[0][c][a]));
        }
        BOOST_TEST_MESSAGE("component " << a << ": relative error " << error / scale << " after "
                           << statistics.maxIterations << " sweeps");
        // relTolerance * rho / (1 - rho) with rho = 0.95 is about 2e-3
        BOOST_CHECK_LT(error / scale, 1e-2);
    }
}

// NaN and infinite values in the input are rejected with std::runtime_error, so that the caller can
// fall back to another solver.
BOOST_AUTO_TEST_CASE(NonFiniteValues)
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    const auto throwsWith = [](auto&& solve, const std::string& text) {
        try {
            solve();
        } catch (const std::runtime_error& e) {
            return std::string(e.what()).find(text) != std::string::npos;
        }
        return false;
    };

    // Cell 1 receives flow from cell 0
    auto matrix = makeMatrix<2>({{}, {0}});
    const auto reset = [&matrix] {
        matrix = 0.0;
        for (int c = 0; c < 2; ++c) {
            matrix[c][c][0][0] = 1.0;
            matrix[c][c][1][1] = 1.0;
        }
        matrix[1][0][0][0] = -0.5;
    };
    Vectors<2> rhs(1, typename Vectors<2>::value_type(2));
    Vectors<2> x = rhs;
    rhs[0] = 1.0;
    x[0] = 0.0;

    reset();
    BOOST_CHECK_NO_THROW(Opm::solveReordered(matrix, rhs, x));

    // Off-diagonal block, in an entry that is not on the block diagonal
    reset();
    matrix[1][0][1][0] = nan;
    BOOST_CHECK(throwsWith([&] { Opm::solveReordered(matrix, rhs, x); },
                           "the block of row 1, column 0"));

    // Diagonal block, checked up front since b / inf would be finite
    reset();
    matrix[0][0][0][1] = inf;
    BOOST_CHECK(throwsWith([&] { Opm::solveReordered(matrix, rhs, x); },
                           "diagonal block of row 0"));

    // Right-hand side
    reset();
    auto badRhs = rhs;
    badRhs[0][1][1] = -inf;
    BOOST_CHECK(throwsWith([&] { Opm::solveReordered(matrix, badRhs, x); },
                           "the right-hand side in row 1"));

    // An initial value is only used by Gauss-Seidel, so elsewhere it is ignored
    auto badX = x;
    badX[0][0][0] = nan;
    BOOST_CHECK_NO_THROW(Opm::solveReordered(matrix, rhs, badX));
    BOOST_CHECK(std::isfinite(badX[0][0][0]));

    // A loop of two cells
    auto loop = makeMatrix<2>({{1}, {0}});
    loop = 0.0;
    for (int c = 0; c < 2; ++c) {
        loop[c][c][0][0] = 2.0;
        loop[c][c][1][1] = 2.0;
        loop[c][1 - c][0][0] = -0.5;
        loop[c][1 - c][1][1] = -0.5;
    }
    badX = x;
    badX[0][1][0] = nan;
    BOOST_CHECK_NO_THROW(Opm::solveReordered(loop, rhs, badX, alwaysDirect()));
    // Whichever cell Gauss-Seidel updates first reads the initial value of the other one
    badX = x;
    badX[0] = nan;
    BOOST_CHECK(throwsWith([&] { Opm::solveReordered(loop, rhs, badX, alwaysIterate()); },
                           "an initial value of Gauss-Seidel"));

    // An infinite coupling in a directly solved loop, where elimination gives the finite solution
    // [-1, -0] for this matrix
    auto scalarLoop = makeMatrix<1>({{1}, {0}});
    scalarLoop[0][0] = 1.0;
    scalarLoop[0][1] = -inf;
    scalarLoop[1][0] = -1.0;
    scalarLoop[1][1] = 2.0;
    Vectors<1> scalarRhs(1, typename Vectors<1>::value_type(2));
    Vectors<1> scalarX = scalarRhs;
    scalarRhs[0] = 1.0;
    BOOST_CHECK(throwsWith(
        [&] { Opm::solveReordered(scalarLoop, scalarRhs, scalarX, alwaysDirect()); },
        "the block of row 0, column 1"));
}

// With a zero right-hand side and a nonzero initial guess the Gauss-Seidel iterates go to zero by a
// constant factor per sweep, here 1/16, so the change relative to the values never gets small. The
// absolute tolerance stops it.
BOOST_AUTO_TEST_CASE(SolutionGoingToZeroConverges)
{
    auto matrix = makeMatrix<1>({{1}, {0}});
    matrix[0][0] = 2.0;
    matrix[0][1] = -0.5;
    matrix[1][0] = -0.5;
    matrix[1][1] = 2.0;
    Vectors<1> rhs(1, typename Vectors<1>::value_type(2));
    rhs[0] = 0.0;

    auto settings = alwaysIterate();
    Vectors<1> x = rhs;
    x[0] = 1.0;
    const auto statistics = Opm::solveReordered(matrix, rhs, x, settings);
    BOOST_CHECK_LE(statistics.maxIterations, settings.maxIterations);
    for (int c = 0; c < 2; ++c) {
        BOOST_CHECK_LE(std::abs(x[0][c][0]), 1e-20);
    }

    // Without it the relative test alone never stops
    settings.absTolerance = 0.0;
    x[0] = 1.0;
    BOOST_CHECK_THROW(Opm::solveReordered(matrix, rhs, x, settings), std::runtime_error);
}
