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
#define BOOST_TEST_MODULE TarjanTests

#include <boost/test/unit_test.hpp>

#include <opm/common/utility/SparseTable.hpp>
#include <opm/simulators/utils/Tarjan.hpp>

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{

using Edges = std::vector<std::vector<int>>;
using Table = Opm::SparseTable<int>;

Table
toTable(const Edges& graph)
{
    Table table;
    for (const auto& row : graph) {
        table.appendRow(row.begin(), row.end());
    }
    return table;
}

// Components as sorted sets of sorted vertices, ignoring the order of the components and of the
// vertices inside them.
std::set<std::vector<int>>
asSet(const Opm::StronglyConnectedComponents& scc)
{
    std::set<std::vector<int>> sets;
    for (const auto component : scc.components()) {
        std::vector<int> members(component.begin(), component.end());
        std::ranges::sort(members);
        sets.insert(std::move(members));
    }
    return sets;
}

// Check that the result is consistent and every edge goes to the same or an earlier component.
// Returns the first problem, or an empty string, so callers need only one assertion.
std::string
consistencyProblem(const Edges& graph, const Opm::StronglyConnectedComponents& scc)
{
    const int n = static_cast<int>(graph.size());
    if (static_cast<int>(scc.vertices.size()) != n) {
        return "wrong number of vertices";
    }
    if (static_cast<int>(scc.componentOf.size()) != n) {
        return "wrong size of componentOf";
    }
    if (scc.offsets.size() != static_cast<std::size_t>(scc.numComponents()) + 1
        || scc.offsets.front() != 0 || scc.offsets.back() != n) {
        return "inconsistent offsets";
    }

    std::vector<int> sorted = scc.vertices;
    std::ranges::sort(sorted);
    for (int v = 0; v < n; ++v) {
        if (sorted[v] != v) {
            return "vertex " + std::to_string(v) + " is missing or repeated";
        }
    }

    for (int k = 0; k < scc.numComponents(); ++k) {
        if (scc[k].empty()) {
            return "component " + std::to_string(k) + " is empty";
        }
        for (const int v : scc[k]) {
            if (scc.componentOf[v] != k) {
                return "componentOf[" + std::to_string(v) + "] is not " + std::to_string(k);
            }
        }
    }

    for (int v = 0; v < n; ++v) {
        for (const int w : graph[v]) {
            if (scc.componentOf[w] > scc.componentOf[v]) {
                return "edge " + std::to_string(v) + " -> " + std::to_string(w)
                    + " goes to a later component";
            }
        }
    }
    return {};
}

void
checkConsistent(const Edges& graph, const Opm::StronglyConnectedComponents& scc)
{
    const std::string problem = consistencyProblem(graph, scc);
    BOOST_REQUIRE_MESSAGE(problem.empty(), problem);
}

// Two vertices must be in the same component exactly when each can reach the other. Returns a
// description of the first mismatch, or an empty string.
std::string
reachabilityProblem(const std::vector<std::vector<char>>& reach,
                    const Opm::StronglyConnectedComponents& scc)
{
    const int n = static_cast<int>(reach.size());
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const bool mutuallyReachable = reach[i][j] != 0 && reach[j][i] != 0;
            if ((scc.componentOf[i] == scc.componentOf[j]) != mutuallyReachable) {
                return "vertices " + std::to_string(i) + " and " + std::to_string(j)
                    + (mutuallyReachable ? " are mutually reachable but in different components"
                                         : " are in the same component but not mutually reachable");
            }
        }
    }
    return {};
}

} // namespace

BOOST_AUTO_TEST_CASE(EmptyGraph)
{
    const auto scc = Opm::tarjan(Edges {});
    BOOST_CHECK_EQUAL(scc.numComponents(), 0);
    BOOST_CHECK(scc.vertices.empty());
    BOOST_CHECK(scc.offsets == std::vector<int> {0});

    BOOST_CHECK_EQUAL(Opm::tarjan(Table {}).numComponents(), 0);
}

BOOST_AUTO_TEST_CASE(VerticesWithoutEdgesAreSeparateComponents)
{
    const Edges graph(4);
    const auto scc = Opm::tarjan(graph);
    BOOST_CHECK_EQUAL(scc.numComponents(), 4);
    checkConsistent(graph, scc);
}

BOOST_AUTO_TEST_CASE(ChainIsOrderedWithEndOfEdgesFirst)
{
    // 0 -> 1 -> 2 -> 3, i.e. 3 is the most upstream vertex.
    const Edges graph {{1}, {2}, {3}, {}};
    const auto scc = Opm::tarjan(graph);
    BOOST_REQUIRE_EQUAL(scc.numComponents(), 4);
    BOOST_CHECK(scc.vertices == std::vector<int>({3, 2, 1, 0}));
    checkConsistent(graph, scc);
}

BOOST_AUTO_TEST_CASE(CycleIsOneComponent)
{
    const Edges graph {{1}, {2}, {0}};
    const auto scc = Opm::tarjan(graph);
    BOOST_REQUIRE_EQUAL(scc.numComponents(), 1);
    BOOST_CHECK_EQUAL(scc[0].size(), 3u);
    checkConsistent(graph, scc);
}

BOOST_AUTO_TEST_CASE(SelfLoopsAndDuplicateEdges)
{
    const Edges graph {{0, 0, 1, 1}, {1}};
    const auto scc = Opm::tarjan(graph);
    BOOST_CHECK_EQUAL(scc.numComponents(), 2);
    BOOST_CHECK_LT(scc.componentOf[1], scc.componentOf[0]);
    checkConsistent(graph, scc);
}

BOOST_AUTO_TEST_CASE(SeveralComponentsInOrder)
{
    // Components {0, 1, 2}, {3, 4}, {5, 6} and {7}, with edges from later to earlier ones.
    const Edges graph {
        {1}, // 0 -> 1
        {2}, // 1 -> 2
        {0}, // 2 -> 0
        {1, 2, 4}, // 3 -> 1, 2, 4
        {3, 5}, // 4 -> 3, 5
        {2, 6}, // 5 -> 2, 6
        {5}, // 6 -> 5
        {4, 6, 7}, // 7 -> 4, 6, 7
    };
    const auto scc = Opm::tarjan(graph);

    const std::set<std::vector<int>> expected {{0, 1, 2}, {3, 4}, {5, 6}, {7}};
    BOOST_CHECK(asSet(scc) == expected);
    BOOST_CHECK_EQUAL(scc.componentOf[7], scc.numComponents() - 1);
    BOOST_CHECK_EQUAL(scc.componentOf[0], 0);
    checkConsistent(graph, scc);
}

BOOST_AUTO_TEST_CASE(SparseTableInput)
{
    // A cycle with a tail: 3 -> 0 -> 1 -> 2 -> 0.
    const Table graph {{1}, {2}, {0}, {0}};
    const auto scc = Opm::tarjan(graph);
    BOOST_REQUIRE_EQUAL(scc.numComponents(), 2);
    BOOST_CHECK_EQUAL(scc[0].size(), 3u);
    BOOST_CHECK_EQUAL(scc[1].size(), 1u);
    BOOST_CHECK_EQUAL(scc[1][0], 3);
}

BOOST_AUTO_TEST_CASE(OtherIntegerTypesForTheEdges)
{
    const std::vector<std::vector<unsigned>> graph {{1u}, {0u}, {1u}};
    const auto scc = Opm::tarjan(graph);
    BOOST_REQUIRE_EQUAL(scc.numComponents(), 2);
    BOOST_CHECK_EQUAL(scc[0].size(), 2u);
    BOOST_CHECK_EQUAL(scc[1].size(), 1u);
    BOOST_CHECK_EQUAL(scc[1][0], 2);
}

BOOST_AUTO_TEST_CASE(EdgesOutsideTheGraphAreRejected)
{
    const Edges tooLarge {{1}};
    const Edges negative {{-1}};
    BOOST_CHECK_THROW(Opm::tarjan(tooLarge), std::out_of_range);
    BOOST_CHECK_THROW(Opm::tarjan(negative), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(VeryLongChainDoesNotOverflowTheStack)
{
    // 0 -> 1 -> ... -> n - 1. A recursive implementation overflows an 8 MB stack on a chain of
    // 100 000 to 200 000 vertices, depending on optimization.
    constexpr int n = 250000;
    std::vector<int> targets(n - 1);
    std::iota(targets.begin(), targets.end(), 1);
    std::vector<int> rowSizes(n, 1);
    rowSizes.back() = 0;
    const Table graph(targets.begin(), targets.end(), rowSizes.begin(), rowSizes.end());

    const auto scc = Opm::tarjan(graph);
    BOOST_REQUIRE_EQUAL(scc.numComponents(), n);
    BOOST_CHECK_EQUAL(scc.vertices.front(), n - 1);
    BOOST_CHECK_EQUAL(scc.vertices.back(), 0);
}

BOOST_AUTO_TEST_CASE(RandomGraphsAgreeWithBruteForceReachability)
{
    // A fixed seed is intended: the test must give the same graphs on every run, so that a failure
    // can be reproduced.
    std::mt19937 rng(20261008);

    for (int trial = 0; trial < 300; ++trial) {
        const int n = 1 + static_cast<int>(rng() % 25);
        const int numEdges = static_cast<int>(rng() % (3 * n + 1));

        Edges graph(n);
        std::vector<std::vector<char>> reach(n, std::vector<char>(n, 0));
        for (int v = 0; v < n; ++v) {
            reach[v][v] = 1;
        }
        for (int e = 0; e < numEdges; ++e) {
            const int v = static_cast<int>(rng() % n);
            const int w = static_cast<int>(rng() % n);
            graph[v].push_back(w);
            reach[v][w] = 1;
        }
        for (int k = 0; k < n; ++k) {
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    if (reach[i][k] != 0 && reach[k][j] != 0) {
                        reach[i][j] = 1;
                    }
                }
            }
        }

        const auto scc = Opm::tarjan(graph);
        checkConsistent(graph, scc);

        // The sparse table must give exactly the same result.
        const auto sccFromTable = Opm::tarjan(toTable(graph));
        BOOST_REQUIRE(sccFromTable.vertices == scc.vertices && sccFromTable.offsets == scc.offsets
                      && sccFromTable.componentOf == scc.componentOf);

        const std::string problem = reachabilityProblem(reach, scc);
        BOOST_REQUIRE_MESSAGE(problem.empty(), problem);
    }
}
