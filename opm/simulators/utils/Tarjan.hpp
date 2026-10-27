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

#ifndef OPM_TARJAN_HEADER_INCLUDED
#define OPM_TARJAN_HEADER_INCLUDED

#include <algorithm>
#include <cstddef>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

/// \file
///
/// Tarjan's algorithm for the strongly connected components (SCCs) of a directed graph, see
/// R. E. Tarjan, "Depth-first search and linear graph algorithms", SIAM J. Comput. 1(2), 1972, or
/// https://en.wikipedia.org/wiki/Tarjan's_strongly_connected_components_algorithm
///
/// This is a C++20 rewrite of the C implementation used by the reordering transport solvers in
/// opm-simulators-legacy. The original is Copyright (C) 2012 Jostein R. Natvig and is licensed
/// under the MIT license.
///
/// The intended use is to order the cells of a reservoir model for transport solves that follow the
/// flow. Build the graph so that every cell has an edge to each of its upwind neighbours. The
/// components are then returned with all upwind components first, so that solving them in the
/// returned order only ever needs values that are already final:
///
///  - A component with one cell and no flow loop is solved with a single local solve.
///  - A component with several cells is a flow loop: its cells depend on each other and have to
///    be solved together, either directly or iteratively.

namespace Opm
{

/// \brief The strongly connected components of a directed graph, as computed by tarjan().
///
/// The components are numbered in reverse topological order of the condensation of the graph: if
/// there is an edge from a vertex in component `a` to a vertex in a different component `b`, then
/// `b < a`. For an upwind graph (edges point from a cell to the cells it receives flow from) this
/// means that component 0 is the most upstream one and that every component only depends on
/// components with a lower number.
struct StronglyConnectedComponents
{
    /// All vertices, grouped by component. Component `k` occupies `vertices[offsets[k]]`, ...,
    /// `vertices[offsets[k + 1] - 1]`.
    std::vector<int> vertices;

    /// Start of each component in `vertices`, with the total vertex count as the final entry.
    /// Always has `numComponents() + 1` entries.
    std::vector<int> offsets {0};

    /// The component each vertex belongs to.
    std::vector<int> componentOf;

    /// \brief The number of components.
    ///
    /// \return The number of components, i.e. `offsets.size() - 1`.
    int numComponents() const
    {
        return static_cast<int>(offsets.size()) - 1;
    }

    /// \brief The vertices of a component.
    ///
    /// \param k Component index, `0 <= k < numComponents()`.
    /// \return A view of the vertices of component `k`.
    std::span<const int> operator[](int k) const
    {
        const auto first = static_cast<std::size_t>(offsets[k]);
        const auto last = static_cast<std::size_t>(offsets[k + 1]);
        return std::span<const int> {vertices}.subspan(first, last - first);
    }

    /// \brief All components in order, for use in a range-based for loop.
    ///
    /// \return A view with one `std::span<const int>` per component, see operator[]().
    auto components() const
    {
        return std::views::iota(0, numComponents())
            | std::views::transform([this](int k) { return (*this)[k]; });
    }
};

/// \brief Compute the strongly connected components of a directed graph with Tarjan's algorithm.
///
/// Linear in the number of vertices plus edges; iterative, so long chains do not overflow the call
/// stack. Self-edges and duplicate edges are allowed, and the result is deterministic.
///
/// \tparam Graph Adjacency list, e.g. `Opm::SparseTable<int>` or `std::vector<std::vector<int>>`;
///               `graph[v]` must not return a container by value.
/// \param graph The directed graph. For an upwind ordering, each cell has an edge to each of the
///              cells it receives flow from.
/// \return The components in reverse topological order, see Opm::StronglyConnectedComponents.
/// \throws std::out_of_range if an edge refers to a vertex outside the graph.
template <class Graph>
StronglyConnectedComponents
tarjan(const Graph& graph)
{
    using Row = decltype(graph[std::size_t {}]);
    using Iterator = std::ranges::iterator_t<Row>;
    static_assert(std::ranges::forward_range<Row>, "A graph row must be a forward range");
    static_assert(std::is_integral_v<std::ranges::range_value_t<Row>>,
                  "A graph row must hold integer vertex indices");

    struct Frame
    {
        int vertex;
        Iterator next;
        Iterator end;
    };

    constexpr int unvisited = -1;
    const int numVertices = static_cast<int>(std::ranges::size(graph));

    StronglyConnectedComponents result;
    result.vertices.reserve(numVertices);
    result.componentOf.assign(numVertices, unvisited);

    // order[v]: discovery number of v. low[v]: smallest discovery number reachable from v through
    // its search subtree and at most one edge to a vertex that is still waiting to be assigned to a
    // component.
    std::vector<int> order(numVertices, unvisited);
    std::vector<int> low(numVertices, 0);
    std::vector<int> pending; // Visited but not yet assigned to a component.
    std::vector<Frame> search; // Explicit stack replacing the recursion.
    int discovered = 0;

    const auto discover = [&](int vertex) {
        order[vertex] = low[vertex] = discovered++;
        pending.push_back(vertex);
        auto&& row = graph[static_cast<std::size_t>(vertex)];
        search.push_back({vertex, std::ranges::begin(row), std::ranges::end(row)});
    };

    for (int root = 0; root < numVertices; ++root) {
        if (order[root] != unvisited) {
            continue;
        }
        discover(root);

        while (!search.empty()) {
            // Do not keep a reference into 'search' across discover().
            const int v = search.back().vertex;

            if (search.back().next != search.back().end) {
                const auto target = *search.back().next++;
                if (std::cmp_less(target, 0) || std::cmp_greater_equal(target, numVertices)) {
                    throw std::out_of_range("Opm::tarjan: edge from vertex " + std::to_string(v)
                                            + " to vertex " + std::to_string(target)
                                            + " in a graph with " + std::to_string(numVertices)
                                            + " vertices");
                }
                const int w = static_cast<int>(target);

                if (order[w] == unvisited) {
                    discover(w);
                } else if (result.componentOf[w] == unvisited) {
                    // w is on the 'pending' stack: v can reach back to it.
                    low[v] = std::min(low[v], order[w]);
                }
                continue;
            }

            // All edges of v are handled. If nothing in its subtree reaches a vertex discovered
            // before v, v is the root of a component made up of everything above it on the
            // 'pending' stack.
            if (low[v] == order[v]) {
                const int component = result.numComponents();
                int member;
                do {
                    member = pending.back();
                    pending.pop_back();
                    result.componentOf[member] = component;
                    result.vertices.push_back(member);
                } while (member != v);
                result.offsets.push_back(static_cast<int>(result.vertices.size()));
            }

            search.pop_back();
            if (!search.empty()) {
                const int parent = search.back().vertex;
                low[parent] = std::min(low[parent], low[v]);
            }
        }
    }

    return result;
}

} // namespace Opm

#endif // OPM_TARJAN_HEADER_INCLUDED
