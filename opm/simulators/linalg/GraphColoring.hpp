/*
  Copyright 2018 Equinor

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
#ifndef OPM_GRAPHCOLORING_HEADER_INCLUDED
#define OPM_GRAPHCOLORING_HEADER_INCLUDED

#include <opm/common/TimingMacros.hpp>

#include <opm/grid/utility/SparseTable.hpp>

#include <algorithm>
#include <cstddef>
#include <deque>
#include <limits>
#include <numeric>
#include <queue>
#include <string>
#include <tuple>
#include <vector>

namespace Opm {
namespace Detail {
template <class Graph>
std::size_t colorGraphWelshPowell(const Graph& graph,
                                  std::deque<typename Graph::VertexDescriptor>& orderedVertices,
                                  std::vector<int>& colors,
                                  int color,
                                  int noVertices)
{
    std::vector<int> forbidden(noVertices, false);
    std::size_t noColored = 0;

    for (auto vertex = orderedVertices.begin(),
              vertexEnd = orderedVertices.end(); vertex != vertexEnd; ++vertex) {
        // Skip forbidden vertices
        while (vertex != vertexEnd && forbidden[*vertex])
            ++vertex;
        if (vertex == vertexEnd) {
            break;
        }

        // Color Vertex
        colors[*vertex] = color;
        ++noColored;
        // Forbid neighors
        for (auto edge = graph.beginEdges(*vertex),
                  endEdge = graph.endEdges(*vertex); edge != endEdge; ++edge) {
            forbidden[edge.target()] = true;
        }
    }
    // forbidden vertices will be colored next for coloring
    using Vertex = typename Graph::VertexDescriptor;
    auto newEnd = std::remove_if(orderedVertices.begin(),
                                 orderedVertices.end(),
                                 [&forbidden](const Vertex& vertex) { return !forbidden[vertex]; });
    orderedVertices.resize(newEnd - orderedVertices.begin());
    return noColored;
}

template <class Graph, class Functor>
std::size_t breadthFirstSearch(const Graph& graph,
                               typename Graph::VertexDescriptor root,
                               Functor functor)
{
    std::vector<int> visited(graph.maxVertex() + 1, false);
    using Vertex = typename Graph::VertexDescriptor;
    std::queue<Vertex> nextVertices;
    std::size_t noVisited = 0;
    nextVertices.push(root);
    visited[root] = true; // We do not visit root.

    while (!nextVertices.empty()) {
        auto current = nextVertices.front();
        for (auto edge = graph.beginEdges(current),
                  endEdge = graph.endEdges(current); edge != endEdge; ++edge) {
            if (!visited[edge.target()]) {
                visited[edge.target()] = true;
                nextVertices.push(edge.target());
                functor(edge.target());
                ++noVisited;
            }
        }
        nextVertices.pop();
    }
    return noVisited;
}
} // end namespace Detail


/// \brief Color the vertices of graph.
///
/// It uses the algorithm of Welsh and Powell for this.
/// \param graph The graph to color. Must adhere to the graph interface of dune-istl.
/// \return A pair of a vector with the colors of the vertices and the number of colors
///         assigned
template <class Graph>
std::tuple<std::vector<int>, int, std::vector<std::size_t>>
colorVerticesWelshPowell(const Graph& graph)
{
    using Vertex = typename Graph::VertexDescriptor;
    std::deque<Vertex> orderedVertices;
    auto noVertices = graph.maxVertex() + 1;
    std::vector<int> degrees(noVertices, 0);
    int maxDegree = 0;
    std::ptrdiff_t firstDegreeChange = 0;

    // populate deque
    for (auto vertex = graph.begin(),
              endVertex = graph.end(); vertex != endVertex; ++vertex) {
        auto currentVertex = *vertex;
        auto& degree = degrees[currentVertex];

        for (auto edge = graph.beginEdges(currentVertex),
                  endEdge = graph.endEdges(currentVertex); edge != endEdge; ++edge) {
            ++degree;
        }

        if (degree >= maxDegree) {
            orderedVertices.emplace_front(currentVertex);
            ++firstDegreeChange;
            if (degree > maxDegree) {
                firstDegreeChange = 1;
                maxDegree = degree;
            }
        } else {
            orderedVertices.emplace_back(currentVertex);
        }
    }

    // order deque by descending degree
    std::stable_sort(orderedVertices.begin() + firstDegreeChange,
                     orderedVertices.end(),
                     [&degrees](const Vertex& v1, const Vertex& v2)
                     { return degrees[v1] > degrees[v2]; });

    // Overwrite degree with color
    auto& colors = degrees;
    std::fill(colors.begin(), colors.end(), -1);

    int color = 0;
    std::vector<std::size_t> verticesPerColor;
    verticesPerColor.reserve(10);

    while (!orderedVertices.empty()) {
        verticesPerColor.push_back(Detail::colorGraphWelshPowell(graph, orderedVertices,
                                                                 colors, color++, noVertices));
    }
    return std::make_tuple(colors, color, verticesPerColor);
}

/// \! Reorder colored graph preserving order of vertices with the same color.
template <class Graph>
std::vector<std::size_t>
reorderVerticesPreserving(const std::vector<int>& colors,
                          int noColors,
                          const std::vector<std::size_t>& verticesPerColor,
                          const Graph& graph)
{
    std::vector<std::size_t> colorIndex(noColors, 0);
    std::vector<std::size_t> indices(graph.maxVertex() + 1);
    std::partial_sum(verticesPerColor.begin(),
                     verticesPerColor.begin() + verticesPerColor.size() - 1,
                     colorIndex.begin() + 1);

    for (const auto& vertex : graph) {
        indices[vertex] = colorIndex[colors[vertex]]++;
    }
    return indices;
}

/// \! Reorder Vetrices in spheres
template <class Graph>
std::vector<std::size_t>
reorderVerticesSpheres(const std::vector<int>& colors,
                       int noColors,
                       const std::vector<std::size_t>& verticesPerColor,
                       const Graph& graph,
                       typename Graph::VertexDescriptor root)
{
    std::vector<std::size_t> colorIndex(noColors, 0);
    const auto notVisitedTag = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> indices(graph.maxVertex() + 1, notVisitedTag);
    using Vertex = typename Graph::VertexDescriptor;
    std::partial_sum(verticesPerColor.begin(),
                     verticesPerColor.begin() + verticesPerColor.size() - 1,
                     colorIndex.begin() + 1);
    std::size_t noVisited = 0;
    auto numberer = [&colorIndex, &colors, &indices](Vertex vertex)
    {
        indices[vertex] = colorIndex[colors[vertex]]++;
    };

    while (noVisited < graph.maxVertex() + 1) {
        numberer(root);
        ++noVisited; // root node already visited and not visited in BFS
        noVisited += Detail::breadthFirstSearch(graph, root, numberer);
        if (noVisited < graph.maxVertex() + 1) {
            // Graph is disconnected search for not yet visited node
            for (auto vertex : graph) {
                if (indices[vertex] == notVisitedTag) {
                    // \todo make sure that this is a peripheral node!
                    root = vertex;
                    break;
                }
            }
        }
    }
    return indices;
}

/// \brief Specify coloring type.
/// \details The coloring types have been implemented initially to parallelize DILU
///          preconditioner and parallel sparse triangular solves.
///          Symmetric coloring will create a dependency from row i to j if
///          both element A_ij and A_ji exists.
///          Lower coloring creates a dependency from row i to j (where i < j) if
///          A_ij is nonzero.
///          Upper coloring creates a dependecy from row i to j (where i > j) if A_ij is nonzero.
enum class ColoringType { SYMMETRIC, LOWER, UPPER };

/// This coloring algorithm interprets the sparsity structure of a matrix as a graph. Each row is given a color
/// or level where all the rows in the same level only have dependencies from lower levels. The level computation
/// is done with dynamic programming, and to improve caching the rows in the same level stay in matrix order.
/// \brief Given a matrix and dependecy type, returns a SparseTable grouping the rows by which
///        can be executed in parallel without breaking dependencies
/// \param matrix A dune sparse matrix
/// \param coloringType The coloringtype determines what constitutes a dependency,
///        see ColoringType definition above
/// \return SparseTable with rows of the matrix grouped into least number of groups
///         while dependencies only come from groups with lower index
template <class M>
Opm::SparseTable<std::size_t>
getMatrixRowColoring(const M& matrix, ColoringType coloringType)
{
    OPM_TIMEBLOCK(createMatrix);

    std::vector<std::size_t> color(matrix.N(), 0);
    std::vector<std::size_t> rowIndices(matrix.N(), 0);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (std::size_t i = 0; i < matrix.N(); i++) {
        rowIndices[i] = i;
    }

    std::vector<std::size_t> colorCnt;
    // These dynamic programming computations only rely on the following observation:
    // level[row_i] = 1 + max{level[row_j]} for all j which i depend on
    // This minimizes the level of each row, and every rows dependencies belong to a lower level set.
    if (coloringType == ColoringType::SYMMETRIC) {
        for (auto i = matrix.begin(); i != matrix.end(); ++i) {
            for (auto a_ij = i->begin(); a_ij.index() != i.index(); ++a_ij) {
                auto a_ji = matrix[a_ij.index()].find(i.index());
                if (a_ji != matrix[a_ij.index()].end()) {
                    color[i.index()] = std::max({color[i.index()], color[a_ij.index()] + 1});
                }
            }
            if (color[i.index()] >= colorCnt.size()) {
                colorCnt.push_back(1);
            } else {
                ++colorCnt[color[i.index()]];
            }
        }
    } else if (coloringType == ColoringType::UPPER) {
        for (auto i = matrix.beforeEnd(); i != matrix.beforeBegin(); --i) {
            for (auto a_ij = ++(i->find(i.index())); a_ij != i->end(); ++a_ij) {
                color[i.index()] = std::max({color[i.index()], color[a_ij.index()] + 1});
            }
            if (color[i.index()] >= colorCnt.size()) {
                colorCnt.push_back(1);
            } else {
                ++colorCnt[color[i.index()]];
            }
        }
    } else if (coloringType == ColoringType::LOWER) {
        for (auto i = matrix.begin(); i != matrix.end(); ++i) {
            for (auto a_ij = i->begin(); a_ij.index() != i.index(); ++a_ij) {
                color[i.index()] = std::max({color[i.index()], color[a_ij.index()] + 1});
            }
            if (color[i.index()] >= colorCnt.size()) {
                colorCnt.push_back(1);
            } else {
                ++colorCnt[color[i.index()]];
            }
        }
    }

    std::stable_sort(rowIndices.begin(),
                     rowIndices.end(),
                     [&color](const std::size_t a, const std::size_t b)
                     { return color[a] < color[b]; });

    return {rowIndices.data(), rowIndices.data() + rowIndices.size(),
            colorCnt.data(), colorCnt.data() + colorCnt.size()};
}

} // end namespace Opm

#endif
