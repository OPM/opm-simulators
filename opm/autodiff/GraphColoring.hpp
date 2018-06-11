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

#include <vector>
#include <deque>
#include <tuple>
#include <algorithm>

namespace Opm
{

namespace Detail
{
template<class Graph>
void colorGraphWelshPowell(const Graph& graph,
                           std::deque<typename Graph::VertexDescriptor>& orderedVertices,
                           std::vector<int>& colors,
                           int color, int noVertices)
{
    std::vector<int> forbidden(noVertices, false);

    for(auto vertex = orderedVertices.begin(),
            vertexEnd = orderedVertices.end();
        vertex != vertexEnd; ++vertex)
    {
        // Skip forbidden vertices
        while(vertex != vertexEnd && forbidden[*vertex])
            ++vertex;
        if ( vertex == vertexEnd )
        {
            break;
        }

        // Color Vertex
        colors[*vertex] = color;
        // Forbid neighors
        for(auto edge = graph.beginEdges(*vertex), endEdge = graph.endEdges(*vertex);
            edge != endEdge; ++edge)
        {
            forbidden[edge.target()] = true;
        }
    }
    // forbidden vertices will be colored next for coloring
    using Vertex = typename Graph::VertexDescriptor;
    auto newEnd = std::remove_if(orderedVertices.begin(), orderedVertices.end(),
                                 [&forbidden](const Vertex& vertex)
                                 {
                                     return !forbidden[vertex];
                                 });
    orderedVertices.resize(newEnd-orderedVertices.begin());
}
} // end namespace Detail


/// \brief Color the vertices of graph-
///
/// It uses the algorithm of Welsh and Powell for this.
/// \param graph The graph to color. Must adhere to the graph interface of dune-istl.
/// \return A pair of a vector with the colors of the vertices and the number of colors
///         assigned
template<class Graph>
std::tuple<std::vector<int>, int> colorVerticesWelshPowell(const Graph& graph)
{
    using Vertex = typename Graph::VertexDescriptor;
    std::deque<Vertex> orderedVertices;
    auto noVertices = graph.maxVertex()+1;
    std::vector<int> degrees(noVertices, 0);
    int maxDegree = 0;
    std::ptrdiff_t firstDegreeChange = 0;

    // populate deque
    for( auto vertex = graph.begin(), endVertex = graph.end();
         vertex != endVertex; ++vertex)
    {
        auto currentVertex = *vertex;
        auto& degree = degrees[currentVertex];

        for(auto edge = graph.beginEdges(currentVertex),
                endEdge = graph.endEdges(currentVertex);
            edge != endEdge; ++edge)
        {
            ++degree;
        }


        if( degree >= maxDegree )
        {
            orderedVertices.emplace_front(currentVertex);
            ++firstDegreeChange;
            if(degree > maxDegree)
            {
                firstDegreeChange = 1;
                maxDegree = degree;
            }
        }
        else
        {
            orderedVertices.emplace_back(currentVertex);
        }
    }

    // order deque by descending degree
    std::stable_sort(orderedVertices.begin() + firstDegreeChange,
                     orderedVertices.end(),
              [&degrees](const Vertex& v1, const Vertex& v2)
              {
                  return degrees[v1] > degrees[v2];
              });

    // Overwrite degree with color
    auto& colors = degrees;
    std::fill(colors.begin(), colors.end(), -1);

    int color = 0;

    while(!orderedVertices.empty())
    {
        Detail::colorGraphWelshPowell(graph, orderedVertices, colors, color++,
                                      noVertices);
    }
    return std::make_tuple(colors, color);
}
}
// end namespace Opm
#endif
