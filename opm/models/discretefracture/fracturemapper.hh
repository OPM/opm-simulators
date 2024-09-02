// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::FractureMapper
 */
#ifndef EWOMS_FRACTURE_MAPPER_HH
#define EWOMS_FRACTURE_MAPPER_HH

#include <opm/models/utils/propertysystem.hh>

#include <algorithm>
#include <set>

namespace Opm {

/*!
 * \ingroup DiscreteFractureModel
 * \brief Stores the topology of fractures.
 */
template <class TypeTag>
class FractureMapper
{
    struct FractureEdge
    {
        FractureEdge(unsigned edgeVertex1Idx, unsigned edgeVertex2Idx)
            : i_(std::min(edgeVertex1Idx, edgeVertex2Idx)),
              j_(std::max(edgeVertex1Idx, edgeVertex2Idx))
        {}

        bool operator<(const FractureEdge& e) const
        { return i_ < e.i_ || (i_ == e.i_ && j_ < e.j_); }

        bool operator==(const FractureEdge& e) const
        { return i_ == e.i_ && j_ == e.j_; }

        unsigned i_;
        unsigned j_;
    };

public:
    /*!
     * \brief Constructor
     */
    FractureMapper()
    {}

    /*!
     * \brief Marks an edge as having a fracture.
     *
     * \param vertexIdx1 The index of the edge's first vertex.
     * \param vertexIdx2 The index of the edge's second vertex.
     */
    void addFractureEdge(unsigned vertexIdx1, unsigned vertexIdx2)
    {
        fractureEdges_.insert(FractureEdge(vertexIdx1, vertexIdx2));
        fractureVertices_.insert(vertexIdx1);
        fractureVertices_.insert(vertexIdx2);
    }

    /*!
     * \brief Returns true iff a fracture cuts through a given vertex.
     *
     * \param vertexIdx The index of the vertex.
     */
    bool isFractureVertex(unsigned vertexIdx) const
    { return fractureVertices_.count(vertexIdx) > 0; }

    /*!
     * \brief Returns true iff a fracture is associated with a given edge.
     *
     * \param vertex1Idx The index of the first vertex of the edge.
     * \param vertex2Idx The index of the second vertex of the edge.
     */
    bool isFractureEdge(unsigned vertex1Idx, unsigned vertex2Idx) const
    {
        FractureEdge tmp(vertex1Idx, vertex2Idx);
        return fractureEdges_.count(tmp) > 0;
    }

private:
    std::set<FractureEdge> fractureEdges_;
    std::set<unsigned> fractureVertices_;
};

} // namespace Opm

#endif // EWOMS_FRACTURE_MAPPER_HH
