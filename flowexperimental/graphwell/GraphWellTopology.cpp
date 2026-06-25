/*
  Copyright 2026 SINTEF Digital, Mathematics and Cybernetics.

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
#include <flowexperimental/graphwell/GraphWellTopology.hpp>

#include <opm/input/eclipse/Schedule/MSW/Segment.hpp>
#include <opm/input/eclipse/Schedule/MSW/WellSegments.hpp>
#include <opm/input/eclipse/Schedule/Well/Connection.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

namespace Opm {

template<class Scalar>
void GraphWellTopology<Scalar>::finalize()
{
    const int nseg = numSegments();
    seg_connections_.assign(nseg, {});
    for (int c = 0; c < numConnections(); ++c) {
        const auto& conn = connections_[c];
        if (conn.up != surface_node)
            seg_connections_[conn.up].push_back(c);
        if (conn.down != surface_node)
            seg_connections_[conn.down].push_back(c);
    }

    seg_perforations_.assign(nseg, {});
    for (int p = 0; p < num_perforations_; ++p)
        seg_perforations_[perf_segment_[p]].push_back(p);
}

template<class Scalar>
GraphWellTopology<Scalar>
GraphWellTopology<Scalar>::fromWellSegments(const WellSegments& segs,
                                            const WellConnections& conns)
{
    GraphWellTopology top;
    const int nseg = static_cast<int>(segs.size());
    top.segments_.resize(nseg);
    for (int s = 0; s < nseg; ++s) {
        top.segments_[s].depth = static_cast<Scalar>(segs[s].depth());
        top.segments_[s].volume = static_cast<Scalar>(segs[s].volume());
    }

    // The top segment is the one without an outlet (outletSegment <= 0).
    int top_segment = 0;
    for (int s = 0; s < nseg; ++s) {
        if (segs[s].outletSegment() <= 0) {
            top_segment = s;
            break;
        }
    }

    // Surface connection (carries the control equation). Index 0 by convention.
    Connection surf;
    surf.up = surface_node;
    surf.down = top_segment;
    surf.kind = ConnectionKind::Surface;
    top.connections_.push_back(surf);
    top.surface_connection_ = 0;

    // One internal connection per segment that has an outlet.
    for (int s = 0; s < nseg; ++s) {
        const int outlet_number = segs[s].outletSegment();
        if (outlet_number <= 0)
            continue;
        const int outlet = segs.segmentNumberToIndex(outlet_number);

        Connection conn;
        conn.down = s;          // deeper / further from surface
        conn.up = outlet;       // nearer the surface
        conn.kind = ConnectionKind::Internal;
        conn.length = static_cast<Scalar>(segs[s].totalLength()
                                          - segs[outlet].totalLength());
        conn.diameter = static_cast<Scalar>(segs[s].internalDiameter());
        conn.area = static_cast<Scalar>(segs[s].crossArea());
        conn.roughness = static_cast<Scalar>(segs[s].roughness());
        conn.depth_diff = static_cast<Scalar>(segs[s].depth() - segs[outlet].depth());
        top.connections_.push_back(conn);
    }

    // Perforations: enumerate OPEN connections in deck order (matching the
    // production MSW perforation numbering) and map each to its segment.
    int perf = 0;
    for (std::size_t i = 0; i < conns.size(); ++i) {
        const Opm::Connection& c = conns.get(i);
        if (c.state() != Opm::Connection::State::OPEN)
            continue;
        const int seg = segs.segmentNumberToIndex(c.segment());
        top.perf_segment_.push_back(seg);
        top.perf_depth_diff_.push_back(static_cast<Scalar>(c.depth() - segs[seg].depth()));
        ++perf;
    }
    top.num_perforations_ = perf;

    top.finalize();
    return top;
}

template<class Scalar>
GraphWellTopology<Scalar>
GraphWellTopology<Scalar>::singleSegment(const Well& well)
{
    GraphWellTopology top;
    const Scalar ref_depth = static_cast<Scalar>(well.getRefDepth());
    top.segments_.resize(1);
    top.segments_[0].depth = ref_depth;
    top.segments_[0].volume = Scalar{0};   // accumulation negligible for a node-less well

    Connection surf;
    surf.up = surface_node;
    surf.down = 0;
    surf.kind = ConnectionKind::Surface;
    top.connections_.push_back(surf);
    top.surface_connection_ = 0;

    const auto& conns = well.getConnections();
    int perf = 0;
    for (std::size_t i = 0; i < conns.size(); ++i) {
        const Opm::Connection& c = conns.get(i);
        if (c.state() != Opm::Connection::State::OPEN)
            continue;
        top.perf_segment_.push_back(0);
        top.perf_depth_diff_.push_back(static_cast<Scalar>(c.depth() - ref_depth));
        ++perf;
    }
    top.num_perforations_ = perf;

    top.finalize();
    return top;
}

template class GraphWellTopology<double>;
#if FLOW_INSTANTIATE_FLOAT
template class GraphWellTopology<float>;
#endif

} // namespace Opm
