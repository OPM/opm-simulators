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

#ifndef OPM_GRAPHWELL_TOPOLOGY_HEADER_INCLUDED
#define OPM_GRAPHWELL_TOPOLOGY_HEADER_INCLUDED

#include <cstddef>
#include <vector>

namespace Opm {

class WellSegments;
class WellConnections;
class Well;

/// \brief Simple graph backbone for the GraphWell multisegment-well reformulation.
///
/// A well is described as a graph of \e segments (control volumes, carrying the
/// mass-conservation DOFs) connected by \e connections ("faces", carrying a single
/// total-flux DOF and the momentum / pressure-drop equation). One distinguished
/// \e surface connection attaches to the top segment and carries the well control
/// equation instead of a pressure-drop equation; its \c up node is the sentinel -1
/// (the surface, which has no mass equation).
///
/// Unlike the production MultisegmentWell, the structure here is a general graph
/// (adjacency lists, no single-outlet/tree assumption), so loops are representable.
///
/// Sign / orientation convention (used throughout assembly):
///   * Each connection has an \c up segment (nearer the wellhead/surface) and a
///     \c down segment (further from the surface, typically deeper).
///   * The connection flux DOF \c Q is the surface-volume rate in the \c down ->
///     \c up direction; it is positive for production (flow towards the surface).
///   * The component rate through the connection is q_c = Q * F(upwind), where the
///     upwind segment is \c down when Q >= 0 and \c up otherwise.
///   * The mass balance of segment \c s is written as accumulation plus net
///     outflow: the flux q_c leaves the \c down segment (+q_c in its residual) and
///     enters the \c up segment (-q_c). The helper \c sign() returns the geometric
///     incidence (+1 for up, -1 for down); the assembler uses the outflow signs
///     (the negation of that) directly.
template<class Scalar>
class GraphWellTopology
{
public:
    enum class ConnectionKind { Internal, Surface };

    /// Sentinel used for the surface node on the \c up side of the surface connection.
    static constexpr int surface_node = -1;

    struct Connection
    {
        int up{surface_node};   //!< segment index nearer the surface (surface_node for surface conn)
        int down{0};            //!< segment index further from the surface
        ConnectionKind kind{ConnectionKind::Internal};
        Scalar length{0};       //!< pipe length between the two segment nodes
        Scalar diameter{0};     //!< inner diameter
        Scalar area{0};         //!< cross-sectional area
        Scalar roughness{0};    //!< absolute roughness
        Scalar depth_diff{0};   //!< depth(down) - depth(up) (>= 0 when down is deeper)
    };

    struct SegmentGeom
    {
        Scalar depth{0};
        Scalar volume{0};
    };

    /// Build the graph from a parsed multisegment well (WELSEGS + COMPSEGS).
    static GraphWellTopology fromWellSegments(const WellSegments& segs,
                                              const WellConnections& conns);

    /// Build the trivial single-segment graph for a standard (non-MSW) well:
    /// one segment at the BHP reference depth, one surface connection, all
    /// perforations on that single segment.
    static GraphWellTopology singleSegment(const Well& well);

    int numSegments() const { return static_cast<int>(segments_.size()); }
    int numConnections() const { return static_cast<int>(connections_.size()); }
    int numPerforations() const { return num_perforations_; }

    const SegmentGeom& segment(int s) const { return segments_[s]; }
    const Connection& connection(int c) const { return connections_[c]; }

    //! Connections incident to segment \c s (both as up and as down end).
    const std::vector<int>& connectionsOfSegment(int s) const { return seg_connections_[s]; }

    //! Perforations attached to segment \c s (indices into the well's perforation list).
    const std::vector<int>& perforationsOfSegment(int s) const { return seg_perforations_[s]; }

    //! Segment index a perforation belongs to.
    int segmentOfPerforation(int perf) const { return perf_segment_[perf]; }

    //! depth(perf) - depth(segment) for the perforation.
    Scalar perforationDepthDiff(int perf) const { return perf_depth_diff_[perf]; }

    //! Index of the (unique) surface connection.
    int surfaceConnection() const { return surface_connection_; }

    //! Index of the top segment (the one carrying the surface connection).
    int topSegment() const { return connections_[surface_connection_].down; }

    //! +1 if \c s is the up end of \c c, -1 if it is the down end, 0 otherwise.
    int sign(int s, int c) const
    {
        const auto& conn = connections_[c];
        if (s == conn.up) return 1;
        if (s == conn.down) return -1;
        return 0;
    }

    //! Upwind segment for connection \c c given the (frozen) sign of its flux.
    //! For Q >= 0 flow is down->up so the upwind side is \c down; otherwise \c up.
    //! When the up node is the surface sentinel the only real segment (down) is used.
    int upwindSegment(int c, Scalar Q) const
    {
        const auto& conn = connections_[c];
        if (Q >= Scalar{0}) return conn.down;
        return conn.up == surface_node ? conn.down : conn.up;
    }

private:
    void finalize();   //!< build incidence lists once segments/connections are set

    std::vector<SegmentGeom> segments_;
    std::vector<Connection> connections_;
    std::vector<std::vector<int>> seg_connections_;
    std::vector<std::vector<int>> seg_perforations_;
    std::vector<int> perf_segment_;
    std::vector<Scalar> perf_depth_diff_;
    int num_perforations_{0};
    int surface_connection_{0};
};

} // namespace Opm

#endif // OPM_GRAPHWELL_TOPOLOGY_HEADER_INCLUDED
