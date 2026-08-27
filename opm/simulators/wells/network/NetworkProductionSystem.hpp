/*
  Copyright 2026 Equinor ASA.

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
#ifndef OPM_NETWORK_PRODUCTION_SYSTEM_HEADER_INCLUDED
#define OPM_NETWORK_PRODUCTION_SYSTEM_HEADER_INCLUDED

#include <opm/simulators/wells/network/NetworkSolve.hpp>

#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/simulators/wells/VFPHelpers.hpp>
#include <opm/simulators/wells/VFPProdProperties.hpp>
#include <opm/material/densead/Evaluation.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace Opm::NetworkSolve {

// ---------------------------------------------------------------------------
// Production networks
//
// The same formulation, with one thing different: a rate is three numbers
// instead of one. That turns out to be the whole of it. VFPPROD derives the
// water and gas fractions from the rate triple it is handed, so the fractions
// never become unknowns, and mixing at a node is then just a sum -- linear.
// The only nonlinearity is still the table lookup.
//
//   unknowns    pressure of every non-terminal node                     nP
//               three phase rates through every parent branch          3nP
//               three phase rates and a bhp for every well              4W
//
//   equations   branch drop  p_n - VFPPROD(thp = p_parent, q_n, alq)     nP
//               node balance per phase                                 3nP
//               inflow perf. per phase  q_wp - ipr_p(bhp_w)              3W
//               control      whichever of THP / BHP / ORAT is active      W
//
// Prototype: ALQ comes from the branch as it does today, guide rates and group
// targets are not handled, and re-routing (BRANPROP changing the tree) is not
// modelled -- the tree is taken as given.
template<class Scalar>
class ProductionSystem : public SystemBase<Scalar>
{
public:
    using State = std::vector<Scalar>;
    using ScalarType = Scalar;
    static constexpr int NP = 3;   // water, oil, gas -- the order VFPPROD wants

    /// Tied: on its rate limit and its tubing at once -- see the residual.
    /// Cmpl: all three of a well's own limits as one complementarity row.
    /// Shut: its tubing cannot lift at this node pressure; q = 0.
    enum class Control { Thp, Bhp, OilRate, Grup, Tied, Cmpl, Shut };

    struct Well
    {
        std::string name;
        int node = 0;
        int vfp_table = 0;
        Scalar alq = 0.0;
        /// Per phase, production positive: q_p = ipr_a[p] + ipr_b[p] * bhp.
        std::array<Scalar, NP> ipr_a{};
        std::array<Scalar, NP> ipr_b{};
        Scalar bhp_limit = 0.0;
        Scalar oil_rate_limit = 0.0;
        /// Held by the group, so its rate counts against the target whatever
        /// control it ends on.
        bool in_group = false;
        Scalar guide = 0.0;
        /// WEFAC as it applies to the network: the branch sees efficiency * q.
        Scalar efficiency = 1.0;
        /// Hydrostatic correction between the tubing table's datum and the
        /// well's reference depth: the well's bhp is the table's less this.
        Scalar vfp_dp = 0.0;
        /// Held at oil_rate_limit, whatever the node pressure: a well the
        /// network is not deciding for is a source, and offering it thp lets
        /// it undercut the rate it was given. Sets the control to OilRate.
        bool pinned = false;
        /// A well the well model has at zero rate at this thp is dead, and more
        /// back-pressure cannot revive it: its tubing allows nothing at or
        /// above this pressure. Zero means not dead.
        Scalar dead_above = 0.0;
        /// Gas the well is lifted with. It goes up the tubing and so into the
        /// branch's gas stream, but it is not produced, so it is not in q. Zero
        /// unless the node is set to add it (NODEPROP item 4).
        Scalar lift_gas = 0.0;
        /// The oil rate the well model has now. A tubing table with a loading
        /// hump gives a well two branches, dead and flowing, and which one the
        /// solve lands on depends on where it starts; this says which it is on.
        Scalar q_start = 0.0;
        /// Shut by the adapter: no inflow performance, or shut by the network
        /// earlier in this step. No IPR; the rows hold q = 0, bhp = limit.
        bool shut = false;
        /// Whether the node adds the well's lift gas to the stream, so a
        /// change of alq is a change of lift_gas too.
        bool node_adds_lift_gas = false;
    };

    ProductionSystem(const VFPProdProperties<Scalar>& props, const UnitSystem& units)
        : props_(&props), units_(&units)
    {}

    void addNode(Node n, const Scalar alq)
    {
        nodes_.push_back(std::move(n));
        branch_alq_.push_back(alq);
        node_source_.push_back({});
        node_choke_target_.push_back(Scalar{0});
        node_choked_.push_back(0);
        node_choke_pressure_.push_back(Scalar{0});
    }
    /// A rate that enters at a node without a well behind it -- satellite
    /// production, or lift gas that is not any well's. Water, oil, gas.
    void setNodeSource(const int node, const std::array<Scalar, NP>& q) { node_source_[node] = q; }

    /// Make a node an autochoke: a valve just upstream of it that throttles
    /// the oil collected there to `target`. The node's pressure is then the
    /// group's common thp -- raised above the upstream pressure until the oil
    /// through it is the target, or left at the upstream pressure when even
    /// that does not reach the target (the choke is open). Oil-rate targets
    /// only, for now.
    void setChokeTarget(const int node, const Scalar target) { node_choke_target_[node] = target; }
    Scalar groupTarget() const { return group_target_; }
    Scalar branchAlq(const int node) const { return branch_alq_[node]; }
    const std::array<Scalar, NP>& nodeSource(const int node) const { return node_source_[node]; }
    bool isChoke(const int node) const { return node_choke_target_[node] > Scalar{0}; }
    bool choked(const int node) const { return node_choked_[node] != 0; }
    /// The pressure the last control selection placed a closed choke at.
    Scalar chokePressure(const int node) const { return node_choke_pressure_[node]; }

    /// What a well's own controls allow it at node pressure p: the least of its
    /// bhp limit, its rate limit and its tubing; zero if the tubing cannot lift.
    Scalar wellAllowance(const Well& well, const Scalar p) const
    {
        if (well.pinned) {
            return well.oil_rate_limit;
        }
        if (well.dead_above > Scalar{0} && p >= well.dead_above) {
            return Scalar{0};
        }
        Scalar allow = ipr(well, 1, well.bhp_limit);
        if (well.oil_rate_limit > Scalar{0}) {
            allow = std::min(allow, well.oil_rate_limit);
        }
        if (hasTubing(well)) {
            const Scalar found = cachedThpPotential(well, p);
            allow = (found > Scalar{0}) ? std::min(allow, found) : Scalar{0};
        }
        return allow;
    }

    /// thpPotential() is a scan of the tubing table, and the choke's root-find
    /// asks for it at dozens of pressures per well per control selection. The
    /// answer is a smooth function of p away from the loading cliff, so keep it
    /// on a grid and interpolate; anything within a hair of a grid point is the
    /// grid point. Cleared whenever a well's data change (finish()).
    Scalar cachedThpPotential(const Well& well, const Scalar p) const
    {
        const int w = static_cast<int>(&well - wells_.data());
        auto& grid = potential_grid_[w];
        constexpr Scalar step = Scalar{0.25} * unit::barsa;   // finer than any table feature
        const long key = static_cast<long>(std::floor(p / step));
        const auto lo = grid.find(key);
        const auto hi = grid.find(key + 1);
        const Scalar p_lo = key * step, p_hi = (key + 1) * step;
        const Scalar v_lo = (lo != grid.end()) ? lo->second
                          : grid.emplace(key, thpPotential(well, p_lo)).first->second;
        const Scalar v_hi = (hi != grid.end()) ? hi->second
                          : grid.emplace(key + 1, thpPotential(well, p_hi)).first->second;
        // A cliff is not interpolated across: one side cannot lift or is
        // unbounded, or the two differ by more than the tubing could over a
        // quarter bar -- which is the stable crossing jumping along the hump.
        // Those get the exact answer; the smooth stretches get the grid.
        constexpr Scalar inf = std::numeric_limits<Scalar>::max();
        if (!(v_lo > Scalar{0}) || !(v_hi > Scalar{0}) || v_lo == inf || v_hi == inf
            || std::abs(v_hi - v_lo) > Scalar{0.05} * std::max(v_lo, v_hi)) {
            return thpPotential(well, p);
        }
        const Scalar t = (p - p_lo) / step;
        return v_lo + t * (v_hi - v_lo);
    }

    /// Whether the well has no crossing at this node pressure, and so produces
    /// nothing. Sticky in one direction only: a well seen alive lower down that
    /// loses its crossing as the pressure rises stays shut for the solve. It has
    /// to be -- shutting it lowers the node pressure, which revives it, which
    /// raises the pressure again, and the system has no fixed point. A start
    /// that is merely far too high is not a death.
    bool cmplDead(const int w, const Scalar p) const
    {
        if (cmpl_dead_[w]) { return true; }
        const bool alive = cachedThpPotential(wells_[w], p) > Scalar{0};
        if (alive) {
            cmpl_alive_p_[w] = std::min(cmpl_alive_p_[w], p);
            return false;
        }
        if (p > cmpl_alive_p_[w]) { cmpl_dead_[w] = 1; }
        return true;
    }
    /// Whether anything under the node answers to its pressure: a well the
    /// well model has pinned does not, and a valve over pinned wells only has
    /// nothing to throttle -- its row would have no unknown to act on.
    bool chokeCanAct(const int node) const
    {
        for (const int w : wells_at_[node]) {
            if (!wells_[w].pinned && hasTubing(wells_[w])) { return true; }
        }
        for (const int c : children_[node]) {
            if (chokeCanAct(c)) { return true; }
        }
        return false;
    }
    Scalar cmplAllowance(const Well& well, const Scalar scan) const
    {
        Scalar allow = std::min(scan, ipr(well, 1, well.bhp_limit));
        if (well.oil_rate_limit > Scalar{0}) { allow = std::min(allow, well.oil_rate_limit); }
        return allow;
    }

    /// Oil the node would collect with its valve open at pressure p: the wells'
    /// allowances, the sources, and the children as the iterate has them.
    Scalar chokeDeliverable(const int node, const Scalar p, const State& x) const
    {
        Scalar total = node_source_[node][1];
        for (const int c : children_[node]) {
            total += nodes_[c].efficiency * x[qIdx(c, 1)];
        }
        for (const int w : wells_at_[node]) {
            total += wells_[w].efficiency * wellAllowance(wells_[w], p);
        }
        return total;
    }
    Scalar chokeTarget(const int node) const { return node_choke_target_[node]; }
    void addWell(Well w) { wells_.push_back(std::move(w)); }
    void setTerminalPressure(const Scalar p) { terminal_pressure_ = p; }
    Scalar terminalPressure() const { return terminal_pressure_; }
    void setRateScale(const Scalar s) { rate_scale_ = s; }
    /// Oil rate the group above this network is asked for.
    void setGroupTarget(const Scalar target) { group_target_ = target; }

    bool grouped() const { return group_target_ > 0.0; }

    void finish()
    {
        children_.assign(nodes_.size(), {});
        wells_at_.assign(nodes_.size(), {});
        for (std::size_t n = 1; n < nodes_.size(); ++n) {
            children_[nodes_[n].parent].push_back(static_cast<int>(n));
        }
        for (std::size_t w = 0; w < wells_.size(); ++w) {
            wells_at_[wells_[w].node].push_back(static_cast<int>(w));
        }
        for (auto& w : wells_) {
            if (w.guide <= 0.0) {
                w.guide = std::max(w.oil_rate_limit, Scalar{1.0});
            }
        }
        if (rate_scale_ <= 0.0) {
            Scalar largest = group_target_;
            for (const auto& w : wells_) {
                largest = std::max(largest, w.oil_rate_limit);
            }
            rate_scale_ = std::max(largest * Scalar{0.01},
                                   unit::convert::from(1.0, unit::cubic(unit::meter) / unit::day));
        }
        potential_grid_.assign(wells_.size(), {});
        recent_.assign(wells_.size(), {Control::Thp, Control::Thp});
        cmpl_.assign(wells_.size(), 0);
        cmpl_dead_.assign(wells_.size(), 0);
        cmpl_alive_p_.assign(wells_.size(), std::numeric_limits<Scalar>::max());
        cmpl_decided_ = false;
        controls_.assign(wells_.size(), Control::Thp);
        for (std::size_t w = 0; w < wells_.size(); ++w) {
            if (!hasTubing(wells_[w])) {
                controls_[w] = Control::Bhp;
            }
            if (grouped() && wells_[w].in_group) {
                controls_[w] = Control::Grup;
            }
            if (wells_[w].pinned) {
                controls_[w] = Control::OilRate;
            }
        }
    }

    int numNodes() const { return static_cast<int>(nodes_.size()) - 1; }
    int numWells() const override { return static_cast<int>(wells_.size()); }
    int size() const override { return 4 * numNodes() + 4 * numWells() + 1; }

    const std::vector<Node>& nodes() const { return nodes_; }
    const std::vector<Well>& wells() const { return wells_; }
    Control control(const int w) const { return controls_[w]; }

    char controlLetter(const int w) const override
    {
        switch (controls_[w]) {
        case Control::Thp:     return 'T';
        case Control::Bhp:     return 'B';
        case Control::OilRate: return 'O';
        case Control::Grup:    return 'G';
        case Control::Tied:    return 'C';
        case Control::Cmpl:    return 'M';
        case Control::Shut:    return 'S';
        }
        return '?';
    }

    int pIdx(const int node) const { return node - 1; }
    int qIdx(const int node, const int ph) const { return numNodes() + NP * (node - 1) + ph; }
    int qwIdx(const int w, const int ph) const { return 4 * numNodes() + NP * w + ph; }
    int bhpIdx(const int w) const { return 4 * numNodes() + NP * numWells() + w; }
    int lambdaIdx() const { return 4 * numNodes() + 4 * numWells(); }

    bool hasTable(const Node& n) const { return n.vfp_table != NoTable; }

    static Scalar ipr(const Well& w, const int ph, const Scalar bhp)
    {
        return w.ipr_a[ph] + w.ipr_b[ph] * bhp;
    }

    /// Pressure below a branch carrying these phase rates. Production rates are
    /// positive here and negative to the table, as the relaxed path does.
    Scalar tableBhp(const int table, const Scalar thp,
                    const std::array<Scalar, NP>& q, const Scalar alq) const
    {
        return props_->bhp(table, -q[0], -q[1], -q[2], thp, alq,
                           Scalar{0}, Scalar{0}, /*use_expvfp=*/false);
    }

    /// Whether thp control is even available: a well the deck gives no VFPPROD
    /// table has no tubing curve, so its rate does not answer to the node
    /// pressure and thp is not one of its controls.
    static bool hasTubing(const Well& w) { return w.vfp_table > 0; }

    /// The oil rate thp control allows at this node pressure -- the *crossing*
    /// of the well's inflow performance with its tubing curve. Zero means the
    /// tubing cannot lift here at all; max() means thp does not bind.
    ///
    /// Searched in bhp rather than in rate, because the tubing lookup wants the
    /// whole phase triple and the inflow performance gives it from a bhp
    /// directly, so no fractions have to be guessed at.
    ///
    /// Two shapes to know about, both referred to elsewhere in this file:
    ///   - the *crossing* is where h(bhp) = bhp - tableBhp(...) turns positive;
    ///   - the *liquid-loading hump* makes h non-monotone. At low rates a tubing
    ///     table needs more pressure than at moderate rates, so h can be negative
    ///     at both ends of the bracket and positive between, with two crossings.
    ///     Only the one where h turns positive with rising bhp is an operating
    ///     point; the other is the loading point.
    ///
    /// TODO: this should not be a search. With constant phase fractions the IPR
    /// is linear in FLO and the table piecewise-linear on its own flow axis, so
    /// the crossing is exact per interval -- VFPHelpers::intersectWithIPR, which
    /// estimateStableBhp already uses. That is ~21 lookups against this scan's
    /// ~136, and exact: the scan's resolution is why cachedThpPotential() below
    /// needs a jump guard.

    Scalar thpPotential(const Well& w, const Scalar p_node) const
    {
        if (!hasTubing(w) || !(w.ipr_b[1] < Scalar{0})) {
            return Scalar{0};
        }
        // The bhp at which the *last* phase stops, not the one at which oil
        // does. Water and gas have their own zero crossings, and at the oil
        // one they are still flowing -- so the tubing still has something to
        // lift there, the bracket does not contain the crossing, and the well
        // reads as unable to produce at all.
        Scalar shut = w.bhp_limit;
        for (int ph = 0; ph < NP; ++ph) {
            if (w.ipr_b[ph] < Scalar{0}) {
                shut = std::max(shut, -w.ipr_a[ph] / w.ipr_b[ph]);
            }
        }
        Scalar lo = w.bhp_limit;
        if (!(shut > lo)) {
            return Scalar{0};
        }
        auto rates = [&](const Scalar bhp) {
            std::array<Scalar, NP> q{};
            for (int ph = 0; ph < NP; ++ph) {
                q[ph] = std::max(ipr(w, ph, bhp), Scalar{0});
            }
            return q;
        };
        // Start where the table starts to know the answer. Below that bhp the
        // IPR rate is beyond the flow axis and the extrapolated tubing curve
        // has sign changes of its own, which the scan took for crossings: the
        // same well read 6225, 0 and 9368 m3/d within a tenth of a bar.
        {
            const auto& t = props_->getTable(w.vfp_table);
            auto flo = [&](const Scalar bhp) {
                const auto q = rates(bhp);
                return std::abs(detail::getFlo(t, -q[0], -q[1], -q[2]));
            };
            const Scalar fmax = t.getFloAxis().back();
            const Scalar f0 = flo(lo), f1 = flo(shut);
            if (f0 > fmax && f0 > f1) {
                lo += (shut - lo) * (f0 - fmax) / (f0 - f1);
                lo = std::min(lo, shut);
            }
        }
        auto h = [&](const Scalar bhp) {
            return bhp - (tableBhp(w.vfp_table, p_node, rates(bhp), w.alq) - w.vfp_dp);
        };
        // At the bhp limit the tubing already needs less than the limit, so thp
        // does not hold the well back: it is the bhp limit that binds. Say so by
        // allowing more than the limit does -- reporting exactly what the limit
        // allows makes a tie, the tie goes to thp, and the thp row then settles
        // the bhp *below* its limit.
        if (h(lo) >= Scalar{0}) {
            return std::numeric_limits<Scalar>::max();
        }
        // h is not monotone. At low rates a tubing table typically needs *more*
        // pressure than at moderate rates -- liquid loading -- so h can be
        // negative at both ends of the bracket and positive in between, with
        // two crossings. A bisection on the ends sees "cannot lift" for a well
        // that lifts perfectly well. Scan instead, and take the crossing where
        // h turns positive with rising bhp: more drawdown there means more
        // excess pressure, which is the one the well settles on. The other is
        // the loading point, and not an operating point.
        constexpr int samples = 96;
        Scalar a = lo, b = shut;
        bool found = false;
        Scalar h_prev = h(lo);
        for (int i = 1; i <= samples; ++i) {
            const Scalar bhp_i = lo + (shut - lo) * Scalar(i) / Scalar(samples);
            const Scalar h_i = h(bhp_i);
            if (h_prev < Scalar{0} && h_i >= Scalar{0}) {
                a = lo + (shut - lo) * Scalar(i - 1) / Scalar(samples);
                b = bhp_i;
                found = true;
                break;
            }
            h_prev = h_i;
        }
        if (!found) {
            return Scalar{0};   // nowhere does the reservoir out-push the tubing
        }
        Scalar bhp = b;
        for (int it = 0; it < 40; ++it) {
            bhp = Scalar{0.5} * (a + b);
            (h(bhp) < Scalar{0} ? a : b) = bhp;
        }
        return std::max(ipr(w, 1, bhp), Scalar{0});
    }

    std::vector<Scalar> guides() const
    {
        std::vector<Scalar> g(wells_.size());
        std::transform(wells_.begin(), wells_.end(), g.begin(),
                       [](const Well& w) { return w.guide; });
        return g;
    }

    std::vector<char> inGroup() const
    {
        std::vector<char> in(wells_.size());
        std::transform(wells_.begin(), wells_.end(), in.begin(),
                       [](const Well& w) { return static_cast<char>(w.in_group); });
        return in;
    }

    State residual(const State& x) const override
    {
        const int nodes = numNodes();
        const int wells = numWells();
        State r(size(), 0.0);
        auto pressure = [&](const int n) { return n == 0 ? terminal_pressure_ : x[pIdx(n)]; };
        auto branchRates = [&](const int n) {
            std::array<Scalar, NP> q{};
            for (int ph = 0; ph < NP; ++ph) {
                q[ph] = x[qIdx(n, ph)];
            }
            return q;
        };

        for (int n = 1; n <= nodes; ++n) {
            const auto& node = nodes_[n];
            const Scalar upstream = pressure(node.parent);
            if (isChoke(n) && complementarity_ && analytic_jacobian_ && chokeCanAct(n)) {
                // The valve's drop, on top of whatever the branch itself loses,
                // p - branch(p_up, q) >= 0, and the group's surplus
                // target - q >= 0, one of them zero: open and under target, or
                // throttling at the target. One row, no state to decide.
                const Scalar open_p = hasTable(node)
                    ? tableBhp(node.vfp_table, upstream, branchRates(n), branch_alq_[n]) : upstream;
                const Scalar drop = (x[pIdx(n)] - open_p) / pressure_scale_;
                const Scalar surplus = (node_choke_target_[n] - x[qIdx(n, 1)]) / rate_scale_;
                r[n - 1] = fb(drop, surplus);
            } else if (isChoke(n) && choked(n) && !(complementarity_ && analytic_jacobian_)) {
                // The valve holds the oil through the node at the target; the
                // node pressure is whatever that takes. No table on a choke
                // branch -- the drop *is* the unknown.
                r[n - 1] = (x[qIdx(n, 1)] - node_choke_target_[n]) / rate_scale_;
            } else {
                r[n - 1] = (hasTable(node)
                    ? x[pIdx(n)] - tableBhp(node.vfp_table, upstream, branchRates(n), branch_alq_[n])
                    : x[pIdx(n)] - upstream) / pressure_scale_;
            }

            for (int ph = 0; ph < NP; ++ph) {
                Scalar balance = x[qIdx(n, ph)] - node_source_[n][ph];
                for (const int c : children_[n]) {
                    balance -= nodes_[c].efficiency * x[qIdx(c, ph)];
                }
                for (const int w : wells_at_[n]) {
                    balance -= wells_[w].efficiency * x[qwIdx(w, ph)];
                    if (ph == 2) {
                        balance -= wells_[w].efficiency * wells_[w].lift_gas;
                    }
                }
                r[nodes + NP * (n - 1) + ph] = balance / rate_scale_;
            }
        }

        Scalar produced = 0.0;
        for (int w = 0; w < wells; ++w) {
            const auto& well = wells_[w];
            const Scalar bhp = x[bhpIdx(w)];
            std::array<Scalar, NP> q{};
            for (int ph = 0; ph < NP; ++ph) {
                q[ph] = x[qwIdx(w, ph)];
                r[4 * nodes + NP * w + ph] = (q[ph] - ipr(well, ph, bhp)) / rate_scale_;
            }
            // Every well the group allocated counts against the target, on
            // whatever control it ended on.
            if (well.in_group) {
                produced += q[1];
            }
            Scalar& control = r[4 * nodes + NP * wells + w];
            switch (controls_[w]) {
            case Control::Thp:
                control = (bhp - (tableBhp(well.vfp_table, pressure(well.node), q, well.alq)
                                  - well.vfp_dp)) / pressure_scale_;
                break;
            case Control::Bhp:
                control = (bhp - well.bhp_limit) / pressure_scale_;
                break;
            case Control::OilRate:
                control = (q[1] - well.oil_rate_limit) / rate_scale_;
                break;
            case Control::Grup:
                control = (q[1] - well.guide * x[lambdaIdx()]) / rate_scale_;
                break;
            case Control::Cmpl: {
                // All three of the well's own limits as one row. a, b, c are
                // the slacks on rate, tubing and bhp; at the solution each is
                // >= 0 and at least one is 0, which is exactly fb()'s zero set.
                // Nothing is selected here, so nothing can flip.
                // Whether the tubing lifts at all is thpPotential()'s answer,
                // not b's: b is negative below the liquid-loading hump too,
                // where the well does flow. q >= 0 is kept by limitStep(), so
                // the lookup takes max(q, 0).
                if (cmplDead(w, pressure(well.node))) {
                    control = q[1] / rate_scale_;
                    break;
                }
                std::array<Scalar, NP> qp{};
                for (int ph = 0; ph < NP; ++ph) { qp[ph] = std::max(q[ph], Scalar{0}); }
                const Scalar a = (well.oil_rate_limit > Scalar{0})
                    ? (well.oil_rate_limit - q[1]) / rate_scale_ : Scalar{1e6};
                const Scalar b = (bhp - (tableBhp(well.vfp_table, pressure(well.node), qp, well.alq)
                                         - well.vfp_dp)) / pressure_scale_;
                const Scalar c = (bhp - well.bhp_limit) / pressure_scale_;
                control = fb(fb(a, b), c);
                break;
            }
            case Control::Shut:
                control = (well.ipr_b[1] < Scalar{0}) ? q[1] / rate_scale_
                                                      : (bhp - well.bhp_limit) / pressure_scale_;
                break;
            case Control::Tied: {
                // The well is at its rate limit with the tubing only just
                // passing it: rate slack a = limit - q and tubing slack
                // b = bhp - tubing(p_node, q) are both non-negative and one of
                // them is zero. Choosing which by the sign at the iterate is
                // what flips; the Fischer-Burmeister function has the same
                // zero set and is smooth everywhere but the origin.
                const Scalar a = (well.oil_rate_limit - q[1]) / rate_scale_;
                const Scalar b = (bhp - (tableBhp(well.vfp_table, pressure(well.node), q, well.alq)
                                         - well.vfp_dp)) / pressure_scale_;
                control = a + b - std::sqrt(a * a + b * b);
                break;
            }
            }
        }

        // With nobody on group control the multiplier is free, so pin it rather
        // than hand the Newton a singular column.
        const bool any = std::find(controls_.begin(), controls_.end(), Control::Grup)
                       != controls_.end();
        r[lambdaIdx()] = grouped() && any ? (produced - group_target_) / rate_scale_
                                          : (x[lambdaIdx()] - lambda0()) / rate_scale_;
        return r;
    }

    /// Each control names the oil rate it would allow at this node pressure and
    /// the smallest wins, exactly as on the injection side -- a producer is just
    /// held back by a bhp that is too *low* rather than too high, so the bhp
    /// limit's allowance is the rate the inflow gives at it.
    ///
    /// Nothing here reads the iterate's rates, bhp or multiplier: mid-Newton
    /// those are not a consistent well state, on rate control q *is* the limit,
    /// and the multiplier is defined by the very active set being chosen here.
    bool updateControls(const State& x) override
    {
        const int n = numWells();
        constexpr Scalar unbounded = std::numeric_limits<Scalar>::max();

        // A choke closes when the wells behind it could deliver more than the
        // target with the valve open -- at the upstream pressure. Once closed,
        // the wells' controls have to be chosen at the pressure the choke will
        // settle at, not at the iterate's: at the upstream pressure the tubing
        // allows more than each well's own limit, every well picks its limit,
        // and the choke row is then left with no unknown to act on. So find
        // that pressure here, from the allowances alone -- the same cheap
        // model the control rule runs on -- and let the Newton only polish it.
        bool changed = false;
        std::vector<Scalar> choke_pressure(numNodes() + 1, Scalar{0});
        for (int node = 1; node <= numNodes(); ++node) {
            if (!isChoke(node) || (complementarity_ && analytic_jacobian_)) {
                continue;            // the row decides; nothing to pre-solve
            }
            const Scalar p_up = (nodes_[node].parent == 0) ? terminal_pressure_
                                                           : x[pIdx(nodes_[node].parent)];
            auto deliverable = [&](const Scalar p) { return chokeDeliverable(node, p, x); };
            const Scalar target = node_choke_target_[node];
            const char want = (deliverable(p_up) > target) ? 1 : 0;
            changed |= (want != node_choked_[node]);
            node_choked_[node] = want;
            choke_pressure[node] = p_up;
            if (want) {
                // Deliverability falls with pressure; bracket upward, then bisect.
                Scalar lo = p_up, hi = p_up;
                for (int it = 0; it < 40 && deliverable(hi) > target; ++it) {
                    lo = hi;
                    hi = hi * Scalar{1.25} + unit::barsa;
                }
                for (int it = 0; it < 50; ++it) {
                    const Scalar mid = Scalar{0.5} * (lo + hi);
                    (deliverable(mid) > target ? lo : hi) = mid;
                }
                choke_pressure[node] = Scalar{0.5} * (lo + hi);
            }
            node_choke_pressure_[node] = choke_pressure[node];
        }

        std::vector<Scalar> own(n), thp(n, unbounded);
        for (int w = 0; w < n; ++w) {
            const auto& well = wells_[w];
            const Scalar p_node = (well.node == 0) ? terminal_pressure_
                : (isChoke(well.node) && choked(well.node) && !(complementarity_ && analytic_jacobian_))
                    ? choke_pressure[well.node]
                : x[pIdx(well.node)];
            // Zero from thpPotential() means the well cannot lift against this
            // node pressure at all -- its table does not reach that high, or the
            // inflow cannot feed the tubing. That makes thp *unavailable*, not a
            // control that allows nothing: taken as an allowance of zero it wins
            // "most restrictive" every time, and the thp row it then imposes
            // says nothing about a rate, so the well produces whatever the
            // tubing crossing happens to be.
            if (hasTubing(well) && !well.pinned
                && !(well.dead_above > Scalar{0} && p_node >= well.dead_above)) {
                const Scalar found = cachedThpPotential(well, p_node);
                if (found > Scalar{0}) {
                    thp[w] = found;
                }
            }
            own[w] = std::min(thp[w], ipr(well, 1, well.bhp_limit));
            if (well.oil_rate_limit > Scalar{0}) {
                own[w] = std::min(own[w], well.oil_rate_limit);
            }
            if (well.pinned) {
                own[w] = well.oil_rate_limit;
            }
        }

        const auto share = shareByGuide(guides(), inGroup(), own, group_target_);

        for (int w = 0; w < n; ++w) {
            const auto& well = wells_[w];
            if (well.pinned) {
                changed |= (controls_[w] != Control::OilRate);
                controls_[w] = Control::OilRate;
                continue;
            }
            // A bhp limit at or above the shut-in pressure produces nothing;
            // the Bhp row would sit the well there and report injection.
            if (!(ipr(well, 1, well.bhp_limit) > Scalar{0})) {
                changed |= (controls_[w] != Control::Shut);
                controls_[w] = Control::Shut;
                continue;
            }
            auto wanted = (thp[w] < unbounded) ? Control::Thp : Control::Bhp;
            Scalar smallest = thp[w];
            Scalar current_allows = unbounded;
            auto consider = [&](const Control c, const Scalar allows) {
                if (allows < smallest) {
                    smallest = allows;
                    wanted = c;
                }
                if (c == controls_[w]) {
                    current_allows = allows;
                }
            };
            if (controls_[w] == Control::Thp) {
                current_allows = thp[w];
            }
            consider(Control::Bhp, ipr(well, 1, well.bhp_limit));
            if (well.oil_rate_limit > Scalar{0}) {
                consider(Control::OilRate, well.oil_rate_limit);
            }
            if (grouped() && well.in_group) {
                consider(Control::Grup, share[w]);
            }
            // A control is overtaken by a margin, not by rounding. A well that
            // sits exactly where its tubing passes its own rate limit -- which
            // is where a choke puts the marginal well -- would otherwise flip
            // between the two every iteration.
            if (wanted != controls_[w] && current_allows < unbounded
                && current_allows <= smallest * (Scalar{1} + Scalar{1e-3})) {
                wanted = controls_[w];
            }
            if (complementarity_ && analytic_jacobian_) {
                // Decided once per solve, from the starting point, and kept:
                // deciding it from the iterate makes the well switch between
                // the row and the active set as the node pressure moves, which
                // is the switching the row exists to remove. A well that turns
                // out unable to lift leaves the row unsatisfied and the solve
                // is handed back, honestly, rather than re-decided mid-way.
                if (!cmpl_decided_) {
                    const bool dead = well.dead_above > Scalar{0}
                        && ((well.node == 0 ? terminal_pressure_ : x[pIdx(well.node)]) >= well.dead_above);
                    // Not conditioned on a finite tubing allowance at the start:
                    // the row itself covers "thp does not bind" (b large) and
                    // "cannot lift" (the q = 0 branch), and a start pressure
                    // far from the answer must not decide a well's row.
                    // A well the well model has at zero rate at this pressure
                    // is shut, not handed to the active set, which would put
                    // it on bhp or rate control and produce through tubing
                    // that cannot lift (three dead wells at 2500 m3/d each
                    // under a 6000 choke target: dump_prod_1112).
                    cmpl_[w] = well.shut ? 2
                        : hasTubing(well) && !well.pinned && !(grouped() && well.in_group)
                        ? (dead ? 2 : 1) : 0;
                }
            }
            if (complementarity_ && analytic_jacobian_ && cmpl_[w] == 2) {
                changed |= (controls_[w] != Control::Shut);
                controls_[w] = Control::Shut;
                continue;
            }
            if (complementarity_ && analytic_jacobian_ && cmpl_[w]) {
                wanted = Control::Cmpl;           // the row decides; nothing to switch
            } else if (controls_[w] == Control::Tied) {
                wanted = Control::Tied;           // sticky: the row decides
            } else if (analytic_jacobian_ && wanted != controls_[w] && wanted == recent_[w][0]
                       && (wanted == Control::OilRate || controls_[w] == Control::OilRate)
                       && well.oil_rate_limit > Scalar{0} && hasTubing(well)) {
                // Only with the assembled Jacobian: the row is not smooth at
                // the origin, and a difference that straddles the kink is not
                // a derivative of anything.
                wanted = Control::Tied;
            }
            recent_[w] = {recent_[w][1], controls_[w]};
            changed |= (wanted != controls_[w]);
            controls_[w] = wanted;
        }
        cmpl_decided_ = true;
        return changed;
    }

    State start(const State& node_pressure) const override
    {
        State x(size(), 0.0);
        for (int n = 1; n <= numNodes(); ++n) {
            x[pIdx(n)] = node_pressure[n];
        }
        for (int w = 0; w < numWells(); ++w) {
            const auto& well = wells_[w];
            // Start each well at the oil rate its own controls allow at the
            // guessed node pressure -- the rate the control rule will pick --
            // and the bhp that gives it. Opening at the bhp limit instead puts
            // a well whose limit is the 1 atm default at a rate off the end of
            // every table, and the Newton never comes back from there.
            Scalar q_oil = ipr(well, 1, well.bhp_limit);
            if (well.oil_rate_limit > Scalar{0}) {
                q_oil = std::min(q_oil, well.oil_rate_limit);
            }
            if (well.pinned) {
                q_oil = well.oil_rate_limit;
            } else if (complementarity_ && well.q_start > Scalar{0}) {
                q_oil = well.q_start;
            } else if (hasTubing(well)) {
                const Scalar p = node_pressure[well.node];
                const Scalar found = thpPotential(well, p);
                if (found > Scalar{0} && found < q_oil) {
                    q_oil = found;
                }
            }
            q_oil = std::max(q_oil, Scalar{0});
            const Scalar bhp = (well.ipr_b[1] < Scalar{0})
                ? std::max((q_oil - well.ipr_a[1]) / well.ipr_b[1], well.bhp_limit)
                : std::max(well.bhp_limit, node_pressure[well.node]);
            x[bhpIdx(w)] = bhp;
            for (int ph = 0; ph < NP; ++ph) {
                x[qwIdx(w, ph)] = std::max(ipr(well, ph, bhp), Scalar{0});
            }
        }
        for (int n = numNodes(); n >= 1; --n) {
            for (int ph = 0; ph < NP; ++ph) {
                Scalar q = node_source_[n][ph];
                for (const int w : wells_at_[n]) {
                    q += wells_[w].efficiency
                       * (x[qwIdx(w, ph)] + (ph == 2 ? wells_[w].lift_gas : Scalar{0}));
                }
                for (const int c : children_[n]) {
                    q += nodes_[c].efficiency * x[qIdx(c, ph)];
                }
                x[qIdx(n, ph)] = q;
            }
        }
        x[lambdaIdx()] = lambda0();
        // A choke starts closed if the guess already carries more oil than the
        // target; updateControls() settles it from the wells' allowances.
        for (int n = 1; n <= numNodes(); ++n) {
            if (isChoke(n)) {
                node_choked_[n] = (x[qIdx(n, 1)] > node_choke_target_[n]) ? 1 : 0;
            }
        }
        return x;
    }

    State pressures(const State& x) const override
    {
        State p(nodes_.size(), terminal_pressure_);
        for (int n = 1; n <= numNodes(); ++n) {
            p[n] = x[pIdx(n)];
        }
        return p;
    }

    /// Every phase of every well, and every well's bhp, at a state.
    std::vector<std::array<Scalar, NP>> wellPhaseRates(const State& x) const override
    {
        std::vector<std::array<Scalar, NP>> q(wells_.size());
        for (int w = 0; w < numWells(); ++w) {
            for (int ph = 0; ph < NP; ++ph) {
                q[w][ph] = x[qwIdx(w, ph)];
            }
        }
        return q;
    }
    State wellBhps(const State& x) const override
    {
        State b(wells_.size());
        for (int w = 0; w < numWells(); ++w) {
            b[w] = x[bhpIdx(w)];
        }
        return b;
    }
    int wellIndex(const std::string& name) const
    {
        for (int w = 0; w < numWells(); ++w) {
            if (wells_[w].name == name) {
                return w;
            }
        }
        return -1;
    }
    /// Give a well a different lift-gas rate -- a trial the optimiser asks
    /// about. Its tubing sees the new alq; its node sees the gas if it adds it.
    void setWellAlq(const int w, const Scalar alq)
    {
        potential_grid_[w].clear();
        // A different lift gas is a different tubing curve: whatever tie the
        // well was on is gone with it, and a trial that inherits Control::Tied
        // is asked for a rate on a curve that no longer passes it.
        if (!controls_.empty() && controls_[w] == Control::Tied) {
            controls_[w] = Control::Thp;
        }
        if (!recent_.empty()) {
            recent_[w] = {Control::Thp, Control::Thp};
        }
        wells_[w].alq = alq;
        if (wells_[w].node_adds_lift_gas) {
            wells_[w].lift_gas = alq;
        }
    }

    /// What a well could flow at node pressure p, ignoring its own rate limit
    /// and any group share: the crossing of its inflow with its tubing, as the
    /// whole phase triple. The gas lift optimiser asks for this potential and
    /// applies the limits itself. Empty when the tubing cannot lift there.
    std::optional<std::array<Scalar, NP + 1>> potentialAt(const int w, const Scalar p) const
    {
        const auto& well = wells_[w];
        if (!hasTubing(well) || !(well.ipr_b[1] < Scalar{0})) {
            return {};
        }
        auto free_well = well;
        free_well.oil_rate_limit = Scalar{0};
        free_well.pinned = false;
        free_well.dead_above = Scalar{0};
        const Scalar q_oil = thpPotential(free_well, p);
        if (!(q_oil > Scalar{0}) || q_oil == std::numeric_limits<Scalar>::max()) {
            return {};
        }
        const Scalar bhp = (q_oil - well.ipr_a[1]) / well.ipr_b[1];
        std::array<Scalar, NP + 1> out{};
        for (int ph = 0; ph < NP; ++ph) {
            out[ph] = std::max(ipr(well, ph, bhp), Scalar{0});
        }
        out[NP] = bhp;
        return out;
    }

    /// Oil rate per well, which is what a caller usually wants back.
    State wellRates(const State& x) const override
    {
        State q(wells_.size());
        for (int w = 0; w < numWells(); ++w) {
            q[w] = x[qwIdx(w, 1)];
        }
        return q;
    }

    Scalar columnScale(const int i) const override
    {
        const bool is_pressure = (i < numNodes())
                              || (i >= bhpIdx(0) && i < lambdaIdx());
        return is_pressure ? pressure_scale_ : rate_scale_;
    }

    /// Bounds by projection, as the injection system has: no node pressure
    /// below an atmosphere, no well bhp below its limit, and no single step
    /// of more than 50 bar on any pressure. The active-set rows are linear
    /// enough in the pressures to do without; the complementarity row is not,
    /// and its first full step from a poor start ran a choke node to minus
    /// three thousand bar.
    State limitStep(const State& x, const State& dx) const override
    {
        Scalar alpha = Scalar{1};
        const Scalar floor = unit::atm;
        const Scalar cap = Scalar{50} * unit::barsa;
        for (int n = 1; n <= numNodes(); ++n) {
            const Scalar d = dx[pIdx(n)];
            if (std::abs(d) > cap) { alpha = std::min(alpha, cap / std::abs(d)); }
            if (d < Scalar{0} && x[pIdx(n)] + alpha * d < floor) {
                alpha = std::min(alpha, (x[pIdx(n)] - floor) / (-d));
            }
        }
        for (int w = 0; w < numWells(); ++w) {
            const Scalar d = dx[bhpIdx(w)];
            if (std::abs(d) > cap) { alpha = std::min(alpha, cap / std::abs(d)); }
            const Scalar lower = std::max(wells_[w].bhp_limit, floor);
            if (d < Scalar{0} && x[bhpIdx(w)] + alpha * d < lower && x[bhpIdx(w)] > lower) {
                alpha = std::min(alpha, (x[bhpIdx(w)] - lower) / (-d));
            }
        }
        alpha = std::max(alpha, Scalar{1e-3});
        if (complementarity_ && analytic_jacobian_) {
            // Clamp the pressures one by one and leave the rates their full
            // step. One alpha for everything let a choke row's pressure
            // demand cut the step to its floor while the shut rows still had
            // their whole rate to remove -- 0.1 % of it per iteration.
            State out = dx;
            for (int n = 1; n <= numNodes(); ++n) {
                Scalar& d = out[pIdx(n)];
                d = std::clamp(d, -cap, cap);
                d = std::max(d, floor - x[pIdx(n)]);
            }
            for (int w = 0; w < numWells(); ++w) {
                Scalar& d = out[bhpIdx(w)];
                d = std::clamp(d, -cap, cap);
                const Scalar lower = std::max(wells_[w].bhp_limit, floor);
                if (x[bhpIdx(w)] > lower) { d = std::max(d, lower - x[bhpIdx(w)]); }
            }
            alpha = Scalar{1};
            return limitCmplRates(x, out, alpha);
        }
        return limitCmplRates(x, dx, alpha);
    }

    State limitCmplRates(const State& x, const State& dx, const Scalar alpha) const
    {
        State d = dx;
        // Two ways a complementarity row proposes a useless oil step: below the
        // liquid-loading hump it linearises to the wrong sign, and a near-tie
        // between two slacks scales badly. So cap the step at what the well
        // could actually deliver, and put a well pushed to zero back on its
        // crossing rather than letting it leave the flowing branch.
        for (int w = 0; w < numWells(); ++w) {
            if (controls_[w] != Control::Cmpl) { continue; }
            const auto& well = wells_[w];
            const int i = qwIdx(w, 1);
            const Scalar p = well.node == 0 ? terminal_pressure_ : x[pIdx(well.node)];
            const Scalar scan = cachedThpPotential(well, p);
            if (cmplDead(w, p)) { continue; }   // the dead row handles it
            const Scalar allow = cmplAllowance(well, scan);
            const Scalar cap = std::max(std::abs(x[i]), allow);
            d[i] = std::clamp(d[i], -cap, cap);
            if (x[i] + alpha * d[i] <= Scalar{0} && well.ipr_b[1] < Scalar{0}) {
                // Re-seed the whole well at the allowance, on its IPR: a rate
                // alone leaves the bhp where a step toward a negative rate put it.
                const Scalar bhp = (allow - well.ipr_a[1]) / well.ipr_b[1];
                for (int ph = 0; ph < NP; ++ph) {
                    d[qwIdx(w, ph)] = (std::max(ipr(well, ph, bhp), Scalar{0}) - x[qwIdx(w, ph)]) / alpha;
                }
                d[bhpIdx(w)] = (bhp - x[bhpIdx(w)]) / alpha;
            }
        }
        static const bool trace = std::getenv("OPM_NETWORK_CM_TRACE") != nullptr;
        if (trace) {
            std::size_t worst = 0;
            for (std::size_t i = 0; i < dx.size(); ++i) { if (std::abs(dx[i]) > std::abs(dx[worst])) { worst = i; } }
            std::fprintf(stderr, "alpha %.3g (max |dx| %.3g at %zu, np %d) |", double(alpha), double(dx[worst]), worst, numNodes());
            for (int n = 1; n <= numNodes(); ++n) {
                std::fprintf(stderr, " %s p %.2f q %.0f tgt %.0f", nodes_[n].name.c_str(), double(x[pIdx(n)] / unit::barsa),
                             double(x[qIdx(n, 1)] * 86400), isChoke(n) ? double(node_choke_target_[n] * 86400) : -1.0);
            }
            for (int w = 0; w < numWells(); ++w) {
                const auto& well = wells_[w];
                const Scalar pw = well.node == 0 ? terminal_pressure_ : x[pIdx(well.node)];
                const Scalar sc = controls_[w] == Control::Cmpl ? cachedThpPotential(well, pw) : Scalar{-1};
                std::fprintf(stderr, " | %s[%c%s] q %.1f d %.1f bhp %.2f scan %.0f", well.name.c_str(), controlLetter(w),
                             cmpl_dead_[w] ? "/dead" : "",
                             double(x[qwIdx(w, 1)] * 86400), double(dx[qwIdx(w, 1)] * 86400), double(x[bhpIdx(w)] / unit::barsa),
                             std::min(double(sc * 86400), 9e9));
            }
            std::fprintf(stderr, "\n");
        }
        State out(d.size());
        for (std::size_t i = 0; i < d.size(); ++i) { out[i] = alpha * d[i]; }
        return out;
    }

    void setAnalyticJacobian(const bool on) { analytic_jacobian_ = on; }
    bool usesAnalyticJacobian() const override { return analytic_jacobian_; }

    /// Close every well that has its own limits to choose between with one
    /// complementarity row instead of an active set: of the rate slack
    /// limit - q, the tubing slack bhp - tubing(p, q) and the bhp slack
    /// bhp - bhp_limit, all are non-negative and at least one is zero.
    /// psi = fb(fb(a, b), c) with fb(u, v) = u + v - sqrt(u^2 + v^2) has exactly
    /// that zero set and is smooth away from the origin, so there is nothing to
    /// switch and nothing to flip -- two wells on ties at once included. Needs
    /// the assembled Jacobian. Pinned, group-held and dead wells keep their rows.
    void setComplementarity(const bool on) { complementarity_ = on; }
    bool usesComplementarity() const { return complementarity_; }

    static Scalar fb(const Scalar u, const Scalar v) { return u + v - std::sqrt(u * u + v * v); }
    /// d fb / du and d fb / dv; the generalised derivative at the origin.
    static std::array<Scalar, 2> dfb(const Scalar u, const Scalar v)
    {
        const Scalar n = std::sqrt(u * u + v * v);
        if (!(n > Scalar{0})) {
            const Scalar g = Scalar{1} - Scalar{1} / std::sqrt(Scalar{2});
            return {g, g};
        }
        return {Scalar{1} - u / n, Scalar{1} - v / n};
    }

    /// A table lookup with its derivatives: the three rate derivatives by
    /// automatic differentiation of the same lookup, the thp derivative by one
    /// more lookup -- the table is piecewise linear in thp, so a small forward
    /// difference inside an interval is the derivative.
    struct Lookup { Scalar value; std::array<Scalar, NP> dq; Scalar dthp; };

    Lookup tableLookup(const int table, const Scalar thp,
                       const std::array<Scalar, NP>& q, const Scalar alq) const
    {
        using Eval = DenseAd::Evaluation<Scalar, NP>;
        // production rates are negative to the table, as in tableBhp()
        const Eval aqua    = Eval::createVariable(-q[0], 0);
        const Eval liquid  = Eval::createVariable(-q[1], 1);
        const Eval vapour  = Eval::createVariable(-q[2], 2);
        const Eval bhp = props_->bhp(table, aqua, liquid, vapour, thp, alq,
                                     Scalar{0}, Scalar{0}, /*use_expvfp=*/false);
        Lookup out;
        out.value = bhp.value();
        for (int ph = 0; ph < NP; ++ph) {
            out.dq[ph] = -bhp.derivative(ph);     // d/dq = -d/d(-q)
        }
        const Scalar h = Scalar{0.01} * unit::barsa;
        out.dthp = (tableBhp(table, thp + h, q, alq) - out.value) / h;
        return out;
    }

    DenseMatrix<Scalar> jacobian(const State& x) const override
    {
        const int nodes = numNodes();
        const int wells = numWells();
        DenseMatrix<Scalar> J(size());
        auto pressure = [&](const int n) { return n == 0 ? terminal_pressure_ : x[pIdx(n)]; };
        auto branchRates = [&](const int n) {
            std::array<Scalar, NP> q{};
            for (int ph = 0; ph < NP; ++ph) { q[ph] = x[qIdx(n, ph)]; }
            return q;
        };
        auto add = [&](const int row, const int col, const Scalar value, const Scalar scale) {
            J(row, col) += value / scale;
        };

        for (int n = 1; n <= nodes; ++n) {
            const auto& node = nodes_[n];
            const int row = n - 1;
            if (isChoke(n) && complementarity_ && analytic_jacobian_ && chokeCanAct(n)) {
                const Scalar upstream = pressure(node.parent);
                Scalar open_p = upstream;
                std::optional<Lookup> branch;
                if (hasTable(node)) {
                    branch = tableLookup(node.vfp_table, upstream, branchRates(n), branch_alq_[n]);
                    open_p = branch->value;
                }
                const Scalar drop = (x[pIdx(n)] - open_p) / pressure_scale_;
                const Scalar surplus = (node_choke_target_[n] - x[qIdx(n, 1)]) / rate_scale_;
                const auto g = dfb(drop, surplus);
                add(row, pIdx(n), g[0], pressure_scale_);
                if (node.parent != 0) {
                    add(row, pIdx(node.parent), -g[0] * (branch ? branch->dthp : Scalar{1}), pressure_scale_);
                }
                if (branch) {
                    for (int ph = 0; ph < NP; ++ph) { add(row, qIdx(n, ph), -g[0] * branch->dq[ph], pressure_scale_); }
                }
                add(row, qIdx(n, 1), -g[1], rate_scale_);
            } else if (isChoke(n) && choked(n) && !(complementarity_ && analytic_jacobian_)) {
                add(row, qIdx(n, 1), 1.0, rate_scale_);
            } else {
                add(row, pIdx(n), 1.0, pressure_scale_);
                if (hasTable(node)) {
                    const auto e = tableLookup(node.vfp_table, pressure(node.parent),
                                               branchRates(n), branch_alq_[n]);
                    if (node.parent != 0) {
                        add(row, pIdx(node.parent), -e.dthp, pressure_scale_);
                    }
                    for (int ph = 0; ph < NP; ++ph) {
                        add(row, qIdx(n, ph), -e.dq[ph], pressure_scale_);
                    }
                } else if (node.parent != 0) {
                    add(row, pIdx(node.parent), -1.0, pressure_scale_);
                }
            }
            for (int ph = 0; ph < NP; ++ph) {
                const int balance = nodes + NP * (n - 1) + ph;
                add(balance, qIdx(n, ph), 1.0, rate_scale_);
                for (const int c : children_[n]) {
                    add(balance, qIdx(c, ph), -nodes_[c].efficiency, rate_scale_);
                }
                for (const int w : wells_at_[n]) {
                    add(balance, qwIdx(w, ph), -wells_[w].efficiency, rate_scale_);
                }
            }
        }

        const bool any_grup = std::find(controls_.begin(), controls_.end(), Control::Grup)
                            != controls_.end();
        for (int w = 0; w < wells; ++w) {
            const auto& well = wells_[w];
            for (int ph = 0; ph < NP; ++ph) {
                const int ipr_row = 4 * nodes + NP * w + ph;
                add(ipr_row, qwIdx(w, ph), 1.0, rate_scale_);
                add(ipr_row, bhpIdx(w), -well.ipr_b[ph], rate_scale_);
            }
            const int row = 4 * nodes + NP * wells + w;
            switch (controls_[w]) {
            case Control::Thp: {
                std::array<Scalar, NP> q{};
                for (int ph = 0; ph < NP; ++ph) { q[ph] = x[qwIdx(w, ph)]; }
                const auto e = tableLookup(well.vfp_table, pressure(well.node), q, well.alq);
                add(row, bhpIdx(w), 1.0, pressure_scale_);
                if (well.node != 0) {
                    add(row, pIdx(well.node), -e.dthp, pressure_scale_);
                }
                for (int ph = 0; ph < NP; ++ph) {
                    add(row, qwIdx(w, ph), -e.dq[ph], pressure_scale_);
                }
                break;
            }
            case Control::Bhp:
                add(row, bhpIdx(w), 1.0, pressure_scale_);
                break;
            case Control::OilRate:
                add(row, qwIdx(w, 1), 1.0, rate_scale_);
                break;
            case Control::Grup:
                add(row, qwIdx(w, 1), 1.0, rate_scale_);
                add(row, lambdaIdx(), -well.guide, rate_scale_);
                break;
            case Control::Cmpl: {
                if (cmplDead(w, pressure(well.node))) {
                    add(row, qwIdx(w, 1), 1.0, rate_scale_);
                    break;
                }
                std::array<Scalar, NP> q{}, qp{};
                for (int ph = 0; ph < NP; ++ph) { q[ph] = x[qwIdx(w, ph)]; qp[ph] = std::max(q[ph], Scalar{0}); }
                const auto e = tableLookup(well.vfp_table, pressure(well.node), qp, well.alq);
                const bool has_rate = well.oil_rate_limit > Scalar{0};
                const Scalar a = has_rate ? (well.oil_rate_limit - q[1]) / rate_scale_ : Scalar{1e6};
                const Scalar b = (x[bhpIdx(w)] - (e.value - well.vfp_dp)) / pressure_scale_;
                const Scalar c = (x[bhpIdx(w)] - well.bhp_limit) / pressure_scale_;
                const auto g_ab = dfb(a, b);                 // d fb(a,b) / da, db
                const auto g_oc = dfb(fb(a, b), c);          // d psi / d fb(a,b), dc
                const Scalar da = g_oc[0] * g_ab[0], db = g_oc[0] * g_ab[1], dc = g_oc[1];
                if (has_rate) { add(row, qwIdx(w, 1), -da, rate_scale_); }
                add(row, bhpIdx(w), db + dc, pressure_scale_);
                if (well.node != 0) { add(row, pIdx(well.node), -db * e.dthp, pressure_scale_); }
                for (int ph = 0; ph < NP; ++ph) {
                    if (q[ph] > Scalar{0}) { add(row, qwIdx(w, ph), -db * e.dq[ph], pressure_scale_); }
                }
                break;
            }
            case Control::Shut:
                if (well.ipr_b[1] < Scalar{0}) { add(row, qwIdx(w, 1), 1.0, rate_scale_); }
                else { add(row, bhpIdx(w), 1.0, pressure_scale_); }
                break;
            case Control::Tied: {
                std::array<Scalar, NP> q{};
                for (int ph = 0; ph < NP; ++ph) { q[ph] = x[qwIdx(w, ph)]; }
                const auto e = tableLookup(well.vfp_table, pressure(well.node), q, well.alq);
                const Scalar a = (well.oil_rate_limit - q[1]) / rate_scale_;
                const Scalar b = (x[bhpIdx(w)] - (e.value - well.vfp_dp)) / pressure_scale_;
                const Scalar norm = std::sqrt(a * a + b * b);
                // at the origin the generalised Jacobian; (1 - 1/sqrt 2) on both
                const Scalar da = (norm > Scalar{0}) ? Scalar{1} - a / norm : Scalar{1} - Scalar{1} / std::sqrt(Scalar{2});
                const Scalar db = (norm > Scalar{0}) ? Scalar{1} - b / norm : Scalar{1} - Scalar{1} / std::sqrt(Scalar{2});
                add(row, qwIdx(w, 1), -da, rate_scale_);
                add(row, bhpIdx(w), db, pressure_scale_);
                if (well.node != 0) {
                    add(row, pIdx(well.node), -db * e.dthp, pressure_scale_);
                }
                for (int ph = 0; ph < NP; ++ph) {
                    add(row, qwIdx(w, ph), -db * e.dq[ph], pressure_scale_);
                }
                break;
            }
            }
            if (grouped() && any_grup && well.in_group) {
                add(lambdaIdx(), qwIdx(w, 1), 1.0, rate_scale_);
            }
        }
        if (!(grouped() && any_grup)) {
            add(lambdaIdx(), lambdaIdx(), 1.0, rate_scale_);
        }
        return J;
    }

private:
    bool analytic_jacobian_ = false;
    bool complementarity_ = false;

    /// The multiplier an even guide-rate split would imply, which is what the
    /// row pins it to while nobody is on group control.
    Scalar lambda0() const
    {
        Scalar guides = 0.0;
        for (const auto& w : wells_) {
            if (w.in_group) {
                guides += w.guide;
            }
        }
        return guides > Scalar{0} ? group_target_ / guides : Scalar{0};
    }

    const VFPProdProperties<Scalar>* props_;
    const UnitSystem* units_;
    std::vector<Node> nodes_;
    std::vector<Scalar> branch_alq_;
    std::vector<std::array<Scalar, NP>> node_source_;
    std::vector<Scalar> node_choke_target_;
    mutable std::vector<char> node_choked_;
    std::vector<Scalar> node_choke_pressure_;
    mutable std::vector<std::map<long, Scalar>> potential_grid_;
    /// The last two controls each well was on. A well that returns to the
    /// control it left two selections ago, across its rate limit, is on a tie
    /// and goes to Control::Tied for the rest of the solve.
    std::vector<std::array<Control, 2>> recent_;
    std::vector<char> cmpl_;            // on the complementarity row this solve
    bool cmpl_decided_ = false;
    /// Died during this solve -- the scan found no crossing at some iterate
    /// -- and stays shut for the rest of it. Reviving with the pressure the
    /// shut-in lowers makes a system with no fixed point.
    mutable std::vector<char> cmpl_dead_;
    mutable std::vector<Scalar> cmpl_alive_p_;   // lowest node pressure seen alive at
    std::vector<Well> wells_;
    std::vector<std::vector<int>> children_;
    std::vector<std::vector<int>> wells_at_;
    std::vector<Control> controls_;

    Scalar terminal_pressure_ = 0.0;
    Scalar group_target_ = 0.0;
    Scalar rate_scale_ = 0.0;
    Scalar pressure_scale_ = unit::barsa;
};

/// Write a production system and its starting pressures, for replay in the
/// bench against the same tables. Everything the solve depends on is here.
template<class Scalar>
void write(const ProductionSystem<Scalar>& system, const std::vector<Scalar>& guess, std::ostream& os)
{
    os.precision(17);
    os << "production\n"
       << "terminal " << system.terminalPressure() << '\n'
       << "group_target " << system.groupTarget() << '\n'
       << "analytic_jacobian " << system.usesAnalyticJacobian() << '\n';
    for (int n = 0; n < static_cast<int>(system.nodes().size()); ++n) {
        const auto& node = system.nodes()[n];
        const auto src = system.nodeSource(n);
        os << "node " << node.name << ' ' << node.parent << ' ' << node.vfp_table << ' '
           << node.efficiency << ' ' << system.branchAlq(n) << ' ' << src[0] << ' ' << src[1] << ' '
           << src[2] << ' ' << system.chokeTarget(n) << '\n';
    }
    for (const auto& w : system.wells()) {
        os << "well " << w.name << ' ' << w.node << ' ' << w.vfp_table << ' ' << w.alq << ' '
           << w.ipr_a[0] << ' ' << w.ipr_a[1] << ' ' << w.ipr_a[2] << ' '
           << w.ipr_b[0] << ' ' << w.ipr_b[1] << ' ' << w.ipr_b[2] << ' '
           << w.bhp_limit << ' ' << w.oil_rate_limit << ' ' << w.in_group << ' ' << w.guide << ' '
           << w.efficiency << ' ' << w.vfp_dp << ' ' << w.pinned << ' ' << w.dead_above << ' '
           << w.lift_gas << ' ' << w.node_adds_lift_gas << ' ' << w.q_start << ' ' << w.shut << '\n';
    }
    os << "guess";
    for (const auto p : guess) { os << ' ' << p; }
    os << '\n';
}

template<class Scalar>
std::pair<ProductionSystem<Scalar>, std::vector<Scalar>>
readProduction(std::istream& is, const VFPProdProperties<Scalar>& props, const UnitSystem& units)
{
    ProductionSystem<Scalar> system(props, units);
    std::vector<Scalar> guess;
    std::string line;
    while (std::getline(is, line)) {
        std::istringstream in(line);
        std::string tag;
        if (!(in >> tag)) { continue; }
        if (tag == "terminal") { Scalar v; in >> v; system.setTerminalPressure(v); }
        else if (tag == "group_target") { Scalar v; in >> v; system.setGroupTarget(v); }
        else if (tag == "analytic_jacobian") { int v; in >> v; system.setAnalyticJacobian(v != 0); }
        else if (tag == "node") {
            Node n; Scalar alq, s0, s1, s2, choke;
            in >> n.name >> n.parent >> n.vfp_table >> n.efficiency >> alq >> s0 >> s1 >> s2 >> choke;
            system.addNode(n, alq);
            const int idx = static_cast<int>(system.nodes().size()) - 1;
            system.setNodeSource(idx, {s0, s1, s2});
            if (choke > Scalar{0}) { system.setChokeTarget(idx, choke); }
        } else if (tag == "well") {
            typename ProductionSystem<Scalar>::Well w; int in_group, pinned, adds;
            in >> w.name >> w.node >> w.vfp_table >> w.alq
               >> w.ipr_a[0] >> w.ipr_a[1] >> w.ipr_a[2] >> w.ipr_b[0] >> w.ipr_b[1] >> w.ipr_b[2]
               >> w.bhp_limit >> w.oil_rate_limit >> in_group >> w.guide >> w.efficiency >> w.vfp_dp
               >> pinned >> w.dead_above >> w.lift_gas >> adds;
            w.in_group = in_group != 0; w.pinned = pinned != 0; w.node_adds_lift_gas = adds != 0;
            in >> w.q_start;             // older dumps: stays 0
            int shut = 0; in >> shut; w.shut = shut != 0;
            system.addWell(std::move(w));
        } else if (tag == "guess") { Scalar v; while (in >> v) { guess.push_back(v); } }
    }
    system.finish();
    return {std::move(system), std::move(guess)};
}

} // namespace Opm::NetworkSolve

#endif // OPM_NETWORK_PRODUCTION_SYSTEM_HEADER_INCLUDED
