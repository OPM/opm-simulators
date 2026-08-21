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
#ifndef OPM_NETWORK_SYSTEM_HEADER_INCLUDED
#define OPM_NETWORK_SYSTEM_HEADER_INCLUDED

#include <opm/input/eclipse/EclipseState/Phase.hpp>
#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/simulators/wells/VFPHelpers.hpp>
#include <opm/simulators/wells/VFPInjProperties.hpp>
#include <opm/simulators/wells/VFPProdProperties.hpp>

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

/// An injection network solved simultaneously in its pressures and its rates.
///
/// The unknowns are the pressure of every non-terminal node, the rate through
/// every node's parent branch, each well's (rate, bhp), and a group multiplier
/// when a target is active. The equations are the branch pressure drops, the
/// node mass balances, each well's inflow performance and one control equation
/// per well, plus the group target.
///
/// The alternative is to eliminate the rates and iterate on the node pressures
/// alone, which is what the fixed-point and bracketing methods do. That residual
/// is only piecewise differentiable -- the control limits put kinks in it -- and
/// needs globalising. This one holds its controls fixed while a step is taken,
/// so it is smooth within an active set and a plain Newton suffices.
///
/// Both the simulator and tests/test_networksolve.cpp fill this in; the wells
/// differ (the simulator's inflow performance comes from the well Jacobian, the
/// bench's from a reference operating point) but the system does not.

constexpr int NoTable = 9999;   // GNETINJE's "no table": pressure passes through

struct Node
{
    std::string name;
    int parent = -1;            // -1 only for the terminal
    int vfp_table = NoTable;
};

template<class Scalar>
struct Well
{
    std::string name;
    int node = 0;
    int vfp_table = 0;
    /// Inflow performance, q = ipr_a + ipr_b * bhp. The simulator's implicit
    /// IPR stores it as q = b*bhp - a, so ipr_a is the negated one.
    Scalar ipr_a = 0.0;
    Scalar ipr_b = 0.0;
    Scalar bhp_limit = 0.0;
    Scalar rate_limit = 0.0;
    /// Whether the group allocated this well. It counts against the group's
    /// target whatever control it ends up on -- a well that runs into its own
    /// bhp or rate limit still injects, and the wells that scale with the
    /// multiplier have to make up the remainder, not the whole target.
    bool in_group = false;
    Scalar guide = 0.0;         // share of a group target
    /// Rate to start the solve from. Zero means work one out from the tables,
    /// which is all the bench can do; the simulator knows what the well is
    /// actually doing and should say so, or the first control selection is made
    /// on a rate that has nothing to do with the current state.
    Scalar q_start = 0.0;
};

/// Which equation closes a well.
enum class Control { Thp, Bhp, Rate, Grup };

template<class Scalar>
struct Result
{
    bool converged = false;
    int iterations = 0;
    std::vector<Scalar> node_pressure;   // every node, terminal included
    std::vector<Scalar> well_rate;       // per well, in the order they were added

    /// Why it stopped, for the caller to report. A converged solve leaves these
    /// at the values that satisfied the test.
    Scalar residual = 0.0;          // max norm of the scaled residual
    bool controls_moving = false;   // an active-set change on the last iteration
    bool guides_moving = false;     // group shares still settling
    /// One letter per well per iteration for the last few iterations, so a
    /// cycling active set can be read off: T thp, B bhp, R rate, G group.
    std::string control_trace;
};

/// Dense square system. The networks this solves have tens of unknowns, so
/// Gaussian elimination with partial pivoting is the whole story.
template<class Scalar>
class DenseMatrix
{
public:
    explicit DenseMatrix(const int n) : n_(n), a_(n * n, 0.0) {}

    Scalar& operator()(const int i, const int j) { return a_[i * n_ + j]; }
    Scalar operator()(const int i, const int j) const { return a_[i * n_ + j]; }

    /// Solves A y = b. False if A is singular to working precision.
    bool solve(std::vector<Scalar> b, std::vector<Scalar>& y) const
    {
        auto a = a_;
        y.assign(n_, 0.0);
        for (int k = 0; k < n_; ++k) {
            int pivot = k;
            for (int i = k + 1; i < n_; ++i) {
                if (std::abs(a[i * n_ + k]) > std::abs(a[pivot * n_ + k])) {
                    pivot = i;
                }
            }
            if (std::abs(a[pivot * n_ + k]) < 1e-300) {
                return false;
            }
            if (pivot != k) {
                for (int j = 0; j < n_; ++j) {
                    std::swap(a[k * n_ + j], a[pivot * n_ + j]);
                }
                std::swap(b[k], b[pivot]);
            }
            for (int i = k + 1; i < n_; ++i) {
                const Scalar f = a[i * n_ + k] / a[k * n_ + k];
                for (int j = k; j < n_; ++j) {
                    a[i * n_ + j] -= f * a[k * n_ + j];
                }
                b[i] -= f * b[k];
            }
        }
        for (int i = n_ - 1; i >= 0; --i) {
            Scalar sum = b[i];
            for (int j = i + 1; j < n_; ++j) {
                sum -= a[i * n_ + j] * y[j];
            }
            y[i] = sum / a[i * n_ + i];
        }
        return true;
    }

private:
    int n_;
    std::vector<Scalar> a_;
};

template<class Scalar>
class System
{
public:
    using State = std::vector<Scalar>;
    using ScalarType = Scalar;

    System(const VFPInjProperties<Scalar>& props, const Phase phase)
        : props_(&props), phase_(phase)
    {}

    void addNode(Node n) { nodes_.push_back(std::move(n)); }
    void addWell(Well<Scalar> w) { wells_.push_back(std::move(w)); }
    void setTerminalPressure(const Scalar p) { terminal_pressure_ = p; }
    /// The target of the one group this system can carry.
    ///
    /// A well belongs to that group by being under group control, not by where
    /// it sits in the network: the group tree and the network tree are
    /// independent and share only their leaves, and nothing here assumes
    /// otherwise.
    ///
    /// What it does assume is a **single** constraining group. Two groups
    /// binding different subsets of these wells would be summed into one target,
    /// which is wrong, and nested groups need a multiplier each with a well's
    /// share the product down its chain. Neither is modelled.
    void setGroupTarget(const Scalar target) { group_target_ = target; }

    /// Residual scale for the rate rows. Without one, rate and pressure rows
    /// differ by several decades and no single tolerance means anything. The
    /// default from finish() is a hundredth of the largest rate in play, which
    /// is why it has to be called after the wells are in.
    void setRateScale(const Scalar s) { rate_scale_ = s; }

    /// Resolve the tree and the defaults. Call once everything is added.
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
                w.guide = std::max(w.rate_limit, Scalar{1.0});
            }
        }
        if (rate_scale_ <= 0.0) {
            Scalar largest = group_target_;
            for (const auto& w : wells_) {
                largest = std::max(largest, w.rate_limit);
            }
            rate_scale_ = std::max(largest * Scalar{0.01},
                                   unit::convert::from(1.0, unit::cubic(unit::meter) / unit::day));
        }
        controls_.assign(wells_.size(), Control::Thp);
        for (std::size_t w = 0; w < wells_.size(); ++w) {
            if (grouped() && wells_[w].in_group) {
                controls_[w] = Control::Grup;
            }
        }
    }

    int numNodes() const { return static_cast<int>(nodes_.size()) - 1; }
    int numWells() const { return static_cast<int>(wells_.size()); }
    bool grouped() const { return group_target_ > 0.0; }
    int size() const { return 2 * numNodes() + 2 * numWells() + (grouped() ? 1 : 0); }

    Phase phase() const { return phase_; }
    Scalar terminalPressure() const { return terminal_pressure_; }
    Scalar groupTarget() const { return group_target_; }
    const std::vector<Node>& nodes() const { return nodes_; }
    const std::vector<Well<Scalar>>& wells() const { return wells_; }
    Control control(const int w) const { return controls_[w]; }

    /// One letter for the trace a failed solve reports.
    char controlLetter(const int w) const
    {
        switch (controls_[w]) {
        case Control::Thp:  return 'T';
        case Control::Bhp:  return 'B';
        case Control::Rate: return 'R';
        case Control::Grup: return 'G';
        }
        return '?';
    }

    int pIdx(const int node) const { return node - 1; }
    int qIdx(const int node) const { return numNodes() + node - 1; }
    int qwIdx(const int w) const { return 2 * numNodes() + w; }
    int bhpIdx(const int w) const { return 2 * numNodes() + numWells() + w; }
    int lambdaIdx() const { return 2 * numNodes() + 2 * numWells(); }

    bool hasTable(const Node& n) const { return n.vfp_table != NoTable; }

    static Scalar ipr(const Well<Scalar>& w, const Scalar bhp) { return w.ipr_a + w.ipr_b * bhp; }

    /// Clamp table lookups to the axes, as the fixed-point pressure computation
    /// does. Leave this off for a Newton: outside the box the residual then goes
    /// flat and the Jacobian is singular in the rates, so there is nothing to
    /// descend. limitStep() is the treatment that works. It exists here only so
    /// that the comparison can be made -- see test_networksolve.cpp.
    void setClampToAxes(const bool on) { clamp_to_axes_ = on; }

    /// A table lookup with the two derivatives the Jacobian needs. They come
    /// free with the interpolation and are otherwise thrown away.
    struct Lookup
    {
        Scalar value = 0.0;
        Scalar dthp = 0.0;    // d(bhp)/d(thp)
        Scalar dflo = 0.0;    // d(bhp)/d(rate)
    };

    Lookup tableLookup(const int table, const Scalar thp, const Scalar q_in) const
    {
        Scalar q = q_in;
        Scalar p = thp;
        const auto& t = props_->getTable(table);
        bool clamped_flo = false;
        bool clamped_thp = false;
        if (clamp_to_axes_) {
            const Scalar lo = t.getFloAxis().front(), hi = t.getFloAxis().back();
            const Scalar plo = t.getTHPAxis().front(), phi = t.getTHPAxis().back();
            clamped_flo = (q < lo) || (q > hi);
            clamped_thp = (p < plo) || (p > phi);
            q = std::clamp(q, lo, hi);
            p = std::clamp(p, plo, phi);
        }
        const Scalar aqua = (phase_ == Phase::WATER) ? q : Scalar{0};
        const Scalar vapour = (phase_ == Phase::GAS) ? q : Scalar{0};
        const auto e = VFPHelpers<Scalar>::bhp(t, aqua, liquid_, vapour, p);
        // Where the lookup was clamped the value no longer moves with the input,
        // which is exactly the flat residual that makes clamping a bad idea for
        // a Newton -- but the derivative has to report it honestly.
        return {e.value, clamped_thp ? Scalar{0} : e.dthp, clamped_flo ? Scalar{0} : e.dflo};
    }

    /// Downstream pressure of a branch, or a well's bhp: the same table lookup.
    Scalar tableBhp(const int table, const Scalar thp, const Scalar q_in) const
    {
        Scalar q = q_in;
        Scalar p = thp;
        if (clamp_to_axes_) {
            const auto& t = props_->getTable(table);
            q = std::clamp(q, t.getFloAxis().front(), t.getFloAxis().back());
            p = std::clamp(p, t.getTHPAxis().front(), t.getTHPAxis().back());
        }
        const Scalar aqua = (phase_ == Phase::WATER) ? q : Scalar{0};
        const Scalar vapour = (phase_ == Phase::GAS) ? q : Scalar{0};
        return props_->bhp(table, aqua, Scalar{0}, vapour, p);
    }

    /// Rate this well would take on THP control at a given node pressure: its
    /// inflow performance met with its tubing curve, then its own limits.
    ///
    /// This is the well's capability at a network pressure, which is what a
    /// share of a group target should be proportional to. Its current rate is
    /// not: that is the split one is trying to decide, so using it as the guide
    /// makes the allocation reproduce whatever it already was.
    Scalar thpPotential(const Well<Scalar>& w, const Scalar p_node) const
    {
        const auto& t = props_->getTable(w.vfp_table);
        const Scalar lo = t.getFloAxis().front();
        const Scalar hi = std::min(w.rate_limit, t.getFloAxis().back());
        if (!(hi > lo)) {
            return Scalar{0};
        }
        // bhp falls with rate at fixed thp in these tables, so f is decreasing.
        const auto f = [&](const Scalar q) { return ipr(w, tableBhp(w.vfp_table, p_node, q)) - q; };
        if (f(lo) <= Scalar{0}) {
            return Scalar{0};
        }
        if (f(hi) >= Scalar{0}) {
            return hi;
        }
        Scalar a = lo, b = hi, q = hi;
        for (int it = 0; it < 60; ++it) {
            q = Scalar{0.5} * (a + b);
            (f(q) > Scalar{0} ? a : b) = q;
        }
        return std::clamp(q, Scalar{0}, w.rate_limit);
    }

    /// Take the guide rates from thpPotential() at the current node pressures
    /// instead of whatever the caller supplied. Only meaningful with a group
    /// target, and only when the caller has no better guide of its own.
    void setGuidesFromPotential(const bool on) { guides_from_potential_ = on; }

    /// Recompute the guides from the current iterate. Returns the largest
    /// relative change, so the caller can tell when they have settled.
    Scalar refreshGuides(const State& x)
    {
        if (!guides_from_potential_ || !grouped()) {
            return Scalar{0};
        }
        Scalar moved = 0.0;
        for (auto& w : wells_) {
            const Scalar p = (w.node == 0) ? terminal_pressure_ : x[pIdx(w.node)];
            const Scalar potential = thpPotential(w, p);
            if (potential > Scalar{0}) {
                moved = std::max(moved, std::abs(potential - w.guide) / std::max(w.guide, potential));
                w.guide = potential;
            }
        }
        return moved;
    }

    /// Largest rate the table describes. Past it the cells are zero-filled and
    /// the interpolation runs away, so this is the edge of the feasible set.
    Scalar maxFlow(const int table) const { return props_->getTable(table).getFloAxis().back(); }

    State residual(const State& x) const
    {
        const int nodes = numNodes();
        const int wells = numWells();
        State r(size(), 0.0);

        auto pressure = [&](const int n) { return n == 0 ? terminal_pressure_ : x[pIdx(n)]; };

        for (int n = 1; n <= nodes; ++n) {
            const auto& node = nodes_[n];
            const Scalar upstream = pressure(node.parent);
            r[n - 1] = hasTable(node)
                ? x[pIdx(n)] - tableBhp(node.vfp_table, upstream, x[qIdx(n)])
                : x[pIdx(n)] - upstream;

            Scalar balance = x[qIdx(n)];
            for (const int c : children_[n]) {
                balance -= x[qIdx(c)];
            }
            for (const int w : wells_at_[n]) {
                balance -= x[qwIdx(w)];
            }
            r[nodes + n - 1] = balance;
        }

        Scalar injected = 0.0;
        for (int w = 0; w < wells; ++w) {
            const auto& well = wells_[w];
            const Scalar q = x[qwIdx(w)];
            const Scalar bhp = x[bhpIdx(w)];
            // Every well the group allocated counts against the target, on
            // whatever control it ended up. Counting only those still on group
            // control asks the rest to deliver the whole target while a limited
            // well injects on top of it, and the group over-delivers by exactly
            // that well's rate. Counting wells the group never allocated is the
            // opposite error and cannot be satisfied at all.
            if (well.in_group) {
                injected += q;
            }

            r[2 * nodes + w] = (q - ipr(well, bhp)) / rate_scale_;

            Scalar& control = r[2 * nodes + wells + w];
            switch (controls_[w]) {
            case Control::Thp:
                control = (bhp - tableBhp(well.vfp_table, pressure(well.node), q)) / pressure_scale_;
                break;
            case Control::Bhp:
                control = (bhp - well.bhp_limit) / pressure_scale_;
                break;
            case Control::Rate:
                control = (q - well.rate_limit) / rate_scale_;
                break;
            case Control::Grup:
                control = (q - well.guide * x[lambdaIdx()]) / rate_scale_;
                break;
            }
        }

        if (grouped()) {
            // With nobody on group control the multiplier is free, so pin it
            // rather than hand the Newton a singular column.
            const bool any = std::find(controls_.begin(), controls_.end(), Control::Grup)
                           != controls_.end();
            r[lambdaIdx()] = any ? (injected - group_target_) / rate_scale_
                                 : (x[lambdaIdx()] - lambda0()) / rate_scale_;
        }

        for (int n = 0; n < nodes; ++n) {
            r[n] /= pressure_scale_;
            r[nodes + n] /= rate_scale_;
        }
        return r;
    }

    /// Take the last well out of the group, so a test can build a network whose
    /// group does not hold every well on it.
    void dropLastFromGroup()
    {
        if (!wells_.empty()) {
            wells_.back().in_group = false;
        }
    }

    /// Reselect each well's control: the most restrictive violated limit wins,
    /// the same rule a clamp would apply. Choosing by a fixed priority instead
    /// makes the active set chatter and the Newton never terminates. Returns
    /// true if anything moved -- converging with a control still moving is not
    /// converged.
    bool updateControls(const State& x)
    {
        bool changed = false;
        for (int w = 0; w < numWells(); ++w) {
            const auto& well = wells_[w];
            const Scalar q = x[qwIdx(w)];

            // Given the node pressure, every control determines the well
            // completely and so names the rate it would allow. The binding one
            // is simply the smallest: nothing is "violated", and a control stops
            // binding by being overtaken, which is how a well leaves one.
            //
            // Nothing here reads the iterate's q or bhp. Mid-Newton those are
            // not a consistent well state, and on rate control q *is* the limit,
            // so a test against them has no stable answer.
            const Scalar p_node = (well.node == 0) ? terminal_pressure_ : x[pIdx(well.node)];

            auto wanted = Control::Thp;
            Scalar smallest = thpPotential(well, p_node);
            auto consider = [&](const Control c, const Scalar allows) {
                if (allows < smallest) {
                    smallest = allows;
                    wanted = c;
                }
            };
            consider(Control::Bhp, ipr(well, well.bhp_limit));
            consider(Control::Rate, well.rate_limit);
            if (grouped() && well.in_group) {
                consider(Control::Grup, well.guide * x[lambdaIdx()]);
            }

            changed |= (wanted != controls_[w]);
            controls_[w] = wanted;
        }
        return changed;
    }

    /// A starting point derived from a guess at every node's pressure.
    State start(const State& node_pressure) const
    {
        State x(size(), 0.0);
        std::vector<Scalar> well_rate(numWells());
        for (int n = 1; n <= numNodes(); ++n) {
            x[pIdx(n)] = node_pressure[n];
        }
        for (int w = 0; w < numWells(); ++w) {
            const auto& well = wells_[w];
            const Scalar guess = well.q_start > Scalar{0}
                ? well.q_start : std::max(well.rate_limit * Scalar{0.1}, rate_scale_);
            x[bhpIdx(w)] = tableBhp(well.vfp_table, node_pressure[well.node], guess);
            // Deliberately just inside the limit, never exactly on it: the
            // control tests below are inclusive, so opening on the limit would
            // put every well on rate control before the solve has begun.
            const Scalar most = Scalar{0.999} * well.rate_limit;
            x[qwIdx(w)] = well.q_start > Scalar{0}
                ? std::min(well.q_start, most)
                : std::clamp(ipr(well, x[bhpIdx(w)]), Scalar{0}, most);
            well_rate[w] = x[qwIdx(w)];
        }
        for (int n = numNodes(); n >= 1; --n) {
            Scalar q = 0.0;
            for (const int w : wells_at_[n]) {
                q += well_rate[w];
            }
            for (const int c : children_[n]) {
                q += x[qIdx(c)];
            }
            x[qIdx(n)] = q;
        }
        if (grouped()) {
            x[lambdaIdx()] = lambda0();
        }
        return x;
    }

    /// Rate of every well, in the order they were added.
    State wellRates(const State& x) const
    {
        State q(wells_.size());
        for (int w = 0; w < numWells(); ++w) {
            q[w] = x[qwIdx(w)];
        }
        return q;
    }

    /// Pressure at every node, terminal included.
    State pressures(const State& x) const
    {
        State p(nodes_.size(), terminal_pressure_);
        for (int n = 1; n <= numNodes(); ++n) {
            p[n] = x[pIdx(n)];
        }
        return p;
    }

    /// Natural magnitude of unknown i, so one step cap or trust radius can apply
    /// to a vector holding both pressures and rates.
    Scalar columnScale(const int i) const
    {
        const bool is_pressure = (i < numNodes()) || (i >= bhpIdx(0) && i < lambdaIdx());
        return is_pressure ? pressure_scale_ : rate_scale_;
    }

    /// Keep the branch flows inside the box the tables describe, by projecting
    /// the offending components rather than scaling all of the step -- one
    /// binding rate should not throttle the pressure updates too. Only the
    /// branch flows: a well's own rate limit already has a control equation, and
    /// bounding it would stop that control ever activating.
    State limitStep(const State& x, const State& dx) const
    {
        State limited = dx;
        for (int n = 1; n <= numNodes(); ++n) {
            const auto& node = nodes_[n];
            if (!hasTable(node)) {
                continue;
            }
            const int i = qIdx(n);
            const Scalar hi = maxFlow(node.vfp_table);
            if (x[i] <= hi && x[i] + limited[i] > hi) {
                limited[i] = hi - x[i];
            }
            if (x[i] >= 0.0 && x[i] + limited[i] < 0.0) {
                limited[i] = -x[i];
            }
        }
        return limited;
    }

    Scalar pressureScale() const { return pressure_scale_; }

    /// Assemble the Jacobian from the table derivatives instead of differencing
    /// the residual. Everything but the two branch/tubing lookups is constant,
    /// so this is n+1 residual evaluations replaced by one pass.
    void setAnalyticJacobian(const bool on) { analytic_jacobian_ = on; }
    bool usesAnalyticJacobian() const { return analytic_jacobian_; }

    /// The Jacobian of residual() at x, entry by entry.
    DenseMatrix<Scalar> jacobian(const State& x) const
    {
        const int nodes = numNodes();
        const int wells = numWells();
        DenseMatrix<Scalar> J(size());

        auto pressure = [&](const int n) { return n == 0 ? terminal_pressure_ : x[pIdx(n)]; };
        // Row i is divided by scale(i) in residual(), so its derivatives are too.
        auto add = [&](const int row, const int col, const Scalar value, const Scalar scale) {
            J(row, col) += value / scale;
        };

        for (int n = 1; n <= nodes; ++n) {
            const auto& node = nodes_[n];
            const int row = n - 1;
            add(row, pIdx(n), 1.0, pressure_scale_);
            if (hasTable(node)) {
                const auto e = tableLookup(node.vfp_table, pressure(node.parent), x[qIdx(n)]);
                if (node.parent != 0) {
                    add(row, pIdx(node.parent), -e.dthp, pressure_scale_);
                }
                add(row, qIdx(n), -e.dflo, pressure_scale_);
            } else if (node.parent != 0) {
                add(row, pIdx(node.parent), -1.0, pressure_scale_);
            }

            const int balance = nodes + n - 1;
            add(balance, qIdx(n), 1.0, rate_scale_);
            for (const int c : children_[n]) {
                add(balance, qIdx(c), -1.0, rate_scale_);
            }
            for (const int w : wells_at_[n]) {
                add(balance, qwIdx(w), -1.0, rate_scale_);
            }
        }

        for (int w = 0; w < wells; ++w) {
            const auto& well = wells_[w];
            const int ipr_row = 2 * nodes + w;
            add(ipr_row, qwIdx(w), 1.0, rate_scale_);
            add(ipr_row, bhpIdx(w), -well.ipr_b, rate_scale_);

            const int row = 2 * nodes + wells + w;
            switch (controls_[w]) {
            case Control::Thp: {
                const auto e = tableLookup(well.vfp_table, pressure(well.node), x[qwIdx(w)]);
                add(row, bhpIdx(w), 1.0, pressure_scale_);
                if (well.node != 0) {
                    add(row, pIdx(well.node), -e.dthp, pressure_scale_);
                }
                add(row, qwIdx(w), -e.dflo, pressure_scale_);
                break;
            }
            case Control::Bhp:
                add(row, bhpIdx(w), 1.0, pressure_scale_);
                break;
            case Control::Rate:
                add(row, qwIdx(w), 1.0, rate_scale_);
                break;
            case Control::Grup:
                add(row, qwIdx(w), 1.0, rate_scale_);
                add(row, lambdaIdx(), -well.guide, rate_scale_);
                break;
            }
        }

        if (grouped()) {
            const bool any = std::find(controls_.begin(), controls_.end(), Control::Grup)
                           != controls_.end();
            if (any) {
                // The group's own wells, matching the residual. Differentiating
                // every well instead is a wrong derivative that no bench case
                // could see, because there every well is in the group.
                for (int w = 0; w < wells; ++w) {
                    if (wells_[w].in_group) {
                        add(lambdaIdx(), qwIdx(w), 1.0, rate_scale_);
                    }
                }
            } else {
                add(lambdaIdx(), lambdaIdx(), 1.0, rate_scale_);
            }
        }
        return J;
    }

private:
    Scalar lambda0() const
    {
        Scalar guides = 0.0;
        for (const auto& w : wells_) {
            guides += w.guide;
        }
        return guides > 0.0 ? group_target_ / guides : Scalar{0};
    }

    const VFPInjProperties<Scalar>* props_;
    Phase phase_;
    std::vector<Node> nodes_;
    std::vector<Well<Scalar>> wells_;
    std::vector<std::vector<int>> children_;
    std::vector<std::vector<int>> wells_at_;
    std::vector<Control> controls_;

    Scalar terminal_pressure_ = 0.0;
    Scalar group_target_ = 0.0;
    Scalar rate_scale_ = 0.0;
    Scalar liquid_ = 0.0;
    bool clamp_to_axes_ = false;
    bool analytic_jacobian_ = false;
    bool guides_from_potential_ = false;
    Scalar pressure_scale_ = unit::barsa;
};

/// Write everything the solve works from, so a failure can be replayed offline.
/// The VFP tables are not included -- the reader supplies them from the deck.
template<class Scalar>
void write(const System<Scalar>& system, const std::vector<Scalar>& guess, std::ostream& os)
{
    os << "phase " << (system.phase() == Phase::GAS ? "GAS" : "WATER") << '\n'
       << "terminal " << system.terminalPressure() << '\n'
       << "group_target " << system.groupTarget() << '\n';
    for (const auto& n : system.nodes()) {
        os << "node " << n.name << ' ' << n.parent << ' ' << n.vfp_table << '\n';
    }
    for (const auto& w : system.wells()) {
        os << "well " << w.name << ' ' << w.node << ' ' << w.vfp_table << ' '
           << w.ipr_a << ' ' << w.ipr_b << ' ' << w.bhp_limit << ' '
           << w.rate_limit << ' ' << w.guide << ' ' << w.q_start << '\n';
    }
    os << "guess";
    for (const auto p : guess) {
        os << ' ' << p;
    }
    os << '\n';
}

/// Rebuild a written system against tables the caller already has. Returns the
/// system and the starting pressures it was given.
template<class Scalar>
std::pair<System<Scalar>, std::vector<Scalar>>
read(std::istream& is, const VFPInjProperties<Scalar>& props)
{
    std::string tag;
    Phase phase = Phase::GAS;
    Scalar terminal = 0.0, target = 0.0;
    std::vector<Node> nodes;
    std::vector<Well<Scalar>> wells;
    std::vector<Scalar> guess;

    std::string line;
    while (std::getline(is, line)) {
        std::istringstream in(line);
        if (!(in >> tag)) {
            continue;
        }
        if (tag == "phase") {
            std::string name;
            in >> name;
            phase = (name == "GAS") ? Phase::GAS : Phase::WATER;
        } else if (tag == "terminal") {
            in >> terminal;
        } else if (tag == "group_target") {
            in >> target;
        } else if (tag == "node") {
            Node n;
            in >> n.name >> n.parent >> n.vfp_table;
            nodes.push_back(std::move(n));
        } else if (tag == "well") {
            Well<Scalar> w;
            in >> w.name >> w.node >> w.vfp_table >> w.ipr_a >> w.ipr_b
               >> w.bhp_limit >> w.rate_limit >> w.guide >> w.q_start;
            wells.push_back(std::move(w));
        } else if (tag == "guess") {
            Scalar p;
            while (in >> p) {
                guess.push_back(p);
            }
        }
    }

    System<Scalar> system(props, phase);
    system.setTerminalPressure(terminal);
    system.setGroupTarget(target);
    for (auto& n : nodes) {
        system.addNode(std::move(n));
    }
    for (auto& w : wells) {
        system.addWell(std::move(w));
    }
    system.finish();
    return {std::move(system), std::move(guess)};
}


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
class ProductionSystem
{
public:
    using State = std::vector<Scalar>;
    using ScalarType = Scalar;
    static constexpr int NP = 3;   // water, oil, gas -- the order VFPPROD wants

    enum class Control { Thp, Bhp, OilRate };

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
    };

    ProductionSystem(const VFPProdProperties<Scalar>& props, const UnitSystem& units)
        : props_(&props), units_(&units)
    {}

    void addNode(Node n, const Scalar alq = 0.0) { nodes_.push_back(std::move(n)); branch_alq_.push_back(alq); }
    void addWell(Well w) { wells_.push_back(std::move(w)); }
    void setTerminalPressure(const Scalar p) { terminal_pressure_ = p; }
    void setRateScale(const Scalar s) { rate_scale_ = s; }

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
        if (rate_scale_ <= 0.0) {
            Scalar largest = 0.0;
            for (const auto& w : wells_) {
                largest = std::max(largest, w.oil_rate_limit);
            }
            rate_scale_ = std::max(largest * Scalar{0.01},
                                   unit::convert::from(1.0, unit::cubic(unit::meter) / unit::day));
        }
        controls_.assign(wells_.size(), Control::Thp);
    }

    int numNodes() const { return static_cast<int>(nodes_.size()) - 1; }
    int numWells() const { return static_cast<int>(wells_.size()); }
    int size() const { return 4 * numNodes() + 4 * numWells(); }

    const std::vector<Node>& nodes() const { return nodes_; }
    const std::vector<Well>& wells() const { return wells_; }
    Control control(const int w) const { return controls_[w]; }

    char controlLetter(const int w) const
    {
        switch (controls_[w]) {
        case Control::Thp:     return 'T';
        case Control::Bhp:     return 'B';
        case Control::OilRate: return 'O';
        }
        return '?';
    }

    int pIdx(const int node) const { return node - 1; }
    int qIdx(const int node, const int ph) const { return numNodes() + NP * (node - 1) + ph; }
    int qwIdx(const int w, const int ph) const { return 4 * numNodes() + NP * w + ph; }
    int bhpIdx(const int w) const { return 4 * numNodes() + NP * numWells() + w; }

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

    State residual(const State& x) const
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
            r[n - 1] = (hasTable(node)
                ? x[pIdx(n)] - tableBhp(node.vfp_table, upstream, branchRates(n), branch_alq_[n])
                : x[pIdx(n)] - upstream) / pressure_scale_;

            for (int ph = 0; ph < NP; ++ph) {
                Scalar balance = x[qIdx(n, ph)];
                for (const int c : children_[n]) {
                    balance -= x[qIdx(c, ph)];
                }
                for (const int w : wells_at_[n]) {
                    balance -= x[qwIdx(w, ph)];
                }
                r[nodes + NP * (n - 1) + ph] = balance / rate_scale_;
            }
        }

        for (int w = 0; w < wells; ++w) {
            const auto& well = wells_[w];
            const Scalar bhp = x[bhpIdx(w)];
            std::array<Scalar, NP> q{};
            for (int ph = 0; ph < NP; ++ph) {
                q[ph] = x[qwIdx(w, ph)];
                r[4 * nodes + NP * w + ph] = (q[ph] - ipr(well, ph, bhp)) / rate_scale_;
            }
            Scalar& control = r[4 * nodes + NP * wells + w];
            switch (controls_[w]) {
            case Control::Thp:
                control = (bhp - tableBhp(well.vfp_table, pressure(well.node), q, well.alq))
                        / pressure_scale_;
                break;
            case Control::Bhp:
                control = (bhp - well.bhp_limit) / pressure_scale_;
                break;
            case Control::OilRate:
                control = (q[1] - well.oil_rate_limit) / rate_scale_;
                break;
            }
        }
        return r;
    }

    /// Most restrictive wins, as for injection -- but a producer is limited by a
    /// bhp that is too *low*, not too high.
    bool updateControls(const State& x)
    {
        bool changed = false;
        for (int w = 0; w < numWells(); ++w) {
            const auto& well = wells_[w];
            auto wanted = Control::Thp;
            Scalar smallest = std::numeric_limits<Scalar>::max();
            auto consider = [&](const bool violated, const Scalar implied, const Control c) {
                if (violated && implied < smallest) {
                    smallest = implied;
                    wanted = c;
                }
            };
            consider(x[bhpIdx(w)] <= well.bhp_limit * (1.0 + 1e-9),
                     ipr(well, 1, well.bhp_limit), Control::Bhp);
            if (well.oil_rate_limit > Scalar{0}) {
                consider(x[qwIdx(w, 1)] > well.oil_rate_limit, well.oil_rate_limit,
                         Control::OilRate);
            }
            changed |= (wanted != controls_[w]);
            controls_[w] = wanted;
        }
        return changed;
    }

    State start(const State& node_pressure) const
    {
        State x(size(), 0.0);
        for (int n = 1; n <= numNodes(); ++n) {
            x[pIdx(n)] = node_pressure[n];
        }
        for (int w = 0; w < numWells(); ++w) {
            const auto& well = wells_[w];
            // Open a little above the bhp limit so the control test does not latch.
            x[bhpIdx(w)] = std::max(well.bhp_limit * Scalar{1.05}, node_pressure[well.node]);
            for (int ph = 0; ph < NP; ++ph) {
                x[qwIdx(w, ph)] = std::max(ipr(well, ph, x[bhpIdx(w)]), Scalar{0});
            }
        }
        for (int n = numNodes(); n >= 1; --n) {
            for (int ph = 0; ph < NP; ++ph) {
                Scalar q = 0.0;
                for (const int w : wells_at_[n]) {
                    q += x[qwIdx(w, ph)];
                }
                for (const int c : children_[n]) {
                    q += x[qIdx(c, ph)];
                }
                x[qIdx(n, ph)] = q;
            }
        }
        return x;
    }

    State pressures(const State& x) const
    {
        State p(nodes_.size(), terminal_pressure_);
        for (int n = 1; n <= numNodes(); ++n) {
            p[n] = x[pIdx(n)];
        }
        return p;
    }

    /// Oil rate per well, which is what a caller usually wants back.
    State wellRates(const State& x) const
    {
        State q(wells_.size());
        for (int w = 0; w < numWells(); ++w) {
            q[w] = x[qwIdx(w, 1)];
        }
        return q;
    }

    Scalar columnScale(const int i) const
    {
        const bool is_pressure = (i < numNodes()) || (i >= bhpIdx(0));
        return is_pressure ? pressure_scale_ : rate_scale_;
    }

    State limitStep(const State&, const State& dx) const { return dx; }

private:
    const VFPProdProperties<Scalar>* props_;
    const UnitSystem* units_;
    std::vector<Node> nodes_;
    std::vector<Scalar> branch_alq_;
    std::vector<Well> wells_;
    std::vector<std::vector<int>> children_;
    std::vector<std::vector<int>> wells_at_;
    std::vector<Control> controls_;

    Scalar terminal_pressure_ = 0.0;
    Scalar rate_scale_ = 0.0;
    Scalar pressure_scale_ = unit::barsa;
};

/// Take the Newton step as it comes. This is what the full system wants: it has
/// no kinks within an active set, so there is nothing for a globalisation to fix.
struct FullStep
{
    template<class Sys, class State>
    State accept(const Sys&, const State& x, const State&, const State& dx) const
    {
        State next(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) {
            next[i] = x[i] + dx[i];
        }
        return next;
    }
};

/// Backtrack until the residual norm drops. Not needed on the full system, but
/// useful when the residual is not smooth.
struct LineSearch
{
    int max_halvings = 12;

    template<class Sys, class State>
    State accept(const Sys& system, const State& x, const State& r, const State& dx) const
    {
        auto norm2 = [](const State& v) {
            double s = 0.0;
            for (const auto e : v) {
                s += e * e;
            }
            return std::sqrt(s);
        };
        const double f0 = norm2(r);
        double lambda = 1.0;
        State trial(x.size());
        for (int k = 0; k < max_halvings; ++k) {
            for (std::size_t i = 0; i < x.size(); ++i) {
                trial[i] = x[i] + lambda * dx[i];
            }
            if (norm2(system.residual(trial)) < f0) {
                return trial;
            }
            lambda *= 0.5;
        }
        return trial;
    }
};

/// Solve the system from a guess at the node pressures. The tolerance is on the
/// scaled residual, so it reads as bar on the pressure rows.
/// Ask a system whether it assembles its own Jacobian, without requiring that
/// every system knows how.
template<class Sys>
bool systemUsesAnalytic(const Sys& system)
{
    if constexpr (requires { system.usesAnalyticJacobian(); }) {
        return system.usesAnalyticJacobian();
    } else {
        return false;
    }
}

template<class Sys, class State>
auto systemJacobian(const Sys& system, const State& x)
{
    if constexpr (requires { system.jacobian(x); }) {
        return system.jacobian(x);
    } else {
        return DenseMatrix<typename Sys::ScalarType>(system.size());
    }
}

/// Solve any of the systems in this file. They differ in what a rate is -- one
/// number for an injection network, three for a production one -- but not in how
/// the Newton, the active set or the bounds work.
template<class Sys, class Globalisation = FullStep>
Result<typename Sys::ScalarType>
solve(Sys& system,
      const std::vector<typename Sys::ScalarType>& node_pressure_guess,
      const typename Sys::ScalarType tolerance = 1e-2,
      const int max_iterations = 50,
      Globalisation globalisation = {})
{
    using Scalar = typename Sys::ScalarType;
    auto x = system.start(node_pressure_guess);
    const int n = system.size();
    Result<Scalar> last;
    std::vector<std::string> trace;
    auto joined = [&trace] {
        std::string out;
        for (const auto& e : trace) {
            out += (out.empty() ? "" : " ") + e;
        }
        return out;
    };

    // Guide rates are explicit: the simulator sets them once per timestep, and
    // this follows that. Refreshing them inside the Newton makes each well's
    // share a moving target while its rate is chasing it, and the active set
    // then cycles between group and thp control instead of settling.
    if constexpr (requires { system.refreshGuides(x); }) {
        system.refreshGuides(x);
    }

    for (int it = 1; it <= max_iterations; ++it) {
        const bool controls_moved = system.updateControls(x);
        const auto r = system.residual(x);

        Scalar worst = 0.0;
        for (const auto e : r) {
            worst = std::max(worst, std::abs(e));
        }
        {   // remember the active set, so a cycle can be seen in the report
            std::string set;
            for (int w = 0; w < system.numWells(); ++w) {
                set += system.controlLetter(w);
            }
            trace.push_back(set);
            if (trace.size() > 8) {
                trace.erase(trace.begin());
            }
        }
        const bool settled = !controls_moved;
        if (worst < tolerance && settled) {
            return {true, it, system.pressures(x), system.wellRates(x), worst, false, false, {}};
        }
        last = {false, it, {}, {}, worst, controls_moved, false, joined()};

        // A system that can hand over an assembled Jacobian does; the rest are
        // differenced. The production prototype has no analytic one yet.
        DenseMatrix<Scalar> J = systemUsesAnalytic(system)
            ? systemJacobian(system, x)
            : [&] {
                  DenseMatrix<Scalar> fd(n);
                  for (int j = 0; j < n; ++j) {
                      auto shifted = x;
                      const Scalar h = 1e-2 * system.columnScale(j);
                      shifted[j] += h;
                      const auto rj = system.residual(shifted);
                      for (int i = 0; i < n; ++i) {
                          fd(i, j) = (rj[i] - r[i]) / h;
                      }
                  }
                  return fd;
              }();

        std::vector<Scalar> negative(n), dx;
        for (int i = 0; i < n; ++i) {
            negative[i] = -r[i];
        }
        if (!J.solve(negative, dx)) {
            last.node_pressure = system.pressures(x);
            last.well_rate = system.wellRates(x);
            return last;
        }

        dx = system.limitStep(x, dx);
        // The residual jumps when a control switches, and that jump is not a
        // failure to make progress. Letting a globalisation veto it stalls the
        // active set instead of resolving it.
        if (controls_moved) {
            for (int i = 0; i < n; ++i) {
                x[i] += dx[i];
            }
        } else {
            x = globalisation.accept(system, x, r, dx);
        }
    }
    last.iterations = max_iterations + 1;
    last.node_pressure = system.pressures(x);
    last.well_rate = system.wellRates(x);
    return last;
}

} // namespace Opm::NetworkSolve

#endif // OPM_NETWORK_SYSTEM_HEADER_INCLUDED
