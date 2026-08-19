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

#include <opm/simulators/wells/VFPInjProperties.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
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
    Scalar guide = 0.0;         // share of a group target
};

/// Which equation closes a well.
enum class Control { Thp, Bhp, Rate, Grup };

template<class Scalar>
struct Result
{
    bool converged = false;
    int iterations = 0;
    std::vector<Scalar> node_pressure;   // every node, terminal included
};

/// Dense square system. The networks this solves have tens of unknowns, so
/// Gaussian elimination with partial pivoting is the whole story.
template<class Scalar>
class DenseMatrix
{
public:
    explicit DenseMatrix(const int n) : n_(n), a_(n * n, 0.0) {}

    Scalar& operator()(const int i, const int j) { return a_[i * n_ + j]; }

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

    System(const VFPInjProperties<Scalar>& props, const Phase phase)
        : props_(&props), phase_(phase)
    {}

    void addNode(Node n) { nodes_.push_back(std::move(n)); }
    void addWell(Well<Scalar> w) { wells_.push_back(std::move(w)); }
    void setTerminalPressure(const Scalar p) { terminal_pressure_ = p; }
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
        controls_.assign(wells_.size(), group_target_ > 0.0 ? Control::Grup : Control::Thp);
    }

    int numNodes() const { return static_cast<int>(nodes_.size()) - 1; }
    int numWells() const { return static_cast<int>(wells_.size()); }
    bool grouped() const { return group_target_ > 0.0; }
    int size() const { return 2 * numNodes() + 2 * numWells() + (grouped() ? 1 : 0); }

    const std::vector<Node>& nodes() const { return nodes_; }
    const std::vector<Well<Scalar>>& wells() const { return wells_; }
    Control control(const int w) const { return controls_[w]; }

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
            injected += q;

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

            auto wanted = Control::Thp;
            Scalar smallest = std::numeric_limits<Scalar>::max();
            auto consider = [&](const bool violated, const Scalar implied, const Control c) {
                if (violated && implied < smallest) {
                    smallest = implied;
                    wanted = c;
                }
            };
            consider(x[bhpIdx(w)] > well.bhp_limit, ipr(well, well.bhp_limit), Control::Bhp);
            consider(q > well.rate_limit, well.rate_limit, Control::Rate);
            if (grouped()) {
                // Inclusive: at the solution the rate equals the share exactly,
                // and a strict test would flip the control every iteration.
                const Scalar share = well.guide * x[lambdaIdx()];
                consider(q >= share * (1.0 - 1e-9), share, Control::Grup);
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
            const Scalar guess = std::max(well.rate_limit * Scalar{0.1}, rate_scale_);
            x[bhpIdx(w)] = tableBhp(well.vfp_table, node_pressure[well.node], guess);
            x[qwIdx(w)] = std::clamp(ipr(well, x[bhpIdx(w)]), Scalar{0}, well.rate_limit);
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
    bool clamp_to_axes_ = false;
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
template<class Scalar, class Globalisation = FullStep>
Result<Scalar> solve(System<Scalar>& system,
                     const std::vector<Scalar>& node_pressure_guess,
                     const Scalar tolerance = 1e-2,
                     const int max_iterations = 50,
                     Globalisation globalisation = {})
{
    auto x = system.start(node_pressure_guess);
    const int n = system.size();

    for (int it = 1; it <= max_iterations; ++it) {
        const bool controls_moved = system.updateControls(x);
        const auto r = system.residual(x);

        Scalar worst = 0.0;
        for (const auto e : r) {
            worst = std::max(worst, std::abs(e));
        }
        if (worst < tolerance && !controls_moved) {
            return {true, it, system.pressures(x)};
        }

        DenseMatrix<Scalar> J(n);
        for (int j = 0; j < n; ++j) {
            auto shifted = x;
            const Scalar h = 1e-2 * system.columnScale(j);
            shifted[j] += h;
            const auto rj = system.residual(shifted);
            for (int i = 0; i < n; ++i) {
                J(i, j) = (rj[i] - r[i]) / h;
            }
        }

        std::vector<Scalar> negative(n), dx;
        for (int i = 0; i < n; ++i) {
            negative[i] = -r[i];
        }
        if (!J.solve(negative, dx)) {
            return {false, it, system.pressures(x)};
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
    return {false, max_iterations + 1, system.pressures(x)};
}

} // namespace Opm::NetworkSolve

#endif // OPM_NETWORK_SYSTEM_HEADER_INCLUDED
