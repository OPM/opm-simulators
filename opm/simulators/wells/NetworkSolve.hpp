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
#ifndef OPM_NETWORK_SOLVE_HEADER_INCLUDED
#define OPM_NETWORK_SOLVE_HEADER_INCLUDED

#include <opm/input/eclipse/Units/Units.hpp>

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/fmatrix.hh>

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

constexpr int NoTable = 9999;   // GNETINJE's "no table": pressure passes through

struct Node
{
    std::string name;
    int parent = -1;            // -1 only for the terminal
    int vfp_table = NoTable;
    /// NEFAC: what this node passes on of what it collects.
    double efficiency = 1.0;
};

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
    /// Iterations on which some well changed control. A solve that has to move
    /// the active set a few times is working; one that keeps moving it is not.
    int switches = 0;
    /// Production only: every well's water/oil/gas and bhp at the solution.
    std::vector<std::array<Scalar, 3>> well_phase_rates;
    std::vector<Scalar> well_bhp;
};

/// Dense square system, solved by Dune. The networks this solves have tens of
/// unknowns, so a dense direct solve is the whole story. Wrapped only to keep
/// (i,j) indexing and to answer "singular" with false instead of an exception --
/// the caller hands the network back to the relaxed update rather than throwing
/// out of the well model.
template<class Scalar>
class DenseMatrix
{
public:
    explicit DenseMatrix(const int n) : a_(n, n, Scalar{0}) {}

    Scalar& operator()(const int i, const int j) { return a_[i][j]; }
    Scalar operator()(const int i, const int j) const { return a_[i][j]; }

    /// Solves A y = b. False if A is singular to working precision.
    bool solve(const std::vector<Scalar>& b, std::vector<Scalar>& y) const
    {
        const auto n = a_.N();
        Dune::DynamicVector<Scalar> rhs(n), x(n, Scalar{0});
        std::copy(b.begin(), b.end(), rhs.begin());
        try {
            a_.solve(x, rhs);
        } catch (const Dune::FMatrixError&) {
            return false;
        }
        y.assign(x.begin(), x.end());
        return true;
    }

private:
    Dune::DynamicMatrix<Scalar> a_;
};
/// Divide a group target by guide rate, take out the wells whose own limits keep
/// them below their share, and re-divide the rest among those that can take it.
/// A well that is out gets no share at all, so its own control binds.
///
/// This is the fixed point an active set would otherwise have to find by
/// iterating, and finding it here is what stops it cycling: a multiplier read
/// from the iterate means "the even split" while nobody is on group control and
/// "the remainder after the others' rates" while somebody is, and each of those
/// two numbers selects the state that produces the other.
template<class Scalar>
std::vector<Scalar> shareByGuide(const std::vector<Scalar>& guide,
                                 const std::vector<char>& in_group,
                                 const std::vector<Scalar>& own,
                                 const Scalar target)
{
    const int n = static_cast<int>(guide.size());
    std::vector<Scalar> share(n, std::numeric_limits<Scalar>::max());
    if (!(target > Scalar{0})) {
        return share;
    }

    std::vector<char> pooled(in_group);
    Scalar guides = 0.0;
    for (int w = 0; w < n; ++w) {
        if (pooled[w] != 0) {
            guides += guide[w];
        }
    }

    Scalar remaining = target;
    for (int pass = 0; pass <= n; ++pass) {
        if (!(guides > Scalar{0})) {
            break;
        }
        int drop = -1;
        Scalar worst = 0.0;
        for (int w = 0; w < n; ++w) {
            if (pooled[w] == 0) {
                continue;
            }
            share[w] = guide[w] / guides * std::max(remaining, Scalar{0});
            if (share[w] - own[w] > worst) {
                worst = share[w] - own[w];
                drop = w;
            }
        }
        if (drop < 0) {
            break;
        }
        pooled[drop] = 0;
        guides -= guide[drop];
        remaining -= own[drop];
        share[drop] = std::numeric_limits<Scalar>::max();
    }
    return share;
}


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

/// Convergence settings for solve(). No defaults: a caller states what it wants.
template<class Scalar>
struct Parameters
{
    /// Max norm of the scaled residual at which the system is converged.
    Scalar tolerance;
    int max_iterations;
};

/// Solve an InjectionSystem or a ProductionSystem by Newton-Raphson, choosing
/// each well's control by an active set as it goes. Takes the node pressures to
/// start from; returns the converged pressures and rates, or a Result with
/// converged false and the reason it stopped.
template<class Sys, class Globalisation>
Result<typename Sys::ScalarType>
solve(Sys& system,
      const std::vector<typename Sys::ScalarType>& node_pressure_guess,
      const Parameters<typename Sys::ScalarType> params,
      Globalisation globalisation)
{
    const auto tolerance = params.tolerance;
    const int max_iterations = params.max_iterations;
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

    // Only under --network-group-control, where the network places the group's
    // split itself: the share is each well's potential at the node pressure, so
    // it cannot be known before the starting pressures. Explicit, like OPM's own
    // guide rates, and for the same reason -- recomputing it every iteration
    // makes each share a moving target while its rate chases it, and the active
    // set cycles between group and thp control. Making it implicit (the share an
    // unknown of the system) or iterating it to a fixed point in an outer loop
    // are the ways past that; both are open.
    if constexpr (requires { system.refreshGuides(x); }) {
        system.refreshGuides(x);
    }

    int switches = 0;
    bool enforcing = false;
    for (int it = 1; it <= max_iterations; ++it) {
        const bool controls_moved = system.updateControls(x);
        switches += controls_moved ? 1 : 0;
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
            // Converged, but possibly with a well parked above its own rate
            // limit -- thp's allowance is capped by that limit while the solve
            // is moving, and a capped allowance ties with it. Now that the
            // iterate is not transient, drop the cap and carry on from here;
            // whoever is over the line goes on rate control and the rest take
            // it up. The `enforcing` flag below is what keeps this to one pass,
            // so a well that is still over the line afterwards converges as it
            // is rather than dropping the cap again.
            if constexpr (requires { system.setEnforceRateLimits(true); }) {
                if (!enforcing && system.rateLimitsViolated(x)) {
                    system.setEnforceRateLimits(true);
                    enforcing = true;
                    continue;
                }
            }
            Result<Scalar> done{true, it, system.pressures(x), system.wellRates(x), worst,
                                false, false, {}, switches};
            if constexpr (requires { system.wellPhaseRates(x); system.wellBhps(x); }) {
                done.well_phase_rates = system.wellPhaseRates(x);
                done.well_bhp = system.wellBhps(x);
            }
            return done;
        }
        last = {false, it, {}, {}, worst, controls_moved, false, joined(), switches};

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

#endif // OPM_NETWORK_SOLVE_HEADER_INCLUDED
