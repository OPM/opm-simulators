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

#ifndef OPM_NETWORK_NODE_PRESSURE_UPDATER_HPP
#define OPM_NETWORK_NODE_PRESSURE_UPDATER_HPP

#include <opm/input/eclipse/Units/Units.hpp>

#include <algorithm>
#include <cmath>
#include <optional>
#include <utility>

namespace Opm {

/// Drives one network node's pressure towards the fixed point of
/// r(P) = P_computed(P) - P, where P_computed is the pressure the network gives for
/// the rates the wells produce/inject with P applied as their THP.
///
/// One of these per node, held by the caller in a map and stepped once per node per
/// network sub-iteration: that is what "per-node" means here. Each node keeps its own
/// bracket and moves on its own evidence, which is also the limitation -- the nodes
/// are coupled through the wells below them, so a node's r changes when a sibling
/// moves, and a bracket end has to be dropped when the new evidence contradicts it
/// (below). Solving every node at once instead is what NetworkSolve::solve() does.
///
/// The well response makes r a decreasing but strongly nonlinear function of P: flat
/// while the wells are on group/rate control (r' = -1) and steep while they are on
/// THP control (r' ~ -10 for gas injectors, since BHP is pinned by the reservoir).
/// A plain damped fixed-point iteration either crawls or, once the THP branch is hit,
/// falls into a limit cycle. Hence: keep a sign-change bracket [lo (r>0), hi (r<0)]
/// per node and take Illinois regula-falsi steps inside it; before a bracket exists,
/// take a secant step from the last two evaluations, clipped to [P, P_computed]
/// (which contains the fixed point when r is decreasing); if neither is possible,
/// fall back to the damped, capped update.
template<typename Scalar>
class NodePressureUpdater
{
public:
    /// Next pressure to apply, given the pressure that was applied and what the
    /// network computed for the resulting rates. `valid` is false when the network
    /// had no VFP solution (rate beyond what the branches can deliver): the pressure
    /// is then treated as too high, r = 1 atm - P.
    ///
    /// `plateau_floor` is the lowest wellhead pressure among the node's wells that are
    /// currently not on THP control: above it the rates do not depend on the node
    /// pressure, below it the wells become THP-limited and the response changes
    /// abruptly. An unbracketed downward step is not taken past it (it stops just below,
    /// so the next evaluation is on the THP branch close to the kink).
    Scalar next(const Scalar applied,
                const Scalar computed,
                const bool valid,
                const Scalar damping_factor,
                const Scalar max_update,
                const std::optional<Scalar>& plateau_floor = std::nullopt)
    {
        const Scalar residual = valid ? computed - applied : Scalar{unit::atm} - applied;

        // Maintain the sign-change bracket. r is decreasing in P, so a point with r > 0
        // must lie below every point with r < 0: an evaluation that contradicts a
        // stored end (the wells changed state, e.g. stopped/reopened, or other nodes
        // moved) invalidates that end.
        const Scalar eps = 1e-3 * unit::barsa;
        if (residual > 0.0) {
            if (hi_ && applied >= hi_->first - eps) {
                hi_.reset();
                same_side_ = 0;
            }
            lo_ = std::make_pair(applied, residual);
            same_side_ = (last_side_ == +1) ? same_side_ + 1 : 1;
            last_side_ = +1;
            r_lo_scale_ = 1.0;
        } else if (residual < 0.0) {
            if (lo_ && applied <= lo_->first + eps) {
                lo_.reset();
                same_side_ = 0;
            }
            hi_ = std::make_pair(applied, residual);
            same_side_ = (last_side_ == -1) ? same_side_ + 1 : 1;
            last_side_ = -1;
            r_hi_scale_ = 1.0;
        }

        // A bracket that has collapsed, or an end whose weight has been halved away,
        // means the retained end is stale: drop it and continue from the newest point.
        if (lo_ && hi_) {
            const bool collapsed = hi_->first - lo_->first < 0.01 * unit::barsa;
            if (collapsed || r_lo_scale_ < 1.0 / 32.0 || r_hi_scale_ < 1.0 / 32.0) {
                if (last_side_ == +1) {
                    hi_.reset();
                } else {
                    lo_.reset();
                }
                same_side_ = 0;
                r_lo_scale_ = r_hi_scale_ = 1.0;
            }
        }

        Scalar target;
        if (lo_ && hi_) {
            // Illinois regula falsi: secant through the bracket ends; every time an end
            // is retained while the other moves, its residual weight is halved so the
            // iteration cannot stall on it.
            if (same_side_ >= 2) {
                if (last_side_ == +1) {
                    r_hi_scale_ *= 0.5;
                } else {
                    r_lo_scale_ *= 0.5;
                }
            }
            const Scalar r_lo = lo_->second * r_lo_scale_;
            const Scalar r_hi = hi_->second * r_hi_scale_;
            const Scalar p_lo = lo_->first;
            const Scalar p_hi = hi_->first;
            target = p_lo + r_lo * (p_hi - p_lo) / (r_lo - r_hi);
            // Stay strictly inside the bracket.
            const Scalar margin = 0.02 * (p_hi - p_lo);
            target = std::clamp(target, p_lo + margin, p_hi - margin);
        } else {
            bool have_secant = false;
            if (last_ && valid) {
                const auto [prev_applied, prev_residual] = *last_;
                const Scalar dp = applied - prev_applied;
                const Scalar dr = residual - prev_residual;
                if (std::abs(dp) > 1e-3 * unit::barsa && dr * dp < 0.0) {
                    const Scalar lo = std::min(applied, computed);
                    const Scalar hi = std::max(applied, computed);
                    target = std::clamp(applied - residual / (dr / dp), lo, hi);
                    have_secant = true;
                }
            }
            if (!have_secant) {
                target = damped(applied, residual, damping_factor, max_update);
            }
            // Unbracketed steps are guesses about a response that may be a step function
            // (wells switching control, stopping); keep them moderate, and do not jump
            // from the group-controlled plateau past the wells' own THP.
            if (plateau_floor && valid && applied > *plateau_floor && computed < *plateau_floor) {
                // Flat down to the floor: go straight to just below the kink.
                target = *plateau_floor - Scalar{0.5 * unit::barsa};
            } else {
                const Scalar cap = std::max(Scalar{0.25} * std::abs(applied), Scalar{10.0 * unit::barsa});
                target = std::clamp(target, applied - cap, applied + cap);
            }
        }
        last_ = std::make_pair(applied, residual);
        last_residual_ = std::abs(residual);
        return target;
    }

    /// How far the fixed point may still be from the last target: the bracket width
    /// while a bracket exists, else the last residual (r' ~ -1 at best, so |r| bounds
    /// the pressure error from below). Use this as the convergence measure instead
    /// of the raw residual, which is amplified by the well response slope.
    Scalar error() const
    {
        if (lo_ && hi_) {
            return hi_->first - lo_->first;
        }
        return last_residual_;
    }

    static Scalar damped(const Scalar applied, const Scalar residual,
                         const Scalar damping_factor, const Scalar max_update)
    {
        const Scalar change = std::min(damping_factor * std::abs(residual), max_update);
        return applied + (residual > 0 ? change : -change);
    }

private:
    std::optional<std::pair<Scalar, Scalar>> lo_;   // applied pressure with r > 0
    std::optional<std::pair<Scalar, Scalar>> hi_;   // applied pressure with r < 0
    std::optional<std::pair<Scalar, Scalar>> last_;
    int last_side_{0};   // +1: lo_ updated last, -1: hi_ updated last
    int same_side_{0};   // consecutive updates of the same bracket end
    Scalar r_lo_scale_{1.0};  // Illinois weights of the retained ends
    Scalar r_hi_scale_{1.0};
    Scalar last_residual_{0.0};
};

} // namespace Opm

#endif // OPM_NETWORK_NODE_PRESSURE_UPDATER_HPP
