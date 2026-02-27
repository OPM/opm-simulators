/*
  Copyright 2026 SINTEF AS.

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

#ifndef OPM_NEWTON_ITERATION_CONTEXT_HPP
#define OPM_NEWTON_ITERATION_CONTEXT_HPP

#include <cassert>

namespace Opm {

/// \brief Context for iteration-dependent decisions in the Newton solver.
///
/// Provides explicit state for iteration-dependent behavior:
///
/// - Global iteration: The top-level Newton iteration count (0-based).
///   Used for NUPCOL checks, group controls, gas lift timing.
///
/// - Local iteration: Iteration count within a nested domain solve.
///   Starts at 0 for each domain.
///
/// - Timestep initialization: One-time setup at the start of each timestep,
///   independent of iteration number.
///
/// For domain-local solves, `forLocalSolve()` creates a context that
/// preserves the global iteration but resets local state.
///
struct NewtonIterationContext {

    NewtonIterationContext() = default;

    /// Getters

    /// Current iteration number (0-based).
    /// Returns the local iteration during domain-local solves,
    /// and the global iteration otherwise.
    int iteration() const
    {
        return inLocalSolve_ ? localIteration_ : globalIteration_;
    }

    /// Whether we are inside a domain-local solve (NLDD).
    bool inLocalSolve() const { return inLocalSolve_; }

    /// Whether timestep initialization has been performed.
    bool timestepInitialized() const { return timestepInitialized_; }

    /// Semantic Queries

    /// Should timestep initialization run?
    /// True only on first global iteration before initialization is done.
    bool needsTimestepInit() const
    {
        return !timestepInitialized_ && !inLocalSolve_;
    }

    /// Is this the first iteration of the global solve (not a local solve)?
    /// Use for one-time-per-timestep logic
    bool isFirstGlobalIteration() const
    {
        return globalIteration_ == 0 && !inLocalSolve_;
    }

    /// Are we within the NUPCOL iteration window?
    /// Always uses global iteration regardless of local solve state.
    bool withinNupcol(int nupcol) const
    {
        return globalIteration_ < nupcol;
    }

    /// Should tolerances be relaxed based on iteration count?
    /// Uses local iteration for local solves, global otherwise.
    bool shouldRelax(int strictIterations) const
    {
        return iteration() >= strictIterations;
    }

    /// Whether inner well iterations should run.
    /// Returns false during local solves, otherwise true while below maxIter.
    bool shouldRunInnerWellIterations(int maxIter) const
    {
        if (inLocalSolve_) return false;
        return globalIteration_ < maxIter;
    }

    /// State Mutations

    /// Mark timestep initialization as complete.
    void markTimestepInitialized()
    {
        timestepInitialized_ = true;
    }

    /// Advance the current iteration counter.
    /// Advances the local iteration during domain-local solves,
    /// and the global iteration otherwise.
    void advanceIteration()
    {
        if (inLocalSolve_) {
            localIteration_++;
        } else {
            globalIteration_++;
        }
    }

    /// Reset all state for a new timestep.
    void resetForNewTimestep()
    {
        globalIteration_ = 0;
        localIteration_ = 0;
        timestepInitialized_ = false;
        inLocalSolve_ = false;
    }

    /// Create a context for a domain-local solve.
    /// Preserves global iteration, resets local iteration to 0.
    NewtonIterationContext forLocalSolve() const
    {
        return {globalIteration_, 0, /*inLocalSolve=*/true, /*timestepInitialized=*/true};
    }

private:
    NewtonIterationContext(int globalIter, int localIter, bool localSolve, bool tsInit)
        : globalIteration_(globalIter), localIteration_(localIter)
        , inLocalSolve_(localSolve), timestepInitialized_(tsInit)
    {}
    int globalIteration_ = 0;
    int localIteration_ = 0;
    bool inLocalSolve_ = false;
    bool timestepInitialized_ = false;
};

/// RAII guard for NLDD domain-local iteration context.
/// Saves the current context on the problem, installs a local-solve
/// context, and restores the original on destruction.
template<class Problem>
class LocalContextGuard {
public:
    LocalContextGuard(Problem& problem)
        : problem_(problem)
        , saved_(problem.iterationContext())
    {
        problem_.mutableIterationContext() = saved_.forLocalSolve();
    }

    ~LocalContextGuard()
    {
        problem_.mutableIterationContext() = saved_;
    }

    LocalContextGuard(const LocalContextGuard&) = delete;
    LocalContextGuard& operator=(const LocalContextGuard&) = delete;
    LocalContextGuard(LocalContextGuard&&) = delete;
    LocalContextGuard& operator=(LocalContextGuard&&) = delete;

    NewtonIterationContext& context() { return problem_.mutableIterationContext(); }
    const NewtonIterationContext& context() const { return problem_.iterationContext(); }

private:
    Problem& problem_;
    NewtonIterationContext saved_;
};

} // namespace Opm

#endif // OPM_NEWTON_ITERATION_CONTEXT_HPP
