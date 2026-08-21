// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/**
 * \file
 *
 * \brief The ODE integrator and phase-pressure function used to solve the
 *        hydrostatic equilibrium problem, shared by the black-oil and the
 *        compositional equilibration facilities.
 */
#ifndef OPM_EQUIL_PRESSURE_FUNCTION_HPP
#define OPM_EQUIL_PRESSURE_FUNCTION_HPP

#include <array>
#include <cassert>
#include <memory>
#include <vector>

namespace Opm {
namespace EQUIL {
namespace Details {

/// Fourth-order Runge-Kutta integrator for the initial-value problem
///   y'(x) = f(x, y),  y(span[0]) = y0
/// on N equidistant steps covering span, with O(h^3) dense output.
template <class Scalar, class RHS>
class RK4IVP
{
public:
    RK4IVP(const RHS& f,
           const std::array<Scalar,2>& span,
           const Scalar y0,
           const int N)
        : N_(N)
        , span_(span)
    {
        const Scalar h = stepsize();
        const Scalar h2 = h / 2;
        const Scalar h6 = h / 6;

        y_.reserve(N + 1);
        f_.reserve(N + 1);

        y_.push_back(y0);
        f_.push_back(f(span_[0], y0));

        for (int i = 0; i < N; ++i) {
            const Scalar x = span_[0] + i*h;
            const Scalar y = y_.back();

            const Scalar k1 = f_[i];
            const Scalar k2 = f(x + h2, y + h2*k1);
            const Scalar k3 = f(x + h2, y + h2*k2);
            const Scalar k4 = f(x + h, y + h*k3);

            y_.push_back(y + h6*(k1 + 2*(k2 + k3) + k4));
            f_.push_back(f(x + h, y_.back()));
        }

        assert (y_.size() == typename std::vector<Scalar>::size_type(N + 1));
    }

    Scalar operator()(const Scalar x) const
    {
        // Dense output (O(h**3)) according to Shampine
        // (Hermite interpolation)
        const Scalar h = stepsize();
        int i = (x - span_[0]) / h;
        // Crude handling of evaluation point outside "span_";
        if (i  <  0) { i = 0;      }
        if (N_ <= i) { i = N_ - 1; }

        // Relative to the interval actually used, so that a point at the end of
        // the span lands on t = 1 of the final interval rather than t = 0.
        const Scalar t = (x - (span_[0] + i*h)) / h;

        const Scalar y0 = y_[i], y1 = y_[i + 1];
        const Scalar f0 = f_[i], f1 = f_[i + 1];

        Scalar u = (1 - 2*t) * (y1 - y0);
        u += h * ((t - 1)*f0 + t*f1);
        u *= t * (t - 1);
        u += (1 - t)*y0 + t*y1;

        return u;
    }

private:
    int N_;
    std::array<Scalar,2> span_;
    std::vector<Scalar>  y_;
    std::vector<Scalar>  f_;

    Scalar stepsize() const
    { return (span_[1] - span_[0]) / N_; }
};

/// Phase pressure as a function of depth, obtained by integrating the
/// hydrostatic ODE
///   dp/ddepth = ODE(depth, p)
/// upwards and downwards from an initial condition (depth, pressure).
template <class Scalar, class ODE>
class PressureFunction
{
public:
    using VSpan = std::array<Scalar, 2>;

    struct InitCond {
        Scalar depth;
        Scalar pressure;
    };

    explicit PressureFunction(const ODE&      ode,
                              const InitCond& ic,
                              const int       nsample,
                              const VSpan&    span)
        : initial_(ic)
    {
        this->value_[Direction::Up] = std::make_unique<Distribution>
            (ode, VSpan {{ ic.depth, span[0] }}, ic.pressure, nsample);

        this->value_[Direction::Down] = std::make_unique<Distribution>
            (ode, VSpan {{ ic.depth, span[1] }}, ic.pressure, nsample);
    }

    PressureFunction(const PressureFunction& rhs)
        : initial_(rhs.initial_)
    {
        this->value_[Direction::Up] =
            std::make_unique<Distribution>(*rhs.value_[Direction::Up]);

        this->value_[Direction::Down] =
            std::make_unique<Distribution>(*rhs.value_[Direction::Down]);
    }

    PressureFunction(PressureFunction&& rhs) = default;

    PressureFunction& operator=(const PressureFunction& rhs)
    {
        this->initial_ = rhs.initial_;

        this->value_[Direction::Up] =
            std::make_unique<Distribution>(*rhs.value_[Direction::Up]);

        this->value_[Direction::Down] =
            std::make_unique<Distribution>(*rhs.value_[Direction::Down]);

        return *this;
    }

    PressureFunction& operator=(PressureFunction&& rhs)
    {
        this->initial_ = rhs.initial_;
        this->value_   = std::move(rhs.value_);

        return *this;
    }

    Scalar value(const Scalar depth) const
    {
        if (depth < this->initial_.depth) {
            // Value above initial condition depth.
            return (*this->value_[Direction::Up])(depth);
        }
        else if (depth > this->initial_.depth) {
            // Value below initial condition depth.
            return (*this->value_[Direction::Down])(depth);
        }
        else {
            // Value *at* initial condition depth.
            return this->initial_.pressure;
        }
    }

private:
    enum Direction : std::size_t { Up, Down, NumDir };

    using Distribution = RK4IVP<Scalar, ODE>;
    using DistrPtr = std::unique_ptr<Distribution>;

    InitCond initial_;
    std::array<DistrPtr, Direction::NumDir> value_;
};

} // namespace Details
} // namespace EQUIL
} // namespace Opm

#endif // OPM_EQUIL_PRESSURE_FUNCTION_HPP
