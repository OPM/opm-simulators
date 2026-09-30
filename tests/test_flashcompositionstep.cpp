/*
  Copyright 2026 SINTEF Digital

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM. If not, see <http://www.gnu.org/licenses/>.
*/

#include "config.h"

#define BOOST_TEST_MODULE FlashCompositionStep
#include <boost/test/unit_test.hpp>

#include <opm/models/ptflash/flashcompositionstep.hh>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <random>

namespace {

constexpr double compositionFloor = 1.0e-8;
constexpr double hydrocarbonFloor = 1.0e-8;
constexpr double maxAmountChange = 0.2;

template <std::size_t n>
struct State
{
    std::array<double, n> z;
    double sw;
};

template <std::size_t n>
State<n> step(State<n> state, const std::array<double, n>& dz, double dSw)
{
    Opm::applyFlashCompositionStep(state.z, state.sw, dz, dSw,
                                   compositionFloor, hydrocarbonFloor, maxAmountChange);
    return state;
}

double hydrocarbonShare(double sw)
{
    return std::max(1.0 - sw, hydrocarbonFloor);
}

// Largest change of the hydrocarbon share and of the component amounts.
template <std::size_t n>
double amountChange(const State<n>& oldState, const State<n>& newState)
{
    const double oldShare = hydrocarbonShare(oldState.sw);
    const double newShare = hydrocarbonShare(newState.sw);
    double change = std::abs(newShare - oldShare);
    for (std::size_t c = 0; c < n; ++c) {
        change = std::max(std::abs(newShare * newState.z[c] - oldShare * oldState.z[c]), change);
    }
    return change;
}

template <std::size_t n>
void checkBounds(const State<n>& oldState, const State<n>& newState)
{
    BOOST_CHECK_SMALL(std::accumulate(newState.z.begin(), newState.z.end(), -1.0), 1.0e-12);
    for (const double zc : newState.z) {
        BOOST_CHECK_GE(zc, compositionFloor * (1.0 - 1.0e-6));
    }
    BOOST_CHECK_GE(newState.sw, 0.0);
    BOOST_CHECK_LE(newState.sw, 1.0);
    BOOST_CHECK_LE(amountChange(oldState, newState), maxAmountChange * (1.0 + 1.0e-12));
}

} // namespace

BOOST_AUTO_TEST_CASE(SmallStepIsTakenInFull)
{
    const State<3> oldState{{0.5, 0.3, 0.2}, 0.3};
    const std::array dz{0.01, -0.02, 0.01};
    const auto newState = step(oldState, dz, 0.05);

    // The amounts h z follow the linearized step, so z moves by h_old/h_new of dz.
    BOOST_CHECK_CLOSE(newState.sw, 0.35, 1.0e-12);
    for (std::size_t c = 0; c < 3; ++c) {
        BOOST_CHECK_CLOSE(newState.z[c], oldState.z[c] + 0.7 / 0.65 * dz[c], 1.0e-10);
    }
}

BOOST_AUTO_TEST_CASE(LargeStepIsDampedInTheAmounts)
{
    // The hydrocarbon share would drop by 0.4, so the step is halved.
    const State<3> oldState{{0.4, 0.3, 0.3}, 0.2};
    const auto newState = step(oldState, {0.4, -0.2, -0.2}, 0.4);

    BOOST_CHECK_CLOSE(newState.sw, 0.4, 1.0e-10);
    BOOST_CHECK_CLOSE(newState.z[0], 2.0 / 3.0, 1.0e-10);
    BOOST_CHECK_CLOSE(newState.z[1], 1.0 / 6.0, 1.0e-10);
    BOOST_CHECK_CLOSE(newState.z[2], 1.0 / 6.0, 1.0e-10);
    BOOST_CHECK_CLOSE(amountChange(oldState, newState), maxAmountChange, 1.0e-10);
}

BOOST_AUTO_TEST_CASE(CompositionFloorDoesNotLengthenTheStep)
{
    // The candidate step takes the second component below zero. Keeping it at the floor
    // must not move the other amounts by more than the limit.
    const State<3> oldState{{0.98, 0.01, 0.01}, 0.0};
    const auto newState = step(oldState, {-0.004, -0.198, 0.202}, 0.2);

    checkBounds(oldState, newState);
    BOOST_CHECK_CLOSE(amountChange(oldState, newState), maxAmountChange, 1.0e-8);
    BOOST_CHECK_LT(newState.z[1], oldState.z[1]);
    BOOST_CHECK_GT(newState.sw, 0.0);
}

BOOST_AUTO_TEST_CASE(CompositionFloorDoesNotLengthenTheStepWithoutWater)
{
    const State<4> oldState{{0.62, 0.08, 0.02, 0.28}, 0.0};
    const auto newState = step(oldState, {0.12, 0.16, -0.08, -0.20}, 0.0);

    checkBounds(oldState, newState);
    BOOST_CHECK_CLOSE(amountChange(oldState, newState), maxAmountChange, 1.0e-8);
    BOOST_CHECK_EQUAL(newState.sw, 0.0);
}

BOOST_AUTO_TEST_CASE(SaturationBoundDoesNotLengthenTheStep)
{
    // Sw is clamped at zero, so the hydrocarbon share grows less than the damping
    // assumed and no longer offsets the decline of the third component.
    const State<3> oldState{{0.117, 0.021, 0.862}, 0.009};
    const auto newState = step(oldState, {0.257, 0.263, -0.52}, -0.288);

    checkBounds(oldState, newState);
    BOOST_CHECK_CLOSE(amountChange(oldState, newState), maxAmountChange, 1.0e-8);
    BOOST_CHECK_LT(newState.sw, oldState.sw);
}

BOOST_AUTO_TEST_CASE(InitiallyAbsentComponentReachesTheFloor)
{
    // An initial composition may hold zeros. Shortening the step must not leave the
    // absent component between zero and the floor.
    const State<4> oldState{{0.99, 0.0, 0.005, 0.005}, 0.0};
    const auto newState = step(oldState, {-0.2, -0.2, 0.2, 0.2}, 0.0);

    checkBounds(oldState, newState);
    BOOST_CHECK_CLOSE(newState.z[1], compositionFloor, 1.0e-6);
}

BOOST_AUTO_TEST_CASE(SaturationFollowsTheFlooredShare)
{
    // The proposed step reaches the hydrocarbon floor. After shortening, composition
    // and Sw must use the same final share so the accepted amount changes stay within the limit.
    const State<3> oldState{{1.0 - 2.0e-8, 1.0e-8, 1.0e-8}, 0.7999999905};
    const auto newState = step(oldState, {-2.0e-8, 4.0e-8, -2.0e-8}, 0.2);

    checkBounds(oldState, newState);
}

BOOST_AUTO_TEST_CASE(WaterFilledCellTakesTheInflowComposition)
{
    // The share h sits at its floor. A linearized balance that puts 0.1 of the pore
    // space with composition zIn into the cell asks for a step in z of order 1/h.
    const State<3> oldState{{0.2, 0.3, 0.5}, 1.0};
    const std::array zIn{0.6, 0.3, 0.1};
    const double dShare = 0.1;
    std::array<double, 3> dz{};
    for (std::size_t c = 0; c < 3; ++c) {
        dz[c] = (hydrocarbonFloor + dShare) * (zIn[c] - oldState.z[c]) / hydrocarbonFloor;
    }
    const auto newState = step(oldState, dz, -dShare);

    BOOST_CHECK_CLOSE(newState.sw, 1.0 - dShare, 1.0e-10);
    for (std::size_t c = 0; c < 3; ++c) {
        BOOST_CHECK_CLOSE(newState.z[c], zIn[c], 1.0e-4);
    }
}

BOOST_AUTO_TEST_CASE(ComponentAtTheFloorDoesNotHoldBackTheOthers)
{
    const State<3> oldState{{compositionFloor, 0.5, 0.5 - compositionFloor}, 0.0};
    const std::array dz{-1.0e-9, 0.1, -0.1 + 1.0e-9};
    const auto newState = step(oldState, dz, 0.0);

    checkBounds(oldState, newState);
    BOOST_CHECK_CLOSE(newState.z[0], compositionFloor, 1.0e-6);
    BOOST_CHECK_CLOSE(newState.z[1], oldState.z[1] + dz[1], 1.0e-6);
    BOOST_CHECK_CLOSE(newState.z[2], oldState.z[2] + dz[2], 1.0e-6);
}

BOOST_AUTO_TEST_CASE(RandomStepsKeepTheBounds)
{
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    std::uniform_int_distribution<std::size_t> pick(0, 3);
    constexpr std::array stepSizes{1.0e-8, 0.3, 3.0, 1.0e3};
    constexpr std::array exponents{1.0, 4.0, 20.0};
    // Includes a share just above maxAmountChange, which a damped step can take to its floor.
    constexpr std::array saturations{0.0, 1.0 - maxAmountChange - hydrocarbonFloor,
                                     1.0 - 1.0e-9, 1.0};
    constexpr std::size_t n = 4;

    for (int sample = 0; sample < 10000; ++sample) {
        State<n> oldState{};
        const double exponent = exponents[pick(gen) % exponents.size()];
        for (auto& zc : oldState.z) {
            zc = std::pow(uniform(gen), exponent);
        }
        oldState.z[pick(gen)] = 0.0;
        // Half of the compositions keep their zeros, as initial compositions may.
        const double zMin = pick(gen) < 2 ? compositionFloor : 0.0;
        const double zSum = std::accumulate(oldState.z.begin(), oldState.z.end(), 0.0);
        for (auto& zc : oldState.z) {
            zc = zMin + (1.0 - n * zMin) * zc / zSum;
        }
        oldState.sw = pick(gen) == 0 ? uniform(gen) : saturations[pick(gen)];

        std::normal_distribution<double> zStep(0.0, stepSizes[pick(gen)]);
        std::array<double, n> dz{};
        for (auto& dzc : dz) {
            dzc = zStep(gen);
        }
        const double dzMean = std::accumulate(dz.begin(), dz.end(), 0.0) / n;
        for (auto& dzc : dz) {
            dzc -= dzMean;
        }
        std::normal_distribution<double> swStep(0.0, stepSizes[pick(gen) % 3]);

        checkBounds(oldState, step(oldState, dz, swStep(gen)));
    }
}
