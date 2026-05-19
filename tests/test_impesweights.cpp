/*
  Copyright 2025 Equinor.

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

#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <algorithm>
#include <array>
#include <cmath>

#define BOOST_TEST_MODULE ImpesWeightsTest
#include <boost/test/unit_test.hpp>

namespace {

constexpr double tolerance = 1.0e-12;

struct Weights {
    double water = 0.0;
    double oil = 0.0;
    double gas = 0.0;
};

Weights analyticBlackoilWeights(const double Bw,
                                const double Bo,
                                const double Bg,
                                const double rs,
                                const double rv,
                                const bool waterActive = true,
                                const bool oilActive = true,
                                const bool gasActive = true)
{
    Weights weights;
    const double denominator = (oilActive && gasActive) ? (1.0 - rs * rv) : 1.0;

    if (waterActive) {
        weights.water = Bw;
    }
    if (oilActive) {
        weights.oil = (Bo - rs * Bg) / denominator;
    }
    if (gasActive) {
        weights.gas = (Bg - rv * Bo) / denominator;
    }

    return weights;
}

Weights coatsBlackoilWeights(const double Bw,
                             const double Bo,
                             const double Bg,
                             const double rs,
                             const double rv,
                             const double rsw,
                             const double rvw,
                             const bool waterActive = true,
                             const bool oilActive = true,
                             const bool gasActive = true)
{
    Weights weights;
    const double denominator = 1.0 - rs * rv - rvw * rsw;

    if (waterActive) {
        weights.water = (Bw * (1.0 - rs * rv) - rsw * (Bg - rv * Bo)) / denominator;
    }
    if (oilActive) {
        weights.oil = (Bo * (1.0 - rvw * rsw) - rs * (Bg - rvw * Bw)) / denominator;
    }
    if (gasActive) {
        weights.gas = (Bg - rvw * Bw - rv * Bo) / denominator;
    }

    return weights;
}

double absMax(const Weights& weights)
{
    return std::max({std::fabs(weights.water), std::fabs(weights.oil), std::fabs(weights.gas)});
}

Weights normalizeByAbsMax(Weights weights)
{
    const double scale = absMax(weights);
    if (scale > 0.0) {
        weights.water /= scale;
        weights.oil /= scale;
        weights.gas /= scale;
    }

    return weights;
}

void checkStandardBlackoilBalance(const Weights& weights,
                                  const double Bw,
                                  const double Bo,
                                  const double Bg,
                                  const double rs,
                                  const double rv,
                                  const bool waterActive = true,
                                  const bool oilActive = true,
                                  const bool gasActive = true)
{
    if (waterActive) {
        BOOST_CHECK_CLOSE_FRACTION(weights.water / Bw, 1.0, tolerance);
    }
    if (oilActive) {
        BOOST_CHECK_CLOSE_FRACTION((weights.oil + rs * weights.gas) / Bo, 1.0, tolerance);
    }
    if (gasActive) {
        BOOST_CHECK_CLOSE_FRACTION((rv * weights.oil + weights.gas) / Bg, 1.0, tolerance);
    }
}

void checkExtendedBlackoilBalance(const Weights& weights,
                                  const double Bw,
                                  const double Bo,
                                  const double Bg,
                                  const double rs,
                                  const double rv,
                                  const double rsw,
                                  const double rvw,
                                  const bool waterActive = true,
                                  const bool oilActive = true,
                                  const bool gasActive = true)
{
    if (waterActive) {
        BOOST_CHECK_CLOSE_FRACTION((weights.water + rsw * weights.gas) / Bw, 1.0, tolerance);
    }
    if (oilActive) {
        BOOST_CHECK_CLOSE_FRACTION((weights.oil + rs * weights.gas) / Bo, 1.0, tolerance);
    }
    if (gasActive) {
        BOOST_CHECK_CLOSE_FRACTION((rvw * weights.water + rv * weights.oil + weights.gas) / Bg,
                                   1.0,
                                   tolerance);
    }
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(StandardBlackoilWeightsSatisfyVolumeBalance)
{
    const double Bw = 1.05;
    const double Bo = 1.20;
    const double Bg = 0.85;
    const double rs = 0.25;
    const double rv = 0.15;

    const auto weights = analyticBlackoilWeights(Bw, Bo, Bg, rs, rv);

    checkStandardBlackoilBalance(weights, Bw, Bo, Bg, rs, rv);
}

BOOST_AUTO_TEST_CASE(UndersaturatedLimitsReduceToExpectedFormulas)
{
    const double Bw = 1.02;
    const double Bo = 1.18;
    const double Bg = 0.81;

    {
        const double rs = 0.0;
        const double rv = 0.12;
        const auto weights = analyticBlackoilWeights(Bw, Bo, Bg, rs, rv);

        BOOST_CHECK_CLOSE_FRACTION(weights.water, Bw, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(weights.oil, Bo, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(weights.gas, Bg - rv * Bo, tolerance);
        checkStandardBlackoilBalance(weights, Bw, Bo, Bg, rs, rv);
    }

    {
        const double rs = 0.18;
        const double rv = 0.0;
        const auto weights = analyticBlackoilWeights(Bw, Bo, Bg, rs, rv);

        BOOST_CHECK_CLOSE_FRACTION(weights.water, Bw, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(weights.oil, Bo - rs * Bg, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(weights.gas, Bg, tolerance);
        checkStandardBlackoilBalance(weights, Bw, Bo, Bg, rs, rv);
    }
}

BOOST_AUTO_TEST_CASE(AnalyticAndCoatsMatchWithoutWaterGasCrossTerms)
{
    const double Bw = 1.08;
    const double Bo = 1.24;
    const double Bg = 0.79;
    const double rs = 0.22;
    const double rv = 0.11;

    const auto analytic = analyticBlackoilWeights(Bw, Bo, Bg, rs, rv);
    const auto coats = coatsBlackoilWeights(Bw, Bo, Bg, rs, rv, 0.0, 0.0);
    const auto normalizedCoats = normalizeByAbsMax(coats);
    const auto normalizedAnalytic = normalizeByAbsMax(analytic);

    BOOST_CHECK_CLOSE_FRACTION(coats.water, analytic.water, tolerance);
    BOOST_CHECK_CLOSE_FRACTION(coats.oil, analytic.oil, tolerance);
    BOOST_CHECK_CLOSE_FRACTION(coats.gas, analytic.gas, tolerance);

    BOOST_CHECK_CLOSE_FRACTION(normalizedCoats.water, normalizedAnalytic.water, tolerance);
    BOOST_CHECK_CLOSE_FRACTION(normalizedCoats.oil, normalizedAnalytic.oil, tolerance);
    BOOST_CHECK_CLOSE_FRACTION(normalizedCoats.gas, normalizedAnalytic.gas, tolerance);
}

BOOST_AUTO_TEST_CASE(ExtendedCoatsWeightsSatisfyVolumeBalance)
{
    const double Bw = 1.04;
    const double Bo = 1.19;
    const double Bg = 0.83;
    const double rs = 0.21;
    const double rv = 0.09;
    const double rsw = 0.04;
    const double rvw = 0.03;

    const auto weights = coatsBlackoilWeights(Bw, Bo, Bg, rs, rv, rsw, rvw);

    checkExtendedBlackoilBalance(weights, Bw, Bo, Bg, rs, rv, rsw, rvw);
}

BOOST_AUTO_TEST_CASE(SingleAndTwoPhaseLimitsRemainConsistent)
{
    {
        const double Bw = 1.07;
        const auto weights = coatsBlackoilWeights(Bw, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                  true, false, false);

        BOOST_CHECK_CLOSE_FRACTION(weights.water, Bw, tolerance);
        BOOST_CHECK_SMALL(weights.oil, tolerance);
        BOOST_CHECK_SMALL(weights.gas, tolerance);
        checkExtendedBlackoilBalance(weights, Bw, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                                     true, false, false);
    }

    {
        const double Bo = 1.23;
        const double Bg = 0.82;
        const double rs = 0.19;
        const double rv = 0.08;
        const auto weights = analyticBlackoilWeights(0.0, Bo, Bg, rs, rv,
                                                     false, true, true);

        BOOST_CHECK_SMALL(weights.water, tolerance);
        checkStandardBlackoilBalance(weights, 1.0, Bo, Bg, rs, rv,
                                     false, true, true);
    }
}