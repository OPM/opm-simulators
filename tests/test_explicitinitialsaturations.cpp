/*
  Copyright 2026, SINTEF Digital

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

#include "config.h"

#define BOOST_TEST_MODULE ExplicitInitialSaturations
#include <boost/test/unit_test.hpp>

#include <opm/simulators/flow/ExplicitInitialSaturations.hpp>

#include <array>
#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>

namespace
{

using Saturations = std::array<double, 3>;
using SuppliedSaturations = std::array<std::optional<double>, 3>;
using ActivePhases = std::array<bool, 3>;

constexpr ActivePhases threePhases {true, true, true};

void
checkSaturations(const SuppliedSaturations& supplied,
                 const ActivePhases& active,
                 const Saturations& expected)
{
    const auto result = Opm::completeExplicitSaturations(supplied, active);
    const auto& actual = result.saturations;
    for (std::size_t phase = 0; phase < actual.size(); ++phase) {
        BOOST_CHECK_SMALL(actual[phase] - expected[phase], 1.0e-12);
        BOOST_CHECK_GE(actual[phase], 0.0);
        BOOST_CHECK_LE(actual[phase], 1.0);
    }
    BOOST_CHECK_SMALL(actual[0] + actual[1] + actual[2] - 1.0, 1.0e-12);
}

} // namespace

BOOST_AUTO_TEST_CASE(AnyPairDefinesTheSameThreePhaseState)
{
    // Omitting a phase must preserve the supplied saturations, rather than
    // scaling those two values to fill the entire pore volume.
    const Saturations expected {0.2, 0.3, 0.5};
    checkSaturations({0.2, 0.3, std::nullopt}, threePhases, expected);
    checkSaturations({0.2, std::nullopt, 0.5}, threePhases, expected);
    checkSaturations({std::nullopt, 0.3, 0.5}, threePhases, expected);
}

BOOST_AUTO_TEST_CASE(AbsentAndPurePhases)
{
    checkSaturations({0.0, 1.0, std::nullopt}, threePhases, {0.0, 1.0, 0.0});
    checkSaturations({0.0, std::nullopt, 0.0}, threePhases, {0.0, 1.0, 0.0});
    checkSaturations({std::nullopt, 0.0, 0.0}, threePhases, {1.0, 0.0, 0.0});
}

BOOST_AUTO_TEST_CASE(AllThreeAreScaledWithinTheTolerance)
{
    const Saturations expected {0.2, 0.3, 0.5};
    checkSaturations({0.2, 0.3, 0.5}, threePhases, expected);
    checkSaturations({0.206, 0.309, 0.515}, threePhases, expected);
    checkSaturations({0.194, 0.291, 0.485}, threePhases, expected);

    // Include both tolerance boundaries despite binary roundoff in their sum.
    checkSaturations({0.21, 0.315, 0.525}, threePhases, expected);
    checkSaturations({0.19, 0.285, 0.475}, threePhases, expected);
}

BOOST_AUTO_TEST_CASE(AllThreeMustBeCloseToUnity)
{
    BOOST_CHECK_THROW(Opm::completeExplicitSaturations({0.21, 0.315, 0.5251}, threePhases),
                      std::invalid_argument);
    BOOST_CHECK_THROW(Opm::completeExplicitSaturations({0.1899, 0.285, 0.475}, threePhases),
                      std::invalid_argument);
    BOOST_CHECK_THROW(Opm::completeExplicitSaturations({0.0, 0.0, 0.0}, threePhases),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(OnlyNontrivialNormalizationNeedsAWarning)
{
    const auto exact = Opm::completeExplicitSaturations({0.2, 0.3, 0.5}, threePhases);
    BOOST_CHECK_EQUAL(exact.normalizationDeviation, 0.0);

    const auto rounded = Opm::completeExplicitSaturations(
        {double {0.2f}, double {0.3f}, double {0.5f}}, threePhases);
    BOOST_CHECK_EQUAL(rounded.normalizationDeviation, 0.0);
    const auto smallExcess = Opm::completeExplicitSaturations({1.0000005, 0.0, 0.0}, threePhases);
    BOOST_CHECK_EQUAL(smallExcess.normalizationDeviation, 0.0);

    const auto pair = Opm::completeExplicitSaturations({double {0.2f}, double {0.8f}, std::nullopt},
                                                       threePhases);
    BOOST_CHECK_EQUAL(pair.normalizationDeviation, 0.0);
    constexpr ActivePhases oilAndGas {false, true, true};
    const auto dry
        = Opm::completeExplicitSaturations({std::nullopt, double {0.4f}, double {0.6f}}, oilAndGas);
    BOOST_CHECK_EQUAL(dry.normalizationDeviation, 0.0);
}

BOOST_AUTO_TEST_CASE(NormalizationDeviationReportsTheInputClosureError)
{
    const auto excess = Opm::completeExplicitSaturations({0.206, 0.309, 0.515}, threePhases);
    BOOST_CHECK_SMALL(excess.normalizationDeviation - 0.03, 1.0e-12);
    const auto deficit = Opm::completeExplicitSaturations({0.194, 0.291, 0.485}, threePhases);
    BOOST_CHECK_SMALL(deficit.normalizationDeviation - 0.03, 1.0e-12);

    const auto boundary = Opm::completeExplicitSaturations({0.21, 0.315, 0.525}, threePhases);
    BOOST_CHECK_SMALL(boundary.normalizationDeviation - 0.05, 1.0e-12);
    const auto smallCorrection
        = Opm::completeExplicitSaturations({1.000002, 0.0, 0.0}, threePhases);
    BOOST_CHECK_SMALL(smallCorrection.normalizationDeviation - 2.0e-6, 1.0e-12);
}

BOOST_AUTO_TEST_CASE(TwoSuppliedPhasesCannotExceedThePoreVolumeBeyondInputPrecision)
{
    BOOST_CHECK_THROW(Opm::completeExplicitSaturations({0.6, 0.41, std::nullopt}, threePhases),
                      std::invalid_argument);
    BOOST_CHECK_THROW(Opm::completeExplicitSaturations({0.6, std::nullopt, 0.41}, threePhases),
                      std::invalid_argument);
    BOOST_CHECK_THROW(Opm::completeExplicitSaturations({std::nullopt, 0.6, 0.41}, threePhases),
                      std::invalid_argument);
    BOOST_CHECK_THROW(Opm::completeExplicitSaturations({0.6, 0.4000011, std::nullopt}, threePhases),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(RoundedPairInputsPreserveTheirSharesAndLeaveNoNegativePhase)
{
    // These are the double values obtained when reading single-precision
    // restart saturations. Their sums can exceed one by a few times 1e-8.
    checkSaturations({double {0.2f}, double {0.8f}, std::nullopt}, threePhases, {0.2, 0.8, 0.0});
    checkSaturations({double {0.2f}, std::nullopt, double {0.8f}}, threePhases, {0.2, 0.0, 0.8});
    checkSaturations({std::nullopt, double {0.4f}, double {0.6f}},
                     threePhases,
                     {0.0, 0.3999999940395357, 0.6000000059604643});
    checkSaturations({0.6, 0.40000000000001, std::nullopt}, threePhases, {0.6, 0.4, 0.0});

    // Include the input-tolerance boundary without introducing additional
    // rounding through a sum of two nonzero decimal values.
    checkSaturations({1.000001, 0.0, std::nullopt}, threePhases, {1.0, 0.0, 0.0});
}

BOOST_AUTO_TEST_CASE(ASmallPositiveRemainderIsPreserved)
{
    // The input tolerance corrects overfilled pore volumes; it does not erase
    // an explicitly implied positive oil or water inventory.
    checkSaturations({0.2, std::nullopt, 0.799999999}, threePhases, {0.2, 1.0e-9, 0.799999999});
    checkSaturations({std::nullopt, 0.2, 0.799999999}, threePhases, {1.0e-9, 0.2, 0.799999999});
}

BOOST_AUTO_TEST_CASE(InvalidActiveSaturationsAreRejected)
{
    for (const double invalid : {-0.01,
                                 std::numeric_limits<double>::infinity(),
                                 -std::numeric_limits<double>::infinity(),
                                 std::numeric_limits<double>::quiet_NaN()}) {
        BOOST_CHECK_THROW(Opm::completeExplicitSaturations({invalid, 0.3, 0.5}, threePhases),
                          std::invalid_argument);
        BOOST_CHECK_THROW(Opm::completeExplicitSaturations({0.2, invalid, 0.5}, threePhases),
                          std::invalid_argument);
        BOOST_CHECK_THROW(Opm::completeExplicitSaturations({0.2, 0.3, invalid}, threePhases),
                          std::invalid_argument);
    }

    // Individually finite values must also be rejected if their sum overflows.
    const double largest = std::numeric_limits<double>::max();
    BOOST_CHECK_THROW(
        Opm::completeExplicitSaturations({largest, largest, std::nullopt}, threePhases),
        std::invalid_argument);
    BOOST_CHECK_THROW(Opm::completeExplicitSaturations({largest, largest, largest}, threePhases),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(AtLeastTwoInputsAreNeededForThreePhases)
{
    BOOST_CHECK_THROW(Opm::completeExplicitSaturations({}, threePhases), std::invalid_argument);
    BOOST_CHECK_THROW(
        Opm::completeExplicitSaturations({0.2, std::nullopt, std::nullopt}, threePhases),
        std::invalid_argument);
    BOOST_CHECK_THROW(
        Opm::completeExplicitSaturations({std::nullopt, 0.3, std::nullopt}, threePhases),
        std::invalid_argument);
    BOOST_CHECK_THROW(
        Opm::completeExplicitSaturations({std::nullopt, std::nullopt, 0.5}, threePhases),
        std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(OneInputIsEnoughForTwoActivePhases)
{
    constexpr ActivePhases oilAndGas {false, true, true};
    checkSaturations({std::nullopt, 0.4, std::nullopt}, oilAndGas, {0.0, 0.4, 0.6});
    checkSaturations({std::nullopt, std::nullopt, 0.6}, oilAndGas, {0.0, 0.4, 0.6});
    checkSaturations({std::nullopt, 0.4, 0.6}, oilAndGas, {0.0, 0.4, 0.6});
    checkSaturations({std::nullopt, 0.384, 0.576}, oilAndGas, {0.0, 0.4, 0.6});
    checkSaturations({std::nullopt, 0.42, 0.63}, oilAndGas, {0.0, 0.4, 0.6});
    BOOST_CHECK_THROW(Opm::completeExplicitSaturations({}, oilAndGas), std::invalid_argument);
    BOOST_CHECK_THROW(
        Opm::completeExplicitSaturations({0.2, std::nullopt, std::nullopt}, oilAndGas),
        std::invalid_argument);
    BOOST_CHECK_THROW(
        Opm::completeExplicitSaturations({std::nullopt, 1.02, std::nullopt}, oilAndGas),
        std::invalid_argument);

    constexpr ActivePhases waterAndGas {true, false, true};
    checkSaturations({0.2, std::nullopt, std::nullopt}, waterAndGas, {0.2, 0.0, 0.8});
}

BOOST_AUTO_TEST_CASE(RoundedDryInputsPreserveTheirShares)
{
    constexpr ActivePhases oilAndGas {false, true, true};
    checkSaturations({std::nullopt, double {0.4f}, double {0.6f}},
                     oilAndGas,
                     {0.0, 0.3999999940395357, 0.6000000059604643});
    checkSaturations({std::nullopt, 0.333333, 0.666666}, oilAndGas, {0.0, 1.0 / 3.0, 2.0 / 3.0});
    checkSaturations({std::nullopt, 0.0, 1.000001}, oilAndGas, {0.0, 0.0, 1.0});
}

BOOST_AUTO_TEST_CASE(AllTwoPhasePairingsAreScaledWithinTheTolerance)
{
    // A fully supplied two-phase state uses the same closure tolerance as a
    // fully supplied three-phase state, regardless of which phase is inactive.
    constexpr std::array<ActivePhases, 3> pairs {
        {{false, true, true}, {true, false, true}, {true, true, false}}};
    constexpr std::array<std::array<double, 2>, 5> inputs {
        {{0.4, 0.6}, {0.388, 0.582}, {0.412, 0.618}, {0.38, 0.57}, {0.42, 0.63}}};

    for (const auto& active : pairs) {
        for (const auto& values : inputs) {
            SuppliedSaturations supplied;
            Saturations expected {};
            std::size_t input = 0;
            for (std::size_t phase = 0; phase < active.size(); ++phase) {
                if (active[phase]) {
                    supplied[phase] = values[input];
                    expected[phase] = input == 0 ? 0.4 : 0.6;
                    ++input;
                }
            }
            checkSaturations(supplied, active, expected);
        }
    }
}

BOOST_AUTO_TEST_CASE(AllTwoPhasePairingsRejectTotalsOutsideTheTolerance)
{
    constexpr std::array<ActivePhases, 3> pairs {
        {{false, true, true}, {true, false, true}, {true, true, false}}};
    constexpr std::array<std::array<double, 2>, 2> inputs {{{0.3799, 0.57}, {0.42, 0.6301}}};

    for (const auto& active : pairs) {
        for (const auto& values : inputs) {
            SuppliedSaturations supplied;
            std::size_t input = 0;
            for (std::size_t phase = 0; phase < active.size(); ++phase) {
                if (active[phase]) {
                    supplied[phase] = values[input++];
                }
            }
            BOOST_CHECK_THROW(Opm::completeExplicitSaturations(supplied, active),
                              std::invalid_argument);
        }
    }
}

BOOST_AUTO_TEST_CASE(TwoPhaseNormalizationReportsTheInputClosureError)
{
    const auto oilAndGas
        = Opm::completeExplicitSaturations({std::nullopt, 0.388, 0.582}, {false, true, true});
    BOOST_CHECK_SMALL(oilAndGas.normalizationDeviation - 0.03, 1.0e-12);
    const auto waterAndGas
        = Opm::completeExplicitSaturations({0.412, std::nullopt, 0.618}, {true, false, true});
    BOOST_CHECK_SMALL(waterAndGas.normalizationDeviation - 0.03, 1.0e-12);
    const auto waterAndOil
        = Opm::completeExplicitSaturations({0.42, 0.63, std::nullopt}, {true, true, false});
    BOOST_CHECK_SMALL(waterAndOil.normalizationDeviation - 0.05, 1.0e-12);
}

BOOST_AUTO_TEST_CASE(InactivePhasesDoNotAffectTheResult)
{
    constexpr ActivePhases oilAndGas {false, true, true};
    checkSaturations({0.8, 0.4, std::nullopt}, oilAndGas, {0.0, 0.4, 0.6});
    checkSaturations({-1.0, 0.4, std::nullopt}, oilAndGas, {0.0, 0.4, 0.6});
    checkSaturations(
        {std::numeric_limits<double>::quiet_NaN(), 0.4, std::nullopt}, oilAndGas, {0.0, 0.4, 0.6});
}

BOOST_AUTO_TEST_CASE(AOnePhaseCellIsFilledByThatPhase)
{
    constexpr ActivePhases waterOnly {true, false, false};
    checkSaturations({}, waterOnly, {1.0, 0.0, 0.0});
    checkSaturations({1.0, 0.1, 0.2}, waterOnly, {1.0, 0.0, 0.0});
    checkSaturations({0.999999, std::nullopt, std::nullopt}, waterOnly, {1.0, 0.0, 0.0});
    checkSaturations({1.000001, std::nullopt, std::nullopt}, waterOnly, {1.0, 0.0, 0.0});
    BOOST_CHECK_THROW(
        Opm::completeExplicitSaturations({0.9999989, std::nullopt, std::nullopt}, waterOnly),
        std::invalid_argument);
    BOOST_CHECK_THROW(
        Opm::completeExplicitSaturations({1.0000011, std::nullopt, std::nullopt}, waterOnly),
        std::invalid_argument);
    BOOST_CHECK_THROW(
        Opm::completeExplicitSaturations({0.9, std::nullopt, std::nullopt}, waterOnly),
        std::invalid_argument);
}
