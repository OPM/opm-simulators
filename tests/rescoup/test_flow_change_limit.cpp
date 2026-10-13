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

#include <config.h>

#define BOOST_TEST_MODULE RescoupFlowChangeLimit
#include <boost/test/unit_test.hpp>

#include <opm/input/eclipse/Parser/ParserKeywords/G.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCoupling.hpp>

#include <cmath>
#include <limits>
#include <vector>

using Opm::ReservoirCoupling::flowChangeTimeStepLimit;

namespace {

constexpr double day = 86400.0;
constexpr double infinity = std::numeric_limits<double>::infinity();
// GRUPMAST item 4 not given
constexpr double no_limit = Opm::ParserKeywords::GRUPMAST::LIMITING_FRACTION::defaultValue;

// Flows are: reservoir production rate, then oil, gas and water injection rates

} // namespace

BOOST_AUTO_TEST_CASE(LimitFromLargestChange)
{
    // The production flow changes by 1%, the gas injection by 2%. The step length times
    // the fraction divided by the largest change: 10 days * 0.005 / 0.02 = 2.5 days.
    const auto result = flowChangeTimeStepLimit(0.005, 10 * day, 0.1 * day,
                                                std::vector<double>{100.0, 0.0, 50.0, 0.0},
                                                std::vector<double>{101.0, 0.0, 51.0, 0.0});
    BOOST_CHECK_CLOSE(result.max_change, 0.02, 1e-8);
    BOOST_CHECK_CLOSE(result.time_step_limit, 2.5 * day, 1e-8);
}

BOOST_AUTO_TEST_CASE(DecreaseCountsAsChange)
{
    const auto result = flowChangeTimeStepLimit(0.01, 10 * day, 0.1 * day,
                                                std::vector<double>{100.0, 0.0, 0.0, 0.0},
                                                std::vector<double>{98.0, 0.0, 0.0, 0.0});
    BOOST_CHECK_CLOSE(result.max_change, 0.02, 1e-8);
    BOOST_CHECK_CLOSE(result.time_step_limit, 5 * day, 1e-8);
}

BOOST_AUTO_TEST_CASE(NoChangeDoesNotLimit)
{
    const auto result = flowChangeTimeStepLimit(0.005, 10 * day, 0.1 * day,
                                                std::vector<double>{100.0, 0.0, 50.0, 0.0},
                                                std::vector<double>{100.0, 0.0, 50.0, 0.0});
    BOOST_CHECK_EQUAL(result.max_change, 0.0);
    BOOST_CHECK_EQUAL(result.time_step_limit, infinity);
}

BOOST_AUTO_TEST_CASE(NeverBelowMinimumStep)
{
    // A change of 100% gives 10 days * 0.005 / 1 = 0.05 days, below the minimum of 0.5 days
    const auto result = flowChangeTimeStepLimit(0.005, 10 * day, 0.5 * day,
                                                std::vector<double>{100.0, 0.0, 0.0, 0.0},
                                                std::vector<double>{200.0, 0.0, 0.0, 0.0});
    BOOST_CHECK_EQUAL(result.time_step_limit, 0.5 * day);
}

BOOST_AUTO_TEST_CASE(MinimumStepIsAtLeastOneSecond)
{
    const auto result = flowChangeTimeStepLimit(1e-9, 10 * day, 0.0,
                                                std::vector<double>{100.0, 0.0, 0.0, 0.0},
                                                std::vector<double>{200.0, 0.0, 0.0, 0.0});
    BOOST_CHECK_EQUAL(result.time_step_limit, 1.0);
}

BOOST_AUTO_TEST_CASE(LimitIsWholeSeconds)
{
    // 1000 s * 0.3 / 0.7 = 428.57 s
    const auto result = flowChangeTimeStepLimit(0.3, 1000.0, 1.0,
                                                std::vector<double>{100.0, 0.0, 0.0, 0.0},
                                                std::vector<double>{170.0, 0.0, 0.0, 0.0});
    BOOST_CHECK_EQUAL(result.time_step_limit, 429.0);
}

BOOST_AUTO_TEST_CASE(DefaultedFractionDoesNotLimit)
{
    const auto result = flowChangeTimeStepLimit(no_limit, 10 * day, 0.1 * day,
                                                std::vector<double>{100.0, 0.0, 0.0, 0.0},
                                                std::vector<double>{200.0, 0.0, 0.0, 0.0});
    BOOST_CHECK_EQUAL(result.time_step_limit, infinity);
}

BOOST_AUTO_TEST_CASE(RateStartingFromZeroGivesMinimumStep)
{
    // The fractional change of a flow that starts from zero is infinite
    const auto result = flowChangeTimeStepLimit(0.005, 10 * day, 0.1 * day,
                                                std::vector<double>{0.0, 0.0, 0.0, 0.0},
                                                std::vector<double>{100.0, 0.0, 0.0, 0.0});
    BOOST_CHECK_EQUAL(result.max_change, infinity);
    BOOST_CHECK_EQUAL(result.time_step_limit, 0.1 * day);
}

BOOST_AUTO_TEST_CASE(RateStartingFromZeroWithDefaultedFractionDoesNotLimit)
{
    // This used to give the minimum step: a defaulted item 4 was not recognised as such
    const auto result = flowChangeTimeStepLimit(no_limit, 10 * day, 0.1 * day,
                                                std::vector<double>{0.0, 0.0, 0.0, 0.0},
                                                std::vector<double>{100.0, 0.0, 0.0, 0.0});
    BOOST_CHECK_EQUAL(result.time_step_limit, infinity);
}

BOOST_AUTO_TEST_CASE(NoStepLengthDoesNotLimit)
{
    const auto result = flowChangeTimeStepLimit(0.005, 0.0, 0.1 * day,
                                                std::vector<double>{100.0, 0.0, 0.0, 0.0},
                                                std::vector<double>{200.0, 0.0, 0.0, 0.0});
    BOOST_CHECK_EQUAL(result.time_step_limit, infinity);
}
