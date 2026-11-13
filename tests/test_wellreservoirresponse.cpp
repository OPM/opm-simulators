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

#define BOOST_TEST_MODULE TestWellReservoirResponse

#include <boost/test/unit_test.hpp>

#include <opm/simulators/wells/WellReservoirResponse.hpp>

#include <vector>

namespace
{

using Response = Opm::WellReservoirResponse<double>;

// Two components.  For the stopped well their cell pressures are 200 and 210, which
// average to 202 with the weights b.
const std::vector<double> b {2.0, 0.5};

// a_c with the cell pressures lowered by drop.
std::vector<double>
ipr_a(const double drop)
{
    return {b[0] * (200.0 - drop), b[1] * (210.0 - drop)};
}

// Inflow at a bhp with the cell pressures lowered by drop.
std::vector<double>
inflow(const double drop, const double bhp)
{
    const auto a = ipr_a(drop);
    return {a[0] - b[0] * bhp, a[1] - b[1] * bhp};
}

const std::vector<double> stopped_rates {0.0, 0.0};

// At bhp 190 with the cells 3 lower, the rates are 14 and 8.5: a drawdown of
// 22.5 / 2.5 = 9.  The ratio is 3 / 9.
constexpr double flowing_drop = 3.0;
constexpr double flowing_bhp = 190.0;
constexpr double ratio = 1.0 / 3.0;

} // Anonymous namespace

BOOST_AUTO_TEST_CASE(RatioFromStoppedToFlowing)
{
    Response response;
    response.update(1, true, ipr_a(0.0), b, stopped_rates);
    BOOST_CHECK(!response.ratio().has_value());
    response.update(2, false, ipr_a(flowing_drop), b, inflow(flowing_drop, flowing_bhp));
    BOOST_REQUIRE(response.ratio().has_value());
    BOOST_CHECK_CLOSE(*response.ratio(), ratio, 1.0e-10);
    BOOST_CHECK_CLOSE(response.stoppedPressureRise(), flowing_drop, 1.0e-10);
}

BOOST_AUTO_TEST_CASE(RatioFromFlowingToStopped)
{
    Response response;
    response.update(1, false, ipr_a(flowing_drop), b, inflow(flowing_drop, flowing_bhp));
    response.update(2, true, ipr_a(0.0), b, stopped_rates);
    BOOST_REQUIRE(response.ratio().has_value());
    BOOST_CHECK_CLOSE(*response.ratio(), ratio, 1.0e-10);
    BOOST_CHECK_EQUAL(response.stoppedPressureRise(), 0.0);
}

BOOST_AUTO_TEST_CASE(NoRatioWithoutChangeOfStatus)
{
    Response flowing;
    flowing.update(1, false, ipr_a(flowing_drop), b, inflow(flowing_drop, flowing_bhp));
    flowing.update(2, false, ipr_a(2.0 * flowing_drop), b, inflow(2.0 * flowing_drop, 180.0));
    BOOST_CHECK(!flowing.ratio().has_value());

    Response stopped;
    stopped.update(1, true, ipr_a(0.0), b, stopped_rates);
    stopped.update(2, true, ipr_a(-1.0), b, stopped_rates);
    BOOST_CHECK(!stopped.ratio().has_value());
}

BOOST_AUTO_TEST_CASE(NoRatioFromFirstIteration)
{
    // The cell pressures of iteration 0 are those of the previous time step.
    Response response;
    response.update(0, true, ipr_a(0.0), b, stopped_rates);
    response.update(1, false, ipr_a(flowing_drop), b, inflow(flowing_drop, flowing_bhp));
    BOOST_CHECK(!response.ratio().has_value());
}

BOOST_AUTO_TEST_CASE(NoRatioFromNonConsecutiveIterations)
{
    Response response;
    response.update(1, true, ipr_a(0.0), b, stopped_rates);
    response.update(3, false, ipr_a(flowing_drop), b, inflow(flowing_drop, flowing_bhp));
    BOOST_CHECK(!response.ratio().has_value());
}

BOOST_AUTO_TEST_CASE(NoRatioIfStoppedCellsAreLower)
{
    Response response;
    response.update(1, true, ipr_a(flowing_drop), b, stopped_rates);
    response.update(2, false, ipr_a(0.0), b, inflow(0.0, flowing_bhp));
    BOOST_CHECK(!response.ratio().has_value());
}

BOOST_AUTO_TEST_CASE(NoRatioWithoutInflow)
{
    Response response;
    response.update(1, true, ipr_a(0.0), b, stopped_rates);
    response.update(2, false, ipr_a(flowing_drop), b, inflow(flowing_drop, flowing_bhp));
    response.update(3, false, {0.0, 0.0}, {0.0, 0.0}, stopped_rates);
    BOOST_CHECK(!response.ratio().has_value());
}

BOOST_AUTO_TEST_CASE(FirstChangeOfStatusOnly)
{
    Response response;
    response.update(1, true, ipr_a(0.0), b, stopped_rates);
    response.update(2, false, ipr_a(flowing_drop), b, inflow(flowing_drop, flowing_bhp));
    response.update(3, true, ipr_a(1.0), b, stopped_rates);
    BOOST_REQUIRE(response.ratio().has_value());
    BOOST_CHECK_CLOSE(*response.ratio(), ratio, 1.0e-10);
}

BOOST_AUTO_TEST_CASE(NoDropAtMeasuredState)
{
    Response response;
    response.update(1, true, ipr_a(0.0), b, stopped_rates);
    response.update(2, false, ipr_a(flowing_drop), b, inflow(flowing_drop, flowing_bhp));
    BOOST_CHECK_SMALL(response.pressureDrop(flowing_bhp), 1.0e-10);
}

BOOST_AUTO_TEST_CASE(CellPressuresIndependentOfLastUpdate)
{
    // The cells are 0.25 * (202 - 195) = 1.75 lower than for the stopped well.
    const double bhp = 195.0;

    Response last_stopped;
    last_stopped.update(1, false, ipr_a(flowing_drop), b, inflow(flowing_drop, flowing_bhp));
    last_stopped.update(2, true, ipr_a(0.0), b, stopped_rates);
    BOOST_CHECK_CLOSE(last_stopped.pressureDrop(bhp), 1.75, 1.0e-10);

    Response last_flowing;
    last_flowing.update(1, true, ipr_a(0.0), b, stopped_rates);
    last_flowing.update(2, false, ipr_a(flowing_drop), b, inflow(flowing_drop, flowing_bhp));
    BOOST_CHECK_CLOSE(last_flowing.pressureDrop(bhp), 1.75 - flowing_drop, 1.0e-10);
}

BOOST_AUTO_TEST_CASE(DropFollowsDrawdown)
{
    Response response;
    response.update(1, true, ipr_a(0.0), b, stopped_rates);
    response.update(2, false, ipr_a(flowing_drop), b, inflow(flowing_drop, flowing_bhp));
    for (const double bhp : {185.0, 195.0, 201.0}) {
        // The inflow with the cells lowered by the drop is the one at a bhp higher by it.
        const auto q = inflow(flowing_drop, bhp + response.pressureDrop(bhp));
        const double drawdown = (q[0] + q[1]) / (b[0] + b[1]);
        BOOST_CHECK_CLOSE(flowing_drop + response.pressureDrop(bhp), ratio * drawdown, 1.0e-8);
    }
}

BOOST_AUTO_TEST_CASE(NoInflowAtStoppedPressure)
{
    Response response;
    response.update(1, true, ipr_a(0.0), b, stopped_rates);
    response.update(2, false, ipr_a(flowing_drop), b, inflow(flowing_drop, flowing_bhp));
    const double bhp = 202.0;
    const auto q = inflow(flowing_drop, bhp + response.pressureDrop(bhp));
    BOOST_CHECK_SMALL(q[0] + q[1], 1.0e-10);
}
