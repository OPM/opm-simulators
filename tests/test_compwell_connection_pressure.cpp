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

#define BOOST_TEST_MODULE CompWellConnectionPressure
#include <boost/test/unit_test.hpp>

#include <opm/simulators/flow/FlowGenericProblem_impl.hpp>

#include <flowexperimental/comp/flow_comp.hpp>

#include "RunDeckHelper.hpp"

#include <algorithm>
#include <cstddef>
#include <initializer_list>
#include <iterator>
#include <string>
#include <vector>

namespace {

// Three 10 m layers without flow between them: gas, oil and gas again. A
// producer is connected to all three, with its reference depth at the top
// connection.
const std::string deck = R"(RUNSPEC
METRIC
TABDIMS
/
OIL
GAS
DIMENS
1 1 3 /
COMPS
3 /
FULLIMP
WELLDIMS
1 3 1 1 /
START
1 'JAN' 2016 /
UNIFOUT

GRID
DX
3*1000 /
DY
3*1 /
DZ
3*10 /
TOPS
0 /
PERMX
3*100 /
PERMY
3*100 /
PERMZ
3*0 /
PORO
3*0.1 /

PROPS
CNAMES
CO2
METHANE
DECANE
/
ROCK
68.9476 0 /
EOS
PR /
BIC
0
0
0
/
ACF
0.22394 0.01142 0.4884 /
PCRIT
73.773 45.992 21.03 /
TCRIT
304.128 190.564 617.7 /
MW
44.00 16.04 142.28 /
VCRIT
0.09412 0.09863 0.60980 /
SGOF
0.0 0.00 1.00 0.0
0.2 0.04 0.64 0.0
0.4 0.16 0.36 0.0
0.6 0.36 0.16 0.0
0.8 0.64 0.04 0.0
1.0 1.00 0.00 0.0
/
STCOND
15.0 1.0 /

SOLUTION
PRESSURE
3*100 /
SGAS
1 0 1 /
TEMPI
3*150 /
ZMF
0.5 0.05 0.5
0.5 0.05 0.5
0.0 0.90 0.0
/

SUMMARY
WBHP
/
CPR
'PROD' /
/

SCHEDULE
WELSPECS
PROD FIELD 1 1 1* GAS /
/
COMPDAT
PROD 1 1 1 3 OPEN 2* 0.0151 /
/
)";

// The producer at its BHP limit.
const std::string bhp_limit = R"(WCONPROD
PROD OPEN BHP 2* 100000 2* 70 /
/
TSTEP
5*0.1 /
END
)";

// The producer at 70 bar.
const std::string bhp_fixed = R"(WCONPROD
PROD OPEN BHP 5* 70 /
/
TSTEP
4*0.1 /
END
)";

// The BHP drops from 70 to 60 bar after three time steps.
const std::string bhp_drop = R"(WCONPROD
PROD OPEN BHP 5* 70 /
/
TSTEP
3*0.1 /
WCONPROD
PROD OPEN BHP 5* 60 /
/
TSTEP
8*0.1 /
END
)";

// The reference depth moves from the top of the reservoir to the top
// connection after three time steps.
const std::string ref_depth_change = R"(WELSPECS
PROD FIELD 1 1 0 GAS /
/
WCONPROD
PROD OPEN BHP 5* 70 /
/
TSTEP
3*0.1 /
WELSPECS
PROD FIELD 1 1 5 GAS /
/
TSTEP
3*0.1 /
END
)";

// The producer is connected to the top two layers, and after three time steps
// to the top and bottom ones.
const std::string connection_move = R"(COMPDAT
PROD 1 1 3 3 SHUT 2* 0.0151 /
/
WCONPROD
PROD OPEN BHP 5* 70 /
/
TSTEP
3*0.1 /
COMPDAT
PROD 1 1 2 2 SHUT 2* 0.0151 /
PROD 1 1 3 3 OPEN 2* 0.0151 /
/
TSTEP
3*0.1 /
END
)";

// the density [kg/m^3] of the wellbore column between two connections
float columnDensity(const float upper, const float lower)
{
    constexpr float gravity = 9.80665f;
    constexpr float layer_thickness = 10.0f;
    return (lower - upper) * 1.0e5f / (gravity * layer_thickness);
}

// the parameter that damps the density of each wellbore section instead of the
// offset of each connection
std::string dampSections(const bool sections)
{
    return std::string("--comp-well-damp-section-densities=") + (sections ? "true" : "false");
}

// The deck has no water phase.
using TypeTag = Opm::Properties::TTag::FlowCompProblem<3, false>;

} // namespace

using RunDeckFixture = Opm::RunDeckFixture;
BOOST_GLOBAL_FIXTURE(RunDeckFixture);

// The wellbore between two connections holds what flows up past the lower
// one. Below the oil layer that is the gas of the bottom layer alone, and
// above it that gas together with the oil, so the column there is heavier.
BOOST_AUTO_TEST_CASE(WellboreColumnHoldsTheFluidFlowingPastIt)
{
    const Opm::TemporaryDirectory dir("compwell_connection_pressure");
    const auto smry = Opm::runDeck<TypeTag>(dir.path(), "HYDROSTATICS", deck + bhp_limit);

    const auto& time = smry.get("TIME");
    const auto& bhp = smry.get("WBHP:PROD");
    const auto& top = smry.get("CPR:PROD:1,1,1");
    const auto& middle = smry.get("CPR:PROD:1,1,2");
    const auto& bottom = smry.get("CPR:PROD:1,1,3");
    BOOST_REQUIRE_GT(time.size(), 0u);
    for (std::size_t step = 0; step < time.size(); ++step) {
        BOOST_TEST_CONTEXT("producer at day " << time[step]) {
            BOOST_CHECK_SMALL(top[step] - bhp[step], 1.0e-4f);
            const float gas_column = columnDensity(middle[step], bottom[step]);
            const float mixed_column = columnDensity(top[step], middle[step]);
            BOOST_CHECK_GT(gas_column, 50.0f);
            BOOST_CHECK_LT(gas_column, 100.0f);
            BOOST_CHECK_GT(mixed_column, 1.15f * gas_column);
        }
    }
}

// By default, the offset of each connection is damped: in each time step it
// moves a third of the way from its last value to the new offset of the
// connection above it plus the weight of the column between them. While the
// gas below the oil keeps its density, the offset of the bottom connection
// lags behind the falling offset of the middle one.
BOOST_AUTO_TEST_CASE(OffsetApproachesTheOneAboveGradually)
{
    const Opm::TemporaryDirectory dir("compwell_connection_pressure_offset_damping");
    const auto smry = Opm::runDeck<TypeTag>(dir.path(), "OFFSETDAMPING", deck + bhp_fixed);

    const auto& top = smry.get("CPR:PROD:1,1,1");
    const auto& middle = smry.get("CPR:PROD:1,1,2");
    const auto& bottom = smry.get("CPR:PROD:1,1,3");
    BOOST_REQUIRE_EQUAL(top.size(), 4u);
    // the first time step is not damped, so it gives the weight of the gas column
    const float gas_column = bottom[0] - middle[0];
    for (std::size_t step = 1; step < top.size(); ++step) {
        BOOST_TEST_CONTEXT("time step " << step) {
            const float middle_offset = middle[step] - top[step];
            const float last_offset = bottom[step - 1] - top[step - 1];
            const float target = middle_offset + gas_column;
            BOOST_CHECK_CLOSE(bottom[step] - top[step],
                              2.0f / 3.0f * last_offset + 1.0f / 3.0f * target, 0.1f);
        }
    }
}

// With the section densities damped, the density of each column moves a third
// of the way from its last value to the one of the flow past it in each time
// step. After the BHP drops, the gas below the oil gets lighter, and each change
// of its column density is two thirds of the one before.
BOOST_AUTO_TEST_CASE(ColumnDensityApproachesItsNewValueGradually)
{
    const Opm::TemporaryDirectory dir("compwell_connection_pressure_damping");
    const auto smry = Opm::runDeck<TypeTag>(dir.path(), "DAMPING", deck + bhp_drop,
                                            {dampSections(true)});

    const auto& middle = smry.get("CPR:PROD:1,1,2");
    const auto& bottom = smry.get("CPR:PROD:1,1,3");
    std::vector<float> changes;
    for (std::size_t step = 1; step < middle.size(); ++step) {
        changes.push_back(columnDensity(middle[step], bottom[step])
                          - columnDensity(middle[step - 1], bottom[step - 1]));
    }
    // The densities are taken at the start of a time step, so they first see
    // the lower BHP in the second time step at 60 bar.
    const auto first = std::ranges::find_if(changes, [](const float change)
                                            { return change < -1.0f; });
    BOOST_REQUIRE_GE(std::distance(first, changes.end()), 4);
    for (auto it = first + 1; it != first + 4; ++it) {
        BOOST_TEST_CONTEXT("time step " << std::distance(first, it) << " after the first change") {
            BOOST_CHECK_CLOSE(*it / *(it - 1), 2.0f / 3.0f, 5.0f);
        }
    }
}

// The damped densities of the wellbore do not depend on the reference depth,
// and damped offsets from the old one count no more, so either way the
// connection at the new reference depth is at the BHP from the first time step
// on.
BOOST_AUTO_TEST_CASE(OffsetsFollowNewReferenceDepth)
{
    const Opm::TemporaryDirectory dir("compwell_connection_pressure_ref_depth");
    for (const bool sections : {false, true}) {
        BOOST_TEST_CONTEXT("section densities damped: " << sections) {
            const auto smry = Opm::runDeck<TypeTag>(dir.path(), "REFDEPTH",
                                                    deck + ref_depth_change,
                                                    {dampSections(sections)});

            const auto& time = smry.get("TIME");
            const auto& bhp = smry.get("WBHP:PROD");
            const auto& top = smry.get("CPR:PROD:1,1,1");
            std::size_t moved = 0;
            for (std::size_t step = 0; step < time.size(); ++step) {
                BOOST_TEST_CONTEXT("producer at day " << time[step]) {
                    if (time[step] < 0.35f) {
                        // 5 m of wellbore fluid above the top connection
                        BOOST_CHECK_GT(top[step] - bhp[step], 0.01f);
                    } else {
                        BOOST_CHECK_SMALL(top[step] - bhp[step], 1.0e-4f);
                        ++moved;
                    }
                }
            }
            BOOST_CHECK_EQUAL(moved, 3u);
        }
    }
}

// The densities and offsets of the last time step count only for the
// connections they were taken at. After the lower connection moves from the oil
// layer to the gas layer below it, the well has as many connections as before,
// but the column above that connection holds the gas from the first time step
// on, not a blend with the oil of the old section.
BOOST_AUTO_TEST_CASE(DensitiesRestartWhenConnectionsMove)
{
    const Opm::TemporaryDirectory dir("compwell_connection_pressure_connection_move");
    for (const bool sections : {false, true}) {
        BOOST_TEST_CONTEXT("section densities damped: " << sections) {
            const auto smry = Opm::runDeck<TypeTag>(dir.path(), "CONNECTIONMOVE",
                                                    deck + connection_move,
                                                    {dampSections(sections)});

            const auto& time = smry.get("TIME");
            const auto& top = smry.get("CPR:PROD:1,1,1");
            const auto& middle = smry.get("CPR:PROD:1,1,2");
            const auto& bottom = smry.get("CPR:PROD:1,1,3");
            std::size_t moved = 0;
            for (std::size_t step = 0; step < time.size(); ++step) {
                BOOST_TEST_CONTEXT("producer at day " << time[step]) {
                    if (time[step] < 0.35f) {
                        BOOST_CHECK_GT(columnDensity(top[step], middle[step]), 200.0f);
                    } else {
                        // the connections are two layers apart
                        const float gas_column
                            = columnDensity(top[step], bottom[step]) / 2.0f;
                        BOOST_CHECK_GT(gas_column, 50.0f);
                        BOOST_CHECK_LT(gas_column, 100.0f);
                        ++moved;
                    }
                }
            }
            BOOST_CHECK_EQUAL(moved, 3u);
        }
    }
}
