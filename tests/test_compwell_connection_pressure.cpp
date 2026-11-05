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

#include <cstddef>
#include <string>

namespace {

// Three 10 m layers without flow between them: gas, oil and gas again. A
// producer at its BHP limit is connected to all three, with its reference
// depth at the top connection.
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
WCONPROD
PROD OPEN BHP 2* 100000 2* 70 /
/
TSTEP
5*0.1 /
END
)";

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
    const auto smry = Opm::runDeck<TypeTag>(dir.path(), "HYDROSTATICS", deck);

    constexpr float gravity = 9.80665f;
    constexpr float layer_thickness = 10.0f;
    // the density [kg/m^3] of the column between two connections
    const auto density = [&](const float upper, const float lower)
    { return (lower - upper) * 1.0e5f / (gravity * layer_thickness); };

    const auto& time = smry.get("TIME");
    const auto& bhp = smry.get("WBHP:PROD");
    const auto& top = smry.get("CPR:PROD:1,1,1");
    const auto& middle = smry.get("CPR:PROD:1,1,2");
    const auto& bottom = smry.get("CPR:PROD:1,1,3");
    BOOST_REQUIRE_GT(time.size(), 0u);
    for (std::size_t step = 0; step < time.size(); ++step) {
        BOOST_TEST_CONTEXT("producer at day " << time[step]) {
            BOOST_CHECK_SMALL(top[step] - bhp[step], 1.0e-4f);
            const float gas_column = density(middle[step], bottom[step]);
            const float mixed_column = density(top[step], middle[step]);
            BOOST_CHECK_GT(gas_column, 50.0f);
            BOOST_CHECK_LT(gas_column, 100.0f);
            BOOST_CHECK_GT(mixed_column, 1.15f * gas_column);
        }
    }
}
