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

#define BOOST_TEST_MODULE CompWellBhpLimit
#include <boost/test/unit_test.hpp>

#include <opm/simulators/flow/FlowGenericProblem_impl.hpp>

#include <flowexperimental/comp/flow_comp.hpp>

#include "RunDeckHelper.hpp"

#include <opm/output/eclipse/WStat.hpp>

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

namespace {

// A 30-cell row with a gas injector at 150 bar and a producer at 50 bar.
// After four report steps, the appended schedule puts one of the wells
// beyond its BHP limit.
std::string deck(const std::string& schedule)
{
    return R"(RUNSPEC
METRIC
TABDIMS
/
OIL
GAS
DIMENS
30 1 1 /
COMPS
3 /
FULLIMP
WELLDIMS
3 1 1 3 /
START
1 'JAN' 2016 /
UNIFOUT

GRID
DX
30*10 /
DY
30*1 /
DZ
30*1 /
TOPS
30*0 /
PERMX
30*100 /
PERMY
30*100 /
PERMZ
30*100 /
PORO
30*0.1 /

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
30*75 /
SGAS
30*1 /
TEMPI
30*150 /
ZMF
30*0.1
30*0.3
30*0.6
/

SUMMARY
WBHP
/
WSTAT
/
WOPR
/
WGPR
/
WOIR
/
WGIR
/
BPR
1 1 1 /
30 1 1 /
/

SCHEDULE
WELSPECS
INJ  FIELD  1 1 1* GAS /
PROD FIELD 30 1 1* GAS /
/
COMPDAT
INJ   1 1 1 1 OPEN 2* 0.0151 /
PROD 30 1 1 1 OPEN 2* 0.0151 /
/
WELLSTRE
ISTR 1.0 0.0 0.0 /
/
WINJGAS
INJ STREAM ISTR /
/
WCONINJE
INJ GAS OPEN BHP 100000 1* 150 /
/
WCONPROD
PROD OPEN BHP 2* 100000 2* 50 /
/
TSTEP
0.01 0.1 0.4 1 /
)" + schedule + "END\n";
}

// Index of the first report step after the appended schedule takes effect.
constexpr std::size_t first_step = 4;

const std::vector<std::string> rates = {"WOPR", "WGPR", "WOIR", "WGIR"};

constexpr auto wstat_prod = static_cast<float>(Opm::WStat::numeric::PROD);
constexpr auto wstat_stop = static_cast<float>(Opm::WStat::numeric::STOP);

// The deck has no water phase.
using TypeTag = Opm::Properties::TTag::FlowCompProblem<3, false>;

} // namespace

using RunDeckFixture = Opm::RunDeckFixture;
BOOST_GLOBAL_FIXTURE(RunDeckFixture);

// Raise the producer's minimum BHP above its cell pressure. It should stop
// rather than inject, then reopen when the injector raises the pressure enough.
BOOST_AUTO_TEST_CASE(ProducerStopsAndReopensAtItsBhpLimit)
{
    // Give each run a separate output directory.
    const Opm::TemporaryDirectory dir("compwell_bhp_limit_producer");

    const auto smry = Opm::runDeck<TypeTag>(dir.path(), "PRODUCER", deck(R"(WCONPROD
PROD OPEN BHP 2* 100000 2* 120 /
/
TSTEP
8*0.25 4*1 /
)"));
    const auto& time = smry.get("TIME");
    const auto& bhp = smry.get("WBHP:PROD");
    const auto& cell_pressure = smry.get("BPR:30,1,1");

    // Find the first report step at which the producer reaches its BHP limit again.
    auto reopened = first_step;
    while (reopened < bhp.size() && std::abs(bhp[reopened] - 120.0f) > 1.0e-3f) {
        ++reopened;
    }
    BOOST_CHECK_GT(reopened, first_step);
    BOOST_CHECK_LT(reopened, bhp.size());

    for (auto step = first_step; step < reopened; ++step) {
        BOOST_TEST_CONTEXT("stopped producer at day " << time[step]) {
            for (const auto& rate : rates) {
                BOOST_CHECK_SMALL(smry.get(rate + ":PROD")[step], 1.0e-6f);
            }
            BOOST_CHECK_CLOSE(bhp[step], cell_pressure[step], 0.05f);
            BOOST_CHECK_GT(bhp[step], bhp[step - 1]);
            BOOST_CHECK_EQUAL(smry.get("WSTAT:PROD")[step], wstat_stop);
        }
    }
    for (auto step = reopened; step < bhp.size(); ++step) {
        BOOST_TEST_CONTEXT("reopened producer at day " << time[step]) {
            // It reopens as its cell pressure passes the limit, so by the end of
            // that report step it may produce little.
            const bool first = step == reopened;
            BOOST_CHECK_GT(smry.get("WOPR:PROD")[step], first ? 0.0f : 0.1f);
            BOOST_CHECK_GT(smry.get("WGPR:PROD")[step], first ? 0.0f : 1.0f);
            BOOST_CHECK_SMALL(smry.get("WOIR:PROD")[step], 1.0e-6f);
            BOOST_CHECK_SMALL(smry.get("WGIR:PROD")[step], 1.0e-6f);
            BOOST_CHECK_EQUAL(smry.get("WSTAT:PROD")[step], wstat_prod);
        }
    }
}

// Lower the injector's maximum BHP below its cell pressure. It should stop
// rather than produce. Its BHP falls with the cell pressure as the producer
// depletes the row, but the injector remains stopped.
BOOST_AUTO_TEST_CASE(InjectorStopsAtItsBhpLimit)
{
    // Give each run a separate output directory.
    const Opm::TemporaryDirectory dir("compwell_bhp_limit_injector");

    const auto smry = Opm::runDeck<TypeTag>(dir.path(), "INJECTOR", deck(R"(WCONINJE
INJ GAS OPEN BHP 100000 1* 40 /
/
TSTEP
8*0.25 4*1 /
)"));
    const auto& time = smry.get("TIME");
    const auto& bhp = smry.get("WBHP:INJ");
    const auto& cell_pressure = smry.get("BPR:1,1,1");

    BOOST_CHECK_GT(time.size(), first_step);
    for (auto step = first_step; step < time.size(); ++step) {
        BOOST_TEST_CONTEXT("stopped injector at day " << time[step]) {
            for (const auto& rate : rates) {
                BOOST_CHECK_SMALL(smry.get(rate + ":INJ")[step], 1.0e-6f);
            }
            BOOST_CHECK_CLOSE(bhp[step], cell_pressure[step], 0.05f);
            BOOST_CHECK_LT(bhp[step], bhp[step - 1]);
            BOOST_CHECK_EQUAL(smry.get("WSTAT:INJ")[step], wstat_stop);
            BOOST_CHECK_GT(smry.get("WGPR:PROD")[step], 1.0f);
        }
    }
}

// A second producer next to the first one, at a lower BHP limit, pulls the
// first one's cell below its limit during a time step. The first producer
// should stop in that time step rather than inject.
BOOST_AUTO_TEST_CASE(ProducerStopsInTheTimeStepItsCellFallsBelowItsLimit)
{
    const Opm::TemporaryDirectory dir("compwell_bhp_limit_neighbour");

    const auto smry = Opm::runDeck<TypeTag>(dir.path(), "NEIGHBOUR", deck(R"(WELSPECS
PROD2 FIELD 29 1 1* GAS /
/
COMPDAT
PROD2 29 1 1 1 OPEN 2* 0.0151 /
/
WCONPROD
PROD2 OPEN BHP 2* 100000 2* 10 /
/
TSTEP
10*0.1 /
)"));
    const auto& time = smry.get("TIME");
    const auto& cell_pressure = smry.get("BPR:30,1,1");

    // The producer flows until its neighbour opens, and its cell then falls
    // below its limit.
    BOOST_REQUIRE_GT(time.size(), first_step);
    BOOST_CHECK_GT(smry.get("WGPR:PROD")[first_step - 1], 1.0f);
    BOOST_CHECK_GT(cell_pressure[first_step - 1], 50.0f);
    BOOST_CHECK_LT(cell_pressure.back(), 50.0f);

    for (auto step = first_step; step < time.size(); ++step) {
        BOOST_TEST_CONTEXT("producer at day " << time[step]) {
            BOOST_CHECK_SMALL(smry.get("WOIR:PROD")[step], 1.0e-6f);
            BOOST_CHECK_SMALL(smry.get("WGIR:PROD")[step], 1.0e-6f);
            if (cell_pressure[step] < 50.0f) {
                BOOST_CHECK_SMALL(smry.get("WOPR:PROD")[step], 1.0e-6f);
                BOOST_CHECK_SMALL(smry.get("WGPR:PROD")[step], 1.0e-6f);
                BOOST_CHECK_EQUAL(smry.get("WSTAT:PROD")[step], wstat_stop);
            }
        }
    }
}
