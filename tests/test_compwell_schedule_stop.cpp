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

#define BOOST_TEST_MODULE CompWellScheduleStop
#include <boost/test/unit_test.hpp>

#include <opm/simulators/flow/FlowGenericProblem_impl.hpp>

#include <flowexperimental/comp/flow_comp.hpp>

#include "RunDeckHelper.hpp"

#include <opm/output/eclipse/WStat.hpp>

#include <cstddef>
#include <string>

namespace {

const std::string props = R"(
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
)";

// Two layers without flow between them, at 100 and 60 bar, and a producer
// connected to both that the schedule stops from the start.
std::string layersDeck(const std::string& crossflow)
{
    return R"(RUNSPEC
METRIC
TABDIMS
/
OIL
GAS
DIMENS
1 1 2 /
COMPS
3 /
FULLIMP
WELLDIMS
1 2 1 1 /
START
1 'JAN' 2016 /
UNIFOUT

GRID
DX
2*1000 /
DY
2*1 /
DZ
2*1 /
TOPS
0 /
PERMX
2*100 /
PERMY
2*100 /
PERMZ
2*0 /
PORO
2*0.1 /
)" + props + R"(
SOLUTION
PRESSURE
100 60 /
SGAS
2*1 /
TEMPI
2*150 /
ZMF
2*0.1
2*0.3
2*0.6
/

SUMMARY
WOPR
/
WGPR
/
WSTAT
/
BPR
1 1 1 /
1 1 2 /
/

SCHEDULE
WELSPECS
PROD FIELD 1 1 1* GAS 3* )" + crossflow + R"( /
/
COMPDAT
PROD 1 1 1 2 OPEN 2* 0.0151 /
/
WCONPROD
PROD STOP BHP 2* 100000 2* 70 /
/
TSTEP
5*0.1 /
END
)";
}

// A 30-cell row with a gas injector at 150 bar and a producer at 50 bar.
// The schedule stops both wells after two report steps and opens them again
// two report steps later.
const std::string row_deck = R"(RUNSPEC
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
2 1 1 2 /
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
)" + props + R"(
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
2*0.1 /
WELOPEN
INJ STOP /
PROD STOP /
/
TSTEP
2*0.1 /
WELOPEN
INJ OPEN /
PROD OPEN /
/
TSTEP
2*0.1 /
END
)";

constexpr auto wstat_prod = static_cast<float>(Opm::WStat::numeric::PROD);
constexpr auto wstat_inj = static_cast<float>(Opm::WStat::numeric::INJ);
constexpr auto wstat_stop = static_cast<float>(Opm::WStat::numeric::STOP);

// The decks have no water phase.
using TypeTag = Opm::Properties::TTag::FlowCompProblem<3, false>;

} // namespace

using RunDeckFixture = Opm::RunDeckFixture;
BOOST_GLOBAL_FIXTURE(RunDeckFixture);

// A stopped producer has no surface rate. With crossflow, the layers exchange
// fluid through its wellbore; without, it has no flow at all.
BOOST_AUTO_TEST_CASE(StoppedProducerCrossflowsOnlyWhenAllowed)
{
    const Opm::TemporaryDirectory dir("compwell_schedule_stop_layers");

    const auto without = Opm::runDeck<TypeTag>(dir.path(), "WITHOUT_CROSSFLOW", layersDeck("NO"));
    const auto with = Opm::runDeck<TypeTag>(dir.path(), "WITH_CROSSFLOW", layersDeck("YES"));
    for (const auto* smry : {&without, &with}) {
        const auto& time = smry->get("TIME");
        BOOST_CHECK_GT(time.size(), 0u);
        for (std::size_t step = 0; step < time.size(); ++step) {
            BOOST_TEST_CONTEXT((smry == &with ? "with" : "without")
                               << " crossflow at day " << time[step]) {
                BOOST_CHECK_SMALL(smry->get("WOPR:PROD")[step], 1.0e-6f);
                BOOST_CHECK_SMALL(smry->get("WGPR:PROD")[step], 1.0e-6f);
                BOOST_CHECK_EQUAL(smry->get("WSTAT:PROD")[step], wstat_stop);
            }
        }
    }

    for (std::size_t step = 0; step < without.get("TIME").size(); ++step) {
        BOOST_TEST_CONTEXT("without crossflow at day " << without.get("TIME")[step]) {
            BOOST_CHECK_SMALL(without.get("BPR:1,1,1")[step] - 100.0f, 1.0e-3f);
            BOOST_CHECK_SMALL(without.get("BPR:1,1,2")[step] - 60.0f, 1.0e-3f);
        }
    }
    BOOST_CHECK_LT(with.get("BPR:1,1,1").back(), 99.0f);
    BOOST_CHECK_GT(with.get("BPR:1,1,2").back(), 61.0f);
}

// WELOPEN stops two flowing wells and opens them again. While stopped, they
// have no surface rate, and their BHP follows their cell pressure.
BOOST_AUTO_TEST_CASE(ScheduleStopsAndOpensWells)
{
    const Opm::TemporaryDirectory dir("compwell_schedule_stop_row");

    const auto smry = Opm::runDeck<TypeTag>(dir.path(), "ROW", row_deck);
    const auto& time = smry.get("TIME");
    BOOST_REQUIRE_EQUAL(time.size(), 6u);
    for (std::size_t step = 0; step < time.size(); ++step) {
        BOOST_TEST_CONTEXT("day " << time[step]) {
            if (step == 2 || step == 3) {
                BOOST_CHECK_SMALL(smry.get("WGIR:INJ")[step], 1.0e-6f);
                BOOST_CHECK_SMALL(smry.get("WOPR:PROD")[step], 1.0e-6f);
                BOOST_CHECK_SMALL(smry.get("WGPR:PROD")[step], 1.0e-6f);
                BOOST_CHECK_EQUAL(smry.get("WSTAT:INJ")[step], wstat_stop);
                BOOST_CHECK_EQUAL(smry.get("WSTAT:PROD")[step], wstat_stop);
                BOOST_CHECK_CLOSE(smry.get("WBHP:INJ")[step], smry.get("BPR:1,1,1")[step], 0.05f);
                BOOST_CHECK_CLOSE(smry.get("WBHP:PROD")[step], smry.get("BPR:30,1,1")[step], 0.05f);
            } else {
                BOOST_CHECK_GT(smry.get("WGIR:INJ")[step], 1.0f);
                BOOST_CHECK_GT(smry.get("WGPR:PROD")[step], 1.0f);
                BOOST_CHECK_EQUAL(smry.get("WSTAT:INJ")[step], wstat_inj);
                BOOST_CHECK_EQUAL(smry.get("WSTAT:PROD")[step], wstat_prod);
            }
        }
    }
}
