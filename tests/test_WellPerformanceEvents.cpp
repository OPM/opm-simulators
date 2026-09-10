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

#define BOOST_TEST_MODULE TestWellPerformanceEvents

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <boost/test/unit_test.hpp>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <opm/simulators/wells/WellPerformanceEvents.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Python/Python.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Action/Actions.hpp>
#include <opm/input/eclipse/Schedule/Action/ActionX.hpp>
#include <opm/input/eclipse/Schedule/Well/WellEnums.hpp>

#include <opm/output/data/Wells.hpp>

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

// Four wells with four connections each.  P1 closes connections one at a
// time (CON), P2 closes to the bottom (+CON), P3 inherits +CON from WECON,
// I1 injects.  At the second report step I1 becomes a producer and P3 an
// injector.
std::string deckString()
{
    return R"(RUNSPEC
TITLE
 WPWE unit test
DIMENS
 4 1 4 /
OIL
WATER
GAS
METRIC
START
 1 'JAN' 2020 /
WELLDIMS
 4 4 2 4 /
GRID
DXV
 4*100.0 /
DYV
 100.0 /
DZV
 4*10.0 /
TOPS
 4*1000.0 /
PORO
 16*0.3 /
PERMX
 16*100.0 /
PERMY
 16*100.0 /
PERMZ
 16*10.0 /
SCHEDULE
WELSPECS
 'P1' 'G1' 1 1 1* 'OIL' /
 'P2' 'G1' 2 1 1* 'OIL' /
 'P3' 'G1' 3 1 1* 'OIL' /
 'I1' 'G1' 4 1 1* 'WATER' /
/
COMPDAT
 'P1' 1 1 1 4 'OPEN' 1* 1* 0.2 /
 'P2' 2 1 1 4 'OPEN' 1* 1* 0.2 /
 'P3' 3 1 1 4 'OPEN' 1* 1* 0.2 /
 'I1' 4 1 1 4 'OPEN' 1* 1* 0.2 /
/
WCONPROD
 'P1' 'OPEN' 'ORAT' 100.0 4* 50.0 /
 'P2' 'OPEN' 'ORAT' 100.0 4* 50.0 /
 'P3' 'OPEN' 'ORAT' 100.0 4* 50.0 /
/
WCONINJE
 'I1' 'WATER' 'OPEN' 'RATE' 100.0 1* 400.0 /
/
CECON
 'P1' 2* 1 4 0.6 2* 'CON' /
 'P2' 2* 1 4 0.6 2* '+CON' /
/
WECON
 'P3' 2* 0.9 2* '+CON' /
/
DATES
 1 'FEB' 2020 /
/
WCONPROD
 'I1' 'OPEN' 'ORAT' 100.0 4* 50.0 /
/
WCONINJE
 'P3' 'WATER' 'OPEN' 'RATE' 100.0 1* 400.0 /
/
DATES
 1 'MAR' 2020 /
/
END
)";
}

struct Setup
{
    Setup()
        : deck  { Opm::Parser{}.parseString(deckString()) }
        , es    { deck }
        , sched { deck, es, std::make_shared<Opm::Python>() }
    {}

    Opm::Deck deck;
    Opm::EclipseState es;
    Opm::Schedule sched;
};

using Entry = Opm::WellStatusSnapshot::Entry;

Opm::WellStatusSnapshot allOpen()
{
    auto snap = Opm::WellStatusSnapshot{};

    for (const auto* well : { "P1", "P2", "P3", "I1" }) {
        snap.wells[well] = Entry { Opm::WellStatus::OPEN, {1, 2, 3, 4} };
    }

    return snap;
}

bool noEvents(const Opm::data::WellEvents& events)
{
    return events == Opm::data::WellEvents{};
}

} // Anonymous namespace

BOOST_FIXTURE_TEST_CASE(NoChange, Setup)
{
    auto tracker = Opm::WellPerformanceEvents{};

    tracker.beginTimeStep(sched, 0, allOpen());
    tracker.accumulate(allOpen());

    for (const auto* well : { "P1", "P2", "P3", "I1" }) {
        BOOST_CHECK(noEvents(tracker.events(well)));
    }

    // Unknown wells report no events rather than failing.
    BOOST_CHECK(noEvents(tracker.events("NO-SUCH-WELL")));
}

BOOST_FIXTURE_TEST_CASE(ConnectionsClosed_CON, Setup)
{
    auto tracker = Opm::WellPerformanceEvents{};
    tracker.beginTimeStep(sched, 0, allOpen());

    auto now = allOpen();
    now.wells["P1"].openCompletions = {1, 2, 3};
    tracker.accumulate(now);

    {
        const auto& ev = tracker.events("P1");
        BOOST_CHECK_EQUAL(ev.connsClosed, 1);
        BOOST_CHECK_EQUAL(ev.closedToBottom, 0);
        BOOST_CHECK_EQUAL(ev.connsOpened, 0);
        BOOST_CHECK_EQUAL(ev.shut, 0);
    }

    // Events accumulate within the time step.
    now.wells["P1"].openCompletions = {1, 2};
    tracker.accumulate(now);

    BOOST_CHECK_EQUAL(tracker.events("P1").connsClosed, 2);
    BOOST_CHECK_EQUAL(tracker.events("P1").closedToBottom, 0);

    BOOST_CHECK(noEvents(tracker.events("P2")));
}

BOOST_FIXTURE_TEST_CASE(LastConnectionClosed_CON, Setup)
{
    // A CON workover closing the last open connection leaves the well
    // closed to the bottom (WPWE3) in addition to the connection count
    // (WPWE2), and the well is shut (WPWE7).
    auto tracker = Opm::WellPerformanceEvents{};

    auto before = allOpen();
    before.wells["P1"].openCompletions = {1};
    tracker.beginTimeStep(sched, 0, before);

    auto now = before;
    now.wells["P1"] = Entry { Opm::WellStatus::SHUT, {}, {}, true };
    tracker.accumulate(now);

    const auto& ev = tracker.events("P1");
    BOOST_CHECK_EQUAL(ev.connsClosed, 1);
    BOOST_CHECK_EQUAL(ev.closedToBottom, 1);
    BOOST_CHECK_EQUAL(ev.shut, 1);
    BOOST_CHECK_EQUAL(ev.stopped, 0);
}

BOOST_FIXTURE_TEST_CASE(LastConnectionClosedWithoutWorkover, Setup)
{
    // WPWE3's first route is restricted to CON workovers.  A well left with
    // no flowing connection by anything else -- the deck, a whole well
    // workover, a well test -- is not reported as closed to the bottom.
    auto tracker = Opm::WellPerformanceEvents{};

    auto before = allOpen();
    before.wells["P1"].openCompletions = {1};
    tracker.beginTimeStep(sched, 0, before);

    auto now = before;
    now.wells["P1"] = Entry { Opm::WellStatus::SHUT, {}, {}, false };
    tracker.accumulate(now);

    const auto& ev = tracker.events("P1");
    BOOST_CHECK_EQUAL(ev.connsClosed, 1);
    BOOST_CHECK_EQUAL(ev.closedToBottom, 0);
    BOOST_CHECK_EQUAL(ev.shut, 1);
}

BOOST_FIXTURE_TEST_CASE(ClosedToBottom_PlusCON, Setup)
{
    // A +CON workover is reported through WPWE3 alone.
    auto tracker = Opm::WellPerformanceEvents{};
    tracker.beginTimeStep(sched, 0, allOpen());

    auto now = allOpen();
    now.wells["P2"].openCompletions = {1, 2};
    now.wells["P2"].closedByConPlus = {3, 4};
    tracker.accumulate(now);

    const auto& ev = tracker.events("P2");
    BOOST_CHECK_EQUAL(ev.closedToBottom, 1);
    BOOST_CHECK_EQUAL(ev.connsClosed, 0);
    BOOST_CHECK_EQUAL(ev.shut, 0);
}

BOOST_FIXTURE_TEST_CASE(ClosedToBottom_WECON, Setup)
{
    // No CECON on P3: the well level WECON procedure (+CON) applies.
    auto tracker = Opm::WellPerformanceEvents{};
    tracker.beginTimeStep(sched, 0, allOpen());

    auto now = allOpen();
    now.wells["P3"].openCompletions = {1, 2, 3};
    now.wells["P3"].closedByConPlus = {4};
    tracker.accumulate(now);

    const auto& ev = tracker.events("P3");
    BOOST_CHECK_EQUAL(ev.closedToBottom, 1);
    BOOST_CHECK_EQUAL(ev.connsClosed, 0);
}

BOOST_FIXTURE_TEST_CASE(MixedCONAndPlusCON, Setup)
{
    auto tracker = Opm::WellPerformanceEvents{};
    tracker.beginTimeStep(sched, 0, allOpen());
    auto now = allOpen();
    now.wells["P2"].openCompletions = {1, 2};
    // Completion 3 closed by CON and 4 by +CON in the same update.
    now.wells["P2"].closedByConPlus = {4};
    tracker.accumulate(now);
    BOOST_CHECK_EQUAL(tracker.events("P2").connsClosed, 1);
    BOOST_CHECK_EQUAL(tracker.events("P2").closedToBottom, 1);
}

BOOST_FIXTURE_TEST_CASE(LumpedCompletionCountsConnections, Setup)
{
    // COMPLUMP maps several connections onto one completion number, so a
    // snapshot carries that number once per connection.  WPWE1 and WPWE2 count
    // connections, so closing or opening a lumped completion must count every
    // connection it carries rather than the completion itself.
    auto tracker = Opm::WellPerformanceEvents{};

    auto before = allOpen();
    before.wells["P1"].openCompletions = {1, 1, 2, 2};
    tracker.beginTimeStep(sched, 0, before);

    auto now = before;
    now.wells["P1"].openCompletions = {1, 1};
    tracker.accumulate(now);
    BOOST_CHECK_EQUAL(tracker.events("P1").connsClosed, 2);

    tracker.beginTimeStep(sched, 0, now);
    tracker.accumulate(before);
    BOOST_CHECK_EQUAL(tracker.events("P1").connsOpened, 2);
}

BOOST_FIXTURE_TEST_CASE(ActualWorkoverOverridesConfiguredLimits, Setup)
{
    auto tracker = Opm::WellPerformanceEvents{};
    tracker.beginTimeStep(sched, 0, allOpen());
    auto now = allOpen();
    // An inactive CECON +CON setting must not override an actual CON closure.
    now.wells["P2"].openCompletions = {1, 2, 3};
    // Conversely, WECON +CON can close connections with CECON CON configured.
    now.wells["P1"].openCompletions = {1, 2};
    now.wells["P1"].closedByConPlus = {3, 4};
    tracker.accumulate(now);
    BOOST_CHECK_EQUAL(tracker.events("P2").connsClosed, 1);
    BOOST_CHECK_EQUAL(tracker.events("P2").closedToBottom, 0);
    BOOST_CHECK_EQUAL(tracker.events("P1").connsClosed, 0);
    BOOST_CHECK_EQUAL(tracker.events("P1").closedToBottom, 1);

    tracker.commitTimeStep(sched, 0);
    tracker.beginTimeStep(sched, 0, now);
    now.wells["P1"].openCompletions = {1};
    tracker.accumulate(now);
    // Previously recorded +CON closures are not new events, and must not
    // affect the count for this timestep's CON closure.
    BOOST_CHECK_EQUAL(tracker.events("P1").connsClosed, 1);
    BOOST_CHECK_EQUAL(tracker.events("P1").closedToBottom, 0);
}

BOOST_FIXTURE_TEST_CASE(WellShutAndStopped, Setup)
{
    auto tracker = Opm::WellPerformanceEvents{};
    tracker.beginTimeStep(sched, 0, allOpen());

    auto now = allOpen();
    now.wells["P1"].status = Opm::WellStatus::SHUT;   // e.g., a WELL workover
    now.wells["P2"].status = Opm::WellStatus::STOP;
    tracker.accumulate(now);

    BOOST_CHECK_EQUAL(tracker.events("P1").shut, 1);
    BOOST_CHECK_EQUAL(tracker.events("P1").stopped, 0);
    BOOST_CHECK_EQUAL(tracker.events("P1").connsClosed, 0);
    BOOST_CHECK_EQUAL(tracker.events("P1").closedToBottom, 0);

    BOOST_CHECK_EQUAL(tracker.events("P2").stopped, 1);
    BOOST_CHECK_EQUAL(tracker.events("P2").shut, 0);

    // Staying shut is not a new event.
    tracker.beginTimeStep(sched, 0, now);
    tracker.accumulate(now);
    BOOST_CHECK(noEvents(tracker.events("P1")));
    BOOST_CHECK(noEvents(tracker.events("P2")));
}

BOOST_FIXTURE_TEST_CASE(ConnectionsOpened, Setup)
{
    auto tracker = Opm::WellPerformanceEvents{};

    auto before = allOpen();
    before.wells["P1"].openCompletions = {1, 2};
    before.wells["P2"] = Entry { Opm::WellStatus::SHUT, {} };
    before.wells["P3"] = Entry { Opm::WellStatus::STOP, {} };
    tracker.beginTimeStep(sched, 0, before);

    auto now = allOpen();
    now.wells["P2"].status = Opm::WellStatus::SHUT;
    now.wells["P3"].status = Opm::WellStatus::STOP;
    tracker.accumulate(now);

    BOOST_CHECK_EQUAL(tracker.events("P1").connsOpened, 2);
    BOOST_CHECK_EQUAL(tracker.events("P1").connsClosed, 0);

    // WPWE1 is only reported for a well that is neither shut nor stopped.
    BOOST_CHECK_EQUAL(tracker.events("P2").connsOpened, 0);
    BOOST_CHECK_EQUAL(tracker.events("P3").connsOpened, 0);
}

BOOST_FIXTURE_TEST_CASE(SynchroniseAbsorbsDeckChanges, Setup)
{
    auto tracker = Opm::WellPerformanceEvents{};
    tracker.beginTimeStep(sched, 0, allOpen());

    auto now = allOpen();
    now.wells["P1"].openCompletions = {1, 2};
    now.wells["P2"].status = Opm::WellStatus::SHUT;

    tracker.synchronise(now);
    tracker.accumulate(now);

    BOOST_CHECK(noEvents(tracker.events("P1")));
    BOOST_CHECK(noEvents(tracker.events("P2")));
}

BOOST_FIXTURE_TEST_CASE(BeginTimeStepResetsEvents, Setup)
{
    auto tracker = Opm::WellPerformanceEvents{};
    tracker.beginTimeStep(sched, 0, allOpen());

    auto now = allOpen();
    now.wells["P1"].openCompletions = {1, 2, 3};
    tracker.accumulate(now);
    BOOST_CHECK_EQUAL(tracker.events("P1").connsClosed, 1);

    tracker.beginTimeStep(sched, 0, now);
    BOOST_CHECK(noEvents(tracker.events("P1")));
}

BOOST_FIXTURE_TEST_CASE(UnknownWellSkipped, Setup)
{
    // A well without a reference status, e.g., one that appeared through
    // ACTIONX, is ignored rather than compared against nothing.
    auto tracker = Opm::WellPerformanceEvents{};

    auto before = allOpen();
    before.wells.erase("P1");
    tracker.beginTimeStep(sched, 0, before);

    auto now = allOpen();
    now.wells["P1"] = Entry { Opm::WellStatus::SHUT, {} };
    BOOST_CHECK_NO_THROW(tracker.accumulate(now));
    BOOST_CHECK(noEvents(tracker.events("P1")));
}

BOOST_FIXTURE_TEST_CASE(TypeSwitch, Setup)
{
    auto tracker = Opm::WellPerformanceEvents{};

    // First sight of the wells: nothing to compare against.
    tracker.beginTimeStep(sched, 0, allOpen());
    for (const auto* well : { "P1", "P2", "P3", "I1" }) {
        BOOST_CHECK(noEvents(tracker.events(well)));
    }

    tracker.beginTimeStep(sched, 1, allOpen());

    BOOST_CHECK_EQUAL(tracker.events("I1").injectorToProducer, 1);
    BOOST_CHECK_EQUAL(tracker.events("I1").producerToInjector, 0);

    BOOST_CHECK_EQUAL(tracker.events("P3").producerToInjector, 1);
    BOOST_CHECK_EQUAL(tracker.events("P3").injectorToProducer, 0);

    BOOST_CHECK(noEvents(tracker.events("P1")));
    BOOST_CHECK(noEvents(tracker.events("P2")));

    tracker.commitTimeStep(sched, 1);

    // The switch is not repeated on the next time step of the same report
    // step, nor on a later report step without a change.
    tracker.beginTimeStep(sched, 1, allOpen());
    BOOST_CHECK(noEvents(tracker.events("I1")));

    tracker.beginTimeStep(sched, 2, allOpen());
    BOOST_CHECK(noEvents(tracker.events("I1")));
    BOOST_CHECK(noEvents(tracker.events("P3")));
}

BOOST_FIXTURE_TEST_CASE(TypeSwitchSurvivesRetry, Setup)
{
    auto tracker = Opm::WellPerformanceEvents{};
    tracker.beginTimeStep(sched, 0, allOpen());
    tracker.commitTimeStep(sched, 0);

    for (int attempt = 0; attempt < 3; ++attempt) {
        tracker.beginTimeStep(sched, 1, allOpen());
        BOOST_CHECK_EQUAL(tracker.events("I1").injectorToProducer, 1);
        BOOST_CHECK_EQUAL(tracker.events("P3").producerToInjector, 1);
    }

    tracker.commitTimeStep(sched, 1);
    tracker.beginTimeStep(sched, 2, allOpen());
    BOOST_CHECK(noEvents(tracker.events("I1")));
    BOOST_CHECK(noEvents(tracker.events("P3")));
}

BOOST_AUTO_TEST_CASE(ActionTypeSwitchWithinReportStep)
{
    auto text = deckString();
    text.insert(text.find("DATES"), R"(
ACTIONX
 'CONVERT' 1 /
 TIME > 1 /
/
WCONINJE
 'P1' 'WATER' 'OPEN' 'RATE' 100.0 1* 400.0 /
/
WCONPROD
 'I1' 'OPEN' 'ORAT' 100.0 4* 50.0 /
/
ENDACTIO
)");
    const auto deck = Opm::Parser{}.parseString(text);
    const Opm::EclipseState es{deck};
    Opm::Schedule schedule{deck, es, std::make_shared<Opm::Python>()};
    auto tracker = Opm::WellPerformanceEvents{};
    tracker.beginTimeStep(schedule, 0, allOpen());
    tracker.commitTimeStep(schedule, 0);

    const auto action = schedule[0].actions.get()["CONVERT"];
    schedule.applyAction(0, action, Opm::Action::Result{true}.matches(),
                         std::unordered_map<std::string, double>{}, true);
    BOOST_REQUIRE(schedule.getWell("P1", 0).isInjector());
    BOOST_REQUIRE(schedule.getWell("I1", 0).isProducer());

    // Both the first attempt and its retry must see the ACTIONX conversion,
    // without waiting for a new report interval.
    for (int attempt = 0; attempt < 2; ++attempt) {
        tracker.beginTimeStep(schedule, 0, allOpen());
        BOOST_CHECK_EQUAL(tracker.events("P1").producerToInjector, 1);
        BOOST_CHECK_EQUAL(tracker.events("I1").injectorToProducer, 1);
    }

    tracker.commitTimeStep(schedule, 0);
    tracker.beginTimeStep(schedule, 0, allOpen());
    BOOST_CHECK(noEvents(tracker.events("P1")));
    BOOST_CHECK(noEvents(tracker.events("I1")));
}

BOOST_AUTO_TEST_CASE(SerializationTestObject)
{
    const auto obj = Opm::WellPerformanceEvents::serializationTestObject();

    BOOST_CHECK(obj == obj);
    BOOST_CHECK(!(obj == Opm::WellPerformanceEvents{}));
    BOOST_CHECK(!noEvents(obj.events("W1")));
}
