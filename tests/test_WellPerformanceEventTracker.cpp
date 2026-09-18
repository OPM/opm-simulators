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

#define BOOST_TEST_MODULE TestWellPerformanceEventTracker

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <boost/test/unit_test.hpp>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <opm/simulators/wells/WellPerformanceEventTracker.hpp>

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
 5 4 2 5 /
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

// Cell index identifying the connection that carries completion 'c' in the
// one-connection-per-completion wells most of these tests use.
constexpr std::size_t cell(const int c)
{
    return 9 + static_cast<std::size_t>(c);
}

// Connections able to flow, one per completion.
void setOpen(Entry& e, std::initializer_list<int> completions)
{
    e.openConnections.clear();
    e.openCompletions.clear();

    for (const int c : completions) {
        e.openConnections.push_back(cell(c));
        e.openCompletions.push_back(c);
    }
}

// Connections able to flow, given as (cell index, completion number).  Used
// where COMPLUMP puts more than one connection on a completion.
void setOpenLumped(Entry& e, std::initializer_list<std::pair<std::size_t, int>> conns)
{
    e.openConnections.clear();
    e.openCompletions.clear();

    for (const auto& [cellIx, complnum] : conns) {
        e.openConnections.push_back(cellIx);
        e.openCompletions.push_back(complnum);
    }
}

Opm::WellStatusSnapshot allOpen()
{
    auto snap = Opm::WellStatusSnapshot{};

    for (const auto* well : { "P1", "P2", "P3", "I1" }) {
        snap.wells[well] = Entry { Opm::WellStatus::OPEN, {}, {}, {}, 1 };
        setOpen(snap.wells[well], {1, 2, 3, 4});
    }

    return snap;
}

bool noEvents(const Opm::data::WellPerformanceEvents& events)
{
    return events == Opm::data::WellPerformanceEvents{};
}

} // Anonymous namespace

BOOST_FIXTURE_TEST_CASE(NoChange, Setup)
{
    auto tracker = Opm::WellPerformanceEventTracker{};

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
    auto tracker = Opm::WellPerformanceEventTracker{};
    tracker.beginTimeStep(sched, 0, allOpen());

    auto now = allOpen();
    setOpen(now.wells["P1"], {1, 2, 3});
    tracker.accumulate(now);

    {
        const auto& ev = tracker.events("P1");
        BOOST_CHECK_EQUAL(ev.connsClosed, 1);
        BOOST_CHECK_EQUAL(ev.closedToBottom, 0);
        BOOST_CHECK_EQUAL(ev.connsOpened, 0);
        BOOST_CHECK_EQUAL(ev.shut, 0);
    }

    // Events accumulate within the time step.
    setOpen(now.wells["P1"], {1, 2});
    tracker.accumulate(now);

    BOOST_CHECK_EQUAL(tracker.events("P1").connsClosed, 2);
    BOOST_CHECK_EQUAL(tracker.events("P1").closedToBottom, 0);

    BOOST_CHECK(noEvents(tracker.events("P2")));
}

BOOST_FIXTURE_TEST_CASE(AllConnectionsClosedByConWorkovers, Setup)
{
    // Closing the last open connection sets WPWE3 when every connection was
    // closed by an individual CON workover.  The well is also shut (WPWE7).
    auto tracker = Opm::WellPerformanceEventTracker{};

    auto before = allOpen();
    setOpen(before.wells["P1"], {1});
    tracker.beginTimeStep(sched, 0, before);

    auto now = before;
    now.wells["P1"] = Entry { Opm::WellStatus::SHUT, {}, {}, {}, 1 };
    tracker.accumulate(now);

    const auto& ev = tracker.events("P1");
    BOOST_CHECK_EQUAL(ev.connsClosed, 1);
    BOOST_CHECK_EQUAL(ev.closedToBottom, 1);
    BOOST_CHECK_EQUAL(ev.shut, 1);
    BOOST_CHECK_EQUAL(ev.stopped, 0);
}

BOOST_FIXTURE_TEST_CASE(IndividualClosuresLeaveOneCompletion, Setup)
{
    // Closures that leave the well flowing over its topmost completion alone
    // have closed it to the bottom, and count towards WPWE2 as well.  Measured
    // in WPWE-PLUSCON-CLOSURE-COUNT at day 40 and
    // WPWE-DECK-CONNECTION-CHANGES at day 311.
    auto tracker = Opm::WellPerformanceEventTracker{};
    tracker.beginTimeStep(sched, 0, allOpen());

    auto now = allOpen();
    setOpen(now.wells["P2"], {1});
    tracker.accumulate(now);

    const auto& ev = tracker.events("P2");
    BOOST_CHECK_EQUAL(ev.closedToBottom, 1);
    BOOST_CHECK_EQUAL(ev.connsClosed, 3);
    BOOST_CHECK_EQUAL(ev.shut, 0);
}

BOOST_FIXTURE_TEST_CASE(DeckClosingLastConnectionSetsWPWE3, Setup)
{
    // A deck change counts towards WPWE2 and, unlike a workover, never sets
    // WPWE7: the well status the deck chooses is not an event.  It does leave
    // the well closed to the bottom, and WPWE3 does not ask what closed the
    // connections -- WPWE-MIXED-CLOSURE-CAUSES reports it for a well whose
    // other three connections the deck had shut.
    //
    // This pins the interpretation that a deck change closing the final
    // connection can set WPWE3.
    auto tracker = Opm::WellPerformanceEventTracker{};

    auto before = allOpen();
    setOpen(before.wells["P1"], {1});
    tracker.beginTimeStep(sched, 0, before);

    auto now = before;
    setOpen(now.wells["P1"], {});
    now.wells["P1"].status = Opm::WellStatus::SHUT;
    tracker.applyDeckChanges(now);

    const auto& ev = tracker.events("P1");
    BOOST_CHECK_EQUAL(ev.connsClosed, 1);
    BOOST_CHECK_EQUAL(ev.closedToBottom, 1);
    BOOST_CHECK_EQUAL(ev.shut, 0);
}

BOOST_FIXTURE_TEST_CASE(PlusConReachIsNotCounted, Setup)
{
    // A '+CON' workover closes the completion whose limit was violated and
    // reaches past it to the rest of the well below.  Only the violation
    // counts towards WPWE2, and reaching past it closed the well to the
    // bottom however much is left open above.
    auto tracker = Opm::WellPerformanceEventTracker{};
    tracker.beginTimeStep(sched, 0, allOpen());

    auto now = allOpen();
    setOpen(now.wells["P1"], {1, 2});
    now.wells["P1"].closedBelowOffender = {cell(4)};  // 3 offended, 4 lies below it
    tracker.accumulate(now);

    const auto& ev = tracker.events("P1");
    BOOST_CHECK_EQUAL(ev.connsClosed, 1);
    BOOST_CHECK_EQUAL(ev.closedToBottom, 1);
}

BOOST_FIXTURE_TEST_CASE(PlusConReachesPastAlreadyClosedConnections, Setup)
{
    // The deck has already shut completions 4 and 5 when a '+CON' fires on 3.
    // Its reach covers them whether or not they were still able to flow, so
    // the well is closed to the bottom exactly as it would be had they been
    // open.  Completions 1 and 2 keep flowing, so the topmost-completion route
    // does not apply and only the reach can report this.
    auto tracker = Opm::WellPerformanceEventTracker{};

    auto before = allOpen();
    setOpen(before.wells["P1"], {1, 2, 3});
    tracker.beginTimeStep(sched, 0, before);

    auto now = before;
    setOpen(now.wells["P1"], {1, 2});
    now.wells["P1"].closedBelowOffender = {cell(4), cell(5)};
    tracker.accumulate(now);

    const auto& ev = tracker.events("P1");
    BOOST_CHECK_EQUAL(ev.connsClosed, 1);      // the offender alone
    BOOST_CHECK_EQUAL(ev.closedToBottom, 1);
}

BOOST_FIXTURE_TEST_CASE(RepeatedPlusConCountsEveryOffender, Setup)
{
    // If multiple '+CON' limits fail in one cascade, each offending
    // connection counts towards WPWE2.  Connections reached from either
    // offender do not.
    auto tracker = Opm::WellPerformanceEventTracker{};
    tracker.beginTimeStep(sched, 0, allOpen());

    auto now = allOpen();
    setOpen(now.wells["P1"], {});
    now.wells["P1"].closedBelowOffender = {cell(2), cell(4)};
    tracker.accumulate(now);

    const auto& ev = tracker.events("P1");
    BOOST_CHECK_EQUAL(ev.connsClosed, 2);
    BOOST_CHECK_EQUAL(ev.closedToBottom, 1);
}

BOOST_FIXTURE_TEST_CASE(RemainingCompletionMustBeTheTopmost, Setup)
{
    // One completion left able to flow is not enough: closures above it mean
    // the wellbore was not closed from the bottom.  Here the deck has shut
    // completions 1 and 2 and a workover closes 4, leaving the well flowing
    // on 3 alone.  Measured in WPWE-MIXED-CLOSURE-CAUSES, where C-1H keeps its
    // fourth connection with three shut above it and reports no WPWE3.
    auto tracker = Opm::WellPerformanceEventTracker{};

    auto before = allOpen();
    setOpen(before.wells["P1"], {3, 4});
    tracker.beginTimeStep(sched, 0, before);

    auto now = before;
    setOpen(now.wells["P1"], {3});
    tracker.accumulate(now);

    const auto& ev = tracker.events("P1");
    BOOST_CHECK_EQUAL(ev.connsClosed, 1);
    BOOST_CHECK_EQUAL(ev.closedToBottom, 0);

    // Closing that last one does leave the well closed to the bottom.
    tracker.commitTimeStep(sched, 0);
    tracker.beginTimeStep(sched, 0, now);
    setOpen(now.wells["P1"], {});
    tracker.accumulate(now);
    BOOST_CHECK_EQUAL(tracker.events("P1").closedToBottom, 1);
}

BOOST_FIXTURE_TEST_CASE(ClosedTailIsNotClosedToBottom, Setup)
{
    // Closing the bottom connections is not enough on its own: two
    // completions are still able to flow.
    auto tracker = Opm::WellPerformanceEventTracker{};
    tracker.beginTimeStep(sched, 0, allOpen());

    auto now = allOpen();
    setOpen(now.wells["P3"], {1, 2});
    tracker.accumulate(now);

    const auto& ev = tracker.events("P3");
    BOOST_CHECK_EQUAL(ev.closedToBottom, 0);
    BOOST_CHECK_EQUAL(ev.connsClosed, 2);
}

BOOST_FIXTURE_TEST_CASE(LumpedCompletionCountsConnections, Setup)
{
    // COMPLUMP maps several connections onto one completion number, so a
    // snapshot carries that number once per connection.  WPWE1 and WPWE2 count
    // connections, so closing or opening a lumped completion must count every
    // connection it carries rather than the completion itself.  WPWE3 instead
    // works on the distinct numbers: one lumped completion left able to flow
    // is still closed to the bottom.  Measured in
    // WPWE-COMPLUMP-CONNECTION-COUNT at day 13, where the reference reports
    // WPWE2=2 with WPWE3=1.
    auto tracker = Opm::WellPerformanceEventTracker{};

    auto before = allOpen();
    setOpenLumped(before.wells["P1"], {{10, 1}, {11, 1}, {12, 2}, {13, 2}});
    tracker.beginTimeStep(sched, 0, before);

    auto now = before;
    setOpenLumped(now.wells["P1"], {{10, 1}, {11, 1}});
    tracker.accumulate(now);
    BOOST_CHECK_EQUAL(tracker.events("P1").connsClosed, 2);
    BOOST_CHECK_EQUAL(tracker.events("P1").closedToBottom, 1);

    tracker.commitTimeStep(sched, 0);
    tracker.beginTimeStep(sched, 0, now);
    tracker.accumulate(before);
    BOOST_CHECK_EQUAL(tracker.events("P1").connsOpened, 2);
}

BOOST_FIXTURE_TEST_CASE(RenumberingCompletionsIsNotAnEvent, Setup)
{
    // COMPLUMP may put the same connections on different completion numbers
    // partway through a run.  Nothing opens or closes, so nothing is reported:
    // the indicators follow the connections, not the numbering.
    auto tracker = Opm::WellPerformanceEventTracker{};

    auto before = allOpen();
    setOpenLumped(before.wells["P1"], {{10, 1}, {11, 2}, {12, 3}, {13, 4}});
    tracker.beginTimeStep(sched, 0, before);

    auto now = before;
    setOpenLumped(now.wells["P1"], {{10, 1}, {11, 1}, {12, 1}, {13, 1}});
    tracker.accumulate(now);

    BOOST_CHECK(noEvents(tracker.events("P1")));
}

BOOST_FIXTURE_TEST_CASE(ConnectionsSwappedWithinOneCompletion, Setup)
{
    // One connection of a lumped completion closes while another opens.  The
    // set of open completion numbers never changes, but a connection did close
    // and a connection did open, and both count.  Completion 2 carries the
    // lump, so the well is not left on its topmost completion and WPWE3 stays
    // silent.
    auto tracker = Opm::WellPerformanceEventTracker{};

    auto before = allOpen();
    setOpenLumped(before.wells["P1"], {{10, 1}, {11, 2}, {12, 2}});
    tracker.beginTimeStep(sched, 0, before);

    auto now = before;
    setOpenLumped(now.wells["P1"], {{10, 1}, {12, 2}, {13, 2}});
    tracker.accumulate(now);

    BOOST_CHECK_EQUAL(tracker.events("P1").connsOpened, 1);
    BOOST_CHECK_EQUAL(tracker.events("P1").connsClosed, 1);
    BOOST_CHECK_EQUAL(tracker.events("P1").closedToBottom, 0);
}

BOOST_FIXTURE_TEST_CASE(ClosuresAreNotRepeated, Setup)
{
    auto tracker = Opm::WellPerformanceEventTracker{};
    tracker.beginTimeStep(sched, 0, allOpen());

    auto now = allOpen();
    setOpen(now.wells["P1"], {1, 2});
    tracker.accumulate(now);
    BOOST_CHECK_EQUAL(tracker.events("P1").connsClosed, 2);
    BOOST_CHECK_EQUAL(tracker.events("P1").closedToBottom, 0);

    // Connections closed in an earlier step count once, and must not add to
    // the count for the closure this step makes.
    tracker.commitTimeStep(sched, 0);
    tracker.beginTimeStep(sched, 0, now);
    setOpen(now.wells["P1"], {1});
    tracker.accumulate(now);
    BOOST_CHECK_EQUAL(tracker.events("P1").connsClosed, 1);
    BOOST_CHECK_EQUAL(tracker.events("P1").closedToBottom, 1);
}

BOOST_FIXTURE_TEST_CASE(WellShutAndStopped, Setup)
{
    auto tracker = Opm::WellPerformanceEventTracker{};
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
    tracker.commitTimeStep(sched, 0);
    tracker.beginTimeStep(sched, 0, now);
    tracker.accumulate(now);
    BOOST_CHECK(noEvents(tracker.events("P1")));
    BOOST_CHECK(noEvents(tracker.events("P2")));
}

BOOST_FIXTURE_TEST_CASE(ConnectionsOpened, Setup)
{
    auto tracker = Opm::WellPerformanceEventTracker{};

    auto before = allOpen();
    setOpen(before.wells["P1"], {1, 2});
    before.wells["P2"] = Entry { Opm::WellStatus::SHUT, {}, {}, {}, 1 };
    before.wells["P3"] = Entry { Opm::WellStatus::STOP, {}, {}, {}, 1 };
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

BOOST_FIXTURE_TEST_CASE(DeckChangesConnectionsButNotStatus, Setup)
{
    // The deck closing a connection is a WPWE2 event; the deck shutting a
    // well is not a WPWE7 one.
    auto tracker = Opm::WellPerformanceEventTracker{};
    tracker.beginTimeStep(sched, 0, allOpen());

    auto now = allOpen();
    setOpen(now.wells["P1"], {1, 2});
    now.wells["P2"].status = Opm::WellStatus::SHUT;

    tracker.applyDeckChanges(now);
    tracker.accumulate(now);

    BOOST_CHECK_EQUAL(tracker.events("P1").connsClosed, 2);
    BOOST_CHECK_EQUAL(tracker.events("P1").closedToBottom, 0);
    BOOST_CHECK(noEvents(tracker.events("P2")));
}

BOOST_FIXTURE_TEST_CASE(BeginTimeStepResetsEvents, Setup)
{
    auto tracker = Opm::WellPerformanceEventTracker{};
    tracker.beginTimeStep(sched, 0, allOpen());

    auto now = allOpen();
    setOpen(now.wells["P1"], {1, 2, 3});
    tracker.accumulate(now);
    BOOST_CHECK_EQUAL(tracker.events("P1").connsClosed, 1);

    tracker.commitTimeStep(sched, 0);
    tracker.beginTimeStep(sched, 0, now);
    BOOST_CHECK(noEvents(tracker.events("P1")));
}

BOOST_FIXTURE_TEST_CASE(ConnectionChangesSurviveRetry, Setup)
{
    auto tracker = Opm::WellPerformanceEventTracker{};
    tracker.beginTimeStep(sched, 0, allOpen());
    tracker.commitTimeStep(sched, 0);

    auto changed = allOpen();
    setOpen(changed.wells["P1"], {1, 2, 3});

    // Do not commit between attempts: each beginTimeStep() represents a retry
    // after the previous attempt failed.
    for (int attempt = 0; attempt < 3; ++attempt) {
        tracker.beginTimeStep(sched, 0, changed);
        BOOST_CHECK_EQUAL(tracker.events("P1").connsClosed, 1);
    }
}

BOOST_FIXTURE_TEST_CASE(WtestReopenSurvivesRetry, Setup)
{
    auto tracker = Opm::WellPerformanceEventTracker{};

    auto closed = allOpen();
    setOpen(closed.wells["P1"], {1});
    tracker.beginTimeStep(sched, 0, closed);
    tracker.commitTimeStep(sched, 0);

    auto reopened = closed;
    setOpen(reopened.wells["P1"], {1, 2});

    for (int attempt = 0; attempt < 3; ++attempt) {
        tracker.beginTimeStep(sched, 0, closed);
        tracker.accumulate(reopened);

        BOOST_CHECK_EQUAL(tracker.events("P1").connsOpened, 1);
        BOOST_CHECK_EQUAL(tracker.events("P1").connsClosed, 0);
    }
}

BOOST_FIXTURE_TEST_CASE(ConvergenceClosureSurvivesRetry, Setup)
{
    auto tracker = Opm::WellPerformanceEventTracker{};
    tracker.beginTimeStep(sched, 0, allOpen());
    tracker.commitTimeStep(sched, 0);

    auto closed = allOpen();
    closed.wells["P1"].status = Opm::WellStatus::STOP;
    closed.wells["P2"].status = Opm::WellStatus::SHUT;

    // AdaptiveTimeStepping persists these simulator decisions before retrying
    // the time step, so beginTimeStep() sees them in its starting snapshot.
    tracker.recordWellStatusChangeForRetry("P1", Opm::WellStatus::STOP);
    tracker.recordWellStatusChangeForRetry("P2", Opm::WellStatus::SHUT);

    for (int attempt = 0; attempt < 3; ++attempt) {
        tracker.beginTimeStep(sched, 0, closed);
        BOOST_CHECK_EQUAL(tracker.events("P1").stopped, 1);
        BOOST_CHECK_EQUAL(tracker.events("P1").shut, 0);
        BOOST_CHECK_EQUAL(tracker.events("P2").stopped, 0);
        BOOST_CHECK_EQUAL(tracker.events("P2").shut, 1);
    }

    tracker.commitTimeStep(sched, 0);
    tracker.beginTimeStep(sched, 0, closed);
    BOOST_CHECK(noEvents(tracker.events("P1")));
    BOOST_CHECK(noEvents(tracker.events("P2")));
}

BOOST_FIXTURE_TEST_CASE(UnknownWellSkipped, Setup)
{
    // A well without a reference status, e.g., one that appeared through
    // ACTIONX, is ignored rather than compared against nothing.
    auto tracker = Opm::WellPerformanceEventTracker{};

    auto before = allOpen();
    before.wells.erase("P1");
    tracker.beginTimeStep(sched, 0, before);

    auto now = allOpen();
    now.wells["P1"] = Entry { Opm::WellStatus::SHUT, {}, {}, {}, 1 };
    BOOST_CHECK_NO_THROW(tracker.accumulate(now));
    BOOST_CHECK(noEvents(tracker.events("P1")));
}

BOOST_FIXTURE_TEST_CASE(TypeSwitch, Setup)
{
    auto tracker = Opm::WellPerformanceEventTracker{};

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
    auto tracker = Opm::WellPerformanceEventTracker{};
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
    auto tracker = Opm::WellPerformanceEventTracker{};
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

BOOST_AUTO_TEST_CASE(WellEnteringTheScheduleIsNotDrilled)
{
    // Flow has no drilling queue, and a well introduced through WELSPECS while
    // the run is under way is not a WPWE0 event either.
    auto text = deckString();
    text.insert(text.find("DATES\n 1 'MAR' 2020 /"), R"(
WELSPECS
 'P4' 'G1' 1 1 1* 'OIL' /
/
COMPDAT
 'P4' 1 1 1 4 'OPEN' 1* 1* 0.2 /
/
WCONPROD
 'P4' 'OPEN' 'ORAT' 100.0 4* 50.0 /
/
)");
    const auto deck = Opm::Parser{}.parseString(text);
    const Opm::EclipseState es{deck};
    Opm::Schedule schedule{deck, es, std::make_shared<Opm::Python>()};

    auto tracker = Opm::WellPerformanceEventTracker{};

    tracker.beginTimeStep(schedule, 0, allOpen());
    for (const auto* well : { "P1", "P2", "P3", "I1" }) {
        BOOST_CHECK(noEvents(tracker.events(well)));
    }
    tracker.commitTimeStep(schedule, 0);

    // P4 enters the schedule at the second report step.
    tracker.beginTimeStep(schedule, 1, allOpen());
    BOOST_CHECK(noEvents(tracker.events("P4")));
    BOOST_CHECK(noEvents(tracker.events("P1")));

    tracker.commitTimeStep(schedule, 1);
    tracker.beginTimeStep(schedule, 1, allOpen());
    BOOST_CHECK(noEvents(tracker.events("P4")));
}

BOOST_AUTO_TEST_CASE(SerializationTestObject)
{
    const auto obj = Opm::WellPerformanceEventTracker::serializationTestObject();

    BOOST_CHECK(obj == obj);
    BOOST_CHECK(!(obj == Opm::WellPerformanceEventTracker{}));
    BOOST_CHECK(!noEvents(obj.events("W1")));
}
