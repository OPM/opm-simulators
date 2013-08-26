/* Copyright 2013 Uni Research AS
 * This file is licensed under GPL3, see http://www.opm-project.org/
*/
#include <config.h>

/* --- Boost.Test boilerplate --- */
#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif

#define NVERBOSE  // Suppress own messages when throw()ing

#define BOOST_TEST_MODULE EventTest
#include <boost/test/unit_test.hpp>

/* --- our own headers --- */
#include <opm/core/utility/Event.hpp>

using namespace std;
using namespace Opm;

// idiomatic implementation of generator and receiver
struct EventGenerator {
    EventSource eventSource_;
    Event& event () { return eventSource_; }
    void action () { eventSource_.signal (); }
};

struct EventReceiver {
    int numOfCalls;
    EventReceiver () : numOfCalls (0) { }
    void handler () { ++numOfCalls; }
private:
    // make sure bind() doesn't implement copy constructor
    EventReceiver (EventReceiver&);
};

// declare a generator, a receiver and connect them
struct EventFixture {
    EventGenerator gen;
    EventReceiver recv;
    void register_handler () {
        gen.event().add<EventReceiver,&EventReceiver::handler>(recv);
    }

    EventFixture () {
        register_handler();
    }
};

BOOST_FIXTURE_TEST_SUITE(EventTest, EventFixture)

BOOST_AUTO_TEST_CASE(none)
{
    BOOST_REQUIRE_EQUAL (recv.numOfCalls, 0);
}

BOOST_AUTO_TEST_CASE(once)
{
    gen.action();
    BOOST_REQUIRE_EQUAL (recv.numOfCalls, 1);
}

BOOST_AUTO_TEST_CASE(twice)
{
    gen.action();
    gen.action();
    BOOST_REQUIRE_EQUAL (recv.numOfCalls, 2);
}

BOOST_AUTO_TEST_CASE(reg_twice)
{
    register_handler();
    gen.action();
    BOOST_REQUIRE_EQUAL (recv.numOfCalls, 2);
}

BOOST_AUTO_TEST_SUITE_END()
