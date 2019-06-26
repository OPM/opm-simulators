/*
  Copyright 2018 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2018 Equinor.

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

#define BOOST_TEST_MODULE TestDeferredLogger

#include <boost/test/unit_test.hpp>

#include <opm/simulators/utils/DeferredLogger.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/OpmLog/LogBackend.hpp>
#include <opm/common/OpmLog/CounterLog.hpp>
#include <opm/common/OpmLog/TimerLog.hpp>
#include <opm/common/OpmLog/StreamLog.hpp>
#include <opm/common/OpmLog/LogUtil.hpp>

using namespace Opm;

void initLogger(std::ostringstream& log_stream) {
    OpmLog::removeAllBackends();
    std::shared_ptr<CounterLog> counter = std::make_shared<CounterLog>();
    std::shared_ptr<StreamLog> streamLog = std::make_shared<StreamLog>( log_stream , Log::DefaultMessageTypes);

    OpmLog::addBackend("COUNTER" , counter);
    OpmLog::addBackend("STREAM" , streamLog);
    BOOST_CHECK_EQUAL( true , OpmLog::hasBackend("COUNTER"));
    BOOST_CHECK_EQUAL( true , OpmLog::hasBackend("STREAM"));

    streamLog->setMessageFormatter(std::make_shared<SimpleMessageFormatter>(true, false));
    streamLog->setMessageLimiter(std::make_shared<MessageLimiter>(2));
}


BOOST_AUTO_TEST_CASE(deferredlogger)
{

    const std::string expected = Log::prefixMessage(Log::MessageType::Info, "info 1") + "\n"
        + Log::prefixMessage(Log::MessageType::Warning, "warning 1") + "\n"
        + Log::prefixMessage(Log::MessageType::Error, "error 1") + "\n"
        + Log::prefixMessage(Log::MessageType::Error, "error 2") + "\n"
        + Log::prefixMessage(Log::MessageType::Problem, "problem 1") + "\n"
        + Log::prefixMessage(Log::MessageType::Bug, "bug 1") + "\n"
        + Log::prefixMessage(Log::MessageType::Debug, "debug 1") + "\n"
        + Log::prefixMessage(Log::MessageType::Note, "note 1") + "\n"
        + Log::prefixMessage(Log::MessageType::Note, "note 2") + "\n"
        + Log::prefixMessage(Log::MessageType::Note, "note 3") + "\n"
        + Log::prefixMessage(Log::MessageType::Note, "Message limit reached for message tag: tagme") + "\n";

    std::ostringstream log_stream;
    initLogger(log_stream);
    auto deferred_logger = Opm::DeferredLogger();
    deferred_logger.info("info 1");
    deferred_logger.warning("warning 1");
    deferred_logger.error("error 1");
    deferred_logger.error("error 2");
    deferred_logger.problem("problem 1");
    deferred_logger.bug("bug 1");
    deferred_logger.debug("debug 1");
    deferred_logger.note("note 1");
    deferred_logger.note("tagme", "note 2");
    deferred_logger.note("tagme", "note 3");
    deferred_logger.note("tagme", "note 3");
    deferred_logger.note("tagme", "note 3");
    deferred_logger.note("tagme", "note 3");
    deferred_logger.note("tagme", "note 3");
    deferred_logger.note("tagme", "note 3");

    deferred_logger.logMessages();

    auto counter = OpmLog::getBackend<CounterLog>("COUNTER");
    BOOST_CHECK_EQUAL( 1 , counter->numMessages(Log::MessageType::Warning) );
    BOOST_CHECK_EQUAL( 1 , counter->numMessages(Log::MessageType::Info) );
    BOOST_CHECK_EQUAL( 2 , counter->numMessages(Log::MessageType::Error) );
    BOOST_CHECK_EQUAL( 1 , counter->numMessages(Log::MessageType::Problem) );
    BOOST_CHECK_EQUAL( 1 , counter->numMessages(Log::MessageType::Bug) );
    BOOST_CHECK_EQUAL( 1 , counter->numMessages(Log::MessageType::Debug) );
    BOOST_CHECK_EQUAL( 8 , counter->numMessages(Log::MessageType::Note) );

    BOOST_CHECK_EQUAL(log_stream.str(), expected);

}
