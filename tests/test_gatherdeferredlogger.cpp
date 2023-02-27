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

#define BOOST_TEST_MODULE TestGatherDeferredLogger
#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>

#include <opm/simulators/utils/gatherDeferredLogger.hpp>
#include <dune/common/parallel/mpihelper.hh>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/OpmLog/LogBackend.hpp>
#include <opm/common/OpmLog/CounterLog.hpp>
#include <opm/common/OpmLog/TimerLog.hpp>
#include <opm/common/OpmLog/StreamLog.hpp>
#include <opm/common/OpmLog/LogUtil.hpp>

using namespace Opm;

#if HAVE_MPI
struct MPIError
{
    MPIError(std::string s, int e) : errorstring(std::move(s)), errorcode(e){}
    std::string errorstring;
    int errorcode;
};

void MPI_err_handler(MPI_Comm*, int* err_code, ...)
{
    std::vector<char> err_string(MPI_MAX_ERROR_STRING);
    int err_length;
    MPI_Error_string(*err_code, err_string.data(), &err_length);
    std::string s(err_string.data(), err_length);
    std::cerr << "An MPI Error ocurred:" << std::endl << s << std::endl;
    throw MPIError(s, *err_code);
}
#endif

bool
init_unit_test_func()
{
    return true;
}

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



BOOST_AUTO_TEST_CASE(NoMessages)
{
    const auto& cc = Dune::MPIHelper::getCommunication();
    std::ostringstream log_stream;
    initLogger(log_stream);

    Opm::DeferredLogger local_deferredlogger;

    Opm::DeferredLogger global_deferredlogger = gatherDeferredLogger(local_deferredlogger, cc);

    if (cc.rank() == 0) {

        global_deferredlogger.logMessages();

        auto counter = OpmLog::getBackend<CounterLog>("COUNTER");
        BOOST_CHECK_EQUAL( 0 , counter->numMessages(Log::MessageType::Info) );

        std::string expected;
        BOOST_CHECK_EQUAL(log_stream.str(), expected);
    }
}

BOOST_AUTO_TEST_CASE(VariableNumberOfMessages)
{
    const auto& cc = Dune::MPIHelper::getCommunication();

    std::ostringstream log_stream;
    initLogger(log_stream);

    Opm::DeferredLogger local_deferredlogger;
    if (cc.rank() == 1) {
        local_deferredlogger.info("info from rank " + std::to_string(cc.rank()));
        local_deferredlogger.warning("warning from rank " + std::to_string(cc.rank()));
    } else if (cc.rank() == 2) {
        local_deferredlogger.bug("tagme", "bug from rank " + std::to_string(cc.rank()));
        local_deferredlogger.bug("tagme", "bug from rank " + std::to_string(cc.rank()));
        local_deferredlogger.bug("tagme", "bug from rank " + std::to_string(cc.rank()));
        local_deferredlogger.bug("tagme", "bug from rank " + std::to_string(cc.rank()));
    }

    Opm::DeferredLogger global_deferredlogger = gatherDeferredLogger(local_deferredlogger, cc);

    if (cc.rank() == 0) {

        global_deferredlogger.logMessages();

        auto counter = OpmLog::getBackend<CounterLog>("COUNTER");
        BOOST_CHECK_EQUAL( 1 , counter->numMessages(Log::MessageType::Info) );
        BOOST_CHECK_EQUAL( 1 , counter->numMessages(Log::MessageType::Warning) );
        BOOST_CHECK_EQUAL( 4 , counter->numMessages(Log::MessageType::Bug) );

        const std::string expected = Log::prefixMessage(Log::MessageType::Info, "info from rank 1") + "\n"
            + Log::prefixMessage(Log::MessageType::Warning, "warning from rank 1") + "\n"
            + Log::prefixMessage(Log::MessageType::Bug, "bug from rank 2") + "\n"
            + Log::prefixMessage(Log::MessageType::Bug, "bug from rank 2") + "\n"
            + Log::prefixMessage(Log::MessageType::Bug, "Message limit reached for message tag: tagme") + "\n";
        BOOST_CHECK_EQUAL(log_stream.str(), expected);
    }
}

BOOST_AUTO_TEST_CASE(AllHaveOneMessage)
{
    const auto& cc = Dune::MPIHelper::getCommunication();

    std::ostringstream log_stream;
    initLogger(log_stream);

    Opm::DeferredLogger local_deferredlogger;
    local_deferredlogger.info("info from rank " + std::to_string(cc.rank()));

    Opm::DeferredLogger global_deferredlogger = gatherDeferredLogger(local_deferredlogger, cc);

    if (cc.rank() == 0) {

        global_deferredlogger.logMessages();

        auto counter = OpmLog::getBackend<CounterLog>("COUNTER");
        BOOST_CHECK_EQUAL( cc.size() , counter->numMessages(Log::MessageType::Info) );

        std::string expected;
        for (int i=0; i<cc.size(); i++) {
            expected += Log::prefixMessage(Log::MessageType::Info, "info from rank "+std::to_string(i)) + "\n";
        }
        BOOST_CHECK_EQUAL(log_stream.str(), expected);
    }
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
#if HAVE_MPI
    // register a throwing error handler to allow for
    // debugging with "catch throw" in gdb
    MPI_Errhandler handler;
    MPI_Comm_create_errhandler(MPI_err_handler, &handler);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, handler);
#endif
    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
