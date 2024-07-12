// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \brief This file serves as an example of how to use the tasklet mechanism for
 *        asynchronous work, especially for tasklets that fail.
 */
#define BOOST_TEST_MODULE TASKLETS_FAILURE

//#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>
#include <chrono>
#include <iostream>
#include <mutex>
#include <thread>
#include <sys/wait.h>
#include <unistd.h>
#include <cassert>

#include "config.h"
#include <opm/models/parallel/tasklets.hh>

std::mutex outputMutex;

Opm::TaskletRunner *runner;

class SleepTasklet : public Opm::TaskletInterface
{
public:
    SleepTasklet(int mseconds, int id)
        : mseconds_(mseconds),
          id_(id)
    {}

    void run() override
    {
        assert(0 <= runner->workerThreadIndex() && runner->workerThreadIndex() < runner->numWorkerThreads());
        std::cout << "Sleep tasklet " << id_ << " of " << mseconds_ << " ms starting sleep on worker thread " << runner->workerThreadIndex() << std::endl;
        std::this_thread::sleep_for(std::chrono::milliseconds(mseconds_));
        outputMutex.lock();
        std::cout << "Sleep tasklet " << id_ << " of " << mseconds_ << " ms completed by worker thread " << runner->workerThreadIndex() << std::endl;
        outputMutex.unlock();
    }

private:
    int mseconds_;
    int id_;
};

class FailingSleepTasklet : public Opm::TaskletInterface
{
public:
    FailingSleepTasklet(int mseconds)
        : mseconds_(mseconds)
    {}
    void run() override
    {
        std::this_thread::sleep_for(std::chrono::milliseconds(mseconds_));
        outputMutex.lock();
        std::cout << "Failing sleep tasklet of " << mseconds_ << " ms failing now, on work thread " << runner->workerThreadIndex() << std::endl;
        outputMutex.unlock();
        throw std::logic_error("Intentional failure for testing");
    }

private:
    int mseconds_;
};

void execute () {
        int numWorkers = 2;
        runner = new Opm::TaskletRunner(numWorkers);

        // the master thread is not a worker thread
        assert(runner->workerThreadIndex() < 0);
        assert(runner->numWorkerThreads() == numWorkers);

        // Dispatch some successful tasklets
        for (int i = 0; i < 5; ++i) {
            runner->barrier();

            if (runner->failure()) {
                exit(EXIT_FAILURE);
            }
            auto st = std::make_shared<SleepTasklet>(10,i);
            runner->dispatch(st);
        }

        runner->barrier();
        if (runner->failure()) {
            exit(EXIT_FAILURE);
        }
        // Dispatch a failing tasklet
        auto failingSleepTasklet = std::make_shared<FailingSleepTasklet>(100);
        runner->dispatch(failingSleepTasklet);

        // Dispatch more successful tasklets
        for (int i = 5; i < 10; ++i) {
            runner->barrier();

            if (runner->failure()) {
                exit(EXIT_FAILURE);
            }
            auto st = std::make_shared<SleepTasklet>(10,i);
            runner->dispatch(st);
        }

        std::cout << "before barrier" << std::endl;
        runner->barrier();
    }
BOOST_AUTO_TEST_SUITE(Tasklets)
BOOST_AUTO_TEST_CASE(TASKLETS_FAILURE) {
    pid_t pid = fork(); // Create a new process, such that this child process can call exit(EXIT_FAILURE)
    if (pid == -1) {
        BOOST_FAIL("Fork failed");
    } else if (pid == 0) {
        // Child process
        execute();
        _exit(0);  // Should never reach here
    } else {
        // Parent process
        std::cout << "Checking failure of child process with parent process, process id " << pid << std::endl;
        int status;
        waitpid(pid, &status, 0);
        BOOST_CHECK(WIFEXITED(status));  // Check if the child process exited
        BOOST_CHECK_EQUAL(WEXITSTATUS(status), EXIT_FAILURE);  // Check if the exit status is EXIT_FAILURE
    }
}
BOOST_AUTO_TEST_SUITE_END()
