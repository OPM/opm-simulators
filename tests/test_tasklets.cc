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
 *        asynchronous work.
 */
#include "config.h"

#include <opm/models/parallel/tasklets.hh>

#include <chrono>
#include <iostream>

std::mutex outputMutex;

Opm::TaskletRunner *runner;

class SleepTasklet : public Opm::TaskletInterface
{
public:
    SleepTasklet(int mseconds)
        : mseconds_(mseconds)
    {
        n_ = numInstantiated_;
        ++ numInstantiated_;
    }

    void run()
    {
        assert(0 <= runner->workerThreadIndex() && runner->workerThreadIndex() < runner->numWorkerThreads());
        std::this_thread::sleep_for(std::chrono::milliseconds(mseconds_));
        outputMutex.lock();
        std::cout << "Sleep tasklet " << n_ << " of " << mseconds_ << " ms completed by worker thread " << runner->workerThreadIndex() << std::endl;
        outputMutex.unlock();
    }

private:
    static int numInstantiated_;
    int n_;
    int mseconds_;
};

void sleepAndPrintFunction();
void sleepAndPrintFunction()
{
    int ms = 100;
    std::this_thread::sleep_for(std::chrono::milliseconds(ms));
    outputMutex.lock();
    std::cout << "Sleep completed by worker thread " << runner->workerThreadIndex() << std::endl;
    outputMutex.unlock();
}

int SleepTasklet::numInstantiated_ = 0;

int main()
{
    int numWorkers = 2;
    runner = new Opm::TaskletRunner(numWorkers);

    // the master thread is not a worker thread
    assert(runner->workerThreadIndex() < 0);
    assert(runner->numWorkerThreads() == numWorkers);

    for (int i = 0; i < 5; ++ i) {
        //auto st = std::make_shared<SleepTasklet>((i + 1)*1000);
        auto st = std::make_shared<SleepTasklet>(100);
        runner->dispatch(st);
    }

    std::cout << "before barrier" << std::endl;
    runner->barrier();
    std::cout << "after barrier" << std::endl;

    runner->dispatchFunction(sleepAndPrintFunction);
    runner->dispatchFunction(sleepAndPrintFunction, /*numInvokations=*/6);

    delete runner;

    return 0;
}

