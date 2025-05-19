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
#include <config.h>
#include <opm/models/parallel/tasklets.hpp>

#include <atomic>
#include <condition_variable>
#include <iostream>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <thread>

namespace Opm {

thread_local TaskletRunner* TaskletRunner::taskletRunner_ = nullptr;
thread_local int TaskletRunner::workerThreadIndex_ = -1;

TaskletRunner::BarrierTasklet::BarrierTasklet(unsigned numWorkers)
    : TaskletInterface(/*refCount=*/numWorkers)
{
    numWorkers_ = numWorkers;
    numWaiting_ = 0;
}

void TaskletRunner::BarrierTasklet::run()
{
    wait();
}

void TaskletRunner::BarrierTasklet::wait()
{
    std::unique_lock<std::mutex> lock(barrierMutex_);

    numWaiting_ += 1;
    if (numWaiting_ >= numWorkers_ + 1) {
        lock.unlock();
        barrierCondition_.notify_all();
    }
    else {
        const auto& areAllWaiting =
            [this]() -> bool
            { return this->numWaiting_ >= this->numWorkers_ + 1; };

        barrierCondition_.wait(lock, /*predicate=*/areAllWaiting);
    }
}

TaskletRunner::TaskletRunner(unsigned numWorkers)
{
    threads_.resize(numWorkers);
    for (unsigned i = 0; i < numWorkers; ++i)
        // create a worker thread
        threads_[i].reset(new std::thread(startWorkerThread_, this, i));
}

TaskletRunner::~TaskletRunner()
{
    if (threads_.size() > 0) {
        // dispatch a tasklet which will terminate the worker thread
        dispatch(std::make_shared<TerminateThreadTasklet>());

        // wait until all worker threads have terminated
        for (auto& thread : threads_)
            thread->join();
    }
}

bool TaskletRunner::failure() const
{
    return this->failureFlag_.load(std::memory_order_relaxed);
}

int TaskletRunner::workerThreadIndex() const
{
    if (TaskletRunner::taskletRunner_ != this)
        return -1;
    return TaskletRunner::workerThreadIndex_;
}

void TaskletRunner::dispatch(std::shared_ptr<TaskletInterface> tasklet)
{
    if (threads_.empty()) {
        // run the tasklet immediately in synchronous mode.
        while (tasklet->referenceCount() > 0) {
            tasklet->dereference();
            try {
                tasklet->run();
            }
            catch (const std::exception& e) {
                std::cerr << "ERROR: Uncaught std::exception when running tasklet: " << e.what()
                          << ". Trying to continue.\n";
                failureFlag_.store(true, std::memory_order_relaxed);
            }
            catch (...) {
                std::cerr << "ERROR: Uncaught exception (general type) when running tasklet. Trying to continue.\n";
                failureFlag_.store(true, std::memory_order_relaxed);
            }
        }
    }
    else {
        // lock mutex for the tasklet queue to make sure that nobody messes with the
        // task queue
        taskletQueueMutex_.lock();

        // add the tasklet to the queue
        taskletQueue_.push(tasklet);

        taskletQueueMutex_.unlock();

        workAvailableCondition_.notify_all();
    }
}

void TaskletRunner::barrier()
{
    unsigned numWorkers = threads_.size();
    if (numWorkers == 0)
        // nothing needs to be done to implement a barrier in synchronous mode
        return;

    // dispatch a barrier tasklet and wait until it has been run by the worker thread
    auto barrierTasklet = std::make_shared<BarrierTasklet>(numWorkers);
    dispatch(barrierTasklet);

    barrierTasklet->wait();
}

void TaskletRunner::startWorkerThread_(TaskletRunner* taskletRunner, int workerThreadIndex)
{
    TaskletRunner::taskletRunner_ = taskletRunner;
    TaskletRunner::workerThreadIndex_ = workerThreadIndex;

    taskletRunner->run_();
}

void TaskletRunner::run_()
{
    while (true) {

        // wait until tasklets have been pushed to the queue. first we need to lock
        // mutex for access to taskletQueue_
        std::unique_lock<std::mutex> lock(taskletQueueMutex_);

        const auto& workIsAvailable =
            [this]() -> bool
            { return !taskletQueue_.empty(); };

        if (!workIsAvailable())
            workAvailableCondition_.wait(lock, /*predicate=*/workIsAvailable);

        // remove tasklet from queue
        std::shared_ptr<TaskletInterface> tasklet = taskletQueue_.front();

        // if tasklet is an end marker, terminate the thread and DO NOT remove the
        // tasklet.
        if (tasklet->isEndMarker()) {
            if(taskletQueue_.size() > 1)
                throw std::logic_error("TaskletRunner: Not all queued tasklets were executed");
            taskletQueueMutex_.unlock();
            return;
        }

        tasklet->dereference();
        if (tasklet->referenceCount() == 0)
            // remove tasklets from the queue as soon as their reference count
            // reaches zero, i.e. the tasklet has been run often enough.
            taskletQueue_.pop();
        lock.unlock();

        // execute tasklet
        try {
            tasklet->run();
        }
        catch (const std::exception& e) {
            std::cerr << "ERROR: Uncaught std::exception when running tasklet: " << e.what() << ".\n";
            failureFlag_.store(true, std::memory_order_relaxed);
        }
        catch (...) {
            std::cerr << "ERROR: Uncaught exception when running tasklet.\n";
            failureFlag_.store(true, std::memory_order_relaxed);
        }
    }
}

} // end namespace Opm
