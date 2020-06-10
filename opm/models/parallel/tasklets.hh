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
 * \brief Provides a mechanism to dispatch work to separate threads
 */
#ifndef EWOMS_TASKLETS_HH
#define EWOMS_TASKLETS_HH

#include <stdexcept>
#include <cassert>
#include <thread>
#include <queue>
#include <mutex>
#include <iostream>
#include <condition_variable>

namespace Opm {

/*!
 * \brief The base class for tasklets.
 *
 * Tasklets are a generic mechanism for potentially running work in a separate thread.
 */
class TaskletInterface
{
public:
    TaskletInterface(int refCount = 1)
        : referenceCount_(refCount)
    {}
    virtual ~TaskletInterface() {}
    virtual void run() = 0;
    virtual bool isEndMarker () const { return false; }

    void dereference()
    { -- referenceCount_; }

    int referenceCount() const
    { return referenceCount_; }

private:
    int referenceCount_;
};

/*!
 * \brief A simple tasklet that runs a function that returns void and does not take any
 *        arguments a given number of times.
 */
template <class Fn>
class FunctionRunnerTasklet : public TaskletInterface
{
public:
    FunctionRunnerTasklet(const FunctionRunnerTasklet&) = default;
    FunctionRunnerTasklet(int numInvocations, const Fn& fn)
        : TaskletInterface(numInvocations)
        , fn_(fn)
    {}
    void run() override
    { fn_(); }

private:
    const Fn& fn_;
};

class TaskletRunner;

// this class stores the thread local static attributes for the TaskletRunner class. we
// cannot put them directly into TaskletRunner because defining static members for
// non-template classes in headers leads the linker to choke in case multiple compile
// units are used.
template <class Dummy = void>
struct TaskletRunnerHelper_
{
    static thread_local TaskletRunner* taskletRunner_;
    static thread_local int workerThreadIndex_;
};

template <class Dummy>
thread_local TaskletRunner* TaskletRunnerHelper_<Dummy>::taskletRunner_ = nullptr;

template <class Dummy>
thread_local int TaskletRunnerHelper_<Dummy>::workerThreadIndex_ = -1;

/*!
 * \brief Handles where a given tasklet is run.
 *
 * Depending on the number of worker threads, a tasklet can either be run in a separate
 * worker thread or by the main thread.
 */
class TaskletRunner
{
    /// \brief Implements a barrier. This class can only be used in the asynchronous case.
    class BarrierTasklet : public TaskletInterface
    {
    public:
        BarrierTasklet(unsigned numWorkers)
            : TaskletInterface(/*refCount=*/numWorkers)
        {
            numWorkers_ = numWorkers;
            numWaiting_ = 0;
        }

        void run()
        { wait(); }

        void wait()
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

    private:
        unsigned numWorkers_;
        unsigned numWaiting_;

        std::condition_variable barrierCondition_;
        std::mutex barrierMutex_;
    };

    /// \brief TerminateThreadTasklet class
    /// Empty tasklet marking thread termination.
    class TerminateThreadTasklet : public TaskletInterface
    {
    public:
        void run()
        { }

        bool isEndMarker() const
        { return true; }
    };

public:
    // prohibit copying of tasklet runners
    TaskletRunner(const TaskletRunner&) = delete;

    /*!
     * \brief Creates a tasklet runner with numWorkers underling threads for doing work.
     *
     * The number of worker threads may be 0. In this case, all work is done by the main
     * thread (synchronous mode).
     */
    TaskletRunner(unsigned numWorkers)
    {
        threads_.resize(numWorkers);
        for (unsigned i = 0; i < numWorkers; ++i)
            // create a worker thread
            threads_[i].reset(new std::thread(startWorkerThread_, this, i));
    }

    /*!
     * \brief Destructor
     *
     * If worker threads were created to run the tasklets, this method waits until all
     * worker threads have been terminated, i.e. all scheduled tasklets are guaranteed to
     * be completed.
     */
    ~TaskletRunner()
    {
        if (threads_.size() > 0) {
            // dispatch a tasklet which will terminate the worker thread
            dispatch(std::make_shared<TerminateThreadTasklet>());

            // wait until all worker threads have terminated
            for (auto& thread : threads_)
                thread->join();
        }
    }

    /*!
     * \brief Returns the index of the current worker thread.
     *
     * If the current thread is not a worker thread, -1 is returned.
     */
    int workerThreadIndex() const
    {
        if (TaskletRunnerHelper_<void>::taskletRunner_ != this)
            return -1;
        return TaskletRunnerHelper_<void>::workerThreadIndex_;
    }

    /*!
     * \brief Returns the number of worker threads for the tasklet runner.
     */
    int numWorkerThreads() const
    { return threads_.size(); }

    /*!
     * \brief Add a new tasklet.
     *
     * The tasklet is either run immediately or deferred to a separate thread.
     */
    void dispatch(std::shared_ptr<TaskletInterface> tasklet)
    {
        if (threads_.empty()) {
            // run the tasklet immediately in synchronous mode.
            while (tasklet->referenceCount() > 0) {
                tasklet->dereference();
                try {
                    tasklet->run();
                }
                catch (const std::exception& e) {
                    std::cerr << "ERROR: Uncaught std::exception when running tasklet: " << e.what() << ". Trying to continue.\n";
                }
                catch (...) {
                    std::cerr << "ERROR: Uncaught exception (general type) when running tasklet. Trying to continue.\n";
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

    /*!
     * \brief Convenience method to construct a new function runner tasklet and dispatch it immediately.
     */
    template <class Fn>
    std::shared_ptr<FunctionRunnerTasklet<Fn> > dispatchFunction(Fn &fn, int numInvocations=1)
    {
        using Tasklet = FunctionRunnerTasklet<Fn>;
        auto tasklet = std::make_shared<Tasklet>(numInvocations, fn);
        this->dispatch(tasklet);
        return tasklet;
    }

    /*!
     * \brief Make sure that all tasklets have been completed after this method has been called
     */
    void barrier()
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

protected:
    // main function of the worker thread
    static void startWorkerThread_(TaskletRunner* taskletRunner, int workerThreadIndex)
    {
        TaskletRunnerHelper_<void>::taskletRunner_ = taskletRunner;
        TaskletRunnerHelper_<void>::workerThreadIndex_ = workerThreadIndex;

        taskletRunner->run_();
    }

    //! do the work until the queue received an end tasklet
    void run_()
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
                std::cerr << "ERROR: Uncaught std::exception when running tasklet: " << e.what() << ". Trying to continue.\n";
            }
            catch (...) {
                std::cerr << "ERROR: Uncaught exception when running tasklet. Trying to continue.\n";
            }
        }
    }

    std::vector<std::unique_ptr<std::thread> > threads_;
    std::queue<std::shared_ptr<TaskletInterface> > taskletQueue_;
    std::mutex taskletQueueMutex_;
    std::condition_variable workAvailableCondition_;
};

} // end namespace Opm
#endif
