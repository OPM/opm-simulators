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
#ifndef OPM_TASKLETS_HPP
#define OPM_TASKLETS_HPP

#include <atomic>
#include <condition_variable>
#include <memory>
#include <mutex>
#include <queue>
#include <thread>

namespace Opm {

/*!
 * \brief The base class for tasklets.
 *
 * Tasklets are a generic mechanism for potentially running work in a separate thread.
 */
class TaskletInterface
{
public:
    explicit TaskletInterface(int refCount = 1)
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
        explicit BarrierTasklet(unsigned numWorkers);

        void run() override;

        void wait();

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
        void run() override
        {}

        bool isEndMarker() const override
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
    explicit TaskletRunner(unsigned numWorkers);

    /*!
     * \brief Destructor
     *
     * If worker threads were created to run the tasklets, this method waits until all
     * worker threads have been terminated, i.e. all scheduled tasklets are guaranteed to
     * be completed.
     */
    ~TaskletRunner();

    bool failure() const;

    /*!
     * \brief Returns the index of the current worker thread.
     *
     * If the current thread is not a worker thread, -1 is returned.
     */
    int workerThreadIndex() const;

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
    void dispatch(std::shared_ptr<TaskletInterface> tasklet);

    /*!
     * \brief Convenience method to construct a new function runner tasklet and dispatch it immediately.
     */
    template <class Fn>
    std::shared_ptr<FunctionRunnerTasklet<Fn> > dispatchFunction(Fn &fn, int numInvocations = 1)
    {
        using Tasklet = FunctionRunnerTasklet<Fn>;
        auto tasklet = std::make_shared<Tasklet>(numInvocations, fn);
        this->dispatch(tasklet);
        return tasklet;
    }

    /*!
     * \brief Make sure that all tasklets have been completed after this method has been called
     */
    void barrier();

private:
    // Atomic flag that is set to failure if any of the tasklets run by the TaskletRunner fails.
    // This flag is checked before new tasklets run or get dispatched and in case it is true, the
    // thread execution will be stopped / no new tasklets will be started and the program will abort.
    // To set the flag and load the flag, we use std::memory_order_relaxed.
    // Atomic operations tagged memory_order_relaxed are not synchronization operations; they do not
    // impose an order among concurrent memory accesses. They guarantee atomicity and modification order
    // consistency. This is the right choice for the setting here, since it is enough to broadcast failure
    // before new run or get dispatched.
    std::atomic<bool> failureFlag_ = false;

protected:
    // main function of the worker thread
    static void startWorkerThread_(TaskletRunner* taskletRunner, int workerThreadIndex);

    //! do the work until the queue received an end tasklet
    void run_();

    std::vector<std::unique_ptr<std::thread> > threads_;
    std::queue<std::shared_ptr<TaskletInterface> > taskletQueue_;
    std::mutex taskletQueueMutex_;
    std::condition_variable workAvailableCondition_;

    static thread_local TaskletRunner* taskletRunner_;
    static thread_local int workerThreadIndex_;
};

} // end namespace Opm

#endif // OPM_TASKLETS_HPP
