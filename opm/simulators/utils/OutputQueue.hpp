/*
  Copyright 2022 Equinor ASA.

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

#ifndef OUTPUT_QUEUE_HPP
#define OUTPUT_QUEUE_HPP

#include <condition_variable>
#include <functional>
#include <mutex>
#include <vector>

namespace Opm {
/// Communication channel between thread creating output requests and
/// consumer thread writing those requests to a file.
///
/// Output thread has access to internal state.  Producer thread uses public
/// interface.  Producer thread creates an object of this type and launches
/// the output thread with a reference to that queue object.
template<class RequestType>
class OutputQueue
{
public:
    class Thread
    {
    public:
        Thread(OutputQueue& queue)
            : m_queue(queue)
        {}

        virtual void write(const std::vector<RequestType>& requests) = 0;
        virtual bool finalRequestWritten() const = 0;

        void writeASynchronous()
        {
            // This is the main function of the convergence output thread.  It runs
            // for the duration of the process, although mostly in an idle/waiting
            // state.  Implementation from Microsoft's "GoingNative" show, episode
            // 53, on threading and parallelism in the STL.

            // Note: Loop terminates only when final request written.
            for (auto localReq = std::vector<RequestType>{} ; ; localReq.clear()) {
                std::unique_lock<std::mutex> lock { m_queue.mtx_ };
                m_queue.cv_.wait(lock, [this]() { return !m_queue.requests_.empty(); });

                // Capture all pending output requests, relinquish lock and write
                // file output outside of the critical section.
                m_queue.requests_.swap(localReq);

                lock.unlock();

                this->write(localReq);

                if (this->finalRequestWritten()) {
                    // No more output coming.  Shut down thread.
                    return;
                }
            }
        }

    private:
        OutputQueue& m_queue;
    };

    /// Push sequence of output requests, typically all substeps whether
    /// converged or not, of a single report step.
    ///
    /// \param[in] requests Output request sequence.  Queue takes ownership.
    void enqueue(std::vector<RequestType>&& requests)
    {
        // Signal output thread if we're going from "no work" to "some work".
        // We don't need to signal if we're going from "some work" to "more
        // work".
        auto must_notify = false;

        {
          std::lock_guard<std::mutex> guard{ this->mtx_ };
          must_notify = this->requests_.empty();
          this->requests_.insert(this->requests_.end(),
                                 std::make_move_iterator(requests.begin()),
                                 std::make_move_iterator(requests.end()));
        }

        if (must_notify) {
            this->cv_.notify_one();
        }
    }

    /// Signal end of output request stream.
    ///
    /// No additional requests should be added to queue following a call to
    /// this member function.  Output thread detects this signal, completes
    /// any pending output requests, and shuts down afterwards.
    void signalLastOutputRequest()
    {
        // Empty request signals end of production.
        this->enqueue(std::vector<RequestType>(1));
    }

private:
    /// Mutex for critical sections protecting 'requests_'.
    std::mutex mtx_{};

    /// Condition variable for threads waiting on changes to 'requests_'.
    std::condition_variable cv_{};

    /// Pending convergence output requests.
    std::vector<RequestType> requests_{};
};

} // namespace Opm

#endif // OUTPUT_QUEUE_HPP
