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

#ifndef EXTRA_CONVERGENCE_OUTPUT_THREAD_HPP
#define EXTRA_CONVERGENCE_OUTPUT_THREAD_HPP

#include <opm/simulators/timestepping/ConvergenceReport.hpp>

#include <condition_variable>
#include <functional>
#include <memory>
#include <mutex>
#include <string_view>
#include <vector>

namespace Opm {
    class ConvergenceOutputConfiguration;
} // namespace Opm

/// \file System for outputting additional convergence information, such as
/// material balance and CNV values, at each non-linear iteration.
///
/// Supports an asynchronous protocol that assumes there is a single thread
/// dedicated to per-iteration file output.  Synchronous file output is
/// available for debugging and development purposes.

namespace Opm
{

/// Forward declaration so that Queue may declare this type its 'friend'.
class ConvergenceOutputThread;

/// Communication channel between thread creating output requests and
/// consumer thread writing those requests to a file.
///
/// Output thread has access to internal state.  Producer thread uses public
/// interface.  Producer thread creates an object of this type and launches
/// the output thread with a reference to that queue object.
class ConvergenceReportQueue
{
public:
    /// Single output request.
    ///
    /// Associates non-linear iteration convergence information to single
    /// report and timestep.
    struct OutputRequest
    {
        /// Current report step
        int reportStep{-1};

        /// Current timestep within \c reportStep.  Expected to be a small
        /// integer.
        int currentStep{-1};

        /// Convergence metrics for each non-linear ieration in the \c
        /// currentStep.
        std::vector<ConvergenceReport> reports{};
    };

    /// Push sequence of output requests, typically all substeps whether
    /// converged or not, of a single report step.
    ///
    /// \param[in] requests Output request sequence.  Queue takes ownership.
    void enqueue(std::vector<OutputRequest>&& requests);

    /// Signal end of output request stream.
    ///
    /// No additional requests should be added to queue following a call to
    /// this member function.  Output thread detects this signal, completes
    /// any pending output requests, and shuts down afterwards.
    void signalLastOutputRequest();

    friend class ConvergenceOutputThread;

private:
    /// Mutex for critical sections protecting 'requests_'.
    std::mutex mtx_{};

    /// Condition variable for threads waiting on changes to 'requests_'.
    std::condition_variable cv_{};

    /// Pending convergence output requests.
    std::vector<OutputRequest> requests_{};
};

/// Encapsulating object for thread processing producer's convergence output
/// requests.
class ConvergenceOutputThread
{
public:
    /// Protocol for converting a phase/component ID into a human readable
    /// phase/component name.
    using ComponentToPhaseName = std::function<std::string_view(int)>;

    /// Protocol for converting an SI elapsed time value into an equivalent
    /// time value in the run's output conventions.
    ///
    /// Will typically use \code UnitSystem::from_si() \endcode.
    using ConvertToTimeUnits = std::function<double(double)>;

    /// Constructor.
    ///
    /// \param[in] outputDir    -- Name of run's output directory.  Any file
    ///                            output will be written to this directory.
    ///
    /// \param[in] baseName     -- Run's base name.  Output files will have this
    ///                            name and a type-specific file extension.
    ///
    /// \param[in] getPhaseName -- Callable object for converting component
    ///                            indices into human readable component names.
    ///
    /// \param[in] convertTime  -- Callable object for converting SI elapsed
    ///                            time values into equivalent elapsed time
    ///                            values using run's time conventions.
    ///
    /// \param[in] config       -- Convergence output configuration options.
    ///                            Determines whether to output additional
    ///                            convergence information and, if so, what
    ///                            information to output.
    ///
    /// \param[in] queue        -- Communication channel between producer thread
    ///                            and this output thread.  User must form a
    ///                            valid queue prior to creating the output
    ///                            thread object.
    explicit ConvergenceOutputThread(std::string_view               outputDir,
                                     std::string_view               baseName,
                                     ComponentToPhaseName           getPhaseName,
                                     ConvertToTimeUnits             convertTime,
                                     ConvergenceOutputConfiguration config,
                                     ConvergenceReportQueue&        queue);

    /// Deleted copy constructor.
    ConvergenceOutputThread(const ConvergenceOutputThread& src) = delete;

    /// Move constructor.
    ConvergenceOutputThread(ConvergenceOutputThread&& src);

    /// Deleted assignment operator.
    ConvergenceOutputThread& operator=(const ConvergenceOutputThread& src) = delete;

    /// Deleted move-assignment operator.
    ConvergenceOutputThread& operator=(ConvergenceOutputThread&& src) = delete;

    /// Destructor.
    ///
    /// Needed for pimpl idiom.
    ~ConvergenceOutputThread();

    /// Perform synchronous file output of a sequence of requests.
    ///
    /// Mostly for development and debugging purposes.
    ///
    /// \param[in] requests Output request sequence.  Thread takes ownership.
    void writeSynchronous(std::vector<ConvergenceReportQueue::OutputRequest>&& requests);

    /// Output thread worker function
    ///
    /// This is the endpoint that users should associate to a \code
    /// std::thread \endcode object.
    ///
    /// Returns once last pending output request is written (cf. \code
    /// ConvergenceReportQueue::signalLastOutputRequest() \endcode.)
    void writeASynchronous();

private:
    /// Private implementation class.
    class Impl;

    /// Pointer to implementation.
    std::unique_ptr<Impl> pImpl_;
};

} // namespace Opm

#endif //  EXTRA_CONVERGENCE_OUTPUT_THREAD_HPP
