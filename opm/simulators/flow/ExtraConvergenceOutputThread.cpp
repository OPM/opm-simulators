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

#include <opm/simulators/flow/ExtraConvergenceOutputThread.hpp>

#include <opm/simulators/flow/ConvergenceOutputConfiguration.hpp>

#include <opm/simulators/timestepping/ConvergenceReport.hpp>

#include <algorithm>
#include <cassert>
#include <condition_variable>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iterator>
#include <mutex>
#include <numeric>
#include <optional>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {
    std::string to_string(const Opm::ConvergenceReport::ReservoirFailure::Type t)
    {
        using Type = Opm::ConvergenceReport::ReservoirFailure::Type;

        const auto type_strings = std::unordered_map<Type, std::string> {
            { Type::Invalid    , "Invalid" },
            { Type::MassBalance, "MB"      },
            { Type::Cnv        , "CNV"     },
        };

        auto strPos = type_strings.find(t);
        assert ((strPos != type_strings.end()) &&
                "Unsupported convergence metric type");

        return strPos->second;
    }

    std::string
    formatMetricColumn(const Opm::ConvergenceOutputThread::ComponentToPhaseName& getPhaseName,
                       const Opm::ConvergenceReport::ReservoirConvergenceMetric& metric)
    {
        std::ostringstream os;
        os << to_string(metric.type()) << '_' << getPhaseName(metric.phase());

        return os.str();
    }

    std::string::size_type
    maxColHeaderSize(const std::string::size_type                                           minColSize,
                     const Opm::ConvergenceOutputThread::ComponentToPhaseName&              getPhaseName,
                     const std::vector<Opm::ConvergenceReport::ReservoirConvergenceMetric>& cols)
    {
        return std::accumulate(cols.begin(), cols.end(), minColSize,
            [&getPhaseName](const std::string::size_type                              maxChar,
                            const Opm::ConvergenceReport::ReservoirConvergenceMetric& metric)
        {
            return std::max(maxChar, formatMetricColumn(getPhaseName, metric).size());
        });
    }

    std::string::size_type
    writeConvergenceHeader(std::ostream&                                             os,
                           const Opm::ConvergenceOutputThread::ComponentToPhaseName& getPhaseName,
                           const Opm::ConvergenceReportQueue::OutputRequest&         firstRequest)
    {
        const auto minColSize = std::string::size_type{11};

        os << std::right << std::setw(minColSize) << "ReportStep" << ' '
           << std::right << std::setw(minColSize) << "TimeStep"   << ' '
           << std::right << std::setw(minColSize) << "Time"       << ' '
           << std::right << std::setw(minColSize) << "Iteration";

        const auto& metrics = firstRequest.reports.front().reservoirConvergence();
        const auto  maxChar = maxColHeaderSize(minColSize, getPhaseName, metrics);

        for (const auto& metric : metrics) {
            os << std::right << std::setw(maxChar + 1)
               << formatMetricColumn(getPhaseName, metric);
        }

        // Note: Newline character intentionally placed in separate output
        // request to not influence right-justification of column header.
        os << std::right << std::setw(maxChar + 1) << "WellStatus" << '\n';

        return maxChar;
    }


    std::string to_string(const Opm::ConvergenceReport::Severity s)
    {
        using S = Opm::ConvergenceReport::Severity;
        switch (s) {
        case S::None:
            return "None";
        case S::Normal:
            return "Normal";
        case S::TooLarge:
            return "TooLarge";
        case S::NotANumber:
            return "NotANumber";
        }
        throw std::logic_error("Unknown ConvergenceReport::Severity");
    }


    std::string to_string(const Opm::ConvergenceReport::WellFailure::Type t)
    {
        using T = Opm::ConvergenceReport::WellFailure::Type;
        switch (t) {
        case T::Invalid:
            return "Invalid";
        case T::MassBalance:
            return "MassBalance";
        case T::Pressure:
            return "Pressure";
        case T::ControlBHP:
            return "ControlBHP";
        case T::ControlTHP:
            return "ControlTHP";
        case T::ControlRate:
            return "ControlRate";
        case T::Unsolvable:
            return "Unsolvable";
        case T::WrongFlowDirection:
            return "WrongFlowDirection";
        }
        throw std::logic_error("Unknown ConvergenceReport::WellFailure::Type");
    }


    void writeConvergenceRequest(std::ostream&                                           os,
                                 const Opm::ConvergenceOutputThread::ConvertToTimeUnits& convertTime,
                                 std::string::size_type                                  colSize,
                                 const Opm::ConvergenceReportQueue::OutputRequest&       request)
    {
        const auto firstColSize = std::string::size_type{11};

        os.setf(std::ios_base::scientific);

        auto iter = 0;
        for (const auto& report : request.reports) {
            os << std::setw(firstColSize)          << request.reportStep  << ' '
               << std::setw(firstColSize)          << request.currentStep << ' '
               << std::setprecision(4)             << std::setw(firstColSize)
               << convertTime(report.reportTime()) << ' '
               << std::setw(firstColSize)          << iter;

            for (const auto& metric : report.reservoirConvergence()) {
                os << std::setprecision(4) << std::setw(colSize + 1) << metric.value();
            }

            os << std::right << std::setw(colSize + 1)
               << (report.wellFailed() ? "FAIL" : "CONV");
            if (report.wellFailed()) {
                for (const auto& wf : report.wellFailures()) {
                    os << " { "
                       << wf.wellName() << ' ' << to_string(wf.type());
                    if (wf.type() == Opm::ConvergenceReport::WellFailure::Type::MassBalance) {
                        os << " Severity=" << to_string(wf.severity())
                           << " Phase=" << wf.phase();
                    }
                    os << " }";
                }
            }
            os << '\n';

            ++iter;
        }
    }
} // Anonymous namespace

// ---------------------------------------------------------------------------

class Opm::ConvergenceOutputThread::Impl
{
public:
    explicit Impl(std::string_view               outputDir,
                  std::string_view               baseName,
                  ComponentToPhaseName           getPhaseName,
                  ConvertToTimeUnits             convertTime,
                  ConvergenceOutputConfiguration config,
                  ConvergenceReportQueue&        queue);

    ConvergenceReportQueue& queue()
    {
        return this->queue_;
    }

    void write(const std::vector<ConvergenceReportQueue::OutputRequest>& requests);
    bool finalRequestWritten() const
    {
        return this->finalRequestWritten_;
    }

private:
    std::reference_wrapper<ConvergenceReportQueue> queue_;
    ComponentToPhaseName getPhaseName_{};
    ConvertToTimeUnits convertTime_{};
    std::optional<std::ofstream> infoIter_{};
    std::string::size_type colSize_{0};
    bool haveOutputIterHeader_{false};
    bool finalRequestWritten_{false};

    void writeIterInfo(const std::vector<ConvergenceReportQueue::OutputRequest>& requests);
};

Opm::ConvergenceOutputThread::Impl::Impl(std::string_view               outputDir,
                                         std::string_view               baseName,
                                         ComponentToPhaseName           getPhaseName,
                                         ConvertToTimeUnits             convertTime,
                                         ConvergenceOutputConfiguration config,
                                         ConvergenceReportQueue&        queue)
    : queue_        { std::ref(queue) }
    , getPhaseName_ { std::move(getPhaseName) }
    , convertTime_  { std::move(convertTime) }
{
    if (config.want(ConvergenceOutputConfiguration::Option::Iterations)) {
        this->infoIter_.emplace
            (std::filesystem::path { outputDir } /
             std::filesystem::path { baseName }.concat(".INFOITER"));
    }
}

void
Opm::ConvergenceOutputThread::Impl::
write(const std::vector<ConvergenceReportQueue::OutputRequest>& requests)
{
    assert (! requests.empty() &&
            "Internal logic error in forming convergence output request");

    this->writeIterInfo(requests);
}

void
Opm::ConvergenceOutputThread::Impl::
writeIterInfo(const std::vector<ConvergenceReportQueue::OutputRequest>& requests)
{
    if (! this->infoIter_.has_value()) {
        return;
    }

    if (! this->haveOutputIterHeader_) {
        this->colSize_ =
            writeConvergenceHeader(this->infoIter_.value(),
                                   this->getPhaseName_,
                                   requests.front());
        this->haveOutputIterHeader_ = true;
    }

    for (const auto& request : requests) {
        writeConvergenceRequest(this->infoIter_.value(),
                                this->convertTime_,
                                this->colSize_,
                                request);

        if (request.reports.empty()) {
            this->finalRequestWritten_ = true;
            break;
        }
    }

    this->infoIter_.value().flush();
}

// ===========================================================================
// Public Interface Below Separator
// ===========================================================================

void Opm::ConvergenceReportQueue::enqueue(std::vector<OutputRequest>&& requests)
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

void Opm::ConvergenceReportQueue::signalLastOutputRequest()
{
    // Empty request signals end of production.
    this->enqueue(std::vector<OutputRequest>(1));
}

// ---------------------------------------------------------------------------

Opm::ConvergenceOutputThread::
ConvergenceOutputThread(std::string_view               outputDir,
                        std::string_view               baseName,
                        ComponentToPhaseName           getPhaseName,
                        ConvertToTimeUnits             convertTime,
                        ConvergenceOutputConfiguration config,
                        ConvergenceReportQueue&        queue)
    : pImpl_ { std::make_unique<Impl>(outputDir,
                                      baseName,
                                      getPhaseName,
                                      convertTime,
                                      config,
                                      queue) }
{}

Opm::ConvergenceOutputThread::~ConvergenceOutputThread() = default;

Opm::ConvergenceOutputThread::ConvergenceOutputThread(ConvergenceOutputThread&& src)
    : pImpl_ { std::move(src.pImpl_) }
{}

void
Opm::ConvergenceOutputThread::
writeSynchronous(std::vector<ConvergenceReportQueue::OutputRequest>&& requests)
{
    this->pImpl_->write(requests);
}

void Opm::ConvergenceOutputThread::writeASynchronous()
{
    // This is the main function of the convergence output thread.  It runs
    // for the duration of the process, although mostly in an idle/waiting
    // state.  Implementation from Microsoft's "GoingNative" show, episode
    // 53, on threading and parallelism in the STL.

    auto& queue = this->pImpl_->queue();

    // Note: Loop terminates only when final request written.
    for (auto localReq = std::vector<ConvergenceReportQueue::OutputRequest>{} ; ; localReq.clear()) {
        std::unique_lock<std::mutex> lock { queue.mtx_ };
        queue.cv_.wait(lock, [&queue]() { return !queue.requests_.empty(); });

        // Capture all pending output requests, relinquish lock and write
        // file output outside of the critical section.
        queue.requests_.swap(localReq);

        lock.unlock();

        this->pImpl_->write(localReq);

        if (this->pImpl_->finalRequestWritten()) {
            // No more output coming.  Shut down thread.
            return;
        }
    }
}
