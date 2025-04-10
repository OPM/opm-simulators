/*
  Copyright 2013, 2015, 2020 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2015 Andreas Lauser
  Copyright 2017 IRIS

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

#include <opm/simulators/flow/SimulatorConvergenceOutput.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/IOConfig/IOConfig.hpp>

#include <opm/simulators/flow/ConvergenceOutputConfiguration.hpp>

#include <algorithm>
#include <iterator>
#include <string_view>
#include <utility>
#include <vector>

namespace Opm {

void SimulatorConvergenceOutput::
startThread(const EclipseState&                           eclState,
            std::string_view                              convOutputOptions,
            std::string_view                              optionName,
            ConvergenceOutputThread::ComponentToPhaseName getPhaseName)
{
    const auto config = ConvergenceOutputConfiguration {
        convOutputOptions, optionName
    };

    if (! config.want(ConvergenceOutputConfiguration::Option::Iterations)) {
        return;
    }

    auto convertTime = ConvergenceOutputThread::ConvertToTimeUnits {
        [usys = eclState.getUnits()](const double time)
        { return usys.from_si(UnitSystem::measure::time, time); }
    };

    this->convergenceOutputQueue_.emplace();
    this->convergenceOutputObject_.emplace
        (eclState.getIOConfig().getOutputDir(),
         eclState.getIOConfig().getBaseName(),
         std::move(getPhaseName),
         std::move(convertTime),
         config, *this->convergenceOutputQueue_);

    this->convergenceOutputThread_
        .emplace(&ConvergenceOutputThread::writeASynchronous,
                 &this->convergenceOutputObject_.value());
}

void SimulatorConvergenceOutput::
write(const std::vector<StepReport>& reports)
{
    if (! this->convergenceOutputThread_.has_value() ||
        (reports.size() == this->alreadyReportedSteps_))
    {
        // Convergence output not requested or we've already written all
        // known reports.  Nothing to do.
        return;
    }

    const auto begin = reports.begin() + this->alreadyReportedSteps_;
    const auto end = reports.end();

    auto requests = std::vector<ConvergenceReportQueue::OutputRequest>{};
    requests.reserve(std::distance(begin, end));

    std::transform(begin, end, std::back_inserter(requests),
                   [](const StepReport& report) {
                       return ConvergenceReportQueue::OutputRequest {
                           report.report_step, report.current_step, report.report
                       };
                   });

    this->alreadyReportedSteps_ = reports.size();

    this->convergenceOutputQueue_->enqueue(std::move(requests));
}

void SimulatorConvergenceOutput::endThread()
{
    if (! this->convergenceOutputThread_.has_value()) {
        return;
    }

    this->convergenceOutputQueue_->signalLastOutputRequest();
    this->convergenceOutputThread_->join();
}

} // namespace Opm
