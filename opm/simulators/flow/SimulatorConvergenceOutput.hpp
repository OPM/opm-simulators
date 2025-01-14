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

#ifndef OPM_SIMULATOR_CONVERGENCE_OUTPUT_HEADER_INCLUDED
#define OPM_SIMULATOR_CONVERGENCE_OUTPUT_HEADER_INCLUDED

#include <opm/simulators/flow/ExtraConvergenceOutputThread.hpp>

#include <cstddef>
#include <optional>
#include <string_view>
#include <thread>
#include <vector>

namespace Opm {

class EclipseState;
struct StepReport;

/// Class handling convergence history output for a simulator.
class SimulatorConvergenceOutput
{
public:
    explicit SimulatorConvergenceOutput(const EclipseState& eclState)
        : eclState_(eclState)
    {}

    void startThread(std::string_view convOutputOptions,
                     std::string_view optionName,
                     ConvergenceOutputThread::ComponentToPhaseName getPhaseName);

    void write(const std::vector<StepReport>& reports);

    void endThread();

private:
    const EclipseState& eclState_;

    std::size_t already_reported_steps_ = 0;

    std::optional<ConvergenceReportQueue> convergenceOutputQueue_{};
    std::optional<ConvergenceOutputThread> convergenceOutputObject_{};
    std::optional<std::thread> convergenceOutputThread_{};
};

} // namespace Opm

#endif // OPM_SIMULATOR_CONVERGENCE_OUTPUT_HEADER_INCLUDED
