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

} // namespace Opm

namespace Opm {

/// Class handling convergence history output for a simulator.
class SimulatorConvergenceOutput
{
public:
    /// Start convergence output thread.
    ///
    /// Thread starts only if explicitly requested in runtime user options.
    ///
    /// \param[in] eclState Static model object.  Needed to determine run's
    /// output directory, base name, and unit conventions.
    ///
    /// \param[in] convOutputOptions Comma separated option string as
    /// required by class ConvergenceOutputConfiguration.
    ///
    /// \param[in] optionName Name of command line option whose value is \p
    /// convOutputOptions.  Used as diagnostic information only.
    ///
    /// \param[in] getPhaseName Callable object for converting component
    /// indices into human readable component names.
    void startThread(const EclipseState& eclState,
                     std::string_view    convOutputOptions,
                     std::string_view    optionName,
                     ConvergenceOutputThread::ComponentToPhaseName getPhaseName);

    /// Create convergence output requests.
    ///
    /// Must not be called before startThread() or after endThread().
    ///
    /// Only those reports which have not previously been emitted will be
    /// written, and only if the run actually requests convergence output at
    /// the non-linear iteration level.
    ///
    /// \param[in] reports All step reports generated in the simulation run
    /// so far.  Class SimulatorConvergenceOutput maintains a record of
    /// which reports have been written.
    void write(const std::vector<StepReport>& reports);

    /// Request that convergence output thread be shut down.
    ///
    /// No additional output requests should be submitted to write() once
    /// this function is called.
    void endThread();

private:
    /// Number of step reports which have already been written.
    std::vector<StepReport>::size_type alreadyReportedSteps_ = 0;

    /// Communication channel between main thread creating output requests
    /// and output thread actually writing those requests to file on disk.
    ///
    /// Nullopt unless convergence output has been requested at the
    /// non-linear iteration level.
    std::optional<ConvergenceReportQueue> convergenceOutputQueue_{};

    /// Output request controller object.
    ///
    /// Nullopt unless convergence output has been requested at the
    /// non-linear iteration level.
    std::optional<ConvergenceOutputThread> convergenceOutputObject_{};

    /// Output request thread.
    ///
    /// Calls output writing member functions on the
    /// convergenceOutputObject_.
    ///
    /// Nullopt unless convergence output has been requested at the
    /// non-linear iteration level.
    std::optional<std::thread> convergenceOutputThread_{};
};

} // namespace Opm

#endif // OPM_SIMULATOR_CONVERGENCE_OUTPUT_HEADER_INCLUDED
