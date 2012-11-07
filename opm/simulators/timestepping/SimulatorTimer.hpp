/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_SIMULATORTIMER_HEADER_INCLUDED
#define OPM_SIMULATORTIMER_HEADER_INCLUDED

#include <iosfwd>
#include <vector>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>

namespace Opm
{

    namespace parameter { class ParameterGroup; }
    class EclipseGridParser;


    class SimulatorTimer
    {
    public:
        /// Default constructor.
        SimulatorTimer();

        /// Initialize from parameters. Accepts the following:
        ///    num_psteps    (default 1)
        ///    stepsize_days (default 1)
        void init(const parameter::ParameterGroup& param);

        /// Initialize from TSTEP field.
        void init(const EclipseGridParser& deck);

        /// Total number of steps.
        int numSteps() const;

        /// Current step number.
        int currentStepNum() const;

        /// Set current step number.
        void setCurrentStepNum(int step);

        /// Current step length.
        /// Note: if done(), it is an error to call currentStepLength().
        double currentStepLength() const;

        /// Current time.
        double currentTime() const;

        boost::posix_time::ptime currentDateTime() const;

        /// Total time.
        double totalTime() const;

        /// Set total time.
        /// This is primarily intended for multi-epoch schedules,
        /// where a timer for a given epoch does not have
        /// access to later timesteps.
        void setTotalTime(double time);

        /// Print a report with current and total time etc.
        /// Note: if done(), it is an error to call report().
        void report(std::ostream& os) const;

        /// Next step.
        SimulatorTimer& operator++();

        /// Return true if op++() has been called numSteps() times.
        bool done() const;

    private:
        std::vector<double> timesteps_;
        int current_step_;
        double current_time_;
        double total_time_;
        boost::gregorian::date start_date_;
    };


} // namespace Opm

#endif // OPM_SIMULATORTIMER_HEADER_INCLUDED
