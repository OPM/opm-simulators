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

#include "config.h"
#include <opm/simulators/timestepping/SimulatorTimer.hpp>
#include <opm/common/utility/parameters/ParameterGroup.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Units/Units.hpp>

#include <boost/date_time/gregorian/gregorian_types.hpp>
#include <boost/date_time/posix_time/conversion.hpp>

#include <cstddef>
#include <numeric>
#include <ostream>

namespace Opm
{

    /// Default constructor.
    SimulatorTimer::SimulatorTimer()
        : current_step_(0),
          current_time_(0.0),
          start_date_(2012,1,1)    // A really arbitrary default starting value?!
    {
    }

    SimulatorTimer SimulatorTimer::serializationTestObject()
    {
        SimulatorTimer res;
        res.timesteps_ = {1.0, 2.0, 3.0};
        res.current_step_ = 4;
        res.current_time_ = 5.0;
        res.total_time_ = 6.0;
        res.start_date_ = boost::gregorian::date(2023,1,31);

        return res;
    }

    /// Initialize from parameters. Accepts the following:
    ///    num_psteps    (default 1)
    ///    stepsize_days (default 1)
    void SimulatorTimer::init(const ParameterGroup& param)
    {
        const int num_psteps = param.getDefault("num_psteps", 1);
        const double stepsize_days = param.getDefault("stepsize_days", 1.0);
        const double stepsize = Opm::unit::convert::from(stepsize_days, Opm::unit::day);
        timesteps_.clear();
        timesteps_.resize(num_psteps, stepsize);
        total_time_ = num_psteps*stepsize;
    }

    /// Use the SimulatorTimer as a shim around opm-parser's Opm::TimeMap
    void SimulatorTimer::init(const Schedule& schedule, std::size_t report_step)
    {
        total_time_ = schedule.seconds( schedule.size() - 1 );
        timesteps_.resize(schedule.size() - 1);
        for (std::size_t i = 0; i < schedule.size() - 1; ++i) {
            timesteps_[i] = schedule.stepLength(i);
        }

        setCurrentStepNum(report_step);
        start_date_ = boost::posix_time::from_time_t(schedule.getStartTime()).date();
    }

    /// Whether the current step is the first step.
    bool SimulatorTimer::initialStep() const
    {
        return (current_step_ == 0);
    }

    /// Total number of steps.
    int SimulatorTimer::numSteps() const
    {
        return timesteps_.size();
    }

    /// Current step number.
    int SimulatorTimer::currentStepNum() const
    {
        return current_step_;
    }

    /// Set current step number.
    void SimulatorTimer::setCurrentStepNum(int step)
    {
        current_step_ = step;
        current_time_ = std::accumulate(timesteps_.begin(), timesteps_.begin() + step, 0.0);
    }


    /// Current step length.
    double SimulatorTimer::currentStepLength() const
    {
        assert(!done());
        return timesteps_[current_step_];
    }

    double SimulatorTimer::stepLengthTaken() const
    {
        assert(current_step_ > 0);
        return timesteps_[current_step_ - 1];
    }

    /// time elapsed since the start of the simulation [s].
    double SimulatorTimer::simulationTimeElapsed() const
    {
        return current_time_;
    }

    boost::posix_time::ptime SimulatorTimer::startDateTime() const
    {
        return boost::posix_time::ptime(start_date_);
    }


    /// Total time.
    double SimulatorTimer::totalTime() const
    {
        return total_time_;
    }

    /// Set total time.
    /// This is primarily intended for multi-epoch schedules,
    /// where a timer for a given epoch does not have
    /// access to later timesteps.
    void SimulatorTimer::setTotalTime(double time)
    {
        total_time_ = time;
    }

    /// Print a report with current and total time etc.
    void SimulatorTimer::report(std::ostream& os) const
    {
        os << "\n\n---------------    Simulation step number " << currentStepNum() << "    ---------------"
           << "\n      Current time (days)     " << Opm::unit::convert::to(simulationTimeElapsed(), Opm::unit::day)
           << "\n      Current stepsize (days) " << Opm::unit::convert::to(currentStepLength(), Opm::unit::day)
           << "\n      Total time (days)       " << Opm::unit::convert::to(totalTime(), Opm::unit::day)
           << "\n" << std::endl;
    }

    /// Next step.
    SimulatorTimer& SimulatorTimer::operator++()
    {
        assert(!done());
        current_time_ += timesteps_[current_step_];
        ++current_step_;
        return *this;
    }

    /// Return true if op++() has been called numSteps() times.
    bool SimulatorTimer::done() const
    {
        return int(timesteps_.size()) == current_step_;
    }

    /// return copy of object
    std::unique_ptr< SimulatorTimerInterface >
    SimulatorTimer::clone() const
    {
       return std::make_unique<SimulatorTimer>(*this);
    }

    bool SimulatorTimer::operator==(const SimulatorTimer& rhs) const
    {
        return this->timesteps_ == rhs.timesteps_ &&
               this->current_step_ == rhs.current_step_ &&
               this->current_time_ == rhs.current_time_ &&
               this->total_time_ == rhs.total_time_ &&
               this->start_date_ == rhs.start_date_;
    }

} // namespace Opm
