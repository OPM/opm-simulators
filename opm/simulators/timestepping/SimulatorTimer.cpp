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
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <ostream>
#include <numeric>

namespace Opm
{

    /// Default constructor.
    SimulatorTimer::SimulatorTimer()
        : current_step_(0),
          current_time_(0.0),
          start_date_(2012,1,1)    // A really arbitrary default starting value?!
    {
    }

    /// Initialize from parameters. Accepts the following:
    ///    num_psteps    (default 1)
    ///    stepsize_days (default 1)
    void SimulatorTimer::init(const parameter::ParameterGroup& param)
    {
        const int num_psteps = param.getDefault("num_psteps", 1);
        const double stepsize_days = param.getDefault("stepsize_days", 1.0);
        const double stepsize = Opm::unit::convert::from(stepsize_days, Opm::unit::day);
        timesteps_.clear();
        timesteps_.resize(num_psteps, stepsize);
        total_time_ = num_psteps*stepsize;
    }

    /// Initialize from TSTEP field.
    void SimulatorTimer::init(const EclipseGridParser& deck)
    {
        timesteps_  = deck.getTSTEP().tstep_;
        total_time_ = std::accumulate(timesteps_.begin(), timesteps_.end(), 0.0);
        start_date_ = deck.getStartDate();
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
        if (current_step_ < 0 || current_step_ > int(timesteps_.size())) {
            // Note that we do allow current_step_ == timesteps_.size(),
            // that is the done() state.
            THROW("Trying to set invalid step number: " << step);
        }
        current_step_ = step;
        current_time_ = std::accumulate(timesteps_.begin(), timesteps_.begin() + step, 0.0);
    }


    /// Current step length.
    double SimulatorTimer::currentStepLength() const
    {
        ASSERT(!done());
        return timesteps_[current_step_];
    }

    /// Current time.
    double SimulatorTimer::currentTime() const
    {
        return current_time_;
    }


    boost::posix_time::ptime SimulatorTimer::currentDateTime() const
    {
      return boost::posix_time::ptime(start_date_) + boost::posix_time::seconds( (int) current_time_ );
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
           << "\n      Current time (days)     " << Opm::unit::convert::to(currentTime(), Opm::unit::day)
           << "\n      Current stepsize (days) " << Opm::unit::convert::to(currentStepLength(), Opm::unit::day)
           << "\n      Total time (days)       " << Opm::unit::convert::to(totalTime(), Opm::unit::day)
           << "\n" << std::endl;
    }

    /// Next step.
    SimulatorTimer& SimulatorTimer::operator++()
    {
        ASSERT(!done());
        current_time_ += timesteps_[current_step_];
        ++current_step_;
        return *this;
    }

    /// Return true if op++() has been called numSteps() times.
    bool SimulatorTimer::done() const
    {
        return int(timesteps_.size()) == current_step_;
    }


} // namespace Opm
