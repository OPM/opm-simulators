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

    /// Use the SimulatorTimer as a shim around opm-parser's Opm::TimeMap
    void SimulatorTimer::init(Opm::TimeMapConstPtr timeMap,
                              size_t beginReportStepIdx,
                              size_t endReportStepIdx)
    {
        timeMap_ = timeMap;
        current_step_ = 0;
        beginReportStepIdx_ = beginReportStepIdx;
        endReportStepIdx_ = std::min(timeMap_->numTimesteps() + 1, endReportStepIdx);
    }

    /// Total number of steps.
    int SimulatorTimer::numSteps() const
    {
        if (timeMap_)
            return endReportStepIdx_ - beginReportStepIdx_;
        else
            return timesteps_.size();
    }

    /// Index of the first considered simulation episode
    size_t SimulatorTimer::beginReportStepIndex() const
    {
        if (!timeMap_) {
            OPM_THROW(std::runtime_error, "indexFirstEpisode() is only implemented "
                      "for simulation timers which are based on Opm::TimeMap");
        }

        return beginReportStepIdx_;
    }

    /// Index of the last considered simulation episode
    size_t SimulatorTimer::endReportStepIndex() const
    {
        if (!timeMap_) {
            OPM_THROW(std::runtime_error, "indexLastEpisode() is only implemented "
                      "for simulation timers which are based on Opm::TimeMap");
        }

        return endReportStepIdx_;
    }

    /// Current step number.
    int SimulatorTimer::currentStepNum() const
    {
        return current_step_ + beginReportStepIdx_;
    }

    /// Set current step number.
    void SimulatorTimer::setCurrentStepNum(int step)
    {
        if (current_step_ < 0 || current_step_ > int(numSteps())) {
            // Note that we do allow current_step_ == timesteps_.size(),
            // that is the done() state.
            OPM_THROW(std::runtime_error, "Trying to set invalid step number: " << step);
        }
        current_step_ = step;
        if (timeMap_)
            current_time_ = std::accumulate(timesteps_.begin(), timesteps_.begin() + step, 0.0);
    }


    /// Current step length.
    double SimulatorTimer::currentStepLength() const
    {
        assert(!done());
        if (timeMap_)
            return timeMap_->getTimeStepLength(beginReportStepIdx_ + current_step_);
        else
            return timesteps_[current_step_];
    }

    double SimulatorTimer::stepLengthTaken() const
    {
        assert(current_step_ > 0);
        if (timeMap_)
            return timeMap_->getTimeStepLength(beginReportStepIdx_ + current_step_ - 1);
        else
            return timesteps_[current_step_ - 1];
    }

    /// time elapsed since the start of the simulation [s].
    double SimulatorTimer::simulationTimeElapsed() const
    {
        if (timeMap_)
            return
                timeMap_->getTimePassedUntil(beginReportStepIdx_ + current_step_);
        else
            return current_time_;
    }

    /// time elapsed since the start of the POSIX epoch (Jan 1st, 1970) [s].
    time_t SimulatorTimer::currentPosixTime() const
    {
        tm t = boost::posix_time::to_tm(currentDateTime());
        return std::mktime(&t);
    }

    boost::posix_time::ptime SimulatorTimer::currentDateTime() const
    {
        if (timeMap_)
            return timeMap_->getStartTime(beginReportStepIdx_ + current_step_);
        else
            return boost::posix_time::ptime(start_date_) + boost::posix_time::seconds( (int) current_time_ );
    }



    /// Total time.
    double SimulatorTimer::totalTime() const
    {
        if (timeMap_)
            return
                timeMap_->getTotalTime();
        else
            return total_time_;
    }

    /// Set total time.
    /// This is primarily intended for multi-epoch schedules,
    /// where a timer for a given epoch does not have
    /// access to later timesteps.
    void SimulatorTimer::setTotalTime(double time)
    {
        if (timeMap_) {
            // well, what can we do if we use opm-parser's TimeMap?
            OPM_THROW(std::logic_error,
                      "Not implemented: SimulatorTimer::setTotalTime() if using a TimeMap.");
        }
        else
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
        if (!timeMap_)
            current_time_ += timesteps_[current_step_];
        ++current_step_;
        return *this;
    }

    /// Return true if op++() has been called numSteps() times.
    bool SimulatorTimer::done() const
    {
        if (timeMap_)
            return current_step_ > int(endReportStepIdx_ - beginReportStepIdx_ - 1);
        else
            return int(timesteps_.size()) == current_step_;
    }


} // namespace Opm
