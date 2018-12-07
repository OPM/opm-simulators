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
#include <opm/parser/eclipse/Units/Units.hpp>
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
    void SimulatorTimer::init(const ParameterGroup& param)
    {
        const int num_psteps = param.getDefault("num_psteps", 1);
        const double stepsize_days = param.getDefault("stepsize_days", 1.0);
        const double stepsize = Opm::unit::convert::from(stepsize_days, Opm::unit::day);
        timesteps_.clear();
        timesteps_.resize(num_psteps, stepsize);
        total_time_ = num_psteps*stepsize;
        for ( size_t i = 0; i < num_psteps; ++i ) {
            report_stepindx_[i] = i;
        }
    }
    void SimulatorTimer::init(std::vector<double> time_steps,std::vector<int> report_stepindx){
        timesteps_.clear();
        timesteps_ = time_steps;
        report_stepindx_ = report_stepindx;
        total_time_ = 0;
        for (auto& dt : timesteps_){
            total_time_ += dt;
       }
    }

    /// Use the SimulatorTimer as a shim around opm-parser's Opm::TimeMap
    void SimulatorTimer::init(const TimeMap& timeMap, size_t report_step)
    {
        total_time_ = timeMap.getTotalTime();
        timesteps_.resize(timeMap.numTimesteps());
        report_stepindx_.resize(timeMap.numTimesteps());
        for ( size_t i = 0; i < timeMap.numTimesteps(); ++i ) {
            timesteps_[i] = timeMap.getTimeStepLength(i);
            report_stepindx_[i] = i;
        }

        setCurrentStepNum(report_step);
        start_date_ = boost::posix_time::from_time_t( timeMap.getStartTime(0)).date();
    }

    /// Whether the current step is the first step.
    bool SimulatorTimer::initialStep() const
    {
        return (current_step_ == 0);
    }

    /// Total number of steps.
    int SimulatorTimer::numSteps() const
    {
        return int(timesteps_.size());
    }

    /// Current step number.
    int SimulatorTimer::currentStepNum() const
    {
        return current_step_;
    }
    int SimulatorTimer::prevReportStepNum() const{
        return report_stepindx_[current_step_-1];
    }

    int SimulatorTimer::reportStepNum() const
    {
       assert(current_step_ <= report_stepindx_.size());
       int report_step = -100;
       if(current_step_ == report_stepindx_.size()){
           // we are at the end of the last time step
            report_step =  int(report_stepindx_.size());
       }else{
            report_step = report_stepindx_[current_step_];
       }
       return report_step;
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
        //assert(!done());
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
        double ctime = 0.0;
        if(std::abs(current_time_) < 1e-7){
            ctime = 0.0;// to avoid rounding error for file names
        }else{
            ctime = current_time_;
        }
        return ctime;
    }

    boost::posix_time::ptime SimulatorTimer::startDateTime() const
    {
        return boost::posix_time::ptime(start_date_);
    }


    boost::posix_time::ptime SimulatorTimer::currentDateTime() const
    {
        // Boost uses only 32 bit long for seconds, but 64 bit for milliseconds.
        // As a workaround for very large times we just use milliseconds.
        // The cast is necessary because boost::posix_time::milliseconds requires
        // an integer argument.
        return startDateTime() + boost::posix_time::milliseconds(static_cast<long long>(simulationTimeElapsed() / Opm::prefix::milli));
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
        os << "\n\n---------------    Simulation report step number " << reportStepNum() << "    ---------------";
        os << "\n\n---------------    Simulation step number " << currentStepNum() <<           "    ---------------";
        os   << "\n      Current time (days)     " << Opm::unit::convert::to(simulationTimeElapsed(), Opm::unit::day);
        if(done()){
           os << "\n      Current stepsize (days) " << " LAST STEP ";
        }else{
           os << "\n      Current stepsize (days) " << Opm::unit::convert::to(currentStepLength(), Opm::unit::day);
        }
        os   << "\n      Total time (days)       " << Opm::unit::convert::to(totalTime(), Opm::unit::day)
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
     /// Previous step
    SimulatorTimer& SimulatorTimer::operator--()
    {
        assert(! initialStep() );
         --current_step_;
        current_time_ -= timesteps_[current_step_];
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
       return std::unique_ptr< SimulatorTimerInterface > (new SimulatorTimer( *this ));
    }



} // namespace Opm
