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

#include <opm/parser/eclipse/EclipseState/Schedule/TimeMap.hpp>

#include <iosfwd>
#include <vector>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>

namespace Opm
{

    namespace parameter { class ParameterGroup; }

    class SimulatorTimer
    {
    public:
        /// Default constructor.
        SimulatorTimer();

        /// Initialize from parameters. Accepts the following:
        ///    num_psteps    (default 1)
        ///    stepsize_days (default 1)
        void init(const parameter::ParameterGroup& param);

        /// Use the SimulatorTimer as a shim around opm-parser's Opm::TimeMap
        void init(TimeMapConstPtr timeMap);

        /// Total number of steps.
        int numSteps() const;

        /// Current step number. This is the number of timesteps that
        /// has been completed from the start of the run. The time
        /// after initialization but before the simulation has started
        /// is timestep number zero.
        int currentStepNum() const;

        /// Set current step number.
        void setCurrentStepNum(int step);

        /// Current step length. This is the length of the step
        /// the simulator will take in the next iteration.
        ///
        /// @note if done(), it is an error to call currentStepLength().
        double currentStepLength() const;

        /// Previous step length. This is the length of the step that
        /// was taken to arrive at this time.
        ///
        /// @note if no increments have been done (i.e. the timer is
        /// still in its constructed state and currentStepNum() == 0),
        /// it is an error to call stepLengthTaken().
        double stepLengthTaken () const;

        /// Time elapsed since the start of the POSIX epoch (Jan 1st,
        /// 1970) until the current time step begins [s].
        time_t currentPosixTime() const;

        /// Time elapsed since the start of the simulation until the
        /// beginning of the current time step [s].
        double simulationTimeElapsed() const;

        /// Return the current time as a posix time object.
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

    class SubStepSimulatorTimer
    {
    protected:
        const SimulatorTimer& simulator_timer_;
        double dt_;
        double current_time_;
        const double total_time_;
        int sub_step_;

        std::vector< double > steps_;
    public:
        SubStepSimulatorTimer( const SimulatorTimer& st, const double lastDt )
            : simulator_timer_( st )
            , dt_( lastDt )
            , current_time_( simulator_timer_.simulationTimeElapsed() )
            , total_time_( current_time_ + simulator_timer_.currentStepLength() )
            , sub_step_( 0 )
            , steps_()
        {
            steps_.reserve( 10 );
        }

        void next( const double new_dt )
        {
            ++sub_step_;
            current_time_ += dt_;
            // store used time step sizes
            steps_.push_back( dt_ );

            double remaining = total_time_ - current_time_;
            // set new time step (depending on remaining time)
            dt_ = ( new_dt > remaining && remaining > 0 ) ? remaining : new_dt ;
        }

        int currentStepNum () const { return sub_step_; }

        double currentStepLength () const
        {
            assert( ! done () );
            return dt_;
        }

        double totalTime() const { return total_time_; }

        double simulationTimeElapsed() const { return current_time_; }

        bool done () const { return (current_time_ >= total_time_) ; }

        void report(std::ostream& os) const
        {
            os << "Sub steps started at time = " << simulator_timer_.simulationTimeElapsed() << std::endl;
            for( size_t i=0; i<steps_.size(); ++i )
            {
                os << " step[ " << i << " ] = " << steps_[ i ] << std::endl;
            }
            std::cout << "sub steps end time = " << simulationTimeElapsed() << std::endl;
        }
    };


} // namespace Opm

#endif // OPM_SIMULATORTIMER_HEADER_INCLUDED
