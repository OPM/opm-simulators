/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.
  Copyright 2014 IRIS AS

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
#include <algorithm>
#include <numeric>

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
        double suggestedMax_;
        double suggestedAverage_;

        double computeInitialTimeStep( const double lastDt ) const 
        {
            const double maxTimeStep = simulator_timer_.currentStepLength();
            const double fraction = (lastDt / maxTimeStep);
            // when lastDt and maxTimeStep are close together, choose the max time step
            if( fraction > 0.95 ) return       maxTimeStep;
            // if lastDt is still pretty large, choose half step size since we have to
            // do two steps anyway
            // if( fraction > 0.85 ) return 0.5 * maxTimeStep;

            // otherwise choose lastDt
            return std::min( lastDt, maxTimeStep );
        }
    public:
        SubStepSimulatorTimer( const SimulatorTimer& st, const double lastDt )
            : simulator_timer_( st )
            , dt_( computeInitialTimeStep( lastDt ) )
            , current_time_( simulator_timer_.simulationTimeElapsed() )
            , total_time_( current_time_ + simulator_timer_.currentStepLength() )
            , sub_step_( 0 )
            , steps_()
            , suggestedMax_( 0.0 )
            , suggestedAverage_( 0.0 )
        {
            steps_.reserve( 10 );
        }

        void next( const double new_dt )
        {
            ++sub_step_;
            current_time_ += dt_;
            // store used time step sizes
            steps_.push_back( dt_ );

            // store some information about the time steps suggested 
            suggestedMax_      = std::max( new_dt, suggestedMax_ );
            suggestedAverage_ += new_dt;

            double remaining = (total_time_ - current_time_);

            if( remaining > 0 ) {

                // set new time step (depending on remaining time)
                if( 1.5 * new_dt > remaining ) {
                    dt_ = remaining; 
                    return ;
                }

                // check for half interval step to avoid very small step at the end
                // remaining *= 0.5;

                if( 2.25 * new_dt > remaining ) {
                    dt_ = 0.5 * remaining ;
                    return ;
                }
            }

            // otherwise set new_dt as is
            dt_ = new_dt;
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

        double averageStepLength() const
        {
            const int size = steps_.size();
            if( size == 0 ) return 0.0;

            const double sum = std::accumulate(steps_.begin(), steps_.end(), 0.0);
            return sum / double(size);
        }

        double maxStepLength () const
        {
            if( steps_.size() == 0 ) return 0.0;
            return *(std::max_element( steps_.begin(), steps_.end() ));
        }

        double minStepLength () const
        {
            if( steps_.size() == 0 ) return 0.0;
            return *(std::min_element( steps_.begin(), steps_.end() ));
        }

        double suggestedMax () const { return suggestedMax_; }
        double suggestedAverage () const 
        { 
            const int size = steps_.size();
            return (size > 0 ) ? (suggestedAverage_ / double(size)) : suggestedAverage_; 
        }

        void report(std::ostream& os) const
        {
            const double factor = 24.0 * 3600.0;
            os << "Sub steps started at time = " << simulator_timer_.simulationTimeElapsed()/factor << " (days)" << std::endl;
            for( size_t i=0; i<steps_.size(); ++i )
            {
                os << " step[ " << i << " ] = " << steps_[ i ]/factor << " (days)" << std::endl;
            }
            std::cout << "sub steps end time = " << simulationTimeElapsed()/factor << " (days)" << std::endl;
        }
    };


} // namespace Opm

#endif // OPM_SIMULATORTIMER_HEADER_INCLUDED
