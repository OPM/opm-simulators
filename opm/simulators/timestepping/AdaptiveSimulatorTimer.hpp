/*
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

#ifndef OPM_ADAPTIVESIMULATORTIMER_HEADER_INCLUDED
#define OPM_ADAPTIVESIMULATORTIMER_HEADER_INCLUDED

#include <cassert>
#include <iostream>
#include <vector>

#include <algorithm>
#include <numeric>

namespace Opm
{

    /////////////////////////////////////////////////////////
    ///
    /// \brief Simulation timer for adaptive time stepping 
    /// 
    /////////////////////////////////////////////////////////
    class AdaptiveSimulatorTimer
    {
    protected:
        const double start_time_;
        const double total_time_;
        double current_time_;
        double dt_;
        int current_step_;

        std::vector< double > steps_;
        double suggestedMax_;
        double suggestedAverage_;

        double computeInitialTimeStep( const double lastDt ) const 
        {
            const double maxTimeStep = total_time_ - start_time_;
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
        /// \brief constructor taking a simulator timer to determine start and end time 
        ///  \param start_time     start time of timer
        ///  \param total_time     total time of timer
        ///  \param lastDt         last suggested length of time step interval
        AdaptiveSimulatorTimer( const double start_time, const double total_time, const double lastDt )
            : start_time_( start_time )
            , total_time_( total_time )
            , current_time_( start_time_ )
            , dt_( computeInitialTimeStep( lastDt ) )
            , current_step_( 0 )
            , steps_()
            , suggestedMax_( 0.0 )
            , suggestedAverage_( 0.0 )
        {
            // reserve memory for sub steps
            steps_.reserve( 10 );
        }

        /// \brief decrease current time step for factor of two
        void halfTimeStep() { dt_ *= 0.5; } 

        /// \brief advance time by currentStepLength and set new step lenght
        void advance( const double new_dt )
        {
            ++current_step_;
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

        /// \brief \copydoc SimulationTimer::currentStepNum 
        int currentStepNum () const { return current_step_; }

        /// \brief \copydoc SimulationTimer::currentStepLength 
        double currentStepLength () const
        {
            assert( ! done () );
            return dt_;
        }

        /// \brief \copydoc SimulationTimer::totalTime
        double totalTime() const { return total_time_; }

        /// \brief \copydoc SimulationTimer::simulationTimeElapsed
        double simulationTimeElapsed() const { return current_time_; }

        /// \brief \copydoc SimulationTimer::done
        bool done () const { return (current_time_ >= total_time_) ; }

        /// \brief return average step length used so far
        double averageStepLength() const
        {
            const int size = steps_.size();
            if( size == 0 ) return 0.0;

            const double sum = std::accumulate(steps_.begin(), steps_.end(), 0.0);
            return sum / double(size);
        }

        /// \brief return max step length used so far
        double maxStepLength () const
        {
            if( steps_.size() == 0 ) return 0.0;
            return *(std::max_element( steps_.begin(), steps_.end() ));
        }

        /// \brief return min step length used so far
        double minStepLength () const
        {
            if( steps_.size() == 0 ) return 0.0;
            return *(std::min_element( steps_.begin(), steps_.end() ));
        }

        /// \brief return max suggested step length
        double suggestedMax () const { return suggestedMax_; }

        /// \brief return average suggested step length
        double suggestedAverage () const 
        { 
            const int size = steps_.size();
            return (size > 0 ) ? (suggestedAverage_ / double(size)) : suggestedAverage_; 
        }

        /// \brief report start and end time as well as used steps so far
        void report(std::ostream& os) const
        {
            const double factor = 24.0 * 3600.0;
            os << "Sub steps started at time = " << start_time_/factor << " (days)" << std::endl;
            for( size_t i=0; i<steps_.size(); ++i )
            {
                os << " step[ " << i << " ] = " << steps_[ i ]/factor << " (days)" << std::endl;
            }
            std::cout << "sub steps end time = " << simulationTimeElapsed()/factor << " (days)" << std::endl;
        }
    };


} // namespace Opm

#endif // OPM_SIMULATORTIMER_HEADER_INCLUDED
