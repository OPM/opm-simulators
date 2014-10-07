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

#include <cassert>
#include <iostream>
#include <vector>

#include <algorithm>
#include <numeric>

#include <opm/core/utility/Units.hpp>
#include <opm/core/simulator/AdaptiveSimulatorTimer.hpp>

namespace Opm
{
    AdaptiveSimulatorTimer::
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

    AdaptiveSimulatorTimer& AdaptiveSimulatorTimer::operator++ ()
    {
        ++current_step_;
        current_time_ += dt_;
        // store used time step sizes
        steps_.push_back( dt_ );
        return *this;
    }

    void AdaptiveSimulatorTimer::
    provideTimeStepEstimate( const double dt_estimate )
    {
        // store some information about the time steps suggested
        suggestedMax_      = std::max( dt_estimate, suggestedMax_ );
        suggestedAverage_ += dt_estimate;

        double remaining = (total_time_ - current_time_);

        if( remaining > 0 ) {

            // set new time step (depending on remaining time)
            if( 1.5 * dt_estimate > remaining ) {
                dt_ = remaining;
                return ;
            }

            // check for half interval step to avoid very small step at the end
            // remaining *= 0.5;

            if( 2.25 * dt_estimate > remaining ) {
                dt_ = 0.5 * remaining;
                return ;
            }
        }

        // otherwise set dt_estimate as is
        dt_ = dt_estimate;
    }

    int AdaptiveSimulatorTimer::
    currentStepNum () const { return current_step_; }

    double AdaptiveSimulatorTimer::currentStepLength () const
    {
        assert( ! done () );
        return dt_;
    }

    double AdaptiveSimulatorTimer::totalTime() const { return total_time_; }

    double AdaptiveSimulatorTimer::simulationTimeElapsed() const { return current_time_; }

    bool AdaptiveSimulatorTimer::done () const { return (current_time_ >= total_time_) ; }

    double AdaptiveSimulatorTimer::averageStepLength() const
    {
        const int size = steps_.size();
        if( size == 0 ) return 0.0;

        const double sum = std::accumulate(steps_.begin(), steps_.end(), 0.0);
        return sum / double(size);
    }

    /// \brief return max step length used so far
    double AdaptiveSimulatorTimer::maxStepLength () const
    {
        if( steps_.size() == 0 ) return 0.0;
        return *(std::max_element( steps_.begin(), steps_.end() ));
    }

    /// \brief return min step length used so far
    double AdaptiveSimulatorTimer::minStepLength () const
    {
        if( steps_.size() == 0 ) return 0.0;
        return *(std::min_element( steps_.begin(), steps_.end() ));
    }

    /// \brief return max suggested step length
    double AdaptiveSimulatorTimer::suggestedMax () const { return suggestedMax_; }

    /// \brief return average suggested step length
    double AdaptiveSimulatorTimer::suggestedAverage () const
    {
        const int size = steps_.size();
        return (size > 0 ) ? (suggestedAverage_ / double(size)) : suggestedAverage_;
    }

    /// \brief report start and end time as well as used steps so far
    void AdaptiveSimulatorTimer::
    report(std::ostream& os) const
    {
        os << "Sub steps started at time = " <<  unit::convert::to( start_time_, unit::day ) << " (days)" << std::endl;
        for( size_t i=0; i<steps_.size(); ++i )
        {
            os << " step[ " << i << " ] = " << unit::convert::to( steps_[ i ], unit::day ) << " (days)" << std::endl;
        }
        std::cout << "sub steps end time = " << unit::convert::to( simulationTimeElapsed(), unit::day ) << " (days)" << std::endl;
    }

    double AdaptiveSimulatorTimer::
    computeInitialTimeStep( const double lastDt ) const
    {
        const double maxTimeStep = total_time_ - start_time_;
        const double fraction = (lastDt / maxTimeStep);
        // when lastDt and maxTimeStep are close together, choose the max time step
        if( fraction > 0.95 ) return       maxTimeStep;

        // otherwise choose lastDt
        return std::min( lastDt, maxTimeStep );
    }

} // namespace Opm
