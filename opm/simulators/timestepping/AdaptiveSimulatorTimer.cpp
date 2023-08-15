/*
  Copyright (c) 2014 IRIS AS

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

#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <ostream>
#include <numeric>
#include <vector>

#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/simulators/timestepping/AdaptiveSimulatorTimer.hpp>

#include <boost/date_time/posix_time/posix_time.hpp>

namespace Opm
{
    AdaptiveSimulatorTimer::
    AdaptiveSimulatorTimer( const SimulatorTimerInterface& timer,
                            const double lastStepTaken,
                            const double maxTimeStep )
        : start_date_time_( std::make_shared<boost::posix_time::ptime>(timer.startDateTime()) )
        , start_time_( timer.simulationTimeElapsed() )
        , total_time_( start_time_ + timer.currentStepLength() )
        , report_step_( timer.reportStepNum() )
        , max_time_step_( maxTimeStep )
        , current_time_( start_time_ )
        , dt_( 0.0 )
        , current_step_( 0 )
        , steps_()
        , lastStepFailed_( false )
    {
        // reserve memory for sub steps
        steps_.reserve( 10 );

        // set appropriate value for dt_
        provideTimeStepEstimate( lastStepTaken );
    }

    bool AdaptiveSimulatorTimer::initialStep () const
    {
        return ( report_step_ == 0 ) && ( current_step_ == 0 );
    }

AdaptiveSimulatorTimer& AdaptiveSimulatorTimer::operator++ ()
    {
        ++current_step_;
        current_time_ += dt_;
        assert(dt_ > 0);
        // store used time step sizes
        steps_.push_back( dt_ );
        return *this;
    }

    void AdaptiveSimulatorTimer::
    provideTimeStepEstimate( const double dt_estimate )
    {
        double remaining = (total_time_ - current_time_);
        // apply max time step if it was set
        dt_ = std::min( dt_estimate, max_time_step_ );
        assert(dt_ > 0);
        if( remaining > 0 ) {

            // set new time step (depending on remaining time)
            if( 1.05 * dt_ > remaining ) {
                dt_ = remaining;
                // check max time step again and use half remaining if too large
                if( dt_ > max_time_step_ ) {
                    dt_ = 0.5 * remaining;
                }
                assert(dt_ > 0);
                return;
            }

            // check for half interval step to avoid very small step at the end
            // remaining *= 0.5;

            if( 1.5 * dt_ > remaining ) {
                dt_ = 0.5 * remaining;
                assert(dt_ > 0);
                return;
            }
        }
    }

    int AdaptiveSimulatorTimer::
    currentStepNum () const { return current_step_; }

    int AdaptiveSimulatorTimer::
    reportStepNum () const { return report_step_; }

    double AdaptiveSimulatorTimer::currentStepLength () const
    {
      assert(dt_ > 0);
        return dt_;
    }

    double AdaptiveSimulatorTimer::stepLengthTaken() const
    {
        assert( ! steps_.empty() );
        return steps_.back();
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
        if( steps_.empty() ) return 0.0;
        return *(std::max_element( steps_.begin(), steps_.end() ));
    }

    /// \brief return min step length used so far
    double AdaptiveSimulatorTimer::minStepLength () const
    {
        if( steps_.empty() ) return 0.0;
        return *(std::min_element( steps_.begin(), steps_.end() ));
    }

    /// \brief report start and end time as well as used steps so far
    void AdaptiveSimulatorTimer::
    report(std::ostream& os) const
    {
        os << "Sub steps started at time = " <<  unit::convert::to( start_time_, unit::day ) << " (days)" << std::endl;
        for (std::size_t i = 0; i < steps_.size(); ++i)
        {
            os << " step[ " << i << " ] = " << unit::convert::to( steps_[ i ], unit::day ) << " (days)" << std::endl;
        }
        os << "sub steps end time = " << unit::convert::to( simulationTimeElapsed(), unit::day ) << " (days)" << std::endl;
    }

    boost::posix_time::ptime AdaptiveSimulatorTimer::startDateTime() const
    {
        return *start_date_time_;
    }

    /// return copy of object
    std::unique_ptr< SimulatorTimerInterface >
    AdaptiveSimulatorTimer::clone() const
    {
        return std::make_unique<AdaptiveSimulatorTimer>(*this);
    }



} // namespace Opm
