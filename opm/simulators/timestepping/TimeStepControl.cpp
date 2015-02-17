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
#include <cmath>
#include <iostream>
#include <stdexcept>

#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/simulator/TimeStepControl.hpp>

namespace Opm
{
    ////////////////////////////////////////////////////////
    //
    //  InterationCountTimeStepControl Implementation
    //
    ////////////////////////////////////////////////////////

    SimpleIterationCountTimeStepControl::
    SimpleIterationCountTimeStepControl( const int target_iterations,
                                         const double decayrate,
                                         const double growthrate,
                                         const bool verbose)
        : target_iterations_( target_iterations )
        , decayrate_( decayrate )
        , growthrate_( growthrate )
        , verbose_( verbose )
    {
        if( decayrate_  > 1.0 ) {
            OPM_THROW(std::runtime_error,"SimpleIterationCountTimeStepControl: decay should be <= 1 " << decayrate_ );
        }
        if( growthrate_ < 1.0 ) {
            OPM_THROW(std::runtime_error,"SimpleIterationCountTimeStepControl: growth should be >= 1 " << growthrate_ );
        }
    }

    double SimpleIterationCountTimeStepControl::
    computeTimeStepSize( const double dt, const int iterations, const SimulatorState& /* state */ ) const
    {
        double dtEstimate = dt ;

        // reduce the time step size if we exceed the number of target iterations
        if( iterations > target_iterations_ )
        {
            // scale dtEstimate down with a given rate
            dtEstimate *= decayrate_;
        }
        // increase the time step size if we are below the number of target iterations
        else if ( iterations < target_iterations_-1 )
        {
            dtEstimate *= growthrate_;
        }

        return dtEstimate;
    }



    ////////////////////////////////////////////////////////
    //
    //  PIDTimeStepControl Implementation
    //
    ////////////////////////////////////////////////////////

    PIDTimeStepControl::PIDTimeStepControl( const double tol, const bool verbose )
        : p0_()
        , sat0_()
        , tol_( tol )
        , errors_( 3, tol_ )
        , verbose_( verbose )
    {}

    void PIDTimeStepControl::initialize( const SimulatorState& state )
    {
        // store current state for later time step computation
        p0_   = state.pressure();
        sat0_ = state.saturation();
    }

    double PIDTimeStepControl::
    computeTimeStepSize( const double dt, const int /* iterations */, const SimulatorState& state ) const
    {
        const std::size_t pSize = p0_.size();
        assert( state.pressure().size() == pSize );
        const std::size_t satSize = sat0_.size();
        assert( state.saturation().size() == satSize );

        // compute u^n - u^n+1
        for( std::size_t i=0; i<pSize; ++i ) {
            p0_[ i ]   -= state.pressure()[ i ];
        }

        for( std::size_t i=0; i<satSize; ++i ) {
            sat0_[ i ] -= state.saturation()[ i ];
        }

        // compute || u^n - u^n+1 ||
        const double stateOld  = euclidianNormSquared( p0_.begin(),   p0_.end() ) +
                                 euclidianNormSquared( sat0_.begin(), sat0_.end() );

        // compute || u^n+1 ||
        const double stateNew  = euclidianNormSquared( state.pressure().begin(),   state.pressure().end()   ) +
                                 euclidianNormSquared( state.saturation().begin(), state.saturation().end() );

        // shift errors
        for( int i=0; i<2; ++i ) {
          errors_[ i ] = errors_[i+1];
        }

        // store new error
        const double error = stateOld / stateNew;
        errors_[ 2 ] =  error ;

        if( error > tol_ )
        {
            // adjust dt by given tolerance
            const double newDt = dt * tol_ / error;
            if( verbose_ )
                std::cout << "Computed step size (tol): " << unit::convert::to( newDt, unit::day ) << " (days)" << std::endl;
            return newDt;
        }
        else
        {
            // values taking from turek time stepping paper
            const double kP = 0.075 ;
            const double kI = 0.175 ;
            const double kD = 0.01 ;
            const double newDt = (dt * std::pow( errors_[ 1 ] / errors_[ 2 ], kP ) *
                                 std::pow( tol_         / errors_[ 2 ], kI ) *
                                 std::pow( errors_[0]*errors_[0]/errors_[ 1 ]/errors_[ 2 ], kD ));
            if( verbose_ )
                std::cout << "Computed step size (pow): " << unit::convert::to( newDt, unit::day ) << " (days)" << std::endl;
            return newDt;
        }
    }



    ////////////////////////////////////////////////////////////
    //
    //  PIDAndIterationCountTimeStepControl  Implementation
    //
    ////////////////////////////////////////////////////////////

    PIDAndIterationCountTimeStepControl::
    PIDAndIterationCountTimeStepControl( const int target_iterations,
                                         const double tol,
                                         const double maxgrowth,
                                         const bool verbose)
        : BaseType( tol, verbose )
        , target_iterations_( target_iterations )
        , maxgrowth_( maxgrowth )
    {}

    double PIDAndIterationCountTimeStepControl::
    computeTimeStepSize( const double dt, const int iterations, const SimulatorState& state ) const
    {
        double dtEstimate = BaseType :: computeTimeStepSize( dt, iterations, state );

        // further reduce step size if to many iterations were used
        if( iterations > target_iterations_ )
        {
            // if iterations was the same or dts were the same, do some magic
            dtEstimate *= double( target_iterations_ ) / double(iterations);
        }

        // limit the growth of the timestep size by the growth factor
        dtEstimate = std::min( dtEstimate, double(maxgrowth_ * dt) );

        return dtEstimate;
    }

} // end namespace Opm
