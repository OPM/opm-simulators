/*
  Copyright 2014 IRIS AS
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 Statoil AS

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
#include <config.h>
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <fstream>
#include <iostream>

#include <opm/common/ErrorMacros.hpp>
#include <opm/parser/eclipse/Units/Units.hpp>
#include <opm/simulators/timestepping/TimeStepControl.hpp>

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
    computeTimeStepSize( const double dt, const int iterations, const RelativeChangeInterface& /* relativeChange */, const double /*simulationTimeElapsed */) const
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
    //  HardcodedTimeStepControl Implementation
    //
    ////////////////////////////////////////////////////////

    HardcodedTimeStepControl::
    HardcodedTimeStepControl( const std::string& filename)
    {
        std::ifstream infile (filename);
        if (!infile.is_open()) {
            OPM_THROW(std::runtime_error,"Incorrect or no filename is provided to the hardcodedTimeStep. Use timestep.control.filename=your_file_name");
        }
        std::string::size_type sz;
        std::string line;
        while ( std::getline(infile, line)) {
            if( line[0] != '-') { // ignore lines starting with '-'
                const double time = std::stod(line,&sz); // read the first number i.e. the actual substep time
                subStepTime_.push_back( time * unit::day );
            }

        }
    }

    double HardcodedTimeStepControl::
    computeTimeStepSize( const double /*dt */, const int /*iterations */, const RelativeChangeInterface& /* relativeChange */ , const double simulationTimeElapsed) const
    {
        auto nextTime = std::upper_bound(subStepTime_.begin(), subStepTime_.end(), simulationTimeElapsed);
        return (*nextTime - simulationTimeElapsed);
    }



    ////////////////////////////////////////////////////////
    //
    //  PIDTimeStepControl Implementation
    //
    ////////////////////////////////////////////////////////

    PIDTimeStepControl::PIDTimeStepControl( const double tol,
                                            const bool verbose )
        : tol_( tol )
        , errors_( 3, tol_ )
        , verbose_( verbose )
    {}

    double PIDTimeStepControl::
    computeTimeStepSize( const double dt, const int /* iterations */, const RelativeChangeInterface& relChange, const double /*simulationTimeElapsed */) const
    {
        // shift errors
        for( int i=0; i<2; ++i ) {
            errors_[ i ] = errors_[i+1];
        }

        // store new error
        const double error = relChange.relativeChange();
        errors_[ 2 ] = error;
        for( int i=0; i<2; ++i ) {
            assert(std::isfinite(errors_[i]));
            //assert(errors_[i]>0);
            if(errors_[i]==0){
                return dt*10;
            }
        }

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
                                         const bool verbose)
        : BaseType( tol, verbose )
        , target_iterations_( target_iterations )
    {}

    double PIDAndIterationCountTimeStepControl::
    computeTimeStepSize( const double dt, const int iterations, const RelativeChangeInterface& relChange,  const double simulationTimeElapsed ) const
    {
        double dtEstimatePID = BaseType :: computeTimeStepSize( dt, iterations, relChange, simulationTimeElapsed);

        // adjust timesteps based on target iteration
        double dtEstimateIter;
        if (iterations > target_iterations_) {
            double off_target_fraction = double(iterations - target_iterations_) / target_iterations_;
            dtEstimateIter = dt / (1.0 + off_target_fraction);
        } else {
            double off_target_fraction = double(target_iterations_ - iterations) / target_iterations_;
            // Be a bit more careful when increasing. The 1.2 factor is from ebos.
            dtEstimateIter = dt * (1.0 + off_target_fraction / 1.2);
        }

        return std::min(dtEstimatePID, dtEstimateIter);
    }

} // end namespace Opm
