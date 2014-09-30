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
#ifndef OPM_TIMESTEPCONTROL_HEADER_INCLUDED
#define OPM_TIMESTEPCONTROL_HEADER_INCLUDED

namespace Opm
{

    class TimeStepControlInterface 
    {
    protected:    
        TimeStepControlInterface() {}
    public:
        virtual double computeTimeStepSize( const double dt, const int iterations ) const = 0;
        virtual ~TimeStepControlInterface () {}
    };

    class IterationCountTimeStepControl : public TimeStepControlInterface
    {
    protected:    
        mutable double prevDt_;
        mutable int prevIterations_;
        const int targetIterationCount_;
        const double adjustmentFactor_;

        const int upperTargetIterationCount_;
        const int lowerTargetIterationCount_;

    public:
        IterationCountTimeStepControl() 
            : prevDt_( 0.0 ), prevIterations_( 0 ), 
              targetIterationCount_( 100 ), adjustmentFactor_( 1.25 ),
              upperTargetIterationCount_( 200 ), lowerTargetIterationCount_( 30 )
        {}

        double computeTimeStepSize( const double dt, const int iterations ) const
        {
            // make sure dt is somewhat reliable
            assert( dt > 0 && dt == dt );
            double newDt = dt;
            double derivation = double(std::abs( iterations - targetIterationCount_ )) / double(targetIterationCount_);

            if( derivation > 0.1 )
            {
                if( iterations < targetIterationCount_ ) 
                {
                    newDt = dt * adjustmentFactor_;
                }
                else 
                    newDt = dt / adjustmentFactor_;
            }
                           
            /*
            if( prevDt_ > 0 && std::abs( dt - prevDt_ ) > 1e-12 ) {
                const double dFdt  = double(iterations - prevIterations_) / ( dt - prevDt_ );
                if( std::abs( dFdt ) > 1e-12 )
                    newDt = dt + (targetIterationCount_ - iterations) / dFdt;
                else
                    // if iterations was the same or dts were the same, do some magic
                    newDt = dt * double( targetIterationCount_ ) / double(targetIterationCount_ - iterations);
            }

            if( newDt < 0 || ! (prevDt_ > 0) || ( iterations == prevIterations_) ) 
            {
                if( iterations > upperTargetIterationCount_ )
                    newDt = dt / adjustmentFactor_;
                else if( iterations < lowerTargetIterationCount_ )
                    newDt = dt * adjustmentFactor_;
                else
                    newDt = dt;
            }
            */

            assert( newDt == newDt );

            //std::cout << "dt = " << dt << "  " << prevDt_ << std::endl;

            prevDt_ = dt;
            prevIterations_ = iterations;

            return newDt;
        }
    };

} // end namespace OPM

#endif
