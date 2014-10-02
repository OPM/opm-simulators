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
        virtual void initialize( const SimulatorState& state ) {}
        virtual double computeTimeStepSize( const double dt, const int iterations, const SimulatorState& ) const = 0;
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

        double computeTimeStepSize( const double dt, const int iterations, const SimulatorState& ) const
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

    class PIDTimeStepControl : public TimeStepControlInterface
    {
    protected:   
        mutable std::vector<double> p0_;
        mutable std::vector<double> sat0_;

        mutable double prevDt_;
        mutable int prevIterations_;
        const int targetIterationCount_;
        const double adjustmentFactor_;

        const int upperTargetIterationCount_;
        const int lowerTargetIterationCount_;

        const double tol_;

        mutable std::vector< double > errors_;

    public:
        PIDTimeStepControl( const double tol = 8e-4 ) 
            : p0_(), sat0_(), prevDt_( 0.0 ), prevIterations_( 0 ), 
              targetIterationCount_( 100 ), adjustmentFactor_( 1.25 ),
              upperTargetIterationCount_( 200 ), lowerTargetIterationCount_( 30 ),
              tol_( tol ),
              errors_( 3, tol_ )
        {}

        void initialize( const SimulatorState& state ) 
        {
            // store current state
            p0_   = state.pressure();
            sat0_ = state.saturation();
        }

        template <class Iterator>
        double inner_product( Iterator it, const Iterator end ) const 
        {
            double product = 0.0 ;
            for( ; it != end; ++it ) 
                product += ( *it * *it );
            return product;
        }

        double computeTimeStepSize( const double dt, const int iterations, const SimulatorState& state ) const
        {
            const size_t size = p0_.size();
            assert( state.pressure().size() == size );
            // compute u^n - u^n+1 
            for( size_t i=0; i<size; ++i ) 
            {
                p0_[ i ]   -= state.pressure()[ i ];
                sat0_[ i ] -= state.saturation()[ i ];
            }

            // compute || u^n - u^n+1 || 
            double stateN0  = inner_product( p0_.begin(),   p0_.end() ) +
                              inner_product( sat0_.begin(), sat0_.end() );

            // compute || u^n+1 || 
            double stateN   = inner_product( state.pressure().begin(), state.pressure().end() ) +
                              inner_product( state.saturation().begin(), state.saturation().end() );


            for( int i=0; i<2; ++i )
              errors_[ i ] = errors_[i+1];

            double error = stateN0 / stateN ;
            errors_[ 2 ] =  error ;

            prevDt_ = dt;
            prevIterations_ = iterations;

            if( error > tol_ )
            {
                // adjust dt by given tolerance 
                std::cout << "Computed step size (tol): " << (dt * tol_ / error ) << std::endl;
                return (dt * tol_ / error );
            }
            else
            {
                const double kP = 0.075 ;
                const double kI = 0.175 ;
                const double kD = 0.01 ;
                double newDt = (dt * std::pow( errors_[ 1 ] / errors_[ 2 ], kP ) *
                             std::pow( tol_         / errors_[ 2 ], kI ) *
                             std::pow( (errors_[0]*errors_[0]/errors_[ 1 ]*errors_[ 2 ]), kD ));
                std::cout << "Computed step size (pow): " << newDt << std::endl;
                return newDt;
            }
        }
    };

} // end namespace OPM

#endif
