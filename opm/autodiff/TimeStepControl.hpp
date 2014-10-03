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

    ///////////////////////////////////////////////////////////////////
    ///
    ///  TimeStepControlInterface 
    /// 
    ///////////////////////////////////////////////////////////////////
    class TimeStepControlInterface 
    {
    protected:    
        TimeStepControlInterface() {}
    public:
        /// \param state simulation state before computing update in the solver (default is empty)
        virtual void initialize( const SimulatorState& state ) {}

        /// compute new time step size suggestions based on the PID controller
        /// \param dt          time step size used in the current step
        /// \param iterations  number of iterations used (linear/nonlinear) 
        /// \param state       new solution state
        ///
        /// \return suggested time step size for the next step
        virtual double computeTimeStepSize( const double dt, const int iterations, const SimulatorState& ) const = 0;

        /// virtual destructor (empty)
        virtual ~TimeStepControlInterface () {}
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///
    ///  PID controller based adaptive time step control as suggested in: 
    ///  Turek and Kuzmin. Algebraic Flux Correction III. Incompressible Flow Problems. Uni Dortmund. 
    /// 
    ///  See also: 
    ///  D. Kuzmin and S.Turek. Numerical simulation of turbulent bubbly flows. Techreport Uni Dortmund. 2004
    ///  
    ///  and the original article: 
    ///  Valli, Coutinho, and Carey. Adaptive Control for Time Step Selection in Finite Element
    ///  Simulation of Coupled Viscous Flow and Heat Transfer. Proc of the 10th 
    ///  International Conference on Numerical Methods in Fluids. 1998.
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class PIDTimeStepControl : public TimeStepControlInterface
    {
    protected:   
        mutable std::vector<double> p0_;
        mutable std::vector<double> sat0_;

        const double tol_;
        mutable std::vector< double > errors_;

        const bool verbose_;
    public:
        /// \brief constructor 
        /// \param tol  tolerance for the relative changes of the numerical solution to be accepted 
        ///             in one time step (default is 1e-3)
        PIDTimeStepControl( const double tol = 1e-3, const bool verbose = false ) 
            : p0_()
            , sat0_() 
            , tol_( tol )
            , errors_( 3, tol_ )
            , verbose_( verbose )
        {}

        /// \brief \copydoc TimeStepControlInterface::initialize
        void initialize( const SimulatorState& state ) 
        {
            // store current state for later time step computation
            p0_   = state.pressure();
            sat0_ = state.saturation();
        }

        /// \brief \copydoc TimeStepControlInterface::computeTimeStepSize
        double computeTimeStepSize( const double dt, const int /* iterations */, const SimulatorState& state ) const
        {
            const size_t size = p0_.size();
            assert( state.pressure().size() == size );
            assert( state.saturation().size() == size );
            assert( sat0_.size() == size );

            // compute u^n - u^n+1 
            for( size_t i=0; i<size; ++i ) 
            {
                p0_[ i ]   -= state.pressure()[ i ];
                sat0_[ i ] -= state.saturation()[ i ];
            }

            // compute || u^n - u^n+1 || 
            const double stateOld  = inner_product( p0_.begin(),   p0_.end() ) +
                                     inner_product( sat0_.begin(), sat0_.end() );

            // compute || u^n+1 || 
            const double stateNew  = inner_product( state.pressure().begin(), state.pressure().end() ) +
                                     inner_product( state.saturation().begin(), state.saturation().end() );

            // shift errors
            for( int i=0; i<2; ++i )
              errors_[ i ] = errors_[i+1];

            // store new error
            const double error = stateOld / stateNew;
            errors_[ 2 ] =  error ;

            if( error > tol_ )
            {
                // adjust dt by given tolerance 
                if( verbose_ )
                    std::cout << "Computed step size (tol): " << (dt * tol_ / error )/86400.0 << " (days)" << std::endl;
                return (dt * tol_ / error );
            }
            else
            {
                // values taking from turek time stepping paper
                const double kP = 0.075 ;
                const double kI = 0.175 ;
                const double kD = 0.01 ;
                double newDt = (dt * std::pow( errors_[ 1 ] / errors_[ 2 ], kP ) *
                             std::pow( tol_         / errors_[ 2 ], kI ) *
                             std::pow( (errors_[0]*errors_[0]/errors_[ 1 ]*errors_[ 2 ]), kD ));
                if( verbose_ )
                    std::cout << "Computed step size (pow): " << newDt/86400.0 << " (days)" << std::endl;
                return newDt;
            }
        }

    protected:    
        // return inner product for given container, here std::vector
        template <class Iterator>
        double inner_product( Iterator it, const Iterator end ) const 
        {
            double product = 0.0 ;
            for( ; it != end; ++it ) 
                product += ( *it * *it );
            return product;
        }

    };

    class PIDAndIterationCountTimeStepControl : public PIDTimeStepControl
    {
        typedef PIDTimeStepControl BaseType;
    protected:   
        const int targetIterationCount_;

    public:
        PIDAndIterationCountTimeStepControl( const int target_iterations = 20,
                                             const double tol = 1e-3, 
                                             const bool verbose = false) 
            : BaseType( tol, verbose )
            , targetIterationCount_( target_iterations )
        {}

        double computeTimeStepSize( const double dt, const int iterations, const SimulatorState& state ) const
        {
            double dtEstimate = BaseType :: computeTimeStepSize( dt, iterations, state );

            // further reduce step size if to many iterations were used
            if( iterations > targetIterationCount_ ) 
            {
                // if iterations was the same or dts were the same, do some magic
                dtEstimate *= double( targetIterationCount_ ) / double(iterations);
            }

            return dtEstimate;
        }
    };


} // end namespace OPM
#endif
