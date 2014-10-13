#ifndef OPM_SUBSTEPPING_HEADER_INCLUDED
#define OPM_SUBSTEPPING_HEADER_INCLUDED

#include <iostream>
#include <utility>

#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/simulator/TimeStepControlInterface.hpp>

namespace Opm {

    // AdaptiveTimeStepping 
    //---------------------
    
    class AdaptiveTimeStepping
    {
    public:    
        //! \brief contructor taking parameter object
        AdaptiveTimeStepping( const parameter::ParameterGroup& param );

        /** \brief  step method that acts like the solver::step method
                    in a sub cycle of time steps 
            
            \param  solver      solver object that must implement a method step( dt, state, well_state )
            \param  state       current state of the solution variables 
            \param  well_state  additional well state object
            \param  time        current simulation time
            \param  timestep    current time step length that is to be sub cycled 
        */          
        template <class Solver, class State, class WellState>
        void step( Solver& solver, State& state, WellState& well_state,
                   const double time, const double timestep );

    protected:
        typedef std::unique_ptr< TimeStepControlInterface > TimeStepControlType;

        TimeStepControlType timeStepControl_; //!< time step control object
        const double initial_fraction_;       //!< fraction to take as a guess for initial time interval
        const double restart_factor_;         //!< factor to multiply time step with when solver fails to converge
        const int solver_restart_max_;        //!< how many restart of solver are allowed
        const bool solver_verbose_;           //!< solver verbosity 
        const bool timestep_verbose_;         //!< timestep verbosity 
        double last_timestep_;                //!< size of last timestep
    };
}

#include <opm/core/simulator/AdaptiveTimeStepping_impl.hpp>
#endif
