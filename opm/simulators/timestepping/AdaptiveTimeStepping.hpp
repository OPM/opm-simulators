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
#ifndef OPM_SUBSTEPPING_HEADER_INCLUDED
#define OPM_SUBSTEPPING_HEADER_INCLUDED

#include <iostream>
#include <utility>

#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/simulator/TimeStepControlInterface.hpp>

namespace Opm {


    // AdaptiveTimeStepping
    //---------------------

    class AdaptiveTimeStepping
    {
    public:
        //! \brief contructor taking parameter object
        //! \param param The parameter object
        //! \param pinfo The information about the data distribution
        //!              and communication for a parallel run.
        AdaptiveTimeStepping( const parameter::ParameterGroup& param,
                              const bool terminal_output = true );

        /** \brief  step method that acts like the solver::step method
                    in a sub cycle of time steps

            \param  timer       simulator timer providing time and timestep
            \param  solver      solver object that must implement a method step( dt, state, well_state )
            \param  state       current state of the solution variables
            \param  well_state  additional well state object
        */
        template <class Solver, class State, class WellState>
        void step( const SimulatorTimer& timer,
                   Solver& solver, State& state, WellState& well_state );

        /** \brief  step method that acts like the solver::step method
                    in a sub cycle of time steps

            \param  timer        simulator timer providing time and timestep
            \param  solver       solver object that must implement a method step( dt, state, well_state )
            \param  state        current state of the solution variables
            \param  well_state   additional well state object
            \param  outputWriter writer object to write sub steps
        */
        template <class Solver, class State, class WellState, class Output>
        void step( const SimulatorTimer& timer,
                   Solver& solver, State& state, WellState& well_state,
                   Output& outputWriter );

    protected:
        template <class Solver, class State, class WellState, class Output>
        void stepImpl( const SimulatorTimer& timer,
                       Solver& solver, State& state, WellState& well_state,
                       Output* outputWriter);

        typedef std::unique_ptr< TimeStepControlInterface > TimeStepControlType;

        TimeStepControlType timeStepControl_; //!< time step control object
        const double restart_factor_;         //!< factor to multiply time step with when solver fails to converge
        const double growth_factor_;          //!< factor to multiply time step when solver recovered from failed convergence
        const double max_growth_;             //!< factor that limits the maximum growth of a time step
        const double max_time_step_;          //!< maximal allowed time step size
        const int solver_restart_max_;        //!< how many restart of solver are allowed
        const bool solver_verbose_;           //!< solver verbosity
        const bool timestep_verbose_;         //!< timestep verbosity
        double suggested_next_timestep_;      //!< suggested size of next timestep
        bool full_timestep_initially_;        //!< beginning with the size of the time step from data file
    };
}

#include <opm/core/simulator/AdaptiveTimeStepping_impl.hpp>
#endif
