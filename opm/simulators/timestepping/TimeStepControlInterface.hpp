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
#ifndef OPM_TIMESTEPCONTROLINTERFACE_HEADER_INCLUDED
#define OPM_TIMESTEPCONTROLINTERFACE_HEADER_INCLUDED

#include <opm/core/simulator/SimulatorState.hpp> 

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

}
#endif
