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


namespace Opm
{

    ///////////////////////////////////////////////////////////////////
    ///
    ///  RelativeChangeInterface
    ///
    ///////////////////////////////////////////////////////////////////
    class RelativeChangeInterface
    {
    protected:
        RelativeChangeInterface() {}
    public:
        /// \return || u^n+1 - u^n || / || u^n+1 ||
        virtual double relativeChange() const = 0;

        /// virtual destructor (empty)
        virtual ~RelativeChangeInterface() {}
    };

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
        /// compute new time step size suggestions based on the PID controller
        /// \param dt          time step size used in the current step
        /// \param iterations  number of iterations used (linear/nonlinear)
        /// \param timeError   object to compute || u^n+1 - u^n || / || u^n+1 ||
        ///
        /// \return suggested time step size for the next step
        virtual double computeTimeStepSize( const double dt, const int iterations, const RelativeChangeInterface& relativeChange , const double simulationTimeElapsed) const = 0;

        /// virtual destructor (empty)
        virtual ~TimeStepControlInterface () {}
    };

}
#endif
