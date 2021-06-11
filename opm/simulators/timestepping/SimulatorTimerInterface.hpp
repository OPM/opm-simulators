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

#ifndef OPM_SIMULATORTIMERINTERFACE_HEADER_INCLUDED
#define OPM_SIMULATORTIMERINTERFACE_HEADER_INCLUDED

#include <memory>

namespace boost { namespace posix_time { class ptime; } }

namespace Opm
{

    class ParameterGroup;

    /// Interface class for SimulatorTimer objects, to be improved.
    class SimulatorTimerInterface
    {
    protected:
        /// Default constructor, protected to not allow explicit instances of this class.
        SimulatorTimerInterface() {}

    public:
        /// destructor
        virtual ~SimulatorTimerInterface() {}

        /// Current step number. This is the number of timesteps that
        /// has been completed from the start of the run. The time
        /// after initialization but before the simulation has started
        /// is timestep number zero.
        virtual int currentStepNum() const = 0;

        /// Current report step number. This might differ from currentStepNum in case of sub stepping
        virtual int reportStepNum() const { return currentStepNum(); }

        /// Current step length. This is the length of the step
        /// the simulator will take in the next iteration.
        ///
        /// @note if done(), it is an error to call currentStepLength().
        virtual double currentStepLength() const = 0;

        /// Previous step length. This is the length of the step that
        /// was taken to arrive at this time.
        ///
        /// @note if no increments have been done (i.e. the timer is
        /// still in its constructed state and currentStepNum() == 0),
        /// it is an error to call stepLengthTaken().
        virtual double stepLengthTaken () const = 0;

        /// Previous report step length. This is the length of the step that
        /// was taken to arrive at this report step time.
        ///
        /// @note if no increments have been done (i.e. the timer is
        /// still in its constructed state and reportStepNum() == 0),
        /// it is an error to call stepLengthTaken().
        virtual double reportStepLengthTaken () const { return stepLengthTaken(); }

        /// Time elapsed since the start of the simulation until the
        /// beginning of the current time step [s].
        virtual double simulationTimeElapsed() const = 0;

        /// advance time by currentStepLength
        virtual void advance() = 0 ;

        /// Return true if timer indicates that simulation of timer interval is finished
        virtual bool done() const = 0;

        /// Whether the current step is the first step.
        virtual bool initialStep() const = 0;

        /// Return start date of simulation
        virtual boost::posix_time::ptime startDateTime() const = 0;

        /// Return the current time as a posix time object.
        virtual boost::posix_time::ptime currentDateTime() const;

        /// Time elapsed since the start of the POSIX epoch (Jan 1st,
        /// 1970) until the current time step begins [s].
        virtual time_t currentPosixTime() const;

        /// Return true if last time step failed
        virtual bool lastStepFailed() const = 0;

        /// return copy of current timer instance
        virtual std::unique_ptr< SimulatorTimerInterface > clone () const = 0;
    };


} // namespace Opm

#endif // OPM_SIMULATORTIMER_HEADER_INCLUDED
