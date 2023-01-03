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
#ifndef OPM_ADAPTIVESIMULATORTIMER_HEADER_INCLUDED
#define OPM_ADAPTIVESIMULATORTIMER_HEADER_INCLUDED

#include <cassert>
#include <iosfwd>
#include <vector>
#include <limits>
#include <algorithm>
#include <memory>
#include <numeric>

#include <opm/simulators/timestepping/SimulatorTimerInterface.hpp>

namespace Opm
{

    /////////////////////////////////////////////////////////
    ///
    /// \brief Simulation timer for adaptive time stepping
    ///
    /////////////////////////////////////////////////////////
    class AdaptiveSimulatorTimer : public SimulatorTimerInterface
    {
    public:
        /// \brief constructor taking a simulator timer to determine start and end time
        ///  \param timer          in case of sub stepping this is the outer timer
        ///  \param lastStepTaken  last suggested time step
        ///  \param maxTimeStep    maximum time step allowed
        AdaptiveSimulatorTimer( const SimulatorTimerInterface& timer,
                                const double lastStepTaken,
                                const double maxTimeStep = std::numeric_limits<double>::max() );

        /// \brief advance time by currentStepLength
        AdaptiveSimulatorTimer& operator++ ();

        /// \brief advance time by currentStepLength
        void advance() { this->operator++ (); }

        /// \brief provide and estimate for new time step size
        void provideTimeStepEstimate( const double dt_estimate );

        /// \brief Whether this is the first step
        bool initialStep () const;

        /// \brief \copydoc SimulationTimer::currentStepNum
        int currentStepNum () const;

        /// \brief return current report step
        int reportStepNum() const;

        /// \brief \copydoc SimulationTimer::currentStepLength
        double currentStepLength () const;

        /// \brief \copydoc SimulationTimer::totalTime
        double totalTime() const;

        /// \brief \copydoc SimulationTimer::simulationTimeElapsed
        double simulationTimeElapsed() const;

        /// \brief \copydoc SimulationTimer::done
        bool done () const;

        /// \brief return average step length used so far
        double averageStepLength() const;

        /// \brief return max step length used so far
        double maxStepLength () const;

        /// \brief return min step length used so far
        double minStepLength () const;

        /// \brief Previous step length. This is the length of the step that
        ///        was taken to arrive at this time.
        double stepLengthTaken () const;

        /// \brief report start and end time as well as used steps so far
        void report(std::ostream& os) const;

        /// \brief start date time of simulation
        boost::posix_time::ptime startDateTime() const;

        /// \brief Return true if last time step failed
        bool lastStepFailed() const {return lastStepFailed_;}

        /// \brief tell the timestepper whether timestep failed or not
        void setLastStepFailed(bool lastStepFailed) {lastStepFailed_ = lastStepFailed;}

        /// return copy of object
        virtual std::unique_ptr< SimulatorTimerInterface > clone() const;

    protected:
        std::shared_ptr<boost::posix_time::ptime> start_date_time_;
        const double start_time_;
        const double total_time_;
        const int report_step_;
        const double max_time_step_;

        double current_time_;
        double dt_;
        int current_step_;

        std::vector< double > steps_;
        bool lastStepFailed_;

    };

} // namespace Opm

#endif // OPM_SIMULATORTIMER_HEADER_INCLUDED
