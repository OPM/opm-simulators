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
#include <iostream>
#include <vector>

#include <algorithm>
#include <numeric>

namespace Opm
{

    /////////////////////////////////////////////////////////
    ///
    /// \brief Simulation timer for adaptive time stepping
    ///
    /////////////////////////////////////////////////////////
    class AdaptiveSimulatorTimer
    {
    public:
        /// \brief constructor taking a simulator timer to determine start and end time
        ///  \param start_time     start time of timer
        ///  \param total_time     total time of timer
        ///  \param lastDt         last suggested length of time step interval
        AdaptiveSimulatorTimer( const double start_time, const double total_time, const double lastDt );

        /// \brief advance time by currentStepLength
        AdaptiveSimulatorTimer& operator++ ();

        /// \brief provide and estimate for new time step size
        void provideTimeStepEstimate( const double dt_estimate );

        /// \brief \copydoc SimulationTimer::currentStepNum
        int currentStepNum () const;

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

        /// \brief return max suggested step length
        double suggestedMax () const;

        /// \brief return average suggested step length
        double suggestedAverage () const;

        /// \brief report start and end time as well as used steps so far
        void report(std::ostream& os) const;

    protected:
        const double start_time_;
        const double total_time_;
        double current_time_;
        double dt_;
        int current_step_;

        std::vector< double > steps_;
        double suggestedMax_;
        double suggestedAverage_;

        double computeInitialTimeStep( const double lastDt ) const;
    };

} // namespace Opm

#endif // OPM_SIMULATORTIMER_HEADER_INCLUDED
