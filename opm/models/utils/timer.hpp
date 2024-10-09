// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \copydoc Opm::Timer
 */
#ifndef OPM_TIMER_HPP
#define OPM_TIMER_HPP

#include <chrono>

namespace Opm {

/*!
 * \ingroup Common
 *
 * \brief Provides an encapsulation to measure the system time
 *
 * This means the wall clock time used by the simulation, the CPU time
 * used by all threads of a single process and the CPU time used by
 * the overall simulation. (i.e., the time used by all threads of all
 * involved processes.)
 */
class Timer
{
    struct TimeData
    {
        std::chrono::high_resolution_clock::time_point realtimeData;
        std::clock_t cputimeData;

        // measure the current time and put it into the object.
        void measure();
    };

public:
    Timer();

    /*!
     * \brief Start counting the time resources used by the simulation.
     */
    void start();

    /*!
     * \brief Stop counting the time resources.
     *
     * Returns the wall clock time the timer was active.
     */
    double stop();

    /*!
     * \brief Stop the measurement reset all timing values
     */
    void halt();

    /*!
     * \brief Make the current point in time t=0 but do not change the status of the timer.
     */
    void reset();

    /*!
     * \brief Return the  real time [s] elapsed  during the periods the  timer was active
     *        since the last reset.
     */
    double realTimeElapsed() const;

    /*!
     * \brief This is an alias for realTimeElapsed()
     *
     * Its main purpose is to make the API of the class a superset of Dune::Timer
     */
    double elapsed() const
    { return realTimeElapsed(); }

    /*!
     * \brief Return the CPU time [s] used by all threads of the local process for the
     *        periods the timer was active
     */
    double cpuTimeElapsed() const;

    /*!
     * \brief Return the CPU time [s] used by all threads of the all processes of program
     *
     * The value returned only differs from cpuTimeElapsed() if MPI is used.
     */
    double globalCpuTimeElapsed() const;

    /*!
     * \brief Adds the time of another timer to the current one
     */
    Timer& operator+=(const Timer& other);

private:
    bool isStopped_;
    double cpuTimeElapsed_;
    double realTimeElapsed_;
    TimeData startTime_;
};

} // namespace Opm

#endif // OPM_TIMER_HPP
