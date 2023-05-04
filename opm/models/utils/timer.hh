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
#ifndef EWOMS_TIMER_HH
#define EWOMS_TIMER_HH

#include <chrono>

#if HAVE_MPI
#include <mpi.h>
#endif

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
    };
public:
    Timer()
    { halt(); }

    /*!
     * \brief Start counting the time resources used by the simulation.
     */
    void start()
    {
        isStopped_ = false;
        measure_(startTime_);
    }

    /*!
     * \brief Stop counting the time resources.
     *
     * Returns the wall clock time the timer was active.
     */
    double stop()
    {
        if (!isStopped_) {
            TimeData stopTime;

            measure_(stopTime);

            const auto& t1 = startTime_.realtimeData;
            const auto& t2 = stopTime.realtimeData;
            std::chrono::duration<double> dt =
                std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);

            realTimeElapsed_ += dt.count();
            cpuTimeElapsed_ +=
                static_cast<double>(stopTime.cputimeData
                                    - startTime_.cputimeData)/CLOCKS_PER_SEC;
        }

        isStopped_ = true;

        return realTimeElapsed_;
    }

    /*!
     * \brief Stop the measurement reset all timing values
     */
    void halt()
    {
        isStopped_ = true;
        cpuTimeElapsed_ = 0.0;
        realTimeElapsed_ = 0.0;
    }

    /*!
     * \brief Make the current point in time t=0 but do not change the status of the timer.
     */
    void reset()
    {
        cpuTimeElapsed_ = 0.0;
        realTimeElapsed_ = 0.0;

        measure_(startTime_);
    }

    /*!
     * \brief Return the  real time [s] elapsed  during the periods the  timer was active
     *        since the last reset.
     */
    double realTimeElapsed() const
    {
        if (isStopped_)
            return realTimeElapsed_;

        TimeData stopTime;

        measure_(stopTime);

        const auto& t1 = startTime_.realtimeData;
        const auto& t2 = stopTime.realtimeData;
        std::chrono::duration<double> dt =
            std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);

        return realTimeElapsed_ + dt.count();
    }

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
    double cpuTimeElapsed() const
    {
        if (isStopped_)
            return cpuTimeElapsed_;

        TimeData stopTime;

        measure_(stopTime);

        const auto& t1 = startTime_.cputimeData;
        const auto& t2 = stopTime.cputimeData;

        return cpuTimeElapsed_ + static_cast<double>(t2 - t1)/CLOCKS_PER_SEC;
    }

    /*!
     * \brief Return the CPU time [s] used by all threads of the all processes of program
     *
     * The value returned only differs from cpuTimeElapsed() if MPI is used.
     */
    double globalCpuTimeElapsed() const
    {
        double val = cpuTimeElapsed();
        double globalVal = val;

#if HAVE_MPI
        MPI_Reduce(&val,
                   &globalVal,
                   /*count=*/1,
                   MPI_DOUBLE,
                   MPI_SUM,
                   /*rootRank=*/0,
                   MPI_COMM_WORLD);
#endif

        return globalVal;
    }

    /*!
     * \brief Adds the time of another timer to the current one
     */
    Timer& operator+=(const Timer& other)
    {
        realTimeElapsed_ += other.realTimeElapsed();
        cpuTimeElapsed_ += other.cpuTimeElapsed();

        return *this;
    }

private:
    // measure the current time and put it into the object passed via
    // the argument.
    static void measure_(TimeData& timeData)
    {
        // Note: On Linux -- or rather fully POSIX compliant systems -- using
        // clock_gettime() would be more accurate for the CPU time.
        timeData.realtimeData = std::chrono::high_resolution_clock::now();
        timeData.cputimeData = std::clock();
    }

    bool isStopped_;
    double cpuTimeElapsed_;
    double realTimeElapsed_;
    TimeData startTime_;
};
} // namespace Opm

#endif
