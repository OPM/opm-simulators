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

#include <config.h>
#include <opm/models/utils/timer.hpp>

#if HAVE_MPI
#include <mpi.h>
#endif

namespace Opm {

void Timer::TimeData::measure()
{
    // Note: On Linux -- or rather fully POSIX compliant systems -- using
    // clock_gettime() would be more accurate for the CPU time.
    realtimeData = std::chrono::high_resolution_clock::now();
    cputimeData = std::clock();
}

Timer::Timer()
{
    halt();
}

void Timer::start()
{
    isStopped_ = false;
    startTime_.measure();
}

double Timer::stop()
{
    if (!isStopped_) {
        TimeData stopTime;

        stopTime.measure();

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

void Timer::halt()
{
    isStopped_ = true;
    cpuTimeElapsed_ = 0.0;
    realTimeElapsed_ = 0.0;
}

void Timer::reset()
{
    cpuTimeElapsed_ = 0.0;
    realTimeElapsed_ = 0.0;

    startTime_.measure();
}

double Timer::realTimeElapsed() const
{
    if (isStopped_)
        return realTimeElapsed_;

    TimeData stopTime;

    stopTime.measure();

    const auto& t1 = startTime_.realtimeData;
    const auto& t2 = stopTime.realtimeData;
    std::chrono::duration<double> dt =
        std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);

    return realTimeElapsed_ + dt.count();
}

double Timer::cpuTimeElapsed() const
{
    if (isStopped_)
        return cpuTimeElapsed_;

    TimeData stopTime;

    stopTime.measure();

    const auto& t1 = startTime_.cputimeData;
    const auto& t2 = stopTime.cputimeData;

    return cpuTimeElapsed_ + static_cast<double>(t2 - t1)/CLOCKS_PER_SEC;
}

double Timer::globalCpuTimeElapsed() const
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

Timer& Timer::operator+=(const Timer& other)
{
    realTimeElapsed_ += other.realTimeElapsed();
    cpuTimeElapsed_ += other.cpuTimeElapsed();

    return *this;
}

} // namespace Opm
