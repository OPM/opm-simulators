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

#include <config.h>
#include <opm/simulators/timestepping/SimulatorTimerInterface.hpp>

#include <boost/date_time/posix_time/conversion.hpp>

namespace Opm
{

boost::posix_time::ptime SimulatorTimerInterface::currentDateTime() const
{
    return startDateTime() + boost::posix_time::seconds( (int) simulationTimeElapsed());
    //boost::posix_time::ptime(startDate()) +  boost::posix_time::seconds( (int) simulationTimeElapsed());
}

time_t SimulatorTimerInterface::currentPosixTime() const
{
    tm t = boost::posix_time::to_tm(currentDateTime());
    return std::mktime(&t);
}

} // namespace Opm
