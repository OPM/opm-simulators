/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#include <opm/core/simulator/SimulatorReport.hpp>
#include <iostream>

namespace Opm
{
    SimulatorReport::SimulatorReport()
        : pressure_time(0.0),
          transport_time(0.0),
          total_time(0.0)
    {
    }

    void SimulatorReport::operator+=(const SimulatorReport& sr)
    {
        pressure_time += sr.pressure_time;
        transport_time += sr.transport_time;
        total_time += sr.total_time;
    }

    void SimulatorReport::report(std::ostream& os)
    {
        os << "Total time taken: " << total_time
           << "\n  Pressure time:  " << pressure_time
           << "\n  Transport time: " << transport_time << std::endl;
    }
    void SimulatorReport::reportParam(std::ostream& os)
    {
        os << "/timing/total_time=" << total_time
           << "\n/timing/pressure/total_time=" << pressure_time
           << "\n/timing/transport/total_time=" << transport_time << std::endl;
    }


} // namespace Opm
