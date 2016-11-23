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

#include "config.h"
#include <opm/core/simulator/SimulatorReport.hpp>
#include <ostream>

namespace Opm
{
    SimulatorReport::SimulatorReport(bool verbose)
        : pressure_time(0.0),
          transport_time(0.0),
          total_time(0.0),
          total_well_iterations(0),
          solver_time(0.0),
          assemble_time(0.0),
          linear_solve_time(0.0),
          update_time(0.0),
          output_write_time(0.0),
          total_linearizations( 0 ),
          total_newton_iterations( 0 ),
          total_linear_iterations( 0 ),
          converged(false),
          verbose_(verbose)
    {
    }

    void SimulatorReport::operator+=(const SimulatorReport& sr)
    {
        pressure_time += sr.pressure_time;
        transport_time += sr.transport_time;
        linear_solve_time += sr.linear_solve_time;
        solver_time += sr.solver_time;
        assemble_time += sr.assemble_time;
        update_time += sr.update_time;
        output_write_time += sr.output_write_time;
        total_time += sr.total_time;
        total_well_iterations += sr.total_well_iterations;
        total_linearizations += sr.total_linearizations;
        total_newton_iterations += sr.total_newton_iterations;
        total_linear_iterations += sr.total_linear_iterations;
    }

    void SimulatorReport::report(std::ostream& os)
    {
        if ( verbose_ )
        {
            os << "Total time taken: " << total_time
               << "\n  Pressure time:  " << pressure_time
               << "\n  Transport time: " << transport_time
               << "\n  Overall Newton Iterations:  " << total_newton_iterations
               << "\n  Overall Linear Iterations:  " << total_linear_iterations
               << std::endl;
        }
    }

    void SimulatorReport::reportFullyImplicit(std::ostream& os)
    {
        if ( verbose_ )
        {
            os << "Total time taken (seconds):   " << total_time << "\n"
               << "Solver time (seconds):        " << solver_time << "\n";
            if (assemble_time > 0.0 || linear_solve_time > 0.0) {
                os << " Assembly time (seconds):     " << assemble_time << "\n"
                   << " Linear solve time (seconds): " << linear_solve_time  << "\n"
                   << " Update time (seconds):       " << update_time  << "\n"
                   << " Output write time (seconds): " << output_write_time  << "\n";

            }
            os << "Overall Well Iterations:      " << total_well_iterations << "\n"
               << "Overall Linearizations:       " << total_linearizations << "\n"
               << "Overall Newton Iterations:    " << total_newton_iterations << "\n"
               << "Overall Linear Iterations:    " << total_linear_iterations << "\n";
        }
    }

    void SimulatorReport::reportParam(std::ostream& os)
    {
        if ( verbose_ )
        {
            os << "/timing/total_time=" << total_time
               << "\n/timing/pressure/total_time=" << pressure_time
               << "\n/timing/transport/total_time=" << transport_time
               << "\n/timing/newton/iterations=" << total_newton_iterations
               << "\n/timing/linear/iterations=" << total_linear_iterations
               << std::endl;
        }
    }


} // namespace Opm
