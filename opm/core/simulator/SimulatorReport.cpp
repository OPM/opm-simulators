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

#include <iomanip>
#include <ostream>
#include <sstream>

namespace Opm
{
    SimulatorReport::SimulatorReport(bool verbose)
        : pressure_time(0.0),
          transport_time(0.0),
          total_time(0.0),
          solver_time(0.0),
          assemble_time(0.0),
	  linear_solve_setup_time(0.0),
          linear_solve_time(0.0),
          update_time(0.0),
          output_write_time(0.0),
          total_well_iterations(0),
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
	linear_solve_setup_time += sr.linear_solve_setup_time;
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

    void SimulatorReport::reportStep(std::ostringstream& ss)
    {
        if ( verbose_ )
        {
            ss << "Time step summary: ";
            if (total_well_iterations != 0) {
                ss << "well its = " << std::setw(2) << total_well_iterations << ", ";
            }
            ss << "newton its = " << std::setw(2) << total_newton_iterations << ", "
               << "linearizations = "  << std::setw(2) << total_linearizations
               << " ("  << std::fixed << std::setprecision(3) << std::setw(6) << assemble_time << " sec), "
               << "linear its = " << std::setw(3) << total_linear_iterations
               << " ("  << std::fixed << std::setprecision(3) << std::setw(6) << linear_solve_time << " sec)";
        }
    }

    void SimulatorReport::reportFullyImplicit(std::ostream& os, const SimulatorReport* failureReport)
    {
        if ( verbose_ )
        {
            double t = total_time;
            os << "Total time (seconds):         " << t;
            os << std::endl;

            t = solver_time + (failureReport ? failureReport->solver_time : 0.0);
            os << "Solver time (seconds):        " << t;
            os << std::endl;

            if (assemble_time > 0.0 || linear_solve_time > 0.0) {
                t = assemble_time + (failureReport ? failureReport->assemble_time : 0.0);
                os << " Assembly time (seconds):     " << t;
                if (failureReport) {
                    os << " (Failed: " << failureReport->assemble_time << "; "
                       << 100*failureReport->assemble_time/t << "%)";
                }
                os << std::endl;

                t = linear_solve_time + (failureReport ? failureReport->linear_solve_time : 0.0);
                os << " Linear solve time (seconds): " << t;
                if (failureReport) {
                    os << " (Failed: " << failureReport->linear_solve_time << "; "
                       << 100*failureReport->linear_solve_time/t << "%)";
                }
                os << std::endl;

		t = linear_solve_setup_time + (failureReport ? failureReport->linear_solve_setup_time : 0.0);
                os << " Linear solve setup time (seconds): " << t;
                if (failureReport) {
		  os << " (Failed: " << failureReport->linear_solve_setup_time << "; "
		     << 100*failureReport->linear_solve_setup_time/t << "%)";
                }
                os << std::endl;
		
                t = update_time + (failureReport ? failureReport->update_time : 0.0);
                os << " Update time (seconds):       " << t;
                if (failureReport) {
                    os << " (Failed: " << failureReport->update_time << "; "
                       << 100*failureReport->update_time/t << "%)";
                }
                os << std::endl;

                t = output_write_time + (failureReport ? failureReport->output_write_time : 0.0);
                os << " Output write time (seconds): " << t;
                os << std::endl;

            }

            int n = total_well_iterations + (failureReport ? failureReport->total_well_iterations : 0);
            os << "Overall Well Iterations:      " << n;
            if (failureReport) {
                os << " (Failed: " << failureReport->total_well_iterations << "; "
                   << 100.0*failureReport->total_well_iterations/n << "%)";
            }
            os << std::endl;

            n = total_linearizations + (failureReport ? failureReport->total_linearizations : 0);
            os << "Overall Linearizations:       " << n;
            if (failureReport) {
                os << " (Failed: " << failureReport->total_linearizations << "; "
                   << 100.0*failureReport->total_linearizations/n << "%)";
            }
            os << std::endl;

            n = total_newton_iterations + (failureReport ? failureReport->total_newton_iterations : 0);
            os << "Overall Newton Iterations:    " << n;
            if (failureReport) {
                os << " (Failed: " << failureReport->total_newton_iterations << "; "
                   << 100.0*failureReport->total_newton_iterations/n << "%)";
            }
            os << std::endl;

            n = total_linear_iterations + (failureReport ? failureReport->total_linear_iterations : 0);
            os << "Overall Linear Iterations:    " << n;
            if (failureReport) {
                os << " (Failed: " << failureReport->total_linear_iterations << "; "
                   << 100.0*failureReport->total_linear_iterations/n << "%)";
            }
            os << std::endl;
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
