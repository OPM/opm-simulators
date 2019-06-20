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

#ifndef OPM_SIMULATORREPORT_HEADER_INCLUDED
#define OPM_SIMULATORREPORT_HEADER_INCLUDED

#include <iosfwd>

namespace Opm
{

    /// A struct for returning timing data from a simulator to its caller.
    struct SimulatorReport
    {
        double pressure_time;
        double transport_time;
        double total_time;
        double solver_time;
        double assemble_time;
        double linear_solve_setup_time;
        double linear_solve_time;
        double update_time;
        double output_write_time;

        unsigned int total_well_iterations;
        unsigned int total_linearizations;
        unsigned int total_newton_iterations;
        unsigned int total_linear_iterations;

        bool converged;

        /// Default constructor initializing all times to 0.0.
        explicit SimulatorReport(bool verbose=true);
        /// Copy constructor
        SimulatorReport(const SimulatorReport&) = default;
        /// Increment this report's times by those in sr.
        void operator+=(const SimulatorReport& sr);
        /// Print a report to the given stream.
        void report(std::ostream& os);
        void reportStep(std::ostringstream& os);
        /// Print a report, leaving out the transport time.
        void reportFullyImplicit(std::ostream& os, const SimulatorReport* failedReport = nullptr);
        void reportParam(std::ostream& os);
    private:
        // Whether to print statistics to std::cout
        bool verbose_;
    };

} // namespace Opm

#endif // OPM_SIMULATORREPORT_HEADER_INCLUDED
