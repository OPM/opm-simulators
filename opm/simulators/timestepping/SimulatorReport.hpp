/*
  Copyright 2012, 2020 SINTEF Digital, Mathematics and Cybernetics.

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

#include <cassert>
#include <cstdlib>
#include <iosfwd>
#include <limits>
#include <vector>

namespace Opm
{

    /// A struct for returning timing data from a simulator to its caller.
    struct SimulatorReportSingle
    {
        double pressure_time = 0.0;
        double transport_time = 0.0;
        double total_time = 0.0;
        double solver_time = 0.0;
        double assemble_time = 0.0;
        double pre_post_time = 0.0;
        double assemble_time_well = 0.0;
        double linear_solve_setup_time = 0.0;
        double linear_solve_time = 0.0;
        double update_time = 0.0;
        double output_write_time = 0.0;

        unsigned int total_well_iterations = 0;
        unsigned int total_linearizations = 0;
        unsigned int total_newton_iterations = 0;
        unsigned int total_linear_iterations = 0;
        unsigned int min_linear_iterations = std::numeric_limits<unsigned int>::max();
        unsigned int max_linear_iterations = 0;

        bool converged = false;
        bool well_group_control_changed = false;
        int exit_status = EXIT_SUCCESS;

        double global_time = 0.0;
        double timestep_length = 0.0;

        static SimulatorReportSingle serializationTestObject();

        bool operator==(const SimulatorReportSingle&) const;
        /// Increment this report's times by those in sr.
        void operator+=(const SimulatorReportSingle& sr);
        /// Print a report suitable for a single simulation step.
        void reportStep(std::ostream& os) const;
        /// Print a report suitable for the end of a fully implicit case, leaving out the pressure/transport time.
        void reportFullyImplicit(std::ostream& os, const SimulatorReportSingle* failedReport = nullptr) const;

        template<class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(pressure_time);
            serializer(transport_time);
            serializer(total_time);
            serializer(solver_time);
            serializer(assemble_time);
            serializer(pre_post_time);
            serializer(assemble_time_well);
            serializer(linear_solve_setup_time);
            serializer(linear_solve_time);
            serializer(update_time);
            serializer(output_write_time);
            serializer(total_well_iterations);
            serializer(total_linearizations);
            serializer(total_newton_iterations);
            serializer(total_linear_iterations);
            serializer(min_linear_iterations);
            serializer(max_linear_iterations);
            serializer(converged);
            serializer(well_group_control_changed);
            serializer(exit_status);
            serializer(global_time);
            serializer(timestep_length);
        }
    };

    struct SimulatorReport
    {
        SimulatorReportSingle success;
        SimulatorReportSingle failure;
        std::vector<SimulatorReportSingle> stepreports;

        static SimulatorReport serializationTestObject();

        bool operator==(const SimulatorReport&) const;
        void operator+=(const SimulatorReportSingle& sr);
        void operator+=(const SimulatorReport& sr);
        void reportFullyImplicit(std::ostream& os) const;
        void fullReports(std::ostream& os) const;

        template<class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(success);
            serializer(failure);
            serializer(stepreports);
        }
    };

    } // namespace Opm

#endif // OPM_SIMULATORREPORT_HEADER_INCLUDED
