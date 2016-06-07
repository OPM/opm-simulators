/*
  Copyright 2013, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014, 2015 Statoil ASA.
  Copyright 2014, 2015 Dr. Markus Blatt - HPC-Simulation-Software & Services
  Copyright 2015 NTNU

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

#ifndef OPM_ITERATION_REPORT_HEADER_INCLUDED
#define OPM_ITERATION_REPORT_HEADER_INCLUDED

namespace Opm {
    /// Class used for reporting the outcome of a nonlinearIteration() call.
    struct IterationReport
    {
        bool failed;
        bool converged;
        int linear_iterations;
        int well_iterations;
    };
} // namespace Opm

#endif // OPM_ITERATION_REPORT_HEADER_INCLUDED
