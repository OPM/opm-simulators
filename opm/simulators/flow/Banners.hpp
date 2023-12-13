/*
  Copyright 2013, 2014, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 IRIS AS
  Copyright 2014 STATOIL ASA.

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
#ifndef OPM_FLOW_BANNERS_HEADER_INCLUDED
#define OPM_FLOW_BANNERS_HEADER_INCLUDED

#include <string>
#include <string_view>

namespace Opm {

struct SimulatorReport;
struct SimulatorReportSingle;

// Print an ASCII-art header to the PRT and DEBUG files.
void printPRTHeader(const int nprocs, const int nthreads,
                    const std::string& parameters,
                    std::string_view moduleVersion,
                    std::string_view compileTimestamp);

// Print flow application banner.
void printFlowBanner(int nprocs, int threads, std::string_view moduleVersionName);

// Print flow application trailer.
void printFlowTrailer(int nprocs, int nthreads,
                      const SimulatorReport& report,
                      const SimulatorReportSingle& localsolves_report);

} // namespace Opm

#endif // OPM_FLOW_BANNERS_EBOS_HEADER_INCLUDED
