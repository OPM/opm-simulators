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
#ifndef OPM_FLOW_UTILS_HEADER_INCLUDED
#define OPM_FLOW_UTILS_HEADER_INCLUDED

#include <functional>
#include <set>
#include <string_view>

namespace Opm { struct SimulatorReport; }

namespace Opm::detail {

void checkAllMPIProcesses();

void mergeParallelLogFiles(std::string_view output_dir,
                           std::string_view deck_filename,
                           bool enableLoggingFalloutWarning);

void handleExtraConvergenceOutput(SimulatorReport& report,
                                  std::string_view option,
                                  std::string_view optionName,
                                  std::string_view output_dir,
                                  std::string_view base_name);

//! \brief Hides unused runtime parameters.
template<class Scalar>
void hideUnusedParameters();

int eclPositionalParameter(std::function<void(const std::string&,
                                              const std::string&)> addKey,
                           std::set<std::string>& seenParams,
                           std::string& errorMsg,
                           const char** argv,
                           int paramIdx);

} // namespace Opm::detail

#endif // OPM_FLOW_UTILS_HEADER_INCLUDED
