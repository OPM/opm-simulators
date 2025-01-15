/*
  Copyright 2013, 2015, 2020 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2015 Andreas Lauser
  Copyright 2017 IRIS

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
#include <opm/simulators/flow/SimulatorFullyImplicitBlackoil.hpp>

namespace Opm::detail {

void registerSimulatorParameters()
{
    Parameters::Register<Parameters::EnableTerminalOutput>
        ("Print high-level information about the simulation's progress to the terminal");
    Parameters::Register<Parameters::EnableAdaptiveTimeStepping>
        ("Use adaptive time stepping between report steps");
    Parameters::Register<Parameters::OutputExtraConvergenceInfo>
        ("Provide additional convergence output "
         "files for diagnostic purposes. "
         "\"none\" gives no extra output and "
         "overrides all other options, "
         "\"steps\" generates an INFOSTEP file, "
         "\"iterations\" generates an INFOITER file. "
         "Combine options with commas, e.g., "
         "\"steps,iterations\" for multiple outputs.");
    Parameters::Register<Parameters::SaveStep>
        ("Save serialized state to .OPMRST file. "
         "Either a specific report step, \"all\" to save "
         "all report steps or \":x\" to save every x'th step."
         "Use negative values of \"x\" to keep only the last "
         "written step, or \"last\" to save every step, keeping "
         "only the last.");
    Parameters::Register<Parameters::LoadStep>
        ("Load serialized state from .OPMRST file. "
         "Either a specific report step, or 0 to load last "
         "stored report step.");
    Parameters::Register<Parameters::SaveFile>
        ("FileName for .OPMRST file used for saving serialized state. "
         "If empty, CASENAME.OPMRST is used.");
    Parameters::Hide<Parameters::SaveFile>();
    Parameters::Register<Parameters::LoadFile>
        ("FileName for .OPMRST file used to load serialized state. "
         "If empty, CASENAME.OPMRST is used.");
    Parameters::Hide<Parameters::LoadFile>();
}

} // namespace Opm::detail
