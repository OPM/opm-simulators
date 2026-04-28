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

#include <opm/common/OpmLog/OpmLog.hpp>

#include <fmt/format.h>

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
    Parameters::Register<Parameters::Slave>
        ("Specify if the simulation is a slave simulation in a master-slave simulation");
    Parameters::Hide<Parameters::Slave>();

}

void logTuning(const Tuning& tuning)
{
    const auto msg = fmt::format(fmt::runtime("Tuning values: "
                                 "MB: {:.2e}, CNV: {:.2e}, NEWTMN: {}, NEWTMX: {}"),
                                 tuning.TRGMBE, tuning.TRGCNV, tuning.NEWTMN, tuning.NEWTMX);
    OpmLog::debug(msg);
    if (tuning.TRGTTE_has_value) {
        OpmLog::warning("Tuning item 2-1 (TRGTTE) is not supported.");
    }
    if (tuning.TRGLCV_has_value) {
        OpmLog::warning("Tuning item 2-4 (TRGLCV) is not supported.");
    }
    if (tuning.XXXTTE_has_value) {
        OpmLog::warning("Tuning item 2-5 (XXXTTE) is not supported.");
    }
    if (tuning.XXXLCV_has_value) {
        OpmLog::warning("Tuning item 2-8 (XXXLCV) is not supported.");
    }
    if (tuning.XXXWFL_has_value) {
        OpmLog::warning("Tuning item 2-9 (XXXWFL) is not supported.");
    }
    if (tuning.TRGFIP_has_value) {
        OpmLog::warning("Tuning item 2-10 (TRGFIP) is not supported.");
    }
    if (tuning.TRGSFT_has_value) {
        OpmLog::warning("Tuning item 2-11 (TRGSFT) is not supported.");
    }
    if (tuning.THIONX_has_value) {
        OpmLog::warning("Tuning item 2-12 (THIONX) is not supported.");
    }
    if (tuning.TRWGHT_has_value) {
        OpmLog::warning("Tuning item 2-13 (TRWGHT) is not supported.");
    }
    if (tuning.LITMAX_has_value) {
        OpmLog::warning("Tuning item 3-3 (LITMAX) is not supported.");
    }
    if (tuning.LITMIN_has_value) {
        OpmLog::warning("Tuning item 3-4 (LITMIN) is not supported.");
    }
    if (tuning.MXWSIT_has_value) {
        OpmLog::warning("Tuning item 3-5 (MXWSIT) is not supported.");
    }
    if (tuning.MXWPIT_has_value) {
        OpmLog::warning("Tuning item 3-6 (MXWPIT) is not supported.");
    }
    if (tuning.DDPLIM_has_value) {
        OpmLog::warning("Tuning item 3-7 (DDPLIM) is not supported.");
    }
    if (tuning.DDSLIM_has_value) {
        OpmLog::warning("Tuning item 3-8 (DDSLIM) is not supported.");
    }
    if (tuning.TRGDPR_has_value) {
        OpmLog::warning("Tuning item 3-9 (TRGDPR) is not supported.");
    }
    if (tuning.XXXDPR_has_value) {
        OpmLog::warning("Tuning item 3-10 (XXXDPR) is not supported.");
    }
    if (tuning.MNWRFP_has_value) {
        OpmLog::warning("Tuning item 3-11 (MNWRFP) is not supported.");
    }
}

} // namespace Opm::detail
