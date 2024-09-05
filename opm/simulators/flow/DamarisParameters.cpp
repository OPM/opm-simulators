// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2022 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2023 Inria, Bretagne–Atlantique Research Center

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/

#include <config.h>

#include <opm/simulators/flow/DamarisParameters.hpp>

#include <opm/models/utils/parametersystem.hpp>

namespace Opm {

void registerDamarisParameters()
{
    Parameters::Register<Parameters::DamarisOutputHdfCollective>
        ("Write output via Damaris using parallel HDF5 to "
         "get single file and dataset per timestep instead "
         "of one per Damaris core with multiple datasets.");
    Parameters::Register<Parameters::DamarisSaveToHdf>
        ("Set to false to prevent output to HDF5. "
         "Uses collective output by default or "
         "set --enable-damaris-collective=false to"
         "use file per core (file per Damaris server).");
    Parameters::Register<Parameters::DamarisSaveMeshToHdf>
        ("Saves the mesh data to the HDF5 file (1st iteration only). "
         "Will set  --damaris-output-hdf-collective to false "
         "so will use file per core (file per Damaris server) output "
        "(global sizes and offset values  of mesh variables are not being provided as yet).");
    Parameters::Register<Parameters::DamarisPythonScript>
        ("Set to the path and filename of a Python script to run on "
         "Damaris server resources with access to OPM flow data.");
    Parameters::Register<Parameters::DamarisPythonParaviewScript>
        ("Set to the path and filename of a Paraview Python script "
         "to run on Paraview Catalyst (1 or 2) on Damaris server "
         "resources with access to OPM flow data.");
    Parameters::Register<Parameters::DamarisSimName>
        ("The name of the simulation to be used by Damaris. "
         "If empty (the default) then Damaris uses \"opm-sim-<random-number>\". "
         "This name is used for the Damaris HDF5 file name prefix. "
         "Make unique if writing to the same output directory.");
    Parameters::Register<Parameters::DamarisLogLevel>
        ("The log level for the Damaris logging system (boost log based). "
         "Levels are: [trace, debug, info, warning, error, fatal]. "
         "Currently debug and info are useful. ");
    Parameters::Register<Parameters::DamarisDaskFile>
        ("The name of a Dask json configuration file (if using Dask for processing).");
    Parameters::Register<Parameters::DamarisDedicatedCores>
        ("Set the number of dedicated cores (MPI processes) "
         "that should be used for Damaris processing (per node). "
         "Must divide evenly into the number of simulation ranks (client ranks).");
    Parameters::Register<Parameters::DamarisDedicatedNodes>
        ("Set the number of dedicated nodes (full nodes) "
         "that should be used for Damaris processing (per simulation). "
         "Must divide evenly into the number of simulation nodes.");
    Parameters::Register<Parameters::DamarisSharedMemorySizeBytes>
        ("Set the size of the shared memory buffer used for IPC "
         "between the simulation and the Damaris resources. "
         "Needs to hold all the variables published, possibly over "
         "multiple simulation iterations.");
    Parameters::Register<Parameters::DamarisSharedMemoryName>
        ("The name of the shared memory area to be used by Damaris for the current. "
         "If empty (the default) then Damaris uses \"opm-damaris-<random-string>\". "
         "This name should be unique if multiple simulations are running on "
         "the same node/server as it is used for the Damaris shmem name and by "
         "the Python Dask library to locate sections of variables.");
    Parameters::Register<Parameters::DamarisLimitVariables>
        ("A comma separated list of variable names that a user wants to pass "
         "through via DamarisOutput::DamarisWriter::writeOutput)() to the "
         "damaris_write() call. This can be used to limit the number of "
         "variables being passed to the Damaris plugins (Paraview, Python and HDF5)");
}

} // namespace Opm
