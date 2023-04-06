/*
  Copyright 2022 SINTEF Digital, Mathematics and Cybernetics.

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

#include <string>
#include <Damaris.h>
#include <opm/simulators/utils/ParallelCommunication.hpp>

/*
    Below is the XML file for Damaris that is supported by Damaris.

    The entries in the map below will be filled by corresponding Damaris
    Keywords.
*/


namespace Opm::DamarisOutput
{
 // Initialize an XML file
 std::string initDamarisXmlFile();
 
 // Initialize Damaris by filling in th XML file and storing it in the chosen directory
 void initializeDamaris(MPI_Comm comm, int mpiRank, std::string OutputDir, bool enableDamarisOutputCollective);
 
/** 
* Set up Damaris Parameters for writing e.g., grid size and communicator to output "PRESSURE" field
* This function is used to set the size of the local arrays and to compute the offsets into a
* global view of the variable for use by HDF5 collective writing mode and for defining distributed 
* arrays in Dask. 
* N.B. As we are to use the HDF5 select-file attribute, the offsets are only used to write the 
*      'select-file' data i.e. GLOBAL_CELL_INDEX
*/
 void setupDamarisWritingPars(Parallel::Communication comm, const int n_elements_local_grid);

} // namespace Opm::DamarisOutput
