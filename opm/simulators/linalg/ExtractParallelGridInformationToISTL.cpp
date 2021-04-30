/*
  Copyright 2015 Dr. Markus Blatt - HPC-Simulation-Software & Services
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
#include <config.h>
#include "ExtractParallelGridInformationToISTL.hpp"
#include <opm/simulators/linalg/ParallelIstlInformation.hpp>
#include <dune/common/version.hh>
#include <dune/common/shared_ptr.hh>
#include <opm/grid/CpGrid.hpp>

 namespace Opm
{
#if defined(HAVE_OPM_GRID)
#if defined(HAVE_MPI) && defined(HAVE_DUNE_ISTL)
// Extracts the information about the data decomposition from the grid for dune-istl
void extractParallelGridInformationToISTL(const Dune::CpGrid& grid, std::any& anyComm)
{
    if(grid.comm().size()>1)
    {
        // this is a parallel run with distributed data.
        Dune::CpGrid& mgrid=const_cast<Dune::CpGrid&>(grid);
        Dune::CpGrid::ParallelIndexSet& idx=mgrid.getCellIndexSet();
        Dune::CpGrid::RemoteIndices& ridx=mgrid.getCellRemoteIndices();
        anyComm=std::any(Opm::ParallelISTLInformation(Dune::stackobject_to_shared_ptr(idx),
                                                        Dune::stackobject_to_shared_ptr(ridx),
                                                        grid.comm()));
    }
}
#else
// Missing support for MPI or dune-istl -> do nothing.
void extractParallelGridInformationToISTL(const Dune::CpGrid&, std::any&)
{}
#endif //defined(HAVE_MPI) && defined(HAVE_DUNE_ISTL)
#endif //defined(HAVE_OPM_GRID)
} // end namespace Opm
