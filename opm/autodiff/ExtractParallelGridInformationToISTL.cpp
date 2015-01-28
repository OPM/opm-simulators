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
#include <opm/core/linalg/ParallelIstlInformation.hpp>
#include <dune/common/version.hh>

 namespace Opm
{
#if defined(HAVE_DUNE_CORNERPOINT) && DUNE_VERSION_NEWER(DUNE_GRID, 2, 3)
#if defined(HAVE_MPI) && defined(HAVE_DUNE_ISTL)
// Extracts the information about the data decomposition from the grid for dune-istl
struct NullDeleter
{
    void operator()(void*)
    {}
};

void extractParallelGridInformationToISTL(boost::any& anyComm, const Dune::CpGrid& grid)
{
    if(grid.comm().size()>1)
    {
        // this is a parallel run with distributed data.
        NullDeleter null_delete;
        Dune::CpGrid& mgrid=const_cast<Dune::CpGrid&>(grid);
        Dune::CpGrid::ParallelIndexSet* idx=&(mgrid.getCellIndexSet());
        Dune::CpGrid::RemoteIndices* ridx=&(mgrid.getCellRemoteIndices());
        anyComm=boost::any(Opm::ParallelISTLInformation(std::shared_ptr<Dune::CpGrid::ParallelIndexSet>(idx, null_delete),
                                                        std::shared_ptr<Dune::CpGrid::RemoteIndices>(ridx, null_delete),
                                                        grid.comm()));
    }
}
#else
// Missing support for MPI or dune-istl -> do nothing.
void extractParallelGridInformationToISTL(boost::any&, const Dune::CpGrid&)
{}
#endif //defined(HAVE_MPI) && defined(HAVE_DUNE_ISTL)
#endif //defined(HAVE_DUNE_CORNERPOINT)
} // end namespace Opm
