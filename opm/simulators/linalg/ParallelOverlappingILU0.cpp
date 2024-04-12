/*
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 Statoil AS

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

#include <opm/simulators/linalg/ParallelOverlappingILU0_impl.hpp>

#include <dune/istl/owneroverlapcopy.hh>

namespace Opm
{

#define INSTANTIATE_PAR(T, Dim, ...)                                                     \
  template class ParallelOverlappingILU0<Dune::BCRSMatrix<MatrixBlock<T,Dim,Dim>>,       \
                                         Dune::BlockVector<Dune::FieldVector<T,Dim>>,    \
                                         Dune::BlockVector<Dune::FieldVector<T,Dim>>,    \
                                         __VA_ARGS__>;                                   \
  template class ParallelOverlappingILU0<Dune::BCRSMatrix<Dune::FieldMatrix<T,Dim,Dim>>, \
                                         Dune::BlockVector<Dune::FieldVector<T,Dim>>,    \
                                         Dune::BlockVector<Dune::FieldVector<T,Dim>>,    \
                                         __VA_ARGS__>;

#if HAVE_MPI
#define INSTANTIATE(T,Dim)                                                \
    INSTANTIATE_PAR(T, Dim, Dune::Amg::SequentialInformation)             \
    INSTANTIATE_PAR(T, Dim, Dune::OwnerOverlapCopyCommunication<int,int>)
#else
#define INSTANTIATE(T,Dim) \
    INSTANTIATE_PAR(T, Dim, Dune::Amg::SequentialInformation)
#endif

#define INSTANTIATE_TYPE(T) \
    INSTANTIATE(T,1)        \
    INSTANTIATE(T,2)        \
    INSTANTIATE(T,3)        \
    INSTANTIATE(T,4)        \
    INSTANTIATE(T,5)        \
    INSTANTIATE(T,6)

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} // end namespace Opm
