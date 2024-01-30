/*
  Copyright 2023 Inria, Bretagne–Atlantique Research Center

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
#include <opm/simulators/utils/GridDataOutput_impl.hpp>

#include <opm/grid/CpGrid.hpp>

#include <opm/simulators/utils/DamarisVar.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/common/gridpart2gridview.hh>
#include <ebos/femcpgridcompat.hh>
#endif // HAVE_DUNE_FEM

namespace Opm::GridDataOutput {

template<class T> using DV = DamarisOutput::DamarisVar<T>;

#define INSTANCE(part, ...) \
    template class SimMeshDataAccessor<__VA_ARGS__, part>; \
    template long SimMeshDataAccessor<__VA_ARGS__,part>:: \
        writeGridPoints<DV<double>>(DV<double>&, DV<double>&, DV<double>&) const; \
    template long SimMeshDataAccessor<__VA_ARGS__,part>:: \
        writeConnectivity<DV<int>>(DV<int>&, ConnectivityVertexOrder) const; \
    template long SimMeshDataAccessor<__VA_ARGS__,part>:: \
        writeOffsetsCells<DV<int>>(DV<int>&) const; \
    template long SimMeshDataAccessor<__VA_ARGS__,part>:: \
        writeCellTypes<DV<char>>(DV<char>&) const;

INSTANCE(1, Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>)

#if HAVE_DUNE_FEM
INSTANCE(1, Dune::Fem::GridPart2GridViewImpl<Dune::Fem::AdaptiveLeafGridPart<Dune::CpGrid, (Dune::PartitionIteratorType)4, false> >)
INSTANCE(1, Dune::GridView<Dune::Fem::GridPart2GridViewTraits<Dune::Fem::AdaptiveLeafGridPart<Dune::CpGrid, Dune::PartitionIteratorType(4), false>>>)
#endif

} // namespace Opm::GridDataOutput
