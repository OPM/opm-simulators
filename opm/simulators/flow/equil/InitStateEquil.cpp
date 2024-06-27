// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/

#include <config.h>
#include <opm/simulators/flow/equil/InitStateEquil_impl.hpp>

#include <opm/grid/CpGrid.hpp>

#if HAVE_DUNE_FEM
#include <dune/common/version.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#if !DUNE_VERSION_GTE(DUNE_FEM, 2, 9)
#include <dune/fem/gridpart/common/gridpart2gridview.hh>
#endif
#include <opm/simulators/flow/FemCpGridCompat.hpp>
#endif

namespace Opm {

template<class Scalar>
using MatLaw = EclMaterialLawManager<ThreePhaseMaterialTraits<Scalar,0,1,2>>;

namespace EQUIL {
namespace DeckDependent {

#define INSTANTIATE_COMP(T, GridView, Mapper)                                      \
    template class InitialStateComputer<BlackOilFluidSystem<T>,                    \
                                        Dune::CpGrid,                              \
                                        GridView,                                  \
                                        Mapper,                                    \
                                        Dune::CartesianIndexMapper<Dune::CpGrid>>; \
    template InitialStateComputer<BlackOilFluidSystem<T>,                          \
                                  Dune::CpGrid,                                    \
                                  GridView,                                        \
                                  Mapper,                                          \
                                  Dune::CartesianIndexMapper<Dune::CpGrid>>::      \
             InitialStateComputer(MatLaw<T>&,                                      \
                                  const EclipseState&,                             \
                                  const Dune::CpGrid&,                             \
                                  const GridView&,                                 \
                                  const Dune::CartesianIndexMapper<Dune::CpGrid>&, \
                                  const T,                                         \
                                  const int,                                       \
                                  const bool);

using GridView = Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>;
using Mapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;

INSTANTIATE_COMP(double, GridView, Mapper)
#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_COMP(float, GridView, Mapper)
#endif

#if HAVE_DUNE_FEM
#if DUNE_VERSION_GTE(DUNE_FEM, 2, 9)
using GridViewFem = Dune::Fem::AdaptiveLeafGridPart<Dune::CpGrid,
                                                   (Dune::PartitionIteratorType)4,
                                                    false>;
#else
using GridViewFem = Dune::Fem::GridPart2GridViewImpl<
                                        Dune::Fem::AdaptiveLeafGridPart<
                                            Dune::CpGrid,
                                            Dune::PartitionIteratorType(4),
                                            false>>;
#endif
using MapperFem = Dune::MultipleCodimMultipleGeomTypeMapper<GridViewFem>;

INSTANTIATE_COMP(double, GridViewFem, MapperFem)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_COMP(float, GridViewFem, MapperFem)
#endif

#endif // HAVE_DUNE_FEM

} // namespace DeckDependent

namespace Details {
#define INSTANTIATE_TYPE(T)                                                 \
    template class PressureTable<BlackOilFluidSystem<T>,EquilReg<T>>;       \
    template void verticalExtent(const std::vector<int>&,                   \
                                 const std::vector<std::pair<T,T>>&,        \
                                 const Parallel::Communication&,            \
                                 std::array<T,2>&);                         \
    template class PhaseSaturations<MatLaw<T>,BlackOilFluidSystem<T>,       \
                                    EquilReg<T>,std::size_t>;               \
    template std::pair<T,T> cellZMinMax<T>(const Dune::cpgrid::Entity<0>&);

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif
}

} // namespace EQUIL
} // namespace Opm
