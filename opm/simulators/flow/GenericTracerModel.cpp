// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

#include <opm/simulators/flow/GenericTracerModel_impl.hpp>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#if HAVE_DUNE_FEM
#include <dune/common/version.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#if !DUNE_VERSION_GTE(DUNE_FEM, 2, 9)
#include <dune/fem/gridpart/common/gridpart2gridview.hh>
#endif
#include <opm/simulators/flow/FemCpGridCompat.hpp>
#endif // HAVE_DUNE_FEM

namespace Opm {

#define INSTANTIATE_TYPE(T)                                                                             \
    template class GenericTracerModel<Dune::CpGrid,                                                     \
                                      Dune::GridView<                                                   \
                                          Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>,               \
                                      Dune::MultipleCodimMultipleGeomTypeMapper<                        \
                                          Dune::GridView<                                               \
                                              Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>>,          \
                                      EcfvStencil<T,Dune::GridView<                                     \
                                                        Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>, \
                                                  false,false>,                                         \
                                      BlackOilFluidSystem<T,BlackOilDefaultIndexTraits>,                \
                                      T>;

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

#if HAVE_DUNE_FEM
#if DUNE_VERSION_GTE(DUNE_FEM, 2, 9)
using GV = Dune::Fem::AdaptiveLeafGridPart<Dune::CpGrid,
                                           (Dune::PartitionIteratorType)4,
                                           false>;
#define INSTANTIATE_FEM_TYPE(T)                                                      \
template class GenericTracerModel<Dune::CpGrid,                                      \
                                  GV,                                                \
                                  Dune::MultipleCodimMultipleGeomTypeMapper<GV>,     \
                                  EcfvStencil<T, GV, false, false>,                  \
                                  BlackOilFluidSystem<T,BlackOilDefaultIndexTraits>, \
                                  T>;
#else
#define INSTANTIATE_FEM_TYPE(T)                                                                     \
    template class GenericTracerModel<Dune::CpGrid,                                                 \
                                      Dune::GridView<                                               \
                                          Dune::Fem::GridPart2GridViewTraits<                       \
                                              Dune::Fem::AdaptiveLeafGridPart<                      \
                                                  Dune::CpGrid,                                     \
                                                  Dune::PartitionIteratorType(4), false>>>,         \
                                      Dune::MultipleCodimMultipleGeomTypeMapper<                    \
                                          Dune::GridView<                                           \
                                              Dune::Fem::GridPart2GridViewTraits<                   \
                                                  Dune::Fem::AdaptiveLeafGridPart<                  \
                                                      Dune::CpGrid,                                 \
                                                      Dune::PartitionIteratorType(4), false>>>>,    \
                                      EcfvStencil<T,Dune::GridView<                                 \
                                                        Dune::Fem::GridPart2GridViewTraits<         \
                                                            Dune::Fem::AdaptiveLeafGridPart<        \
                                                                Dune::CpGrid,                       \
                                                                Dune::PartitionIteratorType(4),     \
                                                                false>>>,                           \
                                                        false,false>,                               \
                                      BlackOilFluidSystem<T,BlackOilDefaultIndexTraits>,            \
                                      T>;                                                           \
    template class GenericTracerModel<Dune::CpGrid,                                                 \
                                      Dune::Fem::GridPart2GridViewImpl<                             \
                                          Dune::Fem::AdaptiveLeafGridPart<                          \
                                              Dune::CpGrid,                                         \
                                              (Dune::PartitionIteratorType)4, false> >,             \
                                      Dune::MultipleCodimMultipleGeomTypeMapper<                    \
                                         Dune::Fem::GridPart2GridViewImpl<                          \
                                             Dune::Fem::AdaptiveLeafGridPart<                       \
                                                Dune::CpGrid,                                       \
                                                Dune::PartitionIteratorType(4), false> > >,         \
                                      EcfvStencil<T, Dune::Fem::GridPart2GridViewImpl<              \
                                                        Dune::Fem::AdaptiveLeafGridPart<            \
                                                            Dune::CpGrid,                           \
                                                            Dune::PartitionIteratorType(4), false>>,\
                                                 false, false>,                                     \
                                      BlackOilFluidSystem<T,BlackOilDefaultIndexTraits>,            \
                                      T>;
#endif

INSTANTIATE_FEM_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_FEM_TYPE(float)
#endif

#endif // HAVE_DUNE_FEM

} // namespace Opm
