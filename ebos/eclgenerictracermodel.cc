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
#include "eclgenerictracermodel_impl.hh"

namespace Opm {
#if HAVE_DUNE_FEM
template class EclGenericTracerModel<Dune::CpGrid,
                                     Dune::GridView<Dune::Fem::GridPart2GridViewTraits<Dune::Fem::AdaptiveLeafGridPart<Dune::CpGrid, Dune::PartitionIteratorType(4), false>>>,
                                     Dune::MultipleCodimMultipleGeomTypeMapper<Dune::GridView<Dune::Fem::GridPart2GridViewTraits<Dune::Fem::AdaptiveLeafGridPart<Dune::CpGrid, Dune::PartitionIteratorType(4), false>>>>,
                                     Opm::EcfvStencil<double,Dune::GridView<Dune::Fem::GridPart2GridViewTraits<Dune::Fem::AdaptiveLeafGridPart<Dune::CpGrid, Dune::PartitionIteratorType(4), false>>>,false,false>,
                                     double>;
template class EclGenericTracerModel<Dune::CpGrid,
                                     Dune::Fem::GridPart2GridViewImpl<Dune::Fem::AdaptiveLeafGridPart<Dune::CpGrid, (Dune::PartitionIteratorType)4, false> >,
                                     Dune::MultipleCodimMultipleGeomTypeMapper<
                                         Dune::Fem::GridPart2GridViewImpl<
                                             Dune::Fem::AdaptiveLeafGridPart<Dune::CpGrid, Dune::PartitionIteratorType(4), false> > >,
                                     Opm::EcfvStencil<double, Dune::Fem::GridPart2GridViewImpl<
                                                                  Dune::Fem::AdaptiveLeafGridPart<Dune::CpGrid, Dune::PartitionIteratorType(4), false> >,
                                                      false, false>,
                                     double>;
#if HAVE_DUNE_ALUGRID

#if HAVE_MPI
    using ALUGrid3CN = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming, Dune::ALUGridMPIComm>;
#else
    using ALUGrid3CN = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming, Dune::ALUGridNoComm>;
#endif //HAVE_MPI
                                    
template class EclGenericTracerModel<ALUGrid3CN, Dune::GridView<Dune::Fem::GridPart2GridViewTraits<Dune::Fem::AdaptiveLeafGridPart<ALUGrid3CN, Dune::PartitionIteratorType(4), false>>>, Dune::MultipleCodimMultipleGeomTypeMapper<Dune::GridView<Dune::Fem::GridPart2GridViewTraits<Dune::Fem::AdaptiveLeafGridPart<ALUGrid3CN, Dune::PartitionIteratorType(4), false>>>>, Opm::EcfvStencil<double,Dune::GridView<Dune::Fem::GridPart2GridViewTraits<Dune::Fem::AdaptiveLeafGridPart<ALUGrid3CN, Dune::PartitionIteratorType(4), false>>>,false,false>,
                                     double>;

template class EclGenericTracerModel<ALUGrid3CN,
Dune::Fem::GridPart2GridViewImpl<Dune::Fem::AdaptiveLeafGridPart<ALUGrid3CN, Dune::PartitionIteratorType(4), false> >,
                                     Dune::MultipleCodimMultipleGeomTypeMapper<
                                         Dune::Fem::GridPart2GridViewImpl<
Dune::Fem::AdaptiveLeafGridPart<ALUGrid3CN, Dune::PartitionIteratorType(4), false> > >,
                                     Opm::EcfvStencil<double, Dune::Fem::GridPart2GridViewImpl<
                                                                  Dune::Fem::AdaptiveLeafGridPart<ALUGrid3CN, Dune::PartitionIteratorType(4), false> >,
                                                      false, false>,
                                     double>;                                     
#endif //HAVE_DUNE_ALUGRID                                     
#else
template class EclGenericTracerModel<Dune::CpGrid,
                                     Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>,
                                     Dune::MultipleCodimMultipleGeomTypeMapper<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>>,
                                     Opm::EcfvStencil<double,Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>,false,false>,
                                     double>;
#if HAVE_DUNE_ALUGRID
#if HAVE_MPI
    using ALUGrid3CN = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming, Dune::ALUGridMPIComm>;
#else
    using ALUGrid3CN = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming, Dune::ALUGridNoComm>;
#endif //HAVE_MPI

template class EclGenericTracerModel<ALUGrid3CN,
                                     Dune::GridView<Dune::ALU3dLeafGridViewTraits<const ALUGrid3CN, Dune::PartitionIteratorType(4)>>,
                                     Dune::MultipleCodimMultipleGeomTypeMapper<Dune::GridView<Dune::ALU3dLeafGridViewTraits<const ALUGrid3CN,
                                         Dune::PartitionIteratorType(4)>>>,
                                     Opm::EcfvStencil<double,Dune::GridView<Dune::ALU3dLeafGridViewTraits<const ALUGrid3CN,
                                        Dune::PartitionIteratorType(4)>>,false,false>,
                                     double>;
#endif //HAVE_DUNE_ALUGRID
#endif //HAVE_DUNE_FEM

template class EclGenericTracerModel<Dune::PolyhedralGrid<3,3,double>,
                                     Dune::GridView<Dune::PolyhedralGridViewTraits<3,3,double,Dune::PartitionIteratorType(4)>>,
                                     Dune::MultipleCodimMultipleGeomTypeMapper<Dune::GridView<Dune::PolyhedralGridViewTraits<3,3,double,Dune::PartitionIteratorType(4)>>>,
                                     Opm::EcfvStencil<double, Dune::GridView<Dune::PolyhedralGridViewTraits<3,3,double,Dune::PartitionIteratorType(4)>>,false,false>,
                                     double>;

} // namespace Opm
