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

#include "config.h"
#include "initstateequil.hh"
#include "initstateequil_impl.hh"

namespace Opm {
namespace EQUIL {
namespace DeckDependent {
#if HAVE_DUNE_FEM
using GridView = Dune::Fem::GridPart2GridViewImpl<
                                     Dune::Fem::AdaptiveLeafGridPart<
                                         Dune::CpGrid,
                                         Dune::PartitionIteratorType(4),
                                         false>>;
#else
using GridView = Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>;
#endif

using Mapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
template class InitialStateComputer<BlackOilFluidSystem<double>,
                                    Dune::CpGrid,
                                    GridView,
                                    Mapper,
                                    Dune::CartesianIndexMapper<Dune::CpGrid>>;

using MatLaw = EclMaterialLawManager<ThreePhaseMaterialTraits<double,0,1,2>>;
template InitialStateComputer<BlackOilFluidSystem<double>,
                              Dune::CpGrid,
                              GridView,
                              Mapper,
                              Dune::CartesianIndexMapper<Dune::CpGrid>>::
    InitialStateComputer(MatLaw&,
                         const EclipseState&,
                         const Dune::CpGrid&,
                         const GridView&,
                         const Dune::CartesianIndexMapper<Dune::CpGrid>&,
                         const double,
                         const bool);
#if HAVE_DUNE_ALUGRID
#if HAVE_MPI
using ALUGridComm = Dune::ALUGridMPIComm;
#else
using ALUGridComm = Dune::ALUGridNoComm;
#endif //HAVE_MPI
using ALUGrid3CN = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming, ALUGridComm>;
using ALUGridView = Dune::GridView<Dune::ALU3dLeafGridViewTraits<const ALUGrid3CN, Dune::PartitionIteratorType(4)>>;
using ALUGridMapper = Dune::MultipleCodimMultipleGeomTypeMapper<ALUGridView>;
template class InitialStateComputer<BlackOilFluidSystem<double>,
                                    ALUGrid3CN,
                                    ALUGridView,
                                    ALUGridMapper,
                                    Dune::CartesianIndexMapper<ALUGrid3CN>>;

template InitialStateComputer<BlackOilFluidSystem<double>,
                              ALUGrid3CN,
                              ALUGridView,
                              ALUGridMapper,
                              Dune::CartesianIndexMapper<ALUGrid3CN>>::
    InitialStateComputer(MatLaw&,
                         const EclipseState&,
                         const ALUGrid3CN&,
                         const ALUGridView&,
                         const Dune::CartesianIndexMapper<ALUGrid3CN>&,
                         const double,
                         const bool);
#endif //HAVE_DUNE_ALUGRID


} // namespace DeckDependent

namespace Details {
    template class PressureTable<BlackOilFluidSystem<double>,EquilReg>;
    template void verticalExtent<std::vector<int>,
                                 Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator>>(
                                 const std::vector<int>&,
                                 const std::vector<std::pair<double,double>>&,
                                 const Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator>&,
                                 std::array<double,2>&);

    using MatLaw = EclMaterialLawManager<ThreePhaseMaterialTraits<double,0,1,2>>;
    template class PhaseSaturations<MatLaw,BlackOilFluidSystem<double>,
                                    EquilReg,size_t>;

    template std::pair<double,double> cellZMinMax(const Dune::cpgrid::Entity<0>& element);
}

} // namespace EQUIL
} // namespace Opm
