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

#include <ebos/eclgenericproblem.hh>
#include <ebos/eclgenericproblem_impl.hh>

#include <dune/grid/common/defaultgridview.hh>
#include <dune/grid/common/gridview.hh>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/LookUpData.hh>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/common/gridpart2gridview.hh>
#include <opm/simulators/flow/FemCpGridCompat.hpp>
#endif // HAVE_DUNE_FEM

namespace Opm {

template class EclGenericProblem<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>,
                                 BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,
                                 double>;

#if HAVE_DUNE_FEM
template class EclGenericProblem<Dune::GridView<Dune::Fem::GridPart2GridViewTraits<Dune::Fem::AdaptiveLeafGridPart<Dune::CpGrid, Dune::PartitionIteratorType(4), false>>>,
                                 BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,
                                 double>;
template class EclGenericProblem<Dune::Fem::GridPart2GridViewImpl<
                                    Dune::Fem::AdaptiveLeafGridPart<
                                        Dune::CpGrid,
                                        Dune::PartitionIteratorType(4),
                                        false> >,
                                 BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,
                                 double>;
#endif // HAVE_DUNE_FEM

} // end namespace Opm
