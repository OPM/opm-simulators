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

#include <dune/grid/common/defaultgridview.hh>
#include <dune/grid/common/gridview.hh>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/LookUpData.hh>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidsystems/GenericOilGasWaterFluidSystem.hpp>

#include <opm/simulators/flow/FlowGenericProblem_impl.hpp>

#if HAVE_DUNE_FEM
#include <dune/common/version.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <opm/simulators/flow/FemCpGridCompat.hpp>
#endif // HAVE_DUNE_FEM

namespace Opm {

using CpLeafGridView = Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>;

#define INSTANTIATE_TYPE(T)                                               \
    template class FlowGenericProblem<                                    \
            CpLeafGridView,                                     \
                      BlackOilFluidSystem<T, BlackOilDefaultFluidSystemIndices>>;

#define INSTANTIATE_COMP_TYPE(T, N, W)                                    \
  template class FlowGenericProblem<                                    \
            CpLeafGridView,                                     \
            GenericOilGasWaterFluidSystem<T, N, W>>;

INSTANTIATE_TYPE(double)
INSTANTIATE_COMP_TYPE(double, 2, false)
INSTANTIATE_COMP_TYPE(double, 2, true)
INSTANTIATE_COMP_TYPE(double, 3, false)
INSTANTIATE_COMP_TYPE(double, 3, true)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
INSTANTIATE_COMP_TYPE(float, 2, false)
INSTANTIATE_COMP_TYPE(float, 2, true)
INSTANTIATE_COMP_TYPE(float, 3, false)
INSTANTIATE_COMP_TYPE(float, 3, true)
#endif

#if HAVE_DUNE_FEM
using GV = Dune::Fem::AdaptiveLeafGridPart<Dune::CpGrid,
                                           (Dune::PartitionIteratorType)4,
                                           false>;
template class FlowGenericProblem<GV,
                                  BlackOilFluidSystem<double, BlackOilDefaultFluidSystemIndices>>;
template class FlowGenericProblem<GV,
                  GenericOilGasWaterFluidSystem<double, 2, false>>;
template class FlowGenericProblem<GV,
                  GenericOilGasWaterFluidSystem<double, 2, true>>;
template class FlowGenericProblem<GV,
                  GenericOilGasWaterFluidSystem<double, 3, false>>;
template class FlowGenericProblem<GV,
                  GenericOilGasWaterFluidSystem<double, 3, true>>;
#if FLOW_INSTANTIATE_FLOAT
template class FlowGenericProblem<GV,
                                  BlackOilFluidSystem<float, BlackOilDefaultFluidSystemIndices>>;
template class FlowGenericProblem<GV,
                  GenericOilGasWaterFluidSystem<float, 2, false>>;
template class FlowGenericProblem<GV,
                  GenericOilGasWaterFluidSystem<float, 2, true>>;
template class FlowGenericProblem<GV,
                  GenericOilGasWaterFluidSystem<float, 3, false>>;
template class FlowGenericProblem<GV,
                  GenericOilGasWaterFluidSystem<float, 3, true>>;
#endif

#endif // HAVE_DUNE_FEM

} // end namespace Opm
