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

#define INSTANTIATE_COMP_VARIANTS(T, N)                                   \
  INSTANTIATE_COMP_TYPE(T, N, false)                                    \
  INSTANTIATE_COMP_TYPE(T, N, true)

#define INSTANTIATE_FLOW_GENERIC_COMPONENTS(T)                            \
  INSTANTIATE_COMP_VARIANTS(T, 2)                                       \
  INSTANTIATE_COMP_VARIANTS(T, 3)                                       \
  INSTANTIATE_COMP_VARIANTS(T, 4)                                       \
  INSTANTIATE_COMP_VARIANTS(T, 5)                                       \
  INSTANTIATE_COMP_VARIANTS(T, 6)                                       \
  INSTANTIATE_COMP_VARIANTS(T, 7)                                       \
  INSTANTIATE_COMP_VARIANTS(T, 8)                                       \
  INSTANTIATE_COMP_VARIANTS(T, 9)                                       \
  INSTANTIATE_COMP_VARIANTS(T, 10)

INSTANTIATE_TYPE(double)
INSTANTIATE_FLOW_GENERIC_COMPONENTS(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
INSTANTIATE_FLOW_GENERIC_COMPONENTS(float)
#endif

#if HAVE_DUNE_FEM
using GV = Dune::Fem::AdaptiveLeafGridPart<Dune::CpGrid,
                                           (Dune::PartitionIteratorType)4,
                                           false>;

#define INSTANTIATE_DUNE_FEM_COMP_TYPE(T, N, W)                           \
template class FlowGenericProblem<GV,                                     \
                  GenericOilGasWaterFluidSystem<T, N, W>>;

#define INSTANTIATE_DUNE_FEM_COMP_VARIANTS(T, N)                          \
  INSTANTIATE_DUNE_FEM_COMP_TYPE(T, N, false)                           \
  INSTANTIATE_DUNE_FEM_COMP_TYPE(T, N, true)

#define INSTANTIATE_DUNE_FEM_FLOW_GENERIC_COMPONENTS(T)                   \
  INSTANTIATE_DUNE_FEM_COMP_VARIANTS(T, 2)                              \
  INSTANTIATE_DUNE_FEM_COMP_VARIANTS(T, 3)                              \
  INSTANTIATE_DUNE_FEM_COMP_VARIANTS(T, 4)                              \
  INSTANTIATE_DUNE_FEM_COMP_VARIANTS(T, 5)                              \
  INSTANTIATE_DUNE_FEM_COMP_VARIANTS(T, 6)                              \
  INSTANTIATE_DUNE_FEM_COMP_VARIANTS(T, 7)                              \
  INSTANTIATE_DUNE_FEM_COMP_VARIANTS(T, 8)                              \
  INSTANTIATE_DUNE_FEM_COMP_VARIANTS(T, 9)                              \
  INSTANTIATE_DUNE_FEM_COMP_VARIANTS(T, 10)

template class FlowGenericProblem<GV,
                                  BlackOilFluidSystem<double, BlackOilDefaultFluidSystemIndices>>;
INSTANTIATE_DUNE_FEM_FLOW_GENERIC_COMPONENTS(double)
#if FLOW_INSTANTIATE_FLOAT
template class FlowGenericProblem<GV,
                                  BlackOilFluidSystem<float, BlackOilDefaultFluidSystemIndices>>;
INSTANTIATE_DUNE_FEM_FLOW_GENERIC_COMPONENTS(float)
#endif

#undef INSTANTIATE_DUNE_FEM_FLOW_GENERIC_COMPONENTS
#undef INSTANTIATE_DUNE_FEM_COMP_VARIANTS
#undef INSTANTIATE_DUNE_FEM_COMP_TYPE

#endif // HAVE_DUNE_FEM

#undef INSTANTIATE_FLOW_GENERIC_COMPONENTS
#undef INSTANTIATE_COMP_VARIANTS

} // end namespace Opm
