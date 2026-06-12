/*
  Copyright 2026, SINTEF Digital

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

#include <opm/simulators/utils/ComponentName.hpp>

#include <opm/material/fluidsystems/GenericOilGasWaterFluidSystem.hpp>

#include <opm/models/ptflash/flashindices.hh>

#include <flowexperimental/comp/flowexp_comp.hpp>

#include "ComponentName_impl.hpp"

namespace Opm {

#define INSTANTIATE_FLOWEXP_COMPONENT_NAME(T, N, W)                        \
template class ComponentName<GenericOilGasWaterFluidSystem<T, N, W>,       \
                             FlashIndices<Properties::TTag::FlowExpCompProblem<N, W>, 0>>;

#define INSTANTIATE_FLOWEXP_COMPONENT_NAME_COUNTS(T, W)                    \
    INSTANTIATE_FLOWEXP_COMPONENT_NAME(T, 2, W)                            \
    INSTANTIATE_FLOWEXP_COMPONENT_NAME(T, 3, W)                            \
    INSTANTIATE_FLOWEXP_COMPONENT_NAME(T, 4, W)                            \
    INSTANTIATE_FLOWEXP_COMPONENT_NAME(T, 5, W)                            \
    INSTANTIATE_FLOWEXP_COMPONENT_NAME(T, 6, W)                            \
    INSTANTIATE_FLOWEXP_COMPONENT_NAME(T, 7, W)                            \
    INSTANTIATE_FLOWEXP_COMPONENT_NAME(T, 8, W)                            \
    INSTANTIATE_FLOWEXP_COMPONENT_NAME(T, 9, W)                            \
    INSTANTIATE_FLOWEXP_COMPONENT_NAME(T, 10, W)

INSTANTIATE_FLOWEXP_COMPONENT_NAME_COUNTS(double, false)
INSTANTIATE_FLOWEXP_COMPONENT_NAME_COUNTS(double, true)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_FLOWEXP_COMPONENT_NAME_COUNTS(float, false)
INSTANTIATE_FLOWEXP_COMPONENT_NAME_COUNTS(float, true)
#endif

#undef INSTANTIATE_FLOWEXP_COMPONENT_NAME_COUNTS
#undef INSTANTIATE_FLOWEXP_COMPONENT_NAME

} // namespace Opm
