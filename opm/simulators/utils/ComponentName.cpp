/*
  Copyright 2013, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014, 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2014, 2015 Statoil ASA.
  Copyright 2015 NTNU
  Copyright 2015, 2016, 2017 IRIS AS

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

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidsystems/GenericOilGasWaterFluidSystem.hpp>

#include <opm/models/blackoil/blackoilvariableandequationindices.hh>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>
#include <opm/models/ptflash/flashindices.hh>

#include <flowexperimental/comp/flow_comp.hpp>

#include "ComponentName_impl.hpp"

namespace Opm {

#include <opm/simulators/utils/InstantiationIndicesMacros.hpp>

#define INSTANTIATE_COMPONENT_NAME_TYPE_TAG(T, N, W, TAG)                  \
template class ComponentName<GenericOilGasWaterFluidSystem<T, N, W>,       \
                             FlashIndices<Properties::TTag::TAG<N, W>, 0>>;

#define INSTANTIATE_COMPONENT_NAME_PROBLEMS(T, N, W)                       \
    INSTANTIATE_COMPONENT_NAME_TYPE_TAG(T, N, W, FlowCompProblem)

#define INSTANTIATE_COMPONENT_NAME_COUNTS(T, W)                            \
    INSTANTIATE_COMPONENT_NAME_PROBLEMS(T, 2, W)                           \
    INSTANTIATE_COMPONENT_NAME_PROBLEMS(T, 3, W)                           \
    INSTANTIATE_COMPONENT_NAME_PROBLEMS(T, 4, W)                           \
    INSTANTIATE_COMPONENT_NAME_PROBLEMS(T, 5, W)                           \
    INSTANTIATE_COMPONENT_NAME_PROBLEMS(T, 6, W)                           \
    INSTANTIATE_COMPONENT_NAME_PROBLEMS(T, 7, W)                           \
    INSTANTIATE_COMPONENT_NAME_PROBLEMS(T, 8, W)                           \
    INSTANTIATE_COMPONENT_NAME_PROBLEMS(T, 9, W)                           \
    INSTANTIATE_COMPONENT_NAME_PROBLEMS(T, 10, W)

INSTANTIATE_TYPE_INDICES(ComponentName, double)
INSTANTIATE_COMPONENT_NAME_COUNTS(double, false)
INSTANTIATE_COMPONENT_NAME_COUNTS(double, true)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE_INDICES(ComponentName, float)
INSTANTIATE_COMPONENT_NAME_COUNTS(float, false)
INSTANTIATE_COMPONENT_NAME_COUNTS(float, true)
#endif

#undef INSTANTIATE_COMPONENT_NAME_COUNTS
#undef INSTANTIATE_COMPONENT_NAME_PROBLEMS
#undef INSTANTIATE_COMPONENT_NAME_TYPE_TAG

} // namespace Opm
