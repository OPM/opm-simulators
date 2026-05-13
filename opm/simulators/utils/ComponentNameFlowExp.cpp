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