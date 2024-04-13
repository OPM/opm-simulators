#include <config.h>

#include <opm/simulators/linalg/PreconditionerFactory_impl.hpp>

namespace Opm {

INSTANTIATE_PF(double,3)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_PF(float,3)
#endif

}
