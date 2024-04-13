#include <config.h>

#include <opm/simulators/linalg/PreconditionerFactory_impl.hpp>

namespace Opm {

INSTANTIATE_PF(double,5)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_PF(float,5)
#endif

}
