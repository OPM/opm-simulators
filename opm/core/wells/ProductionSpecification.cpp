#include "config.h"
#include <opm/core/wells/ProductionSpecification.hpp>

namespace Opm
{

    ProductionSpecification::ProductionSpecification()
    : 
      control_mode_(NONE),
      procedure_(NONE_P),
      oil_max_rate_(-1e100),
      water_max_rate_(-1e100),
      gas_max_rate_(-1e100),
      liquid_max_rate_(-1e100),
      reservoir_flow_max_rate_(-1e100),
      BHP_limit_(-1e100),
      guide_rate_(1.0),
      guide_rate_type_(NONE_GRT)      
    {
    }

}
