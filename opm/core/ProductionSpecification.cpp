#include <opm/core/ProductionSpecification.hpp>

namespace Opm
{

    ProductionSpecification::ProductionSpecification()
    : 
      control_mode_(NONE_CM),
      procedure_(NONE_P),
      oil_max_rate_(-1.0),
      water_production_target_(-1.0),
      fluid_volume_max_rate_(-1.0),
      BHP_limit_(-1.0)
    {
    }

}
