#include <opm/core/ProductionSpecification.hpp>

namespace Opm
{

    ProductionSpecification::ProductionSpecification()
    : 
      control_mode_(NONE_CM),
      procedure_(NONE_P),
      oil_max_rate_(1e100),
      water_production_target_(1e100),
      fluid_volume_max_rate_(1e100),
      BHP_limit_(1e100),
            guide_rate_(1e100)
    {
    }

}
